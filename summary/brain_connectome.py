#!/usr/bin/env python3
"""Brain-space visualization of edge heritability (h2) from MASH / twin estimates.

Reads results/summary/mash_twin_wide.csv (produced by compare_mash_twin.R) and,
for each atlas (gordon 352, probaConns 80) and method (twin, AdjHE),
reconstructs the N x N h2 pconn (via np.triu_indices) and produces:

  (A) a cortical-surface NODE-SUMMARY plot (nilearn view_surf / plot_surf_stat_map)
      using the fsLR gifti surfaces + dlabel vertex->parcel map. No MNI coordinates
      are needed for this view. The node value is the `node-pct` percentile (default
      90th) of the h2 values over the node's incident edges.

  (B) a CONNECTOME plot (nilearn plot_connectome) using N x 3 MNI node coordinates.
      Edges are thresholded at the `edge-pct` percentile (default 90th) of |h2|.
      Coordinates are taken from --centroids (CSV) if given, otherwise derived from
      the dlabel + fsLR surfaces via wb_command (or the volume affine for subcortex).

This script must run where the atlas files (dlabel, fsLR gifti surfaces) and the
Python packages nilearn + nibabel are available (i.e. the HPC).

Example
-------
python brain_connectome.py \
  --wide results/summary/mash_twin_wide.csv \
  --atlas gordon --n 352 \
  --dlabel /path/Gordon352.dlabel.nii \
  --surf-l /path/L.white.surf.gii --surf-r /path/R.white.surf.gii \
  --centroids /path/gordon352_MNIcentroids.csv \
  --methods twin AdjHE \
  --node-pct 90 --edge-pct 90 \
  --outdir results/summary/brain
"""
import argparse
import os
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.affines import apply_affine

import nilearn.plotting as niplot


# --------------------------------------------------------------------------
# column mapping: method token -> column name in mash_twin_wide.csv
# --------------------------------------------------------------------------
def column_for(atlas, method):
    m = method.lower()
    if m == "twin":
        return "Twin_h2"
    if m == "adjhe":
        return f"h2_{atlas}_AdjHE_RE"
    if m == "adjhe_fe":
        return f"h2_{atlas}_AdjHE_FE"
    raise ValueError(f"unknown method token: {method}")


# --------------------------------------------------------------------------
# rebuild the N x N pconn from per-edge h2
# --------------------------------------------------------------------------
def rebuild_pconn(phenos, edge_h2, N):
    k = np.empty(len(phenos), dtype=int)
    for ii, p in enumerate(phenos):
        s = str(p)
        k[ii] = int(s[1:]) if s.startswith("o") else int(s)  # "o123" -> 123
    M = np.full((N, N), np.nan)
    i, j = np.triu_indices(N, 1)
    M[i[k], j[k]] = edge_h2
    M = M + M.T
    np.fill_diagonal(M, np.nan)
    return M


def node_summary(M, pct):
    N = M.shape[0]
    out = np.full(N, np.nan)
    for i in range(N):
        vals = M[i, :]
        vals = vals[~np.isnan(vals)]
        if len(vals):
            out[i] = np.percentile(vals, pct)
    return out


# --------------------------------------------------------------------------
# dlabel -> per-structure vertex->label, and parcel centroids
# --------------------------------------------------------------------------
def load_surf(path):
    img = nib.load(path)
    if hasattr(img, "darrays"):  # GiftiImage: select arrays by intent
        coords = faces = None
        for da in img.darrays:
            it = da.intent
            is_point = it in ("NIFTI_INTENT_POINTSET", "POINTSET", 1008)
            is_tri = it in ("NIFTI_INTENT_TRIANGLE", "TRIANGLE", 1009)
            if is_point:
                coords = np.asarray(da.data, dtype=float)
            elif is_tri:
                faces = np.asarray(da.data, dtype=int)
        if coords is None:  # fallback: standard surface order (first=verts, second=tris)
            arrs = [np.asarray(da.data) for da in img.darrays]
            if arrs:
                coords = arrs[0].astype(float)
                faces = arrs[1].astype(int) if len(arrs) > 1 else np.zeros((0, 0), int)
        if faces is None:
            faces = np.zeros((0, 0), dtype=int)
        return coords, faces
    data = img.agg_data()  # non-gifti fallback
    if isinstance(data, tuple):
        coords = data[0]
        faces = data[1] if len(data) > 1 else np.zeros((0, 0), dtype=int)
    else:
        coords, faces = np.asarray(data, dtype=float), np.zeros((0, 0), dtype=int)
    return np.asarray(coords, dtype=float), np.asarray(faces, dtype=int)


def vertex_mni(surf, wb_cmd):
    out = tempfile.NamedTemporaryFile(suffix=".gii", delete=False).name
    subprocess.run([wb_cmd, "-surface-coordinates-to-mni", surf, out], check=True)
    return np.asarray(nib.load(out).agg_data()[0], dtype=float)


def _brain_models(header):
    """Yield (structure, is_surface, vertex, voxel, vol_affine) for each brain
    model in a CIFTI dlabel, tolerating nibabel API differences.

    Strategy:
    1. Reconstruct surface/volume blocks from the BrainModelAxis `nvertices`
       count dict + the concatenated `vertex`/`voxel` arrays (attributes
       confirmed present across nibabel versions; reliable for cortex-only
       dlabels such as Gordon352 / probaConns80).
    2. Fallback to header.get_index_map(1).brain_models() (canonical API).
    """
    ax = header.get_axis(1)
    # Primary: nvertices + concatenated vertex/voxel arrays
    nv = getattr(ax, "nvertices", None)
    if nv:
        try:
            nv = dict(nv)
            vertex_arr = np.asarray(ax.vertex) if getattr(ax, "vertex", None) is not None else None
            voxel_arr = np.asarray(ax.voxel) if getattr(ax, "voxel", None) is not None else None
            pos_v = 0
            pos_x = 0
            for structure, count in nv.items():
                s = str(structure)
                c = int(count)
                if s.startswith("CORTEX"):
                    vertex = vertex_arr[pos_v:pos_v + c] if vertex_arr is not None else None
                    yield (s, True, vertex, None, None)
                    pos_v += c
                else:
                    voxel = voxel_arr[pos_x:pos_x + c] if voxel_arr is not None else None
                    yield (s, False, None, voxel, None)
                    pos_x += c
            return
        except Exception:
            pass
    # Fallback: canonical header-level brain-model list
    try:
        mim = header.get_index_map(1)
        gen = getattr(mim, "brain_models", None)
        bms = list(gen()) if callable(gen) else (list(gen) if gen is not None else [])
        for bm in bms:
            s = str(bm.structure)
            v = np.asarray(bm.vertex) if getattr(bm, "vertex", None) is not None else None
            vx = np.asarray(bm.voxel) if getattr(bm, "voxel", None) is not None else None
            aff = getattr(getattr(bm, "volume", None), "transform", None)
            yield (s, s.startswith("CORTEX"), v, vx, aff)
        return
    except Exception:
        pass
    raise RuntimeError("Could not extract brain models from CIFTI header")


def derive_centroids(N, dlabel, surf_l, surf_r, wb_cmd):
    img = nib.load(dlabel)
    data = img.get_fdata().astype(int).ravel()

    vtx_mni = {}
    if surf_l:
        vtx_mni["CORTEX_LEFT"] = vertex_mni(surf_l, wb_cmd)
    if surf_r:
        vtx_mni["CORTEX_RIGHT"] = vertex_mni(surf_r, wb_cmd)

    acc = defaultdict(list)
    offset = 0
    for structure, is_surface, vertex, voxel, vol_affine in _brain_models(img.header):
        n = len(vertex) if is_surface else (len(voxel) if voxel is not None else 0)
        elem = data[offset:offset + n]
        if is_surface:
            mni = vtx_mni.get(structure)
            if mni is not None:
                for li, v in zip(elem, vertex):
                    acc[int(li)].append(mni[v])
        else:
            if vol_affine is not None and voxel is not None:
                for li, vv in zip(elem, voxel):
                    acc[int(li)].append(apply_affine(vol_affine, vv))
        offset += n

    labels = list(acc.keys())
    off = 1 if (labels and max(labels) == N) else 0
    C = np.full((N, 3), np.nan)
    for li, pts in acc.items():
        idx = li - off
        if 0 <= idx < N:
            C[idx] = np.mean(pts, axis=0)
    if np.isnan(C).any():
        warnings.warn("Some parcels have no coordinate (subcortex missing or "
                      "label/parcel order mismatch). Filling with 0.")
        C = np.nan_to_num(C, nan=0.0)
    return C


def load_centroids(N, dlabel, surf_l, surf_r, wb_cmd, centroids_csv):
    if centroids_csv:
        C = np.loadtxt(centroids_csv, delimiter=",")
        if C.shape[0] != N:
            raise ValueError(f"centroids have {C.shape[0]} rows, expected {N}")
        return C
    return derive_centroids(N, dlabel, surf_l, surf_r, wb_cmd)


def dlabel_parcel_maps(N, dlabel, n_left=None, n_right=None):
    """Return (tex_left, tex_right, offset) holding the parcel index per vertex."""
    img = nib.load(dlabel)
    data = img.get_fdata().astype(int).ravel()

    off = 1 if (int(data.max()) == N) else 0  # 1-based label tables

    # learn surface vertex counts if not supplied
    if n_left is None:
        n_left = 0
        for structure, is_surface, vertex, voxel, _ in _brain_models(img.header):
            if structure == "CORTEX_LEFT" and vertex is not None:
                n_left = max(n_left, int(np.max(vertex)) + 1)
    if n_right is None:
        n_right = 0
        for structure, is_surface, vertex, voxel, _ in _brain_models(img.header):
            if structure == "CORTEX_RIGHT" and vertex is not None:
                n_right = max(n_right, int(np.max(vertex)) + 1)

    tex_left = np.full(n_left, np.nan)
    tex_right = np.full(n_right, np.nan)

    offset = 0
    for structure, is_surface, vertex, voxel, _ in _brain_models(img.header):
        n = len(vertex) if is_surface else (len(voxel) if voxel is not None else 0)
        elem = data[offset:offset + n]
        if structure == "CORTEX_LEFT":
            tex_left[vertex] = elem - off
        elif structure == "CORTEX_RIGHT":
            tex_right[vertex] = elem - off
        offset += n
    return tex_left, tex_right, off


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def plot_surface(node_h2, surf_l, surf_r, tex_left, tex_right, outdir, tag, vmin, vmax):
    # tex_left/right currently hold parcel indices; convert to node_h2 values
    val_l = np.nan_to_num(node_h2[np.nan_to_num(tex_left, nan=0).astype(int)], nan=0.0) \
        if tex_left is not None else None
    val_r = np.nan_to_num(node_h2[np.nan_to_num(tex_right, nan=0).astype(int)], nan=0.0) \
        if tex_right is not None else None

    for side, (surf, tex, val) in enumerate([(surf_l, tex_left, val_l),
                                             (surf_r, tex_right, val_r)]):
        if surf is None or tex is None:
            continue
        name = "left" if side == 0 else "right"
        png = os.path.join(outdir, f"{tag}_surface_{name}.png")
        niplot.plot_surf_stat_map(surf, val, title=f"{tag} ({name})",
                                  output_file=png, cmap="coolwarm",
                                  colorbar=True, threshold=None,
                                  vmin=vmin, vmax=vmax)
        html = os.path.join(outdir, f"{tag}_surface_{name}.html")
        try:
            niplot.view_surf(surf, val, cmap="coolwarm", vmin=vmin, vmax=vmax,
                             title=f"{tag} ({name})").save_as_html(html)
        except Exception as e:  # pragma: no cover
            warnings.warn(f"view_surf html failed for {name}: {e}")
        print(f"  wrote {png}")


def plot_connectome(M, coords, outdir, tag, thr):
    Mp = np.nan_to_num(M, nan=0.0)
    np.fill_diagonal(Mp, 0.0)
    png = os.path.join(outdir, f"{tag}_connectome.png")
    niplot.plot_connectome(Mp, coords, edge_threshold=thr,
                           node_size=8, title=f"{tag}",
                           output_file=png, cmap="coolwarm")
    print(f"  wrote {png}")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    import subprocess  # local import so --centroids-only path avoids it

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wide", required=True, help="mash_twin_wide.csv")
    ap.add_argument("--atlas", required=True, help="gordon | probaConns (or proba)")
    ap.add_argument("--n", type=int, required=True, help="number of nodes (352 or 80)")
    ap.add_argument("--dlabel", help="parcellation dlabel CIFTI")
    ap.add_argument("--surf-l", help="fsLR left white/gray gifti surface")
    ap.add_argument("--surf-r", help="fsLR right white/gray gifti surface")
    ap.add_argument("--centroids", help="optional Nx3 MNI centroid CSV")
    ap.add_argument("--wb-command", default="wb_command",
                    help="wb_command binary (for centroid derivation)")
    ap.add_argument("--methods", nargs="+", default=["twin", "AdjHE"],
                    help="twin | AdjHE | AdjHE_FE")
    ap.add_argument("--node-pct", type=float, default=90,
                    help="percentile for node summary (default 90)")
    ap.add_argument("--edge-pct", type=float, default=90,
                    help="percentile threshold for connectome edges (default 90)")
    ap.add_argument("--outdir", default="results/summary/brain")
    ap.add_argument("--no-surface", action="store_true")
    ap.add_argument("--no-connectome", action="store_true")
    args = ap.parse_args()

    atlas = "probaConns" if args.atlas.lower().startswith("proba") else args.atlas
    N = args.n
    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.wide):
        sys.exit(f"WIDE file not found: {args.wide}\n"
                 "Run summary/compare_mash_twin.R first (or let "
                 "summary/run_brain_viz.sh generate it).")
    df = pd.read_csv(args.wide)
    if "Set" not in df.columns or "Pheno" not in df.columns:
        sys.exit("wide CSV must contain Set and Pheno columns")
    sub = df[df["Set"] == atlas]
    if sub.empty:
        sys.exit(f"No rows for Set={atlas} in {args.wide}")

    # load surface geometry + dlabel parcel maps once
    surf_l = args.surf_l if (args.surf_l and not args.no_surface) else None
    surf_r = args.surf_r if (args.surf_r and not args.no_surface) else None
    tex_left = tex_right = None
    if surf_l or surf_r:
        if not args.dlabel:
            sys.exit("--dlabel is required for the surface view")
        Vl = load_surf(surf_l)[0].shape[0] if surf_l else None
        Vr = load_surf(surf_r)[0].shape[0] if surf_r else None
        tex_left, tex_right, _ = dlabel_parcel_maps(N, args.dlabel, Vl, Vr)

    # ---- compute matrices for all requested methods -----------------------
    matrices = {}
    node_vals = {}
    for method in args.methods:
        col = column_for(atlas, method)
        if col not in sub.columns:
            warnings.warn(f"column {col} missing for method {method}; skipping")
            continue
        s = sub[["Pheno", col]].dropna()
        if s.empty:
            warnings.warn(f"no data for {method} ({col}); skipping")
            continue
        M = rebuild_pconn(s["Pheno"].values, s[col].values.astype(float), N)
        matrices[method] = M
        node_vals[method] = node_summary(M, args.node_pct)
        print(f"[{atlas}_{method}] edges={int(np.isfinite(M).sum() // 2)} nodes={N}")

    if not matrices:
        sys.exit("no methods produced data; nothing to plot")

    # ---- shared colour / threshold scales for fair comparison -------------
    allv = np.concatenate([np.abs(v[~np.isnan(v)]) for v in node_vals.values()])
    vmax = float(np.nanmax(allv)) if allv.size else 1.0
    vmin = -vmax
    alledges = np.concatenate([np.abs(M[np.triu_indices(N, 1)])
                               for M in matrices.values()])
    gthr = float(np.nanpercentile(alledges, args.edge_pct)) if alledges.size else 0.0

    # ---- surface node-summary plots --------------------------------------
    if not args.no_surface and (surf_l or surf_r):
        for method, node_h2 in node_vals.items():
            plot_surface(node_h2, surf_l, surf_r, tex_left, tex_right,
                         args.outdir, f"{atlas}_{method}", vmin, vmax)

    # ---- connectome plots -------------------------------------------------
    if not args.no_connectome:
        if args.centroids is None and not (args.dlabel and (surf_l or surf_r)):
            warnings.warn("connectome skipped: need --centroids or "
                          "--dlabel + surfaces for derivation")
        else:
            try:
                coords = load_centroids(N, args.dlabel, surf_l, surf_r,
                                        args.wb_command, args.centroids)
            except Exception as e:
                warnings.warn(f"connectome skipped: centroid derivation failed: {e}")
                coords = None
            if coords is not None:
                for method, M in matrices.items():
                    plot_connectome(M, coords, args.outdir, f"{atlas}_{method}", gthr)

    print("Done. Outputs in", args.outdir)


if __name__ == "__main__":
    main()
