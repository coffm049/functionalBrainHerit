#!/usr/bin/env python3
"""Brain-space visualization of edge heritability (h2) from MASH / twin estimates.

Reads results/summary/mash_twin_wide.csv (produced by compare_mash_twin.R) and,
for each atlas (gordon 352, probaConns 80), reconstructs the N x N h2 pconn
(via np.triu_indices) and produces, for twin and AdjHE_RE:

  (A) a cortical-surface NODE-SUMMARY plot (nilearn plot_surf_stat_map) using the
       fsLR gifti surfaces + dlabel vertex->parcel map. The node value is the
       `node-pct` percentile (default 90th) of the h2 values over the node's
       incident edges.

  (B) a CIRCULAR connectome (nodes on a ring, edges drawn for |h2| > 0.3, coloured
       by sign: red positive, blue negative) at <atlas>_<method>_circular.png.

This script must run where the atlas files (dlabel, fsLR gifti surfaces) and the
Python packages nilearn + nibabel are available (i.e. the HPC).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import nibabel as nib

import nilearn.plotting as niplot


# --------------------------------------------------------------------------
# column mapping: method token -> column name in mash_twin_wide.csv
# --------------------------------------------------------------------------
def column_for(atlas, method):
    m = method.lower()
    # Set label is "probaConns" but the wide-CSV column prefix is "proba".
    prefix = "proba" if atlas.lower().startswith("proba") else atlas.lower()
    if m == "twin":
        return "Twin_h2"
    if m in ("adjhe", "adjhe_re"):
        return f"h2_{prefix}_AdjHE_RE"
    if m == "adjhe_fe":
        return f"h2_{prefix}_AdjHE_FE"
    raise ValueError(f"unknown method token: {method}")


# --------------------------------------------------------------------------
# rebuild the N x N pconn from per-edge h2
# --------------------------------------------------------------------------
def rebuild_pconn(phenos, edge_h2, N):
    k = np.empty(len(phenos), dtype=int)
    for ii, p in enumerate(phenos):
        s = str(p)
        k[ii] = int(s[1:]) if s.startswith("o") else int(s)  # "o123" -> 123
    M = np.zeros((N, N))
    i, j = np.triu_indices(N, 1)
    M[i[k], j[k]] = edge_h2
    M[j[k], i[k]] = edge_h2
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
# dlabel -> per-hemisphere surface value maps
# --------------------------------------------------------------------------
def dlabel_values(node_h2, dlabel, N):
    """Return (tex_left, tex_right): per-vertex h2 values (NaN = unassigned),
    using the same mapping as the 01-07 scripts:

        vertex_scalar_map[label_data == region] = value

    i.e. build the full greyordinate map from the dlabel's label array, then
    split it into left/right hemisphere textures via the BrainModelAxis.
    """
    img = nib.load(dlabel)
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()

    # 01-07 style: full greyordinate map, one value per parcel.
    full = np.full(label_data.shape, np.nan, dtype=float)
    for p in range(1, N + 1):
        val = node_h2[p - 1]
        if np.isfinite(val):
            full[label_data == p] = val

    ax = img.header.get_axis(1)
    tex_left = tex_right = None

    # Canonical API: iter_structures() yields (structure, slice, brainmodel),
    # with `slice` already positioned correctly in the greyordinate array
    # (handles subcortical voxels without desyncing the surface vertices).
    try:
        for structure, slc, bm in ax.iter_structures():
            s = str(structure).upper()
            v = getattr(bm, "vertex", None)
            if v is None:
                continue
            v = np.asarray(v)
            if "LEFT" in s:
                tex_left = np.full(int(np.max(v)) + 1, np.nan)
                tex_left[v] = full[slc]
            elif "RIGHT" in s:
                tex_right = np.full(int(np.max(v)) + 1, np.nan)
                tex_right[v] = full[slc]
        if tex_left is not None or tex_right is not None:
            return tex_left, tex_right
    except Exception as e:
        print(f"  [dlabel] iter_structures failed ({e}); using nvertices fallback")

    # Fallback (cortex-only dlabels): nvertices dict + concatenated vertex array.
    nv = dict(getattr(ax, "nvertices", {}) or {})
    va = np.asarray(ax.vertex) if getattr(ax, "vertex", None) is not None else None
    pos = 0
    for structure, count in nv.items():
        s = str(structure).upper()
        c = int(count)
        if "LEFT" in s and va is not None:
            blk = va[pos:pos + c]
            tex_left = np.full(int(np.max(blk)) + 1, np.nan)
            tex_left[blk] = full[pos:pos + c]
        elif "RIGHT" in s and va is not None:
            blk = va[pos:pos + c]
            tex_right = np.full(int(np.max(blk)) + 1, np.nan)
            tex_right[blk] = full[pos:pos + c]
        pos += c
    return tex_left, tex_right


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def plot_surface(val_left, val_right, surf_l, surf_r, outdir, tag, vmax):
    # val_left/right are per-vertex h2 value arrays (NaN = unassigned/background).
    vmax = max(float(vmax) if np.isfinite(vmax) else 1e-6, 1e-6)
    surfs, vals = [], []
    if surf_l is not None and val_left is not None:
        surfs.append(surf_l)
        vals.append(val_left)
    if surf_r is not None and val_right is not None:
        surfs.append(surf_r)
        vals.append(val_right)
    if not surfs:
        return
    # Prefer one figure with both hemispheres; fall back to separate L/R PNGs
    # if this nilearn version does not accept list inputs.
    try:
        png = os.path.join(outdir, f"{tag}_surface.png")
        niplot.plot_surf_stat_map(surfs, vals, title=tag, output_file=png,
                                  cmap="coolwarm", colorbar=True, threshold=None,
                                  vmin=0.0, vmax=vmax)
        print(f"  wrote {png}")
    except (ValueError, TypeError):
        for name, surf, val in (("left", surf_l, val_left),
                                ("right", surf_r, val_right)):
            if surf is None or val is None:
                continue
            png = os.path.join(outdir, f"{tag}_surface_{name}.png")
            niplot.plot_surf_stat_map(surf, val, title=f"{tag} ({name})",
                                      output_file=png, cmap="coolwarm",
                                      colorbar=True, threshold=None,
                                      vmin=0.0, vmax=vmax)
            print(f"  wrote {png}")
    assigned = int(sum(np.count_nonzero(~np.isnan(v)) for v in vals))
    if assigned:
        allv = np.concatenate([v[~np.isnan(v)] for v in vals])
        print(f"  [{tag}] surface assigned={assigned} "
              f"min={np.nanmin(allv):.3f} max={np.nanmax(allv):.3f}")
    else:
        print(f"  [{tag}] surface WARNING: no finite values")


def plot_circular(M, outdir, tag, h2_thr=0.3):
    """Circular (chord-style) connectome: nodes evenly spaced on a ring, edges
    drawn for |h2| > h2_thr, coloured by sign (red=+, blue=-)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    N = M.shape[0]
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xy = np.column_stack([np.cos(angles), np.sin(angles)])
    fig, ax = plt.subplots(figsize=(9, 9))
    iu = np.triu_indices(N, 1)
    vals = M[iu]
    mask = np.abs(vals) > h2_thr
    n_edges = int(mask.sum())
    for (a, b), v in zip(zip(iu[0][mask], iu[1][mask]), vals[mask]):
        ax.plot([xy[a, 0], xy[b, 0]], [xy[a, 1], xy[b, 1]],
                color=("red" if v > 0 else "blue"), alpha=0.25, linewidth=0.4)
    ax.scatter(xy[:, 0], xy[:, 1], s=12, c="black", zorder=5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{tag} circular (|h2|>{h2_thr}, {n_edges} edges)")
    png = os.path.join(outdir, f"{tag}_circular.png")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wide", required=True, help="mash_twin_wide.csv")
    ap.add_argument("--atlas", required=True, help="gordon | probaConns (or proba)")
    ap.add_argument("--n", type=int, required=True, help="number of nodes (352 or 80)")
    ap.add_argument("--dlabel", help="parcellation dlabel CIFTI")
    ap.add_argument("--surf-l", help="fsLR left white/gray gifti surface")
    ap.add_argument("--surf-r", help="fsLR right white/gray gifti surface")
    ap.add_argument("--node-pct", type=float, default=90,
                    help="percentile for node summary (default 90)")
    ap.add_argument("--outdir", default="results/summary/brain")
    args = ap.parse_args()

    atlas = "probaConns" if args.atlas.lower().startswith("proba") else args.atlas
    N = args.n
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.wide)
    sub = df[df["Set"] == atlas]
    if sub.empty:
        sys.exit(f"No rows for Set={atlas} in {args.wide}")

    surf_l = args.surf_l
    surf_r = args.surf_r

    # compute matrices for the two methods we care about (twin + AdjHE_RE)
    methods = ["twin", "AdjHE"]
    matrices = {}
    node_vals = {}
    for method in methods:
        col = column_for(atlas, method)
        if col not in sub.columns:
            sys.exit(f"column {col} missing for method {method}")
        s = sub[["Pheno", col]].dropna()
        if s.empty:
            sys.exit(f"no data for {method} ({col})")
        M = rebuild_pconn(s["Pheno"].values, s[col].values.astype(float), N)
        matrices[method] = M
        node_vals[method] = node_summary(M, args.node_pct)
        print(f"[{atlas}_{method}] edges={int(np.isfinite(M).sum() // 2)} nodes={N}")

    # shared colour scale for fair comparison (heritability is non-negative)
    allv = np.concatenate([np.abs(v[~np.isnan(v)]) for v in node_vals.values()])
    vmax = float(np.nanmax(allv)) if allv.size else 1.0

    # ---- cortical-surface node-summary plots -----------------------------
    if surf_l or surf_r:
        if not args.dlabel:
            sys.exit("--dlabel is required for the surface view")
        # diagnostic: confirm the dlabel maps onto the surface (AdjHE method)
        dL, dR = dlabel_values(node_vals["AdjHE"], args.dlabel, N)
        aL = int(np.count_nonzero(~np.isnan(dL))) if dL is not None else 0
        aR = int(np.count_nonzero(~np.isnan(dR))) if dR is not None else 0
        print(f"[dlabel] left assigned={aL} right assigned={aR} (N={N})")
        for method, node_h2 in node_vals.items():
            val_l, val_r = dlabel_values(node_h2, args.dlabel, N)
            plot_surface(val_l, val_r, surf_l, surf_r,
                         args.outdir, f"{atlas}_{method}", vmax)

    # ---- circular connectome (|h2| > 0.3) --------------------------------
    for method, M in matrices.items():
        plot_circular(M, args.outdir, f"{atlas}_{method}", h2_thr=0.3)

    print("Done. Outputs in", args.outdir)


if __name__ == "__main__":
    main()
