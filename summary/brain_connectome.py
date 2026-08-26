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

    Build the full greyordinate map from the dlabel's label array, then split
    it into left/right hemisphere textures via the BrainModelAxis vertex
    indices (ax.vertex is the concatenated within-structure vertex array, in
    greyordinate order; ax.nvertices gives the per-structure counts).
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
    va = getattr(ax, "vertex", None)
    if va is None:
        raise ValueError("dlabel BrainModelAxis has no vertex array")
    va = np.asarray(va, dtype=int)
    if len(va) != len(full):
        raise ValueError(
            f"dlabel vertex length ({len(va)}) != greyordinate count "
            f"({len(full)}); only cortex-only dlabels are supported"
        )
    nv = getattr(ax, "nvertices", {}) or {}

    tex_left = tex_right = None
    pos = 0
    for structure, count in nv.items():
        s = str(structure).upper()
        c = int(count)
        if c == 0:
            continue
        blk = va[pos:pos + c]
        seg = full[pos:pos + c]
        if "LEFT" in s:
            tex_left = np.full(int(blk.max()) + 1, np.nan)
            tex_left[blk] = seg
        elif "RIGHT" in s:
            tex_right = np.full(int(blk.max()) + 1, np.nan)
            tex_right[blk] = seg
        pos += c
    return tex_left, tex_right


def read_network_csv(path, network_col=None):
    """Return dict {dlabel_value(int): network_name(str)} from a CSV mapping
    each dlabel value (the row, or an explicit key column) to a network name.
    If the first column is numeric and equals 1..N or 0..N-1 it is used as an
    explicit key (1-based or 0-based respectively); otherwise rows are
    positional (row i <-> dlabel value i+1)."""
    df = pd.read_csv(path)
    cols = list(df.columns)
    if network_col is None:
        if "network" in cols:
            network_col = "network"
        elif "name" in cols:
            network_col = "name"
        else:
            network_col = cols[-1]
    key_vals = df.iloc[:, 0].to_numpy()
    net_vals = df[network_col].astype(str).to_numpy()
    result = {}
    n = len(df)
    is_num = np.issubdtype(key_vals.dtype, np.number)
    uniq = set(int(k) for k in key_vals)
    if is_num and uniq == set(range(1, n + 1)):
        for k, nv in zip(key_vals, net_vals):
            result[int(k)] = nv
    elif is_num and uniq == set(range(0, n)):
        for k, nv in zip(key_vals, net_vals):
            result[int(k) + 1] = nv
    else:
        for i, nv in enumerate(net_vals):
            result[i + 1] = nv
            result[i] = nv
    return result


def dlabel_networks(dlabel, N, csv_path=None):
    """Per-node network label (array of N strings). Uses the network CSV when
    given; otherwise falls back to the dlabel CIFTI label table (first token of
    each region name, e.g. 'SomMotA_3' -> 'SomMotA')."""
    if csv_path:
        net_by_value = read_network_csv(csv_path)
        return np.array([net_by_value.get(p, "NA") for p in range(1, N + 1)],
                        dtype=object)
    img = nib.load(dlabel)
    lt = getattr(img, "labeltable", None)
    labels = getattr(lt, "labels", {}) if lt is not None else {}
    nets = []
    for p in range(1, N + 1):
        name = None
        if p in labels:
            lab = labels[p]
            name = getattr(lab, "label", None) or getattr(lab, "key", None)
        nets.append(_network_of(name))
    return np.asarray(nets)


def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    return s.split("_")[0] if "_" in s else s


def dlabel_network_values(dlabel, N, net_id):
    """Return (tex_left, tex_right): per-vertex network-index maps
    (NaN = unassigned). Built the same way as dlabel_values but with the
    integer network id per parcel, for a categorical topography surface."""
    img = nib.load(dlabel)
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, -1, dtype=int)
    for p in range(1, N + 1):
        nid = int(net_id[p - 1])
        if nid >= 0:
            full[label_data == p] = nid
    ax = img.header.get_axis(1)
    va = np.asarray(getattr(ax, "vertex", None), dtype=int)
    if va is None or len(va) != len(full):
        raise ValueError("dlabel vertex/structure mismatch for network map")
    nv = getattr(ax, "nvertices", {}) or {}
    tex_left = tex_right = None
    pos = 0
    for structure, count in nv.items():
        s = str(structure).upper()
        c = int(count)
        if c == 0:
            continue
        blk = va[pos:pos + c]
        seg = full[pos:pos + c].astype(float)
        seg[seg < 0] = np.nan
        if "LEFT" in s:
            tex_left = seg
        elif "RIGHT" in s:
            tex_right = seg
        pos += c
    return tex_left, tex_right


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def plot_surface(val_left, val_right, surf_l, surf_r, outdir, tag, vmax):
    # val_left/right are per-vertex h2 value arrays (NaN = unassigned/background),
    # one per hemisphere, length == gifti vertex count for that hemisphere.
    vmax = max(float(vmax) if np.isfinite(vmax) else 1e-6, 1e-6)
    sides = [("left", surf_l, val_left), ("right", surf_r, val_right)]
    sides = [(n, s, v) for n, s, v in sides if s is not None and v is not None]
    if not sides:
        return
    png = os.path.join(outdir, f"{tag}_surface.png")
    # nilearn 0.14 dropped list inputs to plot_surf_stat_map, so draw both
    # hemispheres into one figure by giving each its own 3D axes.
    try:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6 * len(sides), 5))
        for i, (name, surf, val) in enumerate(sides, start=1):
            ax = fig.add_subplot(1, len(sides), i, projection="3d")
            niplot.plot_surf_stat_map(surf, val, hemi=name, axes=ax, figure=fig,
                                      cmap="coolwarm", colorbar=True,
                                      threshold=None, vmin=0.0, vmax=vmax,
                                      title=f"{tag} ({name})")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {png}")
    except Exception:
        for name, surf, val in sides:
            sp = os.path.join(outdir, f"{tag}_surface_{name}.png")
            niplot.plot_surf_stat_map(surf, val, hemi=name,
                                      cmap="coolwarm", colorbar=True,
                                      threshold=None, vmin=0.0, vmax=vmax,
                                      title=f"{tag} ({name})", output_file=sp)
            print(f"  wrote {sp}")
    assigned = int(sum(np.count_nonzero(~np.isnan(v)) for _, _, v in sides))
    if assigned:
        allv = np.concatenate([v[~np.isnan(v)] for _, _, v in sides])
        print(f"  [{tag}] surface assigned={assigned} "
              f"min={np.nanmin(allv):.3f} max={np.nanmax(allv):.3f}")
    else:
        print(f"  [{tag}] surface WARNING: no finite values")


def plot_circular(M, outdir, tag, h2_thr=0.35, networks=None):
    """Circular (chord-style) connectome: nodes evenly spaced on a ring, edges
    drawn for h2 > h2_thr (positive only). When `networks` (array of network
    labels, one per node) is given, nodes are grouped and coloured by network.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    N = M.shape[0]
    if networks is not None and len(networks) == N:
        order = np.argsort(networks, kind="stable")  # group by network
        net_ordered = np.asarray(networks)[order]
    else:
        order = np.arange(N)
        net_ordered = None

    pos = np.empty(N, dtype=int)
    pos[order] = np.arange(N)

    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xy = np.column_stack([np.cos(angles), np.sin(angles)])

    fig, ax = plt.subplots(figsize=(9, 9))
    iu = np.triu_indices(N, 1)
    vals = M[iu]
    mask = vals > h2_thr
    n_edges = int(mask.sum())
    for (a, b), v in zip(zip(iu[0][mask], iu[1][mask]), vals[mask]):
        pa, pb = pos[a], pos[b]
        ax.plot([xy[pa, 0], xy[pb, 0]], [xy[pa, 1], xy[pb, 1]],
                color="red", alpha=0.25, linewidth=0.4)

    if net_ordered is not None:
        uniq = list(dict.fromkeys(net_ordered.tolist()))
        cmap = plt.get_cmap("tab20")
        net_color = {net: cmap(i % 20) for i, net in enumerate(uniq)}
        node_colors = [net_color[net_ordered[i]] for i in range(N)]
        ax.scatter(xy[:, 0], xy[:, 1], s=20, c=node_colors, zorder=5,
                   edgecolors="black", linewidths=0.3)
        handles = [plt.Line2D([0], [0], marker="o", linestyle="",
                              color=net_color[net], label=net) for net in uniq]
        ax.legend(handles=handles, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    else:
        ax.scatter(xy[:, 0], xy[:, 1], s=12, c="black", zorder=5)

    ax.set_aspect("equal")
    ax.axis("off")
    title = f"{tag} circular (h2>{h2_thr}, {n_edges} edges)"
    if net_ordered is not None:
        title += f", {len(uniq)} networks"
    ax.set_title(title)
    png = os.path.join(outdir, f"{tag}_circular.png")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")


def plot_surface_networks(net_left, net_right, net_names, surf_l, surf_r,
                          outdir, tag):
    """Categorical brain-surface map of network topography (one colour per
    network), both hemispheres in one figure."""
    import matplotlib.pyplot as plt
    sides = [("left", surf_l, net_left), ("right", surf_r, net_right)]
    sides = [(n, s, v) for n, s, v in sides if s is not None and v is not None]
    if not sides:
        return
    K = len(net_names)
    cmap = plt.get_cmap("tab20")
    png = os.path.join(outdir, f"{tag}_networks_surface.png")
    try:
        fig = plt.figure(figsize=(6 * len(sides), 5))
        for i, (name, surf, val) in enumerate(sides, start=1):
            ax = fig.add_subplot(1, len(sides), i, projection="3d")
            niplot.plot_surf_stat_map(surf, val, hemi=name, axes=ax, figure=fig,
                                      cmap="tab20", colorbar=False, threshold=None,
                                      vmin=0.0, vmax=max(K - 1, 1),
                                      title=f"{tag} networks ({name})")
        handles = [plt.Line2D([0], [0], marker="o", linestyle="",
                              color=cmap(i / max(K - 1, 1)), label=net_names[i])
                   for i in range(K)]
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7,
                   frameon=False)
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {png}")
    except Exception:
        for name, surf, val in sides:
            sp = os.path.join(outdir, f"{tag}_networks_surface_{name}.png")
            niplot.plot_surf_stat_map(surf, val, hemi=name, cmap="tab20",
                                      colorbar=False, threshold=None,
                                      vmin=0.0, vmax=max(K - 1, 1),
                                      title=f"{tag} networks ({name})",
                                      output_file=sp)
            print(f"  wrote {sp}")


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
    ap.add_argument("--networks-csv", help="CSV mapping dlabel value -> "
                    "network name (row/key matches dlabel value)")
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

    # fixed 0..0.5 scale so every surface heatmap is directly comparable
    vmax = 0.5

    # network grouping (from the network CSV, when provided)
    networks = None
    net_names = None
    net_id = None
    if args.dlabel and args.networks_csv:
        networks = dlabel_networks(args.dlabel, N, args.networks_csv)
        net_names = list(dict.fromkeys(networks.tolist()))
        net_id = np.array([net_names.index(s) for s in networks], dtype=int)
        print(f"[networks] {len(net_names)} groups: {net_names}")
    elif args.dlabel and not args.networks_csv:
        print("[networks] no --networks-csv given; skipping network grouping")

    # ---- cortical-surface node-summary (h2) plots ------------------------
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
        # network topography surface (categorical, one colour per network)
        if networks is not None:
            nl, nr = dlabel_network_values(args.dlabel, N, net_id)
            plot_surface_networks(nl, nr, net_names, surf_l, surf_r,
                                  args.outdir, f"{atlas}_networks")

    # ---- circular connectome (h2 > 0.35, grouped by network) -------------
    for method, M in matrices.items():
        plot_circular(M, args.outdir, f"{atlas}_{method}",
                      h2_thr=0.35, networks=networks)

    print("Done. Outputs in", args.outdir)


if __name__ == "__main__":
    main()
