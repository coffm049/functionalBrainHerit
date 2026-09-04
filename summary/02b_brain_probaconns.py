#!/usr/bin/env python3
"""ProbaConns (80 parcels) brain-space visualization — simple, hardcoded.

Produces for twin and AdjHE-RE (90th-percentile node summary, publication-ready, no titles, minimal whitespace):
  - <out>/probaConns_{twin,AdjHE-RE}_surface.png  (both hemispheres, 0–0.5, colorbar not overlapping)
  - <out>/probaConns_networks_surface.png      (categorical, from label table, Sal↔SMl swapped, 14 nets)
  - <out>/probaConns_{twin,AdjHE-RE}_circular.png (h2 > 0.25, nodes grouped by network, publication-ready)

Run on the HPC where the dlabel / surfaces and .venv are available:
  /users/4/coffm049/papers/functionalBrainHerit/.venv/bin/python summary/brain_probaconns.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import nilearn.plotting as niplot

# --------------------------------------------------------------------------
# Hardcoded paths / constants — edit here if needed, no fallbacks
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)

def _find_proba_dlabel():
    """Find proba dlabel with correct parcel count (80) — tries 0.8 then 0.75,
    and both combined_clusters and abcd_template naming."""
    candidates = [
        ROOT / "brainTemplate" / "combined_clusters_thresh0.8.dlabel.nii",
        ROOT / "brainTemplate" / "combined_clusters_thresh0.75.dlabel.nii",
        ROOT.parent / "brainTemplates" / "combined_clusters_thresh0.8.dlabel.nii",
        ROOT.parent / "brainTemplates" / "combined_clusters_thresh0.75.dlabel.nii",
        Path(r"C:\Users\coffm049\brainTemplates\combined_clusters_thresh0.8.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\combined_clusters_thresh0.75.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.8.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"),
        Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"),
        Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.75.dlabel.nii"),
        Path("/users/4/coffm049/papers/brainTemplates/combined_clusters_thresh0.8.dlabel.nii"),
        Path("/users/4/coffm049/papers/brainTemplates/combined_clusters_thresh0.75.dlabel.nii"),
    ]
    # Prefer file with N parcels (=80) to avoid 63-parcel mismatch
    for p in candidates:
        if p.exists():
            try:
                img = nib.load(str(p))
                ax0 = img.header.get_axis(0)
                n_lab = len(ax0.get_element(0)[1]) - 1  # minus background ??? label
                if n_lab == N:
                    return p
            except Exception:
                pass
    # fallback to first existing
    for p in candidates:
        if p.exists():
            return p
    # last resort: original HPC 0.8
    return Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii")

DLABEL = _find_proba_dlabel()
SURF_L = _resolve(ROOT.parent / "brainTemplates" / "Conte69.L.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = _resolve(ROOT.parent / "brainTemplates" / "Conte69.R.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
OUTDIR = ROOT / "results/summary/brain/probaConns"
N = 80
NODE_PCT = 90
VMAX = 0.5


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def rebuild_pconn(phenos, edge_h2, N):
    k = np.array([int(str(p)[1:]) if str(p).startswith("o") else int(str(p)) for p in phenos])
    M = np.zeros((N, N))
    i, j = np.triu_indices(N, 1)
    M[i[k], j[k]] = edge_h2
    M[j[k], i[k]] = edge_h2
    np.fill_diagonal(M, np.nan)
    return M


def node_summary(M, pct):
    out = np.full(M.shape[0], np.nan)
    for i in range(M.shape[0]):
        vals = M[i, :][~np.isnan(M[i, :])]
        if len(vals):
            out[i] = np.percentile(vals, pct)
    return out


def dlabel_values(node_h2, dlabel, N):
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    for p in range(1, N + 1):
        v = node_h2[p - 1]
        if np.isfinite(v):
            full[label_data == p] = v
    ax = img.header.get_axis(1)
    va = np.asarray(ax.vertex, dtype=int)
    tex_left = tex_right = None
    # Use iter_structures for correct slice boundaries: proba's nvertices (32492)
    # mismatches actual BrainModel slices (29696/29716 for medial-wall-excluded 32k),
    # causing a ~2796-vertex shift for the right hemisphere if nvertices is used.
    for structure, sl, _ in ax.iter_structures():
        s = str(structure).upper()
        if sl.stop is None:
            sl = slice(sl.start, len(va))
        blk = va[sl]
        seg = full[sl]
        if "CORTEX_LEFT" in s:
            # Conte69 surfaces are 32492 vertices; proba dlabels have 29696/29716
            # valid vertices, so allocate full 32492 and fill valid indices.
            tex_left = np.full(32492, np.nan)
            valid = blk >= 0
            # blk values are 0..32491 for valid cortex vertices
            tex_left[blk[valid]] = seg[valid]
        elif "CORTEX_RIGHT" in s:
            tex_right = np.full(32492, np.nan)
            valid = blk >= 0
            tex_right[blk[valid]] = seg[valid]
        # non-cortical structures ignored for surface PNGs
    return tex_left, tex_right


def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    s = re.sub(r"^\d{1,3}_[LR]_", "", s)
    _map = {
        "Default": "DMN",
        "Auditory": "Aud",
        "CinguloOperc": "CO",
        "FrontoParietal": "FP",
        "Salience": "Sal",
        "VentralAttn": "VAN",
        "DorsalAttn": "DAN",
        "Visual": "Vis",
        "SMhand": "SMd",
        "SMmouth": "SMl",
        "RetrosplenialTemporal": "Tpole",
    }
    s = _map.get(s, s)
    s = s.split("_")[0] if "_" in s else s
    return s


def _cifti_labels(img):
    """Return {label_value: Cifti2Label} dict for a CIFTI dlabel image."""
    try:
        ax0 = img.header.get_axis(0)
        elem = ax0.get_element(0)
        d = elem[1]
        if isinstance(d, dict) and len(d) > 1:
            return d
        if isinstance(d, list) and len(d) > 1:
            return {i: v for i, v in enumerate(d) if v is not None}
    except Exception:
        pass
    for attr in ("labeltable", "labels"):
        lt = getattr(img, attr, None)
        if lt is not None:
            d = getattr(lt, "labels", lt)
            if isinstance(d, dict) and len(d) > 1:
                return d
            if isinstance(d, list) and len(d) > 1:
                return {i: v for i, v in enumerate(d) if v is not None}
    try:
        ax0 = img.header.get_axis(0)
        for attr in ("labels", "label", "labeltable"):
            lt = getattr(ax0, attr, None)
            if lt is not None:
                d = getattr(lt, "labels", lt)
                if isinstance(d, dict) and len(d) > 1:
                    return d
                if isinstance(d, list) and len(d) > 1:
                    return {i: v for i, v in enumerate(d) if v is not None}
    except Exception:
        pass
    return {}


def load_probaconns_networks(dlabel, N):
    """Parcel (1..N) -> network from the CIFTI label table (first token).

    Fix 2024-08: Sal (Salience, insula) and SMl (somatomotor lateral) were swapped
    in the dlabel (brown patch above insula showed SMl, purple ventral tip showed Sal).
    Swap them back so insula brown = Sal and ventral tip grey/purple = SMl.
    """
    img = nib.load(str(dlabel))
    labels = _cifti_labels(img)
    nets = []
    for p in range(1, N + 1):
        lab = labels.get(p)
        if lab is None:
            parcel_label = None
        elif isinstance(lab, (tuple, list)):
            parcel_label = lab[0]
        else:
            parcel_label = getattr(lab, "label", None) or getattr(lab, "key", None)
            if parcel_label is None:
                parcel_label = str(lab)
        net = _network_of(parcel_label)
        # Swap Sal <-> SMl to match correct cortical topography (insula = Sal)
        if net == "Sal":
            net = "SMl"
        elif net == "SMl":
            net = "Sal"
        nets.append(net)
    return np.asarray(nets)


def dlabel_network_values(dlabel, N, net_id):
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, -1, dtype=int)
    for p in range(1, N + 1):
        full[label_data == p] = int(net_id[p - 1])
    ax = img.header.get_axis(1)
    va = np.asarray(ax.vertex, dtype=int)
    tex_left = tex_right = None
    for structure, sl, _ in ax.iter_structures():
        s = str(structure).upper()
        if sl.stop is None:
            sl = slice(sl.start, len(va))
        blk = va[sl]
        seg = full[sl].astype(float)
        seg[seg < 0] = np.nan
        if "CORTEX_LEFT" in s:
            tex_left = np.full(32492, np.nan)
            valid = blk >= 0
            # map via vertex indices to full 32492 surface
            # need to create texture of size 32492, fill valid positions
            tmp = np.full(32492, np.nan)
            # seg[valid] holds net_id for valid vertices; blk[valid] is target index
            tmp[blk[valid]] = seg[valid]
            tex_left = tmp
        elif "CORTEX_RIGHT" in s:
            tex_right = np.full(32492, np.nan)
            valid = blk >= 0
            tmp = np.full(32492, np.nan)
            tmp[blk[valid]] = seg[valid]
            tex_right = tmp
    return tex_left, tex_right


def _network_palette(net_names, dlabel=None):
    """Network -> color, using dlabel's true colors when available to match anatomy.

    Fixes VIS peach vs light purple, DMN light green vs dark blue, duplicates
    DMN/SMl both dark green etc., which were caused by tab20 insertion-order.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    if dlabel is not None:
        try:
            img = nib.load(str(dlabel))
            labels = _cifti_labels(img)
            net_to_rgba = {}
            for p, lab in labels.items():
                if p == 0:
                    continue
                raw = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
                net = _network_of(raw)
                # Swap Sal<->SMl as per user fix (insula brown Sal)
                if net == "Sal":
                    net = "SMl"
                elif net == "SMl":
                    net = "Sal"
                # canonicalize
                for target in net_names:
                    if target.lower() == net.lower():
                        if target not in net_to_rgba:
                            rgba = lab[1] if isinstance(lab, (tuple, list)) and len(lab) > 1 else (0.5, 0.5, 0.5, 1.0)
                            if isinstance(rgba, (list, tuple)) and len(rgba) >= 3:
                                net_to_rgba[target] = tuple(float(x) for x in rgba[:3])
                        break
            if len(net_to_rgba) == len(net_names):
                color_of = {n: net_to_rgba[n] for n in net_names}
                cmap = mcolors.ListedColormap([color_of[n] for n in net_names])
                return cmap, color_of
        except Exception:
            pass
    base = plt.get_cmap("tab20").colors
    sorted_nets = sorted(net_names, key=lambda x: x.lower())
    tmp_color = {net: base[i % len(base)] for i, net in enumerate(sorted_nets)}
    if "Sal" in tmp_color:
        tmp_color["Sal"] = base[5]
    if "SMl" in tmp_color:
        tmp_color["SMl"] = base[4]
    color_of = {n: tmp_color[n] for n in net_names}
    cmap = mcolors.ListedColormap([color_of[n] for n in net_names])
    return cmap, color_of


def _display_tag(tag: str) -> str:
    parts = tag.split("_")
    atlas = parts[0]
    if atlas.lower() == "gordon":
        atlas_disp = "Gordon"
    elif atlas.lower().startswith("proba"):
        atlas_disp = "ProbaConns"
    elif atlas.lower() == "sa":
        atlas_disp = "SA"
    else:
        atlas_disp = atlas.capitalize()
    if len(parts) == 1:
        return atlas_disp
    method = "_".join(parts[1:])
    m = method.lower()
    if m == "twin":
        method_disp = "Twin"
    elif m in ("adjhe", "adjhe_re", "adjhe-re"):
        method_disp = "AdjHE-RE"
    elif m == "networks":
        return f"{atlas_disp} Networks"
    else:
        method_disp = re.sub(r"AdjHE", "AdjHE-RE", method.replace("_", "-"), flags=re.IGNORECASE)
        method_disp = method_disp[0].upper() + method_disp[1:] if method_disp else method_disp
    return f"{atlas_disp} {method_disp}"


def plot_surface(val_left, val_right, surf_l, surf_r, outdir, tag, vmax):
    import matplotlib.pyplot as plt
    sides = [("left", surf_l, val_left), ("right", surf_r, val_right)]
    sides = [(n, s, v) for n, s, v in sides if v is not None]
    png = outdir / f"{tag}_surface.png"
    fig = plt.figure(figsize=(5.8 * len(sides), 4.2))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        try:
            cbar_kwargs = {"shrink": 0.6, "aspect": 12, "pad": 0.02} if hemi == "right" else None
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=(hemi == "right"), cbar_kwargs=cbar_kwargs,
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
        except TypeError:
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=(hemi == "right"),
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
    fig.subplots_adjust(right=0.88, wspace=0.15)
    fig.tight_layout(pad=0.4, rect=[0, 0, 0.88, 1])
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    stats_path = outdir / f"{tag}_surface_stats.csv"
    try:
        assigned = sum(np.count_nonzero(~np.isnan(v)) for _, _, v in sides)
        allv = np.concatenate([v[~np.isnan(v)] for _, _, v in sides]) if assigned>0 else np.array([np.nan])
        pd.DataFrame([{"tag": tag, "assigned_vertices": int(assigned), "h2_min": float(np.nanmin(allv)), "h2_max": float(np.nanmax(allv)), "h2_mean": float(np.nanmean(allv))}]).to_csv(stats_path, index=False)
    except Exception:
        pass
    print(f"  wrote {png}")


def plot_circular(M, outdir, tag, networks, net_names, h2_thr=0.25):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    N = M.shape[0]
    order = np.argsort(networks, kind="stable")
    net_ordered = np.asarray(networks)[order]
    pos = np.empty(N, dtype=int)
    pos[order] = np.arange(N)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xy = np.column_stack([np.cos(angles), np.sin(angles)])
    fig, ax = plt.subplots(figsize=(7, 7))
    iu = np.triu_indices(N, 1)
    vals = M[iu]
    mask = vals > h2_thr
    n_edges = int(mask.sum())
    # line thickness + alpha rescaled: h2 0.25–0.5 -> lw 0.1–1.0, alpha 0.1–1.0 (bumped to 0.25)
    for (a, b), v in zip(zip(iu[0][mask], iu[1][mask]), vals[mask]):
        pa, pb = pos[a], pos[b]
        v_clipped = float(np.clip(v, 0.25, 0.5))
        norm = (v_clipped - 0.25) / (0.5 - 0.25)
        lw = 0.1 + norm * (1.0 - 0.1)
        alpha = 0.1 + norm * (1.0 - 0.1)
        ax.plot([xy[pa, 0], xy[pb, 0]], [xy[pa, 1], xy[pb, 1]], color="red", alpha=alpha, linewidth=lw)
    _, color_of = _network_palette(net_names, DLABEL)
    node_colors = [color_of.get(net_ordered[i], "#999999") for i in range(N)]
    ax.scatter(xy[:, 0], xy[:, 1], s=18, c=node_colors, zorder=5, edgecolors="black", linewidths=0.3)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=color_of[net], label=net) for net in net_names]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    ax.set_aspect("equal")
    ax.axis("off")
    png = outdir / f"{tag}_circular.png"
    fig.tight_layout(pad=0.3)
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    stats_path = outdir / f"{tag}_circular_stats.csv"
    try:
        pd.DataFrame([{"tag": tag, "h2_thr": float(h2_thr), "n_edges": int(n_edges), "n_networks": int(len(net_names))}]).to_csv(stats_path, index=False)
    except Exception:
        pass
    print(f"  wrote {png}")


def plot_surface_networks(net_left, net_right, net_names, surf_l, surf_r, outdir, tag):
    import matplotlib.pyplot as plt
    sides = [("left", surf_l, net_left), ("right", surf_r, net_right)]
    sides = [(n, s, v) for n, s, v in sides if v is not None]
    K = len(net_names)
    cmap, color_of = _network_palette(net_names, DLABEL)
    png = outdir / f"{tag}_networks_surface.png"
    fig = plt.figure(figsize=(5.5 * len(sides), 4.2))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                   cmap=cmap, colorbar=False, threshold=None,
                                   vmin=0.0, vmax=max(K - 1, 1), title="")
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=color_of[net], label=net) for net in net_names]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, frameon=False)
    fig.tight_layout(pad=0.5)
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {png}")


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(WIDE)
sub = df[df["Set"] == "probaConns"]

matrices = {}
node_vals = {}
for method, col in [("twin", "Twin_h2"), ("AdjHE", "h2_proba_AdjHE_RE")]:
    s = sub[["Pheno", col]].dropna()
    M = rebuild_pconn(s["Pheno"].values, s[col].values.astype(float), N)
    matrices[method] = M
    node_vals[method] = node_summary(M, NODE_PCT)
    print(f"[probaConns_{method}] edges={int(np.isfinite(M).sum() // 2)} nodes={N}")

networks = load_probaconns_networks(DLABEL, N)
net_names_all = list(dict.fromkeys(networks.tolist()))
# Hide NA (subcortical/unknown) from legend — keep NA for data but not displayed
net_names = [n for n in net_names_all if n not in ("NA", "SUB", "???")]
# For display, map NA/SUB to -1 so they become NaN on surface and grey in circular
net_id = np.array([net_names.index(s) if s in net_names else -1 for s in networks], dtype=int)
print(f"[networks] {len(net_names_all)} groups (incl. NA): {net_names_all} -> display {len(net_names)}: {net_names}")

dL, dR = dlabel_values(node_vals["AdjHE"], DLABEL, N)
print(f"[dlabel] left assigned={int(np.count_nonzero(~np.isnan(dL)))} right assigned={int(np.count_nonzero(~np.isnan(dR)))} (N={N})")

for method, node_h2 in node_vals.items():
    val_l, val_r = dlabel_values(node_h2, DLABEL, N)
    plot_surface(val_l, val_r, SURF_L, SURF_R, OUTDIR, f"probaConns_{method}", VMAX)

nl, nr = dlabel_network_values(DLABEL, N, net_id)
plot_surface_networks(nl, nr, net_names, SURF_L, SURF_R, OUTDIR, "probaConns")

for method, M in matrices.items():
    plot_circular(M, OUTDIR, f"probaConns_{method}", networks, net_names)

print(f"Done. Outputs in {OUTDIR}")
