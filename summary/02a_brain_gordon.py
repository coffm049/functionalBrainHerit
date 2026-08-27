#!/usr/bin/env python3
"""Gordon (352 parcels) brain-space visualization — simple, hardcoded.

Produces for twin and AdjHE (90th-percentile node summary):
  - <out>/gordon_{twin,AdjHE}_surface.png  (both hemispheres, 0–0.5)
  - <out>/gordon_networks_surface.png      (categorical network topography)
  - <out>/gordon_{twin,AdjHE}_circular.png (h2 > 0.35, nodes grouped by network)

Run on the HPC where the dlabel / surfaces and .venv are available:
  bash summary/run_brain_viz.sh
or directly:
  /users/4/coffm049/papers/functionalBrainHerit/.venv/bin/python summary/brain_gordon.py
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
DLABEL = Path("/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii")
SURF_L = Path("/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = Path("/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
NETWORKS_CSV = Path("/users/4/coffm049/papers/brainTemplates/shortDicionary.csv")
OUTDIR = ROOT / "results/summary/brain/gordon"
N = 352
NODE_PCT = 90
VMAX = 0.5  # fixed 0–0.5 heatmap scale


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
    nv = ax.nvertices
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


def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    # Gordon label is like "1_L_Default" -> "Default" (per 03b: re.sub r"^\d{1,3}_._")
    s = re.sub(r"^\d{1,3}_[LR]_", "", s)
    # 03b maps full names to short names used in shortDicionary.csv
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
    return _map.get(s, s)


def _cifti_labels(img):
    """Return {label_value: Cifti2Label} dict for a CIFTI dlabel image.

    The label table lives on axis 0's get_element(0)[1] dict (per 03b:34),
    not on img.labeltable. Handles both dict and list storage.
    """
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


def load_gordon_networks(dlabel, csv_path, N):
    """Parcel (1..N) -> network shortname for Gordon.

    shortDicionary.csv is a 14-row network dictionary: columns name,shortname.
    Parcel->network comes from the dlabel's CIFTI label table (e.g.
    '1_L_Default' -> 'Default' -> 'DMN' via 03b mapping), translated via
    name->shortname.
    """
    df = pd.read_csv(csv_path)
    name_to_short = dict(zip(df["name"].astype(str), df["shortname"].astype(str)))
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
        short = name_to_short.get(net, net)
        if short == net:
            for k, v in name_to_short.items():
                if k.lower() == net.lower():
                    short = v
                    break
        nets.append(short)
    return np.asarray(nets)


def dlabel_network_values(dlabel, N, net_id):
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, -1, dtype=int)
    for p in range(1, N + 1):
        full[label_data == p] = int(net_id[p - 1])
    ax = img.header.get_axis(1)
    va = np.asarray(ax.vertex, dtype=int)
    nv = ax.nvertices
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


def _network_palette(net_names):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    base = plt.get_cmap("tab20").colors
    color_of = {net: base[i % len(base)] for i, net in enumerate(net_names)}
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
    fig = plt.figure(figsize=(6 * len(sides), 5))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                  cmap="coolwarm", colorbar=True, threshold=None,
                                  vmin=0.0, vmax=vmax, title=f"{_display_tag(tag)} ({hemi})")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")
    assigned = sum(np.count_nonzero(~np.isnan(v)) for _, _, v in sides)
    allv = np.concatenate([v[~np.isnan(v)] for _, _, v in sides])
    print(f"  [{tag}] surface assigned={assigned} min={np.nanmin(allv):.3f} max={np.nanmax(allv):.3f}")


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
    fig, ax = plt.subplots(figsize=(9, 9))
    iu = np.triu_indices(N, 1)
    vals = M[iu]
    mask = vals > h2_thr
    n_edges = int(mask.sum())
    # line thickness + alpha rescaled: h2 0.25–0.5 -> lw 0.1–1.0, alpha 0.1–1.0
    for (a, b), v in zip(zip(iu[0][mask], iu[1][mask]), vals[mask]):
        pa, pb = pos[a], pos[b]
        # clip h2 to 0.25–0.5 then rescale to 0.1–1.0
        v_clipped = float(np.clip(v, 0.25, 0.5))
        norm = (v_clipped - 0.25) / (0.5 - 0.25)  # 0→1
        lw = 0.1 + norm * (1.0 - 0.1)
        alpha = 0.1 + norm * (1.0 - 0.1)
        ax.plot([xy[pa, 0], xy[pb, 0]], [xy[pa, 1], xy[pb, 1]], color="red", alpha=alpha, linewidth=lw)
    _, color_of = _network_palette(net_names)
    node_colors = [color_of[net_ordered[i]] for i in range(N)]
    ax.scatter(xy[:, 0], xy[:, 1], s=20, c=node_colors, zorder=5, edgecolors="black", linewidths=0.3)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=color_of[net], label=net) for net in net_names]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{_display_tag(tag)} Circular (h2>{h2_thr}, {n_edges} edges, {len(net_names)} networks)")
    png = outdir / f"{tag}_circular.png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")


def plot_surface_networks(net_left, net_right, net_names, surf_l, surf_r, outdir, tag):
    import matplotlib.pyplot as plt
    sides = [("left", surf_l, net_left), ("right", surf_r, net_right)]
    sides = [(n, s, v) for n, s, v in sides if v is not None]
    K = len(net_names)
    cmap, color_of = _network_palette(net_names)
    png = outdir / f"{tag}_networks_surface.png"
    fig = plt.figure(figsize=(6 * len(sides), 5))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                  cmap=cmap, colorbar=False, threshold=None,
                                  vmin=0.0, vmax=max(K - 1, 1), title=f"{_display_tag(tag)} Networks ({hemi})")
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=color_of[net], label=net) for net in net_names]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(WIDE)
sub = df[df["Set"] == "gordon"]

matrices = {}
node_vals = {}
for method, col in [("twin", "Twin_h2"), ("AdjHE", "h2_gordon_AdjHE_RE")]:
    s = sub[["Pheno", col]].dropna()
    M = rebuild_pconn(s["Pheno"].values, s[col].values.astype(float), N)
    matrices[method] = M
    node_vals[method] = node_summary(M, NODE_PCT)
    print(f"[gordon_{method}] edges={int(np.isfinite(M).sum() // 2)} nodes={N}")

networks = load_gordon_networks(DLABEL, NETWORKS_CSV, N)
net_names = list(dict.fromkeys(networks.tolist()))
net_id = np.array([net_names.index(s) for s in networks], dtype=int)
print(f"[networks] {len(net_names)} groups: {net_names}")

dL, dR = dlabel_values(node_vals["AdjHE"], DLABEL, N)
print(f"[dlabel] left assigned={int(np.count_nonzero(~np.isnan(dL)))} right assigned={int(np.count_nonzero(~np.isnan(dR)))} (N={N})")

for method, node_h2 in node_vals.items():
    val_l, val_r = dlabel_values(node_h2, DLABEL, N)
    plot_surface(val_l, val_r, SURF_L, SURF_R, OUTDIR, f"gordon_{method}", VMAX)

nl, nr = dlabel_network_values(DLABEL, N, net_id)
plot_surface_networks(nl, nr, net_names, SURF_L, SURF_R, OUTDIR, "gordon")

for method, M in matrices.items():
    plot_circular(M, OUTDIR, f"gordon_{method}", networks, net_names)

print(f"Done. Outputs in {OUTDIR}")
