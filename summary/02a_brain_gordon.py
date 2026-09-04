#!/usr/bin/env python3
"""Gordon (352 parcels) brain-space visualization — simple, hardcoded.

Produces for twin and AdjHE-RE (90th-percentile node summary, publication-ready, no titles, minimal whitespace):
  - <out>/gordon_{twin,AdjHE-RE}_surface.png  (both hemispheres, 0–0.5, colorbar not overlapping)
  - <out>/gordon_networks_surface.png      (categorical network topography, 14 nets, subcortical NA hidden)
  - <out>/gordon_{twin,AdjHE-RE}_circular.png (h2 > 0.33, nodes grouped by network, publication-ready, 0.33-1.0 rescale)

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
# Portable: try relative ../brainTemplates first (local Windows), fallback to HPC absolute
def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)
DLABEL = _resolve(ROOT.parent / "brainTemplates" / "Gordon.subcortical.32k_fs_LR.dlabel.nii",
                  "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii")
SURF_L = _resolve(ROOT.parent / "brainTemplates" / "Conte69.L.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = _resolve(ROOT.parent / "brainTemplates" / "Conte69.R.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
NETWORKS_CSV = _resolve(ROOT.parent / "brainTemplates" / "shortDicionary.csv",
                        "/users/4/coffm049/papers/brainTemplates/shortDicionary.csv")
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
    tex_left = tex_right = None
    # Use iter_structures (not nvertices) for correct slice boundaries;
    # proba's nvertices (32492) mismatches actual slices (29696/29716),
    # causing right-hemisphere shift if nvertices is used. Gordon matches
    # but we use robust method for both.
    for structure, sl, _ in ax.iter_structures():
        s = str(structure).upper()
        if sl.stop is None:
            sl = slice(sl.start, len(va))
        blk = va[sl]
        seg = full[sl]
        if "CORTEX_LEFT" in s:
            tex_left = np.full(32492, np.nan)
            valid = blk >= 0
            tex_left[blk[valid]] = seg[valid]
        elif "CORTEX_RIGHT" in s:
            tex_right = np.full(32492, np.nan)
            valid = blk >= 0
            tex_right[blk[valid]] = seg[valid]
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
    s = _map.get(s, s)
    s = s.split("_")[0] if "_" in s else s
    return s


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

    Recreates the correct solution from test_A_LUT (user confirmed): explicit
    LUT Gordon333_LUT.txt ROI_ID -> Network_Name (col 0 -> col 2) is ground truth.
    Falls back to gordon_modules.csv (352 rows) and then CIFTI+shortDicionary.
    Subcortical / unknown kept as NA hidden from legend.
    """
    # Try explicit LUT first (as suggested by user, test_A_LUT was correct) — this is the ground truth per your last comment
    lut_candidates = [
        ROOT / "Gordon333_LUT.txt",
        ROOT / "brainTemplates" / "Gordon333_LUT.txt",
        ROOT.parent / "brainTemplates" / "Gordon333_LUT.txt",
        ROOT.parent / "brainTemplates" / "Gordon_333_LUT.txt",
        Path(r"C:\Users\coffm049\brainTemplates\Gordon333_LUT.txt"),
        Path("/users/4/coffm049/papers/brainTemplates/Gordon333_LUT.txt"),
        Path("/users/4/coffm049/papers/functionalBrainHerit/brainTemplates/Gordon333_LUT.txt"),
        Path("/users/4/coffm049/papers/functionalBrainHerit/Gordon333_LUT.txt"),
        Path("Gordon333_LUT.txt"),
        Path("brainTemplates/Gordon333_LUT.txt"),
        Path("brainTemplate/Gordon333_LUT.txt"),
    ]
    for cand in lut_candidates:
        if cand.exists():
            try:
                # LUT: ROI_ID, ROI_Name, Network_Name, R,G,B,A  (sep whitespace) — user code
                lut = pd.read_csv(cand, sep=r"\s+", header=None, engine="python", comment="#")
                # Heuristic: col 0 = ROI_ID, col 2 = Network_Name (as user code)
                # Handle header if present
                if str(lut.iloc[0,0]).lower() == "roi_id":
                    lut = lut.iloc[1:]
                # Try col 2 first, fallback to col 1 if needed
                try:
                    id_to_net = dict(zip(lut.iloc[:,0].astype(int), lut.iloc[:,2].astype(str)))
                    # Validate that Network_Name looks like a network (contains letters, not numbers)
                    sample_net = list(id_to_net.values())[0]
                    if sample_net.replace("_","").isdigit():
                        raise ValueError("col 2 not network")
                except Exception:
                    id_to_net = dict(zip(lut.iloc[:,0].astype(int), lut.iloc[:,1].astype(str)))
                print(f"  Using LUT {cand} with {len(id_to_net)} entries")
                # Map network names via _network_of + shortDicionary style to match net_names casing
                nets = []
                for p in range(1, N+1):
                    raw = id_to_net.get(p, "NA")
                    net = _network_of(raw)
                    # For Gordon, the LUT's Network_Name is already like "Default", "Visual" etc., map to short
                    # Use shortDicionary if available to canonicalize
                    try:
                        df_lut = pd.read_csv(csv_path)
                        name_to_short = dict(zip(df_lut["name"].astype(str), df_lut["shortname"].astype(str)))
                        short = name_to_short.get(net, net)
                        if short == net:
                            for k,v in name_to_short.items():
                                if k.lower() == net.lower():
                                    short = v
                                    break
                        nets.append(short)
                    except Exception:
                        nets.append(net)
                if len(nets) == N:
                    print(f"  LUT {cand} succeeded: {len(set(nets))} nets")
                    return np.asarray(nets)
            except Exception as e:
                print(f"  LUT {cand} failed: {e}")
                pass
    # Fallback: CIFTI label table + shortDicionary.csv (test_C, also correct per your note that A and C both correct)
    # This uses the dlabel's own RGBA, so VIS occipital = light green etc., correct.
    # Keep gordon_modules.csv as last fallback only for environments without dlabel.
    print("  LUT not found, falling back to CIFTI+shortDict (test_C, also correct)")
    gordon_modules_candidates = [
        Path(__file__).resolve().parents[1].parent / "brainTemplates" / "gordon_modules.csv",
        Path(r"C:\Users\coffm049\brainTemplates\gordon_modules.csv"),
        Path("/users/4/coffm049/papers/brainTemplates/gordon_modules.csv"),
        ROOT.parent / "brainTemplates" / "gordon_modules.csv",
        Path("brainTemplate/gordon_modules.csv"),
        Path("brainTemplates/gordon_modules.csv"),
    ]
    # Actually try CIFTI first (more accurate than gordon_modules.csv tab20)
    # So skip gordon_modules here and go directly to CIFTI as primary fallback
    # (gordon_modules kept only if CIFTI also fails)
    # Fallback: CIFTI label table + shortDicionary.csv (previous behavior)
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
            tmp = np.full(32492, np.nan)
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
    """Return network -> color mapping, using dlabel's true colors when available.

    This fixes the mismatched legend/cortex colors reported (VIS peach vs light purple,
    DMN light green vs dark blue, duplicates DMN/SMl both dark green, etc.) which were
    caused by tab20 insertion-order palette. When a dlabel is provided, we use its
    label-table RGBA for each network (the colors that correctly paint the cortex in
    wb_view), giving distinct, anatomy-matched colors and no overlap with the
    nilearn dark-grey background.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    # Try to use dlabel's true network colors
    if dlabel is not None:
        try:
            img = nib.load(str(dlabel))
            labels = _cifti_labels(img)
            # Map network shortname -> rgba from label table (first parcel with that network)
            net_to_rgba = {}
            # For Gordon via gordon_modules.csv, map via labelDictionary to find a representative parcel
            # For fallback, use _network_of on label strings
            for p, lab in labels.items():
                if p == 0:
                    continue
                raw = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
                net = _network_of(raw)
                # Translate via shortDicionary if possible to match net_names casing
                # net_names are shortnames like DMN, Vis, etc.
                # Find matching net_names entry case-insensitively
                for target in net_names:
                    if target.lower() == net.lower():
                        if target not in net_to_rgba:
                            rgba = lab[1] if isinstance(lab, (tuple, list)) and len(lab) > 1 else (0.5, 0.5, 0.5, 1.0)
                            # lab[1] is (r,g,b,a) in 0..1
                            if isinstance(rgba, (list, tuple)) and len(rgba) >= 3:
                                net_to_rgba[target] = tuple(float(x) for x in rgba[:3])
                        break
            if len(net_to_rgba) == len(net_names):
                # All networks found in dlabel — use those colors directly
                color_of = {n: net_to_rgba[n] for n in net_names}
                cmap = mcolors.ListedColormap([color_of[n] for n in net_names])
                return cmap, color_of
        except Exception:
            pass
    # Fallback: tab20 in canonical sorted order (not insertion order) to avoid shift
    # Sort net_names alphabetically for stable, distinct colors; then override Sal/SMl per user
    base = plt.get_cmap("tab20").colors
    # Use sorted order for color assignment to avoid DAN-first insertion-order shift
    sorted_nets = sorted(net_names, key=lambda x: x.lower())
    tmp_color = {net: base[i % len(base)] for i, net in enumerate(sorted_nets)}
    # Ensure Sal brown, SMl purple as requested (insula = Sal brown)
    if "Sal" in tmp_color:
        tmp_color["Sal"] = base[5]  # brown
    if "SMl" in tmp_color:
        tmp_color["SMl"] = base[4]  # purple
    # Reorder to original net_names order for display but with fixed colors
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
    # Publication-ready: minimal whitespace, no title, colorbar not overlapping (heatmap is good, keep its settings; brain needs extra right margin)
    fig = plt.figure(figsize=(5.8 * len(sides), 4.2))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        # HPC nilearn <0.12 does not support cbar_kwargs; handle both
        try:
            cbar_kwargs = {"shrink": 0.6, "aspect": 12, "pad": 0.02} if hemi == "right" else None
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=(hemi == "right"), cbar_kwargs=cbar_kwargs,
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
        except TypeError:
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=(hemi == "right"),
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
    # Leave extra space on right for colorbar (fix overlap without cbar_kwargs)
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


def plot_circular(M, outdir, tag, networks, net_names, h2_thr=0.33):
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
    # line thickness + alpha rescaled: h2 0.33–1.0 -> lw 0.1–1.0, alpha 0.1–1.0 (per user: span 0.33-1.0)
    for (a, b), v in zip(zip(iu[0][mask], iu[1][mask]), vals[mask]):
        pa, pb = pos[a], pos[b]
        v_clipped = float(np.clip(v, 0.33, 1.0))
        norm = (v_clipped - 0.33) / (1.0 - 0.33)
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
    # No title for publication; stats saved to CSV for quarto
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
net_names_all = list(dict.fromkeys(networks.tolist()))
# Hide subcortical/NA from legend — keep only the 14 cortical Gordon networks (as in test_B, 14 nets)
# This fixes the duplicate legend entries (DMN/SMl both dark green, PON/Sal light green) and background overlap.
CORTICAL_14 = {"DMN","VIS","Vis","FP","DAN","VAN","SAL","Sal","CO","SMD","SMd","SML","SMl","AUD","Aud","Tpole","MTL","PMN","PON"}
net_names = [n for n in net_names_all if n in CORTICAL_14 or n.lower() in {c.lower() for c in CORTICAL_14}]
# For display, map non-cortical (subcortical NA etc.) to -1 so they become NaN on surface and grey in circular
net_id = np.array([net_names.index(s) if s in net_names else -1 for s in networks], dtype=int)
print(f"[networks] {len(net_names_all)} groups (incl. subcortical): {net_names_all} -> display {len(net_names)} cortical: {net_names}")

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
