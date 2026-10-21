#!/usr/bin/env python3
"""SA (network surface area, 17 Gordon networks) brain-surface visualization — simple & hardcoded.

Follows 03b-topoh2Ests2Brain.py: SA heritability (twin vs AdjHE, 20 PCs) is
per-network (not per-edge), so values are mapped directly via the Gordon
networks dlabel — no pconn rebuild / node-summary needed.

Produces:
  - <out>/SA_{twin,AdjHE}_surface.png/.pdf  (left+right lateral + medial, 0–0.5)

Run on the HPC, where the dlabel and surfaces are available:
  conda activate MASH
  python summary/02c_brain_sa.py

This is the wo_total SA surface. The w_total counterpart is 02d_brain_sa_total.py.

Notes (from 03a/03b/05-topoViz + HPC check 2026-08-27, updated for 14-network fix):
  - SA Set in mash_twin_wide.csv has 17 phenos network_surfarea1..17
    (grep "network_surfarea" -> 17 distinct), but 4,6,17 are null per
    Gordon (no parcels; 05-topoViz.R filters to 1,2,3,5,7,8,9,10,11,12,13,14,15,16).
    Current code maps only the 14 non-null networks; 4/6/17 are excluded to avoid
    index shift (see 01-07 labeling: case_when pheno>=7 -> -2, >=5 -> -1).
  - Template dlabel is Gordon.networks.32k_fs_LR.dlabel.nii.
  - Surfaces are Conte69 inflated 32k.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import nilearn.plotting as niplot

# Imaging Neuroscience requires Arial/Helvetica in figures; importing this
# module applies it. See summary/figure_font.py.
import figure_font  # noqa: F401

# Writes every figure as PNG (for the report and figures.pptx) and PDF (the
# vector file Imaging Neuroscience wants). See summary/plotstyle.py.
from plotstyle import savefig_publication


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
    # handle SA w_total prefix
    method = "_".join(parts[1:])
    # SA w_total: tag like SA_w_total_twin -> method w_total_twin
    if method.lower().startswith("w_total"):
        suffix = method[7:].lstrip("_")
        base = "SA w_total"
        if suffix.lower() == "twin":
            return f"{base} Twin"
        if suffix.lower() in ("adjhe", "adjhe_re", "adjhe-re"):
            return f"{base} AdjHE-RE"
        return f"{base} {suffix}"
    m = method.lower()
    if m == "twin":
        return f"{atlas_disp} Twin"
    if m in ("adjhe", "adjhe_re", "adjhe-re"):
        return f"{atlas_disp} AdjHE-RE"
    return f"{atlas_disp} {method.replace('_', ' ')}"

# --------------------------------------------------------------------------
# Hardcoded paths / constants — edit here if needed, no fallbacks
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)
DLABEL = _resolve(ROOT.parent / "brainTemplates" / "Gordon.networks.32k_fs_LR.dlabel.nii",
                  "/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii")
SURF_L = _resolve(ROOT.parent / "brainTemplates" / "Conte69.L.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = _resolve(ROOT.parent / "brainTemplates" / "Conte69.R.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
OUTDIR = ROOT / "results/summary/brain/SA"
VMAX = 0.5  # fixed 0–0.5 heatmap scale, same as FC gordon/probaConns

# SA networks: 14 Gordon networks (excludes 4,6,17 which are null per 01-07 labeling)
# 05-topoViz.R networks = 1,2,3,5,7,8,9,10,11,12,13,14,15,16 (14 after filtering)
# 03b-topoh2Ests2Brain.py used range(1,15) for SA; 05-topoViz.R case_when adjusted
# network_surfarea4/6/17 -> null (no parcels in Gordon.networks dlabel).
# Keep mapping consistent: pheno numbers 4,6,17 will be dropped (null) to avoid index shift.
NETWORKS = {
    1: "DMN", 2: "VIS", 3: "FP", 5: "DAN", 7: "VAN",
    8: "SAL", 9: "CO", 10: "SMD", 11: "SML", 12: "AUD",
    13: "Tpole", 14: "MTL", 15: "PMN", 16: "PON",
}


def _cifti_labels(img):
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
    return {}


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
    return s.split("_")[0] if "_" in s else s


def sa_region_from_pheno(pheno: str) -> int | None:
    """'network_surfarea12' -> 12, '12' -> 12, else None."""
    s = str(pheno)
    if "network_surfarea" in s:
        try:
            return int(s.split("network_surfarea")[-1])
        except ValueError:
            return None
    try:
        return int(s)
    except ValueError:
        return None


def sa_values_to_textures(sa_h2: dict[int, float], dlabel: Path):
    """Map per-network SA h2 (network index -> h2) onto per-vertex textures.

    SA h2 is per Gordon network; per 01-07 labeling (05-topoViz.R, 03b range(1,15))
    only 14 networks are non-null: 1,2,3,5,7,8,9,10,11,12,13,14,15,16. The dlabel
    is parcel-level (333 parcels like '1_L_Default'). For each SA network k
    with value v and short name NETWORKS[k], find all parcels p whose network
    (via label table) equals that short name, then set vertices where dlabel == p to v.
    Regions 4,6,17 are null and skipped to avoid index shift.
    """
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    # parcel index -> network short name via label table
    labels = _cifti_labels(img)
    parcel_to_net = {}
    for p, lab in labels.items():
        if p == 0:
            continue
        lab_str = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
        parcel_to_net[p] = _network_of(lab_str)
    # invert: network short name -> list of parcel indices, keyed case-
    # insensitively. _network_of() returns short names from its own map
    # ("Visual"->"Vis", "Salience"->"Sal", "SMhand"->"SMd", "SMmouth"->"SMl",
    # "Auditory"->"Aud", ...), but NETWORKS above spells several of them in a
    # different case ("VIS", "SAL", "SMD", "SML", "AUD"). The exact-match lookup
    # that used to be here therefore painted nothing for those networks. Lower
    # casing both sides makes the join spelling-agnostic.
    net_to_parcels = {}
    for p, net in parcel_to_net.items():
        net_to_parcels.setdefault(str(net).lower(), []).append(p)
    painted = set()
    for region, value in sa_h2.items():
        if not np.isfinite(value):
            continue
        net_name = NETWORKS.get(region)
        if net_name is None:
            # 4,6,17 are null per Gordon (no parcels) — skip to avoid shift
            continue
        parcels = net_to_parcels.get(str(net_name).lower(), [])
        if parcels:
            painted.add(str(net_name))
        for p in parcels:
            full[label_data == p] = value
    missing = sorted({str(v) for k, v in NETWORKS.items() if k in sa_h2} - painted)
    if missing:
        print(f"  WARNING: SA networks with no parcels in this dlabel: {missing}")
    ax = img.header.get_axis(1)
    va = np.asarray(ax.vertex, dtype=int)
    tex_left = tex_right = None
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


def plot_sa_surface(val_left, val_right, surf_l, surf_r, outdir, tag, vmax):
    import matplotlib.pyplot as plt
    # Four views: lateral and medial of each hemisphere. Both hemispheres are
    # computed and both reach the stats CSV.
    sides = [("left", surf_l, val_left, "lateral"), ("left", surf_l, val_left, "medial"),
             ("right", surf_r, val_right, "lateral"), ("right", surf_r, val_right, "medial")]
    sides = [(h, s, v, w) for h, s, v, w in sides if v is not None]
    png = outdir / f"{tag}_surface.png"
    fig = plt.figure(figsize=(5.8 * len(sides), 4.2))
    last = len(sides)
    for i, (hemi, surf, val, view) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        # Colour bar on the last panel only, with a generous pad so it does not
        # land on the cortex.
        colorbar = (i == last)
        try:
            cbar_kwargs = {"shrink": 0.6, "aspect": 12, "pad": 0.1} if colorbar else None
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, view=view, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=colorbar, cbar_kwargs=cbar_kwargs,
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
        except TypeError:
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, view=view, axes=ax, figure=fig,
                                       cmap="coolwarm", colorbar=colorbar,
                                       threshold=None, vmin=0.0, vmax=vmax, title="")
    fig.subplots_adjust(right=0.86, wspace=0.18)
    fig.tight_layout(pad=0.4, rect=[0, 0, 0.86, 1])
    savefig_publication(fig, png)
    stats_path = outdir / f"{tag}_surface_stats.csv"
    try:
        both = [v for v in (val_left, val_right) if v is not None]
        assigned = sum(np.count_nonzero(~np.isnan(v)) for v in both)
        allv = np.concatenate([v[~np.isnan(v)] for v in both]) if assigned>0 else np.array([np.nan])
        pd.DataFrame([{"tag": tag, "assigned_vertices": int(assigned), "h2_min": float(np.nanmin(allv)), "h2_max": float(np.nanmax(allv)), "h2_mean": float(np.nanmean(allv))}]).to_csv(stats_path, index=False)
    except Exception:
        pass
    print(f"  wrote {png}")


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(WIDE)
sub = df[df["Set"] == "SA"]

# Build per-network SA h2 dicts for twin and AdjHE (20 PCs)
sa = {}
for method, col in [("twin", "Twin_h2"), ("AdjHE", "h2_SA_AdjHE_RE")]:
    s = sub[["Pheno", col]].dropna()
    # Pheno like "network_surfarea7" -> region 7; keep only those in NETWORKS or all
    mapping = {}
    for pheno, h2 in zip(s["Pheno"].values, s[col].values.astype(float)):
        region = sa_region_from_pheno(pheno)
        if region is not None:
            mapping[region] = float(h2)
    sa[method] = mapping
    print(f"[SA_{method}] networks={len(mapping)} values={mapping}")

for method, mapping in sa.items():
    # Use the max network index present as N for texture sizing is handled by dlabel
    tex_l, tex_r = sa_values_to_textures(mapping, DLABEL)
    plot_sa_surface(tex_l, tex_r, SURF_L, SURF_R, OUTDIR, f"SA_{method}", VMAX)

print(f"Done. Outputs in {OUTDIR}")
print("To reproduce the 14-network topo figure (05-topoViz.R), filter NETWORKS to 1,2,3,5,7,8,9,10,11,12,13,14,15,16 before mapping.")
