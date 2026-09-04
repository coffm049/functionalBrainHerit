#!/usr/bin/env python3
"""SA (network surface area, 17 Gordon networks) brain-surface visualization — simple & hardcoded.

Follows 03b-topoh2Ests2Brain.py: SA heritability (twin vs AdjHE, 30 PCs) is
per-network (not per-edge), so values are mapped directly via the Gordon
networks dlabel — no pconn rebuild / node-summary needed.

Produces:
  - <out>/SA_{twin,AdjHE}_surface.png  (both hemispheres, 0–0.5)

Run on the HPC where the dlabel / surfaces and .venv are available:
  /users/4/coffm049/papers/functionalBrainHerit/.venv/bin/python summary/brain_sa.py

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
    # invert: network short name -> list of parcel indices
    net_to_parcels = {}
    for p, net in parcel_to_net.items():
        net_to_parcels.setdefault(net, []).append(p)
    for region, value in sa_h2.items():
        if not np.isfinite(value):
            continue
        net_name = NETWORKS.get(region)
        if net_name is None:
            # 4,6,17 are null per Gordon (no parcels) — skip to avoid shift
            continue
        parcels = net_to_parcels.get(net_name, [])
        for p in parcels:
            full[label_data == p] = value
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
    sides = [("left", surf_l, val_left), ("right", surf_r, val_right)]
    sides = [(h, s, v) for h, s, v in sides if v is not None]
    png = outdir / f"{tag}_surface.png"
    fig = plt.figure(figsize=(5.8 * len(sides), 4.2))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        cbar_kwargs = {"shrink": 0.6, "aspect": 12, "pad": 0.02} if hemi == "right" else None
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                   cmap="coolwarm", colorbar=(hemi == "right"), cbar_kwargs=cbar_kwargs,
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


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(WIDE)
sub = df[df["Set"] == "SA"]

# Build per-network SA h2 dicts for twin and AdjHE (30 PCs)
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
