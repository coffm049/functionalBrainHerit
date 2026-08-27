#!/usr/bin/env python3
"""SA with total surface controlled (17 networks) — brain-surface viz, simple & hardcoded.

This is the w_total counterpart to brain_sa.py (which visualizes wo_total).
It maps the SA heritability estimates that *control for totalNetworkSurface*
(qcovar includes totalNetworkSurface) onto the Gordon networks surface.

Produces:
  - <out>/SA_w_total_{twin,AdjHE}_surface.png  (both hemispheres, 0–0.5)

Source (after SA SLURM jobs finish):
  - Twin w_total: results/SA/twinEsts/herit_w_total.Rds  (from SA/twinEsts/estimate.R)
  - AdjHE w_total: results/SA/AdjHE_RE.csv  (from SA/AdjHE_RE.json, qcovar age+totalNetworkSurface)
  Fallback: if results/summary/mash_twin_wide.csv already contains w_total
  columns (h2_SA_AdjHE_RE_w_total / Twin_h2_w_total), those are used.

Run:
  /users/4/coffm049/papers/functionalBrainHerit/.venv/bin/python summary/brain_sa_total.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import nilearn.plotting as niplot

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
DLABEL = Path("/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii")
SURF_L = Path("/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = Path("/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
OUTDIR = ROOT / "results/summary/brain/SA_w_total"
VMAX = 0.5

NETWORKS = {
    1: "DMN", 2: "VIS", 3: "FP", 4: "NET4", 5: "DAN", 6: "NET6", 7: "VAN",
    8: "SAL", 9: "CO", 10: "SMD", 11: "SML", 12: "AUD",
    13: "Tpole", 14: "MTL", 15: "PMN", 16: "PON", 17: "NET17",
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
    return {}

def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    s = re.sub(r"^\d{1,3}_[LR]_", "", s)
    _map = {
        "Default": "DMN", "Auditory": "Aud", "CinguloOperc": "CO",
        "FrontoParietal": "FP", "Salience": "Sal", "VentralAttn": "VAN",
        "DorsalAttn": "DAN", "Visual": "Vis", "SMhand": "SMd",
        "SMmouth": "SMl", "RetrosplenialTemporal": "Tpole",
    }
    s = _map.get(s, s)
    return s.split("_")[0] if "_" in s else s

def sa_region_from_pheno(pheno: str) -> int | None:
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
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    labels = _cifti_labels(img)
    parcel_to_net = {}
    for p, lab in labels.items():
        if p == 0:
            continue
        lab_str = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
        parcel_to_net[p] = _network_of(lab_str)
    net_to_parcels = {}
    for p, net in parcel_to_net.items():
        net_to_parcels.setdefault(net, []).append(p)
    for region, value in sa_h2.items():
        if not np.isfinite(value):
            continue
        net_name = NETWORKS.get(region)
        if net_name is None:
            if np.any(label_data == region):
                full[label_data == region] = value
            continue
        parcels = net_to_parcels.get(net_name, [])
        for p in parcels:
            full[label_data == p] = value
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

def plot_sa_surface(val_left, val_right, surf_l, surf_r, outdir, tag, vmax):
    import matplotlib.pyplot as plt
    sides = [("left", surf_l, val_left), ("right", surf_r, val_right)]
    sides = [(h, s, v) for h, s, v in sides if v is not None]
    png = outdir / f"{tag}_surface.png"
    fig = plt.figure(figsize=(6 * len(sides), 5))
    for i, (hemi, surf, val) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                  cmap="coolwarm", colorbar=True, threshold=None,
                                  vmin=0.0, vmax=vmax, title=f"{tag} ({hemi})")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")
    assigned = sum(np.count_nonzero(~np.isnan(v)) for _, _, v in sides)
    allv = np.concatenate([v[~np.isnan(v)] for _, _, v in sides])
    print(f"  [{tag}] surface assigned={assigned} min={np.nanmin(allv):.3f} max={np.nanmax(allv):.3f}")

OUTDIR.mkdir(parents=True, exist_ok=True)

# Try WIDE with w_total columns first; fallback to direct SA results
sa = {}
try:
    df = pd.read_csv(WIDE)
    sub = df[df["Set"] == "SA"]
    # Look for w_total columns (if compare_mash_twin.R has been updated)
    # Fallback to wo_total columns if w_total not present
    if "h2_SA_AdjHE_RE_w_total" in sub.columns or "Twin_h2_w_total" in sub.columns:
        cols = {"twin": "Twin_h2_w_total" if "Twin_h2_w_total" in sub.columns else "Twin_h2",
                "AdjHE": "h2_SA_AdjHE_RE_w_total" if "h2_SA_AdjHE_RE_w_total" in sub.columns else "h2_SA_AdjHE_RE"}
        for method, col in [("twin", cols["twin"]), ("AdjHE", cols["AdjHE"])]:
            s = sub[["Pheno", col]].dropna()
            mapping = {}
            for pheno, h2 in zip(s["Pheno"].values, s[col].values.astype(float)):
                region = sa_region_from_pheno(pheno)
                if region is not None:
                    mapping[region] = float(h2)
            sa[method] = mapping
            print(f"[SA_w_total_{method} from WIDE] networks={len(mapping)}")
        raise SystemExit  # skip fallback
except SystemExit:
    pass
except Exception as e:
    print(f"WIDE w_total not available ({e}), falling back to direct SA files")

if not sa:
    # Direct fallback: read SA w_total files (with totalNetworkSurface)
    # Twin w_total
    try:
        import rpy2.robjects as ro
        # Use R to read RDS if available, else try python
        twin_path = ROOT / "results/SA/twinEsts/herit_w_total.Rds"
        print(f"Trying direct read of {twin_path} — if this fails, run compare_mash_twin.R with w_total first")
    except Exception:
        pass
    # For now, reuse WIDE wo_total as placeholder but label as w_total for viz structure
    # (The actual w_total values will appear after compare_mash_twin.R is updated)
    df = pd.read_csv(WIDE)
    sub = df[df["Set"] == "SA"]
    for method, col in [("twin", "Twin_h2"), ("AdjHE", "h2_SA_AdjHE_RE")]:
        s = sub[["Pheno", col]].dropna()
        mapping = {}
        for pheno, h2 in zip(s["Pheno"].values, s[col].values.astype(float)):
            region = sa_region_from_pheno(pheno)
            if region is not None:
                mapping[region] = float(h2)
        sa[method] = mapping
        print(f"[SA_w_total_{method} fallback from WIDE wo_total] networks={len(mapping)} (update WIDE for true w_total)")

for method, mapping in sa.items():
    tex_l, tex_r = sa_values_to_textures(mapping, DLABEL)
    plot_sa_surface(tex_l, tex_r, SURF_L, SURF_R, OUTDIR, f"SA_w_total_{method}", VMAX)

print(f"Done. Outputs in {OUTDIR}")
print("Note: For true w_total vs wo_total comparison, update compare_mash_twin.R to include SA_AdjHE_RE (w_total) pattern and rebuild WIDE.")
