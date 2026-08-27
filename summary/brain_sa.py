#!/usr/bin/env python3
"""SA (network surface area, 14 Gordon networks) brain-surface visualization — draft, simple & hardcoded.

Follows 03b-topoh2Ests2Brain.py: SA heritability (twin vs AdjHE, 30 PCs) is
per-network (not per-edge), so values are mapped directly via the Gordon
networks dlabel (14 networks) — no pconn rebuild / node-summary needed.

Produces:
  - <out>/SA_{twin,AdjHE}_surface.png  (both hemispheres, 0–0.5)

Run on the HPC where the dlabel / surfaces and .venv are available:
  /users/4/coffm049/papers/functionalBrainHerit/.venv/bin/python summary/brain_sa.py

Draft notes (from 03a/03b/05-topoViz):
  - SA Set in mash_twin_wide.csv uses Pheno like "network_surfarea3" (the
    full SA results have 17 phenos 1..17; the topo set used in 05-topoViz.R
    keeps 14 after filtering 1,2,3,5,7,8,9,10,11,12,13,14,15,16). This draft
    keeps the 14 present in mash_twin_wide.csv for Set=="SA" and maps them
    via the integer suffix; if your SA wide has 17, adjust N / the filter.
  - Template dlabel is Gordon.networks.32k_fs_LR.dlabel.nii (14 networks).
  - Surfaces are Conte69 inflated 32k.
"""
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
DLABEL = Path("/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii")
SURF_L = Path("/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = Path("/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")
OUTDIR = ROOT / "results/summary/brain/SA"
VMAX = 0.5  # fixed 0–0.5 heatmap scale, same as FC gordon/probaConns

# SA networks as in 05-topoViz.R (index -> short name, 1-based)
NETWORKS = {
    1: "DMN", 2: "VIS", 3: "FP", 5: "DAN", 7: "VAN",
    8: "SAL", 9: "CO", 10: "SMD", 11: "SML", 12: "AUD",
    13: "Tpole", 14: "MTL", 15: "PMN", 16: "PON",
}
# The SA dlabel used here has 14 distinct network values (1..14 after the
# filtering in 03b). Keep N as the max network index present.


def sa_region_from_pheno(pheno: str) -> int | None:
    """'network_surfarea12' -> 12, '12' -> 12, else None."""
    s = str(pheno)
    # 03b used pheno.str[16:] to strip "network_surfarea"
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
    """Map per-network SA h2 (network index -> h2) onto per-vertex textures
    via 01-07 style: vertex_map[dlabel == region] = value.
    """
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    for region, value in sa_h2.items():
        if np.isfinite(value):
            full[label_data == region] = value
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
print("Note: if your SA wide has 17 phenos 1..17, extend NETWORKS and the dlabel mapping accordingly;"
      " see 03b-topoh2Ests2Brain.py and 05-topoViz.R for the 14-network filtering.")
