#!/usr/bin/env python3
"""Twin common-effect (C) Manhattan panels and brain maps.

Reads results/summary/mash_twin_wide.csv and plots the twin model's
common/shared-environment component C for Gordon, ProbaConns and SA.

Manhattan panels reuse the Sys-Sys ordering from 05_manhattan.py.
Brain maps show the 90th-percentile node C (Gordon/ProbaConns) or the
per-network C (SA).

Run: python summary/05b_common_effect.py
"""
import importlib
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
manh = importlib.import_module("05_manhattan")

from plotstyle import (ALT_COLORS, FS_AXIS_LABEL, FS_AXIS_LABEL_X, FS_TICK,
                       FS_TICK_ROTATED, MANHATTAN_BOTTOM, MANHATTAN_FIGSIZE,
                       SMALL_COLOR, savefig_publication)

WIDE = ROOT / "results/summary/mash_twin_wide.csv"
PLOT_DIR = ROOT / "results/summary/plots"
BRAIN_DIR = ROOT / "results/summary/brain"
PLOT_DIR.mkdir(parents=True, exist_ok=True)
BRAIN_DIR.mkdir(parents=True, exist_ok=True)

TEMPLATE_DIRS = [
    ROOT.parent / "brainTemplates",
    Path(r"C:\Users\coffm049\brainTemplates"),
    Path("/users/4/coffm049/papers/brainTemplates"),
    Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters"),
]

SA_NETWORKS = manh.SA_NETWORKS


def find_template(name):
    for d in TEMPLATE_DIRS:
        p = d / name
        if p.exists():
            return p
    raise FileNotFoundError(name)


def plot_manhattan(df, atlas, shared_order):
    df, layout, large_order, largest_20 = manh.prepare_manhattan_df(
        df, shared_order=shared_order)

    fig, ax = plt.subplots(1, figsize=MANHATTAN_FIGSIZE)
    x_labels, x_pos = [], []
    for name, group in df.groupby("connection", observed=True):
        if name in largest_20:
            col = ALT_COLORS[large_order.index(name) % 2]
            ax.scatter(group["plot_index"], group["h2"], color=col, s=9,
                       alpha=0.9, zorder=10, linewidths=0.3,
                       edgecolors="black")
            mid = (group["plot_index"].iloc[0] + group["plot_index"].iloc[-1]) / 2.0
            x_labels.append(str(name))
            x_pos.append(mid)
        else:
            ax.scatter(group["plot_index"], group["h2"], color=SMALL_COLOR,
                       s=6, alpha=0.35, zorder=10, linewidths=0)

    if x_pos:
        order = np.argsort(x_pos)
        ax.set_xticks([x_pos[i] for i in order])
        ax.set_xticklabels([x_labels[i] for i in order], rotation=45,
                           ha="right", fontsize=FS_TICK_ROTATED,
                           fontweight="bold")
        ax.tick_params(axis="x", pad=2)
        ax.set_xlabel("Sys-Sys connection", fontsize=FS_AXIS_LABEL_X)

    ax.set_xlim(layout["xlim"])
    ax.set_ylim([0, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=FS_TICK)
    ax.set_ylabel(r"Common effect ($C$)", size=FS_AXIS_LABEL)

    fig.subplots_adjust(bottom=MANHATTAN_BOTTOM)
    out = PLOT_DIR / f"manhattan_{atlas}_TwinCommonEffect.png"
    savefig_publication(fig, out)
    print(f"  wrote {out}")


def rebuild_pconn(phenos, values, N):
    k = np.array([int(str(p)[1:]) if str(p).startswith("o") else int(str(p))
                  for p in phenos])
    M = np.zeros((N, N))
    i, j = np.triu_indices(N, 1)
    M[i[k], j[k]] = values
    M[j[k], i[k]] = values
    np.fill_diagonal(M, np.nan)
    return M


def node_summary(M, pct):
    out = np.full(M.shape[0], np.nan)
    for i in range(M.shape[0]):
        vals = M[i, :][~np.isnan(M[i, :])]
        if len(vals):
            out[i] = np.percentile(vals, pct)
    return out


def dlabel_values(node_vals, dlabel, N):
    import nibabel as nib
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    for p in range(1, N + 1):
        v = node_vals[p - 1]
        if np.isfinite(v):
            full[label_data == p] = v
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


def _cifti_labels(img):
    try:
        ax0 = img.header.get_axis(0)
        d = ax0.get_element(0)[1]
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
    mapping = {
        "Default": "DMN", "Auditory": "AUD", "CinguloOperc": "CO",
        "FrontoParietal": "FP", "Salience": "SAL", "VentralAttn": "VAN",
        "DorsalAttn": "DAN", "Visual": "VIS", "SMhand": "SMD",
        "SMmouth": "SML", "RetrosplenialTemporal": "Tpole",
    }
    s = mapping.get(s, s)
    return s.split("_")[0] if "_" in s else s


def sa_region_from_pheno(pheno):
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


def sa_values_to_textures(mapping, dlabel):
    import nibabel as nib
    img = nib.load(str(dlabel))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    full = np.full(label_data.shape, np.nan, dtype=float)
    labels = _cifti_labels(img)
    net_to_parcels = {}
    for p, lab in labels.items():
        if p == 0:
            continue
        lab_str = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
        net_to_parcels.setdefault(str(_network_of(lab_str)).lower(), []).append(p)
    for region, value in mapping.items():
        if not np.isfinite(value):
            continue
        net_name = SA_NETWORKS.get(region)
        if net_name is None:
            continue
        for p in net_to_parcels.get(str(net_name).lower(), []):
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


def plot_surface(val_left, val_right, surf_l, surf_r, out, tag):
    import nilearn.plotting as niplot
    sides = [("left", surf_l, val_left, "lateral"),
             ("left", surf_l, val_left, "medial"),
             ("right", surf_r, val_right, "lateral"),
             ("right", surf_r, val_right, "medial")]
    sides = [(h, s, v, w) for h, s, v, w in sides if v is not None]
    vmax = 0.5

    fig = plt.figure(figsize=(5.8 * len(sides), 4.2))
    last = len(sides)
    for i, (hemi, surf, val, view) in enumerate(sides, start=1):
        ax = fig.add_subplot(1, len(sides), i, projection="3d")
        colorbar = (i == last)
        try:
            cbar_kwargs = {"shrink": 0.6, "aspect": 12, "pad": 0.1} if colorbar else None
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, view=view,
                                      axes=ax, figure=fig, cmap="coolwarm",
                                      colorbar=colorbar, cbar_kwargs=cbar_kwargs,
                                      threshold=None, vmin=0.0, vmax=vmax, title="")
        except TypeError:
            niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, view=view,
                                      axes=ax, figure=fig, cmap="coolwarm",
                                      colorbar=colorbar, threshold=None,
                                      vmin=0.0, vmax=vmax, title="")
    fig.subplots_adjust(right=0.86, wspace=0.18)
    fig.tight_layout(pad=0.4, rect=[0, 0, 0.86, 1])
    savefig_publication(fig, out)
    print(f"  wrote {out}")


def plot_brain(atlas, N, wide):
    try:
        import nibabel  # noqa: F401
        import nilearn  # noqa: F401
    except Exception as e:
        print(f"  skip brain {atlas}: {e}")
        return
    try:
        if atlas == "gordon":
            dlabel = find_template("Gordon.subcortical.32k_fs_LR.dlabel.nii")
        elif atlas == "probaConns":
            dlabel = find_template("combined_clusters_thresh0.8.dlabel.nii")
        else:
            dlabel = find_template("Gordon.networks.32k_fs_LR.dlabel.nii")
        surf_l = find_template("Conte69.L.inflated.32k_fs_LR.surf.gii")
        surf_r = find_template("Conte69.R.inflated.32k_fs_LR.surf.gii")
    except FileNotFoundError as e:
        print(f"  skip brain {atlas}: missing {e}")
        return

    sub = wide[wide["Set"] == atlas][["Pheno", "C"]].dropna()
    tag = f"{atlas}_TwinCommonEffect"

    if atlas == "SA":
        mapping = {}
        for pheno, val in zip(sub["Pheno"], sub["C"]):
            region = sa_region_from_pheno(pheno)
            if region is not None:
                mapping[region] = float(val)
        val_l, val_r = sa_values_to_textures(mapping, dlabel)
    else:
        M = rebuild_pconn(sub["Pheno"].values, sub["C"].values.astype(float), N)
        node_vals = node_summary(M, 90)
        val_l, val_r = dlabel_values(node_vals, dlabel, N)

    outdir = BRAIN_DIR / atlas
    outdir.mkdir(parents=True, exist_ok=True)
    plot_surface(val_l, val_r, surf_l, surf_r,
                 outdir / f"{tag}_surface.png", tag)


def main():
    wide = manh.load_wide()
    for atlas in ["gordon", "probaConns", "SA"]:
        N = manh.set_N[atlas]
        df = manh.build_long_for_set(wide, atlas, N, "C")
        if df.empty:
            print(f"skip {atlas}: no C data")
            continue
        print(f"{atlas}: {len(df)} edges")
        h2_dfs = {}
        for method, col in manh.specs_for(atlas):
            d = manh.build_long_for_set(wide, atlas, N, col)
            if not d.empty:
                h2_dfs[method] = d
        shared_order = manh._get_shared_order(h2_dfs)
        plot_manhattan(df, atlas, shared_order)
        plot_brain(atlas, N, wide)
    print("Done.")


if __name__ == "__main__":
    main()
