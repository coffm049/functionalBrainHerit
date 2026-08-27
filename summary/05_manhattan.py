#!/usr/bin/env python3
"""Pairwise-system Manhattan — grouped by Sys-Sys, sorted by median h2.

For each atlas (Gordon 352, ProbaConns 80, SA 17) and both methods
(Twin, AdjHE-RE) from mash_twin_wide.csv (30 PCs):

  x = edge index grouped by Sys-Sys network pair (e.g. DMN-VIS), ordered by
      largest median h2 within that Sys-Sys group (descending).
  y = h2 (0–1, full heritability range).
  Largest 20 Sys-Sys groups by number of connections are plotted normally
  with alternating colours; remaining groups are compressed 100× on the
  x-axis and shown in grey.

Outputs: results/summary/plots/manhattan_{Set}_{Method}.png
Run: .venv/bin/python summary/05_manhattan.py
"""
import math
import re
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
PLOT_DIR = ROOT / "results/summary/plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

GORDON_DLABEL = Path(
    "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii"
)
PROBA_DLABEL = Path(
    "/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/"
    "ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"
)
SHORT_DICT = Path("/users/4/coffm049/papers/brainTemplates/shortDicionary.csv")
SA_NETWORKS = {
    1: "DMN", 2: "VIS", 3: "FP", 4: "NET4", 5: "DAN", 6: "NET6",
    7: "VAN", 8: "SAL", 9: "CO", 10: "SMD", 11: "SML", 12: "AUD",
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
        "Default": "DMN", "Auditory": "Aud",
        "CinguloOperc": "CO", "FrontoParietal": "FP",
        "Salience": "Sal", "VentralAttn": "VAN",
        "DorsalAttn": "DAN", "Visual": "Vis",
        "SMhand": "SMd", "SMmouth": "SMl",
        "RetrosplenialTemporal": "Tpole",
    }
    s = _map.get(s, s)
    return s.split("_")[0] if "_" in s else s


def load_networks_for_atlas(atlas, N):
    import nibabel as nib
    if atlas == "gordon":
        df = pd.read_csv(SHORT_DICT)
        name_to_short = dict(zip(df["name"].astype(str), df["shortname"].astype(str)))
        img = nib.load(str(GORDON_DLABEL))
        labels = _cifti_labels(img)
        nets = []
        for p in range(1, N + 1):
            lab = labels.get(p)
            s = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab)) if lab is not None else None
            net = _network_of(s)
            nets.append(name_to_short.get(net, net))
            if nets[-1] == net:
                for k, v in name_to_short.items():
                    if k.lower() == net.lower():
                        nets[-1] = v
                        break
        return np.asarray(nets)
    elif atlas == "probaConns":
        import nibabel as nib
        img = nib.load(str(PROBA_DLABEL))
        labels = _cifti_labels(img)
        nets = []
        for p in range(1, N + 1):
            lab = labels.get(p)
            s = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab)) if lab is not None else None
            nets.append(_network_of(s))
        return np.asarray(nets)
    else:
        return np.array([SA_NETWORKS.get(i + 1, f"NET{i+1}") for i in range(N)])


def decode_pheno_to_ij(pheno, N):
    s = str(pheno)
    k = int(s[1:]) if s.startswith("o") else int(s)
    i, j = np.triu_indices(N, 1)
    return int(i[k]), int(j[k])


def build_long_for_set(wide, atlas, N, h2_col):
    sub = wide[wide["Set"] == atlas][["Pheno", h2_col]].dropna()
    sub = sub.copy()
    if atlas == "SA":
        sub["connection"] = sub["Pheno"].apply(
            lambda p: SA_NETWORKS.get(int(str(p).split("network_surfarea")[-1]) if "network_surfarea" in str(p) else 0, str(p))
        )
        sub["h2"] = sub[h2_col].astype(float)
        sub["signif"] = sub["h2"] > 0.2
        return sub[["connection", "h2", "signif"]]
    networks = load_networks_for_atlas(atlas, N)
    rows = []
    for _, r in sub.iterrows():
        pheno, h2 = r["Pheno"], float(r[h2_col])
        try:
            i, j = decode_pheno_to_ij(pheno, N)
            ni = networks[i] if i < len(networks) else "NA"
            nj = networks[j] if j < len(networks) else "NA"
            conn = f"{ni}-{nj}"
        except Exception:
            conn = str(pheno)
        rows.append((conn, h2, h2 > 0.2))
    return pd.DataFrame(rows, columns=["connection", "h2", "signif"])


def manhattan_for_df(df, atlas, method, out_path):
    """
    Single-panel Manhattan y=h2 (0-1) — groups by network-pair.

    Network-network pairs accounting for <1.5% of total ROI-ROI connections
    are rescaled 100-fold on the x-axis and displayed in grey, exactly as
    in 03-matrixManhattan2.py. Large groups (>=1.5%) are shown first,
    ordered by median h2 descending; small groups are appended at the right
    and compressed, so the grey squished tail is on the right with no blank
    gaps between large groups.
    """
    total_edges = len(df)
    small_thresh = total_edges * 0.015  # 1.5% threshold
    
    # Calculate group statistics
    group_stats = df.groupby("connection").agg(
        median_h2=("h2", "median"),
        size=("h2", "size")
    ).sort_values("h2", ascending=False)  # Sort by median h2 descending
    
    # Identify small groups (< 1.5% of total edges)
    group_sizes = df.groupby("connection").size()
    is_small = {name: size < total_edges * 0.015 for name, size in group_sizes.items()}
    
    # connection_order: by median h2 descending (already sorted in group_stats)
    connection_order = group_stats.index.tolist()
    
    df["connection"] = pd.Categorical(df["connection"], categories=connection_order, ordered=True)
    df = df.sort_values("connection").reset_index(drop=True).reset_index(drop=False).rename(columns={"index": "idx"})
    df["index"] = df["idx"]
    grouped = df.groupby("connection", observed=True)
    total_edges = len(df)
    
    # Find divider: first index of first small group in the sorted data
    small_starts = [g["index"].min() for name, g in df.groupby("connection", observed=True) 
                    if (len(g) / total_edges) < 0.015]
    divider = min(small_starts) if small_starts else None
    
    # Top 30 groups get TABLEAU colors, rest grey
    nlabels = 30
    colors_large = np.fromiter(mcolors.TABLEAU_COLORS.values(), dtype="<U7")
    colors_large = np.tile(colors_large, math.ceil(nlabels / len(colors_large)))[:nlabels]
    colors_large = np.append(colors_large, np.repeat("grey", max(0, len(grouped) - nlabels)))
    
    fig, ax = plt.subplots(1, figsize=(12, 3.8))
    
    # divider is start of first small group
    small_starts = [g["index"].min() for name, g in grouped if (len(g) / total) < 0.015]
    divider = min(small_starts) if small_starts else None
    
    # Color mapping: top 30 get TABLEAU colors, rest grey
    color_idx = 0
    for idx, (name, group) in enumerate(grouped):
        is_small = (len(group) / total) < 0.015
        if is_small and divider is not None:
            # Compress 100x on x-axis, place on RIGHT
            group["index"] = (group["index"] - divider) / 100 + divider
            color = "grey"
            alpha, sz = 0.3, 6
        else:
            color = plt.get_cmap("tab20").colors[color_idx % 20]
            color_idx = (color_idx + 1) % 20
            alpha, sz = 0.9, 10
        
        ax.scatter(group["index"], group["h2"], color=color, s=10 if not is_small else 6,
                   alpha=0.9 if not is_small else 0.4, zorder=10, 
                   linewidths=0.3, edgecolors="black")
    
    ax.set_xlabel("")
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    try:
        xmax = max(g["index"].max() for _, g in grouped)
    except Exception:
        xmax = df["index"].max()
    ax.set_xlim([0, xmax * 1.02])
    ax.set_ylim([0, 1])
    disp = {"gordon": "Gordon", "probaConns": "ProbaConns", "SA": "SA"}.get(atlas, atlas)
    ax.set_ylabel(r"Heritability ($h^2$)", size=11)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_title(f"{disp} — {method}", fontsize=12, fontweight="bold")
    
    # Legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:blue", markersize=7, label="Large groups (top 30)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=7, label="Small groups (<1.5%, 100× compressed)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=False)
    
    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    n_small = sum((len(g) / total) < 0.015 for _, g in grouped)
    print(f"  wrote {out_path}  n={len(df)} groups={len(grouped)} small={n_small}")


wide = pd.read_csv(WIDE)
set_N = {"gordon": 352, "probaConns": 80, "SA": 17}
for atlas in ["gordon", "probaConns", "SA"]:
    N = set_N[atlas]
    for method, col in [("Twin", "Twin_h2"), ("AdjHE-RE", f"h2_{atlas}_AdjHE_RE" if atlas != "SA" else "h2_SA_AdjHE_RE")]:
        if col not in wide.columns:
            alt = col.replace("probaConns", "proba")
            if alt in wide.columns:
                col = alt
            else:
                print(f"skip {atlas} {method}: {col} not in wide")
                continue
        sub = wide[wide["Set"] == atlas]
        if sub[col].dropna().empty:
            print(f"skip {atlas} {method}: no data in {col}")
            continue
        df = build_long_for_set(wide, atlas, N, col)
        if df.empty:
            print(f"skip {atlas} {method}: empty df")
            continue
        out = PLOT_DIR / f"manhattan_{atlas}_{method.replace('-','')}.png"
        manhattan_for_df(df, atlas, method, out)

print("Done.")
