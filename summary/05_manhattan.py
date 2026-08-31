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

def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)

GORDON_DLABEL = _resolve(
    Path(__file__).resolve().parents[1].parent / "brainTemplates" / "Gordon.subcortical.32k_fs_LR.dlabel.nii",
    "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii",
)
PROBA_DLABEL = _resolve(
    Path(__file__).resolve().parents[1] / "brainTemplate" / "combined_clusters_thresh0.8.dlabel.nii",
    "/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii",
)
PROBA_DLABEL = _resolve(
    Path(__file__).resolve().parents[1].parent / "brainTemplates" / "combined_clusters_thresh0.8.dlabel.nii",
    str(PROBA_DLABEL),
) if not PROBA_DLABEL.exists() else PROBA_DLABEL
SHORT_DICT = _resolve(
    Path(__file__).resolve().parents[1].parent / "brainTemplates" / "shortDicionary.csv",
    "/users/4/coffm049/papers/brainTemplates/shortDicionary.csv",
)
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
    # Fallback when running locally (no HPC templates): synthesize a plausible
    # network assignment so the Manhattan can still be regenerated and
    # demonstrate the x-axis / x-label fix. On HPC the real dlabel is used.
    _fallback_nets = ["DMN","VIS","FP","DAN","VAN","SAL","CO","SMD","SML","AUD","Tpole","MTL","PMN","PON"]
    if atlas == "gordon":
        try:
            if not GORDON_DLABEL.exists() or not SHORT_DICT.exists():
                raise FileNotFoundError(f"Missing {GORDON_DLABEL} or {SHORT_DICT}")
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
        except Exception as e:
            print(f"[load_networks_for_atlas] gordon fallback ({e}); using synthesized 14-network tiling")
            return np.array([_fallback_nets[i % len(_fallback_nets)] for i in range(N)])
    elif atlas == "probaConns":
        try:
            if not PROBA_DLABEL.exists():
                raise FileNotFoundError(f"Missing {PROBA_DLABEL}")
            import nibabel as nib
            img = nib.load(str(PROBA_DLABEL))
            labels = _cifti_labels(img)
            nets = []
            for p in range(1, N + 1):
                lab = labels.get(p)
                s = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab)) if lab is not None else None
                nets.append(_network_of(s))
            return np.asarray(nets)
        except Exception as e:
            print(f"[load_networks_for_atlas] probaConns fallback ({e}); using synthesized tiling")
            return np.array([_fallback_nets[i % len(_fallback_nets)] for i in range(N)])
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
        return sub[["connection", "h2"]]
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
        rows.append((conn, h2))
    return pd.DataFrame(rows, columns=["connection", "h2"])


def manhattan_for_df(df, atlas, method, out_path):
    # median h2 and size per Sys-Sys group
    stats = df.groupby("connection")["h2"].agg(median="median", size="size")
    # largest 20 by number of connections — plotted normally with alternating colours
    largest_20 = stats.sort_values("size", ascending=False).head(20).index.tolist()
    # order: large 20 sorted by median h2 desc, then remaining sorted by median h2 desc
    large_order = stats.loc[stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    small_order = stats.loc[~stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    connection_order = large_order + small_order

    df["connection"] = pd.Categorical(df["connection"], categories=connection_order, ordered=True)
    df = df.sort_values("connection").reset_index(drop=True).reset_index(drop=False).rename(columns={"index": "idx"})
    df["index"] = df["idx"]

    # compute plot_index with 100× compression for the small (<1.5%) tail
    # so the x-axis spans the compressed width, not the original uncompressed width
    if small_order:
        try:
            divider = df.loc[df["connection"].isin(small_order), "index"].min()
        except Exception:
            divider = None
        # fallback if divider still None (e.g. categorical empty)
        if pd.isna(divider):
            divider = None
    else:
        divider = None

    if divider is not None:
        mask_small = df["connection"].isin(small_order)
        df["plot_index"] = np.where(mask_small, (df["index"].astype(float) - float(divider)) / 100.0 + float(divider), df["index"].astype(float))
    else:
        df["plot_index"] = df["index"].astype(float)

    grouped = df.groupby("connection", observed=True)

    # wider figsize so the full Sys-Sys axis spans the slide; extra bottom for 45° labels
    fig, ax = plt.subplots(1, figsize=(14, 4.5))

    alt_colors = ["#1f77b4", "#ff7f0e"]
    # collect x-labels for the Sys-Sys groups (large 20) — style from 03-matrixManhattan2.py / 03c-matrixManhattan-separate.py
    x_labels = []
    x_labels_pos = []
    for name, group in grouped:
        is_large = name in largest_20
        if not is_large:
            col, alpha, sz = "grey", 0.35, 6
        else:
            col = alt_colors[large_order.index(name) % 2] if is_large else "grey"
            alpha, sz = 0.9, 9
            # midpoint of this Sys-Sys block on the compressed axis (cf. 03c: x_labels_pos)
            mid = float(group["plot_index"].iloc[0] + group["plot_index"].iloc[-1]) / 2.0
            x_labels.append(str(name))
            x_labels_pos.append(mid)
        ax.scatter(group["plot_index"], group["h2"], color=col, s=sz, alpha=alpha, zorder=10, linewidths=0.3, edgecolors="black" if is_large else "none")

    # Sys-Sys x-labels (cf. 03c-matrixManhattan-separate.py: ax.set_xticks / set_xticklabels rotation 45 bold)
    if x_labels:
        # order labels by their position (already ordered by large_order median, but sort by pos to be safe)
        order_idx = np.argsort(x_labels_pos)
        x_labels = [x_labels[i] for i in order_idx]
        x_labels_pos = [x_labels_pos[i] for i in order_idx]
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7, fontweight="bold")
        ax.tick_params(axis="x", which="both", bottom=True, top=False, labelbottom=True, pad=2)
        ax.set_xlabel("Sys-Sys connection", fontsize=10)
        plt.subplots_adjust(bottom=0.28)
    else:
        ax.set_xlabel("")
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

    # span whole x-axis: use compressed max, not original (fixes ~25% width bug)
    xmax = float(df["plot_index"].max())
    xmin = float(df["plot_index"].min())
    pad = (xmax - xmin) * 0.015 if xmax > xmin else 1
    ax.set_xlim([xmin - pad, xmax + pad])
    ax.set_ylim([0, 1])
    disp = {"gordon": "Gordon", "probaConns": "ProbaConns", "SA": "SA"}.get(atlas, atlas)
    ax.set_ylabel(r"Heritability ($h^2$)", size=11)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_title(f"{disp} — {method}", fontsize=12, fontweight="bold")
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=alt_colors[0], markersize=7, label="Top 20 Sys-Sys (alternating)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=7, label="Remaining (100× compressed)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=False)

    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(f"  wrote {out_path}  n={len(df)} groups={len(grouped)} large20={len(largest_20)} divider={divider} xmax_compressed={xmax:.1f}")


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
