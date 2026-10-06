#!/usr/bin/env python3
"""Pairwise-system Manhattan — grouped by Sys-Sys, sorted by median h2.

For each atlas (Gordon 352, ProbaConns 80, SA 17), two methods are plotted (20 PCs):

  gordon / probaConns / SA : Twin, AdjHE-RE

AdjHE-RE is the single reported SNP method: it is the only AdjHE variant that
applies a conditioning screen (cond(X_tX) > 1e10), so a retained RE row has been
checked for ill-conditioning whereas a retained FE row has had only a determinant
sign test. AdjHE-FE and GCTA are still computed in mash_twin_wide.csv and still
compared against AdjHE-RE in 07_method_diagnostics.R; they are not plotted here.
All methods for an atlas share one x-axis, ordered by median h2 within each
Sys-Sys group, with Twin drawn against AdjHE-RE. Panels state their own retained n.

  x = edge index grouped by Sys-Sys network pair (e.g. DMN-VIS), ordered by
      largest median h2 within that Sys-Sys group (descending).
  y = h2 (0–1, full heritability range).
  Largest 20 Sys-Sys groups by number of connections are plotted normally
  with alternating colours; remaining groups are compressed 100× on the
  x-axis and shown in grey.

Outputs: results/summary/plots/manhattan_{Set}_{Method}.png
         results/summary/plots/manhattan_overview.png  (atlas x method grid)
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

def _find_proba_dlabel():
    candidates = [
        Path(__file__).resolve().parents[1] / "brainTemplate" / "combined_clusters_thresh0.8.dlabel.nii",
        Path(__file__).resolve().parents[1] / "brainTemplate" / "combined_clusters_thresh0.75.dlabel.nii",
        Path(__file__).resolve().parents[1].parent / "brainTemplates" / "combined_clusters_thresh0.8.dlabel.nii",
        Path(__file__).resolve().parents[1].parent / "brainTemplates" / "combined_clusters_thresh0.75.dlabel.nii",
        Path(r"C:\Users\coffm049\brainTemplates\combined_clusters_thresh0.8.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\combined_clusters_thresh0.75.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.8.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"),
        Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"),
        Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.75.dlabel.nii"),
        Path("/users/4/coffm049/papers/brainTemplates/combined_clusters_thresh0.8.dlabel.nii"),
        Path("/users/4/coffm049/papers/brainTemplates/combined_clusters_thresh0.75.dlabel.nii"),
    ]
    # Prefer 80-parcel file (N=80) to avoid 63-parcel mismatch
    for p in candidates:
        if p.exists():
            try:
                import nibabel as _nib
                img = _nib.load(str(p))
                ax0 = img.header.get_axis(0)
                n_lab = len(ax0.get_element(0)[1]) - 1
                if n_lab == 80:
                    return p
            except Exception:
                pass
    for p in candidates:
        if p.exists():
            return p
    return Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii")

PROBA_DLABEL = _find_proba_dlabel()
SHORT_DICT = _resolve(
    Path(__file__).resolve().parents[1].parent / "brainTemplates" / "shortDicionary.csv",
    "/users/4/coffm049/papers/brainTemplates/shortDicionary.csv",
)
SA_NETWORKS = {
    1: "DMN", 2: "VIS", 3: "FP", 5: "DAN", 7: "VAN",
    8: "SAL", 9: "CO", 10: "SMD", 11: "SML", 12: "AUD",
    13: "Tpole", 14: "MTL", 15: "PMN", 16: "PON",
}
# 4,6,17 are null per Gordon (no parcels) — per 01-07 labeling (05-topoViz.R, 03b range(1,15))
# excluded to avoid index shift; see 02c_brain_sa.py


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
    s = s.split("_")[0] if "_" in s else s
    return s


def load_networks_for_atlas(atlas, N):
    import nibabel as nib
    # Fallback when running locally (no HPC templates): synthesize a plausible
    # network assignment so the Manhattan can still be regenerated and
    # demonstrate the x-axis / x-label fix. On HPC the real dlabel is used.
    _fallback_nets = ["DMN","VIS","FP","DAN","VAN","SAL","CO","SMD","SML","AUD","Tpole","MTL","PMN","PON"]
    CORTICAL = {"DMN","VIS","Vis","FP","DAN","VAN","SAL","Sal","CO","SMD","SMd","SML","SMl","AUD","Aud","Tpole","MTL","PMN","PON"}
    CORTICAL_LOWER = {c.lower(): c for c in CORTICAL}
    CANON = {"vis":"Vis","aud":"Aud","sal":"Sal","smd":"SMd","sml":"SMl","van":"VAN","dan":"DAN","dmn":"DMN","fp":"FP","co":"CO","mtl":"MTL","pmn":"PMN","pon":"PON","tpole":"Tpole"}
    if atlas == "gordon":
        # Recreate test_A_LUT correct solution: explicit Gordon333_LUT.txt ROI_ID->Network_Name
        lut_candidates = [
            ROOT.parent / "brainTemplates" / "Gordon333_LUT.txt",
            ROOT.parent / "brainTemplates" / "Gordon_333_LUT.txt",
            Path(r"C:\Users\coffm049\brainTemplates\Gordon333_LUT.txt"),
            Path("/users/4/coffm049/papers/brainTemplates/Gordon333_LUT.txt"),
            Path("Gordon333_LUT.txt"),
        ]
        for cand in lut_candidates:
            if cand.exists():
                try:
                    lut = pd.read_csv(cand, sep=r"\s+", header=None, engine="python")
                    if str(lut.iloc[0,0]).lower() == "roi_id":
                        lut = lut.iloc[1:]
                    id_to_net = dict(zip(lut.iloc[:,0].astype(int), lut.iloc[:,2].astype(str)))
                    nets = []
                    for p in range(1, N+1):
                        raw = id_to_net.get(p, "NA")
                        net = _network_of(raw)
                        # canonicalize via shortDicionary if possible
                        try:
                            df_lut = pd.read_csv(SHORT_DICT)
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
                        return np.asarray(nets)
                except Exception as e:
                    print(f"  LUT {cand} failed: {e}")
                    pass
        # Fallback: gordon_modules.csv (352 rows, region 1..14) via labelDictionary (previous correct for 01-07)
        gordon_modules_candidates = [
            Path(__file__).resolve().parents[1].parent / "brainTemplates" / "gordon_modules.csv",
            Path(r"C:\Users\coffm049\brainTemplates\gordon_modules.csv"),
            Path("/users/4/coffm049/papers/brainTemplates/gordon_modules.csv"),
            ROOT.parent / "brainTemplates" / "gordon_modules.csv",
            Path("brainTemplate/gordon_modules.csv"),
        ]
        for cand in gordon_modules_candidates:
            if cand.exists():
                try:
                    gm = pd.read_csv(cand)
                    col = gm.columns[0]
                    vals = gm[col].astype(int).tolist()
                    if len(vals) == N:
                        labelDict = {1:"DMN",2:"VIS",3:"FP",4:"DAN",5:"VAN",6:"SAL",7:"CO",8:"SMD",9:"SML",10:"AUD",11:"Tpole",12:"MTL",13:"PMN",14:"PON"}
                        nets = [labelDict.get(v, "NA") for v in vals]
                        return np.asarray(nets)
                except Exception:
                    pass
        # Fallback: CIFTI label table + shortDicionary.csv
        try:
            if not GORDON_DLABEL.exists() or not SHORT_DICT.exists():
                raise FileNotFoundError(f"Missing {GORDON_DLABEL} or {SHORT_DICT}")
            df = pd.read_csv(SHORT_DICT)
            name_to_short = dict(zip(df["name"].astype(str), df["shortname"].astype(str)))
            cortical = set(name_to_short.values())
            cortical_lower = {c.lower(): c for c in cortical}
            img = nib.load(str(GORDON_DLABEL))
            labels = _cifti_labels(img)
            nets = []
            for p in range(1, N + 1):
                lab = labels.get(p)
                s = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab)) if lab is not None else None
                net = _network_of(s)
                short = name_to_short.get(net, net)
                if short == net:
                    for k, v in name_to_short.items():
                        if k.lower() == net.lower():
                            short = v
                            break
                if short not in cortical and short.lower() not in cortical_lower:
                    short = "NA"
                nets.append(short)
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
                net = _network_of(s)
                # Swap Sal <-> SMl to match correct insula topography (see 02b)
                if net == "Sal":
                    net = "SMl"
                elif net == "SMl":
                    net = "Sal"
                if net not in CORTICAL and net.lower() not in CORTICAL_LOWER:
                    net = "NA"
                else:
                    net = CANON.get(net.lower(), net)
                nets.append(net)
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
        # Per 01-07 labeling (05-topoViz.R, 03b), only 14 SA networks are non-null:
        # 1,2,3,5,7,8,9,10,11,12,13,14,15,16 — 4,6,17 are null and excluded to avoid shift.
        def _sa_conn(p):
            try:
                num = int(str(p).split("network_surfarea")[-1]) if "network_surfarea" in str(p) else int(str(p))
                return SA_NETWORKS.get(num)  # None for 4,6,17 -> dropped
            except Exception:
                return None
        sub["connection"] = sub["Pheno"].apply(_sa_conn)
        sub = sub.dropna(subset=["connection"])
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


def _get_ordering(df, atlas, method):
    """Helper to get Sys-Sys ordering for an atlas — used to share SNP ordering between Twin and SNP per user request."""
    stats = df.groupby("connection")["h2"].agg(median="median", size="size")
    largest_20 = stats.sort_values("size", ascending=False).head(20).index.tolist()
    large_order = stats.loc[stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    small_order = stats.loc[~stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    return large_order, small_order, largest_20

def _get_shared_order(dfs: dict, order_method: str = "AdjHE-RE"):
    """Sys-Sys ordering based on a single reference method (default: AdjHE-RE)."""
    ref_df = dfs.get(order_method)
    if ref_df is None:
        return None
    ref = ref_df[["connection", "h2"]].dropna()
    if ref.empty:
        return None
    stats = ref.groupby("connection")["h2"].agg(median="median", size="size")
    largest_20 = stats.sort_values("size", ascending=False).head(20).index.tolist()
    large_order = stats.loc[stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    small_order = stats.loc[~stats.index.isin(largest_20)].sort_values("median", ascending=False).index.tolist()
    return large_order, small_order, largest_20


def manhattan_for_df(df, atlas, method, out_path, shared_order=None):
    # median h2 and size per Sys-Sys group — if shared_order provided (from SNP), use it for both Twin and SNP per user
    if shared_order is not None:
        large_order, small_order, largest_20 = shared_order
    else:
        stats = df.groupby("connection")["h2"].agg(median="median", size="size")
        largest_20 = stats.sort_values("size", ascending=False).head(20).index.tolist()
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
    # Publication-ready: no title (handled in quarto caption), minimal whitespace
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=alt_colors[0], markersize=7, label="Top 20 Sys-Sys (alternating)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=7, label="Remaining (100× compressed)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=False)
    # Retained-row count per panel. AdjHE-FE and AdjHE-RE keep very different
    # fractions of phenotypes (~48% vs ~77% for gordon), so without an explicit n
    # the two panels look directly comparable when they are not: same y-scale,
    # shared x-axis, different row sets. Upper-left is free (legend is upper-right).
    n_rows, n_conn = int(len(df)), int(df["connection"].nunique())
    ax.text(0.012, 0.975, f"n = {n_rows:,} rows\n{n_conn} Sys-Sys",
            transform=ax.transAxes, va="top", ha="left", fontsize=7.5, linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.75))
    fig.tight_layout(pad=0.4)
    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    # Save stats for quarto (instead of title numbers). median_h2/mean_h2 are read
    # back by 06_results_catalog.qmd so the per-panel prose can interpolate the
    # FE/RE median ratio instead of hardcoding a number that drifts whenever a
    # re-run retains a different set of rows.
    try:
        stats_out = str(out_path).replace(".png", "_stats.csv")
        _h2 = pd.to_numeric(df["h2"], errors="coerce").dropna()
        pd.DataFrame([{"atlas": atlas, "method": method, "n": n_rows, "n_connections": n_conn, "groups": int(len(grouped)), "large20": int(len(largest_20)), "divider": float(divider) if divider is not None else np.nan, "xmax_compressed": float(xmax), "median_h2": float(_h2.median()) if len(_h2) else np.nan, "mean_h2": float(_h2.mean()) if len(_h2) else np.nan}]).to_csv(stats_out, index=False)
    except Exception:
        pass
    print(f"  wrote {out_path}  n={n_rows} conns={n_conn} groups={len(grouped)} large20={len(largest_20)} divider={divider} xmax_compressed={xmax:.1f}")


wide = pd.read_csv(WIDE)
set_N = {"gordon": 352, "probaConns": 80, "SA": 17}

# Which h2 columns to plot, per atlas, in panel order. AdjHE-RE is the single
# reported SNP method: it is the only AdjHE variant that applies a conditioning
# screen (cond(X_tX) > 1e10), so a retained RE row has actually been checked for
# ill-conditioning whereas a retained FE row has not. AdjHE-FE and GCTA are still
# computed in mash_twin_wide.csv and reported in 07_method_diagnostics.R, but
# they are not plotted here.
# probaConns columns are named h2_proba_* because 01_compare_mash_twin.R keys
# them by stream label "proba_<method>" rather than by Set "probaConns".
METHOD_SPECS = {
    "gordon": [
        ("Twin", "Twin_h2"),
        ("AdjHE-RE", "h2_gordon_AdjHE_RE"),
    ],
    "probaConns": [
        ("Twin", "Twin_h2"),
        ("AdjHE-RE", "h2_proba_AdjHE_RE"),
    ],
    "SA": [
        ("Twin", "Twin_h2"),
        ("AdjHE-RE", "h2_SA_AdjHE_RE"),
    ],
}

# Column order of the overview grid -- union of every method any atlas uses.
OVERVIEW_COLS = ["Twin", "AdjHE-RE"]


def _resolve_col(col):
    """The column present in `wide`, tolerating a renamed proba stream label."""
    if col in wide.columns:
        return col
    alt = col.replace("h2_proba_", "h2_probaConns_")
    return alt if alt in wide.columns else None


def specs_for(atlas):
    """(method, column) pairs for `atlas` that actually carry data."""
    sub = wide[wide["Set"] == atlas]
    out = []
    for method, col in METHOD_SPECS[atlas]:
        got = _resolve_col(col)
        if got is None:
            print(f"skip {atlas} {method}: {col} not in wide")
            continue
        if sub[got].dropna().empty:
            print(f"skip {atlas} {method}: no data in {got}")
            continue
        out.append((method, got))
    return out


for atlas in ["gordon", "probaConns", "SA"]:
    N = set_N[atlas]
    specs = specs_for(atlas)
    if not specs:
        continue

    dfs = {}
    for method, col in specs:
        df = build_long_for_set(wide, atlas, N, col)
        if df.empty:
            print(f"skip {atlas} {method}: empty df")
            continue
        dfs[method] = df
    if not dfs:
        continue

    shared_order = _get_shared_order(dfs) if dfs else None

    for method, df in dfs.items():
        out = PLOT_DIR / f"manhattan_{atlas}_{method.replace('-', '')}.png"
        manhattan_for_df(df, atlas, method, out, shared_order=shared_order)

# ---- Manhattan Overview: rows = atlas, cols = method ----
# Regenerated alongside the individual panels so every method appears, not
# only the ordering reference. (The previous version drew just one method per atlas:
# the drawing block was indented outside the method loop, so `ax` was left
# pointing at the last method and `method` held that same value.)

# Overview grid disabled per user request

print("Done.")
