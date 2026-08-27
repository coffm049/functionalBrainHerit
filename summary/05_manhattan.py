#!/usr/bin/env python3
"""Pairwise-system Manhattan plots — Python, faithful to 03-matrixManhattan2.py
and 03c-matrixManhattan-separate.py (groups by network pairs).

Groups edges by network-pair (e.g. DMN-DMN, DMN-VIS) using the parcel->network
mapping from the Gordon / ProbaConns dlabel label tables (via shortDicionary
for Gordon). Produces the same 2-panel figure as the 01-07 examples:

  top:    Manhattan scatter  y = h2  vs  x = edge index grouped by
          network-pair (ordered by twin mean h2 descending, as in
          connectionOrder = df.query(Type==twin).groupby(connection).mean)
  bottom: bar  y = proportion heritable (signif mean) per network-pair,
          first 30 pairs coloured (TABLEAU), tail compressed.

For every phenotype Set (gordon 352, probaConns 80, SA 17) and both methods
(Twin, AdjHE-RE) from mash_twin_wide.csv (30 PCs).

Outputs: results/summary/plots/manhattan_{Set}_{Method}.png

Run on HPC where the dlabels and .venv are available:
  .venv/bin/python summary/05_manhattan.py
"""
import itertools
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

# Gordon parcel->network via label table + shortDicionary (as in 02a_brain_gordon.py)
GORDON_DLABEL = Path(
    "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii"
)
PROBA_DLABEL = Path(
    "/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/"
    "ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"
)
SHORT_DICT = Path(
    "/users/4/coffm049/papers/brainTemplates/shortDicionary.csv"
)

# SA 17-network short names (as in 02c_brain_sa.py / 05-topoViz.R)
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
    """Return array networks[N] (network name per ROI index 0..N-1)."""
    import nibabel as nib
    if atlas == "gordon":
        df = pd.read_csv(SHORT_DICT)
        name_to_short = dict(
            zip(df["name"].astype(str), df["shortname"].astype(str))
        )
        img = nib.load(str(GORDON_DLABEL))
        labels = _cifti_labels(img)
        nets = []
        for p in range(1, N + 1):
            lab = labels.get(p)
            s = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab)) if lab is not None else None
            net = _network_of(s)
            nets.append(name_to_short.get(net, net))
            # fallback case-insensitive
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
    else:  # SA 17 networks are themselves
        return np.array([SA_NETWORKS.get(i + 1, f"NET{i+1}") for i in range(N)])


def decode_pheno_to_ij(pheno, N):
    """Pheno like 'o123' or '123' -> (i,j) via triu_indices(N,1)."""
    s = str(pheno)
    k = int(s[1:]) if s.startswith("o") else int(s)
    i, j = np.triu_indices(N, 1)
    return int(i[k]), int(j[k])


def build_long_for_set(wide, atlas, N, h2_col, twin_col="Twin_h2"):
    """Return long df with columns connection, h2, signif for one atlas/method."""
    sub = wide[wide["Set"] == atlas][["Pheno", h2_col]].dropna()
    # keep only phenos that decode correctly for this N
    sub = sub.copy()
    sub["k"] = sub["Pheno"].apply(
        lambda p: int(str(p)[1:]) if str(p).startswith("o") else int(str(p)) if str(p).isdigit() else None
    )
    # for SA, Pheno is network_surfarea* not oK — handle separately
    if atlas == "SA":
        # SA is per-network, not pairwise; connection is just the network name
        sub["connection"] = sub["Pheno"].apply(
            lambda p: SA_NETWORKS.get(int(str(p).split("network_surfarea")[-1]) if "network_surfarea" in str(p) else 0, str(p))
        )
        sub["h2"] = sub[h2_col].astype(float)
        sub["signif"] = sub["h2"] > 0.2
        return sub[["connection", "h2", "signif"]]

    # FC: pairwise
    networks = load_networks_for_atlas(atlas, N)
    rows = []
    for _, r in sub.iterrows():
        pheno = r["Pheno"]
        h2 = float(r[h2_col])
        try:
            i, j = decode_pheno_to_ij(pheno, N)
            ni = networks[i] if i < len(networks) else "NA"
            nj = networks[j] if j < len(networks) else "NA"
            conn = f"{ni}-{nj}"
            # canonicalise order: e.g. Aud-DMN vs DMN-Aud -> same group; sort?
            # 03 script keeps i-j order as per triu, not sorted, so keep as is
        except Exception:
            conn = str(pheno)
        rows.append((conn, h2, h2 > 0.2))
    df = pd.DataFrame(rows, columns=["connection", "h2", "signif"])
    return df


def manhattan_for_df(df, atlas, method, out_path):
    """Single-panel Manhattan y=h2 (0-1) — groups by network-pair.

    Network-network pairs accounting for <1.5% of total ROI-ROI connections
    are rescaled 100-fold on the x-axis and displayed in grey, exactly as
    in 03-matrixManhattan2.py (original used nlabels=30 + pROIs=0.015).
    """
    connection_order = (
        df.groupby("connection")["h2"].mean().sort_values(ascending=False).index.tolist()
    )
    df["connection"] = pd.Categorical(
        df["connection"], categories=connection_order, ordered=True
    )
    df = df.sort_values("connection").reset_index(drop=True).reset_index(drop=False)
    df = df.rename(columns={"index": "idx"})
    df["index"] = df["idx"]

    grouped = df.groupby("connection", observed=True)
    total = len(df)
    # 1.5% threshold for small groups (as per description)
    small_thresh = 0.015
    is_small = {name: (len(g) / total) < small_thresh for name, g in grouped}
    # Tableau colours for large groups, grey for small
    nlabels = 30
    colors_large = np.fromiter(mcolors.TABLEAU_COLORS.values(), dtype="<U7")
    colors_large = np.tile(colors_large, math.ceil(nlabels / len(colors_large)))[:nlabels]

    fig, ax = plt.subplots(1, figsize=(9, 3.5))

    # divider is the start of the first small group in sorted order
    small_starts = [g["index"].min() for name, g in grouped if is_small.get(name, False)]
    divider = min(small_starts) if small_starts else None

    color_idx = 0
    for name, group in grouped:
        small = is_small.get(name, False)
        if small and divider is not None:
            group["index"] = (group["index"] - divider) / 100 + divider
            col = "grey"
        else:
            col = colors_large[color_idx % len(colors_large)]
            color_idx += 1
        group.plot(kind="scatter", x="index", y="h2", color=col,
                   ax=ax, s=10 if not small else 6, zorder=10, alpha=0.9 if not small else 0.6)

    ax.set_xlabel("")
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    # xlim must cover the compressed tail
    try:
        xmax = df["index"].max()
        # df["index"] was modified in-place for small groups; recompute max from plotted data
        xmax = max(g["index"].max() for _, g in grouped) if len(grouped) else xmax
    except Exception:
        xmax = df["index"].max()
    ax.set_xlim([0, xmax * 1.02])
    ax.set_ylim([0, 1])
    disp = {"gordon": "Gordon", "probaConns": "ProbaConns", "SA": "SA"}.get(atlas, atlas)
    ax.set_ylabel(r"Heritability ($h^2$)", size=11)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_title(f"{disp} — {method}", fontsize=12, fontweight="bold")

    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    n_small = sum(is_small.values())
    print(f"  wrote {out_path}  n={len(df)} groups={len(grouped)} small(<1.5%)={n_small}")


# --------------------------------------------------------------------------
# Main: loop over all Sets and both methods
# --------------------------------------------------------------------------
wide = pd.read_csv(WIDE)

# Map Set -> N and h2 columns
set_N = {"gordon": 352, "probaConns": 80, "SA": 17}
for atlas in ["gordon", "probaConns", "SA"]:
    N = set_N[atlas]
    for method, col in [("Twin", "Twin_h2"),
                        ("AdjHE-RE", f"h2_{atlas}_AdjHE_RE" if atlas != "SA" else "h2_SA_AdjHE_RE")]:
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
