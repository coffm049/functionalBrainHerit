#!/usr/bin/env python3
"""Top-5 Gordon FC follow-up, part 1: distributions + White-subset inputs.

Top 5 by AdjHE-RE SNP h2 (30 PCs) from mash_twin_wide.csv. All 9 edges at the
h2=1.0 boundary are tied (identical var); these 5 break the tie by Twin h2.
All five are within-MTL edges over a tight parcel cluster (7/14/15, 174/176,
297, 139/262) — consistent with a parcel-level artifact rather than biology.

Run on HPC (has pconns.parquet + covars TSV):
    conda activate MASH   # or any env with pandas/pyarrow/matplotlib
    python -u summary/top5_fc_followup.py 2>&1 | tee /tmp/top5.log

Writes results/FCs/gordon/top5/:
    distribution_stats.csv, top5_distributions_raw.png,
    top5_distributions_resid.png
White-stratum moments are included per edge; the White-only GCTA run itself
is handled inside MASH via "covar_filter" (see top5_gcta_white.json).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Top 5: (edge, parcel_a, net_a, parcel_b, net_b, Twin_h2, SNP_h2)
TOP5 = [
    ("o46024", 174, "MTL", 176, "MTL", 0.347649, 1.0),
    ("o4984", 14, "MTL", 176, "MTL", 0.337098, 1.0),
    ("o46145", 174, "MTL", 297, "MTL", 0.259557, 1.0),
    ("o5105", 14, "MTL", 297, "MTL", 0.246904, 1.0),
    ("o39320", 139, "MTL", 262, "MTL", 0.231422, 1.0),
]
EDGES = [e for e, *_ in TOP5]

RANDO = Path("/projects/standard/rando149/coffm049/ABCD")
PHENO = RANDO / "Workflow/02_Phenotypes/FCsTopo/pconns.parquet"
COVAR = RANDO / "Results/02_Phenotypes/Covars_with_totalNetworkSurface.tsv"
EIGEN = RANDO / "Results/01_Gene_QC/filters/filter1/Eigens/no_rels.eigenvec"
OUTDIR = Path("/users/4/coffm049/papers/functionalBrainHerit/results/FCs/gordon/top5")
OUTDIR.mkdir(parents=True, exist_ok=True)


def moments(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    m, s = float(np.mean(x)), float(np.std(x, ddof=1))
    z = (x - m) / s if s > 0 else np.zeros_like(x)
    return {
        "n": n, "mean": m, "sd": s, "min": float(np.min(x)), "max": float(np.max(x)),
        "skew": float(np.mean(z ** 3)), "kurtosis": float(np.mean(z ** 4) - 3),
        "frac_absz_gt4": float(np.mean(np.abs(z) > 4)),
    }


def main():
    ph = pd.read_parquet(PHENO, columns=["IID"] + EDGES)
    print(f"pheno rows={len(ph)} cols={ph.columns.tolist()[:7]}...")
    co = pd.read_table(COVAR, sep="\t")
    print(f"covars shape={co.shape}")
    print("race_ethnicity levels:")
    print(co["race_ethnicity"].value_counts(dropna=False).to_string())

    df = ph.merge(co, on="IID", how="inner")
    print(f"pheno x covar overlap N={len(df)}")

    # ---- White (homogeneous-ancestry) subset ----
    if "White" not in set(co["race_ethnicity"].dropna().unique()):
        print("FATAL: 'White' not a race_ethnicity level; see levels above", file=sys.stderr)
        sys.exit(1)
    white = df[df["race_ethnicity"] == "White"].copy()
    print(f"White subset N={len(white)} (MASH subsets internally via covar_filter)")

    # ---- residualize on RE covariate set + 30 PCs (same as estimation) ----
    pc = pd.read_table(EIGEN, sep=r"\s+", header=None)
    pc = pc.iloc[:, :32]
    pc.columns = ["FID", "IID"] + [f"pc{i}" for i in range(1, 31)]
    df = df.merge(pc[["IID"] + [f"pc{i}" for i in range(1, 31)]], on="IID", how="left")

    rows = []
    resid = {}
    for e, a, na, b, nb, twin, snp in TOP5:
        m_raw = moments(df[e].values)
        m_raw = {"edge": e, "pa": a, "pb": b, "twin_h2": twin, "snp_h2": snp, **{f"raw_{k}": v for k, v in m_raw.items()}}
        # complete-case residualization: age + female + site dummies + pc1..30
        cols = ["age", "female"] + [f"pc{i}" for i in range(1, 31)]
        sub = df[["IID", e, "abcd_site"] + cols].dropna()
        site = pd.get_dummies(sub["abcd_site"], prefix="site", dtype=float)
        X = np.column_stack([np.ones(len(sub)), sub[cols].astype(float).values, site.values])
        y = sub[e].astype(float).values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ beta
        resid[e] = pd.Series(r, index=sub.index)
        m_res = moments(r)
        m_raw.update({f"resid_{k}": v for k, v in m_res.items()})
        m_raw["resid_n"] = len(sub)
        m_white = moments(df.loc[df["race_ethnicity"] == "White", e].values)
        m_raw.update({f"white_{k}": v for k, v in m_white.items()})
        rows.append(m_raw)
        print(f"{e} ({a}-{b} {na}-{nb}): raw skew={m_raw['raw_skew']:.2f} kurt={m_raw['raw_kurtosis']:.2f} "
              f"|z|>4={m_raw['raw_frac_absz_gt4']:.4f} | resid skew={m_res['skew']:.2f} kurt={m_res['kurtosis']:.2f} "
              f"| white skew={m_white['skew']:.2f} kurt={m_white['kurtosis']:.2f}")

    stats = pd.DataFrame(rows)
    stats.to_csv(OUTDIR / "distribution_stats.csv", index=False)
    print(f"Wrote {OUTDIR / 'distribution_stats.csv'}")

    # ---- plots ----
    fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharey=True)
    for ax, (e, *_r) in zip(axes, TOP5):
        v = df[e].dropna().values
        ax.hist(v, bins=60, color="#1f77b4", alpha=0.8)
        ax.set_title(f"{e}\nTwin={_r[-2]:.2f} SNP=1.0", fontsize=9)
        ax.set_xlabel("FC (raw)")
    axes[0].set_ylabel("count")
    fig.tight_layout()
    fig.savefig(OUTDIR / "top5_distributions_raw.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharey=True)
    for ax, (e, *_r) in zip(axes, TOP5):
        v = resid[e].values
        ax.hist(v, bins=60, color="#ff7f0e", alpha=0.8)
        ax.set_title(e, fontsize=9)
        ax.set_xlabel("FC (residualized)")
    axes[0].set_ylabel("count")
    fig.tight_layout()
    fig.savefig(OUTDIR / "top5_distributions_resid.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote distribution PNGs to {OUTDIR}")
    print("White-only GCTA needs no subset files: top5_gcta_white.json uses covar_filter.")


if __name__ == "__main__":
    main()
