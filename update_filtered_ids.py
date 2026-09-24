#!/usr/bin/env python3
"""
Compute the common analysis pool (filtered_ids) across the GRM-derived
filtered_ids, the gordon parquet, and the probaConns parquet, so every SNP
estimate uses one shared subject set. Backs up the current (broader) pool to
filtered_ids_broad.tsv for twin analyses, which keep a larger sample.

Output format matches what MASH expects: headerless, tab-separated, FID IID.
"""
import os
import pyarrow.parquet as pq
import pandas as pd

P = "/projects/standard/rando149/coffm049"
IDS = f"{P}/filtered_ids.tsv"
BROAD = f"{P}/filtered_ids_broad.tsv"
GORDON = f"{P}/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet"
PROBA = f"{P}/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet"

# Every SNP config must be a subset of the GRM it reads. Configs use either the
# full or the no_rels GRM, so intersect with both id files to be safe.
GRM_NO_RELS = f"{P}/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/no_rels/no_rels.grm.id"
GRM_FULL = f"{P}/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full.grm.id"
GRMS = [GRM_NO_RELS, GRM_FULL]


def parquet_iids(path):
    tbl = pq.read_table(path, columns=["IID"])
    return set(tbl.column("IID").to_pylist())


def grm_iids(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing GRM id file {path}")
    df = pd.read_table(path, header=None, names=["FID", "IID"], dtype=str)
    return set(df.IID)


def main():
    if not os.path.isfile(IDS):
        raise FileNotFoundError(f"Missing {IDS}")
    for f in (GORDON, PROBA):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing {f}")
    for f in GRMS:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing GRM id file {f}")

    cur = pd.read_table(IDS, header=None, names=["FID", "IID"], dtype=str)
    cur_ids = set(cur.IID)
    gordon_ids = parquet_iids(GORDON)
    proba_ids = parquet_iids(PROBA)

    grm_id_sets = {os.path.basename(os.path.dirname(g)): grm_iids(g) for g in GRMS}
    common = cur_ids & gordon_ids & proba_ids & set.intersection(*grm_id_sets.values())
    if not common:
        raise ValueError("Empty common pool - check parquet IID columns")

    # Preserve original (GRM) order for the common pool
    keep = cur["IID"].isin(common)
    common_df = cur.loc[keep, ["FID", "IID"]].copy()

    if not os.path.isfile(BROAD):
        cur.to_csv(BROAD, sep="\t", header=False, index=False)
        print(f"Backed up broader pool -> {BROAD} ({len(cur)} subjects)")

    common_df.to_csv(IDS, sep="\t", header=False, index=False)
    print(f"Current pool      : {len(cur_ids)} subjects")
    print(f"Gordon parquet    : {len(gordon_ids)} unique IIDs (in current: {len(cur_ids & gordon_ids)})")
    print(f"ProbaConns parquet: {len(proba_ids)} unique IIDs (in current: {len(cur_ids & proba_ids)})")
    for name, s in grm_id_sets.items():
        print(f"GRM {name:<10s}: {len(s)} IIDs (in current: {len(cur_ids & s)})")
    print(f"COMMON pool       : {len(common)} subjects -> wrote {IDS}")


if __name__ == "__main__":
    main()