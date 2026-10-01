#!/usr/bin/env python3
"""
Rebuild the SNP analysis pool (filtered_ids) as one identical subject set.

The pool is the intersection of the two genomics GRMs (no_rels and full,
because every SNP config reads one of those two and passes this pool as
"ids") with both phenotype parquets (gordon pconns and probaConns). That
intersection is the whole point of the file: it keeps the same subjects in
every atlas and every estimator, so a difference between two methods or two
atlases is a difference in method or atlas and not in sample.

Intersecting from the GRM outward rather than narrowing whatever the
previous pool happened to be means the result can grow as well as shrink --
if probaConns is recompiled with more subjects, so is the pool.

Per-phenotype loading still drops subjects missing covariates, so a subject
in the shared pool may appear in fewer than all four estimators of an
atlas. That residual difference is unavoidable and is reported by the audit.

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

# Every SNP config reads either the full or the no_rels GRM and passes this
# pool as "ids", so a subject must be in both id files to be usable anywhere.
GRM_NO_RELS = f"{P}/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/no_rels/no_rels.grm.id"
GRM_FULL = f"{P}/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full.grm.id"

# Refuse to publish a pool smaller than this: the intersection has run near
# 3,2xx subjects, so anything this small means a parquet column or an id file
# was misparsed rather than genuinely sparse.
MIN_POOL = 1000


def parquet_iids(path):
    tbl = pq.read_table(path, columns=["IID"])
    return set(tbl.column("IID").to_pylist())


def read_ids(path):
    """Read a headerless FID/IID file.

    summary/prep_iids.sh:98 accepts exactly a tab or a space in these files
    and verifies the column counts that way, so match it: take whichever
    delimiter yields two fields on every non-blank line. Anything else fails
    loudly with its field-count histogram instead of handing back a
    shaped-but-wrong frame -- a silent misparse here would republish the
    pool from garbage.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing id file {path}")

    # Default universal-newline reading also splits lone-CR files correctly.
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip() for ln in fh]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValueError(f"{path}: file is empty")

    for split in (lambda s: s.split("\t"), str.split, lambda s: s.split(",")):
        parts = [split(ln) for ln in lines]
        widths = {len(p) for p in parts}
        if min(widths) >= 2:
            # FID and IID are always columns 1 and 2, the same fields
            # prep_iids.sh:140 selects (iid col 2); trailing extras are
            # ignored rather than treated as a parse failure.
            if widths != {2}:
                print(f"  note: {path}: {sorted(widths)} fields/line; "
                      f"using columns 1-2 (FID, IID)")
            return pd.DataFrame(
                {"FID": [p[0] for p in parts], "IID": [p[1] for p in parts]},
                dtype=str,
            )

    hist = {}
    for ln in lines:
        w = len(ln.split())
        hist[w] = hist.get(w, 0) + 1
    bad = next((ln for ln in lines if len(ln.split()) != 2), lines[0])
    raise ValueError(
        f"{path}: expected exactly 2 'FID IID' fields per line.\n"
        f"  {len(lines)} non-blank lines; whitespace field-count histogram "
        f"{dict(sorted(hist.items()))}\n"
        f"  first line          : {lines[0]!r}\n"
        f"  first line that is not 2 fields: {bad!r}"
    )


def main():
    if not os.path.isfile(IDS):
        raise FileNotFoundError(f"Missing {IDS}")
    for f in (GORDON, PROBA):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing {f}")
    for f in (GRM_NO_RELS, GRM_FULL):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing GRM id file {f}")

    cur = read_ids(IDS)
    cur_ids = set(cur.IID)
    gordon_ids = parquet_iids(GORDON)
    proba_ids = parquet_iids(PROBA)
    grm_no_rels_df = read_ids(GRM_NO_RELS)
    grm_full_df = read_ids(GRM_FULL)
    grm_no_rels = set(grm_no_rels_df.IID)
    grm_full = set(grm_full_df.IID)

    # Source of truth: the no_rels cohort, in GRM order, narrowed to subjects
    # usable by every SNP config (present in both GRMs) and phenotyped in both
    # atlases. Base is the GRM rather than the previous pool so the set is
    # recomputed correctly instead of only ever shrinking.
    grm = grm_no_rels_df
    usable = grm_full & gordon_ids & proba_ids
    pool = grm[grm.IID.isin(usable)].copy()
    if pool.empty:
        raise ValueError(
            "Empty pool - no_rels.grm.id shares no subjects with full GRM, "
            f"pconns.parquet, and probaConns.parquet "
            f"(full={len(grm_full)}, gordon={len(gordon_ids)}, "
            f"proba={len(proba_ids)})"
        )
    pool_ids = set(pool.IID)
    if len(pool) < MIN_POOL:
        raise ValueError(
            f"Pool of {len(pool)} subjects is below MIN_POOL={MIN_POOL} - "
            f"no_rels parsed {len(grm_no_rels)}, full {len(grm_full)}, "
            f"gordon {len(gordon_ids)}, proba {len(proba_ids)} - refusing to write"
        )

    # Only worth backing up if the rebuild would actually drop somebody.
    if cur_ids - pool_ids and not os.path.isfile(BROAD):
        cur.to_csv(BROAD, sep="\t", header=False, index=False)
        print(f"Backed up previous pool -> {BROAD} ({len(cur)} subjects)")

    pool.to_csv(IDS, sep="\t", header=False, index=False)

    print(f"Current pool       : {len(cur_ids)} subjects")
    print(f"no_rels GRM        : {len(grm_no_rels)} IIDs")
    print(f"full GRM           : {len(grm_full)} IIDs (no_rels not in full: {len(grm_no_rels - grm_full)})")
    print(f"gordon pconns      : {len(gordon_ids)} IIDs")
    print(f"probaConns         : {len(proba_ids)} IIDs")
    print(f"POOL (4-way shared): {len(pool)} subjects -> wrote {IDS}")
    print(f"  added            : {len(pool_ids - cur_ids)}")
    print(f"  dropped          : {len(cur_ids - pool_ids)}")
    print("Ceilings (how much each constraint alone costs):")
    print(f"  no_rels only             : {len(grm[grm.IID.isin(grm_full)])}")
    print(f"  no_rels x gordon         : {len(grm[grm.IID.isin(grm_full & gordon_ids)])}")
    print(f"  no_rels x proba          : {len(grm[grm.IID.isin(grm_full & proba_ids)])}")
    print("  after covariate load-time drops: see 00_audit_snp.sh")


if __name__ == "__main__":
    main()
