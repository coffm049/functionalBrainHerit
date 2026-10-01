#!/usr/bin/env python3
"""
Rebuild the SNP analysis pool (filtered_ids) from the genomics-derived
no_rels GRM cohort, which is the authoritative sample definition.

The pool used to be narrowed further by intersecting with both phenotype
parquets and by intersecting with whatever the previous pool already
contained. Both were shrinking operations: the result could never grow, and
the parquet intersection capped the pool at probaConns' 5,552 subjects even
though no_rels defines 7,234.

The pool now comes from no_rels.grm.id alone. Each estimate intersects it
with its own phenotype and covariates while loading, so gordon can reach the
full 7,234 ceiling while probaConns stays limited by its own 5,552
phenotypes. Nothing is lost by not pre-intersecting: a subject with no
phenotype is dropped at load time anyway.

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

# Refuse to publish a pool smaller than this: no_rels defines 7,234 subjects,
# so anything this small means the id file was misparsed or mispointed.
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

    # Source of truth: the no_rels cohort, in GRM order.
    grm = grm_no_rels_df
    pool = grm[grm.IID.isin(grm_full)].copy()
    if pool.empty:
        raise ValueError(
            "Empty pool - no_rels.grm.id and full.grm.id share no subjects"
        )
    pool_ids = set(pool.IID)
    if len(pool) < MIN_POOL:
        raise ValueError(
            f"Pool of {len(pool)} subjects is below MIN_POOL={MIN_POOL} for a "
            f"no_rels.grm.id of {len(grm_no_rels)} parsed - refusing to write"
        )

    # Only worth backing up if the rebuild would actually drop somebody.
    if cur_ids - pool_ids and not os.path.isfile(BROAD):
        cur.to_csv(BROAD, sep="\t", header=False, index=False)
        print(f"Backed up previous pool -> {BROAD} ({len(cur)} subjects)")

    pool.to_csv(IDS, sep="\t", header=False, index=False)

    print(f"Current pool       : {len(cur_ids)} subjects")
    print(f"no_rels GRM        : {len(grm_no_rels)} IIDs")
    print(f"full GRM           : {len(grm_full)} IIDs (no_rels not in full: {len(grm_no_rels - grm_full)})")
    print(f"POOL (GRM only)    : {len(pool)} subjects -> wrote {IDS}")
    print(f"  added            : {len(pool_ids - cur_ids)}")
    print(f"  dropped          : {len(cur_ids - pool_ids)}")
    print("Parquet overlap (reported only - MASH intersects per phenotype):")
    print(f"  gordon pconns    : {len(pool_ids & gordon_ids)}")
    print(f"  probaConns       : {len(pool_ids & proba_ids)}")


if __name__ == "__main__":
    main()
