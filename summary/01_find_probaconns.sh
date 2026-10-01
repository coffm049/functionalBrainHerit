#!/bin/bash
# 01_find_probaconns.sh - locate probaConns .pconn.nii files that
# compileProbaCons.py never compiled (read-only)
#
# Usage: bash summary/01_find_probaconns.sh
#
# Why this exists:
#   compileProbaCons.SLURM:15-16 searched only two directories,
#     group1_10minonly_FD0p2 and group2_10minonly_FD0p2,
#   while the gordon equivalent, compilePcons.SLURM:16, searched the whole
#   abcd-hcp-pipeline tree. That asymmetry is why gordon's parquet has 10,123
#   subjects and probaConns' has 5,552, and probaConns is what caps the shared
#   filtered_ids pool. This script answers: are there eligible .pconn.nii
#   files on disk for subjects that proba.files never listed?
#
# Reports and writes nothing:
#   1. every directory group under threshold_template_maps/data with counts of
#      .nii, .pconn.nii, and *10*.pconn.nii - the "*10*" pattern being what
#      compileProbaCons.SLURM selects on (the 10-minute scan)
#   2. line counts of proba.files, pconn.files, FC.files
#   3. IIDs on disk that proba.files never listed (the recoverable set), and
#      listed paths that no longer exist on disk
#   4. how much of the recoverable set is usable in the shared pool - present
#      in filtered_ids.tsv and in gordon's pconns.parquet
#   5. matrix shape of a sample file per group, because compileProbaCons.py:25
#      raises unless the matrix gives exactly 80x80 -> 3160 edges, and
#      compileProbaCons.py:48 catches that and skips the subject silently
#
# Paths default to the production locations and can be overridden in the
# environment, which is how the fixture test exercises this script.
#
# Nothing here submits jobs or writes to results/.

set -uo pipefail
shopt -s nullglob

ROOT="${ROOT:-/users/4/coffm049/papers/functionalBrainHerit}"
MAPS="${MAPS:-/projects/standard/feczk001/shared/projects/ABCD/threshold_template_maps/data}"
FCSTOPO="${FCSTOPO:-/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo}"
POOL="${POOL:-/projects/standard/rando149/coffm049/filtered_ids.tsv}"
GORDON="${GORDON:-$FCSTOPO/pconns.parquet}"
PROBA="${PROBA:-$FCSTOPO/probaConns.parquet}"
# The compiled shape is fixed at N_REGIONS=80 -> 80*79/2 = 3160 edges
# (compileProbaCons.py:7-8).
N_REGIONS="${N_REGIONS:-80}"
want_edges=$(( N_REGIONS * (N_REGIONS - 1) / 2 ))

CONDA_SH="${CONDA_SH:-/users/4/coffm049/miniconda3/etc/profile.d/conda.sh}"

cd "$ROOT" || exit 1

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

hdr() { printf '\n===== %s =====\n' "$1"; }

# compileProbaCons.py:14-19, applied to basename(path).split("_",1)[0]
# (compileProbaCons.py:39). Replicated exactly so the disk/proba.files diff
# compares like with like rather than depending on my own idea of the id.
iid_of() {
  local s="${1##*/}"
  s="${s%%_*}"
  case "$s" in
    sub-*)     s="${s#sub-}" ;;
    NDAR_INV*) s="${s#NDAR_INV}" ;;
    NDARINV*)  s="${s#NDARINV}" ;;
  esac
  printf 'sub-%s\n' "$s"
}

if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
  conda activate MASH >/dev/null 2>&1 || true
fi
HAVE_ARROW=0
HAVE_NIB=0
python -c 'import pyarrow.parquet' >/dev/null 2>&1 && HAVE_ARROW=1
python -c 'import nibabel'         >/dev/null 2>&1 && HAVE_NIB=1

################################################################################
if [ ! -d "$MAPS" ]; then
  echo "MISSING data directory: $MAPS"
  exit 1
fi

hdr "1. threshold_template_maps/data INVENTORY"
# One pass: find every .nii, then classify by basename so the group directory
# name (which itself contains "10minonly") cannot be mistaken for a match.
find "$MAPS" -name '*.nii' 2>/dev/null | sed "s|^$MAPS/||" > "$tmp/nii_rel"
printf '%-42s %7s %7s %7s\n' GROUP NII PCONN PCONN10
awk -F/ -v OFS='\t' '
  { n = $NF
    g = (NF > 1 ? $1 : "(top level)")
    total[g]++
    if (n ~ /\.pconn\.nii$/) pconn[g]++
    if (n ~ /10/ && n ~ /\.pconn\.nii$/) p10[g]++ }
  END { for (g in total) printf "%s\t%d\t%d\t%d\n", g, total[g]+0, pconn[g]+0, p10[g]+0 }' \
  "$tmp/nii_rel" | sort -t$'\t' -k4,4nr -k2,2nr > "$tmp/groups"

if [ ! -s "$tmp/groups" ]; then
  echo "no .nii files found under $MAPS"
fi
while IFS=$'\t' read -r g n p q; do
  printf '%-42s %7d %7d %7d\n' "$g" "$n" "$p" "$q"
done < "$tmp/groups"

tot_nii=$(awk -F'\t' '{ s += $2 } END { print s + 0 }' "$tmp/groups")
tot_p10=$(awk -F'\t' '{ s += $4 } END { print s + 0 }' "$tmp/groups")
printf '%-42s %7d %7s %7d\n' TOTAL "$tot_nii" "-" "$tot_p10"
echo
echo "compileProbaCons.SLURM:15-16 searched only group1_10minonly_FD0p2 and"
echo "group2_10minonly_FD0p2. Any other group with PCONN10 > 0 was never listed."

################################################################################
hdr "2. FILE LISTS"
for f in proba.files pconn.files FC.files; do
  path="$FCSTOPO/$f"
  if [ -f "$path" ]; then
    n=$(grep -c . "$path" || true)
    when=$(date -r "$path" '+%Y-%m-%d %H:%M' 2>/dev/null || echo '?')
    printf '%-14s %8d lines   modified %s\n' "$f" "$n" "$when"
  else
    printf '%-14s %8s\n' "$f" MISSING
  fi
done

################################################################################
hdr "3. ELIGIBLE ON DISK vs proba.files"
if [ ! -f "$FCSTOPO/proba.files" ]; then
  echo "MISSING $FCSTOPO/proba.files - cannot diff"
  exit 1
fi

# The exact pattern compileProbaCons.SLURM passes to find -name, matched
# against the basename only, as find does.
awk '{ n = $0; sub(/^.*\//, "", n)
       if (n ~ /10/ && n ~ /\.pconn\.nii$/) print }' "$tmp/nii_rel" \
  > "$tmp/disk_eligible_rel"
n_disk=$(grep -c . "$tmp/disk_eligible_rel" || true)

: > "$tmp/disk_iids"
while IFS= read -r rel; do
  [ -n "$rel" ] || continue
  iid_of "$rel" >> "$tmp/disk_iids"
done < "$tmp/disk_eligible_rel"
sort -u "$tmp/disk_iids" -o "$tmp/disk_iids"

: > "$tmp/listed_paths"
while IFS= read -r p; do
  [ -n "$p" ] || continue
  printf '%s\n' "$p" >> "$tmp/listed_paths"
done < "$FCSTOPO/proba.files"
sort -u "$tmp/listed_paths" -o "$tmp/listed_paths"

: > "$tmp/listed_iids"
while IFS= read -r p; do
  iid_of "$p" >> "$tmp/listed_iids"
done < "$tmp/listed_paths"
sort -u "$tmp/listed_iids" -o "$tmp/listed_iids"

n_listed=$(grep -c . "$tmp/listed_paths" || true)
n_listed_iid=$(grep -c . "$tmp/listed_iids" || true)

# Stale entries: listed in proba.files but gone from disk.
: > "$tmp/stale"
while IFS= read -r p; do
  [ -f "$p" ] || printf '%s\n' "$p" >> "$tmp/stale"
done < "$tmp/listed_paths"
n_stale=$(grep -c . "$tmp/stale" || true)

comm -23 "$tmp/disk_iids" "$tmp/listed_iids" > "$tmp/new_iids"
comm -13 "$tmp/disk_iids" "$tmp/listed_iids" > "$tmp/listed_not_disk_iids"
n_new=$(grep -c . "$tmp/new_iids" || true)
n_orphan=$(grep -c . "$tmp/listed_not_disk_iids" || true)

printf 'eligible *.pconn.nii paths on disk : %8d\n' "$n_disk"
printf '  distinct IIDs                    : %8d\n' "$(grep -c . "$tmp/disk_iids" || true)"
printf 'paths listed in proba.files        : %8d\n' "$n_listed"
printf '  distinct IIDs                    : %8d\n' "$n_listed_iid"
printf 'listed but missing from disk       : %8d\n' "$n_stale"
printf 'RECOVERABLE (on disk, not listed)  : %8d IIDs\n' "$n_new"
printf 'listed but no eligible file on disk: %8d IIDs\n' "$n_orphan"

if [ "$n_new" -eq 0 ]; then
  echo
  echo "Nothing to recover: proba.files already covers every eligible file."
else
  echo
  echo "sample of recoverable IIDs:"
  head -10 "$tmp/new_iids" | sed 's/^/    /'
fi

################################################################################
hdr "4. IS THE RECOVERABLE SET USABLE? (genotyped + gordon-phenotyped)"
in_pool=0
if [ -f "$POOL" ]; then
  # filtered_ids is FID\tIID (prep_iids.sh:140 keeps column 2).
  cut -f2 "$POOL" | sed 's/[[:space:]]*$//' | sort -u > "$tmp/pool_iids"
  comm -12 "$tmp/new_iids" "$tmp/pool_iids" > "$tmp/new_in_pool"
  in_pool=$(grep -c . "$tmp/new_in_pool" || true)
  printf 'in filtered_ids.tsv (genotyped)    : %8d of %d\n' "$in_pool" "$n_new"
else
  printf '%-34s: %s\n' "filtered_ids.tsv" "MISSING $POOL"
fi

in_gordon=0
if [ "$HAVE_ARROW" -eq 1 ] && [ -f "$GORDON" ]; then
  python - "$GORDON" <<'PY' > "$tmp/gordon_iids" 2>/dev/null || true
import sys, pyarrow.parquet as pq
t = pq.read_table(sys.argv[1], columns=["IID"])
sys.stdout.write("\n".join(t.column("IID").to_pylist()))
PY
  if [ -s "$tmp/gordon_iids" ]; then
    sort -u "$tmp/gordon_iids" -o "$tmp/gordon_iids"
    comm -12 "$tmp/new_iids" "$tmp/gordon_iids" > "$tmp/new_in_gordon"
    in_gordon=$(grep -c . "$tmp/new_in_gordon" || true)
    printf 'in pconns.parquet (gordon)         : %8d of %d\n' "$in_gordon" "$n_new"
  else
    echo "could not read IID column from $GORDON"
  fi
else
  [ "$HAVE_ARROW" -eq 0 ] && echo "pyarrow unavailable - gordon check skipped"
  [ ! -f "$GORDON" ] && echo "MISSING $GORDON - gordon check skipped"
fi

in_proba=0
if [ "$HAVE_ARROW" -eq 1 ] && [ -f "$PROBA" ]; then
  python - "$PROBA" <<'PY' > "$tmp/proba_iids" 2>/dev/null || true
import sys, pyarrow.parquet as pq
t = pq.read_table(sys.argv[1], columns=["IID"])
sys.stdout.write("\n".join(t.column("IID").to_pylist()))
PY
  if [ -s "$tmp/proba_iids" ]; then
    sort -u "$tmp/proba_iids" -o "$tmp/proba_iids"
    comm -12 "$tmp/new_iids" "$tmp/proba_iids" > "$tmp/new_in_proba"
    in_proba=$(grep -c . "$tmp/new_in_proba" || true)
    printf 'already in probaConns.parquet      : %8d of %d\n' "$in_proba" "$n_new"
  fi
fi

################################################################################
hdr "5. MATRIX SHAPE (compileProbaCons.py needs $N_REGIONS x $N_REGIONS = $want_edges edges)"
if [ "$HAVE_NIB" -eq 0 ]; then
  echo "nibabel unavailable - shape check skipped."
  echo "  Run this under the MASH conda env to enable it."
else
  printf '%-42s %-14s %s\n' GROUP SHAPES "verdict"
  while IFS=$'\t' read -r g n p q; do
    [ "${q:-0}" -gt 0 ] || continue
    sample=$(awk -F/ -v g="$g" '
      $1 == g { n = $NF
                if (n ~ /10/ && n ~ /\.pconn\.nii$/) { print; exit } }' \
      "$tmp/nii_rel")
    [ -n "$sample" ] || continue
    shapes=$(python - "$MAPS/$sample" <<'PY' 2>/dev/null
import sys, nibabel as nib, numpy as np
d = nib.load(sys.argv[1]).get_fdata()
v = d[np.triu_indices(d.shape[0], 1)].shape[0]
print(f"{d.shape[0]}x{d.shape[1]}->edges={v}")
PY
)
    if [ -z "$shapes" ]; then
      verdict="READ_FAIL"
    elif [ "$shapes" = "${N_REGIONS}x${N_REGIONS}->edges=$want_edges" ]; then
      verdict=OK
    else
      verdict="WRONG_SHAPE (compileProbaCons.py:25 would raise; :48 swallows it)"
    fi
    printf '%-42s %-14s %s\n' "$g" "${shapes:-?}" "$verdict"
  done < "$tmp/groups"
fi

################################################################################
hdr "6. SUMMARY"
echo "probaConns subjects compiled : $n_listed_iid  (proba.files IIDs)"
echo "eligible on disk             : $(grep -c . "$tmp/disk_iids" || true) IIDs"
echo "recoverable                  : $n_new IIDs"
echo "  of which genotyped         : $in_pool"
echo "  of which gordon-phenotyped : $in_gordon"
echo "  already compiled           : $in_proba"
echo
if [ "$in_pool" -gt 0 ]; then
  echo "Recovering those subjects raises the shared filtered_ids pool by up to"
  echo "  $in_pool subjects, for every atlas and every estimator, because the"
  echo "  pool is the 4-way intersection (update_filtered_ids.py)."
  echo
  echo "Next step if the numbers justify it:"
  echo "  1. widen compileProbaCons.SLURM:15-16 to every group with PCONN10 > 0"
  echo "  2. resubmit sbatch compileProbaCons.SLURM  (rebuilds probaConns.parquet)"
  echo "  3. rerun: python update_filtered_ids.py"
else
  echo "No recoverable subject is in the genotyped pool, so a recompile would"
  echo "not change filtered_ids.tsv. probaConns stays capped by its own data."
fi

printf '\n%s\n' "Sweep complete. This script writes nothing and submits nothing."
