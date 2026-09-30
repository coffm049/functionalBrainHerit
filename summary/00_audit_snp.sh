#!/bin/bash
# 00_audit_snp.sh - audit MASH SNP heritability outputs (read-only)
#
# Usage: bash summary/00_audit_snp.sh
#
# Per method/Set stream reports: file count, row-size histogram, array index
# gaps, row-size anomalies, N used, flag distribution, and a KEEP / INCOMPLETE /
# STALE / ANOMALY / RERUN verdict.
#
# Expected shapes derived from the submit scripts:
#   FCs/gordon/submit.sh:16-27      TOTAL=61776, CHUNK=208 (FE/RE/HEreg), 100 (GCTA)
#   FCs/probaConns/submit.sh:16-27  TOTAL=3160,  CHUNK=10 (FE/GCTA/HEreg), 40 (RE)
#   make_configs.py:26-27           last chunk is short when total % chunk != 0
#   SA/*.json                       17 phenotypes, one CSV per method
#
# Rows per file = chunk phenotypes x len(npc): MASH emits one row per
# (phenotype x npc value), tagged by the 'PCs' column
# (MASH/src/Estimate/estimators/all_estimators.py:455, itertools.product).
# The npc field is read from the template each submit.sh selects:
#   gordon fe/re npc=[0,20] -> x2   proba re  npc=[0,20] -> x2
#   gordon gcta/he npc=[20] -> x1   others   npc=[20]    -> x1
# 01_compare_mash_twin.R:33-36 already collapses these by filtering PCs == npc.
#
# This script never writes to results/ and never submits jobs.

set -uo pipefail
shopt -s nullglob

ROOT=/users/4/coffm049/papers/functionalBrainHerit
POOL=/projects/standard/rando149/coffm049/filtered_ids.tsv
cd "$ROOT" || exit 1

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# key|glob|chunk|total|expected_files|npc_multiplier
STREAMS=(
  "gordon_AdjHE_FE|results/FCs/gordon/pconns.AdjHE.FE.*.csv|208|61776|297|2"
  "gordon_AdjHE_RE|results/FCs/gordon/pconns.AdjHE.RE.*.csv|208|61776|297|2"
  "gordon_GCTA|results/FCs/gordon/pconns.GCTA.GCTA.*.csv|100|61776|618|1"
  "gordon_HEreg|results/FCs/gordon/pconns.HEreg.HEreg.*.csv|208|61776|297|1"
  "proba_AdjHE_FE|results/FCs/probaConns/probaConns.AdjHE.FE.*.csv|10|3160|316|1"
  "proba_AdjHE_RE|results/FCs/probaConns/probaConns.AdjHE.RE.*.csv|40|3160|79|2"
  "proba_GCTA|results/FCs/probaConns/probaConns.GCTA.GCTA.*.csv|10|3160|316|1"
  "proba_HEreg|results/FCs/probaConns/probaConns.HEreg.HEreg.*.csv|10|3160|316|1"
)

# pass stream table to R
: > "$tmp/streams.tsv"
for spec in "${STREAMS[@]}"; do
  IFS='|' read -r key glob chunk total expect mult <<< "$spec"
  printf '%s\t%s\n' "$key" "$glob" >> "$tmp/streams.tsv"
done

################################################################################
hdr() { printf '\n===== %s =====\n' "$1"; }

################################################################################
hdr "1. STREAM INVENTORY"
printf '%-18s %6s %7s %8s  %s\n' STREAM FILES EXPECT EXPROWS "ROW-SIZE HISTOGRAM (count x rows)"
for spec in "${STREAMS[@]}"; do
  IFS='|' read -r key glob chunk total expect mult <<< "$spec"
  files=( $glob )
  n=${#files[@]}
  if [ "$n" -eq 0 ]; then
    printf '%-18s %6d %7d %8s  MISSING\n' "$key" 0 "$expect" "$(( chunk * mult ))"
    continue
  fi
  hist=$(for f in "${files[@]}"; do echo $(( $(wc -l < "$f") - 1 )); done \
         | sort -n | uniq -c | awk '{printf "%dx%s ", $1, $2}')
  printf '%-18s %6d %7d %8d  %s\n' "$key" "$n" "$expect" "$(( chunk * mult ))" "$hist"
done

################################################################################
hdr "2. SA FILES (one CSV per method, 17 phenotypes)"
found_sa=0
for f in results/SA/*.csv; do
  [ -e "$f" ] || continue
  found_sa=1
  printf '%-36s %5d rows  %s\n' "$f" "$(( $(wc -l < "$f") - 1 ))" "$(date -r "$f" '+%Y-%m-%d %H:%M')"
done
[ "$found_sa" -eq 0 ] && echo "no SA csv files found"

################################################################################
hdr "3. INDEX GAPS AND ROW-SIZE ANOMALIES"
for spec in "${STREAMS[@]}"; do
  IFS='|' read -r key glob chunk total expect mult <<< "$spec"
  files=( $glob )
  [ ${#files[@]} -gt 0 ] || continue

  inv="$tmp/$key.inv"
  : > "$inv"
  for f in "${files[@]}"; do
    b=$(basename "$f" .csv)
    printf '%s\t%s\t%s\n' "${b##*.}" "$(( $(wc -l < "$f") - 1 ))" "$f" >> "$inv"
  done
  sort -n -k1,1 "$inv" > "$inv.s"

  exp_rows=$(( chunk * mult ))
  tailrows=$(( (total % chunk) * mult ))

  gaps=$(awk -F'\t' '
    prev != "" && $1+0 != prev+1 { printf "%d-%d ", prev+1, $1-1 }
    { prev = $1+0 }' "$inv.s")
  printf '%s' "$gaps" > "$tmp/$key.gaps"

  if [ -n "$gaps" ]; then
    printf '%-18s GAPS (missing array indices): %s\n' "$key" "$gaps"
  else
    printf '%-18s no gaps\n' "$key"
  fi

  anom=$(awk -F'\t' -v c="$exp_rows" -v t="$tailrows" \
    '$2+0 != c+0 && (t == 0 || $2+0 != t+0)' "$inv.s")
  if [ -n "$anom" ]; then
    printf '%-18s ANOMALY (expect %s phenos x %s npc = %s rows, tail %s):\n' \
      "$key" "$chunk" "$mult" "$exp_rows" "$tailrows"
    acount=$(printf '%s\n' "$anom" | wc -l)
    printf '%s\n' "$anom" | head -n 10 |
      while IFS=$'\t' read -r ix rw path; do
        printf '    idx %-6s %6s rows  %s  %s\n' "$ix" "$rw" "$path" "$(date -r "$path" '+%Y-%m-%d %H:%M')"
      done
    [ "$acount" -gt 10 ] && printf '    ... and %d more\n' "$((acount - 10))"
    printf '%s\n' "$anom" > "$tmp/$key.anom"
  fi
done

################################################################################
hdr "4. N USED AND FLAG DISTRIBUTION (first 5 files per stream)"
if [ -f "$POOL" ]; then
  pool_n=$(wc -l < "$POOL")
  echo "current filtered_ids.tsv N = $pool_n"
else
  pool_n=""
  echo "WARNING: pool not found at $POOL - N comparison will be skipped"
fi
echo

cat > "$tmp/naudit.R" <<'RSCRIPT'
suppressMessages({library(readr); library(dplyr)})
root <- Sys.getenv("AUDIT_ROOT")
sl  <- read.delim(Sys.getenv("AUDIT_TSV"), header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("key", "pattern"))
out <- file(Sys.getenv("AUDIT_OUT"), open = "wt")
for (i in seq_len(nrow(sl))) {
  fs <- Sys.glob(file.path(root, sl$pattern[i]))
  if (!length(fs)) next
  fs <- sort(fs)[seq_len(min(5, length(fs)))]
  d <- tryCatch(suppressMessages(bind_rows(lapply(fs, read_csv, show_col_types = FALSE))),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) {
    cat(sprintf("%-18s READ_ERROR\n", sl$key[i]))
    cat(sprintf("%s\tNA\n", sl$key[i]), file = out)
    next
  }
  nvals <- if ("n" %in% names(d)) {
    u <- unique(na.omit(as.character(d$n))); paste(u, collapse = ",")
  } else "NA"
  fl <- if ("flag" %in% names(d)) {
    tb <- table(d$flag)
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " ")
  } else "no flag column"
  h2med <- if ("h2" %in% names(d)) sprintf("%.3f", median(d$h2, na.rm = TRUE)) else "NA"
  cat(sprintf("%-18s sampled=%4d rows  N=%-10s med_h2=%-7s flags: %s\n",
              sl$key[i], nrow(d), nvals, h2med, fl))
  cat(sprintf("%s\t%s\n", sl$key[i], nvals), file = out)
}
close(out)
RSCRIPT

(
  source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh 2>/dev/null
  conda activate gdc 2>/dev/null
  AUDIT_ROOT="$ROOT" AUDIT_TSV="$tmp/streams.tsv" AUDIT_OUT="$tmp/n_used.tsv" \
    Rscript --vanilla "$tmp/naudit.R"
)
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$tmp/n_used.tsv" ]; then
  echo "R audit unavailable (rc=$rc); N/flag checks skipped."
else
  echo
fi

################################################################################
hdr "5. VERDICTS"
printf '%-18s %-12s %s\n' STREAM VERDICT DETAIL
for spec in "${STREAMS[@]}"; do
  IFS='|' read -r key glob chunk total expect mult <<< "$spec"
  files=( $glob )
  n=${#files[@]}

  if [ "$n" -eq 0 ]; then
    printf '%-18s %-12s %s\n' "$key" "RERUN" "no output files"
    continue
  fi

  gaps=""
  [ -f "$tmp/$key.gaps" ] && gaps=$(cat "$tmp/$key.gaps")

  obs_n=""
  if [ -s "$tmp/n_used.tsv" ]; then
    obs_n=$(awk -F'\t' -v k="$key" '$1 == k { print $2; exit }' "$tmp/n_used.tsv")
  fi

  verdict="KEEP"; detail="$n/$expect files"
  if [ "$n" -ne "$expect" ]; then
    verdict="INCOMPLETE"; detail="$n/$expect files"
  fi
  if [ -n "$gaps" ]; then
    verdict="INCOMPLETE"; detail="$detail, gaps $gaps"
  fi
  if [ -n "$pool_n" ] && [ -n "$obs_n" ] && [ "$obs_n" != "NA" ] && [ "$obs_n" != "$pool_n" ]; then
    verdict="STALE"; detail="$detail, N=$obs_n vs pool=$pool_n"
  fi
  if [ -f "$tmp/$key.anom" ]; then
    acount=$(wc -l < "$tmp/$key.anom")
    detail="$detail, $acount row-size anomalies (see section 3)"
    [ "$verdict" = "KEEP" ] && verdict="ANOMALY"
  fi
  printf '%-18s %-12s %s\n' "$key" "$verdict" "$detail"
done

printf '\n%s\n' "Audit complete. Rerun this after each submission batch."
