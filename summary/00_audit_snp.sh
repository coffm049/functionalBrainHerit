#!/bin/bash
# 00_audit_snp.sh - audit MASH SNP heritability outputs (read-only)
#
# Usage: bash summary/00_audit_snp.sh
#
# Per method/Set stream reports: file count, row-size histogram, array index
# gaps, row-size anomalies, N used, flag distribution, and a KEEP / INCOMPLETE /
# STALE / ANOMALY / RERUN verdict.
#
# KEEP+ means coverage only: section 4 produced no N output, so staleness is
# UNVERIFIED. Do not treat KEEP+ as proof the results were computed on the
# current sample.
#
# Expected shapes derived from the submit scripts:
#   FCs/gordon/submit.sh:16-27      TOTAL=61776, CHUNK=208 (FE/RE/HEreg), 100 (GCTA)
#   FCs/probaConns/submit.sh:16-27  TOTAL=3160,  CHUNK=10 (FE/GCTA/HEreg), 40 (RE)
#   make_configs.py:26-27           last chunk is short when total % chunk != 0
#   SA/*.json                       18 mpheno (17 network_surfarea +
#                                    anthro_height_calc), one CSV per method
#
# Rows per file = chunk phenotypes x len(npc): MASH emits one row per
# (phenotype x npc value), tagged by the 'PCs' column
# (MASH/src/Estimate/estimators/all_estimators.py:455, itertools.product).
# The npc field is read from the template each submit.sh selects. Every active
# template is npc=[20] (one row per phenotype) after the AdjHE-RE fix:
# random_groups was entering the PC projection in _adjhe_3comp, which made
# QSQ = Q D D' Q vanish identically (S is built from the same site indicator
# matrix D), zeroing XtX's middle row and flagging 100% of RE rows
# ill_conditioned at every npc and sample size. See
# MASH tests/test_adjhe_site_projection.py.
# 01_compare_mash_twin.R:33-36 selects PCs == use_npc, which falls back to
# max(avail) = 20.
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
  "gordon_AdjHE_FE|results/FCs/gordon/pconns.AdjHE.FE.*.csv|208|61776|297|1"
  "gordon_AdjHE_RE|results/FCs/gordon/pconns.AdjHE.RE.*.csv|208|61776|297|1"
  "gordon_GCTA|results/FCs/gordon/pconns.GCTA.GCTA.*.csv|100|61776|618|1"
  "gordon_HEreg|results/FCs/gordon/pconns.HEreg.HEreg.*.csv|208|61776|297|1"
  "proba_AdjHE_FE|results/FCs/probaConns/probaConns.AdjHE.FE.*.csv|10|3160|316|1"
  "proba_AdjHE_RE|results/FCs/probaConns/probaConns.AdjHE.RE.*.csv|40|3160|79|1"
  "proba_GCTA|results/FCs/probaConns/probaConns.GCTA.GCTA.*.csv|10|3160|316|1"
  "proba_HEreg|results/FCs/probaConns/probaConns.HEreg.HEreg.*.csv|10|3160|316|1"
)

################################################################################
hdr() { printf '\n===== %s =====\n' "$1"; }

# Rows a SA run must have = len(npc) x len(mpheno), read from its template so
# the check stays correct when a template gains a phenotype. Handles both the
# multi-line ("npc": [ / 20 / ]) and inline ("npc": [20]) spellings; a line
# that only opens the array has to fall through to the digit rule below.
sa_expected() {
  awk '
    /"npc"[[:space:]]*:/    { sec = "npc"
                              s = $0; sub(/.*\[/, "", s); sub(/\].*/, "", s)
                              if (s ~ /[0-9]/) { npc += split(s, a, /,/); sec = "" }
                              next }
    /"mpheno"[[:space:]]*:/ { sec = "ph"
                              s = $0; sub(/^[^:]*:[[:space:]]*/, "", s); sub(/\].*$/, "", s)
                              if (s ~ /"/) { ph += gsub(/"/, "", s) / 2; sec = "" }
                              next }
    /^[[:space:]]*\]/       { sec = "";    next }
    sec == "npc" && /^[[:space:]]*[0-9]/ { npc++ }
    sec == "ph"  && /^[[:space:]]*"/     { ph++ }
    END { print (npc + 0) * (ph + 0) }
  ' "$1"
}

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
hdr "2. SA FILES (expected rows = len(npc) x len(mpheno) from SA/<name>.json)"
printf '%-24s %5s %5s  %-9s  %s\n' FILE ROWS EXP STATUS MODIFIED
: > "$tmp/sa.tsv"
found_sa=0
for f in results/SA/*.csv; do
  [ -e "$f" ] || continue
  found_sa=1
  b=$(basename "$f" .csv)
  rows=$(( $(wc -l < "$f") - 1 ))
  when=$(date -r "$f" '+%Y-%m-%d %H:%M')
  tpl="SA/$b.json"
  exp=""
  [ -f "$tpl" ] && exp=$(sa_expected "$tpl")
  status="OK"
  if [ ! -f "$tpl" ]; then
    status="NO_TPL"
  elif ! [ "${exp:-0}" -gt 0 ] 2>/dev/null; then
    status="BAD_TPL"
  elif [ "$rows" -ne "$exp" ]; then
    status="STALE"
  fi
  printf '%-24s %5d %5s  %-9s  %s\n' "$b" "$rows" "${exp:--}" "$status" "$when"
  printf '%s\t%s\t%s\t%s\n' "$b" "$rows" "${exp:--}" "$status" >> "$tmp/sa.tsv"
done
[ "$found_sa" -eq 0 ] && echo "no SA csv files found"

# A template with no output was never run (or was cleaned out).
for tpl in SA/*.json; do
  [ -e "$tpl" ] || continue
  b=$(basename "$tpl" .json)
  [ -f "results/SA/$b.csv" ] && continue
  exp=$(sa_expected "$tpl")
  printf '%-24s %5s %5s  %-9s\n' "$b" "-" "$exp" "NO_OUTPUT"
  printf '%s\t%s\t%s\t%s\n' "$b" "-" "$exp" "NO_OUTPUT" >> "$tmp/sa.tsv"
done

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
hdr "4. N USED AND FLAG DISTRIBUTION (first 5 FC files per stream, all SA csvs)"
if [ -f "$POOL" ]; then
  pool_n=$(wc -l < "$POOL")
  pool_n="${pool_n// /}"
  pool_mtime=$(date -r "$POOL" '+%Y-%m-%d %H:%M')
  echo "current filtered_ids.tsv: N=$pool_n  modified $pool_mtime"
else
  pool_n=""
  pool_mtime=""
  echo "WARNING: pool not found at $POOL - N comparison will be skipped"
fi
echo

# Sample size column is 'N' (uppercase), set per method at
# MASH/src/Estimate/estimators/all_estimators.py:170,175,193,201; the
# phenotype column is 'pheno'. N is per-phenotype and can vary slightly with
# missingness, so track min/max rather than assuming one value.
# Flags are counted PER npc value (the 'PCs' column) so a rerun that changes
# npc can actually be judged: ill_conditioned at npc=0 but ok at npc=20 is a
# completely different verdict from ill_conditioned at both.
# Deliberately awk and not R: tying this check to a conda env made the STALE
# verdict silently skippable whenever R failed to launch.
flag_summary() {
  local key=$1
  shift
  [ $# -gt 0 ] || return 0
  awk -F',' -v key="$key" -v mp="$tmp/n_used.tsv" -v bp="$tmp/flags.tsv" '
    FNR == 1 {
      if (++fn > 5) exit
      ncol = fcol = pcol = 0
      for (i = 1; i <= NF; i++) {
        h = $i; gsub(/["\r]/, "", h)
        if (h == "N") ncol = i
        else if (h == "flag") fcol = i
        else if (h == "PCs") pcol = i
      }
      next
    }
    {
      rows++
      if (ncol) {
        v = $ncol; gsub(/["\r]/, "", v); v += 0
        if (v > 0) {
          if (nmax == 0 || v > nmax) nmax = v
          if (nmin == 0 || v < nmin) nmin = v
        }
      }
      if (fcol) {
        v = $fcol; gsub(/["\r]/, "", v)
        p = "?"
        if (pcol) {
          p = $pcol; gsub(/["\r]/, "", p)
          if (p == "") p = "?"
        }
        k = p "|" v
        if (k in fv) fv[k]++
        else { fv[k] = 1; ord[++no] = k }
      }
    }
    END {
      ns = "NA"
      if (nmax > 0) ns = (nmin == nmax) ? nmax : nmin "-" nmax
      fs = ""
      bad = 0
      for (i = 1; i <= no; i++) {
        split(ord[i], a, /\|/)
        fs = fs (fs == "" ? "" : ", ") "npc" a[1] ":" a[2] "=" fv[ord[i]]
        # same bad-flag list 01_compare_mash_twin.R:21 drops on
        if (a[2] ~ /(ill_conditioned|singular|nan_solve|h2_gt_1_invalid|nonpos_det|neg_sigma_g)/)
          bad += fv[ord[i]]
      }
      if (fs == "") fs = "no flag column"
      printf "%-20s sampled=%4d rows  N=%-12s flags: %s\n", key, rows, ns, fs
      printf "%s\t%s\n", key, ns >> mp
      printf "%s\t%d\t%d\n", key, bad, rows >> bp
    }' "$@"
}

: > "$tmp/flags.tsv"
for spec in "${STREAMS[@]}"; do
  IFS='|' read -r key glob chunk total expect mult <<< "$spec"
  files=( $glob )
  [ ${#files[@]} -gt 0 ] || continue
  sample=()
  for f in "${files[@]}"; do
    sample+=("$f")
    [ ${#sample[@]} -ge 5 ] && break
  done
  flag_summary "$key" "${sample[@]}"
done

# SA streams are read by 01 too and carry the same flag column.
for f in results/SA/*.csv; do
  [ -e "$f" ] || continue
  flag_summary "SA/$(basename "$f" .csv)" "$f"
done

# Secondary signal: CSVs written before the pool file's last edit cannot have
# used the pool's current contents.
if [ -n "$pool_mtime" ]; then
  n_total=$(find results/FCs results/SA -name '*.csv' 2>/dev/null | wc -l)
  n_pre=$(find results/FCs results/SA -name '*.csv' ! -newer "$POOL" 2>/dev/null | wc -l)
  echo
  echo "CSVs older than filtered_ids.tsv: $n_pre / $n_total  (predate the pool's last edit)"
  if [ "$n_pre" -gt 0 ]; then
    echo "  sample of those:"
    find results/FCs results/SA -name '*.csv' ! -newer "$POOL" 2>/dev/null | head -5 | sed 's/^/    /'
  fi
fi

# The awk loop above writes n_used.tsv; if it produced nothing, staleness
# cannot be evaluated and the verdicts must say so rather than look clean.
rc=0
[ -s "$tmp/n_used.tsv" ] || rc=1
if [ "$rc" -ne 0 ]; then
  echo "N check produced no output: coverage verified, staleness was NOT."
fi

################################################################################
hdr "5. VERDICTS"
printf '%-18s %-12s %s\n' STREAM VERDICT DETAIL

# n_checked is 0 only when section 4 produced no n_used.tsv at all. In that
# case nothing below can detect staleness, so say so rather than printing a
# clean KEEP that implies sample-size verification.
n_checked=0
[ "$rc" -eq 0 ] && [ -s "$tmp/n_used.tsv" ] && n_checked=1

# Does every sampled row carry a flag that 01_compare_mash_twin.R:22 drops?
# A stream can be complete and not stale yet still contribute nothing.
FLAG_TXT=""; FLAG_VERDICT=""
flag_detail() {
  local k=$1 line fbad ftot
  FLAG_TXT=""; FLAG_VERDICT=""
  [ -s "$tmp/flags.tsv" ] || return 0
  line=$(awk -F'\t' -v k="$k" '$1 == k { print $2 "\t" $3; exit }' "$tmp/flags.tsv")
  [ -n "$line" ] || return 0
  fbad=${line%%$'\t'*}
  ftot=${line##*$'\t'}
  case "$fbad$ftot" in *[!0-9]*) return 0 ;; esac
  [ "$ftot" -gt 0 ] || return 0
  if [ "$fbad" -eq "$ftot" ]; then
    FLAG_VERDICT="ALL_FLAGGED"
    FLAG_TXT="all $ftot sampled rows flagged, 01 drops them"
  elif [ "$fbad" -gt 0 ]; then
    FLAG_TXT="$fbad/$ftot sampled rows flagged"
  fi
}

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
  if [ "$n_checked" -eq 1 ]; then
    obs_n=$(awk -F'\t' -v k="$key" '$1 == k { print $2; exit }' "$tmp/n_used.tsv")
  fi

  verdict="KEEP"; detail="$n/$expect files"
  if [ "$n" -ne "$expect" ]; then
    verdict="INCOMPLETE"; detail="$n/$expect files"
  fi
  if [ -n "$gaps" ]; then
    verdict="INCOMPLETE"; detail="$detail, gaps $gaps"
  fi
  # Results are stale only if they used MORE subjects than the current pool
  # contains - that cannot be a subset of the current pool, so the run must
  # predate it. Fewer subjects than the pool is expected (phenotype
  # missingness / per-phenotype drops), so it is reported but not flagged.
  nmax="${obs_n##*-}"
  nmax="${nmax//[^0-9]/}"
  if [ -n "$nmax" ] && [ -n "$pool_n" ]; then
    if [ "$nmax" -gt "$pool_n" ]; then
      verdict="STALE"; detail="$detail, N=$obs_n exceeds pool=$pool_n"
    elif [ "$nmax" -lt "$pool_n" ]; then
      detail="$detail, N=$obs_n < pool=$pool_n"
    else
      detail="$detail, N=$obs_n"
    fi
  elif [ "$n_checked" -eq 1 ]; then
    detail="$detail, N=NA"
  fi
  if [ -f "$tmp/$key.anom" ]; then
    acount=$(wc -l < "$tmp/$key.anom")
    detail="$detail, $acount row-size anomalies (see section 3)"
    [ "$verdict" = "KEEP" ] && verdict="ANOMALY"
  fi
  flag_detail "$key"
  if [ -n "$FLAG_TXT" ]; then
    detail="$detail, $FLAG_TXT"
    [ -n "$FLAG_VERDICT" ] && case "$verdict" in KEEP|ANOMALY) verdict="$FLAG_VERDICT" ;; esac
  fi
  if [ "$n_checked" -ne 1 ]; then
    verdict="$verdict+"
    detail="$detail, N UNVERIFIED (section 4 failed)"
  fi
  printf '%-18s %-12s %s\n' "$key" "$verdict" "$detail"
done

if [ -s "$tmp/sa.tsv" ]; then
  printf '\n%-24s %-12s %s\n' SA_FILE VERDICT DETAIL
  while IFS=$'\t' read -r b rows exp status; do
    case "$status" in
      OK)        v="KEEP";      d="rows=$rows matches template" ;;
      STALE)     v="STALE";     d="rows=$rows but SA/$b.json expects $exp - rerun" ;;
      NO_TPL)    v="UNKNOWN";   d="rows=$rows, no SA/$b.json to check against" ;;
      BAD_TPL)   v="UNKNOWN";   d="could not read npc/mpheno from SA/$b.json" ;;
      NO_OUTPUT) v="NO_OUTPUT";  d="SA/$b.json exists, results/SA/$b.csv absent" ;;
      *)         v="UNKNOWN";   d="status=$status" ;;
    esac
    flag_detail "SA/$b"
    [ -n "$FLAG_TXT" ] && d="$d; $FLAG_TXT"
    [ -n "$FLAG_VERDICT" ] && case "$v" in KEEP) v="$FLAG_VERDICT" ;; esac
    printf '%-24s %-12s %s\n' "SA/$b" "$v" "$d"
  done < "$tmp/sa.tsv"
fi

if [ "$n_checked" -ne 1 ]; then
  printf '\n%s\n' "NOTE: KEEP+ = coverage only. Staleness (N vs pool) was NOT checked."
  printf '%s\n'       "      Section 4 produced no output; fix it before trusting these verdicts."
fi

printf '\n%s\n' "Audit complete. Rerun this after each submission batch."
