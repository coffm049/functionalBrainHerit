#!/bin/bash
# prep_iids.sh - copy MASH's input files into a temp dir with subject IDs
# canonicalized, leaving MASH itself untouched.
#
# Canonical form (identical to compileProbaCons.py:14-19, twinEsts.R:31 and
# SA/joinIDs.py:7): strip any leading "sub-" / "NDAR_INV" / "NDARINV", then
# prefix "sub-".
#     NDAR_INV00CY2MDM -> sub-00CY2MDM
#     sub-NDARINV00CY  -> sub-00CY
#     sub-00CY         -> sub-00CY      (already canonical, file not copied)
#
# MASH joins every input on exact FID+IID strings (load_data.py:86 and
# Estimation.py:45) and never normalizes prefixes, so one file left in the
# wrong namespace silently yields an empty join.
#
# Usage:
#   bash summary/prep_iids.sh                      # -> /scratch.global/coffm049/abcdEsts
#   OUT=/scratch/x bash summary/prep_iids.sh        # choose destination
#   DRY_RUN=1 bash summary/prep_iids.sh             # report only, write nothing
#   bash summary/prep_iids.sh FCs/gordon/reExample2.json
#     -> also writes <OUT>/<dataset>.<template>.iidfix.json with paths swapped
#        (dataset = basename of the template's directory, so gordon and
#        probaConns templates that share a filename cannot collide)
#
# The destination is created (mkdir -p) before anything is written into it.
#
# Sources are never modified. Every rewritten file is verified before being
# accepted: a copy is discarded unless (a) only the ID columns changed,
# (b) each ID changed by exactly canon(), so row order is proven intact, and
# (c) no non-canonical ID remains.

set -euo pipefail
shopt -s nullglob

# All three roots are overridable so the script can be smoke-tested against a
# fixture tree without touching the real data.
ROOT=${ROOT:-/users/4/coffm049/papers/functionalBrainHerit}
BASE=${BASE:-/projects/standard/rando149/coffm049}
OUT=${OUT:-/scratch.global/coffm049/abcdEsts}
DRY_RUN=${DRY_RUN:-0}

cd "$ROOT"
# create the destination first: everything below writes report.txt/paths.tsv
# into it unconditionally, so a missing directory must abort here, not later
if ! mkdir -p "$OUT" 2>/dev/null || [ ! -w "$OUT" ]; then
  echo "ERROR: cannot create or write $OUT" >&2
  exit 1
fi
mkdir -p "$OUT/grm" || { echo "ERROR: cannot create $OUT/grm" >&2; exit 1; }

# drop artifacts from any previous run so nothing stale survives into paths.tsv
rm -f "$OUT"/ids.tsv "$OUT"/pc_eigenvec.tsv "$OUT"/covar.tsv \
      "$OUT"/pheno_*.parquet "$OUT"/*.iidfix.json "$OUT"/.canon.sed
rm -rf "$OUT/grm"
mkdir -p "$OUT/grm"

REPORT="$OUT/report.txt"
PATHS="$OUT/paths.tsv"
: > "$REPORT"
: > "$PATHS"

say() { printf '%s\n' "$*" | tee -a "$REPORT"; }
hr()  { say "----------------------------------------------------------------"; }

say "prep_iids.sh  out=$OUT  dry_run=$DRY_RUN"
say "canonical rule: strip sub-/NDAR_INV/NDARINV then prefix 'sub-'"
hr

# GNU sed BRE. Anchored to a field boundary (start of line or whitespace) so
# numeric FIDs and unrelated fields are never touched. Three passes handle the
# chained form sub-NDARINVxxx in one go.
SEDRULES="$OUT/.canon.sed"
cat > "$SEDRULES" <<'SED'
s/\(^\|[[:space:]]\)NDAR_INV/\1sub-/g
s/\(^\|[[:space:]]\)NDARINV/\1sub-/g
s/\(^\|[[:space:]]\)sub-\(NDAR_INV\|NDARINV\)/\1sub-/g
SED

# canon()/expect() mirrored as awk so the verifier can recompute the expected
# result. canon() always prefixes "sub-"; expect() only applies it when the
# source value actually carries a namespace prefix, so a numeric FID (or a
# header cell) must come through byte-identical - which is exactly what sed
# does, since its rules are anchored to those prefixes.
AWK_CANON='
function canon(s,   i, p, cr) {
  cr = ""
  if (s ~ /\r$/) { cr = "\r"; sub(/\r$/, "", s) }
  p[1] = "sub-"; p[2] = "NDAR_INV"; p[3] = "NDARINV"
  for (i = 1; i <= 3; i++) if (index(s, p[i]) == 1) s = substr(s, length(p[i]) + 1)
  return "sub-" s cr
}
function expect(s) {
  if (index(s, "sub-") == 1 || index(s, "NDAR_INV") == 1 || index(s, "NDARINV") == 1) return canon(s)
  return s
}
'

delim_of() { if head -n 1 "$1" | grep -q $'\t'; then printf '\t'; else printf ' '; fi; }
nf_hist()   { awk -F"$1" '{ n[NF]++ } END { for (k in n) printf "%s:%s ", k, n[k] }' "$2"; }

# mask_id_cols: blank the ID columns so the remaining bytes must be identical
mask_id_cols() { awk -F"$1" -v OFS="$1" -v c="$2" '{ $1 = "<F>"; $c = "<I>"; print }' "$3"; }

# verify <src> <dst> <col> <delim>  -> prints "ok" or the failure reasons
verify() {
  local src=$1 dst=$2 col=$3 delim=$4
  local why=""
  [ "$(wc -l < "$src")" = "$(wc -l < "$dst")" ] || why="$why line-count"
  [ "$(nf_hist "$delim" "$src")" = "$(nf_hist "$delim" "$dst")" ] || why="$why field-counts"
  [ "$(mask_id_cols "$delim" "$col" "$src" | md5sum | cut -d' ' -f1)" = \
    "$(mask_id_cols "$delim" "$col" "$dst" | md5sum | cut -d' ' -f1)" ] || why="$why non-id-columns"

  # per-line: every ID must equal expect(source ID) -> proves order preserved
  awk -F"$delim" -v c="$col" "$AWK_CANON"'
    FNR == NR { f[FNR] = $1; i[FNR] = $c; last = FNR; next }
    FNR > last { exit }
    {
      # compare with trailing CR ignored: some seds translate CRLF to LF on
      # read, and either line ending parses fine in pandas
      a = expect(f[FNR]); b = $1; sub(/\r$/, "", a); sub(/\r$/, "", b); if (a != b) bf++
      a = expect(i[FNR]); b = $c; sub(/\r$/, "", a); sub(/\r$/, "", b); if (a != b) bi++
    }
    END { if (bf + bi > 0) print "id-not-exactly-canon(" bf "/" bi ")" }
  ' "$src" "$dst"

  # nothing still carrying a foreign prefix (header cells and numeric ids are
  # not foreign, so they are correctly not counted here)
  awk -F"$delim" -v c="$col" '{ s = $c; sub(/\r$/, "", s);
                    if (s ~ /^(sub-)?(NDAR_INV|NDARINV)/) n++ }
                  END { if (n + 0 > 0) print "non-canonical-left(" n ")" }' "$dst" |
    tr '\n' ' '
}

printf 'name\tstatus\tsrc\tdst\n' > "$PATHS"

# name|source path|iid col (2, or HEADER)|kind (text|grm)
MANIFEST=(
  "ids|$BASE/filtered_ids.tsv|2|text"
  "grm_no_rels|$BASE/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/no_rels/no_rels.grm.id|2|grm"
  "grm_full|$BASE/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full.grm.id|2|grm"
  "pc_eigenvec|$BASE/ABCD/Results/01_Gene_QC/filters/filter1/Eigens/no_rels.eigenvec|2|text"
  "covar|$BASE/ABCD/Results/02_Phenotypes/Covars_with_totalNetworkSurface.tsv|HEADER|text"
)

for entry in "${MANIFEST[@]}"; do
  IFS='|' read -r name src col kind <<< "$entry"
  if [ ! -f "$src" ]; then
    say "SKIP   $name  missing: $src"
    printf '%s\t%s\t%s\t%s\n' "$name" "MISSING" "$src" "-" >> "$PATHS"
    continue
  fi

  delim=$(delim_of "$src")

  if [ "$col" = "HEADER" ]; then
    col=$(awk -F"$delim" -v OFS="$delim" '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          h = $i; gsub(/["\r]/, "", h); gsub(/^#/, "", h)
          if (tolower(h) ~ /^(iid|subject_id|src_subject_id|subid)$/) { print i; exit }
        }
        exit 1
      }' "$src") || {
      say "SKIP   $name  no IID-like column; header = $(head -n 1 "$src")"
      printf '%s\t%s\t%s\t%s\n' "$name" "NO_IID_COL" "$src" "-" >> "$PATHS"
      continue
    }
  fi

  # destination: GRM outputs keep their .grm.id basename so prefix works
  if [ "$kind" = "grm" ]; then
    dst="$OUT/grm/$(basename "$src")"
    base_noext=${src%.grm.id}
  else
    dst="$OUT/$name.tsv"
  fi

  before_n=$(wc -l < "$src")
  before_foreign=$(awk -F"$delim" -v c="$col" '{ s = $c; sub(/\r$/, "", s);
                       if (s ~ /^(sub-)?(NDAR_INV|NDARINV)/) n++ }
                       END { print n + 0 }' "$src")

  if [ "$DRY_RUN" = 1 ]; then
    say "DRY    $name  col=$col lines=$before_n non-canonical=$before_foreign"
    printf '%s\t%s\t%s\t%s\n' "$name" "DRY_RUN" "$src" "-" >> "$PATHS"
    continue
  fi

  if [ "$before_foreign" -eq 0 ]; then
    # nothing to write - and remove any copy left by an earlier run so a stale
    # file can never be picked up as this file's destination
    rm -f "$dst"
    say "OK     $name  already canonical ($before_n ids) - using ORIGINAL"
    printf '%s\t%s\t%s\t%s\n' "$name" "UNCHANGED" "$src" "$src" >> "$PATHS"
    [ "$kind" = "grm" ] && printf 'grm_prefix_%s\tUNCHANGED\t%s\t%s\n' "$name" "$base_noext" "$base_noext" >> "$PATHS"
    continue
  fi

  sed -f "$SEDRULES" "$src" > "$dst"

  errs=$(verify "$src" "$dst" "$col" "$delim")
  if [ -n "$(printf '%s' "$errs" | tr -d ' ')" ]; then
    say "FAIL   $name  ${errs}(source left untouched, destination removed)"
    rm -f "$dst"
    printf '%s\t%s\t%s\t%s\n' "$name" "FAIL" "$src" "-" >> "$PATHS"
    continue
  fi

  if [ "$kind" = "grm" ]; then
    # MASH reads <prefix>.grm.id AND <prefix>.grm.bin - symlink the binary
    [ -f "$base_noext.grm.bin" ] && ln -sf "$base_noext.grm.bin" "$OUT/grm/$(basename "$base_noext").grm.bin"
    say "OK     $name  rewrote $before_foreign/$before_n ids -> $dst (+ symlinked .grm.bin)"
    printf '%s\t%s\t%s\t%s\n' "$name" "REWRITTEN" "$src" "$dst" >> "$PATHS"
    printf 'grm_prefix_%s\tREWRITTEN\t%s\t%s\n' "$name" "$base_noext" "$OUT/grm/$(basename "$base_noext")" >> "$PATHS"
  else
    say "OK     $name  rewrote $before_foreign/$before_n ids -> $dst"
    printf '%s\t%s\t%s\t%s\n' "$name" "REWRITTEN" "$src" "$dst" >> "$PATHS"
  fi
done

hr
#--- parquet inputs: awk/sed cannot touch these, use pyarrow ------------------
for p in "$BASE/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet" \
         "$BASE/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet"; do
  name="pheno_$(basename "$p" .parquet)"
  if [ ! -f "$p" ]; then
    say "SKIP   $name  missing: $p"
    printf '%s\t%s\t%s\t%s\n' "$name" "MISSING" "$p" "-" >> "$PATHS"
    continue
  fi
  if ! command -v python3 >/dev/null 2>&1 || ! python3 -c "import pyarrow" >/dev/null 2>&1; then
    say "SKIP   $name  python3/pyarrow unavailable (run under 'conda activate MASH')"
    printf '%s\t%s\t%s\t%s\n' "$name" "NO_PYARROW" "$p" "-" >> "$PATHS"
    continue
  fi
  res=$(python3 - "$p" "$OUT/$name.parquet" <<'PY' 2>&1 || true
import sys, os, hashlib
import pyarrow as pa, pyarrow.parquet as pq

def canon(v):
    s = str(v)
    for p in ("sub-", "NDAR_INV", "NDARINV"):
        if s.startswith(p):
            s = s[len(p):]
    return "sub-" + s

src, dst = sys.argv[1], sys.argv[2]
t = pq.read_table(src)
if "IID" not in t.column_names:
    print("ERR no IID column (columns: %s)" % ",".join(t.column_names)); sys.exit(0)
old = t.column("IID").to_pylist()
new = [canon(x) for x in old]
if len(set(new)) != len(set(old)):
    print("ERR canonicalizing merges two distinct ids"); sys.exit(0)
changed = sum(1 for a, b in zip(old, new) if a != b)
if changed == 0:
    print("UNCHANGED %d" % len(old)); sys.exit(0)

rest = t.drop(["IID"])
sig_before = hashlib.sha256(rest.to_pandas().to_csv(index=False).encode()).hexdigest()
idx = t.schema.get_field_index("IID")
t2 = t.set_column(idx, pa.field("IID", pa.string()), pa.array(new))
sig_after = hashlib.sha256(t2.drop(["IID"]).to_pandas().to_csv(index=False).encode()).hexdigest()
if sig_before != sig_after:
    print("ERR non-IID columns changed"); sys.exit(0)
if t2.num_rows != t.num_rows:
    print("ERR row count changed"); sys.exit(0)
pq.write_table(t2, dst)
print("CHANGED %d" % changed)
PY
)
  case "$res" in
    CHANGED\ *)
      say "OK     $name  rewrote ${res#CHANGED } ids -> $OUT/$name.parquet"
      printf '%s\t%s\t%s\t%s\n' "$name" "REWRITTEN" "$p" "$OUT/$name.parquet" >> "$PATHS" ;;
    UNCHANGED*)
      say "OK     $name  already canonical (${res#UNCHANGED }) - using ORIGINAL"
      printf '%s\t%s\t%s\t%s\n' "$name" "UNCHANGED" "$p" "$p" >> "$PATHS" ;;
    *)
      say "FAIL   $name  ${res:-no output}"
      printf '%s\t%s\t%s\t%s\n' "$name" "FAIL" "$p" "-" >> "$PATHS" ;;
  esac
done

#--- optional: emit configs whose paths point at the normalized copies --------
# Before touching any config, decide whether this run is safe to submit from.
# A copy that failed verification is always fatal. A parquet we could not copy
# is fatal only when OTHER inputs were rewritten: that combination would leave
# MASH joining normalized ids against un-normalized phenotypes, which is
# exactly the silent empty join this script exists to prevent. (A plain
# MISSING file fails loudly downstream, so it stays a warning.)
hard=$(awk -F'\t' 'NR>1 && ($2=="FAIL" || $2=="NO_IID_COL") {c++} END{print c+0}' "$PATHS")
pyskip=$(awk -F'\t' 'NR>1 && $2=="NO_PYARROW" {c++} END{print c+0}' "$PATHS")
rew=$(awk -F'\t' 'NR>1 && $2=="REWRITTEN" {c++} END{print c+0}' "$PATHS")
if [ "$hard" -gt 0 ] || { [ "$pyskip" -gt 0 ] && [ "$rew" -gt 0 ]; }; then
  say "ERROR: inputs failed verification (hard=$hard no_pyarrow=$pyskip rewritten=$rew)"
  say "       not writing configs - fix the inputs and rerun"
  exit 1
fi

# Pure sed over paths.tsv: every row whose destination differs from its source
# becomes one s|old|new|g expression, so no JSON parsing (and no python) is
# needed and only path strings can change.
for tpl in "$@"; do
  base=$(basename "$tpl")
  ds=$(basename "$(dirname "$tpl")")
  [ "$ds" = "." ] && ds=local
  if [ ! -f "$tpl" ]; then say "SKIP   config $ds/$base  missing: $tpl"; continue; fi
  if [ "$DRY_RUN" = 1 ]; then say "DRY    config $ds/$base"; continue; fi
  outjson="$OUT/$ds.${base%.json}.iidfix.json"
  args=(); n=0
  while IFS=$'\t' read -r nm st src dst; do
    [ "$nm" = "name" ] && continue                    # header row
    [ -z "$dst" ] && continue
    [ "$dst" = "-" ] && continue                      # file failed or missing
    [ "$dst" = "$src" ] && continue                   # unchanged, nothing to do
    args+=(-e "s|$src|$dst|g")
    n=$((n + 1))
  done < "$PATHS"
  if [ "$n" -eq 0 ]; then
    cp "$tpl" "$outjson"
    say "OK     config $ds/$base -> $outjson (no path changes needed)"
    continue
  fi
  sed "${args[@]}" "$tpl" > "$outjson"
  if cmp -s "$tpl" "$outjson"; then
    say "WARN   config $ds/$base unchanged - none of the rewritten paths appear in it"
    rm -f "$outjson"
  else
    say "OK     config $ds/$base -> $outjson ($n path substitution(s))"
  fi
done

hr
say "Summary:"
awk -F'\t' 'NR > 1 && $1 !~ /^grm_prefix/ { s[$2]++ }
            END { for (k in s) printf "  %-12s %d\n", k, s[k] }' "$PATHS" | tee -a "$REPORT"
say ""
say "Paths MASH should read:"
awk -F'\t' 'NR > 1 && $1 !~ /^grm_prefix/ && $4 != "-" { printf "  %-14s %s\n", $1, $4 }' "$PATHS" | tee -a "$REPORT"
awk -F'\t' '$1 ~ /^grm_prefix/ && $4 != $3 { printf "  %-14s %s\n", $1, $4 }' "$PATHS" | tee -a "$REPORT"
say ""
say "Smoke-test one chunk with a rewritten config:"
say "  conda activate MASH && MASH --argfile $OUT/<dataset>.<template>.iidfix.json"
say "Full report: $REPORT"
