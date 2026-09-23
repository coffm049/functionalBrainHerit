#!/bin/bash
# Submit twin estimates for Gordon (array job)
# Usage: bash submit_twins.sh [--chunk N] [--time HH:MM:SS] [--mem Ng] [--partition NAME] [--dry-run]

set -euo pipefail

REPO=/users/4/coffm049/papers/functionalBrainHerit
CHUNK=528
TOTAL=61776
TIME=6:00:00
MEM=32g
PARTITION=""
DRY_RUN=false

while [ $# -gt 0 ]; do
  case "$1" in
    --chunk) CHUNK=$2; shift 2 ;;
    --time) TIME=$2; shift 2 ;;
    --mem) MEM=$2; shift 2 ;;
    --partition) PARTITION=$2; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

NJOBS=$(( (TOTAL + CHUNK - 1) / CHUNK ))
echo "Submitting Gordon twin estimates: $NJOBS jobs (chunk=$CHUNK, total=$TOTAL)"

cd "$REPO/FCs/gordon"
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/logs

ARGS=(--time="$TIME" --mem="$MEM" --array=0-$((NJOBS - 1)))
[ -n "$PARTITION" ] && ARGS+=(-p "$PARTITION")
ARGS+=(--export=ALL)
ARGS+=(--job-name="gordon_twin")
ARGS+=(--output="/users/4/coffm049/papers/functionalBrainHerit/logs/gordon_twin_%A_%a.out")
ARGS+=(--error="/users/4/coffm049/papers/functionalBrainHerit/logs/gordon_twin_%A_%a.err")

if [ "$DRY_RUN" = true ]; then
  echo "Would run: sbatch ${ARGS[@]} twinEsts.slurm"
else
  sbatch "${ARGS[@]}" twinEsts.slurm
fi