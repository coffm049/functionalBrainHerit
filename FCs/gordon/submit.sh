#!/bin/bash
# Generate per-chunk MASH configs and submit a SLURM array to estimate them.
# usage: bash submit.sh <Method> <FE|RE> [--chunk N] [--time HH:MM:SS] [--mem Ng] [--partition NAME]
set -euo pipefail

METHOD=$1
KIND=$2
WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/gordon
cd "$WORKINGDIRECTORY"

# MASH writes logs into the parent of 'out'; create them first
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/results/FCs/gordon \
         /users/4/coffm049/papers/functionalBrainHerit/results/FCs/gordon/pconn.0

# gordon (Gordon parcellation): 61776 FC phenotypes (o0..o61775)
CHUNK=208
TOTAL=61776
PREFIX=pconns
TIME=3:00:00
MEM=32g
PARTITION=

case "$KIND" in
  FE) TEMPLATE=feExample2.json ;;
  RE) TEMPLATE=reExample2.json ;;
  *) echo "KIND must be FE or RE"; exit 1 ;;
esac

shift 2
while [ $# -gt 0 ]; do
  case "$1" in
    --chunk) CHUNK=$2; shift 2 ;;
    --time) TIME=$2; shift 2 ;;
    --mem) MEM=$2; shift 2 ;;
    --partition) PARTITION=$2; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

NJOBS=$(python ../make_configs.py --template "$TEMPLATE" --method "$METHOD" --kind "$KIND" \
  --chunk "$CHUNK" --total "$TOTAL" --prefix "$PREFIX" --all)

ARGS=(--time="$TIME" --mem="$MEM" --array=0-$((NJOBS - 1)))
[ -n "$PARTITION" ] && ARGS+=(-p "$PARTITION")
ARGS+=(--export=ALL,METHOD="$METHOD",KIND="$KIND",WORKINGDIRECTORY="$WORKINGDIRECTORY")
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/logs
ARGS+=(--job-name="FCgordon_${METHOD}_${KIND}")
ARGS+=(--output="/users/4/coffm049/papers/functionalBrainHerit/logs/FCgordon_${METHOD}_${KIND}_%A_%a.out")
ARGS+=(--error="/users/4/coffm049/papers/functionalBrainHerit/logs/FCgordon_${METHOD}_${KIND}_%A_%a.err")

echo "Submitting $NJOBS jobs for $METHOD $KIND (chunk=$CHUNK phenos):"
sbatch "${ARGS[@]}" run.slurm