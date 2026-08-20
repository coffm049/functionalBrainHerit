#!/bin/bash
# Generate per-chunk MASH configs and submit a SLURM array to estimate them.
# usage: bash submit.sh <Method> <FE|RE> [--chunk N] [--time HH:MM:SS] [--mem Ng] [--partition NAME]
set -euo pipefail

METHOD=$1
KIND=$2
WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/probaConns
cd "$WORKINGDIRECTORY"

# probabilistic connectomes: 3160 FC phenotypes (o0..o3159)
CHUNK=10
TOTAL=3160
PREFIX=probaConns
TIME=3:00:00
MEM=32g
PARTITION=

case "$KIND" in
  FE) CHUNK=10; TEMPLATE=feExample2.json ;;
  RE) CHUNK=40; TEMPLATE=reExample2.json ;;
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

echo "Submitting $NJOBS jobs for $METHOD $KIND (chunk=$CHUNK phenos):"
sbatch "${ARGS[@]}" run.slurm