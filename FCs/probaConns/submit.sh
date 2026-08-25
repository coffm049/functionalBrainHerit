#!/bin/bash
# Generate per-chunk MASH configs and submit a SLURM array to estimate them.
# usage: bash submit.sh <Method> <FE|RE|GCTA|HEreg> [--chunk N] [--time HH:MM:SS] [--mem Ng] [--partition NAME] [--total N]
set -euo pipefail

METHOD=$1
KIND=$2
WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/probaConns
cd "$WORKINGDIRECTORY"

# MASH writes logs into the parent of 'out'; create them first
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/results/FCs/probaConns \
         /users/4/coffm049/papers/functionalBrainHerit/results/FCs/probaConns/probaConns

# probabilistic connectomes: 3160 FC phenotypes (o0..o3159)
CHUNK=10
TOTAL=3160
PREFIX=probaConns
TIME=3:00:00
MEM=64g
PARTITION=

case "$KIND" in
  FE) CHUNK=10; TEMPLATE=feExample2.json ;;
  RE) CHUNK=40; TEMPLATE=reExample2.json ;;
  GCTA) TEMPLATE=gctaExample.json ;;
  HEreg) TEMPLATE=heExample.json ;;
  *) echo "KIND must be FE, RE, GCTA or HEreg"; exit 1 ;;
esac

shift 2
while [ $# -gt 0 ]; do
  case "$1" in
    --chunk) CHUNK=$2; shift 2 ;;
    --time) TIME=$2; shift 2 ;;
    --mem) MEM=$2; shift 2 ;;
    --partition) PARTITION=$2; shift 2 ;;
    --total) TOTAL=$2; shift 2 ;;
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
ARGS+=(--job-name="FCproba_${METHOD}_${KIND}")
ARGS+=(--output="/users/4/coffm049/papers/functionalBrainHerit/logs/FCproba_${METHOD}_${KIND}_%A_%a.out")
ARGS+=(--error="/users/4/coffm049/papers/functionalBrainHerit/logs/FCproba_${METHOD}_${KIND}_%A_%a.err")

echo "Submitting $NJOBS jobs for $METHOD $KIND (chunk=$CHUNK phenos):"
sbatch "${ARGS[@]}" run.slurm