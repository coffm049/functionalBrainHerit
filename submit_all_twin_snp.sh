#!/bin/bash
# Master script to submit all SNP heritability and twin estimates
# Usage: bash submit_all_snp.sh [--methods METHOD_LIST] [--partition NAME] [--dry-run]

set -euo pipefail

REPO=/users/4/coffm049/papers/functionalBrainHerit
METHODS="AdjHE GCTA HEreg"
PARTITION=""
DRY_RUN=false

while [ $# -gt 0 ]; do
  case "$1" in
    --methods) METHODS=$2; shift 2 ;;
    --partition) PARTITION=$2; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

cd "$REPO"
mkdir -p logs

echo "=== Submitting SNP estimates ==="
for METHOD in $METHODS; do
  echo "--- Submitting $METHOD for Gordon ---"
  if [ "$DRY_RUN" = true ]; then
    echo "  Would run: bash FCs/gordon/submit.sh $METHOD FE"
    echo "  Would run: bash FCs/gordon/submit.sh $METHOD RE"
  else
    bash FCs/gordon/submit.sh "$METHOD" FE
    bash FCs/gordon/submit.sh "$METHOD" RE
  fi

  echo "--- Submitting $METHOD for ProbaConns ---"
  if [ "$DRY_RUN" = true ]; then
    echo "  Would run: bash FCs/probaConns/submit.sh $METHOD FE"
    echo "  Would run: bash FCs/probaConns/submit.sh $METHOD RE"
  else
    bash FCs/probaConns/submit.sh "$METHOD" FE
    bash FCs/probaConns/submit.sh "$METHOD" RE
  fi
done

echo "--- Submitting SA estimates ---"
for M in $METHODS; do
  case $M in
    AdjHE) echo "Submitting SA AdjHE-RE..."; [ "$DRY_RUN" = true ] && echo "  Would run: sbatch SA/SLURM_AdjHE_RE_with_total.sh" || sbatch SA/SLURM_AdjHE_RE_with_total.sh ;;
    GCTA) echo "Submitting SA GCTA..."; [ "$DRY_RUN" = true ] && echo "  Would run: sbatch SA/GCTA_height.SLURM" || sbatch SA/GCTA_height.SLURM ;;
    HEreg) echo "Submitting SA HEreg..."; [ "$DRY_RUN" = true ] && echo "  Would run: sbatch SA/SLURM_HEreg.sh" || sbatch SA/SLURM_HEreg.sh ;;
  esac
done

echo ""
echo "=== Submitting Twin estimates (conda env: twinEst) ==="
if [ "$DRY_RUN" = true ]; then
  echo "  Would run: sbatch FCs/gordon/twinEsts.SLURM"
  echo "  Would run: sbatch FCs/probaConns/twinEsts.SLURM"
  echo "  Would run: sbatch SA/twinEsts/estimate.SLURM"
else
  echo "Submitting Gordon twins..."
  sbatch FCs/gordon/twinEsts.SLURM
  echo "Submitting ProbaConns twins..."
  sbatch FCs/probaConns/twinEsts.SLURM
  echo "Submitting SA twins..."
  sbatch SA/twinEsts/estimate.SLURM
fi

echo ""
echo "Done. Check job queue with: squeue -u \$USER"