#!/bin/bash
# Master script to submit all SNP heritability and twin estimates
# Usage: bash submit_all_twin_snp.sh [--methods METHOD_LIST] [--partition NAME] [--dry-run]

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
  if [ "$METHOD" = "AdjHE" ]; then
    # AdjHE has FE and RE variants
    echo "--- Submitting AdjHE-FE for Gordon ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/gordon/submit.sh AdjHE FE"
    else
      bash FCs/gordon/submit.sh AdjHE FE
    fi

    echo "--- Submitting AdjHE-RE for Gordon ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/gordon/submit.sh AdjHE RE"
    else
      bash FCs/gordon/submit.sh AdjHE RE
    fi

    echo "--- Submitting AdjHE-FE for ProbaConns ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/probaConns/submit.sh AdjHE FE"
    else
      bash FCs/probaConns/submit.sh AdjHE FE
    fi

    echo "--- Submitting AdjHE-RE for ProbaConns ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/probaConns/submit.sh AdjHE RE"
    else
      bash FCs/probaConns/submit.sh AdjHE RE
    fi
  elif [ "$METHOD" = "GCTA" ]; then
    # GCTA has no FE/RE variants - pass GCTA as KIND
    echo "--- Submitting GCTA for Gordon ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/gordon/submit.sh GCTA GCTA"
    else
      bash FCs/gordon/submit.sh GCTA GCTA
    fi

    echo "--- Submitting GCTA for ProbaConns ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/probaConns/submit.sh GCTA GCTA"
    else
      bash FCs/probaConns/submit.sh GCTA GCTA
    fi
  elif [ "$METHOD" = "HEreg" ]; then
    # HEreg has no FE/RE variants - pass HEreg as KIND
    echo "--- Submitting HEreg for Gordon ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/gordon/submit.sh HEreg HEreg"
    else
      bash FCs/gordon/submit.sh HEreg HEreg
    fi

    echo "--- Submitting HEreg for ProbaConns ---"
    if [ "$DRY_RUN" = true ]; then
      echo "  Would run: bash FCs/probaConns/submit.sh HEreg HEreg"
    else
      bash FCs/probaConns/submit.sh HEreg HEreg
    fi
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