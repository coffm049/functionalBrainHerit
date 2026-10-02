#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_AdjHE_FE_total
#SBATCH --output=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_FE_total_%j.out
#SBATCH --error=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_FE_total_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA AdjHE-FE with total surface area controlled (qcovar includes totalNetworkSurface,
# abcd_site a fixed discrete covariate since there is no random group in the FE spec).
# Uses SA/AdjHE_FE.json  (qcovar: age + totalNetworkSurface, 17 network_surfarea* phenos)
# Compare to wo_total: SA/AdjHE_FE_wo_total.json (qcovar: age only)
# Run: sbatch SA/SLURM_AdjHE_FE_with_total.sh

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

mkdir -p /users/4/coffm049/papers/functionalBrainHerit/results/SA
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/logs

# Prefer a pre-normalized copy of this config if summary/prep_iids.sh made one
ARGFILE=/users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_FE.json
IIDFIX=${IIDFIX:-/scratch.global/coffm049/abcdEsts}
NORM="$IIDFIX/SA.$(basename "${ARGFILE%.json}").iidfix.json"
if [ -f "$NORM" ]; then ARGFILE=$NORM; echo "Using pre-normalized inputs: $ARGFILE"; fi
MASH --argfile "$ARGFILE"
