#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --time=12:00:00
#SBATCH --job-name=FCtop5_GCTA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu
# GCTA (via MASH) on the 5 boundary h2=1.0 Gordon FCs: full sample + White subset.
# The White subset is handled inside MASH via covar_filter (race_ethnicity).
# Run AFTER summary/top5_fc_followup.py (for the distribution plots/stats).
# Usage: sbatch FCs/gordon/top5.slurm   (from the functionalBrainHerit root on HPC)
set -euo pipefail
cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/gordon
echo "=== full sample ==="
MASH --argfile "$WORKINGDIRECTORY/top5_gcta_full.json"
echo "=== White subset ==="
MASH --argfile "$WORKINGDIRECTORY/top5_gcta_white.json"
echo "Done"
