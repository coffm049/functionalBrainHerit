#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_GCTA
#SBATCH --output=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_GCTA_%j.out
#SBATCH --error=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_GCTA_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA GCTA (17 network_surfarea* + anthro_height_calc)
# Uses SA/GCTA.json (out: results/SA/GCTA -> results/SA/GCTA.csv)

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

mkdir -p /users/4/coffm049/papers/functionalBrainHerit/results/SA
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/logs

# Prefer a pre-normalized copy of this config if summary/prep_iids.sh made one
ARGFILE=/users/4/coffm049/papers/functionalBrainHerit/SA/GCTA.json
IIDFIX=${IIDFIX:-/scratch.global/coffm049/abcdEsts}
NORM="$IIDFIX/SA.$(basename "${ARGFILE%.json}").iidfix.json"
if [ -f "$NORM" ]; then ARGFILE=$NORM; echo "Using pre-normalized inputs: $ARGFILE"; fi
MASH --argfile "$ARGFILE"