#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_AdjHE_RE_wo_total
#SBATCH --output=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_RE_wo_total_%j.out
#SBATCH --error=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_RE_wo_total_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA AdjHE-RE without totalNetworkSurface control (qcovar: age)
# Uses SA/AdjHE_RE_wo_total.json (out: results/SA/AdjHE_RE_wo_total -> results/SA/AdjHE_RE_wo_total.csv)

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

mkdir -p /users/4/coffm049/papers/functionalBrainHerit/results/SA
mkdir -p /users/4/coffm049/papers/functionalBrainHerit/logs

MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE_wo_total.json