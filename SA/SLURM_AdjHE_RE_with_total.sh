#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_AdjHE_RE_total
#SBATCH --output=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_RE_total_%j.out
#SBATCH --error=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_AdjHE_RE_total_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA AdjHE_RE with total surface area controlled (qcovar includes totalNetworkSurface)
# Uses SA/AdjHE_RE.json  (qcovar: age + totalNetworkSurface, 17 network_surfarea* phenos)
# Compare to wo_total: SA/AdjHE_RE_wo_total.json (qcovar: age only)
# Run: sbatch SA/SLURM_AdjHE_RE_with_total.sh

cd ~/software/MASH
source /standard/projects/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

mkdir -p /standard/projects/coffm049/papers/functionalBrainHerit/results/SA
mkdir -p /standard/projects/coffm049/papers/functionalBrainHerit/logs

MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE.json
