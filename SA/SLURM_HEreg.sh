#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_HEreg
#SBATCH --output=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_HEreg_%j.out
#SBATCH --error=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_HEreg_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA HEreg with total surface area controlled
# Uses SA/HE.json (qcovar: age + totalNetworkSurface, 17 network_surfarea* phenos)

cd ~/software/MASH
source /standard/projects/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

mkdir -p /standard/projects/coffm049/papers/functionalBrainHerit/results/SA
mkdir -p /standard/projects/coffm049/papers/functionalBrainHerit/logs

MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/HE.json