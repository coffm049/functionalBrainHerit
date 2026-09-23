#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --job-name=SA_MASH
#SBATCH --output=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_MASH_%j.out
#SBATCH --error=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_MASH_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source /standard/projects/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

# MASH writes a .log into the parent of each config's 'out' dir; create them first
mkdir -p /standard/projects/coffm049/papers/functionalBrainHerit/results/SA \
         /standard/projects/coffm049/papers/functionalBrainHerit/results/SA/GCTA_wo_total \
         /standard/projects/coffm049/papers/functionalBrainHerit/results/SA/AdjHE_RE_wo_total

# Organize the data for proper use
#python /standard/projects/coffm049/papers/functionalBrainHerit/SA/joinIDs.py

# Estimates on each
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/HE.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/AdjHE_FE.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/AdjHE_Naive.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/GCTA.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/GCTA_Naive.json
#MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/GCTA_wo_total.json
MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/GCTA_wo_total.json
MASH --argfile /standard/projects/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE_wo_total.json

#rm ~/Temppheno.csv