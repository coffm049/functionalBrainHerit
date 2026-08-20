#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

# Organize the data for proper use
#python /users/4/coffm049/papers/functionalBrainHerit/SA/joinIDs.py

# Estimates on each
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/HE.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_FE.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_Naive.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/GCTA.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/GCTA_Naive.json
#MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/GCTA_wo_total.json
MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/GCTA_wo_total.json
MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/AdjHE_RE_wo_total.json

#rm ~/Temppheno.csv