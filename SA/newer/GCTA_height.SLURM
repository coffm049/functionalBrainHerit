#!/bin/bash -l
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

# Estimates on each
MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/newer/GCTA_height.json