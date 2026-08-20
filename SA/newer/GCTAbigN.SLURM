#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

MASH --argfile /users/4/coffm049/papers/functionalBrainHerit/SA/newer/GCTAbigN.json

#rm ~/Temppheno.csv