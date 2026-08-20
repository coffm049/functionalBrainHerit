#!/bin/bash -l
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu


WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/gordon
cd $WORKINGDIRECTORY
cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

# Estimates on each
MASH --argfile $WORKINGDIRECTORY/AdjHE-17-RE.json