#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

MASH --argfile "$WORKINGDIRECTORY/temp/${METHOD}.${KIND}.${SLURM_ARRAY_TASK_ID}.json"