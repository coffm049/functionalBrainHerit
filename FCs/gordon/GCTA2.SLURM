#!/bin/bash -l        
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-296

METHOD=GCTA
# SLURM_ARRAY_TASK_ID=1
# 0-100
WORKINGDIRECTORY=/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/gordon
cd $WORKINGDIRECTORY
source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate 
# Create temp directory if it doesn't exist
mkdir -p $WORKINGDIRECTORY/temp

python makeFE.py -i $SLURM_ARRAY_TASK_ID -m $METHOD

cd ~/software/MASH

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/$METHOD.FE.${SLURM_ARRAY_TASK_ID}.json
