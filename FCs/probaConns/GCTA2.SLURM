#!/bin/bash -l        
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-315
#SBATCH --account=saonli

METHOD=GCTA
# SLURM_ARRAY_TASK_ID=1
# 0-100
WORKINGDIRECTORY=/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns
cd $WORKINGDIRECTORY
source ~/software/MASH/.venv/bin/activate

# Create temp directory if it doesn't exist
mkdir -p $WORKINGDIRECTORY/temp

python makeFE.py -i $SLURM_ARRAY_TASK_ID -m $METHOD

cd ~/software/MASH

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/$METHOD.FE.${SLURM_ARRAY_TASK_ID}.json
