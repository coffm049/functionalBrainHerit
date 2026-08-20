#!/bin/bash -l        
#SBATCH --time=3:30:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-296

module load python/3.9.3_anaconda2021.11_mamba

METHOD=AdjHE
# SLURM_ARRAY_TASK_ID=1
# 0-100

WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/gordon
cd $WORKINGDIRECTORY
source ~/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate 
# Create temp directory if it doesn't exist
mkdir -p $WORKINGDIRECTORY/temp

python makeRE.py -i $SLURM_ARRAY_TASK_ID -m $METHOD

cd ~/tools/MASH

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/$METHOD.RE.${SLURM_ARRAY_TASK_ID}.json
