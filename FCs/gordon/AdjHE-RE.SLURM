#!/bin/bash -l        
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-100


WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/gordon
cd $WORKINGDIRECTORY
conda activate Herit_ests
cd ~/tools/MASH2

# run ./makejsons.sh referencing the job array id to make json file
$WORKINGDIRECTORY/makeRE.sh $SLURM_ARRAY_TASK_ID AdjHE $WORKINGDIRECTORY

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/AdjHE.RE.${SLURM_ARRAY_TASK_ID}.json
