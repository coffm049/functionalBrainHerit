#!/bin/bash -l        
#SBATCH --time=3:30:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-9


WORKINGDIRECTORY=/home/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns
cd $WORKINGDIRECTORY
source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate
cd ~/tools/MASH

# run ./makejsons.sh referencing the job array id to make json file
$WORKINGDIRECTORY/makeFE.sh $SLURM_ARRAY_TASK_ID AdjHE $WORKINGDIRECTORY

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/AdjHE.FE_0.${SLURM_ARRAY_TASK_ID}.json
