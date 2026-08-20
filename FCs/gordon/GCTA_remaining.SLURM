#!/bin/bash -l        
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH -p msismall
#SBATCH --array=39,59,78,98


WORKINGDIRECTORY=/home/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs
cd $WORKINGDIRECTORY
conda activate Herit_ests
cd ~/tools/MASH2

# run ./makejsons.sh referencing the job array id to make json file
/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/makeFE.sh $SLURM_ARRAY_TASK_ID GCTA $WORKINGDIRECTORY

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/GCTA.FE.${SLURM_ARRAY_TASK_ID}.json

