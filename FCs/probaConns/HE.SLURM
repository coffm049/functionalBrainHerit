#!/bin/bash -l 
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-296


source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate
WORKINGDIRECTORY=/home/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns
cd $WORKINGDIRECTORY
mkdir -p $WORKINGDIRECTORY/temp
METHOD=HEreg
# run ./makejsons.sh referencing the job array id to make json file
python makeFE.py -i $SLURM_ARRAY_TASK_ID -m $METHOD

cd ~/software/MASH
# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/${METHOD}.FE.${SLURM_ARRAY_TASK_ID}.json
 
