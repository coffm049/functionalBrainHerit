#!/bin/bash -l        
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-9

WORKINGDIRECTORY=/home/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns
source .venv/bin/activate

cd $WORKINGDIRECTORY

# run ./makejsons.sh referencing the job array id to make json file
./makeFE.sh $SLURM_ARRAY_TASK_ID GCTA $WORKINGDIRECTORY


cd ~/software/MASH

# Estimates on each 
python -m Estimate --argfile $WORKINGDIRECTORY/temp/GCTA.FE.${SLURM_ARRAY_TASK_ID}.json
"/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/data/probaConns/part.1.parquet"
