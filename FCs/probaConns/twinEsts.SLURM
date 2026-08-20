#!/bin/bash -l        
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=0-158


# srun -N 1  --mem=32gb -t 1:00:00 -p interactive --pty bash 
# 0-117

WORKINGDIRECTORY=~/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns
cd $WORKINGDIRECTORY
module load R/4.4.0-openblas-rocky8

Rscript twinEsts.R --vanilla ${SLURM_ARRAY_TASK_ID}
