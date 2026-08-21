#!/bin/bash -l        
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --job-name=probaConns_twin
#SBATCH --output=/users/4/coffm049/papers/functionalBrainHerit/logs/probaConns_twin_%A_%a.out
#SBATCH --error=/users/4/coffm049/papers/functionalBrainHerit/logs/probaConns_twin_%A_%a.err
#SBATCH --array=0-157


# srun -N 1  --mem=32gb -t 1:00:00 -p interactive --pty bash 
# 0-117

WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/FCs/probaConns
cd $WORKINGDIRECTORY
module load R/4.4.0-openblas-rocky8

Rscript twinEsts.R --vanilla ${SLURM_ARRAY_TASK_ID}
