#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --job-name=SA_twin
#SBATCH --output=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_twin_%j.out
#SBATCH --error=/standard/projects/coffm049/papers/functionalBrainHerit/logs/SA_twin_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

source /standard/projects/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate twinEst
# srun --nodes=1 --ntasks-per-node=1 --time=01:00:00 --pty bash -i

WORKINGDIRECTORY=/standard/projects/coffm049/papers/functionalBrainHerit/SA/twinEsts
cd $WORKINGDIRECTORY
mkdir -p $WORKINGDIRECTORY

Rscript estimate.R --vanilla
Rscript estimateWOtot.R --vanilla