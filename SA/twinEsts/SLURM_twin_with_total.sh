#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --job-name=SA_twin_total
#SBATCH --output=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_twin_total_%j.out
#SBATCH --error=/users/4/coffm049/papers/functionalBrainHerit/logs/SA_twin_total_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu

# SA twin ACE with total surface area controlled
# Model: value ~ site_id_l + age + female + household.income + high.educ + totalNetworkSurface
# Script: SA/twinEsts/estimate.R  -> results/SA/twinEsts/herit_w_total.Rds
# Compare to wo_total: SA/twinEsts/estimateWOtot.R (no totalNetworkSurface) -> herit_wo_total.Rds
# Run: sbatch SA/twinEsts/SLURM_twin_with_total.sh

module load R/4.4.0-openblas-rocky8

WORKINGDIRECTORY=/users/4/coffm049/papers/functionalBrainHerit/SA/twinEsts
cd $WORKINGDIRECTORY

Rscript estimate.R --vanilla
