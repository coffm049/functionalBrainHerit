#!/bin/bash -l        
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu


module load python3
cd ~/tools/MASH

# Estimates on each 
python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/GCTA_height.json
