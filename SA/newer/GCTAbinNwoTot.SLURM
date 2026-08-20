#!/bin/bash -l        
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu

cd ~/software/MASH
source .venv/bin/activate

python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/newer/GCTAbigNWOtot.json
#rm ~/Temppheno.csv
