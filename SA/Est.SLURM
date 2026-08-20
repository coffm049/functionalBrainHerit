#!/bin/bash -l        
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu

source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate
cd ~/software/MASH

# Organize the data for proper use 
#python /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/joinIDs.py

# Estimates on each 
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/HE.json
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/AdjHE_FE.json
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/AdjHE_RE.json
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/AdjHE_Naive.json
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/GCTA.json
#python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/GCTA_Naive.json
python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/GCTA_wo_total.json
python -m Estimate --argfile /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/AdjHE_RE_wo_total.json

#rm ~/Temppheno.csv
