#!/bin/bash -l         
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu


cd ~/ABCD/Workflow/02_Phenotypes/FCsTopo 
source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate

python compilePcons.py
#python filterPconns.py