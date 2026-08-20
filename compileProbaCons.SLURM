#!/bin/bash -l         
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
find /home/feczk001/shared/projects/ABCD/threshold_template_maps/data/group1_10minonly_FD0p2 -name "*10*.pconn.nii" > proba.files
find /home/feczk001/shared/projects/ABCD/threshold_template_maps/data/group2_10minonly_FD0p2 -name "*10*.pconn.nii" >> proba.files


cd ~/ABCD/Workflow/02_Phenotypes/FCsTopo 
source /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/venv/bin/activate

python compileProbaCons.py
#python filterProbaconns.py
