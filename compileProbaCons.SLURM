#!/bin/bash -l         
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu

REPO=/users/4/coffm049/papers/functionalBrainHerit
FCSTOPO=/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo

cd "$REPO"
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

find /projects/standard/feczk001/shared/projects/ABCD/threshold_template_maps/data/group1_10minonly_FD0p2 -name "*10*.pconn.nii" > "$FCSTOPO/proba.files"
find /projects/standard/feczk001/shared/projects/ABCD/threshold_template_maps/data/group2_10minonly_FD0p2 -name "*10*.pconn.nii" >> "$FCSTOPO/proba.files"

python compileProbaCons.py
#python filterProbaconns.py