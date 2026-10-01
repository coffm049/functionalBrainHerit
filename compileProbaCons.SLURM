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

# Only group1 and group2 are searched here, while compilePcons.SLURM:16
# searches the whole abcd-hcp-pipeline tree. That asymmetry is why probaConns
# has 5,552 subjects and gordon has 10,123, and it is what caps the shared
# filtered_ids pool. Run summary/01_find_probaconns.sh to see whether other
# groups hold *10*.pconn.nii files that were never listed.
find /projects/standard/feczk001/shared/projects/ABCD/threshold_template_maps/data/group1_10minonly_FD0p2 -name "*10*.pconn.nii" > "$FCSTOPO/proba.files"
find /projects/standard/feczk001/shared/projects/ABCD/threshold_template_maps/data/group2_10minonly_FD0p2 -name "*10*.pconn.nii" >> "$FCSTOPO/proba.files"

python compileProbaCons.py
#python filterProbaconns.py