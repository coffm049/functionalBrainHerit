#!/bin/bash -l         
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu

REPO=/users/4/coffm049/papers/functionalBrainHerit
FCSTOPO=/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo
SRC=/spaces/regulated/abcd-ref-data/abcc/derivatives/abcd-hcp-pipeline_v0.1.4

cd "$REPO"
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

find "$SRC" -name "*atlas-Gordon2014FreeSurferSubcortical*thresh-fd0p2mm_censor-10min*conndata-network_connectivity.pconn.nii" > "$FCSTOPO/pconn.files"

python compilePcons.py
#python filterPconns.py