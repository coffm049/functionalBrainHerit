import os
import re

import h5py
import nibabel as nib
import nilearn
import numpy as np
import pandas as pd
from nilearn import plotting

cifti_template = nib.load(
    "brainTemplate/example.pscalar.nii"
)  # Replace with the actual path to your template


def reformat2nii(estFile):

    # Load the scalar values from the CSV
    df = pd.read_csv(estFile)
    # if estFile has twin in its name, drop the first row
    # if "twin" in estFile:
    #    df = df.drop(0)

    # remove leading "network_surfarea" from pheno column
    df["Region"] = df["pheno"].str[16:]
    df["Value"] = df.h2
    df.loc[df["Value"] == 0, "Value"] = 1e-3
    # dropping empty parcels
    df["Region"] = range(1, 15)

    # load in the template matrix
    netlabel = pd.read_csv("brainTemplate/gordon_modules.csv")
    netlabel.columns = ["Region"]

    df = pd.merge(netlabel, df, how="left", on="Region")

    vertex_scalar_map = np.array(df.h2).reshape(1, 352)

    # Create a new CIFTI-2 image with your scalar data
    new_cifti_img = nib.Cifti2Image(
        vertex_scalar_map,
        header=cifti_template.header,
        nifti_header=cifti_template.nifti_header,
    )

    # Save to a new .dscalar.nii file
    output_file = f"{os.path.basename(estFile)}HeatMap.pscalar.nii"
    nib.save(new_cifti_img, output_file)

    return None


for estimator in os.listdir("../results/FCs100824/Topography/SA"):
    if ("Naive" not in estimator) and (".R" not in estimator):
        print(estimator)
        reformat2nii(f"../results/FCs100824/Topography/SA/{estimator}")
