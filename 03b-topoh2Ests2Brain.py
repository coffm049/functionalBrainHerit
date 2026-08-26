import os
import re

import h5py
import nibabel as nib
import nilearn
import numpy as np
import pandas as pd
from nilearn import plotting

# Load an existing CIFTI-2 template file to get the correct structure (like the HCP template)
# cifti_template = nib.load(
#     "brainTemplate/ones_final.dscalar.nii"
# )  # Replace with the actual path to your template
# cifti_template = nib.load(
#     "brainTemplate/single_map.dscalar.nii"
# )  # Replace with the actual path to your template
# cifti_template = nib.load(
#     "brainTemplate/Gordon.32k_fs_LR.dscalar.nii"
# )  # Replace with the actual path to your template
# dlabelFile = nib.load("brainTemplate/Gordon.networks.32k_fs_LR.dlabel.nii")
# dlabelFile = nib.load("brainTemplate/abcd0.75.dlabel.nii")
cifti_template = nib.load(
    "brainTemplate/single_map.dscalar.nii"
)  # Replace with the actual path to your template


# dlabelFile = nib.load("formattedData/networkSurfLabels.dscalar.nii")
dlabelFile = nib.load(
    "brainTemplate/abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"
)
dlabelFile = nib.load("brainTemplate/Gordon.networks.32k_fs_LR.dlabel.nii")

brainColors = dlabelFile.header.get_axis(0).get_element(0)[1]
brainColors.pop(0)
rois = {label: val for (label, val) in brainColors.values()}.keys()
rois = [re.sub(r"^\d{1,3}_._", "", item) for item in rois]

# map names from rois to a set of new names "Default" -> DMN, Auditory -> AUD
rois = [re.sub(r"Default", "DMN", item) for item in rois]
rois = [re.sub(r"Auditory", "AUD", item) for item in rois]
rois = [re.sub(r"CinguloOperc", "CO", item) for item in rois]
rois = [re.sub(r"FrontoParietal", "FP", item) for item in rois]
rois = [re.sub(r"Salience", "SAL", item) for item in rois]
rois = [re.sub(r"VentralAttn", "VAN", item) for item in rois]
rois = [re.sub(r"Visual", "VIS", item) for item in rois]
rois = [re.sub(r"DorsalAttn", "DAN", item) for item in rois]
array(
    ["", "RetrosplenialTemporal", "SMhand", "SMmouth", "", "", "Visual"],
    dtype="<U21",
)


# rois = [re.sub(r"_LABEL_\d{1,3}$", "", item) for item in rois]
# replace rois with the cumulative number of unique observations in the list
# rois = np.unique(rois, return_counts=True)
labelDictionary = pd.DataFrame(
    {
        "index": list(range(1, 15)),
        "name": [
            "SMD",
            "SML",
            "Tpole",
            "MTL",
            "PMN",
            "PON",
        ],
    }
)


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
    netlabel = pd.read_csv("brainTemplate/GRP1_template_parcel.csv")
    netlabel.columns = ["Region"]

    df = pd.merge(netlabel, df, how="left", on="Region")

    vertex_scalar_map = np.zeros(cifti_template.shape)
    for region, value in zip(df.Region, df.Value):
        vertex_scalar_map[dlabelFile.get_fdata() == region] = value

    # Create a new CIFTI-2 image with your scalar data
    new_cifti_img = nib.Cifti2Image(
        vertex_scalar_map,
        header=cifti_template.header,
        nifti_header=cifti_template.nifti_header,
    )

    # Save to a new .dscalar.nii file
    output_file = f"{os.path.basename(estFile)}HeatMap.dscalar.nii"
    nib.save(new_cifti_img, output_file)

    return None


for estimator in os.listdir("../results/FCs100824/Topography/SA"):
    if ("Naive" not in estimator) and (".R" not in estimator):
        print(estimator)
        reformat2nii(f"../results/FCs100824/Topography/SA/{estimator}")
