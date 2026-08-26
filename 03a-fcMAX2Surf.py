import itertools
import os

import nibabel as nib
import numpy as np
import pandas as pd

# import h5py
# import nilearn
# from nilearn import plotting

# Load an existing CIFTI-2 template file to get the correct structure (like the HCP template)
# cifti_template = nib.load(
#     "brainTemplate/91282_Greyordinates.dscalar.nii"
# )  # Replace with the actual path to your template
# dlabelFile = nib.load(
#     "brainTemplate/abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"
# )
cifti_template = nib.load(
    "brainTemplate/single_map.dscalar.nii"
)  # Replace with the actual path to your template


def reformat2nii(method, pheno):
    estFile = (
        pd.read_csv("data/allEstimatesLabeled.csv")
        .query(f"phenoClass == '{pheno}'")
        .query(f"Type == '{method}'")
    )
    if pheno == "Gordon":
        # [ ] convert dlablel to dsclaar
        # - FEZ sent a script
        # - sex effects are likely batch effects
        # [ ] use workbench to make sure conversion looks well dscalar shoould look same as dlabel
        labels = pd.read_csv("brainTemplate/gordon_modules.csv")
        labels.columns = ["label"]
        # dscalarFile = nib.load("brainTemplate/Gordon.subcortical.32k_fs_LR.dlabel.nii")
        # cifti_template = nib.load(
        #     "brainTemplate/Gordon.subcortical.32k_fs_LR.dscalar.nii"
        # )
        cifti_template = nib.load(
            "brainTemplate/single_map.dscalar.nii"
        ) 
        dscalarFile = nib.load("brainTemplate/Gordon.networks.32k_fs_LR.dlabel.nii")

    elif pheno == "PFN":
        labels = pd.read_csv("brainTemplate/GRP1_template_parcel.csv")
        labels.columns = ["label"]
        # dlabelFile = nib.load("brainTemplate/91282_Greyordinates.dscalar.nii")
        # dlabelFile = nib.load("formattedData/networkSurfLabels.dscalar.nii")
        dscalarFile = nib.load("formattedData/networkSurfLabels.dscalar.nii")
        cifti_template = dscalarFile
        # dlabelFile = nib.load("brainTemplate/91282_Greyordinates.dscalar.nii")
        # replace dlabel with "networksurfLabels.dlabel.nii
        # might need to convert this to dsclar first

    else:
        print("File needs to specify either PFN or Gordon as the atlas used")

    # label_data = dlabelFile.get_fdata().astype(
    #     np.float32
    # )  # shape: (n_grayordinates, n_maps)
    # # axis0 = dlabelFile.header.get_axis(0)  # Time/Map axis — likely length 1

    # axis0 = nib.cifti2.ScalarAxis(["converted_from_label"])
    # axis1 = dlabelFile.header.get_axis(1)  # BrainModel axis

    # new_hdr = nib.cifti2.Cifti2Header.from_axes((axis0, axis1))  # Same axes, new intent

    # # Step 4: Create new image with scalar intent
    # dscalar_img = nib.cifti2.Cifti2Image(dataobj=label_data, header=new_hdr)

    # # Step 5: Save it as a true .dscalar.nii
    # nib.save(dscalar_img, "formattedData/networkSurfLabels.dscalar.nii")
    # print("✅ Saved as true dscalar: converted.dscalar.nii")

    # create new column containing first value from connection separated by -
    # for each row
    estFile[["l1", "l2"]] = estFile["connection"].str.split("-", n=1, expand=True)

    # summarize by roi
    # df = estFile.groupby(["phenoClass", "Type", "r1"]).agg({"h2": np.nanmedian})
    df = estFile.groupby(["r1"]).agg({"h2": np.nanmax})

    ############# LEFT OFF HERE

    vertex_scalar_map = np.zeros_like(cifti_template.get_fdata(), dtype=float)

    vertex_label_map = np.zeros_like(dscalarFile.get_fdata(), dtype=int)
    for region, value in zip(df.index + 1, df.h2):
        # print(region)
        vertex_scalar_map[dscalarFile.get_fdata() == region] = value
        vertex_label_map[dscalarFile.get_fdata() == region] = int(labels.label[region])

    # Create a new CIFTI-2 image with your scalar data
    new_cifti_img = nib.Cifti2Image(
        vertex_scalar_map,  # dlabelFile
        header=cifti_template.header,
        nifti_header=cifti_template.nifti_header,
    )
    if method == "GCTA":
        method = "SNP"
    # Save to a new .dscalar.nii file
    output_file = f"formattedData/roiMaxh2{method}-{pheno}.dscalar.nii"
    nib.save(new_cifti_img, output_file)

    return None


for me, ph in itertools.product(["GCTA", "twin"], ["Gordon", "PFN"]):
    print(me, ph)
    reformat2nii(method=me, pheno=ph)
