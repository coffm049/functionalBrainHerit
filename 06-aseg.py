import os

import nibabel as nib
import numpy as np
import pandas as pd

estFile = pd.read_csv("formattedData/twinVol.csv")
dlabelFile = nib.load("brainTemplate/Atlas_ROIs.2.nii")
dscalarFile = nib.load("brainTemplate/GordonBackground.dscalar.nii")        # cifti_template = nib.load(
labels = pd.read_csv("brainTemplate/gordon_modules.csv")
labels.columns = ["label"]

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
    df = estFile.groupby(["r1"]).agg({"h2": np.nanmedian})

    ############# LEFT OFF HERE

    vertex_scalar_map = np.zeros_like(dlabelFile.get_fdata(), dtype=float)

    vertex_label_map = np.zeros_like(dlabelFile.get_fdata(), dtype=int)
    for region, value in zip(df.index + 1, df.h2):
        print(region)
        vertex_scalar_map[dlabelFile.get_fdata() == region] = value
        vertex_label_map[dlabelFile.get_fdata() == region] = int(labels.label[region])

    # Create a new CIFTI-2 image with your scalar data
    new_cifti_img = nib.Cifti2Image(
        vertex_scalar_map,  # dlabelFile
        header=dlabelFile.header,
        nifti_header=dlabelFile.nifti_header,
    )
    if method == "GCTA":
        method = "SNP"
    # Save to a new .dscalar.nii file
    output_file = f"formattedData/roiMedian{method}-{pheno}.dscalar.nii"
    nib.save(new_cifti_img, output_file)

    return None


for estimator in os.listdir("formattedData"):
    if (
        ("roiMedianHerit" in estimator)
        and (".nii" not in estimator)
        and ("padded" not in estimator)
        and ("RE" not in estimator)
        and (".csv" in estimator)
    ):

        me, ph = estimator[14:].split("_")
        reformat2nii(method=me, pheno=ph)

reformat2nii(method=me, pheno=ph)
