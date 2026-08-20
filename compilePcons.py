import numpy as np
import nibabel as nib
import pandas as pd
import dask.array as da
import dask.dataframe as dd
import tqdm

# 352 labels confirming dimensions
#gordon = nib.load("/spaces/ngdr/ref-data/abcd/nda-3165-2020-09/Gordon2014FreeSurferSubcortical_dparc.dlabel.nii")

# srun -N 1 -t 1:00:00 -p interactive --pty bash --mem=8G
# srun -N 1 --mem-per-cpu=32gb -t 1:00:00 -p interactive --pty bash 
# /spaces/regulated/abcd-ref-data/abcc/derivatives/abcd-hcp-pipeline_v0.1.4/sub-003RTV85/ses-00A/func/sub-003RTV85_ses-00A_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries_thresh-fd0p2mm_censor-10min_conndata-network_connectivity.pconn.nii
def extractPconn(subject) :
    basepath = "/spaces/ngdr/ref-data/abcd/nda-3165-2020-09/derivatives/abcd-hcp-pipeline-v0.1.3/sub-"
    filepath = f"{basepath}{subject}/ses-baselineYear1Arm1/func/sub-{subject}_ses-baselineYear1Arm1_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii_10_minutes_of_data_at_FD_0.2.pconn.nii"
    data = nib.load(filepath).get_fdata()
    # assert data.shape == (352, 352), f"Matrix is not 352x352, it is {data.shape}"
    # There are only 352 x 352 conectsion (so diagonal might already be removed, however there are 352 values in gordon)
    # replace all values of data with abs > 3 with np.nan
    #data = np.where(abs(data) > 3, np.nan, data) 
    # Add the one to not include the diagonal
    return data[np.triu_indices(data.shape[0], 1)]

#removeSubjects = pd.read_csv("derivatives_subjects_to_avoid_minus_dwi.csv", names=  ["IID", "Year"], sep = ",")
# Filter Year to just "baselineYear1Arm1"
#removeSubjects = removeSubjects[removeSubjects["Year"] == "baselineYear1Arm1"]
#Subjects= pd.read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/pconns/pconn_subs.files", names = ["IID"], header=None)
#Subjects = pd.read_csv("FC.files", names = ["IID"], header= None)
Subjects = pd.read_csv("/panfs/jay/groups/31/rando149/coffm049/Temppheno.csv")
# Remove the "_" in the IID column
Subjects["IID"] = Subjects["IID"].str.replace("_", "")
#Subjects["IID"] = Subjects["IID"].str.replace("sub-", "")

# Drop the subjects that have their IID appear in removeSubjects
# Subjects = Subjects[~Subjects["IID"].isin(removeSubjects["IID"])]

# rename the index when reading the file
FID = pd.read_table("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/IDs/IDs.txt", names = ["FID", "IID"], sep = " ")
FID["IID"] = FID["IID"].str.replace("_", "")

Subjects = pd.merge(Subjects, FID, how = "left", on = "IID").dropna()
Subjects["FID"] = Subjects["FID"].astype(int)
# Drop duplicates from Subjects
#Subjects = Subjects.drop_duplicates(subset = "FID")
subjects = Subjects.IID

# (352 * 352 - 352) / 2 = 61776
pconns = np.empty(shape = (len(subjects),61776), dtype = float)

pconns = pd.DataFrame(pconns)
pconns.index= subjects
####### #N#@ Code
deletes = []
for subj in tqdm.tqdm(subjects) :
    try :
        pconns.loc[subj] = extractPconn(subj)
    except FileNotFoundError :
        #deletes.append(ind)
        pconns.loc[subj] = np.nan
    except KeyError:
        pconns.loc[subj] = np.nan

# append letter o to each column
pconns.columns = pconns.columns.map(lambda x: f"o{x}")

pconns.to_parquet("pconns.parquet", index= True)
#pconns = np.delete(pconns, deletes, axis =1)
# Drop the rows denoted in the "deletes" list which contains the indices of the subjects that were not found
#Subjects = Subjects.drop(Subjects.index[deletes])
#pconns2 = dd.from_array(pconns, columns = Subjects.IID, chunksize = 620)

# append "o" to each row index value
#pconns2.index = pconns2.index.map(lambda x: f"o{x}")


# remove any rows of pconns that have more than 50% of their values as np.nan
#pconns = pconns.dropna(thresh=int(0.5*pconns.shape[1]))
#np.savetxt("pconns.csv", pconns, delimiter=",")

#dd.to_parquet(pconns2, "data/pconns", write_index= True)
#dd.to_csv(pconns2, "data/pconns", write_index= True)


