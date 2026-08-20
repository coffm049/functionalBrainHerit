import os
import numpy as np
import nibabel as nib
import pandas as pd


def extractPconn(filepath) :
    data = nib.load(filepath).get_fdata()
    # Add the one to not include the diagonal
    return data[np.triu_indices(data.shape[0], 1)]

# removeSubjects = pd.read_csv("derivatives_subjects_to_avoid_minus_dwi.csv", names=  ["IID", "Year"], sep = ",")
# Filter Year to just "baselineYear1Arm1"
# removeSubjects = removeSubjects[removeSubjects["Year"] == "baselineYear1Arm1"]
Subjects = pd.read_csv("proba.files", names = ["IID"], header= None)
filepaths = Subjects.copy()
Subjects["IID"] = Subjects["IID"].apply(os.path.basename)
# Remove the "_" in the IID column
Subjects["IID"] = Subjects["IID"].apply(lambda x : x.split("_", 1)[0].split("-", 1)[1])

# Drop the subjects that have their IID appear in removeSubjects
#Subjects = Subjects[~Subjects["IID"].isin(removeSubjects["IID"])]
# rename the index when reading the file
FID = pd.read_table("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/IDs/IDs.txt", names = ["FID", "IID"], sep = " ")
FID["IID"] = FID["IID"].str.replace("_", "")

Subjects = pd.merge(Subjects, FID, how = "left", on = "IID").dropna()
Subjects["FID"] = Subjects["FID"].astype(int)
# Drop duplicates from Subjects
#Subjects = Subjects.drop_duplicates(subset = "FID")
subjects = Subjects.IID


pconns = np.empty(shape = (len(subjects), 3160), dtype = float)


deletes = []
for ind, subj in enumerate(subjects) :
    print(ind)
    try :
        pconns[ind, :] = extractPconn(filepaths["IID"][ind])
    except FileNotFoundError :
        pconns.loc[subj] = np.nan
    except KeyError:
        pconns.loc[subj] = np.nan

#pconns = np.delete(pconns, deletes, axis =1)
# Drop the rows denoted in the "deletes" list which contains the indices of the subjects that were not found
#Subjects = Subjects.drop(Subjects.index[deletes])
# pconns = dd.from_array(pconns, columns = Subjects.IID, chunksize = 316)

# append "o" to each row index value
pconns = pd.DataFrame(pconns)
pconns.index = subjects
pconns.columns = pconns.columns.map(lambda x: f"o{x}")

# replace any value within pconns with abs(data) > 3 by np.nan
#pconns = pconns.map_partitions(lambda df: df.applymap(lambda x: np.nan if abs(x) > 3 else x))

# remove any rows of pconns that have more than 50% of their values as np.nan
#pconns = pconns.dropna(thresh=int(0.5*pconns.shape[1]))

pconns.to_parquet("probaConns.parquet", index= True)

from pyarrow.parquet import ParquetDataset
a = ParquetDataset("probaConns.parquet")

for piece in a.fragments:
    nrows += piece.get_metadata().num_rows
