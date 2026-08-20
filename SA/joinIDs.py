import pandas as pd


def Adding_IDs(saIDs, FIDs, SA, output):
    saIDs = pd.read_csv(saIDs, header=None)
    saIDs.columns = ['IID']
    saIDs['IID'] = saIDs['IID'].str.replace('sub-NDARINV', 'NDAR_INV')
    FIDs = pd.read_csv(FIDs, header=None, sep="\s+")
    FIDs.columns = ["FID", "IID"]
    saIDs = pd.merge(saIDs, FIDs, on="IID", how="left")
    SA = pd.read_csv(SA, header=0)
    pd.concat([saIDs.iloc[:, 0:2], SA], axis=1).to_csv(output, index=False, header=True, sep=" ")


Adding_IDs(saIDs="/panfs/jay/groups/4/miran045/shared/projects/ABCD_topography/data/Topography/IDs/ABCD_group1IDs.csv",
           FIDs="/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
           SA="/panfs/jay/groups/4/miran045/shared/projects/ABCD_topography/data/Topography/Surfacearea/surfarea.csv",
           output="~/Temppheno.csv")