import itertools
import math

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from numpy._core.defchararray import index

# df = pd.read_csv("data/allEstimates2.csv")  # .fillna(0)
df = pd.read_csv("data/allEstimates.csv")  # .fillna(0)
labelDictionary = pd.DataFrame(
    {
        "index": list(range(1, 15)),
        "name": [
            "DMN",
            "VIS",
            "FP",
            "DAN",
            "VAN",
            "SAL",
            "CO",
            "SMD",
            "SML",
            "AUD",
            "Tpole",
            "MTL",
            "PMN",
            "PON",
        ],
    }
)
# get label names
labels1 = pd.read_csv("brainTemplate/GRP1_template_parcel.csv")
labels1["name"] = [labelDictionary.loc[i - 1, "name"] for i in labels1.region]
nameMat = np.array([f"{i}-{j}" for i in labels1.name for j in labels1.name]).reshape(
    80, 80
)

test = np.triu_indices_from(nameMat, 1)
longNames1 = pd.DataFrame(nameMat[test], columns=["connection"])
longNames1["phenoClass"] = "PFN"
longNames1["pheno"] = range(longNames1.shape[0])
longNames1["r1"] = test[0]
longNames1["r2"] = test[1]

labels2 = pd.read_csv("brainTemplate/gordon_modules.csv")
labels2["name"] = [labelDictionary.loc[i - 1, "name"] for i in labels2.region]
nameMat = np.array([f"{i}-{j}" for i in labels2.name for j in labels2.name]).reshape(
    labels2.shape[0], labels2.shape[0]
)

test = np.triu_indices_from(nameMat, 1)
longNames2 = pd.DataFrame(nameMat[test], columns=["connection"])
longNames2["phenoClass"] = "Gordon"
longNames2["pheno"] = range(longNames2.shape[0])
longNames2["r1"] = test[0]
longNames2["r2"] = test[1]

longNames = pd.concat([longNames1, longNames2]).reset_index(drop=True)
test = pd.merge(df, longNames, how="left", on=["phenoClass", "pheno"])


test.to_csv("data/allEstimatesLabeled.csv", index=False)
