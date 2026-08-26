import re

import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

avoid = pd.read_csv(
    "data/derivatives_subjects_to_avoid_minus_dwi.csv", names=["IID", "year"]
)

# filter to where baseline is in year but isn't the whole string
avoid = avoid[avoid.year.str.contains("base")]

geno = pd.read_table("data/no_rels.grm.id", names=["FID", "IID"])

# subset to only one observation for each FID in geno
geno = geno.groupby("FID").first()

probas = pd.read_csv("data/proba.files", names=["file"])
probas.file = probas.file.str.extract(r"sub-([^ ]*?)_ses")
probas["IID"] = probas.file.str.slice(0, 4) + "_" + probas.file.str.slice(4)
fcFiles = pd.read_csv("data/FC.files", names=["IID"])

Abc = len(set(geno["IID"]) - set(probas["IID"]) - set(fcFiles["IID"]) - set(avoid.IID))
aBc = len(set(probas["IID"]) - set(geno["IID"]) - set(fcFiles["IID"]) - set(avoid.IID))
abC = len(set(fcFiles["IID"]) - set(geno["IID"]) - set(probas["IID"]) - set(avoid.IID))
ABc = len(set(geno["IID"]) & set(probas["IID"]) - set(fcFiles["IID"]) - set(avoid.IID))
AbC = len(set(geno["IID"]) & set(fcFiles["IID"]) - set(probas["IID"]) - set(avoid.IID))
aBC = len(set(probas["IID"]) & set(fcFiles["IID"]) - set(geno["IID"]) - set(avoid.IID))
ABC = len(set(geno["IID"]) & set(probas["IID"]) & set(fcFiles["IID"]) - set(avoid.IID))


matplotlib.use("TkAgg")  # or 'Qt5Agg', 'GTK3Agg', etc.

# (Abc, aBc, ABc, abC, AbC, aBC, ABC)
v = venn3(
    subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
    set_labels=(
        f"Genotyping \n ({geno.shape[0]})",
        f"Topo + PFM\n ({probas.shape[0]})",
        f"Gordon\n ({fcFiles.shape[0]})",
    ),
)


# assign title
plt.title("ABCD (n=11,878)")
plt.show()

# save the figure
plt.savefig("figures/venn.png", dpi=300)
