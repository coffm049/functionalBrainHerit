import os
import numpy as np
import nibabel as nib
import pandas as pd
from tqdm import tqdm

N_REGIONS = 80
N_EDGES = N_REGIONS * (N_REGIONS - 1) // 2

FILES = "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/proba.files"
OUT = "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet"


def normalize_id(s):
    s = str(s)
    for prefix in ("sub-", "NDAR_INV", "NDARINV"):
        if s.startswith(prefix):
            s = s[len(prefix):]
    return s


def extractPconn(filepath):
    data = nib.load(filepath).get_fdata()
    vec = data[np.triu_indices(data.shape[0], 1)]
    if vec.shape[0] != N_EDGES:
        raise ValueError(f"Unexpected edge count {vec.shape[0]}")
    return vec


def main():
    with open(FILES) as f:
        paths = [ln.strip() for ln in f if ln.strip()]

    rows = pd.DataFrame({
        "path": paths,
        "IID": [normalize_id(os.path.basename(p).split("_", 1)[0]) for p in paths],
    }).drop_duplicates(subset="IID").reset_index(drop=True)

    pconns = pd.DataFrame(np.nan, index=rows["IID"],
                          columns=[f"o{i}" for i in range(N_EDGES)], dtype=float)

    for iid, path in tqdm(zip(rows["IID"], rows["path"]), total=len(rows)):
        try:
            pconns.loc[iid] = extractPconn(path)
        except (FileNotFoundError, KeyError, OSError, ValueError):
            pass

    pconns.to_parquet(OUT, index=True)
    print(f"Wrote {pconns.shape[0]} x {pconns.shape[1]} -> {OUT}")


if __name__ == "__main__":
    main()