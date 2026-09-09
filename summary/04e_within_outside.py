#!/usr/bin/env python3
"""
Within vs outside network heritability — 10th, 50th, 90th percentiles.

For each FC atlas (Gordon 352, ProbaConns 80) and each method (Twin, AdjHE-RE, 30 PCs):
  * Load networks per parcel (Gordon: gordon_modules.csv 352->14 via labelDict; Proba: abcd 0.75 via CIFTI + Sal<->SMl swap)
  * For each edge (oK -> i,j via triu_indices), classify as:
      - within: networks[i] == networks[j] and not NA
      - outside: networks[i] != networks[j] and neither NA
  * Collect h2 values for within and outside edges separately
  * Compute 10th, 50th (median), 90th percentiles (and mean, sd, n for context)

Outputs:
  results/summary/within_outside_h2.csv  (atlas, method, group, n, mean, sd, p10, p50, p90, min, max)
  Also per-atlas/per-method detailed CSVs if needed.

Run:
  .venv/bin/python summary/04e_within_outside.py
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import nibabel as nib

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
OUT_CSV = ROOT / "results/summary/within_outside_h2.csv"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)

GORDON_DLABEL = _resolve(ROOT.parent / "brainTemplates" / "Gordon.subcortical.32k_fs_LR.dlabel.nii",
                         "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii")
SHORT_DICT = _resolve(ROOT.parent / "brainTemplates" / "shortDicionary.csv",
                      "/users/4/coffm049/papers/brainTemplates/shortDicionary.csv")

def _cifti_labels(img):
    try:
        ax0 = img.header.get_axis(0)
        d = ax0.get_element(0)[1]
        if isinstance(d, dict) and len(d) > 1:
            return d
        if isinstance(d, list) and len(d) > 1:
            return {i: v for i, v in enumerate(d) if v is not None}
    except: pass
    return {}

def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    s = re.sub(r"^\d{1,3}_[LR]_", "", s)
    _map = {"Default":"DMN","Auditory":"Aud","CinguloOperc":"CO","FrontoParietal":"FP","Salience":"Sal","VentralAttn":"VAN","DorsalAttn":"DAN","Visual":"Vis","SMhand":"SMd","SMmouth":"SMl","RetrosplenialTemporal":"Tpole"}
    s = _map.get(s, s)
    return s.split("_")[0] if "_" in s else s

def load_gordon_networks(dlabel, csv_path, N):
    # Use gordon_modules.csv as ground truth (352 -> 14) per 02a
    for cand in [Path(__file__).resolve().parents[1].parent / "brainTemplates" / "gordon_modules.csv",
                 Path(r"C:\Users\coffm049\brainTemplates\gordon_modules.csv"),
                 Path("/users/4/coffm049/papers/brainTemplates/gordon_modules.csv")]:
        if cand.exists():
            try:
                gm = pd.read_csv(cand)
                vals = gm[gm.columns[0]].astype(int).tolist()
                if len(vals) == N:
                    labelDict = {1:"DMN",2:"VIS",3:"FP",4:"DAN",5:"VAN",6:"SAL",7:"CO",8:"SMD",9:"SML",10:"AUD",11:"Tpole",12:"MTL",13:"PMN",14:"PON"}
                    return np.array([labelDict.get(v, "NA") for v in vals])
            except: pass
    df = pd.read_csv(csv_path)
    name_to_short = dict(zip(df["name"].astype(str), df["shortname"].astype(str)))
    img = nib.load(str(dlabel))
    labels = _cifti_labels(img)
    nets=[]
    for p in range(1,N+1):
        lab=labels.get(p)
        s=lab[0] if isinstance(lab,(tuple,list)) else getattr(lab,"label",str(lab)) if lab else None
        net=_network_of(s)
        short=name_to_short.get(net,net)
        if short==net:
            for k,v in name_to_short.items():
                if k.lower()==net.lower():
                    short=v
                    break
        nets.append(short)
    return np.array(nets)

def find_proba_dlabel(N=80):
    cands=[ROOT.parent / "brainTemplates" / "abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii",
           Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"),
           Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.75.dlabel.nii")]
    for p in cands:
        if p.exists():
            try:
                img=nib.load(str(p))
                if len(img.header.get_axis(0).get_element(0)[1])-1==N:
                    return p
            except: pass
    return cands[0]

def load_proba_networks(dlabel, N):
    img=nib.load(str(dlabel))
    labels=_cifti_labels(img)
    nets=[]
    for p in range(1,N+1):
        lab=labels.get(p)
        s=lab[0] if isinstance(lab,(tuple,list)) else getattr(lab,"label",str(lab)) if lab else None
        net=_network_of(s)
        if net=="Sal": net="SMl"
        elif net=="SMl": net="Sal"
        nets.append(net)
    canon={"vis":"Vis","aud":"Aud","sal":"Sal","smd":"SMd","sml":"SMl","van":"VAN","dan":"DAN","dmn":"DMN","fp":"FP","co":"CO","mtl":"MTL","pmn":"PMN","pon":"PON","tpole":"Tpole"}
    return np.array([canon.get(n.lower(), n) for n in nets])

def pctiles(x, qs=[0.1,0.5,0.9]):
    x=np.array(x)
    x=x[np.isfinite(x)]
    if len(x)==0:
        return {f"p{int(q*100)}": np.nan for q in qs}
    return {f"p{int(q*100)}": float(np.nanpercentile(x, q*100)) for q in qs}

wide=pd.read_csv(WIDE)
results=[]
for atlas,N in [("gordon",352),("probaConns",80)]:
    if atlas=="gordon":
        networks=load_gordon_networks(GORDON_DLABEL, SHORT_DICT, N)
    else:
        dlabel=find_proba_dlabel(N)
        networks=load_proba_networks(dlabel, N)
    # Keep only cortical 14 for Gordon/Proba (filter NA)
    # But for within/outside classification, we need to handle NA: edges involving NA are excluded from both groups
    for method, col in [("Twin","Twin_h2"), ("AdjHE-RE", f"h2_{atlas}_AdjHE_RE")]:
        if col not in wide.columns:
            alt=col.replace("probaConns","proba")
            if alt in wide.columns:
                col=alt
            else:
                continue
        sub=wide[wide["Set"]==atlas][["Pheno",col]].dropna()
        if sub.empty:
            continue
        # Build mapping from pheno -> h2
        pheno_to_h2=dict(zip(sub["Pheno"].astype(str), sub[col].astype(float)))
        # For each edge, classify within vs outside
        within_vals=[]
        outside_vals=[]
        for pheno, h2 in pheno_to_h2.items():
            s=str(pheno)
            k=int(s[1:]) if s.startswith("o") else int(s)
            i,j=np.triu_indices(N,1)
            # k is index into triu
            # Need to decode: i[k],j[k]
            ii=int(i[k]); jj=int(j[k])
            ni=networks[ii] if ii < len(networks) else "NA"
            nj=networks[jj] if jj < len(networks) else "NA"
            if ni=="NA" or nj=="NA" or ni in ("SUB","???") or nj in ("SUB","???"):
                continue
            if ni == nj:
                within_vals.append(float(h2))
            else:
                outside_vals.append(float(h2))
        # Compute stats
        for group, vals in [("within", within_vals), ("outside", outside_vals)]:
            if len(vals)==0:
                continue
            arr=np.array(vals)
            # 10th, 50th, 90th
            p10, p50, p90 = np.nanpercentile(arr, [10,50,90])
            mean=np.nanmean(arr)
            sd=np.nanstd(arr, ddof=1) if len(arr)>1 else np.nan
            n=len(arr)
            results.append({"atlas": atlas, "method": method, "group": group, "n": int(n),
                            "mean": float(mean), "sd": float(sd) if np.isfinite(sd) else np.nan,
                            "p10": float(p10), "p50": float(p50), "p90": float(p90),
                            "min": float(np.nanmin(arr)), "max": float(np.nanmax(arr))})
            print(f"[{atlas} {method} {group}] n={n} mean={mean:.3f} p10={p10:.3f} p50={p50:.3f} p90={p90:.3f}")

out=pd.DataFrame(results)
# Round for display but keep full for CSV
out.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV} with {len(out)} rows")
# Also print formatted
print(out.to_string(index=False))
