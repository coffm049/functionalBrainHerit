#!/usr/bin/env python3
"""
Correlation between FC heritability (per-parcel average / 90th percentile) and SA heritability by system.

For each atlas (Gordon 352, ProbaConns 80) and each method (Twin, AdjHE-RE, 30 PCs):
  * Rebuild the full N x N FC h2 matrix from mash_twin_wide.csv (oK -> triu)
  * For each parcel p, compute:
      - mean_h2[p] = mean of row p (all edges incident to p, nan on diagonal)
      - p90_h2[p]  = 90th percentile of row p
  * Map parcels -> network via the same lookup used for brain/Manhattan:
      Gordon: gordon_modules.csv (352 -> 14) via labelDict, or LUT+shortDict fallback
      Proba:  abcd_template_matching_combined_clusters_thresh0.75 (80) via CIFTI label table + Sal<->SMl swap
  * Aggregate by network: for each network, mean of parcel-level mean_h2 / p90_h2 across parcels in that network
  * Merge with SA h2 per network (mash_twin_wide.csv SA, 14 networks excl. 4/6/17) for same method
  * Compute Pearson r and Spearman rho across the 14 networks

Outputs:
  results/summary/sa_fc_h2_correlation.csv  (atlas, method, metric, n_networks, pearson_r, spearman_rho, p_pearson, etc.)
  Also per-atlas scatter plots if requested.

Run:
  .venv/bin/python summary/04b_h2_SA_correlation.py
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import nibabel as nib

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
OUT_CSV = ROOT / "results/summary/sa_fc_h2_correlation.csv"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

# SA 14 networks (excludes 4,6,17 null)
SA_NETWORKS = {1:"DMN",2:"VIS",3:"FP",5:"DAN",7:"VAN",8:"SAL",9:"CO",10:"SMD",11:"SML",12:"AUD",13:"Tpole",14:"MTL",15:"PMN",16:"PON"}

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
        elem = ax0.get_element(0)
        d = elem[1]
        if isinstance(d, dict) and len(d) > 1:
            return d
        if isinstance(d, list) and len(d) > 1:
            return {i: v for i, v in enumerate(d) if v is not None}
    except Exception:
        pass
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
    # Use gordon_modules.csv as ground truth for Gordon (as in 02a)
    cands = [Path(__file__).resolve().parents[1].parent / "brainTemplates" / "gordon_modules.csv",
             Path(r"C:\Users\coffm049\brainTemplates\gordon_modules.csv"),
             Path("/users/4/coffm049/papers/brainTemplates/gordon_modules.csv"),
             ROOT.parent / "brainTemplates" / "gordon_modules.csv"]
    for cand in cands:
        if cand.exists():
            try:
                gm = pd.read_csv(cand)
                col = gm.columns[0]
                vals = gm[col].astype(int).tolist()
                if len(vals) == N:
                    labelDict = {1:"DMN",2:"VIS",3:"FP",4:"DAN",5:"VAN",6:"SAL",7:"CO",8:"SMD",9:"SML",10:"AUD",11:"Tpole",12:"MTL",13:"PMN",14:"PON"}
                    return np.array([labelDict.get(v, "NA") for v in vals])
            except Exception:
                pass
    # Fallback CIFTI
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
    # canonicalize
    canon={"vis":"Vis","aud":"Aud","sal":"Sal","smd":"SMd","sml":"SMl","van":"VAN","dan":"DAN","dmn":"DMN","fp":"FP","co":"CO","mtl":"MTL","pmn":"PMN","pon":"PON","tpole":"Tpole"}
    return np.array([canon.get(n.lower(), n) for n in nets])

def rebuild_pconn(phenos, edge_h2, N):
    k=np.array([int(str(p)[1:]) if str(p).startswith("o") else int(str(p)) for p in phenos])
    M=np.zeros((N,N))
    i,j=np.triu_indices(N,1)
    M[i[k], j[k]]=edge_h2
    M[j[k], i[k]]=edge_h2
    np.fill_diagonal(M, np.nan)
    return M

# Load SA h2 per network (14)
wide=pd.read_csv(WIDE)
sa_sub=wide[wide["Set"]=="SA"]
# Map pheno -> network via SA_NETWORKS (14) — use method keys "Twin" and "AdjHE-RE" for later lookup
sa_map={}
for _,row in sa_sub.iterrows():
    pheno=row["Pheno"]
    try:
        num=int(str(pheno).split("network_surfarea")[-1]) if "network_surfarea" in str(pheno) else int(str(pheno))
        net=SA_NETWORKS.get(num)
        if net is None:
            continue
        sa_map[net] = {"Twin": float(row["Twin_h2"]) if pd.notna(row["Twin_h2"]) else np.nan,
                       "AdjHE-RE": float(row["h2_SA_AdjHE_RE"]) if pd.notna(row["h2_SA_AdjHE_RE"]) else np.nan}
    except: pass

# Build FC per-parcel summaries and correlate
results=[]
for atlas,N in [("gordon",352),("probaConns",80)]:
    # Load networks for this atlas
    if atlas=="gordon":
        networks=load_gordon_networks(GORDON_DLABEL, SHORT_DICT, N)
    else:
        dlabel=find_proba_dlabel(N)
        networks=load_proba_networks(dlabel, N)
    # Map network -> parcel indices (case-insensitive for SA)
    net_to_idx={}
    for idx, net in enumerate(networks):
        net_to_idx.setdefault(net, []).append(idx)
    # Keep only the 14 SA networks (case-insensitive: VIS vs Vis, SAL vs Sal, etc.)
    sa_lower={k.lower(): k for k in sa_map}
    networks_keep=[]
    for n in net_to_idx:
        # Find matching SA network case-insensitively
        cand=sa_lower.get(n.lower())
        if cand is not None:
            networks_keep.append(n)
    # Also ensure canonical SA name for merging
    print(f"[{atlas}] networks {sorted(set(networks))} -> keep {sorted(networks_keep)}")

    for method, col in [("Twin","Twin_h2"), ("AdjHE-RE", f"h2_{atlas}_AdjHE_RE" if atlas!="SA" else "h2_SA_AdjHE_RE")]:
        if col not in wide.columns:
            alt=col.replace("probaConns","proba")
            if alt in wide.columns:
                col=alt
            else:
                continue
        sub=wide[wide["Set"]==atlas][["Pheno",col]].dropna()
        if sub.empty:
            continue
        M=rebuild_pconn(sub["Pheno"].values, sub[col].values.astype(float), N)
        # Per-parcel mean and p90
        parcel_mean=[]
        parcel_p90=[]
        for p in range(N):
            vals=M[p, :][~np.isnan(M[p,:])]
            if len(vals)==0:
                parcel_mean.append(np.nan)
                parcel_p90.append(np.nan)
            else:
                parcel_mean.append(float(np.nanmean(vals)))
                parcel_p90.append(float(np.nanpercentile(vals, 90)))
        parcel_mean=np.array(parcel_mean)
        parcel_p90=np.array(parcel_p90)

        # Aggregate by network (case-insensitive SA lookup: VIS vs Vis, SAL vs Sal, etc.)
        # Also compute within-system FC h2: mean of edges where both parcels are in same network
        sa_lower2={k.lower(): v for k,v in sa_map.items()}
        rows=[]
        for net in networks_keep:
            idxs=net_to_idx[net]
            # Find SA entry case-insensitively
            sa_key=sa_lower2.get(net.lower())
            if sa_key is None:
                # Try direct
                sa_key=sa_map.get(net)
                if sa_key is None:
                    continue
                sa_h2_entry=sa_key
            else:
                sa_h2_entry=sa_key
            # sa_h2_entry is dict with Twin/AdjHE-RE
            if isinstance(sa_h2_entry, dict):
                sa_h2=sa_h2_entry.get(method, np.nan)
            else:
                sa_h2=np.nan
            # Fallback to direct sa_map[net] if needed
            if not np.isfinite(sa_h2):
                # Try case-insensitive direct
                for k,v in sa_map.items():
                    if k.lower()==net.lower():
                        sa_h2=v.get(method, np.nan)
                        break
            fc_mean=np.nanmean([parcel_mean[i] for i in idxs])
            fc_p90=np.nanmean([parcel_p90[i] for i in idxs])
            # Within-system: edges where both ends in same network
            within_vals=[]
            for ii in idxs:
                for jj in idxs:
                    if ii < jj and not np.isnan(M[ii, jj]):
                        within_vals.append(M[ii, jj])
            fc_within=np.nanmean(within_vals) if len(within_vals)>0 else np.nan
            rows.append((net, sa_h2, fc_mean, fc_p90, fc_within))
        df_net=pd.DataFrame(rows, columns=["network","sa_h2","fc_mean","fc_p90","fc_within"]).dropna()
        if len(df_net) < 3:
            continue
        # Correlations
        for metric in ["fc_mean","fc_p90","fc_within"]:
            x=df_net["sa_h2"].values
            y=df_net[metric].values
            # Pearson and Spearman
            try:
                r=np.corrcoef(x,y)[0,1]
            except: r=np.nan
            try:
                from scipy.stats import spearmanr, pearsonr
                rho,_=spearmanr(x,y)
                # p-value for pearson
                _, pval = pearsonr(x,y)
            except:
                rho=np.nan
                pval=np.nan
            results.append({"atlas": atlas, "method": method, "metric": metric, "n_networks": int(len(df_net)),
                            "pearson_r": float(r) if np.isfinite(r) else np.nan,
                            "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                            "p_pearson": float(pval) if np.isfinite(pval) else np.nan})
            print(f"[{atlas} {method} {metric}] n={len(df_net)} r={r:.3f} rho={rho:.3f}")

out=pd.DataFrame(results)
out.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV} with {len(out)} rows")
print(out.to_string(index=False))
