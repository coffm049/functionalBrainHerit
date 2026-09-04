#!/usr/bin/env python3
"""
Heatmap of FC heritability (h2) ordered by Sys-Sys network labels.

For each atlas (Gordon 352, ProbaConns 80) and each method (Twin, AdjHE-RE, 30 PCs):
  * Rebuild the full N x N h2 matrix from mash_twin_wide.csv
  * Order parcels by network (as in circular plot: argsort(networks))
  * Plot heatmap with networks as blocks, ordered by Sys-Sys labels, red (low) -> yellow (high), 0-1
  * Save as results/summary/plots/heatmap_{atlas}_{method}.png

Publication-ready: minimal whitespace, no title numbers, stats to CSV, AdjHE-RE naming, 0-1 colorbar.

Run:
  .venv/bin/python summary/04c_FC_heatmap.py
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
PLOT_DIR = ROOT / "results/summary/plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

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
    # Use gordon_modules.csv as ground truth for Gordon (as in 02a)
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

def rebuild_pconn(phenos, edge_h2, N):
    k=np.array([int(str(p)[1:]) if str(p).startswith("o") else int(str(p)) for p in phenos])
    M=np.zeros((N,N))
    i,j=np.triu_indices(N,1)
    M[i[k], j[k]]=edge_h2
    M[j[k], i[k]]=edge_h2
    np.fill_diagonal(M, np.nan)
    return M

wide=pd.read_csv(WIDE)
for atlas,N in [("gordon",352),("probaConns",80)]:
    if atlas=="gordon":
        networks=load_gordon_networks(GORDON_DLABEL, SHORT_DICT, N)
    else:
        dlabel=find_proba_dlabel(N)
        networks=load_proba_networks(dlabel, N)
    # Order parcels by network (as in circular)
    order=np.argsort(networks, kind="stable")
    networks_ordered=networks[order]
    # Unique networks in order
    uniq, counts = np.unique(networks_ordered, return_counts=True)
    print(f"[{atlas}] order {networks_ordered[:10]}... uniq {list(zip(uniq, counts))}")

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
        M=rebuild_pconn(sub["Pheno"].values, sub[col].values.astype(float), N)
        # Order matrix
        M_ordered=M[np.ix_(order, order)]
        # Publication-ready heatmap: red (low) -> yellow (high), 0-1 as requested (was 0-0.5)
        fig, ax = plt.subplots(figsize=(6, 6))
        red_yellow = mcolors.LinearSegmentedColormap.from_list("red_yellow", ["#d73027", "#ffff00"])
        im=ax.imshow(M_ordered, cmap=red_yellow, vmin=0, vmax=1.0, interpolation="nearest", aspect="auto")
        # Add network boundaries
        cum=np.cumsum(counts)[:-1]
        for c in cum:
            ax.axhline(c-0.5, color="black", linewidth=0.5, alpha=0.7)
            ax.axvline(c-0.5, color="black", linewidth=0.5, alpha=0.7)
        # Ticks: network labels at middle of each block
        tick_pos=np.cumsum(counts) - counts/2
        # Shorten labels for display
        ax.set_xticks(tick_pos)
        ax.set_yticks(tick_pos)
        ax.set_xticklabels(uniq, rotation=45, ha="right", fontsize=6, fontweight="bold")
        ax.set_yticklabels(uniq, fontsize=6, fontweight="bold")
        ax.tick_params(axis="both", which="both", length=0)
        # No title for publication (handled in quarto)
        cbar=fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label(r"Heritability ($h^2$)", fontsize=9)
        cbar.ax.tick_params(labelsize=7)
        fig.tight_layout(pad=0.5)
        out=PLOT_DIR / f"heatmap_{atlas}_{method.replace('-','')}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        # Save stats
        try:
            pd.DataFrame([{"atlas": atlas, "method": method, "n": int(N), "h2_mean": float(np.nanmean(M)), "h2_max": float(np.nanmax(M))}]).to_csv(str(out).replace(".png","_stats.csv"), index=False)
        except: pass
        print(f"  wrote {out}")

print("Done heatmaps")
