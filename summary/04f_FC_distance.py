#!/usr/bin/env python3
"""
FC heritability vs distance — per edge, both Twin and AdjHE-RE, Euclidean + geodesic.

For each atlas (Gordon 352, ProbaConns 80) and each method (Twin, AdjHE-RE, 30 PCs):
  * Rebuild N x N h2 matrix from mash_twin_wide.csv
  * Compute parcel centroids from Conte69 surfaces + dlabel (as in 02a/b)
  * For each edge (i,j) with h2, compute:
      - euclidean_dist = ||centroid[i] - centroid[j]|| (MNI, mm)
      - geodesic_dist = surface geodesic via wb_command if available, else np.nan
  * Save per-edge CSV and simple scatter plot (h2 vs distance, no GAM/bins).

Outputs:
  results/summary/fc_distance_per_edge.csv  (atlas, method, i, j, h2, euclidean_dist, geodesic_dist, same_network)
  results/summary/plots/fc_distance_{atlas}_{method}.png  (scatter, per edge, Twin vs AdjHE-RE separate)

Run:
  .venv/bin/python summary/04f_FC_distance.py
"""
from pathlib import Path
import re
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib

ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "results/summary/mash_twin_wide.csv"
OUT_CSV = ROOT / "results/summary/fc_distance_per_edge.csv"
PLOT_DIR = ROOT / "results/summary/plots"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
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

def parcel_centroids(dlabel, N, surf_l, surf_r):
    """Compute centroids and centroid vertices for each parcel (1..N) from Conte69 surfaces + dlabel."""
    img=nib.load(str(dlabel))
    label_data=np.asarray(img.get_fdata()).astype(int).ravel()
    ax=img.header.get_axis(1)
    va=np.asarray(ax.vertex, dtype=int)
    # Load surfaces
    try:
        coords_l = nib.load(str(surf_l)).darrays[0].data
        coords_r = nib.load(str(surf_r)).darrays[0].data
    except Exception as e:
        print(f"  centroid: failed to load surfaces {e}, using Euclidean via vertex indices")
        coords_l = coords_r = None
    # Collect coords and vertex indices per parcel
    parcel_coords={p: [] for p in range(1,N+1)}
    parcel_vertex={p: [] for p in range(1,N+1)}
    parcel_hemi={}
    for structure, sl, _ in ax.iter_structures():
        s=str(structure).upper()
        if sl.stop is None:
            sl=slice(sl.start, len(va))
        blk=va[sl]
        seg=label_data[sl]
        hemi="L" if "CORTEX_LEFT" in s else "R" if "CORTEX_RIGHT" in s else None
        if hemi is None:
            continue
        for v,p in zip(blk, seg):
            if p==0 or v<0 or p> N: continue
            coords = coords_l[v] if hemi=="L" and coords_l is not None and v < len(coords_l) else coords_r[v] if hemi=="R" and coords_r is not None and v < len(coords_r) else None
            if coords is not None:
                parcel_coords[p].append(coords)
                parcel_vertex[p].append((hemi, int(v)))
    centroids={}
    centroid_vertex={}
    for p in range(1,N+1):
        coords=parcel_coords[p]
        verts=parcel_vertex[p]
        if len(coords)>0:
            centroids[p]=np.mean(coords, axis=0)
            # Choose vertex closest to centroid as representative
            dists=[np.linalg.norm(np.array(c)-centroids[p]) for c in coords]
            centroid_vertex[p]=verts[int(np.argmin(dists))]
        else:
            centroids[p]=np.array([np.nan, np.nan, np.nan])
            centroid_vertex[p]=(None, -1)
    return centroids, centroid_vertex

def geodesic_distances_via_wb(surf_l, surf_r, centroids, centroid_vertex, tmpdir=Path("/tmp")):
    """Try wb_command geodesic distance via python subprocess; fallback to nan.

    Returns a cache dict and helper to compute geodesic per edge.
    We precompute for each source parcel's centroid vertex via wb_command once,
    caching the full distance map, so we need at most N calls (352), not N^2.
    """
    try:
        result=subprocess.run(["wb_command", "-help"], capture_output=True, timeout=5)
        has_wb=result.returncode==0
    except:
        has_wb=False
    if not has_wb:
        print("  wb_command not found, geodesic will be NaN (Euclidean only)")
        return None
    def _find_midthickness(hemi):
        cands = [
            surf_l.parent / f"Conte69.{hemi}.midthickness.32k_fs_LR.surf.gii",
            Path(f"/users/4/coffm049/papers/brainTemplates/Conte69.{hemi}.midthickness.32k_fs_LR.surf.gii"),
            Path(f"C:/Users/coffm049/brainTemplates/Conte69.{hemi}.midthickness.32k_fs_LR.surf.gii"),
        ]
        for c in cands:
            if c.exists():
                return c
        return surf_l if hemi=="L" else surf_r
    has_left = any(v[0]=="L" for v in centroid_vertex.values() if v[0] is not None)
    has_right = any(v[0]=="R" for v in centroid_vertex.values() if v[0] is not None)
    print(f"  wb_command found — will compute geodesic for same-hemisphere pairs (L:{has_left} R:{has_right}) via centroid vertices")
    tmpdir.mkdir(parents=True, exist_ok=True)
    return {"centroid_vertex": centroid_vertex, "has_wb": True, "surf_l": _find_midthickness("L"), "surf_r": _find_midthickness("R"), "tmpdir": tmpdir, "cache": {}}

# Main
SURF_L = _resolve(ROOT.parent / "brainTemplates" / "Conte69.L.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii")
SURF_R = _resolve(ROOT.parent / "brainTemplates" / "Conte69.R.inflated.32k_fs_LR.surf.gii",
                  "/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii")

def _geodesic_for_pair(atlas, i, j, centroids, centroid_vertex, geo_cache):
    """Compute geodesic distance for parcels i,j (0-indexed) via wb_command if available."""
    if geo_cache is None or not geo_cache.get("has_wb"):
        return np.nan
    # Only same-hemisphere geodesic is meaningful
    hemi_i, vi = centroid_vertex.get(i+1, (None, -1))
    hemi_j, vj = centroid_vertex.get(j+1, (None, -1))
    if hemi_i is None or hemi_j is None or hemi_i != hemi_j or vi < 0 or vj < 0:
        return np.nan
    surf = geo_cache["surf_l"] if hemi_i == "L" else geo_cache["surf_r"]
    # Cache key: (hemi, vi)
    cache_key = (hemi_i, vi)
    if cache_key not in geo_cache["cache"]:
        # Run wb_command for this source vertex
        tmpdir = geo_cache["tmpdir"]
        out_metric = tmpdir / f"geodesic_{atlas}_{hemi_i}_{vi}.func.gii"
        # Also need a dscalar for output? Use metric
        # wb_command -surface-geodesic-distance <surf> <vertex> <out>
        # On some versions, need -borders or -limit
        try:
            # Use midthickness if available, else inflated
            cmd = ["wb_command", "-surface-geodesic-distance", str(surf), str(vi), str(out_metric)]
            # print(f"  wb_command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, timeout=30)
            if result.returncode != 0:
                # Try alternative: -surface-geodesic-distance with -borders
                # For now, just fail and return nan
                geo_cache["cache"][cache_key] = None
                return np.nan
            # Read the metric file (contains distance per vertex)
            # The metric is a func.gii with one column, we can read via nibabel
            try:
                metric = nib.load(str(out_metric))
                # func.gii has darrays
                data = metric.darrays[0].data
                geo_cache["cache"][cache_key] = np.array(data, dtype=float)
            except Exception as e:
                geo_cache["cache"][cache_key] = None
                return np.nan
        except Exception as e:
            geo_cache["cache"][cache_key] = None
            return np.nan
    dist_map = geo_cache["cache"].get(cache_key)
    if dist_map is None or vj >= len(dist_map):
        return np.nan
    return float(dist_map[vj])

wide=pd.read_csv(WIDE)
rows=[]
for atlas,N in [("gordon",352),("probaConns",80)]:
    # Load networks and centroids
    if atlas=="gordon":
        dlabel=GORDON_DLABEL
        networks=load_gordon_networks(dlabel, SHORT_DICT, N)
    else:
        dlabel=find_proba_dlabel(N)
        networks=load_proba_networks(dlabel, N)
    centroids, centroid_vertex = parcel_centroids(dlabel, N, SURF_L, SURF_R)
    # Try geodesic (will fallback to None) — now returns cache dict
    geo_cache=geodesic_distances_via_wb(SURF_L, SURF_R, centroids, centroid_vertex)
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
        # For each edge, compute distances
        for pheno, h2 in zip(sub["Pheno"].astype(str), sub[col].astype(float)):
            s=str(pheno)
            k=int(s[1:]) if s.startswith("o") else int(s)
            i,j=np.triu_indices(N,1)
            ii=int(i[k]); jj=int(j[k])
            # centroids are 1-indexed parcels
            c1=centroids.get(ii+1, np.array([np.nan,np.nan,np.nan]))
            c2=centroids.get(jj+1, np.array([np.nan,np.nan,np.nan]))
            euclid=np.linalg.norm(c1-c2) if np.all(np.isfinite(c1)) and np.all(np.isfinite(c2)) else np.nan
            # geodesic via wb_command if available and same hemisphere
            geod=np.nan
            if geo_cache is not None:
                try:
                    geod=_geodesic_for_pair(atlas, ii, jj, centroids, centroid_vertex, geo_cache)
                except:
                    geod=np.nan
            same=networks[ii]==networks[jj] and networks[ii] not in ("NA","SUB","???")
            rows.append({"atlas": atlas, "method": method, "pheno": pheno, "i": ii, "j": jj,
                         "h2": float(h2), "euclidean_dist": float(euclid) if np.isfinite(euclid) else np.nan,
                         "geodesic_dist": float(geod) if np.isfinite(geod) else np.nan,
                         "same_network": bool(same), "network_i": str(networks[ii]), "network_j": str(networks[jj])})
        print(f"[{atlas} {method}] collected {len([r for r in rows if r['atlas']==atlas and r['method']==method])} edges")

df=pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV} with {len(df)} rows")
# Simple scatter plots per atlas/method
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for atlas in ["gordon","probaConns"]:
    for method in ["Twin","AdjHE-RE"]:
        sub=df[(df["atlas"]==atlas) & (df["method"]==method)]
        if sub.empty:
            continue
        fig, ax=plt.subplots(figsize=(5,4))
        # Sample for speed if large
        plot_df=sub.dropna(subset=["euclidean_dist","h2"])
        if len(plot_df)>20000:
            plot_df=plot_df.sample(20000, random_state=0)
        ax.scatter(plot_df["euclidean_dist"], plot_df["h2"], s=1, alpha=0.2, c="#1f77b4")
        ax.set_xlabel("Euclidean centroid distance (mm)")
        ax.set_ylabel(r"Heritability ($h^2$)")
        ax.set_ylim(0,1)
        # No title per publication
        fig.tight_layout(pad=0.5)
        out=PLOT_DIR / f"fc_distance_{atlas}_{method.replace('-','')}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(f"  wrote {out}")
        # Also geodesic if available (currently all NaN, skip)
        if sub["geodesic_dist"].notna().sum()>10:
            fig, ax=plt.subplots(figsize=(5,4))
            plot_df=sub.dropna(subset=["geodesic_dist","h2"])
            if len(plot_df)>20000:
                plot_df=plot_df.sample(20000, random_state=0)
            ax.scatter(plot_df["geodesic_dist"], plot_df["h2"], s=1, alpha=0.2, c="#ff7f0e")
            ax.set_xlabel("Geodesic distance (mm via wb_command)")
            ax.set_ylabel(r"Heritability ($h^2$)")
            ax.set_ylim(0,1)
            fig.tight_layout(pad=0.5)
            out=PLOT_DIR / f"fc_distance_geodesic_{atlas}_{method.replace('-','')}.png"
            fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.05)
            plt.close(fig)
            print(f"  wrote {out}")

print("Done")
