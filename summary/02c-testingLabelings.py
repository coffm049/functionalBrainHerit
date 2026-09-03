#!/usr/bin/env python3
"""
Gordon 333 / Proba 80 ΓÇö Basic atlas labeling tests (no h2).

Goal: Revisit the Gordon 333 network labeling from the ground up to find the
simple mistake that causes the mismatched legend/cortex colors reported:

  VIS pink vs light green occipital, DMN peach vs bright orange angular,
  MTL light green vs occipital, AUD blue vs brown, Tpole dark brown vs red PON,
  plus duplicate legend entries (DMN/SMl dark green, PON/Sal light green) and
  SUB/NA overlapping background.

This script *does not* use h2 values. It simply visualizes the dlabel itself
via different labeling approaches so we can compare which one paints the cortex
correctly. Outputs go to results/summary/brain/gordon/test_*.png and
results/summary/brain/probaConns/test_*.png for side-by-side inspection.

Approaches tested for Gordon (N=333 cortical, N=352 with subcortical):
  A. Explicit LUT: Gordon333_LUT.txt  ROI_ID -> Network_Name (user suggested)
     lut = pd.read_csv("Gordon333_LUT.txt", sep=r"\s+", header=None)
     id_to_net = dict(zip(lut[0], lut[2]))  # col 0 = ROI_ID, col 2 = Network_Name
  B. gordon_modules.csv (352 rows, region 1..14) via labelDictionary
     as in 02a-manhattanOrder.py: labelDict {1:DMN,2:VIS,3:FP,4:DAN,5:VAN,6:SAL,7:CO,8:SMD,9:SML,10:AUD,11:Tpole,12:MTL,13:PMN,14:PON}
  C. CIFTI label table + shortDicionary.csv (14 rows) via _network_of + name_to_short
     as in previous summary/02a_brain_gordon.py
  D. Direct dlabel parcel ID (no network, just ROI) ΓÇö sanity check for boundaries

For Proba (N=80, 80 parcels):
  A. CIFTI label table via _network_of (DMN_LABEL_1 -> DMN etc.)
  B. GRP1_template_parcel.csv (80 rows, region 1..14) via labelDictionary
  C. Explicit swap Sal<->SMl (user observed brown insula should be Sal)

All approaches use the same surface mapping (iter_structures, 32492) so the
*only* difference is parcel->network, isolating the labeling bug.

Run:
  .venv/bin/python summary/02c-testingLabelings.py
Outputs:
  results/summary/brain/gordon/test_A_LUT.png
  results/summary/brain/gordon/test_B_gordon_modules.png
  results/summary/brain/gordon/test_C_cifti_shortDict.png
  results/summary/brain/gordon/test_D_direct.png
  results/summary/brain/probaConns/test_A_cifti.png
  results/summary/brain/probaConns/test_B_grp1.png
  results/summary/brain/probaConns/test_C_swapped.png
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import nilearn.plotting as niplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = Path(__file__).resolve().parents[1]
def _resolve(rel, abs_path):
    p = Path(rel)
    return p if p.exists() else Path(abs_path)

def _find_gordon(subcortical=True):
    cands = [
        ROOT.parent / "brainTemplates" / ("Gordon.subcortical.32k_fs_LR.dlabel.nii" if subcortical else "Gordon333.32k_fs_LR.dlabel.nii"),
        Path(r"C:\Users\coffm049\brainTemplates\Gordon.subcortical.32k_fs_LR.dlabel.nii") if subcortical else Path(r"C:\Users\coffm049\brainTemplates\Gordon333.32k_fs_LR.dlabel.nii"),
        Path("/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii") if subcortical else Path("/users/4/coffm049/papers/brainTemplates/Gordon333.32k_fs_LR.dlabel.nii"),
        _resolve(ROOT.parent / "brainTemplates" / ("Gordon.subcortical.32k_fs_LR.dlabel.nii" if subcortical else "Gordon333.32k_fs_LR.dlabel.nii"),
                 "/users/4/coffm049/papers/brainTemplates/Gordon.subcortical.32k_fs_LR.dlabel.nii" if subcortical else "/users/4/coffm049/papers/brainTemplates/Gordon333.32k_fs_LR.dlabel.nii"),
    ]
    for c in cands:
        if c.exists() and c.stat().st_size > 1000:
            return c
    return cands[0]

DLABEL_GORDON_SUB = _find_gordon(subcortical=True)
DLABEL_GORDON_333 = _find_gordon(subcortical=False)
DLABEL_GORDON_NET = _resolve(ROOT.parent / "brainTemplates" / "Gordon.networks.32k_fs_LR.dlabel.nii",
                             "/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii")
DLABEL_PROBA = None
for cand in [
    ROOT.parent / "brainTemplates" / "abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii",
    Path(r"C:\Users\coffm049\brainTemplates\abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii"),
    Path("/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.75.dlabel.nii"),
    ROOT.parent / "brainTemplates" / "combined_clusters_thresh0.75.dlabel.nii",
    _resolve(ROOT / "brainTemplate" / "combined_clusters_thresh0.75.dlabel.nii", "/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.75.dlabel.nii"),
    Path(r"C:\Users\coffm049\brainTemplates\combined_clusters_thresh0.75.dlabel.nii"),
]:
    if cand.exists():
        # Check that it has 80 parcels (not 63)
        try:
            img = nib.load(str(cand))
            if len(img.header.get_axis(0).get_element(0)[1]) - 1 == 80:
                DLABEL_PROBA = cand
                break
        except Exception:
            pass
else:
    DLABEL_PROBA = cand

def _find_surf(hemi):
    cands = [
        ROOT.parent / "brainTemplates" / f"Conte69.{hemi}.inflated.32k_fs_LR.surf.gii",
        Path(rf"C:\Users\coffm049\brainTemplates\Conte69.{hemi}.inflated.32k_fs_LR.surf.gii"),
        Path(f"/users/4/coffm049/papers/brainTemplates/Conte69.{hemi}.inflated.32k_fs_LR.surf.gii"),
    ]
    for c in cands:
        if c.exists():
            return c
    return cands[0]
SURF_L = _find_surf("L")
SURF_R = _find_surf("R")

OUT_GORDON = ROOT / "results/summary/brain/gordon"
OUT_PROBA = ROOT / "results/summary/brain/probaConns"
OUT_GORDON.mkdir(parents=True, exist_ok=True)
OUT_PROBA.mkdir(parents=True, exist_ok=True)

# Helpers copied from 02a/02b for consistency
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
    for attr in ("labeltable", "labels"):
        lt = getattr(img, attr, None)
        if lt is not None:
            d = getattr(lt, "labels", lt)
            if isinstance(d, dict) and len(d) > 1:
                return d
            if isinstance(d, list) and len(d) > 1:
                return {i: v for i, v in enumerate(d) if v is not None}
    return {}

def _network_of(name):
    if not name:
        return "NA"
    s = str(name)
    s = re.sub(r"^\d{1,3}_[LR]_", "", s)
    _map = {"Default": "DMN", "Auditory": "Aud", "CinguloOperc": "CO", "FrontoParietal": "FP",
            "Salience": "Sal", "VentralAttn": "VAN", "DorsalAttn": "DAN", "Visual": "Vis",
            "SMhand": "SMd", "SMmouth": "SMl", "RetrosplenialTemporal": "Tpole"}
    s = _map.get(s, s)
    return s.split("_")[0] if "_" in s else s

def dlabel_to_textures(label_data, full, dlabel_path, surf_l=SURF_L, surf_r=SURF_R):
    """Map per-greyordinate 'full' values to Conte69 textures via iter_structures."""
    img = nib.load(str(dlabel_path))
    ax = img.header.get_axis(1)
    va = np.asarray(ax.vertex, dtype=int)
    tex_l = tex_r = None
    for structure, sl, _ in ax.iter_structures():
        s = str(structure).upper()
        if sl.stop is None:
            sl = slice(sl.start, len(va))
        blk = va[sl]
        seg = full[sl]
        if "CORTEX_LEFT" in s:
            tex_l = np.full(32492, np.nan)
            valid = blk >= 0
            # need to handle seg that may be int network ids
            tex_l[blk[valid]] = seg[valid]
        elif "CORTEX_RIGHT" in s:
            tex_r = np.full(32492, np.nan)
            valid = blk >= 0
            tex_r[blk[valid]] = seg[valid]
    return tex_l, tex_r

def plot_network_textures(tex_l, tex_r, net_names, title, out_png, dlabel_for_colors=None):
    """Plot two hemispheres with categorical network colors."""
    # Build palette: try dlabel true colors first, else tab20 sorted
    base = plt.get_cmap("tab20").colors
    if dlabel_for_colors is not None and Path(dlabel_for_colors).exists():
        try:
            img = nib.load(str(dlabel_for_colors))
            labels = _cifti_labels(img)
            net_to_rgba = {}
            for p, lab in labels.items():
                if p == 0:
                    continue
                raw = lab[0] if isinstance(lab, (tuple, list)) else getattr(lab, "label", str(lab))
                net = _network_of(raw)
                # For Gordon 333, the network name via gordon_modules may not match raw; try to map via shortDict
                # Instead, just collect unique nets that are in net_names
                for target in net_names:
                    if target.lower() == net.lower():
                        if target not in net_to_rgba:
                            rgba = lab[1] if isinstance(lab, (tuple, list)) and len(lab) > 1 else (0.5, 0.5, 0.5, 1.0)
                            if isinstance(rgba, (list, tuple)) and len(rgba) >= 3:
                                net_to_rgba[target] = tuple(float(x) for x in rgba[:3])
                        break
            if len(net_to_rgba) == len(net_names):
                color_of = {n: net_to_rgba[n] for n in net_names}
                cmap = mcolors.ListedColormap([color_of[n] for n in net_names])
            else:
                raise ValueError("incomplete dlabel color mapping")
        except Exception as e:
            print(f"  palette fallback (dlabel colors incomplete: {e}) -> tab20 sorted")
            sorted_nets = sorted(net_names, key=lambda x: x.lower())
            tmp = {net: base[i % len(base)] for i, net in enumerate(sorted_nets)}
            if "Sal" in tmp:
                tmp["Sal"] = base[5]
            if "SMl" in tmp:
                tmp["SMl"] = base[4]
            color_of = {n: tmp[n] for n in net_names}
            cmap = mcolors.ListedColormap([color_of[n] for n in net_names])
    else:
        sorted_nets = sorted(net_names, key=lambda x: x.lower())
        tmp = {net: base[i % len(base)] for i, net in enumerate(sorted_nets)}
        if "Sal" in tmp:
            tmp["Sal"] = base[5]
        if "SMl" in tmp:
            tmp["SMl"] = base[4]
        color_of = {n: tmp[n] for n in net_names}
        cmap = mcolors.ListedColormap([color_of[n] for n in net_names])

    fig = plt.figure(figsize=(6*2, 5))
    for i, (hemi, surf, val) in enumerate([("left", SURF_L, tex_l), ("right", SURF_R, tex_r)], start=1):
        if val is None:
            continue
        ax = fig.add_subplot(1, 2, i, projection="3d")
        niplot.plot_surf_stat_map(str(surf), val, hemi=hemi, axes=ax, figure=fig,
                                  cmap=cmap, colorbar=False, threshold=None,
                                  vmin=0, vmax=max(len(net_names)-1,1), title=f"{title} ({hemi})")
    handles = [plt.Line2D([0],[0], marker="o", linestyle="", color=color_of[n], label=n) for n in net_names]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_png} with {len(net_names)} nets: {net_names}")

def approach_A_lut_gordon(dlabel_path):
    """Explicit LUT: Gordon333_LUT.txt ROI_ID -> Network_Name, as user suggested."""
    lut_candidates = [
        ROOT / "Gordon333_LUT.txt",
        ROOT.parent / "brainTemplates" / "Gordon333_LUT.txt",
        Path(r"C:\Users\coffm049\brainTemplates\Gordon333_LUT.txt"),
        Path("/users/4/coffm049/papers/brainTemplates/Gordon333_LUT.txt"),
        Path("Gordon333_LUT.txt"),
    ]
    # Also try to find any LUT-like file in repo
    for cand in lut_candidates:
        if cand.exists():
            try:
                lut = pd.read_csv(cand, sep=r"\s+", header=None)
                # Expect col 0 = ROI_ID, col 2 = Network_Name
                id_to_net = dict(zip(lut[0].astype(int), lut[2].astype(str)))
                print(f"  LUT found {cand} with {len(id_to_net)} entries, sample {list(id_to_net.items())[:3]}")
                return id_to_net
            except Exception as e:
                print(f"  LUT {cand} read failed: {e}")
    # Fallback: use dlabel label table itself as LUT (ROI_ID -> network via _network_of)
    print("  No LUT file found, falling back to dlabel label table as LUT")
    img = nib.load(str(dlabel_path))
    labels = _cifti_labels(img)
    id_to_net = {}
    for k, v in labels.items():
        if k == 0:
            continue
        raw = v[0] if isinstance(v, (tuple, list)) else getattr(v, "label", str(v))
        id_to_net[k] = _network_of(raw)
    return id_to_net

def approach_B_gordon_modules():
    """gordon_modules.csv (352 rows) via labelDictionary 1..14."""
    cands = [
        ROOT.parent / "brainTemplates" / "gordon_modules.csv",
        Path(r"C:\Users\coffm049\brainTemplates\gordon_modules.csv"),
        Path("/users/4/coffm049/papers/brainTemplates/gordon_modules.csv"),
    ]
    for cand in cands:
        if cand.exists():
            gm = pd.read_csv(cand)
            col = gm.columns[0]
            vals = gm[col].astype(int).tolist()
            labelDict = {1:"DMN",2:"VIS",3:"FP",4:"DAN",5:"VAN",6:"SAL",7:"CO",8:"SMD",9:"SML",10:"AUD",11:"Tpole",12:"MTL",13:"PMN",14:"PON"}
            # Map to shortnames as in 02a (VIS->Vis etc.)
            canon = {"VIS":"Vis","SAL":"Sal","SMD":"SMd","SML":"SMl","AUD":"Aud"}
            nets = [canon.get(labelDict.get(v, "NA"), labelDict.get(v, "NA")) for v in vals]
            id_to_net = {i+1: net for i, net in enumerate(nets)}
            print(f"  gordon_modules {cand} -> {len(nets)} nets, unique {sorted(set(nets))}")
            return id_to_net
    return {}

def approach_C_cifti_shortdict(dlabel_path):
    """CIFTI label table + shortDicionary.csv (previous summary/02a)."""
    short_cands = [
        ROOT.parent / "brainTemplates" / "shortDicionary.csv",
        Path(r"C:\Users\coffm049\brainTemplates\shortDicionary.csv"),
    ]
    for cand in short_cands:
        if cand.exists():
            df = pd.read_csv(cand)
            name_to_short = dict(zip(df["name"].astype(str), df["shortname"].astype(str)))
            img = nib.load(str(dlabel_path))
            labels = _cifti_labels(img)
            id_to_net = {}
            for k, v in labels.items():
                if k == 0:
                    continue
                raw = v[0] if isinstance(v, (tuple, list)) else getattr(v, "label", str(v))
                net = _network_of(raw)
                short = name_to_short.get(net, net)
                if short == net:
                    for kk, vv in name_to_short.items():
                        if kk.lower() == net.lower():
                            short = vv
                            break
                id_to_net[k] = short
            return id_to_net
    return {}

if __name__ == "__main__":
    print("=== Gordon 333 testing ===")
    # Use the 333 cortical dlabel if valid, else subcortical
    dlabel_333 = DLABEL_GORDON_333 if DLABEL_GORDON_333.exists() and DLABEL_GORDON_333.stat().st_size > 1000 else DLABEL_GORDON_SUB
    print(f"Using Gordon dlabel: {dlabel_333} (333 exists: {DLABEL_GORDON_333.exists()}, size {DLABEL_GORDON_333.stat().st_size if DLABEL_GORDON_333.exists() else 0})")
    img = nib.load(str(dlabel_333))
    label_data = np.asarray(img.get_fdata()).astype(int).ravel()
    print(f"  label_data shape {label_data.shape}, unique {sorted(np.unique(label_data))[:10]}...{sorted(np.unique(label_data))[-5:]}")

    # A: LUT
    print("\nApproach A: LUT explicit ROI_ID -> Network")
    id_to_net_A = approach_A_lut_gordon(dlabel_333)
    # Build per-vertex full array: map each greyordinate's label_data value via id_to_net
    # For visualization, we need to map network names to integer ids for colormap
    nets_A = sorted(set(id_to_net_A.values()) - {"NA","SUB","???"})
    net_to_idx_A = {n:i for i,n in enumerate(nets_A)}
    full_A = np.full(label_data.shape, np.nan)
    for roi_id, net in id_to_net_A.items():
        if net in net_to_idx_A:
            full_A[label_data == roi_id] = net_to_idx_A[net]
    tex_l, tex_r = dlabel_to_textures(label_data, full_A, dlabel_333)
    plot_network_textures(tex_l, tex_r, nets_A, "Gordon 333 Test A: LUT explicit", OUT_GORDON / "test_A_LUT.png", dlabel_333)

    # B: gordon_modules.csv
    print("\nApproach B: gordon_modules.csv + labelDict")
    id_to_net_B = approach_B_gordon_modules()
    if id_to_net_B:
        nets_B = sorted(set(id_to_net_B.values()) - {"NA","SUB","???"})
        net_to_idx_B = {n:i for i,n in enumerate(nets_B)}
        full_B = np.full(label_data.shape, np.nan)
        for roi_id, net in id_to_net_B.items():
            if net in net_to_idx_B:
                full_B[label_data == roi_id] = net_to_idx_B[net]
        tex_l, tex_r = dlabel_to_textures(label_data, full_B, dlabel_333)
        plot_network_textures(tex_l, tex_r, nets_B, "Gordon 333 Test B: gordon_modules.csv", OUT_GORDON / "test_B_gordon_modules.png", dlabel_333)
    else:
        print("  B skipped (no gordon_modules.csv)")

    # C: CIFTI + shortDicionary
    print("\nApproach C: CIFTI label table + shortDicionary")
    id_to_net_C = approach_C_cifti_shortdict(dlabel_333)
    if id_to_net_C:
        nets_C = sorted(set(id_to_net_C.values()) - {"NA","SUB","???"})
        net_to_idx_C = {n:i for i,n in enumerate(nets_C)}
        full_C = np.full(label_data.shape, np.nan)
        for roi_id, net in id_to_net_C.items():
            if net in net_to_idx_C:
                full_C[label_data == roi_id] = net_to_idx_C[net]
        tex_l, tex_r = dlabel_to_textures(label_data, full_C, dlabel_333)
        plot_network_textures(tex_l, tex_r, nets_C, "Gordon 333 Test C: CIFTI+shortDict", OUT_GORDON / "test_C_cifti_shortDict.png", dlabel_333)

    # D: Direct parcel ID (boundaries only, no network)
    print("\nApproach D: Direct parcel ID (boundaries)")
    # Just visualize label_data % 14 to see boundaries
    full_D = label_data.astype(float)
    full_D[label_data == 0] = np.nan
    # Map to network via modulo for visualization
    tex_l, tex_r = dlabel_to_textures(label_data, full_D % 14, dlabel_333)
    # Use 14 nets for palette
    nets_D = [f"ROI_{i}" for i in range(14)]
    plot_network_textures(tex_l, tex_r, nets_D[:14], "Gordon 333 Test D: Direct parcel", OUT_GORDON / "test_D_direct.png", dlabel_333)

    print("\n=== Proba 80 testing ===")
    if DLABEL_PROBA and Path(DLABEL_PROBA).exists():
        print(f"Using Proba dlabel: {DLABEL_PROBA}")
        img_p = nib.load(str(DLABEL_PROBA))
        label_data_p = np.asarray(img_p.get_fdata()).astype(int).ravel()
        # Proba A: CIFTI
        id_to_net_pA = {}
        labels_p = _cifti_labels(img_p)
        for k,v in labels_p.items():
            if k==0: continue
            raw = v[0] if isinstance(v,(tuple,list)) else getattr(v,"label",str(v))
            id_to_net_pA[k] = _network_of(raw)
        nets_pA = sorted(set(id_to_net_pA.values()) - {"NA","SUB","???"})
        net_to_idx_pA = {n:i for i,n in enumerate(nets_pA)}
        full_pA = np.full(label_data_p.shape, np.nan)
        for roi_id, net in id_to_net_pA.items():
            if net in net_to_idx_pA:
                full_pA[label_data_p == roi_id] = net_to_idx_pA[net]
        tex_l, tex_r = dlabel_to_textures(label_data_p, full_pA, DLABEL_PROBA)
        plot_network_textures(tex_l, tex_r, nets_pA, "Proba 80 Test A: CIFTI", OUT_PROBA / "test_A_cifti.png", DLABEL_PROBA)

        # Proba B: GRP1
        grp = pd.read_csv(ROOT.parent / "brainTemplates" / "GRP1_template_parcel.csv") if (ROOT.parent / "brainTemplates" / "GRP1_template_parcel.csv").exists() else pd.read_csv(r"C:\Users\coffm049\brainTemplates\GRP1_template_parcel.csv")
        labelDict = {1:"DMN",2:"VIS",3:"FP",4:"DAN",5:"VAN",6:"SAL",7:"CO",8:"SMD",9:"SML",10:"AUD",11:"Tpole",12:"MTL",13:"PMN",14:"PON"}
        canon = {"VIS":"Vis","SAL":"Sal","SMD":"SMd","SML":"SMl","AUD":"Aud"}
        vals = grp[grp.columns[0]].astype(int).tolist()
        nets_b = [canon.get(labelDict.get(v,"NA"), labelDict.get(v,"NA")) for v in vals]
        id_to_net_pB = {i+1: net for i, net in enumerate(nets_b)}
        nets_pB = sorted(set(nets_b) - {"NA","SUB","???"})
        net_to_idx_pB = {n:i for i,n in enumerate(nets_pB)}
        full_pB = np.full(label_data_p.shape, np.nan)
        for roi_id, net in id_to_net_pB.items():
            if net in net_to_idx_pB:
                full_pB[label_data_p == roi_id] = net_to_idx_pB[net]
        tex_l, tex_r = dlabel_to_textures(label_data_p, full_pB, DLABEL_PROBA)
        plot_network_textures(tex_l, tex_r, nets_pB, "Proba 80 Test B: GRP1", OUT_PROBA / "test_B_grp1.png", DLABEL_PROBA)

        # Proba C: Swapped Sal<->SMl
        id_to_net_pC = {k: ("SMl" if v=="Sal" else "Sal" if v=="SMl" else v) for k,v in id_to_net_pA.items()}
        nets_pC = sorted(set(id_to_net_pC.values()) - {"NA","SUB","???"})
        net_to_idx_pC = {n:i for i,n in enumerate(nets_pC)}
        full_pC = np.full(label_data_p.shape, np.nan)
        for roi_id, net in id_to_net_pC.items():
            if net in net_to_idx_pC:
                full_pC[label_data_p == roi_id] = net_to_idx_pC[net]
        tex_l, tex_r = dlabel_to_textures(label_data_p, full_pC, DLABEL_PROBA)
        plot_network_textures(tex_l, tex_r, nets_pC, "Proba 80 Test C: CIFTI Swapped Sal/SMl", OUT_PROBA / "test_C_swapped.png", DLABEL_PROBA)
    else:
        print(f"Proba dlabel not found: {DLABEL_PROBA}")

    print("Done. Compare PNGs in results/summary/brain/gordon/test_*.png and probaConns/test_*.png")
