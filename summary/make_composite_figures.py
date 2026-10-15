#!/usr/bin/env python3
"""Composite figures: brain surface + Manhattan, and barplot + two surfaces.

Two layouts, both driven from summary/plotstyle.py so text sizing matches the
standalone panels:

  FC  (gordon, probaConns)   2 rows, full width
                             top    : cortical surface for one method
                             bottom : Sys-Sys Manhattan for the same method

  SA                         3 panels
                             left, spanning both rows : w_total horizontal barplot
                             top right               : AdjHE-RE w_total surface
                             bottom right            : Twin w_total surface

What is re-drawn and what is placed as an image
------------------------------------------------
The Manhattan half is drawn as vector text directly into the composite, so its
axis labels stay crisp and its font sizes are governed by plotstyle. The
surfaces and the barplot are placed as images: rendering them needs nilearn and
brainTemplates (the dlabel and inflated surfaces), which are not available on a
laptop, so the PNGs produced on the HPC are reused. If you re-run this on the
HPC where those templates exist, re-running 02a-02d first will refresh the
images automatically because they are read from disk by path.

Usage
-----
    python summary/make_composite_figures.py              # every layout
    python summary/make_composite_figures.py --list
    python summary/make_composite_figures.py fc           # FC layouts only
    python summary/make_composite_figures.py sa           # SA layout only
    python summary/make_composite_figures.py --from-pull  # read images from
                                                         # results/hpc_pull
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plotstyle import (DPI, FC_COMPOSITE_FIGSIZE, FC_COMPOSITE_ROWS,
                       SA_COMPOSITE_BAR_WIDTH, SA_COMPOSITE_FIGSIZE,
                       SA_COMPOSITE_ROWS, FS_PANEL_LABEL, savefig_publication)

SUMMARY = Path(__file__).resolve().parent
ROOT = SUMMARY.parent
RESULTS = ROOT / "results" / "summary"
PLOTS = RESULTS / "plots"
BRAIN = RESULTS / "brain"
PULL = ROOT / "results" / "hpc_pull"
OUT = RESULTS / "composite"

sys.path.insert(0, str(SUMMARY))
import importlib
_man = importlib.import_module("05_manhattan")

ATLAS_SET_N = {"gordon": 352, "probaConns": 80, "SA": 17}
METHODS = [("Twin", "Twin"), ("AdjHE-RE", "AdjHERE")]
ATLAS_TITLE = {"gordon": "Gordon", "probaConns": "ProbaConns", "SA": "SA"}

SURFACE_PATHS = {
    ("gordon", "AdjHERE"): "gordon/gordon_AdjHERE_surface.png",
    ("gordon", "Twin"): "gordon/gordon_twin_surface.png",
    ("probaConns", "AdjHERE"): "probaConns/probaConns_AdjHERE_surface.png",
    ("probaConns", "Twin"): "probaConns/probaConns_twin_surface.png",
    ("SA", "AdjHERE"): "SA_w_total/SA_w_total_AdjHE_surface.png",
    ("SA", "Twin"): "SA_w_total/SA_w_total_twin_surface.png",
}

SA_BARPLOT = "SA_w_total_horizontal_barplot.png"


def _resolve(rel, use_pull):
    """Find an image in the results tree, or fall back to the pulled snapshot."""
    for base in ([PULL] if use_pull else []) + [BRAIN.parent, PULL]:
        p = base / rel if base is not PULL else PULL / "brain" / rel
        if p.exists():
            return p
    return None


def _resolve_plot(rel, use_pull):
    for base in ([PULL / "plots"] if use_pull else []) + [PLOTS, PULL / "plots"]:
        p = base / rel
        if p.exists():
            return p
    return None


def _trim_whitespace(img, tol=0.02, pad_frac=0.004):
    """Crop near-white margins from a placed image.

    The HPC surface and barplot PNGs are written with bbox_inches="tight" but
    still carry wide white margins, so when they are dropped into a subplot they
    shrink to fill only the middle of the panel. Cropping first lets each panel
    use the space it was allotted.
    """
    import numpy as np
    a = np.asarray(img)
    if a.ndim == 3 and a.shape[2] >= 3:
        rgb = a[..., :3]
        if a.dtype == np.uint8:
            rgb = rgb / 255.0
        content = rgb.min(axis=2) < (1.0 - tol)
    else:
        return img
    if not content.any():
        return img
    rows = np.where(content.any(axis=1))[0]
    cols = np.where(content.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return img
    h, w = a.shape[0], a.shape[1]
    py = int(pad_frac * rows.size) or 1
    px = int(pad_frac * cols.size) or 1
    r0, r1 = max(0, rows[0] - py), min(h, rows[-1] + py + 1)
    c0, c1 = max(0, cols[0] - px), min(w, cols[-1] + px + 1)
    return a[r0:r1, c0:c1]


def _place_image(ax, path, label=None, corner_label=None):
    """Show a PNG filling the axes, preserving aspect, with no frame."""
    img = _trim_whitespace(mpimg.imread(str(path)))
    ax.imshow(img, interpolation="antialiased")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if corner_label:
        ax.text(0.005, 0.985, corner_label, transform=ax.transAxes, ha="left",
                va="top", fontsize=FS_PANEL_LABEL, fontweight="bold")
    if label:
        ax.text(0.5, 1.015, label, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=FS_PANEL_LABEL)


def _wide_axes(fig, row, top, height, left=0.0, width=1.0):
    ax = fig.add_axes([left, top, width, height])
    return ax


def build_fc(wide, atlas, method, tag, use_pull):
    """Two rows: surface on top, Manhattan below."""
    surf = _resolve(SURFACE_PATHS[(atlas, tag)], use_pull)
    if surf is None:
        print(f"  SKIP {atlas} {method}: surface PNG not found")
        return None

    h_surf, h_man = FC_COMPOSITE_ROWS
    fig = plt.figure(figsize=FC_COMPOSITE_FIGSIZE)

    ax_top = _wide_axes(fig, 0, 1.0 - h_surf, h_surf)
    _place_image(ax_top, surf, corner_label="A")

    ax_bot = _wide_axes(fig, 1, 0.0, h_man)
    df = _man.build_long_for_set(wide, atlas, ATLAS_SET_N[atlas],
                                 _man.METHOD_SPECS[atlas][0][1] if method == "Twin"
                                 else _man.METHOD_SPECS[atlas][1][1])
    if df.empty:
        plt.close(fig)
        print(f"  SKIP {atlas} {method}: no rows")
        return None

    dfs = {"Twin": _man.build_long_for_set(wide, atlas, ATLAS_SET_N[atlas],
                                           _man.METHOD_SPECS[atlas][0][1]),
           "AdjHE-RE": _man.build_long_for_set(wide, atlas, ATLAS_SET_N[atlas],
                                               _man.METHOD_SPECS[atlas][1][1])}
    shared = _man._get_shared_order(dfs)
    prepared, layout, large_order, largest_20 = _man.prepare_manhattan_df(df, shared)
    _man.draw_manhattan_panel(ax_bot, prepared, layout, large_order, largest_20,
                              atlas=atlas)
    # reserve room for the rotated Sys-Sys labels under the bottom axes
    fig.subplots_adjust(bottom=0.13, top=0.97, left=0.05, right=0.995)

    out = OUT / f"fc_{atlas}_{tag}.png"
    savefig_publication(fig, out)
    plt.close(fig)
    print(f"  wrote {out.name}  (xlim={layout['xlim']})")
    return out


def build_sa(wide, use_pull):
    """Barplot at left spanning both rows; two w_total surfaces at right."""
    bar = _resolve_plot(SA_BARPLOT, use_pull)
    s_adjhe = _resolve(SURFACE_PATHS[("SA", "AdjHERE")], use_pull)
    s_twin = _resolve(SURFACE_PATHS[("SA", "Twin")], use_pull)
    missing = [n for n, p in (("barplot", bar), ("AdjHE-RE surface", s_adjhe),
                              ("Twin surface", s_twin)) if p is None]
    if missing:
        print(f"  SKIP SA: missing {', '.join(missing)}")
        return None

    r_top_h, r_bot_h = SA_COMPOSITE_ROWS
    fig = plt.figure(figsize=SA_COMPOSITE_FIGSIZE)
    bw = SA_COMPOSITE_BAR_WIDTH
    right_w = 1.0 - bw
    gap = 0.012

    ax_bar = fig.add_axes([0.0, 0.0, bw, 1.0])
    _place_image(ax_bar, bar, corner_label="A")

    ax_top = fig.add_axes([bw + gap, 1.0 - r_top_h - gap, right_w, r_top_h])
    _place_image(ax_top, s_adjhe, corner_label="B")

    ax_bot = fig.add_axes([bw + gap, 0.0, right_w, r_bot_h - gap])
    _place_image(ax_bot, s_twin, corner_label="C")

    fig.subplots_adjust(left=0.0, right=1.0, top=0.97, bottom=0.03)

    out = OUT / "sa_w_total_composite.png"
    savefig_publication(fig, out)
    plt.close(fig)
    print(f"  wrote {out.name}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="?", default="all",
                    choices=["all", "fc", "sa"], help="which layout(s) to build")
    ap.add_argument("--list", action="store_true", help="list what will be built")
    ap.add_argument("--from-pull", action="store_true",
                    help="prefer images under results/hpc_pull/")
    args = ap.parse_args()

    wide_path = RESULTS / "mash_twin_wide.csv"
    if not wide_path.exists():
        wide_path = PULL / "mash_twin_wide.csv"
    if not wide_path.exists():
        print("mash_twin_wide.csv not found in results/summary or results/hpc_pull")
        return 2

    if args.list:
        print("FC composites (2 rows):")
        for atlas in ("gordon", "probaConns"):
            for method, tag in METHODS:
                print(f"  fc_{atlas}_{tag}.png   surface + manhattan  [{method}]")
        print("SA composite (3 panels):")
        print("  sa_w_total_composite.png   barplot | AdjHE-RE + Twin surfaces")
        print(f"\nwide: {wide_path}")
        return 0

    import pandas as pd
    wide = pd.read_csv(wide_path)
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"wide: {wide_path}")
    if args.which in ("all", "fc"):
        for atlas in ("gordon", "probaConns"):
            for method, tag in METHODS:
                build_fc(wide, atlas, method, tag, args.from_pull)
    if args.which in ("all", "sa"):
        build_sa(wide, args.from_pull)
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())