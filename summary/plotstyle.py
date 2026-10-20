#!/usr/bin/env python3
"""Shared plotting style and axis-layout helpers.

Two things live here:

1. Publication defaults (dpi, bbox, padding, palette, font sizes) so every
   figure in the report is written the same way and looks the same size.
2. compute_axis_layout() -- the Sys-Sys x-axis geometry for the Manhattan
   panels, as a pure function with no plotting and no dataframe access.

compute_axis_layout exists because the x-axis span was previously computed
inside the plotting function from data-dependent fallbacks. That made it
impossible to reason about without rendering a figure, and it produced three
successive bugs: panels whose x-axis spanned only their own retained rows,
over-compression on smaller atlases, and a NameError from a partially renamed
variable. The function below takes the ordering as an argument and returns the
geometry, so the same ordering always yields the same axis regardless of which
rows a given method happens to retain.

Usage:
    from plotstyle import compute_axis_layout, savefig_publication
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt

# Imaging Neuroscience requires Arial or Helvetica in figures. Importing this
# module sets that, so every script that imports plotstyle inherits it.
from figure_font import RESOLVED_FONT, apply_font  # noqa: F401

# --- publication output -----------------------------------------------------
DPI = 300
BBOX = "tight"
PAD_INCHES = 0.05

# Every figure is written in each of these formats. PNG keeps the Quarto
# catalog and figures.pptx working; PDF is the vector file the journal wants
# (figure_font sets pdf.fonttype=42, so text is embedded rather than outlined).
# Override for a quick run with FBH_FIGURE_FORMATS="png".
FIGURE_FORMATS = tuple(
    ext.strip().lstrip(".").lower()
    for ext in os.environ.get("FBH_FIGURE_FORMATS", "png,pdf").split(",")
    if ext.strip()
)

# --- palette ----------------------------------------------------------------
ALT_COLORS = ["#1f77b4", "#ff7f0e"]
SMALL_COLOR = "grey"
CMAP = "coolwarm"

# --- font scale (points) ----------------------------------------------------
# One scale for the whole report. The rotated 45 degree Sys-Sys labels are the
# densest text here, so they sit only one step below the tick size; at 7 pt they
# were effectively unreadable.
FS_TITLE = 12
FS_AXIS_LABEL = 11
FS_AXIS_LABEL_X = 11
FS_TICK = 10
FS_TICK_ROTATED = 9

# --- Manhattan geometry -----------------------------------------------------
COMPRESS_FACTOR = 100.0
# The compressed tail is only worth drawing when it is a real fraction of the
# panel. Below this the divider lands too early and squashes the plot.
MIN_SMALL_FRACTION = 0.20
XLIM_PAD_FRACTION = 0.015

# --- figure sizes (inches) --------------------------------------------------
# Wide because a Sys-Sys x-axis carries 20 rotated labels across the full width.
MANHATTAN_FIGSIZE = (14.0, 4.5)
# Bottom margin reserved for those rotated labels.
MANHATTAN_BOTTOM = 0.28


def figure_paths(path, formats=None) -> list[Path]:
    """`path` written in every publication format.

    The extension on `path` names the primary format and stays first, so a
    caller that does `stats = str(path).replace('.png', '_stats.csv')` keeps
    working while the PDF is written next to it.
    """
    fmts = FIGURE_FORMATS if formats is None else tuple(formats)
    stem = Path(path).with_suffix("")
    primary = Path(path).suffix.lstrip(".").lower() or "png"
    ordered = [primary] + [f for f in fmts if f != primary]
    return [stem.with_suffix("." + f) for f in ordered]


def savefig_publication(fig, path, dpi=DPI, bbox=BBOX, pad_inches=PAD_INCHES,
                        formats=None):
    """Write a figure in every publication format, then close it.

    Returns the primary path (the one passed in), so existing callers that
    build a stats filename from it are unaffected.
    """
    for out in figure_paths(path, formats):
        fig.savefig(out, dpi=dpi, bbox_inches=bbox, pad_inches=pad_inches)
    plt.close(fig)
    return path


def compute_axis_layout(connection_order, largest_20, small_order,
                        group_counts=None, compress_factor=COMPRESS_FACTOR,
                        min_small_fraction=MIN_SMALL_FRACTION):
    """Return the x-axis geometry for one Manhattan panel.

    Each Sys-Sys group is allotted a slot whose width is proportional to how
    many edges it holds, and each edge inside the group gets its own distinct
    position across that slot. That is what makes the number of points in a
    group legible: a group with 400 edges is visibly wide, a group with 6 is a
    thin spike. Collapsing a group to a single x position would throw that
    information away.

    Group widths come from ``group_counts``, which must be the POOLED edge
    count across every method for the atlas. Using one method's counts would
    give that method's axis and no other, which is what previously made the
    probaConns Twin panel span a different range from its AdjHE-RE neighbour.

    Parameters
    ----------
    connection_order : list
        Sys-Sys groups in the order they are drawn, large groups first.
    largest_20 : list
        The subset of ``connection_order`` drawn at full resolution.
    small_order : list
        The subset drawn in the compressed tail.
    group_counts : dict, optional
        connection -> pooled edge count. Defaults to one unit per group.

    Returns
    -------
    dict with:
        slot_start : dict  connection -> left edge of its slot (edge units)
        slot_width : dict  connection -> slot width
        divider    : float where the compressed tail starts, or None
        xlim       : tuple (lo, hi) covering the full ordering
        total      : float total un-compressed width
        compressed : bool
    """
    order = list(connection_order)
    counts = group_counts or {}
    slot_start, slot_width = {}, {}
    cursor = 0.0
    for conn in order:
        w = float(counts.get(conn, 0) or 0)
        if w <= 0:
            w = 1.0
        slot_start[conn] = cursor
        slot_width[conn] = w
        cursor += w
    total = cursor

    if not order:
        return {"slot_start": {}, "slot_width": {}, "divider": None,
                "xlim": (0.0, 1.0), "total": 0.0, "compressed": False}

    small_set = set(small_order)
    small_width = sum(slot_width[c] for c in order if c in small_set)

    # Only compress when the tail is a meaningful fraction of the panel.
    compress = small_width > 0 and (small_width / total) > min_small_fraction

    divider = None
    if compress:
        first_small = next((c for c in order if c in small_set), None)
        divider = slot_start[first_small] if first_small is not None else None
        if divider is None:
            compress = False

    if compress and divider is not None:
        xmax = divider + (total - divider) / compress_factor
    else:
        divider = None
        xmax = total

    pad = (xmax * XLIM_PAD_FRACTION) if xmax > 0 else 1.0
    return {"slot_start": slot_start, "slot_width": slot_width,
            "divider": divider, "xlim": (-pad, xmax + pad),
            "total": total, "compressed": compress}


def apply_plot_index(df, layout, small_order, compress_factor=COMPRESS_FACTOR):
    """Add a ``plot_index`` column: one distinct position per edge.

    ``df`` must carry a ``connection`` column. Edges are ranked within their own
    group and laid across that group's slot, so no two edges share an x value.
    """
    import numpy as np

    slot_start = layout["slot_start"]
    slot_width = layout["slot_width"]
    divider = layout["divider"]

    df = df.copy()
    unknown = set(df["connection"]) - set(slot_start)
    if unknown:
        raise ValueError(f"connection_order does not cover: {sorted(unknown)[:5]}")

    df = df.sort_values("connection", kind="stable")
    # rank within group, in row order, so the layout is reproducible
    df["_rank"] = df.groupby("connection", observed=True).cumcount()
    # how many edges THIS method retains in each group
    df["_own"] = df.groupby("connection", observed=True)["connection"].transform("size")
    df["_start"] = df["connection"].map(slot_start).astype(float)
    df["_width"] = df["connection"].map(slot_width).astype(float)

    # Spread this method's own edges across the whole of the group's pooled
    # slot: first edge at slot start + half a step, last at the slot end. A
    # method that retains fewer edges than the pooled count therefore fills its
    # slot at lower density instead of leaving a blank gap before the next pair.
    own = df["_own"].where(df["_own"] > 0, 1.0).astype(float)
    width = df["_width"].where(df["_width"] > 0, 1.0).astype(float)
    df["plot_index"] = df["_start"] + (df["_rank"] + 0.5) / own * width

    if divider is not None:
        mask_small = df["connection"].isin(set(small_order))
        df.loc[mask_small, "plot_index"] = (
            divider + (df.loc[mask_small, "plot_index"] - divider) / compress_factor)

    df = df.drop(columns=["_rank", "_own", "_start", "_width"])
    if np.isnan(df["plot_index"]).any():
        raise ValueError("plot_index contains NaN after layout")
    return df


def check_layout_agrees(layout_a, layout_b, label_a="A", label_b="B"):
    """Assert two layouts for the same ordering agree on the x-axis.

    Used as a plain-assert guard (no pytest) to catch the failure where one
    method's panel spanned only its own retained rows.
    """
    if tuple(layout_a["xlim"]) != tuple(layout_b["xlim"]):
        raise AssertionError(
            f"x-axis mismatch between {label_a} and {label_b}: "
            f"{tuple(layout_a['xlim'])} vs {tuple(layout_b['xlim'])}"
        )
    return True