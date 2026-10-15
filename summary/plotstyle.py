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

import matplotlib.pyplot as plt

# --- publication output -----------------------------------------------------
DPI = 300
BBOX = "tight"
PAD_INCHES = 0.05

# --- palette ----------------------------------------------------------------
ALT_COLORS = ["#1f77b4", "#ff7f0e"]
SMALL_COLOR = "grey"
CMAP = "coolwarm"

# --- font scale (points) ----------------------------------------------------
# One scale for the whole report. A standalone panel and the same panel inside a
# composite must use the same sizes, or the two disagree once they sit side by
# side in a manuscript.
FS_PANEL_LABEL = 11   # short caption above a composite sub-panel
FS_TITLE = 12
FS_AXIS_LABEL = 11
FS_AXIS_LABEL_X = 11
FS_TICK = 10
# Rotated 45 degree Sys-Sys labels are the densest text in the report, so they
# sit one step below the tick size rather than far below it -- at 7 pt they were
# effectively unreadable once the panel was scaled into a composite figure.
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
# Composite layouts. FC stacks a surface above a Manhattan; SA puts a barplot
# at the left spanning both rows with two surfaces stacked at the right.
FC_COMPOSITE_FIGSIZE = (13.0, 9.6)
FC_COMPOSITE_ROWS = (0.44, 0.56)   # (surface, manhattan) height fractions
SA_COMPOSITE_FIGSIZE = (13.6, 8.4)
# The barplot is naturally wide-ish (roughly 1.4:1) while its column is tall and
# narrow, so give it more width than an even split or it letterboxes badly.
SA_COMPOSITE_BAR_WIDTH = 0.40
SA_COMPOSITE_ROWS = (0.5, 0.5)


def savefig_publication(fig, path, dpi=DPI, bbox=BBOX, pad_inches=PAD_INCHES):
    """Write a figure at publication resolution and close it."""
    fig.savefig(path, dpi=dpi, bbox_inches=bbox, pad_inches=pad_inches)
    plt.close(fig)
    return path


def compute_axis_layout(connection_order, largest_20, small_order,
                        n_rows=None, compress_factor=COMPRESS_FACTOR,
                        min_small_fraction=MIN_SMALL_FRACTION):
    """Return the x-axis geometry for one Manhattan panel.

    Parameters
    ----------
    connection_order : list
        Sys-Sys groups in the order they are drawn, large groups first.
    largest_20 : list
        The subset of ``connection_order`` drawn at full resolution.
    small_order : list
        The subset drawn in the compressed tail.
    n_rows : int, optional
        Number of rows the calling method actually retained. Used only to
        decide whether the compressed tail is worth drawing at all.

    Returns
    -------
    dict with:
        conn_pos : dict   connection -> position in connection_order
        divider  : float  where the compressed tail starts, or None
        xlim     : tuple  (lo, hi) covering the full ordering
        compressed : bool
    """
    order = list(connection_order)
    n_total = len(order)
    conn_pos = {conn: i for i, conn in enumerate(order)}

    if n_total == 0:
        return {"conn_pos": {}, "divider": None, "xlim": (0.0, 1.0),
                "compressed": False}

    small_set = set(small_order)
    large_n = len([c for c in order if c not in small_set])

    # Only compress when the tail is a meaningful fraction of the panel.
    small_n = n_total - large_n
    compress = bool(small_n) and small_n > 0
    if n_rows is not None and n_rows > 0:
        compress = compress and (small_n / float(n_total)) > min_small_fraction

    if compress:
        # Divider sits at the first position occupied by a small group. Use the
        # full ordering, not the caller's rows, so every panel for an atlas
        # agrees on where the tail begins.
        small_positions = [conn_pos[c] for c in small_order if c in conn_pos]
        divider = float(min(small_positions)) if small_positions else None
        if divider is None:
            compress = False

    if compress and divider is not None:
        # Tail runs from divider to the end of the ordering, squashed by
        # compress_factor. Large groups keep 1:1 spacing up to the divider.
        tail_len = float(n_total - 1 - divider)
        xmax = divider + tail_len / compress_factor
    else:
        divider = None
        xmax = float(n_total - 1)

    pad = (xmax * XLIM_PAD_FRACTION) if xmax > 0 else 1.0
    xlim = (-pad, xmax + pad)

    return {"conn_pos": conn_pos, "divider": divider, "xlim": xlim,
            "compressed": compress}


def apply_plot_index(df, layout, small_order, compress_factor=COMPRESS_FACTOR):
    """Add a ``plot_index`` column to ``df`` using a computed layout.

    ``df`` must already carry a ``connection`` column and a ``conn_pos``
    column (build it from ``layout['conn_pos']``).
    """
    import numpy as np
    import pandas as pd

    conn_pos = layout["conn_pos"]
    divider = layout["divider"]

    df = df.copy()
    df["conn_pos"] = df["connection"].map(conn_pos)
    if df["conn_pos"].isna().any():
        raise ValueError(
            "connection_order does not cover every row: "
            f"{sorted(df.loc[df['conn_pos'].isna(), 'connection'].unique())[:5]}"
        )

    if divider is not None:
        mask_small = df["connection"].isin(set(small_order))
        df["plot_index"] = np.where(
            mask_small,
            (df["conn_pos"].astype(float) - float(divider)) / compress_factor + float(divider),
            df["conn_pos"].astype(float),
        )
    else:
        df["plot_index"] = df["conn_pos"].astype(float)

    if pd.isna(df["plot_index"]).any():
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