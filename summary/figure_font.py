#!/usr/bin/env python3
"""Set Arial as the figure font for every matplotlib figure in the report.

Imaging Neuroscience's guide for authors asks that figures use Arial or
Helvetica. Matplotlib's default is DejaVu Sans, which does not meet that rule.
Arial is present on this machine, so it is selected here rather than left to the
DejaVu fallback.

Imported for its side effect by the plotting scripts. Call apply_font() when a
font must be applied to an already-imported pyplot.
"""
from __future__ import annotations

import matplotlib

# Imaging Neuroscience: "Please use Arial or Helvetica fonts within figures."
FIGURE_FONT_SANS = ["Arial", "Helvetica", "DejaVu Sans"]


def apply_font():
    """Point matplotlib at Arial, keeping Helvetica/DejaVu as fallbacks.

    Returns the family name matplotlib actually resolved to, which may be a
    fallback if Arial is missing -- callers can warn on that.
    """
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = FIGURE_FONT_SANS
    matplotlib.rcParams["pdf.fonttype"] = 42      # embed a real font in PDFs
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["svg.fonttype"] = "none"  # keep text as text in SVG

    # Work out which family actually wins. findfont() takes a concrete family
    # name, so it must be handed each candidate in turn -- passing the alias
    # "sans-serif" makes it raise ParseException rather than resolve.
    import matplotlib.font_manager as fm
    cache = {f.name for f in fm.fontManager.ttflist}
    for name in FIGURE_FONT_SANS:
        if name in cache:
            return name
    return "unknown"


RESOLVED_FONT = apply_font()

if RESOLVED_FONT != "Arial":
    import warnings
    warnings.warn(
        f"Imaging Neuroscience asks for Arial/Helvetica in figures, but "
        f"matplotlib resolved the sans-serif family to {RESOLVED_FONT!r}. "
        f"Arial is not installed on this machine.",
        RuntimeWarning,
        stacklevel=2,
    )