#!/usr/bin/env python3
"""Trim white margins from generated PNGs.

The brain-surface figures come out of nilearn's 3D axes with wide blank bands
above and below the hemispheres. `bbox_inches="tight"` does not remove them: a
3D axes reports a bounding box larger than the mesh it actually draws, so the
extra space survives into the file.

This crops each image to the bounding box of its non-white pixels. It is pure
pixel work, so it needs neither nilearn nor brainTemplates and can be run on a
pulled copy of the results.

Usage
-----
    python summary/trim_png_margins.py                       # brain figs in results/
    python summary/trim_png_margins.py --dir results/hpc_pull/brain
    python summary/trim_png_margins.py --glob "results/summary/brain/**/*_surface.png"
    python summary/trim_png_margins.py --inplace --backup    # overwrite, keep .orig

By default the trimmed image is written alongside the original as
<name>_trimmed.png, so nothing is destroyed until you have looked at it.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_GLOBS = [
    "results/summary/brain/**/*_surface.png",
    "results/summary/brain/**/*_circular.png",
]


def trim_array(img, tol=0.02, pad_frac=0.004):
    """Return (cropped_image, (w0,h0), (w1,h1)).

    A pixel counts as content when its darkest channel is below 1 - tol, so
    light-grey cortex (about 0.8) is kept while a pure white margin (1.0) is
    dropped. tol is deliberately loose so that very pale values, and any
    anti-aliased edge, survive.
    """
    a = np.asarray(img)
    if a.ndim == 3 and a.shape[2] >= 3:
        rgb = a[..., :3]
        if a.dtype == np.uint8:
            rgb = rgb / 255.0
        content = rgb.min(axis=2) < (1.0 - tol)
    else:
        return a, a.shape[1], a.shape[0]

    if not content.any():
        return a, a.shape[1], a.shape[0]

    rows = np.where(content.any(axis=1))[0]
    cols = np.where(content.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return a, a.shape[1], a.shape[0]

    h, w = a.shape[0], a.shape[1]
    # keep a sliver of margin so the outermost anti-aliased pixel is not clipped
    py = max(1, int(pad_frac * rows.size))
    px = max(1, int(pad_frac * cols.size))
    r0, r1 = max(0, rows[0] - py), min(h, rows[-1] + py + 1)
    c0, c1 = max(0, cols[0] - px), min(w, cols[-1] + px + 1)
    return a[r0:r1, c0:c1], w, h


def trim_file(path: Path, inplace: bool = False, backup: bool = False,
              tol: float = 0.02) -> bool:
    img = mpimg.imread(str(path))
    if img.ndim == 2:
        print(f"  skip {path.name}: greyscale")
        return False
    out, w0, h0 = trim_array(img, tol=tol)
    h1, w1 = out.shape[0], out.shape[1]
    if (w1, h1) == (w0, h0):
        print(f"  no margin to trim: {path.name} ({w0}x{h0})")
        return False

    if inplace:
        if backup and not path.with_suffix(path.suffix + ".orig").exists():
            shutil.copy2(path, path.with_suffix(path.suffix + ".orig"))
        fig = plt.figure(figsize=(w1 / 300, h1 / 300), dpi=300)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(out, interpolation="antialiased")
        ax.set_axis_off()
        fig.savefig(path, dpi=300)
        plt.close(fig)
        dst = path
    else:
        dst = path.with_name(path.stem + "_trimmed.png")
        fig = plt.figure(figsize=(w1 / 300, h1 / 300), dpi=300)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(out, interpolation="antialiased")
        ax.set_axis_off()
        fig.savefig(dst, dpi=300)
        plt.close(fig)

    print(f"  {path.name}: {w0}x{h0} -> {w1}x{h1}  "
          f"({100.0 * (1 - (w1 * h1) / float(w0 * h0)):.1f}% of pixels removed)"
          + ("" if inplace else f"  -> {dst.name}"))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", action="append", default=[],
                    help="directory to scan recursively (repeatable)")
    ap.add_argument("--glob", action="append", default=[],
                    help="glob relative to repo root (repeatable)")
    ap.add_argument("--inplace", action="store_true", help="overwrite originals")
    ap.add_argument("--backup", action="store_true",
                    help="with --inplace, keep a .orig copy")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="white threshold; higher trims less aggressively")
    args = ap.parse_args()

    globs = list(args.glob) or list(DEFAULT_GLOBS)
    files = []
    for d in args.dir:
        files += sorted(Path(d).rglob("*.png"))
    for g in globs:
        files += sorted(ROOT.glob(g))

    # de-duplicate, keep order
    seen, uniq = set(), []
    for f in files:
        s = str(f)
        if s in seen:
            continue
        # Skip our own side-output and backups so a second run is idempotent and
        # does not trim already-trimmed copies of themselves.
        if f.name.endswith("_trimmed.png") or f.name.endswith(".orig"):
            continue
        seen.add(s)
        uniq.append(f)

    if not uniq:
        print("no PNGs matched")
        return 0

    print(f"scanning {len(uniq)} file(s)")
    n = 0
    for f in uniq:
        if trim_file(f, inplace=args.inplace, backup=args.backup, tol=args.tol):
            n += 1
    print(f"\n{n} of {len(uniq)} trimmed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())