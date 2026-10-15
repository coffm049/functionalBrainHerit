#!/usr/bin/env python3
"""Export the figure slides of a PPTX to individual PNGs sized for submission.

Resolutions, not scale factors
------------------------------
The slide is 13.33 x 7.5 in, but no journal page is that wide -- MIT Press
text width is about 7 in, and Imaging Neuroscience asks for figures supplied as
separate files. So the meaningful resolution is the one the figure has *at its
final printed size*, not on the slide.

Each slide is therefore exported at a target width in inches and a target dpi,
and the pixel width is computed as width_in * dpi. Exporting wider does not add
detail: it only upscales the embedded rasters further. The embedded images here
are 300 dpi PNGs, so a slide exported at 3.5 in and 600 dpi holds them well
above their native detail while keeping the file small.

Text-only slides are skipped: they carry no figure, and Imaging Neuroscience
wants figure files, not slides.

Usage
-----
    python summary/export_pptx_figures.py                       # defaults
    python summary/export_pptx_figures.py --width-in 7 --dpi 300
    python summary/export_pptx_figures.py --include-text-only

Requires PowerPoint (COM). The pptx must not be open in PowerPoint.
"""
from __future__ import annotations

import argparse
import re
import struct
import zipfile
from pathlib import Path

import xml.etree.ElementTree as ET

EMU = 914400.0
ROOT = Path(__file__).resolve().parents[1]

NS = {
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def slide_geometry(pptx: Path):
    """Return (slide_w_in, slide_h_in, [(slide_no, n_images, n_text_shapes)])."""
    with zipfile.ZipFile(pptx) as z:
        pres = ET.fromstring(z.read("ppt/presentation.xml"))
        sz = pres.find("p:sldSz", NS)
        sw, sh = int(sz.get("cx")) / EMU, int(sz.get("cy")) / EMU

        names = [n for n in z.namelist()
                 if re.match(r"ppt/slides/slide\d+\.xml$", n)]
        names.sort(key=lambda s: int(re.findall(r"\d+", s)[0]))

        out = []
        for n in names:
            num = int(re.findall(r"\d+", n)[0])
            root = ET.fromstring(z.read(n))
            n_img = sum(1 for _ in root.iter("{%s}pic" % NS["p"]))
            n_txt = sum(1 for e in root.iter()
                        if e.tag.endswith("}sp") and e.tag != "{%s}pic" % NS["p"])
            out.append((num, n_img, n_txt))
        return sw, sh, out


def picture_report(pptx: Path, target_width_in: float):
    """Per-picture printed width and effective dpi at the target print size.

    A picture placed ``pin`` inches wide on a slide that is printed
    ``target_width_in`` wide ends up ``(pin / slide_width) * target_width_in``
    inches wide in the paper. Its real detail is its own pixel count, so the
    effective resolution is ``pixels / printed_width_in``.

    This is the number that matters. Resolution "on the slide" is misleading:
    the slide is 13.33 in wide, which is wider than any journal page, so a
    figure can look under-resolved on the slide and still be far above 300 dpi
    once printed at 3.5 in.
    """
    with zipfile.ZipFile(pptx) as z:
        pres = ET.fromstring(z.read("ppt/presentation.xml"))
        sz = pres.find("p:sldSz", NS)
        slide_w_in = int(sz.get("cx")) / EMU

        dims = {}
        for n in z.namelist():
            if n.startswith("ppt/media/"):
                b = z.read(n)
                if b[:8] == b"\x89PNG\r\n\x1a\n":
                    w, h = struct.unpack(">II", b[16:24])
                    dims[n.split("/")[-1]] = (w, h)

        names = [n for n in z.namelist()
                 if re.match(r"ppt/slides/slide\d+\.xml$", n)]
        names.sort(key=lambda s: int(re.findall(r"\d+", s)[0]))

        rows = []
        for n in names:
            num = int(re.findall(r"\d+", n)[0])
            rels = ET.fromstring(
                z.read("ppt/slides/_rels/slide%d.xml.rels" % num))
            rmap = {c.get("Id"): c.get("Target").split("/")[-1] for c in rels}
            root = ET.fromstring(z.read(n))
            for pic in root.iter("{%s}pic" % NS["p"]):
                blip = pic.find(".//a:blip", NS)
                ext = pic.find(".//a:xfrm/a:ext", NS)
                if blip is None or ext is None:
                    continue
                tgt = rmap.get(blip.get("{%s}embed" % NS["r"]), "?")
                px_wide = int(ext.get("cx")) / EMU
                printed = (px_wide / slide_w_in) * target_width_in
                px = dims.get(tgt, (0, 0))[0]
                rows.append({
                    "slide": num, "image": tgt,
                    "on_slide_in": px_wide, "printed_in": printed,
                    "pixels": px,
                    "effective_dpi": (px / printed) if printed else 0.0,
                })
        return slide_w_in, rows


def export(pptx: Path, outdir: Path, width_in: float, dpi: int,
           include_text_only: bool = False) -> int:
    import win32com.client as win32
    from win32com.client import constants

    outdir.mkdir(parents=True, exist_ok=True)
    px = int(round(width_in * dpi))

    app = win32.Dispatch("PowerPoint.Application")
    pres = None
    written = 0
    try:
        # Open(FileName, ReadOnly, Untitled, WithWindow)
        pres = app.Presentations.Open(str(pptx.resolve()), True, False, False)
        slide_w = pres.PageSetup.SlideWidth / 72.0
        slide_h = pres.PageSetup.SlideHeight / 72.0
        height_px = int(round(px * slide_h / slide_w))

        for i in range(1, pres.Slides.Count + 1):
            slide = pres.Slides(i)
            n_shapes = slide.Shapes.Count
            n_pics = sum(1 for j in range(1, n_shapes + 1)
                         if slide.Shapes(j).Type == 13)  # msoPicture
            if n_pics == 0 and not include_text_only:
                print(f"  slide {i}: skipped (no images)")
                continue

            dst = outdir / f"slide{i:02d}.png"
            # Slide.Export(FileName, FilterName, ScaleWidth, ScaleHeight)
            slide.Export(str(dst), "PNG", px, height_px)
            size = dst.stat().st_size / 1e6
            flag = "  <-- OVER 10Mb" if size > 10 else ""
            print(f"  slide {i}: {n_pics} image(s) -> {dst.name} "
                  f"({px} x {height_px} px, {size:.2f} Mb){flag}")
            written += 1
    finally:
        if pres is not None:
            pres.Close()
        app.Quit()
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pptx", default="figures.pptx", type=Path)
    ap.add_argument("--outdir", default="results/summary/slides", type=Path)
    ap.add_argument("--width-in", type=float, default=3.5,
                    help="final printed width of the figure in inches")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--include-text-only", action="store_true")
    args = ap.parse_args()

    pptx = args.pptx if args.pptx.is_absolute() else ROOT / args.pptx
    outdir = args.outdir if args.outdir.is_absolute() else ROOT / args.outdir

    if not pptx.exists():
        print(f"not found: {pptx}")
        return 2

    sw, sh, slides = slide_geometry(pptx)
    px = int(round(args.width_in * args.dpi))

    print(f"pptx        : {pptx.name}")
    print(f"slide size  : {sw:.2f} x {sh:.2f} in")
    print(f"export size : {args.width_in} in at {args.dpi} dpi "
          f"-> {px} x {int(round(px * sh / sw))} px")
    print()
    print("slides:")
    for num, n_img, n_txt in slides:
        kind = "FIGURE" if n_img else "text-only"
        print(f"  {num}: {kind} ({n_img} image, {n_txt} text shape)")
    print()

    # Effective resolution at the final printed size -- the number that decides
    # whether the figure meets a journal's floor.
    slide_w, rows = picture_report(pptx, args.width_in)
    if rows:
        print(f"effective resolution of each embedded figure at {args.width_in} in wide:")
        for r in rows:
            print(f"  slide {r['slide']}  {r['image']:14} "
                  f"{r['pixels']:5d} px -> prints {r['printed_in']:.2f} in "
                  f"= {r['effective_dpi']:6.0f} dpi")
        worst = min(r["effective_dpi"] for r in rows)
        verdict = "meets 300 dpi" if worst >= 300 else "BELOW 300 dpi"
        print(f"  worst = {worst:.0f} dpi -> {verdict}")
        print()

    try:
        n = export(pptx, outdir, args.width_in, args.dpi, args.include_text_only)
    except ImportError:
        print("pywin32 is required: pip install pywin32")
        return 2
    print(f"\nexported {n} slide(s) to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())