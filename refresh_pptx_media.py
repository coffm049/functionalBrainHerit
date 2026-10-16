"""Refresh the images embedded in figures.pptx with the current figures.

figures.pptx is a hand-assembled deck: the images are placed and sized in
PowerPoint, but the PNGs themselves come from results/. After a figure rerun the
deck still holds the old PNGs, so they have to be swapped in place. Swap by media
path so PowerPoint's placement and extents are preserved.
"""
import os
import shutil
import sys
import zipfile

MAP = {
    "ppt/media/image1.png": "results/summary/brain/gordon/gordon_AdjHERE_surface.png",
    "ppt/media/image2.png": "results/summary/plots/manhattan_gordon_AdjHERE.png",
    "ppt/media/image3.png": "results/summary/brain/gordon/gordon_twin_surface.png",
    "ppt/media/image4.png": "results/summary/plots/manhattan_gordon_Twin.png",
    "ppt/media/image5.png": "results/summary/brain/probaConns/probaConns_AdjHERE_surface.png",
    "ppt/media/image6.png": "results/summary/plots/manhattan_probaConns_AdjHERE.png",
    "ppt/media/image7.png": "results/summary/brain/probaConns/probaConns_twin_surface.png",
    "ppt/media/image8.png": "results/summary/plots/manhattan_probaConns_Twin.png",
    "ppt/media/image9.png": "results/summary/plots/SA_w_total_horizontal_barplot.png",
    "ppt/media/image10.png": "results/summary/brain/SA_w_total/SA_w_total_twin_surface.png",
    "ppt/media/image11.png": "results/summary/brain/SA_w_total/SA_w_total_AdjHE_surface.png",
}

missing = [p for p in MAP.values() if not os.path.exists(p)]
if missing:
    sys.exit("source figures missing:\n  " + "\n  ".join(missing))

shutil.copy2("figures.pptx", "figures.pptx.bak")

src = zipfile.ZipFile("figures.pptx")
names = set(src.namelist())
out = {n: (open(MAP[n], "rb").read() if n in MAP else src.read(n)) for n in src.namelist()}
src.close()

tmp = "figures.pptx.new"
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
    for n, data in out.items():
        dst.writestr(n, data)

# Verify the rebuilt archive before it replaces the deck.
chk = zipfile.ZipFile(tmp)
bad = chk.testzip()
assert bad is None, f"corrupt entry in rebuilt archive: {bad}"
replaced = sum(1 for n in MAP if n in names)
assert all(n in chk.namelist() for n in names), "rebuilt archive lost entries"
chk.close()

shutil.move(tmp, "figures.pptx")
print(f"replaced {replaced}/{len(MAP)} media files; archive verified")

extra = set(MAP) - names
if extra:
    print("not present in deck (harmless):", sorted(extra))