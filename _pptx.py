"""Refresh the images embedded in figures.pptx with the current figures."""
import shutil
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

missing_files = [p for p in MAP.values() if not __import__("os").path.exists(p)]
if missing_files:
    raise SystemExit("source figures missing:\n  " + "\n  ".join(missing_files))

shutil.copy2("figures.pptx", "figures.pptx.bak")

src = zipfile.ZipFile("figures.pptx")
names = set(src.namelist())
changed = []
out = {}
for name in src.namelist():
    data = src.read(name)
    if name in MAP:
        with open(MAP[name], "rb") as fh:
            data = fh.read()
        changed.append(name)
    out[name] = data
src.close()

with zipfile.ZipFile("figures.pptx", "w", zipfile.ZIP_DEFLATED) as dst:
    for name, data in out.items():
        dst.writestr(name, data)

print("replaced %d/%d media files" % (len(changed), len(MAP)))
not_in_deck = set(MAP) - names
if not_in_deck:
    print("not present in deck (harmless):", sorted(not_in_deck))