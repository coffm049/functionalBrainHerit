"""Place the SA barplot on slide 8 correctly.

Identifies the barplot picture by the image it references rather than by its
crop marker, so it can be run repeatedly and repairs the position whether or not
an earlier pass already removed the crop.

Slide is 13.333 x 7.5 in. The surface panels occupy x >= 5.9 in, so the barplot
is fitted to the width to their left, vertically centred on the row of panels.
"""
import re
import shutil
import struct
import zipfile

EMU = 914400
SLIDE_W, SLIDE_H = 13.3333, 7.5
SURFACE_LEFT = 5.9          # surfaces begin here
MARGIN = 0.10
BARPLOT_IMAGE = "image9.png"


def emu(inches):
    return int(round(inches * EMU))


def inches(v):
    return v / EMU


src = "figures.pptx"
z = zipfile.ZipFile(src)
xml = z.read("ppt/slides/slide8.xml").decode("utf-8")
rels = z.read("ppt/slides/_rels/slide8.xml.rels").decode("utf-8")

# map rId -> media target
rid_to_img = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))

target = None
target_idx = None
for i, m in enumerate(re.finditer(r"<p:pic>.*?</p:pic>", xml, re.S)):
    pic = m.group(0)
    rid = re.search(r'r:embed="(rId\d+)"', pic)
    if not rid:
        continue
    t = rid_to_img[rid.group(1)].replace("../", "").split("/")[-1]
    if t == BARPLOT_IMAGE:
        target, target_idx = pic, i
        break

if target is None:
    raise SystemExit(f"{BARPLOT_IMAGE} not found on slide 8")

img = z.read("ppt/media/" + BARPLOT_IMAGE)
iw, ih = struct.unpack(">II", img[16:24])
aspect = iw / ih

# Vertical span of the two surface rows, to centre the barplot against.
# Skip by index: comparing the picture strings with `is` fails, because they come
# from two separate regex searches.
exts = []
for i, m in enumerate(re.finditer(r"<p:pic>.*?</p:pic>", xml, re.S)):
    if i == target_idx:
        continue
    pic = m.group(0)
    e = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"/>', pic)
    o = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"/>', pic)
    if e and o:
        exts.append((int(o.group(2)), int(e.group(2))))
top = min(y for y, _ in exts)
bottom = max(y + h for y, h in exts)
band_centre = (top + bottom) / 2

top_in = inches(top)
bottom_in = inches(bottom)
band_centre_in = (top_in + bottom_in) / 2

# Everything below is in inches; EMU conversion happens only at write time.
# (top/bottom arrive in EMU from the XML, so they are converted here once.)
w = SURFACE_LEFT - 2 * MARGIN
h = w / aspect
x = MARGIN
y = band_centre_in - h / 2

if y < MARGIN:
    # would fall off the top; shrink to the available height instead
    h = top_in - 2 * MARGIN
    w = h * aspect
    x = MARGIN
    y = MARGIN

if x + w > SURFACE_LEFT + 1e-6:
    raise SystemExit(f"barplot would overlap the surface panels "
                     f"(right edge {x + w:.2f} in vs {SURFACE_LEFT} in)")
if y + h > SLIDE_H - MARGIN + 1e-6:
    raise SystemExit(f"barplot would fall off the bottom of the slide "
                     f"(bottom {y + h:.2f} in vs {SLIDE_H} in)")

print(f"image            : {BARPLOT_IMAGE}  {iw} x {ih}  aspect {aspect:.4f}")
print(f"surface band     : y {top_in:.2f} .. {bottom_in:.2f} in"
      f"  (centre {band_centre_in:.2f})")
print(f"barplot frame    : x {x:.2f}  y {y:.2f}  {w:.2f} x {h:.2f} in"
      f"  (aspect {w / h:.4f})")

new_pic = re.sub(r'<a:srcRect[^/>]*/>', "", target)          # ensure no crop
old_o = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"/>', target)
old_e = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"/>', target)
print(f"was              : x {inches(int(old_o.group(1))):.2f}"
      f"  y {inches(int(old_o.group(2))):.2f}"
      f"  {inches(int(old_e.group(1))):.2f} x {inches(int(old_e.group(2))):.2f} in")

new_pic = new_pic.replace(old_o.group(0), f'<a:off x="{emu(x)}" y="{emu(y)}"/>', 1)
new_pic = new_pic.replace(old_e.group(0), f'<a:ext cx="{emu(w)}" cy="{emu(h)}"/>', 1)

new_xml = xml.replace(target, new_pic, 1)
names = z.namelist()
out = {n: z.read(n) for n in names}
z.close()

shutil.copy2(src, "figures.pptx.bak")
with zipfile.ZipFile(src, "w", zipfile.ZIP_DEFLATED) as dst:
    for n in names:
        dst.writestr(n, new_xml.encode("utf-8") if n == "ppt/slides/slide8.xml" else out[n])

chk = zipfile.ZipFile(src)
assert chk.testzip() is None
x8 = chk.read("ppt/slides/slide8.xml").decode("utf-8")
for m in re.finditer(r"<p:pic>.*?</p:pic>", x8, re.S):
    p = m.group(0)
    rid = re.search(r'r:embed="(rId\d+)"', p)
    if rid and rid_to_img[rid.group(1)].replace("../", "").split("/")[-1] == BARPLOT_IMAGE:
        assert "<a:srcRect" not in p, "crop still present"
        e = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"/>', p)
        o = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"/>', p)
        fx, fy, fw, fh = int(o.group(1)), int(o.group(2)), int(e.group(1)), int(e.group(2))
        assert abs(fw / fh - aspect) < 1e-3, "aspect not restored"
        assert 0 <= fx and fy >= 0 and fx + fw <= emu(SURFACE_LEFT) and fy + fh <= emu(SLIDE_H)
        print(f"verified         : x {inches(fx):.2f} y {inches(fy):.2f}"
              f" {inches(fw):.2f} x {inches(fh):.2f} in, no crop, on-slide")
chk.close()