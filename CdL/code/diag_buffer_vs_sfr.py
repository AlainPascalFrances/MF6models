"""
Clarify: STREAM_BUFFER (refinement corridor) vs SFR cells (channel centreline).

  refined cells  = cells whose centroid lies within STREAM_BUFFER (60 m) of the
                   stream -> these are made SMALL by Triangle (cell SIZE control).
  SFR cells      = cells the stream CENTRELINE crosses -> the actual reaches.

SFR is a ~1-cell-wide thread down the spine of the (120 m-wide) corridor, so
refined-cell-count >> SFR-cell-count.  Also writes a zoom figure.
"""
import config
import pickle
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, Patch
from matplotlib.collections import PatchCollection
import geopandas as gpd
from shapely.ops import unary_union, linemerge
from shapely.geometry import LineString, MultiPolygon, Point
from shapely.validation import make_valid
from pyproj import CRS
from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect

WS      = str(config.MODEL)
WS_PKL  = WS + r"\voronoi_grid.pkl"
SFR_PKG = WS + r"\cdl_gwf.sfr"
GPKG    = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
TARGET  = CRS.from_epsg(3763)
STREAM_BUFFER = 60.0      # must match the model

with open(WS_PKL, "rb") as f:
    gridprops_vg, _ = pickle.load(f)
vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
xc = np.array([vgrid.xcellcenters[i] for i in range(ncpl)])
yc = np.array([vgrid.ycellcenters[i] for i in range(ncpl)])
try:
    ix = GridIntersect(vgrid, method="vertex")
except TypeError:
    ix = GridIntersect(vgrid)

# model SFR cells (node = packagedata col 3, 1-based)
model_sfr = set()
with open(SFR_PKG) as f:
    blk = False
    for line in f:
        s = line.strip(); low = s.lower()
        if low.startswith("begin packagedata"): blk = True; continue
        if low.startswith("end packagedata"): break
        if blk and s and not s.startswith("#"):
            model_sfr.add(int(s.split()[2]) - 1)

ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws_poly = unary_union(ws.geometry.apply(lambda g: make_valid(g)).buffer(0).values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(TARGET)
merged = linemerge(unary_union(streams.geometry.apply(make_valid).values))
stream_union = (merged if merged.geom_type == "LineString"
                else unary_union(list(merged.geoms))).intersection(ws_poly)

# refined cells = centroid within STREAM_BUFFER of the stream
buf = stream_union.buffer(STREAM_BUFFER)
refined = {i for i in range(ncpl) if buf.contains(Point(xc[i], yc[i]))}

print(f"grid cells (ncpl)                         : {ncpl}")
print(f"refined cells (centroid within {STREAM_BUFFER:.0f} m)   : {len(refined)}")
print(f"SFR cells (stream centreline crosses)     : {len(model_sfr)}")
print(f"  -> SFR is {100*len(model_sfr)/len(refined):.0f}% of the refined corridor cells; "
      f"the other {len(refined)-len(model_sfr)} refined cells flank the channel (UZF/DRN, not SFR).")

# ---- zoom figure on a stream section ----
cx, cy = stream_union.interpolate(0.5, normalized=True).coords[0]
HW = 220.0
x0, x1, y0, y1 = cx - HW, cx + HW, cy - HW, cy + HW
win = [i for i in range(ncpl) if x0 <= xc[i] <= x1 and y0 <= yc[i] <= y1]
fig, ax = plt.subplots(figsize=(8, 8))
ax.add_collection(PatchCollection([MplPoly(vgrid.get_cell_vertices(i)) for i in win],
                                  facecolor="none", edgecolor="0.7", lw=0.5))
ax.add_collection(PatchCollection([MplPoly(vgrid.get_cell_vertices(i)) for i in win if i in refined],
                                  facecolor="navajowhite", alpha=0.6, edgecolor="none"))
ax.add_collection(PatchCollection([MplPoly(vgrid.get_cell_vertices(i)) for i in win if i in model_sfr],
                                  facecolor="deepskyblue", alpha=0.75, edgecolor="none"))
for ln in ([stream_union] if stream_union.geom_type == "LineString" else stream_union.geoms):
    lx, ly = ln.xy; ax.plot(lx, ly, color="navy", lw=1.4)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
ax.legend(handles=[
    Patch(facecolor="navajowhite", alpha=0.6, label=f"refined corridor (≤{STREAM_BUFFER:.0f} m, small cells)"),
    Patch(facecolor="deepskyblue", alpha=0.75, label="SFR cells (centreline thread)"),
    plt.Line2D([0], [0], color="navy", lw=1.4, label="stream centreline"),
], loc="upper right", fontsize=8)
ax.set_title("STREAM_BUFFER refines cell SIZE; SFR = centreline thread only", fontsize=10)
out = Path(WS) / "diag" / "buffer_vs_sfr.png"
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"wrote {out}")
