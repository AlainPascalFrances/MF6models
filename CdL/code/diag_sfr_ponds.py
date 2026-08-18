"""
Diagnostic: do the on-channel pond cells actually get SFR reaches (LAK inactive)?

Compares, on the CURRENT all-pond grid, three sets per pond:
  (a) footprint cells          — grid cells the pond polygon occupies
  (b) stream-crossed cells     — (a) that the stream centreline crosses  -> SHOULD be SFR
  (c) model SFR cells          — cells actually in the written cdl_gwf.sfr packagedata

If (b) - (c) is non-empty, the MODEL is dropping stream-crossed pond cells from SFR
(a real bug).  If (b) is a subset of (c), every stream-crossed pond cell is an SFR
reach and the "no SFR in ponds" impression on model_map.png is a rendering/zoom effect.

Purely from files on disk (cached grid pkl + the last-written .sfr); no MODFLOW run.
"""
import pickle
import numpy as np
import geopandas as gpd
from shapely.ops import unary_union, linemerge
from shapely.geometry import LineString, MultiPolygon
from shapely.validation import make_valid
from pyproj import CRS

from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect

WS      = r"E:\00code_ws\DRYAD\CdL_model"
WS_PKL  = WS + r"\voronoi_grid.pkl"
SFR_PKG = WS + r"\cdl_gwf.sfr"
GPKG    = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
TARGET  = CRS.from_epsg(3763)

# ---- grid + intersector ----
with open(WS_PKL, "rb") as f:
    gridprops_vg, _ = pickle.load(f)
vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
try:
    ix = GridIntersect(vgrid, method="vertex")
except TypeError:
    try:
        ix = GridIntersect(vgrid, rtree=True)
    except TypeError:
        ix = GridIntersect(vgrid)
print(f"grid: ncpl = {ncpl}")

# ---- model's ACTUAL SFR cells: parse packagedata (rno layer node ...), node is 1-based ----
model_sfr = set()
with open(SFR_PKG) as f:
    inblk = False
    for line in f:
        s = line.strip()
        low = s.lower()
        if low.startswith("begin packagedata"):
            inblk = True; continue
        if low.startswith("end packagedata"):
            break
        if inblk and s and not s.startswith("#"):
            t = s.split()
            model_sfr.add(int(t[2]) - 1)        # node 1-based -> 0-based
print(f"model SFR (cdl_gwf.sfr): {len(model_sfr)} unique cells over its reaches")

# ---- watershed (clip everything, as the model does) ----
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws["geometry"] = ws.geometry.apply(lambda g: make_valid(g)).buffer(0)
ws_poly = unary_union(ws.geometry.values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)

# ---- streams (mirror the model's load/merge/dedup) ----
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(TARGET)
streams["geometry"] = streams.geometry.apply(lambda g: make_valid(g))
merged = linemerge(unary_union(streams.geometry.values))
lines = [merged] if merged.geom_type == "LineString" else list(merged.geoms)
def dedup(l):
    c = list(dict.fromkeys(l.coords)); return LineString(c) if len(c) > 1 else None
lines = [dedup(l) for l in lines if l is not None]
stream_union = unary_union([l for l in lines if l is not None]).intersection(ws_poly)
res_s = ix.intersect(stream_union, geo_dataframe=False)
geo_sfr = {int(c) for c in res_s["cellids"]}
print(f"stream-crossed cells (geometric): {len(geo_sfr)}\n")

# ---- per-pond comparison ----
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
ponds["geometry"] = ponds.geometry.apply(lambda g: make_valid(g))
print(f"{'FID':>3} {'d_str':>6} {'cells':>5} {'xStrm':>5} {'inSFR':>5} {'MISSING':>7}  channel")
print("-" * 60)
tot_missing = 0
for fid in range(len(ponds)):
    g = ponds.geometry.iloc[fid].intersection(ws_poly)
    if g.is_empty:
        continue
    if isinstance(g, MultiPolygon):
        g = max(g.geoms, key=lambda p: p.area)
    res_p = ix.intersect(g, geo_dataframe=False)
    pcells = {int(c) for c in res_p["cellids"]}
    xstrm = pcells & geo_sfr           # stream-crossed pond cells (should be SFR)
    insfr = pcells & model_sfr         # pond cells actually in the model SFR
    missing = xstrm - model_sfr        # stream-crossed but NOT in model SFR  <-- the bug signal
    tot_missing += len(missing)
    d2s = g.distance(stream_union)
    flag = "  <-- DROPPED" if missing else ""
    print(f"{fid:>3} {d2s:6.1f} {len(pcells):5d} {len(xstrm):5d} {len(insfr):5d} {len(missing):7d}  "
          f"{'ON ' if d2s <= 1.0 else 'off'}{flag}")

print(f"\nTOTAL stream-crossed pond cells MISSING from model SFR: {tot_missing}")
if tot_missing == 0:
    print("=> every stream-crossed pond cell IS an SFR reach in the model. "
          "The 'no SFR in ponds' look on model_map.png is a zoom/rendering effect, not a model gap.")
else:
    print("=> the model is NOT placing SFR in some stream-crossed pond cells — real gap to fix.")
