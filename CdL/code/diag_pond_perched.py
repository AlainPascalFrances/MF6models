"""
Quantify how PERCHED each charca is: pond bottom elevation vs the steady-state
water table (no-lake SS heads in spinup_heads.npy).

  perched margin = WT - pond_bottom   (negative => WT BELOW pond bottom = perched)

Segments the 19 ponds into 'in-contact' (LAK behaves) vs 'perched' (the hard ones),
and flags good Phase-0 subset candidates. Files on disk only; no MODFLOW run.
"""
import config
import pickle
import numpy as np
import geopandas as gpd
from shapely.ops import unary_union, linemerge
from shapely.geometry import MultiPolygon
from shapely.validation import make_valid
from pyproj import CRS
from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect

WS      = str(config.MODEL)
WS_PKL  = WS + r"\voronoi_grid.pkl"
HEADS   = WS + r"\spinup_heads.npy"                          # (nlay, ncpl) no-lake SS IC
LAYERS  = WS + r"\conceptual\layers\voronoi_layers.npz"     # has 'top'
GPKG    = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
TARGET  = CRS.from_epsg(3763)
POND_DEPTH = 4.0

with open(WS_PKL, "rb") as f:
    gridprops_vg, _ = pickle.load(f)
vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
try:
    ix = GridIntersect(vgrid, method="vertex")
except TypeError:
    ix = GridIntersect(vgrid)

heads = np.load(HEADS)                      # (4, ncpl)
heads = np.where(np.abs(heads) < 1e29, heads, np.nan)
wt = np.nanmax(heads, axis=0)               # water table = highest saturated head in the column
top = np.load(LAYERS)["top"].astype(float)  # model top (land surface) per cell

ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws_poly = unary_union(ws.geometry.apply(lambda g: make_valid(g)).buffer(0).values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(TARGET)
merged = linemerge(unary_union(streams.geometry.apply(make_valid).values))
stream_union = (merged if merged.geom_type == "LineString"
                else unary_union(list(merged.geoms))).intersection(ws_poly)
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
ponds["geometry"] = ponds.geometry.apply(make_valid)

print(f"{'FID':>3} {'area':>5} {'d_str':>6} {'surf':>6} {'pbot':>6} {'WT':>6} {'margin':>7}  class")
print("-" * 62)
rows = []
for fid in range(len(ponds)):
    g = ponds.geometry.iloc[fid].intersection(ws_poly)
    if g.is_empty:
        continue
    if isinstance(g, MultiPolygon):
        g = max(g.geoms, key=lambda p: p.area)
    cells = [int(c) for c in ix.intersect(g, geo_dataframe=False)["cellids"]]
    surf = float(np.nanmean(top[cells]))
    pbot = surf - POND_DEPTH
    wtl  = float(np.nanmean(wt[cells]))
    margin = wtl - pbot                      # <0 => perched
    d2s  = g.distance(stream_union)
    cls  = "PERCHED" if margin < -0.5 else ("contact" if margin >= 0 else "marginal")
    rows.append((fid, margin, d2s, cls))
    print(f"{fid:>3} {g.area:5.0f} {d2s:6.1f} {surf:6.1f} {pbot:6.1f} {wtl:6.1f} {margin:7.1f}  "
          f"{cls}{'  [on-ch]' if d2s <= 1.0 else ''}")

perched = [r for r in rows if r[3] == "PERCHED"]
contact = [r for r in rows if r[3] == "contact"]
print(f"\n{len(perched)} perched, {len(contact)} in-contact, {len(rows)-len(perched)-len(contact)} marginal.")
print("Phase-0 candidates:")
if contact:
    fid_c = min(contact, key=lambda r: abs(r[1]))[0]
    print(f"  in-contact (easy)  : FID {fid_c}  (margin {dict((r[0],r[1]) for r in rows)[fid_c]:+.1f} m)")
if perched:
    fid_p = min(perched, key=lambda r: r[1])[0]   # most perched
    print(f"  perched (hard)     : FID {fid_p}  (margin {dict((r[0],r[1]) for r in rows)[fid_p]:+.1f} m)")
print("  (pond 6 is the documented killer — include it explicitly.)")
