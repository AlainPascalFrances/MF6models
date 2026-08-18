"""
PREVIEW the smooth-transition Voronoi grid (Daoud et al. 2022 style) WITHOUT touching
the model workspace or its cached grid — builds Triangle+Voronoi to a scratch dir and
plots cell size (sqrt area) zoomed on a pond and a stream reach, to verify the grading
steps smoothly (2.5 -> 6 -> 12 -> 25 -> 60 -> stream/far) before adopting it in the model.
"""
import numpy as np
import geopandas as gpd
from shapely.ops import unary_union, linemerge
from shapely.geometry import LineString, MultiPolygon, Polygon as ShPoly
from shapely.validation import make_valid
from pyproj import CRS
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly
from matplotlib.collections import PatchCollection
import flopy
from flopy.discretization import VertexGrid
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid

GPKG   = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
TRIEXE = r"C:\00MODFLOW\win64\triangle.exe"
SCRATCH = r"C:\Users\ALAIN~1.FRA\AppData\Local\Temp\claude\E--tmp-claude\e6df0a7e-94fe-4a92-af96-7d7d05484f82\scratchpad\tri_preview"
OUT    = r"E:\00code_ws\DRYAD\CdL_model\diag\transition_grid_preview.png"
TARGET = CRS.from_epsg(3763)
import os; os.makedirs(SCRATCH, exist_ok=True)

# --- params (mirror the model) ---
CELL_NEAR_STREAM, CELL_FAR, STREAM_BUFFER, POND_CELL = 40.0, 100.0, 60.0, 5.0   # match the model (user redesign)

# --- geometry (mirror the model) ---
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws["geometry"] = ws.geometry.apply(lambda g: make_valid(g)).buffer(0)
ws_poly = unary_union(ws.geometry.values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(TARGET)
merged = linemerge(unary_union(streams.geometry.apply(make_valid).values))
stream_union = merged if merged.geom_type == "LineString" else unary_union(list(merged.geoms))
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
ponds["geometry"] = ponds.geometry.apply(make_valid)
pond_polys_all = [ (g if not isinstance(g, MultiPolygon) else max(g.geoms, key=lambda p: p.area))
                   for g in (p.intersection(ws_poly) for p in ponds.geometry) if not g.is_empty ]

n_outer = max(int(np.ceil(ws_poly.exterior.length / CELL_FAR)), 50)
domain_xy = [ws_poly.exterior.interpolate(i / n_outer, normalized=True).coords[0] for i in range(n_outer)]
ws_inner = ws_poly.buffer(-STREAM_BUFFER)
if isinstance(ws_inner, MultiPolygon):
    ws_inner = max(ws_inner.geoms, key=lambda p: p.area)
streams_clipped = stream_union.intersection(ws_inner)
stream_buf = unary_union(streams_clipped.buffer(STREAM_BUFFER, cap_style=1, join_style=1))
stream_buf_polys = list(stream_buf.geoms) if isinstance(stream_buf, MultiPolygon) else [stream_buf]

# --- Triangle PSLG (mirror the model) ---
tri = Triangle(model_ws=SCRATCH, angle=20, exe_name=TRIEXE, additional_args=["-j"])
tri.add_polygon(domain_xy)
coarse_pt = ws_poly.difference(unary_union(stream_buf_polys)).representative_point()
tri.add_region((coarse_pt.x, coarse_pt.y), attribute=0, maximum_area=CELL_FAR ** 2)
# pond cores at POND_CELL (exact pond boundaries)
for _pg in pond_polys_all:
    _r = _pg.exterior; _n = max(int(np.ceil(_r.length / POND_CELL)), 8)
    tri.add_polygon([_r.interpolate(i / _n, normalized=True).coords[0] for i in range(_n)])
    _rp = _pg.representative_point(); tri.add_region((_rp.x, _rp.y), attribute=2, maximum_area=POND_CELL ** 2)

# CONTINUOUS nested transition zones from a sizing field:
#   size(x) = min( POND_CELL + GRADE*d_pond , CELL_NEAR_STREAM + GRADE*max(0, d_stream-STREAM_BUFFER) , CELL_FAR )
# The buffer distance for level L grows CONTINUOUSLY with L, so the level sets {size<=L}
# STRICTLY nest -> distinct, non-coincident boundaries (fixes Triangle's topological error).
GRADE  = 0.5
LEVELS = [10.0, 20.0, 40.0, 70.0]   # match the model TRANS_LEVELS (POND_CELL 5 -> ... -> CELL_FAR 100)
ponds_u = unary_union(pond_polys_all)
def _clip_then_buffer(geom, B):
    # clip the FEATURE so its B-buffer reaches the watershed edge only TANGENTIALLY (no
    # shared clip-edge between levels -> avoids Triangle's coincident-segment failure)
    inside = ws_poly.buffer(-B)
    if inside.is_empty:
        return None
    g = geom.intersection(inside)
    return None if g.is_empty else g.buffer(B)
prev = ponds_u
for L in LEVELS:
    parts = []
    gp = _clip_then_buffer(ponds_u, (L - POND_CELL) / GRADE)
    if gp is not None:
        parts.append(gp)
    if L >= CELL_NEAR_STREAM:
        gs = _clip_then_buffer(streams_clipped, STREAM_BUFFER + (L - CELL_NEAR_STREAM) / GRADE)
        if gs is not None:
            parts.append(gs)
    if not parts:
        continue
    G = unary_union(parts)
    if G.is_empty:
        continue
    for poly in (G.geoms if isinstance(G, MultiPolygon) else [G]):
        for ring in [poly.exterior] + list(poly.interiors):
            n = max(int(np.ceil(ring.length / L)), 8)
            tri.add_polygon([ring.interpolate(i / n, normalized=True).coords[0] for i in range(n)])
    annulus = G.difference(prev)
    for ap in (annulus.geoms if isinstance(annulus, MultiPolygon) else [annulus]):
        if ap.is_empty or ap.area < L ** 2:
            continue
        sp = ap.representative_point(); tri.add_region((sp.x, sp.y), attribute=int(L), maximum_area=L ** 2)
    prev = G

print("running triangle …"); tri.build(verbose=False)
vor = VoronoiGrid(tri)
vgrid = VertexGrid(**vor.get_gridprops_vertexgrid(), nlay=1)
ncpl = vgrid.ncpl
xc = np.array([vgrid.xcellcenters[i] for i in range(ncpl)])
yc = np.array([vgrid.ycellcenters[i] for i in range(ncpl)])
verts = [vgrid.get_cell_vertices(i) for i in range(ncpl)]
area = np.array([ShPoly(v).area for v in verts])
width = np.sqrt(area)        # ~ cell width (m)
print(f"ncpl = {ncpl};  cell width min/median/max = {width.min():.1f}/{np.median(width):.1f}/{width.max():.1f} m")

# --- plot: full + zoom on the biggest pond + a stream reach ---
big = max(range(len(pond_polys_all)), key=lambda i: pond_polys_all[i].area)
pc = pond_polys_all[big].centroid
sp = stream_union.interpolate(0.45, normalized=True)
fig, axes = plt.subplots(1, 3, figsize=(18, 7))
for ax, (cx, cy, hw, ttl) in zip(axes, [
        (ws_poly.centroid.x, ws_poly.centroid.y, 4000, f"full grid (ncpl={ncpl})"),
        (pc.x, pc.y, 130, "zoom: pond (2.5->6->12->25->60 m)"),
        (sp.x, sp.y, 260, "zoom: stream (40->75->120->150 m)")]):
    win = [i for i in range(ncpl) if abs(xc[i] - cx) <= hw and abs(yc[i] - cy) <= hw]
    pcoll = PatchCollection([MplPoly(verts[i]) for i in win], cmap="viridis")
    pcoll.set_array(width[win]); pcoll.set_edgecolor("0.4"); pcoll.set_linewidth(0.2)
    ax.add_collection(pcoll)
    gpd.GeoSeries(pond_polys_all, crs=TARGET).boundary.plot(ax=ax, color="red", lw=0.8)
    for ln in ([stream_union] if stream_union.geom_type == "LineString" else stream_union.geoms):
        lx, ly = ln.xy; ax.plot(lx, ly, color="white", lw=0.8)
    ax.set_xlim(cx - hw, cx + hw); ax.set_ylim(cy - hw, cy + hw); ax.set_aspect("equal")
    ax.set_title(ttl, fontsize=10); plt.colorbar(pcoll, ax=ax, label="cell width (m)", shrink=0.6)
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print(f"wrote {OUT}")
