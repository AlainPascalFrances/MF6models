"""
Derive each pond's CONTRIBUTING CATCHMENT AREA from the DEM (for the deep-U3
recharge-fallback: focused recharge = P + captured runoff from the catchment).

Self-contained D8 flow routing (no flow-routing library needed):
  1. DEM decimated to ~20 m, clipped to the watershed (a drainage divide).
  2. Priority-flood pit fill (Barnes et al. 2014) so flow never traps in sinks.
  3. D8 steepest-descent flow direction.
  4. Flow accumulation (topological, high->low filled elevation).
  5. Per pond: contributing area = max accumulation over its footprint x cell area.

-> prints a per-pond table and writes diag\pond_catchments.png (accumulation + ponds).
"""
import config
import heapq
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio import Affine
from rasterio.features import rasterize
import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.validation import make_valid
from shapely.ops import unary_union
from pyproj import CRS
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEM   = (str(config.MODEL) + r"\gis\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif")
GPKG  = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
OUT   = (str(config.MODEL) + r"\diag\pond_catchments.png")
TARGET = CRS.from_epsg(3763)
DS = 40   # 0.5 m DEM -> 20 m working grid (drainage areas don't need finer)
DEEP_U3 = {1, 2, 3, 8, 9, 10, 12, 13, 14}   # the recharge-fallback ponds (L4 connectors)

# ---- 1. DEM decimated + clipped to watershed ----
with rasterio.open(DEM) as src:
    nd = src.nodata; W, H = src.width, src.height; ow, oh = W // DS, H // DS
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.bilinear).astype(float)
    tr = src.transform * Affine.scale(W / ow, H / oh)
cell = tr.a
dem[(dem >= 1e30) | (dem == nd)] = np.nan
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws_poly = unary_union(ws.geometry.apply(lambda g: make_valid(g)).buffer(0).values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
inside = rasterize([(ws_poly, 1)], out_shape=dem.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
dem[~inside] = np.nan
nodata = np.isnan(dem)
h, w = dem.shape
print(f"grid {dem.shape} @ {cell:.0f} m, {int((~nodata).sum())} cells in watershed")

NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
DIST = np.array([2 ** .5, 1, 2 ** .5, 1, 1, 2 ** .5, 1, 2 ** .5])

# ---- 2. priority-flood pit fill ----
filled = dem.copy()
closed = nodata.copy()
pq = []
for i in range(h):
    for j in range(w):
        if nodata[i, j]:
            continue
        if i in (0, h - 1) or j in (0, w - 1) or any(
            0 <= i + di < h and 0 <= j + dj < w and nodata[i + di, j + dj] for di, dj in NB):
            heapq.heappush(pq, (dem[i, j], i, j)); closed[i, j] = True
while pq:
    e, i, j = heapq.heappop(pq)
    for di, dj in NB:
        ni, nj = i + di, j + dj
        if 0 <= ni < h and 0 <= nj < w and not closed[ni, nj] and not nodata[ni, nj]:
            filled[ni, nj] = max(dem[ni, nj], e)
            closed[ni, nj] = True
            heapq.heappush(pq, (filled[ni, nj], ni, nj))

# ---- 3. D8 steepest-descent downstream cell ----
ds_idx = np.full((h, w), -1, dtype=np.int64)   # flat index of downstream cell, -1 = outlet/none
for i in range(h):
    for j in range(w):
        if nodata[i, j]:
            continue
        best, bestslope = -1, 0.0
        for k, (di, dj) in enumerate(NB):
            ni, nj = i + di, j + dj
            if 0 <= ni < h and 0 <= nj < w and not nodata[ni, nj]:
                slope = (filled[i, j] - filled[ni, nj]) / (DIST[k] * cell)
                if slope > bestslope:
                    bestslope, best = slope, ni * w + nj
        ds_idx[i, j] = best

# ---- 4. flow accumulation (process high -> low filled elevation) ----
acc = np.where(nodata, 0.0, 1.0)
fe = np.where(nodata, -np.inf, filled)
order = np.argsort(fe.ravel())[::-1]
accf = acc.ravel(); dsf = ds_idx.ravel()
for idx in order:
    d = dsf[idx]
    if d >= 0:
        accf[d] += accf[idx]
acc = accf.reshape(h, w)

# ---- 5. per-pond contributing area ----
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
ponds["geometry"] = ponds.geometry.apply(make_valid)
cell_area = cell * cell
print(f"\n{'FID':>3} {'pond_m2':>8} {'catch_ha':>9} {'ratio':>7}  group")
print("-" * 46)
results = {}
for fid in range(len(ponds)):
    g = ponds.geometry.iloc[fid].intersection(ws_poly)
    if g.is_empty:
        continue
    if isinstance(g, MultiPolygon):
        g = max(g.geoms, key=lambda p: p.area)
    mask = rasterize([(g.buffer(cell), 1)], out_shape=dem.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
    if not mask.any():
        ci, cj = [int(x) for x in (~tr * (g.centroid.x, g.centroid.y))][::-1]
        mask[ci, cj] = True
    a_contrib = float(np.nanmax(acc[mask])) * cell_area
    a_pond = float(g.area)
    results[fid] = (a_pond, a_contrib)
    grp = "DEEP-U3 (recharge)" if fid in DEEP_U3 else "LAK / other"
    print(f"{fid:>3} {a_pond:8.0f} {a_contrib/1e4:9.2f} {a_contrib/a_pond:7.1f}  {grp}")

# ---- save per-pond areas for the model (recharge fallback reads this) ----
import csv
CSVOUT = (str(config.MODEL) + r"\pond_catchments.csv")
with open(CSVOUT, "w", newline="") as f:
    wtr = csv.writer(f); wtr.writerow(["fid", "pond_area_m2", "catchment_area_m2", "ratio"])
    for fid in sorted(results):
        ap, ac = results[fid]; wtr.writerow([fid, round(ap, 1), round(ac, 1), round(ac / ap, 2)])
print(f"wrote {CSVOUT}")

# ---- map ----
fig, ax = plt.subplots(figsize=(8, 9))
im = ax.imshow(np.log10(np.where(acc > 0, acc, np.nan)), cmap="Blues",
               extent=[tr.c, tr.c + w * cell, tr.f + h * tr.e, tr.f])
ponds.boundary.plot(ax=ax, color="red", lw=0.8)
for fid in DEEP_U3:
    if fid in results:
        c = ponds.geometry.iloc[fid].centroid
        ax.annotate(str(fid), (c.x, c.y), fontsize=7, color="k")
plt.colorbar(im, ax=ax, label="log10(flow accumulation, cells)", shrink=0.6)
ax.set_title("Flow accumulation + pond catchments (deep-U3 ponds labelled)")
ax.set_aspect("equal")
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"\nwrote {OUT}")
