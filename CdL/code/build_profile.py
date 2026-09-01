r"""CdL / DRYAD conceptual model — STEP 1 of 2.
Builds the topographic-section transect  START(SE watershed divide) -> ERT P1 -> P2 -> P3 -> P4
(polyline through the ERT-profile centroids) and samples the 0.5 m LiDAR DEM along it.
Writes  <OUT>/profile_data.npz , consumed by  cross_section.py  (STEP 2).

Run in the flopy env, e.g.:
  & "C:\miniconda3\Scripts\conda.exe" run -p C:\sw\miniconda3\envs\mf6models --no-capture-output python -u build_profile.py
"""
import config
import os
import numpy as np
import geopandas as gpd
import rasterio
from shapely.geometry import LineString, Point

OUT  = (str(config.MODEL) + r"\conceptual")
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
DEM  = (str(config.MODEL) + r"\gis\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif")
STEP = 2.0  # along-transect sampling interval (m)
os.makedirs(OUT, exist_ok=True)

# ---- ERT centroids (EPSG:3763) + watershed boundary ----
ert = gpd.read_file(GPKG, layer="ert_electrodes_cdl")
ws_bnd = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763).geometry.union_all().boundary
cent = {e: np.array([ert[ert["ert"] == e].geometry.x.mean(),
                     ert[ert["ert"] == e].geometry.y.mean()]) for e in (1, 2, 3, 4)}
P1, P2, P3, P4 = cent[1], cent[2], cent[3], cent[4]

# ---- START on the SE watershed limit: extend the P2->P1 direction (ESE) from P1 to the boundary ----
u = (P1 - P2) / np.hypot(*(P1 - P2))
inter = LineString([tuple(P1), tuple(P1 + 8000.0 * u)]).intersection(ws_bnd)
pts = [inter] if isinstance(inter, Point) else [p for p in getattr(inter, "geoms", []) if isinstance(p, Point)]
pts.sort(key=lambda p: Point(*P1).distance(p))      # nearest boundary crossing = first exit (SE)
START = np.array([pts[0].x, pts[0].y])

# ---- END on the NW watershed limit: extend the P3->P4 direction beyond P4 to the boundary ----
u2 = (P4 - P3) / np.hypot(*(P4 - P3))
inter2 = LineString([tuple(P4), tuple(P4 + 8000.0 * u2)]).intersection(ws_bnd)
pts2 = [inter2] if isinstance(inter2, Point) else [p for p in getattr(inter2, "geoms", []) if isinstance(p, Point)]
pts2.sort(key=lambda p: Point(*P4).distance(p))     # nearest boundary crossing past P4 (NW)
END = np.array([pts2[0].x, pts2[0].y])

# ---- transect polyline + cumulative vertex distances ----
verts = [START, P1, P2, P3, P4, END]
vert_names = ["START\n(SE divide)", "P1", "P2", "P3", "P4", "END\n(NW limit)"]
line = LineString([tuple(v) for v in verts])
vert_dist = np.concatenate([[0.0], np.cumsum([np.hypot(*(verts[i] - verts[i - 1]))
                                              for i in range(1, len(verts))])])

# ---- densify + sample DEM ----
dists = np.linspace(0, line.length, int(line.length // STEP) + 1)
sp = [line.interpolate(dd) for dd in dists]
xs = np.array([p.x for p in sp]); ys = np.array([p.y for p in sp])
with rasterio.open(DEM) as src:
    nd = src.nodata
    zz = np.array([v[0] for v in src.sample(zip(xs, ys))], dtype=float)
zz[(zz == nd) | (zz > 1e30)] = np.nan

np.savez(os.path.join(OUT, "profile_data.npz"),
         dists=dists, xs=xs, ys=ys, zz=zz, vert_dist=vert_dist,
         verts=np.array(verts), vert_names=np.array(vert_names, dtype=object))

print(f"START (SE divide) = ({START[0]:.1f}, {START[1]:.1f})")
print(f"transect length   = {line.length:.1f} m   vertices(d): " +
      ", ".join(f"{n.splitlines()[0]}={d:.0f}" for n, d in zip(vert_names, vert_dist)))
print(f"DEM elevation along transect: {np.nanmin(zz):.1f} .. {np.nanmax(zz):.1f} m  (NaN={np.isnan(zz).sum()})")
print("saved", os.path.join(OUT, "profile_data.npz"))
