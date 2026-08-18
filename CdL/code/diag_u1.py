"""Fix test: dive coordinate = distance from the U1 polygon EDGE (~p1), not from U2 (~p2). U1 overlaps U2."""
import os, numpy as np, geopandas as gpd, rasterio
from rasterio import Affine
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import rowcol
from scipy.ndimage import distance_transform_edt, gaussian_filter

OUT = r"E:\00code_ws\DRYAD\CdL_model\conceptual"
GPKG = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
DEM  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif"
HYDRO = "dryad_modelo_nbs__gc_35a_cdl_hydrostrat"; U1_BASE, DS = -10.0, 20

with rasterio.open(DEM) as src:
    nd = src.nodata; W, H = src.width, src.height; ow, oh = W // DS, H // DS
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.bilinear).astype(float)
    tr = src.transform * Affine.scale(W / ow, H / oh); crs = src.crs
cell = tr.a; dem[(dem >= 1e30) | (dem == nd)] = np.nan
g = gpd.read_file(GPKG, layer=HYDRO).to_crs(crs)
def rast(code):
    sub = g[g["Codigo"].astype(str) == code]
    return rasterize([(geom, 1) for geom in sub.geometry], out_shape=dem.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
p1, p2, p3 = rast("1"), rast("2"), rast("3"); p3 = p3 | ~(p1 | p2 | p3)
print(f"U1&U2 overlap cells: {(p1&p2).sum()} ({100*(p1&p2).sum()/p1.sum():.0f}% of U1)")
d1 = distance_transform_edt(p1) * cell                      # depth INTO the U1 polygon from its edge (0 at edge)
d1 = gaussian_filter(d1, sigma=10)                          # smooth the transect weave against the irregular edge
_demf = dem.copy()
if np.isnan(_demf).any():
    _, _j = distance_transform_edt(np.isnan(_demf), return_indices=True); _demf = _demf[_j[0], _j[1]]
dem_smooth = gaussian_filter(_demf, sigma=12)
print(f"U1 d1 (dist from U1 edge) range: 0..{np.nanmax(np.where(p1,d1,np.nan)):.0f} m")

pf = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
tx, ty, sdist = pf["xs"], pf["ys"], pf["dists"]
rr, cc = rowcol(tr, tx, ty); rr = np.clip(rr, 0, dem.shape[0]-1); cc = np.clip(cc, 0, dem.shape[1]-1)
inU1 = p1[rr, cc]
print(f"\nsdist  dem   d1 | bot1@150 thk@150 | bot1@120 thk@120")
for k in range(len(sdist)):
    if inU1[k] and k % 25 == 0:
        row = f"{sdist[k]:5.0f} {dem[rr[k],cc[k]]:4.1f} {d1[rr[k],cc[k]]:4.0f} |"
        for D in (150, 120):
            r1 = min(d1[rr[k], cc[k]] / D, 1.0)
            b = min(U1_BASE + (dem_smooth[rr[k], cc[k]] - U1_BASE) * (1 - r1) ** 2, dem[rr[k], cc[k]])
            row += f"  {b:6.1f} {dem[rr[k],cc[k]]-b:5.1f} |"
        print(row)
