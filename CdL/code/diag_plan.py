"""Plan view: U1/U2/U3 + transect + d23 field + U1's d23-min tip, to understand how the transect crosses U1."""
import os, numpy as np, geopandas as gpd, rasterio
from rasterio import Affine
from rasterio.enums import Resampling
from rasterio.features import rasterize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt

OUT = r"E:\00code_ws\DRYAD\CdL_model\conceptual"
GPKG = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
DEM  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif"
HYDRO = "dryad_modelo_nbs__gc_35a_cdl_hydrostrat"; DS = 20

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
d23 = distance_transform_edt(~p3) * cell

comp = np.full(dem.shape, np.nan); comp[p3] = 0; comp[p2] = 1; comp[p1] = 2
ext = [tr.c, tr.c + tr.a*dem.shape[1], tr.f + tr.e*dem.shape[0], tr.f]
pf = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
tx, ty, sdist = pf["xs"], pf["ys"], pf["dists"]

# zoom to U1 bbox (+ buffer)
ys, xs = np.where(p1); xmin = tr.c + xs.min()*tr.a; xmax = tr.c + xs.max()*tr.a
ymax = tr.f + ys.min()*tr.e; ymin = tr.f + ys.max()*tr.e; buf = 400

fig, ax = plt.subplots(1, 2, figsize=(18, 9))
from matplotlib.colors import ListedColormap
ax[0].imshow(comp, extent=ext, origin="upper", cmap=ListedColormap(["#fdf6a8", "#6b6b6b", "#d9d9d9"]))
ax[0].plot(tx, ty, "r-", lw=2)
# mark sdist ticks along the transect every 200 m
for s in range(0, int(sdist.max()), 200):
    k = int(np.argmin(np.abs(sdist - s))); ax[0].annotate(f"{s}", (tx[k], ty[k]), fontsize=6, color="red")
d23u1 = np.where(p1, d23, np.nan); rmin, cmin = np.unravel_index(np.nanargmin(d23u1), d23u1.shape)
ax[0].plot(tr.c+(cmin+.5)*tr.a, tr.f+(rmin+.5)*tr.e, "b*", ms=18)
ax[0].set_title("units (yellow U3 / grey U2 / light U1) + transect (red, sdist labels) + U1 d23-min tip (blue star)")
ax[0].set_xlim(xmin-buf, xmax+buf); ax[0].set_ylim(ymin-buf, ymax+buf); ax[0].set_aspect("equal")

im = ax[1].imshow(np.where(p1, d23, np.nan), extent=ext, origin="upper", cmap="turbo")
cs = ax[1].contour(np.where(p1|p2, d23, np.nan), levels=range(0, 3500, 250), extent=ext, origin="upper", colors="k", linewidths=0.5)
ax[1].clabel(cs, fontsize=6); ax[1].plot(tx, ty, "r-", lw=2)
ax[1].set_title("d23 (dist from U3 outcrop) on U1, + d23 contours; how does the transect cross the d23 field?")
ax[1].set_xlim(xmin-buf, xmax+buf); ax[1].set_ylim(ymin-buf, ymax+buf); ax[1].set_aspect("equal")
plt.colorbar(im, ax=ax[1], shrink=0.7, label="d23 (m)")
fig.savefig(os.path.join(OUT, "diag_plan.png"), dpi=130, bbox_inches="tight"); print("saved diag_plan.png")
print(f"U1 bbox x {xmin:.0f}..{xmax:.0f} y {ymin:.0f}..{ymax:.0f} ; d23-min tip at x={tr.c+(cmin+.5)*tr.a:.0f} y={tr.f+(rmin+.5)*tr.e:.0f}")
