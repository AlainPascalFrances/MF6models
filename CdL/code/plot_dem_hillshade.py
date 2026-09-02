"""
Improved DEM figure: hypsometric map + HILLSHADE overlay to highlight the slopes
(replaces the flat colour-only dem_envelope_check.png, user 2026-08-19).

The hillshade is blended INTO the elevation colours with matplotlib's LightSource
(soft-light blend), so relief/slopes read clearly while the elevation colour scale is
preserved.  Overlays kept from the original check figure: DEM valid-data envelope (red),
raster bbox (black) and the watershed boundary (blue).

Run in Spyder / the flopy env.  Tune the light (AZDEG/ALTDEG/VERT_EXAG/BLEND_MODE) to taste.
"""
import config
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")                      # Spyder: comment out for an interactive window
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize
import rasterio
from rasterio.enums import Resampling

# ---------------- config (paths + look) ----------------
DEM_TIF = Path(str(config.MODEL) + r"\gis\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif")
GPKG    = Path(str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
WATERSHED_LAYER = "watershed_cdl_fixed"
OUT     = Path(str(config.MODEL) + r"\conceptual\dem_envelope_hillshade.png")  # set to dem_envelope_check.png to replace the old one

CMAP       = "terrain"     # terrain Blues hypsometric colour scale (blue->green->yellow->brown->white)
AZDEG      = 135.0         # illumination azimuth (NW), degrees
ALTDEG     = 5.0          # illumination altitude, degrees
VERT_EXAG  = 5.0           # vertical exaggeration — raise to make gentle slopes stronger
BLEND_MODE = "soft"        # 'soft' (natural) | 'overlay' (stronger relief) | 'hsv'
MAX_PIX    = 5000          # cap the longest raster side (decimate on read) for speed; raise for more detail
CRS        = 3763

# ---------------- 1. read the DEM (decimated to MAX_PIX) ----------------
with rasterio.open(DEM_TIF) as src:
    fw, fh = src.width, src.height
    scale = max(1, int(np.ceil(max(fw, fh) / MAX_PIX)))
    ow, oh = fw // scale, fh // scale
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.bilinear).astype("float32")
    transform = src.transform * src.transform.scale(fw / ow, fh / oh)
    nodata, bounds = src.nodata, src.bounds
print(f">> DEM read {ow}x{oh} (decimation {scale}x); pixel {transform.a:.2f} m")

if nodata is not None:
    dem[dem == nodata] = np.nan
dem[dem < -1e4] = np.nan                    # guard against sentinel nodata
valid = np.isfinite(dem)
vmin, vmax = np.nanpercentile(dem, [0.5, 99.5])
extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
dx, dy = transform.a, -transform.e          # pixel size (m); dy positive
print(f">> elevation {np.nanmin(dem):.1f}–{np.nanmax(dem):.1f} m  (colour {vmin:.1f}–{vmax:.1f})")

# ---------------- 2. hypsometric + hillshade blend ----------------
ls = LightSource(azdeg=AZDEG, altdeg=ALTDEG)
cmap = plt.get_cmap(CMAP)
demf = np.where(valid, dem, vmin)           # fill nodata so shade() has no NaNs
rgb = ls.shade(demf, cmap=cmap, vmin=vmin, vmax=vmax, blend_mode=BLEND_MODE,
               vert_exag=VERT_EXAG, dx=dx, dy=dy)     # -> (H, W, 4) RGBA
rgb[~valid, 3] = 0.0                         # nodata -> transparent

fig, ax = plt.subplots(figsize=(9, 12))
ax.set_facecolor("white")
ax.imshow(rgb, extent=extent, origin="upper", interpolation="nearest")

# colourbar (needs a ScalarMappable since we drew an RGBA image)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin, vmax)); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02, shrink=0.8)
cb.set_label("elevation (m)")

# ---------------- 3. overlays: bbox (black), valid envelope (red), watershed (blue) ----------------
# raster bbox
ax.add_patch(plt.Rectangle((bounds.left, bounds.bottom), bounds.right - bounds.left,
                           bounds.top - bounds.bottom, fill=False, ec="black", lw=1.0, zorder=4))
# DEM valid-data envelope (vectorise the valid mask)
try:
    from rasterio import features
    import shapely
    from shapely.ops import unary_union
    geoms = [shapely.geometry.shape(g) for g, v in
             features.shapes(valid.astype("uint8"), mask=valid, transform=transform) if v == 1]
    env = unary_union(geoms)
    for poly in (env.geoms if env.geom_type == "MultiPolygon" else [env]):
        xs, ys = poly.exterior.xy
        ax.plot(xs, ys, color="red", lw=1.6, zorder=5)
except Exception as e:
    print(f"   (valid-data envelope skipped: {e!r})")
# watershed
try:
    import geopandas as gpd
    ws = gpd.read_file(GPKG, layer=WATERSHED_LAYER).to_crs(CRS)
    ws.boundary.plot(ax=ax, color="blue", lw=1.4, zorder=6)
except Exception as e:
    print(f"   (watershed overlay skipped: {e!r})")

ax.set_xlabel(f"X (m, EPSG:{CRS})"); ax.set_ylabel("Y (m)")
ax.set_title(f"DEM (MDT) hypsometric + hillshade  (az {AZDEG:.0f}\u00b0, alt {ALTDEG:.0f}\u00b0, "
             f"vert.exag {VERT_EXAG:g})\nvalid envelope (red), raster bbox (black), watershed (blue)",
             fontsize=10)
ax.set_aspect("equal")
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f">> wrote {OUT}")
