"""Plan-view map of the 3 hydrostratigraphic units (gpkg layer dryad_modelo_nbs__gc_35a_cdl_hydrostrat,
Codigo 1/2/3) draped over the LiDAR-DTM hillshade — same style as geology_plan.png (ERT P1-P6, streams,
watershed, transect).  -> <OUT>\\hydrostrat_plan.png"""
import config
import os, copy, numpy as np, geopandas as gpd, rasterio
from rasterio import Affine
from rasterio.enums import Resampling
from matplotlib.colors import LightSource
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

OUT  = (str(config.MODEL) + r"\conceptual")
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
DEM  = (str(config.MODEL) + r"\gis\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif")
HYDRO = "dryad_modelo_nbs__gc_35a_cdl_hydrostrat"
halo = [pe.withStroke(linewidth=2.2, foreground="white")]

UCOL = {"1": ("#d9d9d9", "U1 — alluvium · K = 25 m/d"),
        "2": ("#6b6b6b", "U2 — terraces · K ≈ 16–77 m/d"),
        "3": ("#fdf6a8", "U3 — Plio-Miocene · K ≈ 1–2 m/d")}
g = gpd.read_file(GPKG, layer=HYDRO).to_crs(3763)
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
wsg = ws.geometry.union_all()

# transect (from profile_data.npz)
d = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
xs, ys, verts = d["xs"], d["ys"], d["verts"]

# hillshade from the LiDAR DTM (decimated)
DS = 6
with rasterio.open(DEM) as src:
    nd = src.nodata; W, H = src.width, src.height; ow, oh = W // DS, H // DS
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.nearest).astype(float)
    tr = src.transform * Affine.scale(W / ow, H / oh); db = src.bounds
nd_mask = (dem >= 1e30) | (dem == nd); dem[nd_mask] = np.nan
ls = LightSource(azdeg=315, altdeg=45)
hs = ls.hillshade(np.where(np.isnan(dem), np.nanmean(dem), dem), vert_exag=4.0, dx=tr.a, dy=-tr.e)
hs[nd_mask] = np.nan
dem_extent = [db.left, db.right, db.bottom, db.top]
gray = copy.copy(plt.cm.gray); gray.set_bad(alpha=0.0)

fig, ax = plt.subplots(figsize=(11, 12.5)); ax.set_aspect("equal")
ax.imshow(np.ma.masked_invalid(hs), extent=dem_extent, origin="upper", cmap=gray, vmin=0.15, vmax=1.0, zorder=0)
# units: draw U3, then U2, then U1 on top (U1 overlaps U2 in the north)
for c in ["3", "2", "1"]:
    sub = g[g["Codigo"].astype(str) == c]
    if len(sub):
        sub.plot(ax=ax, color=UCOL[c][0], ec="0.4", lw=0.3, alpha=0.55, zorder=1)
        gp = sub.geometry.union_all().intersection(wsg)          # label placed inside the watershed
        rp = (gp if not gp.is_empty else sub.geometry.union_all()).representative_point()
        ax.text(rp.x, rp.y, "U" + c, ha="center", va="center", fontsize=14, fontweight="bold",
                path_effects=halo, zorder=7)
# streams (white casing + light blue)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(3763)
streams.plot(ax=ax, color="white", lw=3.0, alpha=0.7, zorder=2.8)
streams.plot(ax=ax, color="#2BAEEA", lw=1.5, zorder=3.0)
# watershed
ws.boundary.plot(ax=ax, edgecolor="blue", lw=1.3, zorder=4)
# cattle ponds (charcas): navy polygons + small FID labels (positional index = model FID; font < unit labels)
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(3763).reset_index(drop=True)
ponds.plot(ax=ax, color="navy", ec="white", lw=0.4, alpha=0.9, zorder=5.4)
for i, geom in enumerate(ponds.geometry):
    rp = geom.representative_point()
    ax.plot(rp.x, rp.y, marker="o", ms=3.2, color="navy", mec="white", mew=0.5, zorder=5.6)
    ax.annotate(str(i), (rp.x, rp.y), textcoords="offset points", xytext=(3.5, 3.0),
                fontsize=7, fontweight="bold", color="navy", path_effects=halo, zorder=10)
# ERT P1-P6 (red on transect, violet others)
ert = gpd.read_file(GPKG, layer="ert_electrodes_cdl")
for e in sorted(ert["ert"].unique()):
    sub = ert[ert["ert"] == e].sort_values("id_electro"); on = e in (1, 2, 3, 4)
    ax.plot(sub.geometry.x, sub.geometry.y, "-", color=("crimson" if on else "darkviolet"),
            lw=3.2, zorder=6, solid_capstyle="round")
    ax.annotate(f"P{e}", (sub.geometry.x.mean(), sub.geometry.y.mean()), textcoords="offset points",
                xytext=(7, 7), fontsize=12, fontweight="bold",
                color=("crimson" if on else "darkviolet"), path_effects=halo, zorder=10)
# transect
ax.plot(xs, ys, "-", color="black", lw=2.0, zorder=5)
ax.plot(verts[0][0], verts[0][1], "o", color="black", ms=7, zorder=8)
ax.annotate("START\n(SE divide)", (verts[0][0], verts[0][1]), textcoords="offset points",
            xytext=(8, -2), fontsize=10, fontweight="bold", path_effects=halo, zorder=10)
ax.plot(verts[-1][0], verts[-1][1], "o", color="black", ms=7, zorder=8)
ax.annotate("END\n(NW limit)", (verts[-1][0], verts[-1][1]), textcoords="offset points",
            xytext=(9, 4), fontsize=10, fontweight="bold", path_effects=halo, zorder=10)
# legends
leg = ax.legend(handles=[Patch(fc=UCOL[c][0], ec="0.5", label=UCOL[c][1]) for c in ["1", "3"]],
                loc="lower left", fontsize=9, title="Hydrostratigraphic units (U2 labelled on map)")
ax.add_artist(leg)
ax.legend(handles=[Line2D([0], [0], color="crimson", lw=3, label="ERT on transect (P1–P4)"),
                   Line2D([0], [0], color="darkviolet", lw=3, label="ERT other (P5, P6)"),
                   Line2D([0], [0], color="black", lw=2, label="Transect"),
                   Line2D([0], [0], color="#2BAEEA", lw=2, label="Streams"),
                   Line2D([0], [0], color="blue", lw=1.3, label="Watershed"),
                   Line2D([0], [0], marker="o", color="navy", mec="white", mew=0.5, ls="none",
                          ms=6, label="Cattle ponds (charcas) — FID")],
          loc="upper right", fontsize=8.5)
ax.set_title("Hydrostratigraphic units (U1/U2/U3) over LiDAR-DTM hillshade — ERT P1–P6, transect & watershed")
ax.annotate("Hillshade: LiDAR DTM 0.5 m\nillum. 315°/45°, VE 4×\nunits α = 0.55", (0.99, 0.012),
            xycoords="axes fraction", ha="right", va="bottom", fontsize=7.5, color="0.25",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.9))
ax.set_xlabel("X (m, EPSG:3763)"); ax.set_ylabel("Y (m)")
fig.savefig(os.path.join(OUT, "hydrostrat_plan.png"), dpi=160, bbox_inches="tight")
print("saved hydrostrat_plan.png")
