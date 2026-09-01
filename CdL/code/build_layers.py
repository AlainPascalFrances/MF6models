"""Build the top/bottom ELEVATION surfaces (m a.s.l.) + thickness of the 3 hydrostratigraphic units of
the CdL/DRYAD model, from the DTM + the user's modified hydrostrat map (gpkg layer, Codigo 1/2/3, where
Codigo 1/2/3), per the user's recipe (2026-06-23):

  U1 alluvium      : top = DTM ; bottom = DTM (0 thickness) at the U1/U2 contact, dives PARABOLICALLY to
                     U1_BASE (-10 m) over D_DIVE (300 m N from the contact = dist from U2).  thk = top - bottom.
  U2 terraces      : BOTTOM = DTM (0 thickness) at the U2/U3 contact, dives (linear) to a FLAT FLOOR (+1 m)
                     over D_STEEP (500 m).  TOP = DTM where U2 outcrops, = U1 bottom where unit 1 overlaps U2;
                     the U2 bottom is kept <= DTM (always below topography).  thickness = top - bottom.
  U3 Plio-Miocene  : bottom = AQ_BASE (-35 m).  top = DTM where U3 outcrops, else the bottom of the unit
                     directly above it (U2 bottom under the terraces, U1 bottom under the alluvium).

Outputs in <OUT>\\layers\\ : GeoTIFF top/bot/thk per unit ; voronoi_layers.npz/.csv (top + botm[3,ncpl]
at the DISV centroids, MODFLOW-ready) ; layers_maps.png (top/bottom/thickness maps, 3x3) ;
layers_check.png (cross-section validation vs the transect).
"""
import config
import os, csv, pickle
import numpy as np
import geopandas as gpd
import rasterio
from rasterio import Affine
from rasterio.enums import Resampling
from rasterio.features import rasterize, geometry_mask
from rasterio.transform import rowcol
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = (str(config.MODEL) + r"\conceptual")
LAYDIR = os.path.join(OUT, "layers"); os.makedirs(LAYDIR, exist_ok=True)
WS   = str(config.MODEL)
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
DEM  = (str(config.MODEL) + r"\gis\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif")
HYDRO = "dryad_modelo_nbs__gc_35a_cdl_hydrostrat"

# ---- recipe parameters (tunable) — user's 2026-06-23 revision ----
# U1 (alluvium) bottom = DTM (0 thickness) at the U1/U2 contact, dives PARABOLICALLY to U1_BASE over D_DIVE.
# U2 (terraces) bottom = DTM (0 thickness) at the U2/U3 contact, dives LINEARLY to FLOOR over D_STEEP.  U3 basal.
D_DIVE   = 150.0         # m, DEPTH into the U1 body (from its polygon edge) over which U1 dives (parabolic) to U1_BASE; ~150 = half-width of the ~300 m-wide U1, so it reaches −10 m at the valley centre
U1_BASE  = -10.0         # m a.s.l., U1 (alluvium) bottom floor (reached D_DIVE from the U1/U2 contact)
FLOOR_HI = 0.0           # m a.s.l., U2 (terraces) bottom = FLAT FLOOR at +1 m (FLOOR_HI=FLOOR_LO); keeps the
FLOOR_LO = 0.0           #   top of U3 below the DTM under the cover, so U3 never outcrops within the terraces.
D_STEEP  = 900.0         # m, U2 parabolic dive distance: contact(~2100 m) -> +1 m floor reached at ~3000 m
D_FLOOR  = 3000.0        # m, (U2 floor slope, inactive while FLOOR_HI==FLOOR_LO)
AQ_BASE  = -35.0         # m a.s.l., base of unit 3 (model bottom)
DS = 20                  # DEM decimation 0.5 m -> 10 m working grid
KH = {"U1": 25.0, "U2": 40.0, "U3": 1.5}   # m/d (for the grid K columns; terraces repr. of 16-77)

# ---- 10 m DEM grid ----
with rasterio.open(DEM) as src:
    nd = src.nodata; W, H = src.width, src.height; ow, oh = W // DS, H // DS
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.bilinear).astype(float)
    tr = src.transform * Affine.scale(W / ow, H / oh); crs = src.crs
cell = tr.a
dem[(dem >= 1e30) | (dem == nd)] = np.nan
print(f"grid {dem.shape} @ {cell:.1f} m ; DEM {np.nanmin(dem):.1f}..{np.nanmax(dem):.1f} m a.s.l.")

# ---- rasterize the 3 hydrostrat units SEPARATELY (they overlap: U1 over U2 in the north) ----
g = gpd.read_file(GPKG, layer=HYDRO).to_crs(crs)
def rast(code):
    sub = g[g["Codigo"].astype(str) == code]
    return rasterize([(geom, 1) for geom in sub.geometry], out_shape=dem.shape,
                     transform=tr, fill=0, dtype="uint8").astype(bool)
p1, p2, p3 = rast("1"), rast("2"), rast("3")
p3 = p3 | ~(p1 | p2 | p3)                                    # gaps -> formation (basal default)
print(f"   present: U1 {100*p1.mean():.0f}% | U2 {100*p2.mean():.0f}% | U3-outcrop {100*p3.mean():.0f}% "
      f"| U1&U2 overlap {100*(p1 & p2).mean():.1f}%")

# ---- distances: d23 = dist from the U3 outcrop (for the U2 dive) ; d_u1 = DEPTH INTO the U1 polygon from its edge ----
d23  = distance_transform_edt(~p3) * cell                   # m to the nearest U3 outcrop (= U2/U3 contact)
# U1 alluvium OVERLAPS U2 (alluvium over terraces) over ~7% of U1, so the dive measures DEPTH INTO the U1 body from its
# OWN edge — distance_transform_edt(p1) — NOT distance from U2 (which is 0 throughout the overlap and pinched U1 to zero
# thickness across it).  Gaussian-smoothed so the bottom doesn't wiggle where the transect weaves against the ragged edge.
d_u1 = gaussian_filter(distance_transform_edt(p1) * cell, sigma=10)

# ---- U1 (alluvium) bottom: a SMOOTH parabolic wedge — 0 thickness at the U1 polygon edge, diving to U1_BASE a distance
#      D_DIVE into the body — anchored on a SMOOTHED valley-floor DTM (not the jagged terrace edge) so it stays clean. ----
_demf = dem.copy()                                          # fill NaN (nearest valid) so the Gaussian doesn't bleed NaN inward
if np.isnan(_demf).any():
    _, _j = distance_transform_edt(np.isnan(_demf), return_indices=True); _demf = _demf[_j[0], _j[1]]
dem_smooth = gaussian_filter(_demf, sigma=12)              # ~120 m low-pass -> smooth valley-floor anchor (no terrace roughness)
r1 = np.clip(d_u1 / D_DIVE, 0.0, 1.0)                       # 0 at the U1 edge -> 1 a distance D_DIVE into the U1 body
bot1 = np.where(p1, np.minimum(U1_BASE + (dem_smooth - U1_BASE) * (1.0 - r1) ** 2, dem), np.nan)

# ---- U2 (terraces) bottom: DTM (0 thickness) at the U2/U3 contact, dives PARABOLICALLY to FLOOR over D_STEEP ----
# (d23 = distance from the U3 outcrop is computed above, shared with the U1 down-catchment dive)
floor = np.clip(FLOOR_HI - (FLOOR_HI - FLOOR_LO) * (d23 / D_FLOOR), FLOOR_LO, FLOOR_HI)
r2 = np.clip(d23 / D_STEEP, 0.0, 1.0)                        # 0 at the contact -> 1 after D_STEEP
# parabola tangent (flat) at the floor: DTM at the contact (r2=0), eases onto FLOOR at r2=1
bot2 = np.minimum(floor + (dem - floor) * (1.0 - r2) ** 2, dem)

# ---- assemble the unit surfaces (m a.s.l.) ----
top1 = np.where(p1, dem, np.nan)
thk1 = top1 - bot1                                           # U1 thickness = DTM - diving bottom
top2 = np.where(p2, np.where(p1, bot1, dem), np.nan)         # under U1 overlap: = U1 bottom ; else DTM (outcrop)
bot2_u = np.where(p2, np.minimum(bot2, top2), np.nan)        # ensure bottom <= top
thk2 = np.where(p2, np.clip(top2 - bot2_u, 0.0, None), np.nan)
top3 = np.where(p2, bot2_u, np.where(p1, bot1, dem))         # under U2 -> U2 bottom ; under U1 -> U1 bottom ; else DTM
bot3 = np.full_like(dem, AQ_BASE)
thk3 = top3 - bot3
surfaces = {
    "top_U1_alluvium": top1, "bot_U1_alluvium": bot1, "thk_U1_alluvium": thk1,
    "top_U2_terraces": top2, "bot_U2_terraces": bot2_u, "thk_U2_terraces": thk2,
    "top_U3_pliomioc": top3, "bot_U3_pliomioc": bot3, "thk_U3_pliomioc": thk3,
}

# ---- domain mask (watershed) + write GeoTIFFs ----
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(crs)
wsgeom = ws.geometry.union_all()
ws_mask = geometry_mask([wsgeom], out_shape=dem.shape, transform=tr, invert=True)
valid = ws_mask & np.isfinite(dem)
prof = dict(driver="GTiff", height=dem.shape[0], width=dem.shape[1], count=1, dtype="float32",
            crs=crs, transform=tr, nodata=np.float32(np.nan), compress="deflate")
for nm, arr in surfaces.items():
    with rasterio.open(os.path.join(LAYDIR, nm + ".tif"), "w", **prof) as dst:
        dst.write(np.where(valid, arr, np.nan).astype("float32"), 1)
print("wrote", len(surfaces), "GeoTIFFs to", LAYDIR)
for u, nm in [("U1", "thk_U1_alluvium"), ("U2", "thk_U2_terraces"), ("U3", "thk_U3_pliomioc")]:
    a = surfaces[nm][valid]; a = a[np.isfinite(a)]
    if a.size:
        print(f"  {nm}: {a.min():5.1f}..{a.max():5.1f} m (mean {a.mean():4.1f})")

# ---- sample onto the Voronoi (DISV) grid centroids -> MODFLOW top + botm[3,ncpl] ----
gp = pickle.load(open(os.path.join(WS, "voronoi_grid.pkl"), "rb"))[0]
ncpl = int(gp["ncpl"]); cell2d = gp["cell2d"]
xc = np.array([c[1] for c in cell2d], float); yc = np.array([c[2] for c in cell2d], float)
def samp(arr, fill_nan=True):
    rr, cc = rowcol(tr, xc, yc); rr = np.clip(rr, 0, arr.shape[0]-1); cc = np.clip(cc, 0, arr.shape[1]-1)
    v = arr[np.asarray(rr), np.asarray(cc)]
    return v
g_top = samp(dem)
# botm[0]=bot of U1 (= U1 bottom where U1 present, else = top); botm[1]=top of U3; botm[2]=AQ_BASE
_sb1 = samp(bot1)                                           # U1 bottom (nan where no U1)
g_b0 = np.where(np.isfinite(_sb1), _sb1, g_top)             # botm[0] = U1 (parabolic) bottom, or top where no U1
g_b1 = samp(top3)                                            # top of unit 3 = base of the cover
g_b2 = np.full(ncpl, AQ_BASE)
fin = np.isfinite(g_top)
if (~fin).any():
    _, j = cKDTree(np.c_[xc[fin], yc[fin]]).query(np.c_[xc[~fin], yc[~fin]])
    for a in (g_top, g_b0, g_b1): a[~fin] = a[fin][j]
g_b1 = np.minimum(g_b1, g_b0 - 0.01); g_b2 = np.minimum(g_b2, g_b1 - 0.01)   # strictly decreasing
botm = np.vstack([g_b0, g_b1, g_b2])
np.savez(os.path.join(LAYDIR, "voronoi_layers.npz"), ncpl=ncpl, xc=xc, yc=yc, top=g_top, botm=botm,
         names=np.array(["U1_alluvium", "U2_terraces", "U3_pliomioc"]), aq_base=AQ_BASE)
with open(os.path.join(LAYDIR, "voronoi_layers.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["icell", "xc", "yc", "top", "bot_U1", "bot_U2", "bot_U3",
                                   "thk_U1", "thk_U2", "thk_U3"])
    for i in range(ncpl):
        t, b0, b1, b2 = g_top[i], botm[0, i], botm[1, i], botm[2, i]
        w.writerow([i, f"{xc[i]:.2f}", f"{yc[i]:.2f}", f"{t:.2f}", f"{b0:.2f}", f"{b1:.2f}", f"{b2:.2f}",
                    f"{t-b0:.2f}", f"{b0-b1:.2f}", f"{b1-b2:.2f}"])
print(f"sampled {ncpl} DISV cells -> voronoi_layers.npz/.csv ; grid thk U1 {(g_top-botm[0]).max():.1f} | "
      f"U2 {(botm[0]-botm[1]).max():.1f} | U3 {(botm[1]-botm[2]).max():.1f} m")

# ============================ maps: top / bottom / thickness, per unit ============================
ext = [tr.c, tr.c + tr.a*dem.shape[1], tr.f + tr.e*dem.shape[0], tr.f]
gb = g.total_bounds; wsS = gpd.GeoSeries([wsgeom])
pf = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
tx, ty = pf["xs"], pf["ys"]
def panel(ax, arr, title, cmap, kind):
    a = np.where(valid, arr, np.nan)
    im = ax.imshow(a, extent=ext, origin="upper", cmap=cmap)
    wsS.boundary.plot(ax=ax, edgecolor="red", lw=0.8); ax.plot(tx, ty, "k-", lw=0.8)
    ax.set_xlim(gb[0], gb[2]); ax.set_ylim(gb[1], gb[3]); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=9); plt.colorbar(im, ax=ax, shrink=0.78, label="m")
rows = [("U1 alluvium", "top_U1_alluvium", "bot_U1_alluvium", "thk_U1_alluvium"),
        ("U2 terraces", "top_U2_terraces", "bot_U2_terraces", "thk_U2_terraces"),
        ("U3 Plio-Miocene", "top_U3_pliomioc", "bot_U3_pliomioc", "thk_U3_pliomioc")]
fig, axs = plt.subplots(3, 3, figsize=(14, 13))
for r, (lab, tnm, bnm, knm) in enumerate(rows):
    panel(axs[r, 0], surfaces[tnm], f"{lab} — TOP (m a.s.l.)", "viridis", "elev")
    panel(axs[r, 1], surfaces[bnm], f"{lab} — BOTTOM (m a.s.l.)", "viridis", "elev")
    panel(axs[r, 2], surfaces[knm], f"{lab} — THICKNESS (m)", "YlGnBu", "thk")
fig.suptitle("CdL/DRYAD — hydrostratigraphic layers: top / bottom / thickness  (modified U1∩U2; U3 base −35 m)",
             fontsize=13, fontweight="bold")
fig.savefig(os.path.join(OUT, "layers_maps.png"), dpi=140, bbox_inches="tight")
print("saved layers_maps.png")

# ============================ validation: surfaces along the transect ============================
rr, cc = rowcol(tr, tx, ty); rr = np.clip(rr, 0, dem.shape[0]-1); cc = np.clip(cc, 0, dem.shape[1]-1)
sdist = pf["dists"]
p_dem = dem[rr, cc]; p_t3 = top3[rr, cc]; p_b3 = bot3[rr, cc]
p_b1 = np.where(np.isfinite(bot1[rr, cc]), bot1[rr, cc], p_dem)
# same layout/frame as cross_section_geology.png + cross_section_hydrostrat.png: cross-section (6) + a
# thin surface-geology strip (0.6) below, figsize 16.5x9, so the three figures are directly comparable.
fig2 = plt.figure(figsize=(16.5, 9.0))
gs2 = fig2.add_gridspec(2, 1, height_ratios=[6, 0.6], hspace=0.10)
ax = fig2.add_subplot(gs2[0]); axs = fig2.add_subplot(gs2[1], sharex=ax)
ax.fill_between(sdist, p_b3, p_t3, color="#fdf6a8", label="U3 Plio-Miocene")
ax.fill_between(sdist, p_t3, p_b1, color="#6b6b6b", label="U2 terraces")
ax.fill_between(sdist, p_b1, p_dem, color="#d9d9d9", label="U1 alluvium")
ax.plot(sdist, p_dem, "k-", lw=1.4); ax.plot(sdist, p_t3, color="0.3", lw=1.0, label="top of U3 (= base of cover)")
ax.axhline(AQ_BASE, color="0.1", ls=":", lw=1.0)
ax.set_ylim(AQ_BASE, np.nanmax(p_dem) + 8.0); ax.set_ylabel("Elevation (m a.s.l.)")
ax.set_title("Built layer surfaces along the transect (U1 alluvium → −10 m lens, parabolic; U2 terraces → 0 m floor, parabolic; U3 basal −35 m)")
ax.legend(loc="upper right", fontsize=8, ncol=2); ax.tick_params(labelbottom=False)
# --- surface-geology base strip = hydrostrat units (U1/U2/U3) along the transect (same colours as the hydrostrat section) ---
HCOL = {"u1": "#d9d9d9", "u2": "#6b6b6b", "u3": "#fdf6a8"}; ULBL = {"u1": "U1", "u2": "U2", "u3": "U3"}
surf_unit = np.where(p1[rr, cc], "u1", np.where(p2[rr, cc], "u2", "u3"))   # surface outcrop unit along the transect
edges = [0] + [k for k in range(1, len(sdist)) if surf_unit[k] != surf_unit[k - 1]] + [len(sdist)]
axs.set_ylim(0, 1)
for a, b in zip(edges[:-1], edges[1:]):
    d0 = sdist[a]; d1 = sdist[b] if b < len(sdist) else sdist[-1]; u = surf_unit[a]
    axs.fill_betweenx([0, 1], d0, d1, color=HCOL[u])
    if d1 - d0 > 120:
        axs.text(0.5 * (d0 + d1), 0.5, ULBL[u], ha="center", va="center", fontsize=9, fontweight="bold")
axs.set_yticks([]); axs.set_ylabel("surface\ngeology", fontsize=7.5, rotation=0, ha="right", va="center")
axs.set_xlabel("Distance along transect (m)   —   SE → NW", labelpad=8); axs.set_xlim(0, sdist.max())
fig2.savefig(os.path.join(OUT, "layers_check.png"), dpi=140, bbox_inches="tight")
print("saved layers_check.png\nDONE")
