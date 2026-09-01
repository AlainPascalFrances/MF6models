"""Layer diagnostics for the CdL MF6 model, built EXACTLY as the model does:
  - k_layers_map.png   : Kh per MF6 layer, ONE shared K colour ramp + a single colorbar
  - ibound_uzf_map.png : idomain (ibound) per MF6 layer + the UZF surface boundary (iuzfbnd)
The 3 MF6 numerical layers are NOT the 3 geologic units (see titles)."""
import config
import pickle, numpy as np, geopandas as gpd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import flopy
from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect
from shapely.geometry import MultiPolygon
from shapely.validation import make_valid

WS = Path(str(config.MODEL))
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
HYDRO = "dryad_modelo_nbs__gc_35a_cdl_hydrostrat"
VL = WS / "conceptual" / "layers" / "voronoi_layers.npz"
POND_SUBSET = [6, 12, 10, 18]                 # must match the model
POND_DEPTH = 4.0
KH_ALLU, KH_TERR, KH_FORM = 25.0, 40.0, 1.5

# --- grid ---
gp = pickle.load(open(WS / "voronoi_grid.pkl", "rb"))[0]
vg = VertexGrid(**gp, nlay=1); ncpl = vg.ncpl
xc = np.array(vg.xcellcenters); yc = np.array(vg.ycellcenters)
nlay = 3

# --- K arrays (point-in-polygon surfK, as in the model) ---
vl = np.load(VL)
top, botm = np.asarray(vl["top"], float), np.asarray(vl["botm"], float)
thk_cover = np.maximum(top - botm[0], 0.0) + np.maximum(botm[0] - botm[1], 0.0)
_h = gpd.read_file(GPKG, layer=HYDRO).to_crs(3763); _h["Codigo"] = _h["Codigo"].astype(str)
_pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xc, yc), crs=3763)
_sj = gpd.sjoin(_pts, _h[["Codigo", "geometry"]], how="left", predicate="within")
code = _sj[~_sj.index.duplicated(keep="first")].reindex(range(ncpl))["Codigo"].fillna("3").to_numpy()
surfK = np.select([code == "1", code == "2"], [KH_ALLU, KH_TERR], default=KH_FORM).astype(float)
coverK = np.where(thk_cover > POND_DEPTH, surfK, KH_FORM)
k_array = np.vstack([surfK, coverK, np.full(ncpl, KH_FORM)])

# --- lake footprints -> idomain + UZF (as the model: excavate the top layer at lake cells) ---
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
ws_poly = ws.geometry.union_all()
ws_poly = max(ws_poly.geoms, key=lambda p: p.area) if isinstance(ws_poly, MultiPolygon) else ws_poly
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(3763).reset_index(drop=True)
ponds["geometry"] = ponds.geometry.apply(make_valid)
def largest(g): return max(g.geoms, key=lambda p: p.area) if isinstance(g, MultiPolygon) else g
try:
    ix = GridIntersect(vg, method="vertex")
except TypeError:
    try:
        ix = GridIntersect(vg, rtree=True)
    except TypeError:
        ix = GridIntersect(vg)
lake_nodes = []; pond_polys = []
for fid in POND_SUBSET:
    g = largest(largest(ponds.geometry.iloc[fid]).intersection(ws_poly))
    pond_polys.append((fid, g))
    lake_nodes += [int(c) for c in ix.intersect(g)["cellids"]]
lake_nodes = sorted(set(lake_nodes))
idomain = np.ones((nlay, ncpl), int); idomain[0, lake_nodes] = 0     # _nexc=1: only top layer excavated
iuzfbnd = np.ones(ncpl, int); iuzfbnd[lake_nodes] = 0                # UZF on every surface cell except lakes
print(f"lake footprint cells = {len(lake_nodes)} ; idomain L1 inactive = {(idomain[0]==0).sum()} ; "
      f"UZF surface cells = {(iuzfbnd==1).sum()}/{ncpl}")

LAB = ["MF6 Layer 1 — top 4 m (lake host)\nK = outcropping unit (U1/U2/U3)",
       "MF6 Layer 2 — Quaternary cover (U1/U2)\nformation (U3) where cover < 4 m",
       "MF6 Layer 3 — U3 Plio-Miocene formation"]

# ===== Figure 1: Kh per layer, ONE shared ramp + single colorbar =====
fig, axes = plt.subplots(1, 3, figsize=(18, 7.2))
ca = None
for L in range(3):
    pmv = flopy.plot.PlotMapView(modelgrid=vg, ax=axes[L])
    ca = pmv.plot_array(k_array[L], cmap="viridis", vmin=KH_FORM, vmax=KH_TERR)
    ws.boundary.plot(ax=axes[L], color="k", lw=0.8)
    axes[L].set_aspect("equal"); axes[L].set_title(LAB[L], fontsize=11)
    axes[L].set_xlabel("X (m)"); axes[L].set_ylabel("Y (m)" if L == 0 else "")
cb = fig.colorbar(ca, ax=axes, shrink=0.85, pad=0.02)
cb.set_label("Kh (m/d)  —  shared ramp"); cb.set_ticks([1.5, 25, 40])
cb.set_ticklabels(["1.5  (U3 formation)", "25  (U1 alluvium)", "40  (U2 terraces)"])
fig.suptitle("CdL — horizontal K per MF6 layer  (K33 = 0.1·Kh).  The 3 NUMERICAL layers are NOT the 3 "
             "geologic units:\nU1 alluvium + U2 terraces are merged into the 'cover'; U3 is the basal formation.",
             fontsize=12)
fig.savefig(WS / "conceptual" / "layers" / "k_layers_map.png", dpi=150, bbox_inches="tight")
print("saved k_layers_map.png")

# ===== Figure 2: idomain per layer + UZF surface boundary =====
bcmap = ListedColormap(["crimson", "0.85"]); bn = BoundaryNorm([-0.5, 0.5, 1.5], bcmap.N)
fig2, ax2 = plt.subplots(1, 4, figsize=(22, 6.8))
panels = [(idomain[0], "idomain — MF6 Layer 1\n(lakes excavated = inactive)"),
          (idomain[1], "idomain — MF6 Layer 2\n(fully active)"),
          (idomain[2], "idomain — MF6 Layer 3\n(fully active)"),
          (iuzfbnd, "uzf_iuzfbnd (surface)\n(UZF on all cells except lakes)")]
for ax, (arr, ttl) in zip(ax2, panels):
    pmv = flopy.plot.PlotMapView(modelgrid=vg, ax=ax)
    im = pmv.plot_array(arr, cmap=bcmap, norm=bn)
    ws.boundary.plot(ax=ax, color="k", lw=0.8)
    inact = np.where(arr == 0)[0]                                  # lake cells (tiny) -> mark them visibly
    ax.scatter(xc[inact], yc[inact], s=12, c="crimson", edgecolors="none", zorder=5)
    for fid, g in pond_polys:
        gpd.GeoSeries([g], crs=3763).boundary.plot(ax=ax, color="navy", lw=1.3, zorder=6)
        ax.annotate(f"FID {fid}", (g.centroid.x, g.centroid.y), color="navy", fontsize=7.5,
                    fontweight="bold", ha="center", va="bottom", xytext=(0, 4),
                    textcoords="offset points", zorder=7)
    ax.set_aspect("equal"); ax.set_title(ttl, fontsize=11); ax.set_xlabel("X (m)")
cb2 = fig2.colorbar(im, ax=ax2, shrink=0.8, pad=0.02, ticks=[0, 1])
cb2.set_ticklabels(["0 = inactive / no UZF (lake)", "1 = active / UZF"])
fig2.suptitle("CdL — idomain (ibound) per MF6 layer + UZF surface boundary  "
              "(POND_SUBSET = 4 lakes; only the top layer is excavated at lake cells)", fontsize=12)
fig2.savefig(WS / "conceptual" / "layers" / "ibound_uzf_map.png", dpi=150, bbox_inches="tight")
print("saved ibound_uzf_map.png")
