"""PROTOTYPE of the user's chosen layering: 4 MF6 layers = [pond/soil HOST (top POND_DEPTH) + U1 + U2 + U3].
  - HOST (layer 1): top POND_DEPTH m; K = surface OUTCROP unit (point-in-polygon, 25/40/1.5). Lakes excavate THIS.
  - U1/U2/U3 (layers 2-4): CONSTANT K 25/40/1.5; bottom = the unit bottom (capped below the host).  Where a
    unit-layer pinches out (thickness ~0 because the unit is absent or entirely inside the host) it is a
    vertical PASS-THROUGH (idomain = -1).  U3 (basal) is always active.
Produces proto_k_per_unit.png + proto_ibound_uzf.png (4 panels each) for validation, no model run."""
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
POND_SUBSET = [6, 12, 10, 18]; POND_DEPTH = 4.0; MIN_THK = 0.5; MIN_FLOOR = 0.1
KH_ALLU, KH_TERR, KH_FORM = 25.0, 40.0, 1.5
nlay = 4

gp = pickle.load(open(WS / "voronoi_grid.pkl", "rb"))[0]
vg = VertexGrid(**gp, nlay=1); ncpl = vg.ncpl
xc = np.array(vg.xcellcenters); yc = np.array(vg.ycellcenters)
vl = np.load(VL); top, ub = np.asarray(vl["top"], float), np.asarray(vl["botm"], float)   # ub = [botU1, botU2, -35]

# --- surface outcrop unit K (host) by point-in-polygon ---
_h = gpd.read_file(GPKG, layer=HYDRO).to_crs(3763); _h["Codigo"] = _h["Codigo"].astype(str)
_sj = gpd.sjoin(gpd.GeoDataFrame(geometry=gpd.points_from_xy(xc, yc), crs=3763), _h[["Codigo", "geometry"]],
                how="left", predicate="within")
code = _sj[~_sj.index.duplicated(keep="first")].reindex(range(ncpl))["Codigo"].fillna("3").to_numpy()
surfK = np.select([code == "1", code == "2"], [KH_ALLU, KH_TERR], default=KH_FORM).astype(float)

# --- 4-layer botm: host base = top-POND_DEPTH ; unit layers capped below it, strictly decreasing ---
botm = np.zeros((nlay, ncpl))
botm[0] = top - POND_DEPTH                                  # host base
botm[1] = np.minimum(ub[0], botm[0] - MIN_FLOOR)            # U1 layer base
botm[2] = np.minimum(ub[1], botm[1] - MIN_FLOOR)            # U2 layer base
botm[3] = np.minimum(ub[2], botm[2] - MIN_FLOOR)            # U3 layer base (-35)
thk = np.vstack([top - botm[0], botm[0] - botm[1], botm[1] - botm[2], botm[2] - botm[3]])

# --- lake footprints (excavate the HOST only) ---
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763); ws_poly = ws.geometry.union_all()
ws_poly = max(ws_poly.geoms, key=lambda p: p.area) if isinstance(ws_poly, MultiPolygon) else ws_poly
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(3763).reset_index(drop=True)
ponds["geometry"] = ponds.geometry.apply(make_valid)
def largest(g): return max(g.geoms, key=lambda p: p.area) if isinstance(g, MultiPolygon) else g
try: ix = GridIntersect(vg, method="vertex")
except TypeError:
    try: ix = GridIntersect(vg, rtree=True)
    except TypeError: ix = GridIntersect(vg)
lake_nodes, pond_polys = [], []
for fid in POND_SUBSET:
    g = largest(largest(ponds.geometry.iloc[fid]).intersection(ws_poly)); pond_polys.append((fid, g))
    lake_nodes += [int(c) for c in ix.intersect(g)["cellids"]]
lake_nodes = np.array(sorted(set(lake_nodes)))

# --- idomain: host=1 (0 at ponds) ; units = 1 active / -1 pass-through (pinched) ; U3 always active ---
idomain = np.ones((nlay, ncpl), int)
for L in (1, 2):
    idomain[L] = np.where(thk[L] > MIN_THK, 1, -1)
idomain[0, lake_nodes] = 0                                  # lakes excavate ONLY the host
iuzfbnd = np.ones(ncpl, int); iuzfbnd[lake_nodes] = 0       # UZF on the host, except lakes
KH = [surfK, np.full(ncpl, KH_ALLU), np.full(ncpl, KH_TERR), np.full(ncpl, KH_FORM)]
Kmask = np.array([np.where(idomain[L] == 1, KH[L], np.nan) for L in range(nlay)])
for L in range(nlay):
    nm = ["HOST (surface unit)", "U1 (25)", "U2 (40)", "U3 (1.5)"][L]
    print(f"L{L+1} {nm}: active {int((idomain[L]==1).sum())} | pass-through {int((idomain[L]==-1).sum())} "
          f"| pond {int((idomain[L]==0).sum())}")

LAB = ["L1 = HOST (top 4 m)\nK = surface outcrop unit", "L2 = U1 alluvium (K=25)",
       "L3 = U2 terraces (K=40)", "L4 = U3 formation (K=1.5)"]
def overlay(ax):
    for fid, g in pond_polys: gpd.GeoSeries([g], crs=3763).boundary.plot(ax=ax, color="navy", lw=1.1, zorder=6)

fig, ax = plt.subplots(1, 4, figsize=(21, 6.6)); ca = None
for L in range(nlay):
    pmv = flopy.plot.PlotMapView(modelgrid=vg, ax=ax[L]); ca = pmv.plot_array(Kmask[L], cmap="viridis", vmin=1.5, vmax=40)
    ws.boundary.plot(ax=ax[L], color="k", lw=0.7); overlay(ax[L]); ax[L].set_aspect("equal"); ax[L].set_title(LAB[L], fontsize=10)
cb = fig.colorbar(ca, ax=ax, shrink=0.82, pad=0.02, ticks=[1.5, 25, 40]); cb.set_ticklabels(["1.5 (U3)", "25 (U1)", "40 (U2)"]); cb.set_label("Kh (m/d)")
fig.suptitle("PROTOTYPE — 4 layers: HOST + U1 + U2 + U3 (constant K below; white = pass-through / pond)", fontsize=12)
fig.savefig(WS / "conceptual" / "layers" / "proto_k_per_unit.png", dpi=150, bbox_inches="tight"); print("saved proto_k_per_unit.png")

cmap = ListedColormap(["crimson", "#9ecae1", "0.7"]); norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
fig2, a2 = plt.subplots(1, nlay, figsize=(21, 6.6)); im = None
for L in range(nlay):
    pmv = flopy.plot.PlotMapView(modelgrid=vg, ax=a2[L]); im = pmv.plot_array(idomain[L], cmap=cmap, norm=norm)
    ws.boundary.plot(ax=a2[L], color="k", lw=0.7); overlay(a2[L]); a2[L].set_aspect("equal")
    a2[L].set_title("idomain — " + LAB[L].split("\n")[0], fontsize=10)
cb2 = fig2.colorbar(im, ax=a2, shrink=0.82, pad=0.02, ticks=[-1, 0, 1]); cb2.set_ticklabels(["-1 pass-through", "0 pond", "1 active"])
fig2.suptitle("PROTOTYPE — idomain per layer (ponds excavate ONLY the host; UZF/iuzfbnd = host minus lakes)", fontsize=12)
fig2.savefig(WS / "conceptual" / "layers" / "proto_ibound_uzf.png", dpi=150, bbox_inches="tight"); print("saved proto_ibound_uzf.png")
