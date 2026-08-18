"""Produce the per-layer K (NPF) maps from voronoi_layers.npz, building the K arrays EXACTLY as the
MF6 model does (cdl_gwf_model_opusv1.py lines ~499-520).  Run after build_layers.py to inspect the K
layers from the corrected hydrostratigraphy.  -> <OUT>\\layers\\k_layers_map.png"""
import pickle, numpy as np, geopandas as gpd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import flopy
from flopy.discretization import VertexGrid

WS = Path(r"E:\00code_ws\DRYAD\CdL_model")
GPKG = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
VL = WS / "conceptual" / "layers" / "voronoi_layers.npz"
OUT = WS / "conceptual" / "layers" / "k_layers_map.png"

# --- model constants (must match cdl_gwf_model_opusv1.py) ---
POND_DEPTH = 4.0
KH_ALLU, KH_TERR, KH_FORM = 25.0, 40.0, 1.5

vl = np.load(VL); ncpl = int(vl["ncpl"])
top, botm = np.asarray(vl["top"], float), np.asarray(vl["botm"], float)
thk_cover = np.maximum(top - botm[0], 0.0) + np.maximum(botm[0] - botm[1], 0.0)

gp = pickle.load(open(WS / "voronoi_grid.pkl", "rb"))[0]
vg = VertexGrid(**gp, nlay=1)
xc = np.array(vg.xcellcenters); yc = np.array(vg.ycellcenters)

# surfK by POINT-IN-POLYGON on the hydrostrat outcrop (mirrors the model's K block)
_h = gpd.read_file(GPKG, layer="dryad_modelo_nbs__gc_35a_cdl_hydrostrat").to_crs(3763)
_h["Codigo"] = _h["Codigo"].astype(str)
_pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xc, yc), crs=3763)
_sj = gpd.sjoin(_pts, _h[["Codigo", "geometry"]], how="left", predicate="within")
_code = _sj[~_sj.index.duplicated(keep="first")].reindex(range(ncpl))["Codigo"].fillna("3").to_numpy()
allu = _code == "1"; terr = _code == "2"
surfK = np.select([allu, terr], [KH_ALLU, KH_TERR], default=KH_FORM).astype(float)   # L1 = outcrop K
coverK = np.where(thk_cover > POND_DEPTH, surfK, KH_FORM)           # L2 (cover where it extends below the pond)
k_array = np.vstack([surfK, coverK, np.full(ncpl, KH_FORM)])        # L3 = formation
k33 = 0.1 * k_array
unit = np.where(allu, 1, np.where(terr, 2, 3))                      # surface unit code
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
print(f"ncpl={ncpl}; surfK unit counts U1/U2/U3 = "
      f"{int((unit==1).sum())}/{int((unit==2).sum())}/{int((unit==3).sum())}")
print(f"Kh layer1 range {surfK.min():.1f}-{surfK.max():.1f}; layer2 {coverK.min():.1f}-{coverK.max():.1f}")

lname = ["L1 pond host", "L2 cover", "L3 formation"]
fig, axes = plt.subplots(2, 3, figsize=(17, 11), squeeze=False)
for L in range(3):
    for r, (arr, lbl) in enumerate([(k_array[L], "Kh"), (k33[L], "K33")]):
        ax = axes[r, L]
        pmv = flopy.plot.PlotMapView(modelgrid=vg, ax=ax)
        ca = pmv.plot_array(arr, cmap="cividis", vmin=KH_FORM * (0.1 if r else 1), vmax=KH_TERR * (0.1 if r else 1))
        ws.boundary.plot(ax=ax, color="k", lw=0.8)
        fig.colorbar(ca, ax=ax, shrink=0.7)
        ax.set_aspect("equal"); ax.set_title(f"{lname[L]} — {lbl} (m/d)", fontsize=11)
fig.suptitle("CdL — aquifer K per layer from the CORRECTED hydrostratigraphy "
             "(K: alluvium 25 / terraces 40 / formation 1.5 m/d)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"saved {OUT}")
