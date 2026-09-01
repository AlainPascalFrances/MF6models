"""Diagnose the no-lake steady-state water table vs the (qualitative) observed piezometers.
Samples spinup_heads.npy at each obs point and tabulates model vs observed head & depth-to-water."""
import config
import pickle, re, numpy as np, pandas as pd, geopandas as gpd
from pathlib import Path
from flopy.discretization import VertexGrid
from scipy.spatial import cKDTree

WS = Path(str(config.MODEL))
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
XLSX = (str(config.MODEL) + r"\gis\WP3_modeling\3.1.1\modelo_numerico\piezos_qualitative_month_198101_202605.xlsx")

# grid + surfaces + SS heads
with open(WS / "voronoi_grid.pkl", "rb") as f:
    gvg, _ = pickle.load(f)
vg = VertexGrid(**gvg, nlay=1); ncpl = vg.ncpl
xc, yc = np.array(vg.xcellcenters), np.array(vg.ycellcenters)
vl = np.load(WS / "conceptual" / "layers" / "voronoi_layers.npz")
top, ub = np.asarray(vl["top"], float), np.asarray(vl["botm"], float)   # ub = [bot_U1, bot_U2, -35]
head = np.load(WS / "spinup_heads.npy").reshape(-1, ncpl)
head = np.where(np.abs(head) < 1e29, head, np.nan)
# reconstruct the model's 4-layer botm (HOST + U1 + U2 + U3), as in cdl_gwf_model_opusv1.py §5
POND_DEPTH, AQ_BASE, MIN_FLOOR = 4.0, -35.0, 0.1
botm = np.zeros((4, ncpl))
botm[0] = top - POND_DEPTH
botm[1] = np.minimum(ub[0], botm[0] - MIN_FLOOR)
botm[2] = np.minimum(ub[1], botm[1] - MIN_FLOOR)
botm[3] = np.minimum(np.minimum(ub[2], AQ_BASE), botm[2] - MIN_FLOOR)
wt = np.full(ncpl, np.nan)
for k in range(head.shape[0]):
    m = np.isnan(wt) & (head[k] > botm[k]); wt[m] = head[k][m]
wt = np.where(np.isnan(wt), head[0], wt)
tree = cKDTree(np.c_[xc, yc])

# observed (qualitative) targets
xl = pd.ExcelFile(XLSX)
obs = {}
for sh in xl.sheet_names:
    nm = re.match(r"\s*(p\d+)", sh, re.I).group(1).upper()
    d = pd.read_excel(XLSX, sheet_name=sh)
    obs[nm] = (pd.to_numeric(d["head_monthly"], errors="coerce"), pd.to_numeric(d["depth"], errors="coerce"))

pts = gpd.read_file(GPKG, layer="obs_points_cdl").to_crs(3763)
print(f"{'piez':>4} | {'model top':>9} {'mod head':>8} {'mod depth':>9} | "
      f"{'obs head (min-mean-max)':>24} {'obs depth':>11} | {'head misfit':>11}")
rows = []
for _, r in pts.iterrows():
    nm = str(r["Name"]).upper()
    d, n = tree.query([r.geometry.x, r.geometry.y])
    mtop, mh, mdep = top[n], wt[n], top[n] - wt[n]
    if nm in obs:
        oh, od = obs[nm][0].dropna(), obs[nm][1].dropna()
        ohm, odm = oh.mean(), od.mean()
        misfit = mh - ohm
        print(f"{nm:>4} | {mtop:9.1f} {mh:8.1f} {mdep:9.1f} | "
              f"{oh.min():6.1f}/{ohm:5.1f}/{oh.max():5.1f}{'':>6} {od.min():4.1f}-{od.max():<4.1f}  | {misfit:+11.1f}")
        rows.append((nm, mtop, mh, mdep, ohm, odm, misfit))
    else:
        print(f"{nm:>4} | {mtop:9.1f} {mh:8.1f} {mdep:9.1f} |  (no observed sheet)")
df = pd.DataFrame(rows, columns=["piez", "model_top", "model_head", "model_depth", "obs_head", "obs_depth", "head_misfit"])
print(f"\nhead misfit: mean={df.head_misfit.mean():+.1f} m  | mean ABS={df.head_misfit.abs().mean():.1f} m  "
      f"| range {df.head_misfit.min():+.1f}..{df.head_misfit.max():+.1f} m")
print(f"model depth mean={df.model_depth.mean():.1f} m vs observed depth mean={df.obs_depth.mean():.1f} m "
      f"-> model is {df.model_depth.mean()-df.obs_depth.mean():+.1f} m too deep on average")
print(f"\ninterpretation: misfit<0 => model head too LOW (water table too deep / over-drained) at that piezo")
