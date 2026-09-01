"""Scatter of OBSERVED vs COMPUTED steady-state heads at the 7 piezometers (1:1 line + stats).
Reads spinup_heads.npy (the SS water table) + the SYNTHETIC piezo targets (SNIRH-regionalised,
same source as obs_head_timeseries.png).  Saves into the model's input folder (_input/<stamp>/)."""
import config
import pickle, re, numpy as np, pandas as pd, geopandas as gpd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flopy.discretization import VertexGrid
from scipy.spatial import cKDTree

WS = Path(str(config.MODEL))
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
# OBSERVED heads = synthetic water-table elevations regionalised from SNIRH analogs; the
# qualitative Excel is kept only as a fallback if the synthetic csv is unavailable.
SYN_PIEZO_CSV = (str(config.PEST) + r"\snirh_data_availability\cdl_synthetic_piezo.csv")
XLSX = (str(config.MODEL) + r"\gis\WP3_modeling\3.1.1\modelo_numerico\piezos_qualitative_month_198101_202605.xlsx")
stamp = (WS / "last_run_stamp.txt").read_text().strip()
OUT = WS / "_input" / stamp        # match the model's PREPROC_DIR (= _input\<stamp>); NOT the old preproc\
OUT.mkdir(parents=True, exist_ok=True)

with open(WS / "voronoi_grid.pkl", "rb") as f:
    gvg, _ = pickle.load(f)
vg = VertexGrid(**gvg, nlay=1); ncpl = vg.ncpl
xc, yc = np.array(vg.xcellcenters), np.array(vg.ycellcenters)
head = np.load(WS / "spinup_heads.npy").reshape(-1, ncpl)   # (nlay, ncpl); nlay-agnostic via -1
head = np.where(np.abs(head) < 1e29, head, np.nan)
# Water table = highest head in each column (the model's own convention — §5b builds _pond_wt = nanmax over layers).
# FIX (2026-06-27): the old loop did `head[k] > botm[k]`, indexing the 3-row UNIT-bottom array (bU1/bU2/bU3 from the
# npz) by NUMERICAL layer k=0..nlay-1 -> IndexError once k>=3 (already broke at nlay=4; worse at the nlay=5 redesign).
# nanmax needs no per-layer cell bottoms, so it is correct for ANY nlay.
wt = np.nanmax(head, axis=0)
tree = cKDTree(np.c_[xc, yc])

# OBSERVED heads per point: synthetic wl_elev_m (mean = observed value, min/max = seasonal range).
obs = {}
if Path(SYN_PIEZO_CSV).exists():
    syn = pd.read_csv(SYN_PIEZO_CSV)
    for pt, g in syn.groupby("cdl"):
        s = pd.to_numeric(g["wl_elev_m"], errors="coerce").dropna()
        if len(s):
            obs[str(pt).upper()] = s
    _src = "synthetic, SNIRH-regionalised"
else:
    xl = pd.ExcelFile(XLSX)
    for sh in xl.sheet_names:
        nm = re.match(r"\s*(p\d+)", sh, re.I).group(1).upper()
        obs[nm] = pd.to_numeric(pd.read_excel(XLSX, sheet_name=sh)["head_monthly"], errors="coerce").dropna()
    _src = "qualitative xlsx"
pts = gpd.read_file(GPKG, layer="obs_points_cdl").to_crs(3763)

o_mean, o_lo, o_hi, c, names = [], [], [], [], []
for _, r in pts.iterrows():
    nm = str(r["Name"]).upper()
    if nm not in obs: continue
    _, n = tree.query([r.geometry.x, r.geometry.y])
    names.append(nm); c.append(wt[n])
    o_mean.append(obs[nm].mean()); o_lo.append(obs[nm].min()); o_hi.append(obs[nm].max())
o_mean, c = np.array(o_mean), np.array(c)
rmse = float(np.sqrt(np.mean((c - o_mean) ** 2))); bias = float(np.mean(c - o_mean))

fig, ax = plt.subplots(figsize=(7.2, 7))
lo, hi = min(o_mean.min(), c.min()) - 3, max(o_mean.max(), c.max()) + 3
ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
# observed seasonal range as horizontal error bars
ax.errorbar(o_mean, c, xerr=[o_mean - np.array(o_lo), np.array(o_hi) - o_mean], fmt="o",
            color="tab:brown", ms=8, capsize=3, ecolor="0.6", label="piezometers (obs range)")
for nm, ox, cy in zip(names, o_mean, c):
    ax.annotate(nm, (ox, cy), textcoords="offset points", xytext=(7, 4), fontsize=10, fontweight="bold")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel(f"OBSERVED head (m a.s.l., {_src})"); ax.set_ylabel("COMPUTED steady-state head (m a.s.l.)")
ax.set_title(f"CdL — SS observed vs computed heads  (recharge = 20% of P = 101.8 mm/yr)\n"
             f"RMSE = {rmse:.1f} m   bias = {bias:+.1f} m   (n={len(c)})")
ax.grid(alpha=0.3); ax.legend(loc="upper left")
fig.tight_layout()
out = OUT / "ss_obs_vs_computed_heads.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"RMSE={rmse:.2f} m  bias={bias:+.2f} m")
print(f"saved {out}")
