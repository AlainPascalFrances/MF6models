"""
Post-correct the SYNTHETIC piezometer series (user 2026-08-14).  Reads the raw generator
output, applies two corrections, and OVERWRITES cdl_synthetic_piezo.csv (used by make_obs,
postprocess_cdl and ss_obs_vs_computed).  Re-runnable: it backs the raw up to
cdl_synthetic_piezo_raw.csv on first run and ALWAYS corrects from that backup.

Corrections (in this order):
  1. pre-1995 rise: for P0/P5/P6 the pre-1995 depths are reduced (water table raised) by 2 m
     to remove the early offset artefact -> a continuous series.
  2. depth rescale: each point's depth-to-water-table is linearly rescaled to a physically
     sensible [min, max] range (user-specified below).  wl_elev is recomputed = ground - depth
     (ground = wl_elev + depth, constant per point in the raw series).
"""
import shutil
from pathlib import Path
import numpy as np, pandas as pd

CSV = Path(r"E:\00code_ws\DRYAD\CdL_pest\snirh_data_availability\cdl_synthetic_piezo.csv")
RAW = CSV.with_name("cdl_synthetic_piezo_raw.csv")

# depth-to-water-table range [min, max] (m below ground) per piezometer
TARGET = {"P0": (1.0, 4.5), "P1": (1.0, 10.0), "P2": (1.0, 10.0), "P3": (2.0, 6.0),
          "P4": (0.5, 4.0), "P5": (2.0, 6.0), "P6": (1.0, 7.0)}   # 2026-08-15: P6 min 2->1
PRE1995_RISE = {"P0": 1.0, "P5": 2.0, "P6": 2.0}     # pre-1995 depth reduced by this (WT rises).
#   2026-08-15: P0 rise 2->1 (i.e. its pre-1995 WT is 1 m DEEPER than before). P5/P6 unchanged.
CUT = np.datetime64("1995-01-01")

if not RAW.exists():                                  # preserve the raw generator output once
    shutil.copy2(CSV, RAW)
    print(f">> backed up raw -> {RAW.name}")
df = pd.read_csv(RAW, parse_dates=["date"])

# MODEL land-surface elevation at each piezometer cell (from the MF6 grid). The depth targets are
# depth-below-the-MODEL-top (what obs_depth_timeseries.png and the calibration measure), NOT the
# synthetic-DEM ground (which differs by a few m), so recompute wl_elev = model_top - depth. This
# stops the corrected WT crossing the surface where model_top < synthetic ground. Fallback: synthetic.
import flopy, geopandas as gpd
from scipy.spatial import cKDTree
GRB = r"E:\00code_ws\DRYAD\CdL_model\cdl_gwf.disv.grb"
GPKG = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
OBS_LAYER = "obs_points_cdl"
model_ground = {}
try:
    _mg = flopy.mf6.utils.MfGrdFile(GRB).modelgrid
    _top = np.asarray(_mg.top, dtype=float).ravel()
    _tree = cKDTree(np.c_[np.array(_mg.xcellcenters).ravel(), np.array(_mg.ycellcenters).ravel()])
    _pts = gpd.read_file(GPKG, layer=OBS_LAYER).to_crs(3763)
    _ncol = "Name" if "Name" in _pts.columns else _pts.columns[0]
    for _, _r in _pts.iterrows():
        if _r.geometry is not None and not _r.geometry.is_empty:
            _, _idx = _tree.query([_r.geometry.x, _r.geometry.y])
            model_ground[str(_r[_ncol]).strip().upper()] = float(_top[_idx])
    print(">> model land surface at piezos: " + ", ".join(f"{k}={v:.1f}" for k, v in sorted(model_ground.items())))
except Exception as _e:
    print(f">> (model grid top not read: {_e!r}; using synthetic ground)")

out = []
for pt, g in df.groupby("cdl"):
    g = g.copy().sort_values("date")
    P = str(pt).upper()
    ground = model_ground.get(P, float((g["wl_elev_m"] + g["depth_m"]).mean()))   # MODEL top (fallback: synthetic)
    d = g["depth_m"].to_numpy(dtype=float).copy()
    # (1) pre-1995 rise (P0/P5/P6)
    if P in PRE1995_RISE:
        pre = g["date"].to_numpy() < CUT
        d[pre] -= PRE1995_RISE[P]
        print(f"   {P}: raised {int(pre.sum())} pre-1995 obs by {PRE1995_RISE[P]:.1f} m")
    # (2) linear rescale to the target depth range
    lo, hi = TARGET.get(P, (float(d.min()), float(d.max())))
    dmin, dmax = float(d.min()), float(d.max())
    d = lo + (hi - lo) * (d - dmin) / (dmax - dmin) if dmax > dmin else np.full_like(d, 0.5 * (lo + hi))
    g["depth_m"] = d
    g["wl_elev_m"] = ground - d
    out.append(g)

res = pd.concat(out).sort_values(["cdl", "date"])
res.to_csv(CSV, index=False)
print(f">> wrote {CSV.name}")
for pt, g in res.groupby("cdl"):
    print(f"   {pt}: depth {g.depth_m.min():.2f}-{g.depth_m.max():.2f} m | "
          f"wl_elev {g.wl_elev_m.min():.1f}-{g.wl_elev_m.max():.1f} m a.s.l.")
