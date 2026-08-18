"""
MODIS MOD16A2GF actual-ET -> land-cover ZONE monthly means (mm), the RS target for
the PEST++ 'aet' observation group.  Aggregated to the SAME zones and area-weighting
as the model side (model_aet_zonal.py), so the two are directly comparable.

Reads HDF4 via rasterio/GDAL (needs libgdal-hdf4).  GDAL supplies the MODIS sinusoidal
geolocation, so we transform the Voronoi-cell centroids into that CRS, sample the ET
grid at each cell, area-weight per zone, and convert 8-day composites to monthly totals.
Output: CdL_pest\mod16_aet_zonal.csv
"""
import glob, re, pickle
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from rasterio.transform import rowcol
from pyproj import Transformer
from flopy.discretization import VertexGrid
from shapely.geometry import Polygon

WS   = Path(r"E:\00code_ws\DRYAD\CdL_model")
PEST = Path(r"E:\00code_ws\DRYAD\CdL_pest")
MODDIR = Path(r"Y:\RS\MODIS4CDL")
GRID = "MOD_Grid_MOD16A2"
ZONE_NAMES = {1: "oak_broadleaf", 2: "pine", 3: "matos", 4: "grass_crops", 5: "bare_built_water"}
def sub(f):                                   # ET_500m subdataset path for a granule
    return f"HDF4_EOS:EOS_GRID:{f}:{GRID}:ET_500m"

# --- grid, zones, cell areas -------------------------------------------------
gp = pickle.load(open(WS / "voronoi_grid.pkl", "rb"))[0]
vgrid = VertexGrid(**gp, nlay=1); ncpl = vgrid.ncpl
xc = np.array(vgrid.xcellcenters).ravel(); yc = np.array(vgrid.ycellcenters).ravel()
area = np.array([Polygon(vgrid.get_cell_vertices(i)).area for i in range(ncpl)])
zone_id = np.load(PEST / "cos_zones.npz")["zone_id"]
zones = sorted(z for z in np.unique(zone_id) if z != 0)

granules = sorted(glob.glob(str(MODDIR / "MOD16A2GF*.hdf")))
print(f">> {len(granules)} granules")

# geolocation (identical for all granules of this tile) -> cell row/col + read window
with rasterio.open(sub(granules[0])) as r:
    tr, crs = r.transform, r.crs
rows, cols = rowcol(tr, *Transformer.from_crs(3763, crs.to_wkt(), always_xy=True).transform(xc, yc))
rows, cols = np.asarray(rows), np.asarray(cols)
inb = (rows >= 0) & (rows < 2400) & (cols >= 0) & (cols < 2400)
print(f"   {inb.sum()}/{ncpl} cells inside tile; row {rows[inb].min()}-{rows[inb].max()}, "
      f"col {cols[inb].min()}-{cols[inb].max()}")
r0, c0 = int(rows[inb].min()), int(cols[inb].min())
nr, nc = int(rows[inb].max()) - r0 + 1, int(cols[inb].max()) - c0 + 1

recs = []
for f in granules:
    mm = re.search(r"\.A(\d{4})(\d{3})\.", f); yr, doy = int(mm.group(1)), int(mm.group(2))
    start = datetime(yr, 1, 1) + timedelta(days=doy - 1)
    ndays = min(8, (datetime(yr, 12, 31) - start).days + 1)      # last composite of year is short
    with rasterio.open(sub(f)) as r:
        a = r.read(1, window=Window(c0, r0, nc, nr)).astype(float)
    a[(a < -32767) | (a > 32700)] = np.nan                       # mask fill/QC codes
    a *= 0.1                                                      # -> mm/8day
    etc = np.full(ncpl, np.nan)
    etc[inb] = a[rows[inb] - r0, cols[inb] - c0]
    zvals = {}
    for z in zones:
        msk = (zone_id == z) & np.isfinite(etc)
        if msk.sum() > 0:
            zvals[ZONE_NAMES[z]] = float(np.average(etc[msk], weights=area[msk]))
    recs.append((start, ndays, zvals))

# --- 8-day composites -> daily rate -> monthly totals ------------------------
zcols = [ZONE_NAMES[z] for z in zones]
day0 = min(r[0] for r in recs); day1 = max(r[0] + timedelta(days=r[1]) for r in recs)
days = pd.date_range(day0, day1 - timedelta(days=1), freq="D")
daily = pd.DataFrame(0.0, index=days, columns=zcols); cov = pd.Series(False, index=days)
for start, ndays, zvals in recs:
    idx = pd.date_range(start, periods=ndays, freq="D"); cov.loc[idx] = True
    for zc, val in zvals.items():
        daily.loc[idx, zc] += val / ndays                        # mm/day
monthly = daily.resample("MS").sum()
full = cov.resample("MS").sum().eq(cov.resample("MS").size())    # drop partial months
monthly[~full] = np.nan
monthly.index.name = "date"
monthly.to_csv(PEST / "mod16_aet_zonal.csv")
print(f">> wrote mod16_aet_zonal.csv ({monthly.shape[0]} months "
      f"{monthly.index.min().date()}..{monthly.index.max().date()})")
print("annual MODIS AET (mm/yr) per zone:")
print((monthly.resample("YS").sum().replace(0, np.nan).mean().round(0)).to_string())
