"""
Land-cover ZONE map for the CdL PEST++ calibration — groups the COS-2025 majority
class per Voronoi cell into a handful of vegetation zones.  Both the MODEL actual-ET
(UZF) and the MODIS/Copernicus RS-ET are later aggregated to these zones (monthly
zone-mean AET) to form the ET observation group.

Outputs (in CdL_pest): cos_zones.npz (cos_class, zone_id per cell), .csv, .png
"""
import pickle, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter
import rasterio
from rasterio.mask import mask as rmask
from shapely.geometry import Polygon, Point
from flopy.discretization import VertexGrid

WORKSPACE = Path(r"E:\00code_ws\DRYAD\CdL_model")
PEST_DIR  = Path(r"E:\00code_ws\DRYAD\CdL_pest"); PEST_DIR.mkdir(parents=True, exist_ok=True)
COS_TIF   = r"E:\ArcGis_Data\WorkSpace\COS\COSc2025\COSc_2025_N3_v0_TM06.tif"

# COS N3 class -> (zone_id, zone_name).  Grouping mirrors the extdp rooting-depth classes.
ZONE_DEF = {
    311: (1, "oak_broadleaf"), 312: (1, "oak_broadleaf"), 313: (1, "oak_broadleaf"),
    321: (2, "pine"), 322: (2, "pine"), 323: (2, "pine"),
    410: (3, "matos"),
    420: (4, "grass_crops"), 211: (4, "grass_crops"), 212: (4, "grass_crops"), 213: (4, "grass_crops"),
    500: (5, "bare_built_water"), 100: (5, "bare_built_water"),
    610: (5, "bare_built_water"), 620: (5, "bare_built_water"),
}
ZONE_DEFAULT = (4, "grass_crops")     # unmapped -> grass-like
ZONE_NAMES = {0: "unclassified", 1: "oak_broadleaf", 2: "pine", 3: "matos",
              4: "grass_crops", 5: "bare_built_water"}

gp = pickle.load(open(WORKSPACE / "voronoi_grid.pkl", "rb"))[0]
vgrid = VertexGrid(**gp, nlay=1)
ncpl = vgrid.ncpl
xc = np.array(vgrid.xcellcenters).ravel(); yc = np.array(vgrid.ycellcenters).ravel()

cos_class = np.zeros(ncpl, dtype=int)
with rasterio.open(COS_TIF) as src:
    cnd = src.nodata
    for i in range(ncpl):
        poly = Polygon(vgrid.get_cell_vertices(i))
        try:
            out, _ = rmask(src, [poly], crop=True, all_touched=True, nodata=0)
            v = out[0].ravel(); v = v[(v != 0) & (v != (cnd or 0))]
            cos_class[i] = Counter(v.tolist()).most_common(1)[0][0] if v.size else 0
        except Exception:
            cos_class[i] = 0
    miss = np.where(cos_class == 0)[0]
    if miss.size:
        for v, i in zip(src.sample([(xc[i], yc[i]) for i in miss]), miss):
            cos_class[i] = int(v[0]) if v[0] not in (0, cnd) else 0

zone_id = np.array([ZONE_DEF.get(int(c), ZONE_DEFAULT)[0] for c in cos_class], dtype=int)

np.savez(PEST_DIR / "cos_zones.npz", cos_class=cos_class, zone_id=zone_id,
         zone_names=np.array([f"{k}:{v}" for k, v in ZONE_NAMES.items()]))
import pandas as pd
pd.DataFrame({"cell": np.arange(ncpl), "cos_class": cos_class, "zone_id": zone_id,
              "zone_name": [ZONE_NAMES[z] for z in zone_id]}).to_csv(PEST_DIR / "cos_zones.csv", index=False)

print("Zone cell counts:")
for z, n in sorted(Counter(zone_id.tolist()).items()):
    print(f"   {z} {ZONE_NAMES[z]:18s} {n:5d} cells ({100*n/ncpl:.1f}%)")

# map
fig, ax = plt.subplots(figsize=(8, 9))
cmap = matplotlib.colors.ListedColormap(["0.7", "#1b7837", "#7fbf7b", "#d9a066", "#f7f4b8", "#b0b0b0"])
pc = ax.scatter(xc, yc, c=zone_id, cmap=cmap, vmin=-0.5, vmax=5.5, s=6)
ax.set_aspect("equal"); ax.set_title("CdL land-cover zones (COS-2025 majority per cell) — for AET aggregation")
cb = fig.colorbar(pc, ticks=range(6)); cb.ax.set_yticklabels([ZONE_NAMES[i] for i in range(6)])
fig.tight_layout(); fig.savefig(PEST_DIR / "cos_zones.png", dpi=130, bbox_inches="tight")
print("wrote", PEST_DIR / "cos_zones.npz", "/ .csv / .png")
