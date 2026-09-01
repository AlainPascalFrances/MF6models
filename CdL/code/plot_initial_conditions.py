"""Map the INITIAL CONDITIONS the transient lake run starts from = the no-lake steady-state
spin-up heads (spinup_heads.npy).  Plots water-table elevation, depth to water table, and the
per-layer head field on the Voronoi grid, and quantifies the lake-stage vs surrounding-water-table
mismatch that makes SP1 stiff.  -> WORKSPACE\\initial_conditions_map.png"""
import config
import pickle, numpy as np, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import flopy
from flopy.discretization import VertexGrid
from shapely.geometry import Point
from shapely.validation import make_valid

WORKSPACE = Path(str(config.MODEL))
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
VLAYERS = WORKSPACE / "conceptual" / "layers" / "voronoi_layers.npz"
HEADS = WORKSPACE / "spinup_heads.npy"
POND_SUBSET = [6, 12]            # the 2 modelled ponds (positional FID)
POND_DEPTH = 4.0                 # m (matches the model: lake bottom = ground - POND_DEPTH)

# --- grid ---
with open(WORKSPACE / "voronoi_grid.pkl", "rb") as f:
    gridprops_vg, gridprops_disv = pickle.load(f)
vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
xc, yc = np.array(vgrid.xcellcenters), np.array(vgrid.ycellcenters)

# --- surfaces + heads ---
vl = np.load(VLAYERS)
top = np.asarray(vl["top"], float)               # ground surface (ncpl,)
botm = np.asarray(vl["botm"], float)             # (3, ncpl) U1/U2/U3 bottoms
head = np.load(HEADS).reshape(-1, ncpl)          # (nlay, ncpl)
nlay = head.shape[0]
head = np.where(np.abs(head) < 1e29, head, np.nan)

# --- water table = head in the topmost SATURATED layer (head above that layer's bottom) ---
wt = np.full(ncpl, np.nan)
for k in range(nlay):
    wet = np.isnan(wt) & (head[k] > botm[k])
    wt[wet] = head[k][wet]
wt = np.where(np.isnan(wt), head[0], wt)         # fallback: top-layer head
depth = top - wt                                  # >0 below ground; <0 head above ground (seepage)
print(f"water table: {np.nanmin(wt):.1f}..{np.nanmax(wt):.1f} m a.s.l. | "
      f"depth: {np.nanmin(depth):.1f}..{np.nanmax(depth):.1f} m "
      f"(median {np.nanmedian(depth):.2f} m; {100*np.mean(depth<0.05):.0f}% within 5 cm of surface)")

# --- context layers (EPSG:3763) ---
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(3763)
ponds_all = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(3763).reset_index(drop=True)

# --- lake-stage vs surrounding water-table mismatch (why SP1 shocks) ---
cells = gpd.GeoDataFrame(geometry=[Point(x, y) for x, y in zip(xc, yc)], crs=3763)
print("\n   pond | ground(mean) | init stage (g-0.5) | water table under pond | mismatch (stage - WT)")
pond_info = []
for fid in POND_SUBSET:
    poly = make_valid(ponds_all.geometry.iloc[fid])
    inside = cells.within(poly).values
    if not inside.any():                          # tiny pond: fall back to nearest cell
        inside = np.zeros(ncpl, bool); inside[np.argmin((xc - poly.centroid.x) ** 2 + (yc - poly.centroid.y) ** 2)] = True
    g = float(np.mean(top[inside])); wtp = float(np.nanmean(wt[inside])); stage = g - 0.5
    pond_info.append((fid, poly, g, wtp, stage))
    print(f"   FID {fid:>2} | {g:11.2f} | {stage:17.2f} | {wtp:21.2f} | {stage - wtp:+.2f} m")

# --- figure: 2x2 (WT elevation, depth-to-water, head L0, head bottom layer) ---
panels = [("Water-table elevation (m a.s.l.)", wt, "viridis", None),
          ("Depth to water table (m below ground)", depth, "Spectral_r",
           (np.nanpercentile(depth, 2), np.nanpercentile(depth, 98))),
          (f"Head — layer 1 (top/soil) (m)", head[0], "viridis", None),
          (f"Head — layer {nlay} (formation) (m)", head[-1], "viridis", None)]
fig, axes = plt.subplots(2, 2, figsize=(15, 15))
for ax, (ttl, arr, cmap, lim) in zip(axes.ravel(), panels):
    pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
    kw = dict(cmap=cmap)
    if lim: kw.update(vmin=lim[0], vmax=lim[1])
    ca = pmv.plot_array(arr, **kw)
    ws.boundary.plot(ax=ax, color="k", lw=0.9, zorder=4)
    streams.plot(ax=ax, color="#2BAEEA", lw=0.8, zorder=3)
    for fid, poly, *_ in pond_info:
        gpd.GeoSeries([poly], crs=3763).boundary.plot(ax=ax, color="navy", lw=1.6, zorder=5)
        ax.annotate(f"FID {fid}", (poly.centroid.x, poly.centroid.y), fontsize=8, fontweight="bold",
                    color="navy", ha="center", va="bottom", zorder=6)
    fig.colorbar(ca, ax=ax, shrink=0.72)
    ax.set_title(ttl, fontsize=11); ax.set_aspect("equal")
    ax.set_xlabel("X (m, EPSG:3763)"); ax.set_ylabel("Y (m)")
fig.suptitle("CdL — INITIAL CONDITIONS for the transient lake run\n"
             "(no-lake steady-state spin-up heads = spinup_heads.npy; transient SP1 starts from this)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = WORKSPACE / "initial_conditions_map.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nsaved {out}")
