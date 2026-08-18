"""Final CdL/DRYAD topographic cross-section with ERT P1-P4 and geological unit boundaries (Codigo).
v2: map shows all ERT P1-P6 + one label per geological unit; profile labels de-overlapped."""
import os, copy, numpy as np, geopandas as gpd, pandas as pd, rasterio
from rasterio import Affine
from rasterio.enums import Resampling
from shapely.geometry import LineString, Point, MultiPolygon
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LightSource
import matplotlib.patheffects as pe

OUT = r"E:\00code_ws\DRYAD\CdL_model\conceptual"
GPKG = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
DEM  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif"

# ---- transect + topography ----
d = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
dists, xs, ys, zz = d["dists"], d["xs"], d["ys"], d["zz"]
vert_dist = d["vert_dist"]; verts = d["verts"]
line = LineString([tuple(v) for v in verts])
ert_dist = {"P1": vert_dist[1], "P2": vert_dist[2], "P3": vert_dist[3], "P4": vert_dist[4]}

# ---- geology ----
g = gpd.read_file(GPKG, layer="gc_35a_cdl").to_crs(3763)
desc_map = g[["Codigo", "NomeOrigin"]].drop_duplicates().set_index("Codigo")["NomeOrigin"].to_dict()
diss = g.dissolve(by="Codigo", as_index=False)

# stratigraphic order: most-recent -> oldest (epoch from Idade4, tie-break by CodigoLU)
def _epoch_rank(v):
    v = str(v).lower()
    if "pleist" in v and "holoc" in v: return 2
    if "holoc" in v: return 1
    if "pleist" in v: return 3
    if "plioc" in v: return 4
    if "mioc" in v: return 5
    return 99
_age = g[["Codigo", "CodigoLU", "Idade4"]].drop_duplicates().copy()
_age["r"] = _age["Idade4"].map(_epoch_rank)
# within the same epoch: explicit (user-confirmed) order first, then CodigoLU.
# Pliocénico: PSA is younger / on top of PSM (confirmed by A. Frances, 2026-06-19).
STRAT_HINT = {"PSA": 0, "PSM": 1}
_age["h"] = _age["Codigo"].map(lambda c: STRAT_HINT.get(c, 0))
age_order = list(dict.fromkeys(_age.sort_values(["r", "h", "CodigoLU"])["Codigo"]))  # youngest -> oldest (deduped: 'a' merges 2 rows)

segs = []
for _, row in diss.iterrows():
    inter = line.intersection(row.geometry)
    if inter.is_empty:
        continue
    parts = list(inter.geoms) if inter.geom_type.startswith("Multi") else [inter]
    for p in parts:
        if p.geom_type != "LineString" or p.length == 0:
            continue
        ds = line.project(Point(p.coords[0])); de = line.project(Point(p.coords[-1]))
        d0, d1 = sorted([ds, de]); segs.append([d0, d1, row["Codigo"]])
segs.sort()
merged = []
for s in segs:
    if merged and merged[-1][2] == s[2] and abs(merged[-1][1] - s[0]) < 1.0:
        merged[-1][1] = s[1]
    else:
        merged.append(list(s))
segs = merged
boundaries = [((segs[i-1][1] + segs[i][0]) / 2.0, segs[i-1][2], segs[i][2]) for i in range(1, len(segs))]

# AML geology codes + official AML 100k symbology RGB (user-provided)
UNIT_COLORS = {"a": "#e1e1e1", "Qdae": "#d1d1c7", "Qt": "#bbbdb8", "Qi'": "#a6aba7",
               "PSA": "#ffffb7", "PSM": "#fff091", "MAT": "#ffff00"}
def col(c): return UNIT_COLORS.get(c, "#dddddd")
# short names to keep the in-box profile legend compact
SHORT = {"a": "Aluviões/aterros", "Qdae": "Dunas/areias eólicas", "Qt": "Terraços f./m.",
         "Qi'": "Cascalheiras", "PSA": "F. Serra de Almeirim", "PSM": "F. Santa Marta",
         "MAT": "F. Alcoentre–Tomar"}
units_order = list(dict.fromkeys(s[2] for s in segs))
halo = [pe.withStroke(linewidth=2.2, foreground="white")]

# ---- conceptual hydrostratigraphy (hydrostrigraphy.docx; Zbyszewski & Ferreira 1968 + ERT/piezos) ----
LAYER1_THK = 4.0           # uniform low-K top (soil + fine sands); doc: 3-5 m
THK3_MAX   = 25.0          # Layer-3 terrace wedge thickness at the NW end; pinches out to 0 at the SE contact
AQ_BASE    = -45.0         # flat base of the modelled unconfined aquifer (>=40-50 m, not reached)
WT_DEP_SE, WT_DEP_NW = 8.5, 1.5    # water-table depth: ~8-9 m (SE) -> near-surface (NW stream)
FORMATION  = {"PSM", "PSA"}        # low-K Plio-Miocene outcrop (SE) -> terrace pinch-out point
contact_d  = next((bd for bd, ca, cb in boundaries if ca in FORMATION and cb not in FORMATION),
                  dists.max() * 0.4)
HS = {"l1": "#e9e0c8", "l2": "#b7bd8e", "l3": "#f2c45a"}   # Layer1 top / Layer2 Plio-Mio (basal) / Layer3 terraces

BASE_T = -30.0             # plot floor for the blank template
BASE_H = -52.0             # plot floor for the interpretation (below the aquifer base)
zmin, zmax = float(np.nanmin(zz)), float(np.nanmax(zz))
def topo_at(dd): return float(np.interp(dd, dists, zz))

# vertically-stacked layer geometry (units overlap vertically; oldest/deepest at the base):
#   Layer 1 (top, uniform)  >  Layer 3 (terrace wedge, NW only)  >  Layer 2 (Plio-Miocene, basal everywhere)
top1 = zz - LAYER1_THK                                                    # base of Layer 1
thk3 = np.where(dists >= contact_d,
                THK3_MAX * np.clip((dists - contact_d) / (dists.max() - contact_d), 0, 1), 0.0)
top2 = top1 - thk3                                                        # base of Layer 3 = top of Layer 2 (deepens NW)
xs_units = [c for c in age_order if c in {s[2] for s in segs}]


def draw_strip(axs):
    axs.set_ylim(0, 1)
    for d0, d1, c in segs:
        axs.fill_betweenx([0, 1], d0, d1, color=col(c))
        if (d1 - d0) >= 110:
            axs.text(0.5 * (d0 + d1), 0.5, c, ha="center", va="center", fontsize=9, fontweight="bold")
    for bd, _, _ in boundaries:
        axs.axvline(bd, color="0.25", lw=0.8)
        axs.annotate(f"{bd:.0f}", (bd, 0), textcoords="offset points", xytext=(0, -3), rotation=90,
                     ha="center", va="top", fontsize=7, color="0.3", annotation_clip=False)
    axs.set_yticks([]); axs.set_ylabel("surface\ngeology", fontsize=7.5, rotation=0, ha="right", va="center")
    axs.set_xlabel("Distance along transect (m)   —   SE → NW", labelpad=14)


def draw_common(ax, base):
    for nm, dd in ert_dist.items():
        ee = topo_at(dd)
        ax.plot([dd, dd], [base, ee], color="crimson", lw=1.0, alpha=0.5, zorder=4)
        ax.plot(dd, ee, "v", color="crimson", ms=11, mec="black", mew=0.7, zorder=9)
        ax.annotate(f"ERT {nm}", (dd, ee), textcoords="offset points", xytext=(0, 13), ha="center",
                    fontsize=11, fontweight="bold", color="crimson", path_effects=halo, zorder=10)
    for dd, lab, ha_, xo in [(0, "START\n(SE divide)", "left", 8), (dists.max(), "END\n(NW limit)", "right", 2)]:
        ax.plot(dd, topo_at(dd), "o", color="black", ms=7, zorder=9)
        ax.annotate(lab, (dd, topo_at(dd)), textcoords="offset points", xytext=(xo, -14),
                    va="top", ha=ha_, fontsize=10, fontweight="bold", path_effects=halo, zorder=10)
    ax.plot(dists, zz, color="black", lw=1.9, zorder=6)
    ax.set_ylabel("Elevation (m a.s.l.)"); ax.grid(axis="y", alpha=0.25)
    ax.annotate("SE", (0.006, 0.985), xycoords="axes fraction", fontsize=13, fontweight="bold", va="top")
    ax.annotate("NW", (0.99, 0.985), xycoords="axes fraction", fontsize=13, fontweight="bold", ha="right", va="top")


def add_ve(fig, ax, base, ytop):
    fig.canvas.draw(); pos = ax.get_position()
    ve = (dists.max() / (pos.width * fig.get_figwidth())) / ((ytop - base) / (pos.height * fig.get_figheight()))
    ax.annotate(f"V.E. ≈ {ve:.0f}×", (0.5, 0.985), xycoords="axes fraction",
                ha="center", va="top", fontsize=9, style="italic", color="0.35")


# ===================== FIGURE 1a: blank topographic section (kept for hand-drawing) =====================
fig = plt.figure(figsize=(16.5, 9.0))
gs = fig.add_gridspec(2, 1, height_ratios=[6, 0.6], hspace=0.10)
ax = fig.add_subplot(gs[0]); axs = fig.add_subplot(gs[1], sharex=ax)
YTOP = zmax + 24
ax.set_ylim(BASE_T, YTOP); ax.set_xlim(0, dists.max())
for bd, _, _ in boundaries:                                              # faint surface-contact guides
    ax.plot([bd, bd], [BASE_T, topo_at(bd)], color="0.6", ls="--", lw=0.8, zorder=3)
draw_common(ax, BASE_T)
ax.set_title("CdL / DRYAD — Topographic cross-section (SE→NW): ERT & surface geology — blank section to draw the layers",
             fontsize=12.5, fontweight="bold")
handles = [Patch(fc=col(c), ec="0.6", label=f"{c} — {SHORT.get(c, desc_map.get(c, ''))}") for c in xs_units]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.012, 0.965), fontsize=8, framealpha=0.9,
          title="Surface geology (recent→old)", title_fontsize=8, labelspacing=0.25, borderpad=0.4, handlelength=1.2)
draw_strip(axs)
add_ve(fig, ax, BASE_T, YTOP)
fig.savefig(os.path.join(OUT, "cross_section_topography.png"), dpi=170, bbox_inches="tight")
print("saved cross_section_topography.png (blank section)")

# ---- calibrate the pixel<->section transform of THIS exact image (inject 2 markers, detect) ----
from PIL import Image as _PILImage
_cal = [(dists.max() * 0.08, BASE_T + 8.0), (dists.max() * 0.92, YTOP - 8.0)]   # (dist, elev) interior pts
for _xd, _yd in _cal:
    ax.plot(_xd, _yd, "s", ms=4, color=(1, 0, 1), mec=(1, 0, 1), zorder=30)
_tmp = os.path.join(OUT, "_cal_tmp.png")
fig.savefig(_tmp, dpi=170, bbox_inches="tight")
_im = np.asarray(_PILImage.open(_tmp).convert("RGB"))
_msk = (_im[:, :, 0] > 200) & (_im[:, :, 1] < 100) & (_im[:, :, 2] > 180)
_ry, _rx = np.where(_msk)
_mid = 0.5 * (_rx.min() + _rx.max()); _Lm = _rx < _mid
colL, rowL, colR, rowR = _rx[_Lm].mean(), _ry[_Lm].mean(), _rx[~_Lm].mean(), _ry[~_Lm].mean()
(_d1, _e1), (_d2, _e2) = _cal
_a = (colR - colL) / (_d2 - _d1); _b = colL - _a * _d1          # col  = a*dist + b
_c = (rowR - rowL) / (_e2 - _e1); _dd = rowL - _c * _e1         # row  = c*elev + d
np.savez(os.path.join(OUT, "section_transform.npz"), a=_a, b=_b, c=_c, d=_dd,
         BASE_T=BASE_T, YTOP=YTOP, dmax=float(dists.max()), imgW=_im.shape[1], imgH=_im.shape[0])
os.remove(_tmp); plt.close(fig)
print(f"transform: col={_a:.3f}*dist{_b:+.1f} ; row={_c:.4f}*elev{_dd:+.1f} ; cal img {_im.shape[1]}x{_im.shape[0]}")

# (FIGURE 1b - obsolete parametric hydrostratigraphy - REMOVED. geo_section.py now
#  produces cross_section_hydrostrat.png from the user's digitized shapefile, so this
#  block must not run here or it clobbers that file.)

# ======================= FIGURE 2: plan-view geology over hillshade =======================
# hillshade from the LiDAR DTM (decimated for a plot background)
DS = 6
with rasterio.open(DEM) as src:
    nd = src.nodata; W, H = src.width, src.height
    ow, oh = W // DS, H // DS
    dem = src.read(1, out_shape=(oh, ow), resampling=Resampling.nearest).astype(float)
    tr = src.transform * Affine.scale(W / ow, H / oh)
    db = src.bounds
nd_mask = (dem >= 1e30) | (dem == nd)
dem[nd_mask] = np.nan
dx, dy = tr.a, -tr.e
ls = LightSource(azdeg=315, altdeg=45)
hs = ls.hillshade(np.where(np.isnan(dem), np.nanmean(dem), dem), vert_exag=4.0, dx=dx, dy=dy)
hs[nd_mask] = np.nan
dem_extent = [db.left, db.right, db.bottom, db.top]
gray = copy.copy(plt.cm.gray); gray.set_bad(alpha=0.0)

fig2, ax2 = plt.subplots(figsize=(11, 12.5))
ax2.set_aspect("equal")
ax2.imshow(np.ma.masked_invalid(hs), extent=dem_extent, origin="upper", cmap=gray,
           vmin=0.15, vmax=1.0, zorder=0)
for c in sorted(g["Codigo"].unique()):
    g[g["Codigo"] == c].plot(ax=ax2, color=col(c), ec="0.4", lw=0.3, alpha=0.55, zorder=1)
# streams (light blue, white casing for legibility over the drape)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(3763)
streams.plot(ax=ax2, color="white", lw=3.0, alpha=0.7, zorder=2.8)
streams.plot(ax=ax2, color="#2BAEEA", lw=1.5, zorder=3.0)
# watershed
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
ws.boundary.plot(ax=ax2, edgecolor="blue", lw=1.3, zorder=4)
# cattle ponds (charcas): navy polygons + small FID labels (positional index = model FID; font < unit labels)
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(3763).reset_index(drop=True)
ponds.plot(ax=ax2, color="navy", ec="white", lw=0.4, alpha=0.9, zorder=5.4)
for i, geom in enumerate(ponds.geometry):
    rp = geom.representative_point()
    ax2.plot(rp.x, rp.y, marker="o", ms=3.2, color="navy", mec="white", mew=0.5, zorder=5.6)
    ax2.annotate(str(i), (rp.x, rp.y), textcoords="offset points", xytext=(3.5, 3.0),
                 fontsize=7, fontweight="bold", color="navy", path_effects=halo, zorder=10)
# ALL ERT profiles P1-P6 (red=on transect, violet=others) + labels
ert = gpd.read_file(GPKG, layer="ert_electrodes_cdl")
for e in sorted(ert["ert"].unique()):
    sub = ert[ert["ert"] == e].sort_values("id_electro")
    on = e in (1, 2, 3, 4)
    ax2.plot(sub.geometry.x, sub.geometry.y, "-", color=("crimson" if on else "darkviolet"),
             lw=3.2, zorder=6, solid_capstyle="round")
    cx, cy = sub.geometry.x.mean(), sub.geometry.y.mean()
    ax2.annotate(f"P{e}", (cx, cy), textcoords="offset points", xytext=(7, 7), fontsize=12,
                 fontweight="bold", color=("crimson" if on else "darkviolet"),
                 path_effects=halo, zorder=10)
# transect
ax2.plot(xs, ys, "-", color="black", lw=2.0, zorder=5)
ax2.plot(verts[0][0], verts[0][1], "o", color="black", ms=7, zorder=8)
ax2.annotate("START\n(SE divide)", (verts[0][0], verts[0][1]), textcoords="offset points",
             xytext=(8, -2), fontsize=10, fontweight="bold", path_effects=halo, zorder=10)
ax2.plot(verts[-1][0], verts[-1][1], "o", color="black", ms=7, zorder=8)
ax2.annotate("END\n(NW limit)", (verts[-1][0], verts[-1][1]), textcoords="offset points",
             xytext=(9, 4), fontsize=10, fontweight="bold", path_effects=halo, zorder=10)

# legends (kept inside the map)
handles2 = [Patch(fc=col(c), ec="0.5", label=f"{c} — {desc_map.get(c, '')}") for c in age_order]
leg2 = ax2.legend(handles=handles2, loc="lower left", fontsize=8.5,
                  title="Codigo — recent → old (top → bottom)")
ax2.add_artist(leg2)
leg_lines = ax2.legend(handles=[Line2D([0], [0], color="crimson", lw=3, label="ERT on transect (P1–P4)"),
                                Line2D([0], [0], color="darkviolet", lw=3, label="ERT other (P5, P6)"),
                                Line2D([0], [0], color="black", lw=2, label="Transect"),
                                Line2D([0], [0], color="#2BAEEA", lw=2, label="Streams"),
                                Line2D([0], [0], color="blue", lw=1.3, label="Watershed"),
                                Line2D([0], [0], marker="o", color="navy", mec="white", mew=0.5,
                                       ls="none", ms=6, label="Cattle ponds (charcas) — FID")],
                       loc="upper right", fontsize=8.5)
ax2.set_title("Geology gc_35a_cdl over LiDAR-DTM hillshade — ERT P1–P6, transect & watershed")
# hillshade credit: compact boxed note in the bottom-right corner (kept clear of the legend)
hs_ann = ax2.annotate("Hillshade: LiDAR DTM 0.5 m\nillum. 315°/45°, VE 4×\ngeology α = 0.55",
                      (0.99, 0.012), xycoords="axes fraction", ha="right", va="bottom", fontsize=7.5,
                      color="0.25", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.9))

# one label per unit — most-isolated of its largest polygons, kept INSIDE the map frame and
# clear of the two legends and the hillshade note
fig2.canvas.draw()
def _box(artist):
    bb = artist.get_window_extent(); inv = ax2.transData.inverted()
    (x0, y0) = inv.transform((bb.x0, bb.y0)); (x1, y1) = inv.transform((bb.x1, bb.y1))
    return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
keep_out = [_box(leg2), _box(leg_lines), _box(hs_ann)]
def _inbox(x, y, b, pad=80.0):
    return (b[0] - pad) <= x <= (b[2] + pad) and (b[1] - pad) <= y <= (b[3] + pad)
def _cands(c):
    geom = diss.loc[diss["Codigo"] == c, "geometry"].iloc[0]
    polys = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    polys = sorted(polys, key=lambda p: p.area, reverse=True)[:6]
    return [(p.representative_point().x, p.representative_point().y, p.area) for p in polys]
M = 230.0                                            # keep labels this far inside the map frame
xlo, xhi, ylo, yhi = db.left + M, db.right - M, db.bottom + M, db.top - M
def _clamp(x, y): return (min(max(x, xlo), xhi), min(max(y, ylo), yhi))
placed = []
ymid = (db.bottom + db.top) / 2.0
for c in age_order:
    best, best_s = None, -1e18
    for (x, y, ar) in _cands(c):
        xc, yc = _clamp(x, y)
        if any(_inbox(xc, yc, b) for b in keep_out):
            continue
        md = min([((xc - px) ** 2 + (yc - py) ** 2) ** 0.5 for px, py in placed], default=1e9)
        s = min(md, 1500) + 0.005 * ar ** 0.5      # isolation first, larger polygon as tie-break
        if s > best_s:
            best_s, best = s, (xc, yc)
    if best is None:                                # all candidates inside a keep-out -> nudge out
        xc, yc = _clamp(*_cands(c)[0][:2])
        for b in sorted(keep_out, key=lambda bb: bb[1]):
            if _inbox(xc, yc, b):
                yc = b[3] + 260 if (b[1] + b[3]) / 2 < ymid else b[1] - 260
        best = _clamp(xc, yc)
    placed.append(best)
    ax2.text(best[0], best[1], c, ha="center", va="center", fontsize=11, fontweight="bold",
             color="black", path_effects=halo, zorder=9)
ax2.set_xlabel("X (m, EPSG:3763)"); ax2.set_ylabel("Y (m)")
fig2.savefig(os.path.join(OUT, "geology_plan.png"), dpi=160, bbox_inches="tight")
print("saved geology_plan.png")

# ======================= boundary CSV =======================
rows = [{"dist_m": round(db, 1), "elev_m": round(topo_at(db), 2), "from_code": ca, "to_code": cb,
         "x": round(line.interpolate(db).x, 1), "y": round(line.interpolate(db).y, 1),
         "from_unit": desc_map.get(ca, ""), "to_unit": desc_map.get(cb, "")}
        for db, ca, cb in boundaries]
pd.DataFrame(rows).to_csv(os.path.join(OUT, "profile_geology_boundaries.csv"), index=False, encoding="utf-8-sig")
print("DONE")
