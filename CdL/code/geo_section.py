"""Build the geological & hydrostratigraphic cross-sections from the user's digitized shapefile
(geology_crosssection_cdl.shp, in QGIS pixel coords of cross_section_topography.png).

Pipeline: pixel -> section (dist, elev) via section_transform.npz (calibrated in cross_section.py);
clip polygons to below the topography; resolve overlaps (younger wins); fill the rest with PSM;
smooth (Chaikin) the boundaries; make the Qi/PSM contact undulated (erosive). Then:
  - cross_section_geology.png   : units coloured by `codigo`
  - cross_section_hydrostrat.png: merged -> unit 1 alluvium (Qat-NW, K=50), unit 2 terraces
                                  (Qi+Qae+small Qat), unit 3 Plio-Miocene (PSM+PSA)
"""
import config
import os
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.ops import unary_union
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patheffects as pe

OUT  = (str(config.MODEL) + r"\conceptual")
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
SHP  = os.path.join(OUT, "geology_crosssection_cdl.shp")
halo = [pe.withStroke(linewidth=2.2, foreground="white")]

# AML geology codes + official AML 100k symbology RGB (user-provided)
UNIT_COLORS = {"a": "#e1e1e1", "Qdae": "#d1d1c7", "Qt": "#bbbdb8", "Qi'": "#a6aba7",
               "PSA": "#ffffb7", "PSM": "#fff091", "MAT": "#ffff00"}
DESC = {"a": "Aluviões e/ou aterros", "Qdae": "Dunas e areias eólicas", "Qt": "Terraços f./m.",
        "Qi'": "Cascalheiras (terraço)", "PSA": "F. Serra de Almeirim", "PSM": "F. Santa Marta",
        "MAT": "F. Alcoentre–Tomar"}
RANK = {"a": 1, "Qdae": 2, "Qt": 3, "Qi'": 3, "PSA": 4, "PSM": 5}            # young -> old (MAT = basal background)
AGE_ORDER = ["a", "Qdae", "Qt", "Qi'", "PSA", "PSM", "MAT"]                  # recent -> old
BASE_GEO = -35.0

# ---------------- transform + transect ----------------
t = np.load(os.path.join(OUT, "section_transform.npz"))
a, b, c, d = float(t["a"]), float(t["b"]), float(t["c"]), float(t["d"])
def px_to_sec(X, Y):              # X=col, Y=-row (QGIS pixel world); col=a*dist+b, row=c*elev+d
    return (X - b) / a, (-Y - d) / c

pf = np.load(os.path.join(OUT, "profile_data.npz"), allow_pickle=True)
dists, zz, vert_dist = pf["dists"], pf["zz"], pf["vert_dist"]
dmax = float(dists.max())
def topo_at(x): return float(np.interp(x, dists, zz))
ert_dist = {f"P{i}": float(vert_dist[i]) for i in (1, 2, 3, 4)}
zmax = float(np.nanmax(zz)); YTOP = zmax + 8.0

# surface-geology strip (from the map layer)
gmap = gpd.read_file(GPKG, layer="gc_35a_cdl").to_crs(3763)
line = LineString([tuple(v) for v in pf["verts"]])
segs = []
for _, row in gmap.dissolve(by="Codigo", as_index=False).iterrows():
    inter = line.intersection(row.geometry)
    if inter.is_empty:
        continue
    for p in (inter.geoms if inter.geom_type.startswith("Multi") else [inter]):
        if p.geom_type == "LineString" and p.length:
            ds, de = line.project(Point(p.coords[0])), line.project(Point(p.coords[-1]))
            segs.append([min(ds, de), max(ds, de), row["Codigo"]])
segs.sort()

# ---------------- load + transform digitized polygons ----------------
def tr_poly(geom):
    def ring(co): return [px_to_sec(x, y) for x, y in co]
    if geom.geom_type == "Polygon":
        return Polygon(ring(geom.exterior.coords), [ring(r.coords) for r in geom.interiors])
    return MultiPolygon([tr_poly(p) for p in geom.geoms])

shp = gpd.read_file(SHP)
shp["sec"] = shp.geometry.apply(tr_poly)

# ---------------- erosive undulation for the Qi base (NO spline smoothing -> contacts stay EXACT) ----------------
def undulate_qi_base(geom, amp=1.0, taper=320.0, step=10.0):
    """densify + perturb Qi's LOWER boundary by a small irregular sine (erosive), tapered to 0 at both ends
       so the left contact stays progressive and the right (alluvium) contact stays on its vertex."""
    if geom.geom_type != "Polygon":
        return MultiPolygon([undulate_qi_base(p, amp, taper, step) for p in geom.geoms])
    co = list(geom.exterior.coords); dens = []
    for i in range(len(co) - 1):
        p, q = np.array(co[i]), np.array(co[i + 1]); n = max(1, int(np.hypot(*(q - p)) / step))
        for k in range(n):
            dens.append(p + (q - p) * k / n)
    dens.append(np.array(co[-1])); e = np.array(dens); x, y = e[:, 0], e[:, 1]
    ymin, ymax = y.min(), y.max(); xmin, xmax = x.min(), x.max()
    low = np.clip((ymin + 0.55 * (ymax - ymin) - y) / (0.55 * (ymax - ymin) + 1e-9), 0, 1)
    tap = np.clip((x - xmin) / taper, 0, 1) * np.clip((xmax - x) / taper, 0, 1)
    wave = np.sin(2*np.pi*x/250) + 0.55*np.sin(2*np.pi*x/119 + 1.3) + 0.35*np.sin(2*np.pi*x/61 + 2.1)
    out = e.copy(); out[:, 1] = y + amp * low * tap * wave
    return Polygon(out, geom.interiors).buffer(0)

# below-topography clip region
below = Polygon([(dists[0], BASE_GEO)] + list(zip(dists, zz)) + [(dists[-1], BASE_GEO)])

# EXACT digitized contacts (no smoothing, no undulation); MAT is the basal background; clip to below-topo
feats = []
for _, r in shp.iterrows():
    if r["codigo"] == "MAT":
        continue                                     # MAT = basal background (complement of the others)
    g = r["sec"].buffer(0).intersection(below)
    parts = [p for p in (g.geoms if g.geom_type.startswith("Multi") else [g]) if (not p.is_empty and p.area > 4)]
    if parts:
        feats.append([r["codigo"], unary_union(parts), r.geometry.centroid.x])   # drop sliver artifacts

# resolve overlaps (younger wins) + MAT basal background
feats.sort(key=lambda f: RANK.get(f[0], 9))
placed, finals = None, []
for cod, g, cx in feats:
    gg = g if placed is None else g.difference(placed)
    if not gg.is_empty:
        finals.append([cod, gg, cx])
    placed = g if placed is None else placed.union(g)
mat = below.difference(placed)
print("processed units:", [(c, round(g.area)) for c, g, _ in finals], "| MAT (basal) area", round(mat.area))
import pickle
with open(os.path.join(OUT, "_geo.pkl"), "wb") as _f:
    pickle.dump({"finals": [(c, g) for c, g, _ in finals], "psm": mat,
                 "raw": [(r["codigo"], r["sec"]) for _, r in shp.iterrows()]}, _f)

# ---- QA overlay: processed contacts (red) + raw digitized vertices (black) on the reference image ----
ref = os.path.join(OUT, "cross_section_geology - Copy.png")
if os.path.exists(ref):
    from PIL import Image as _IM
    im = np.asarray(_IM.open(ref).convert("RGB")); Hh, Ww = im.shape[:2]
    ext = [(0 - b) / a, (Ww - b) / a, (Hh - d) / c, (0 - d) / c]
    zooms = [(1300, 2900, 18, 42), (2850, 4350, -2, 26), (4300, 5460, -14, 12)]
    figq, axqs = plt.subplots(len(zooms), 1, figsize=(16, 13))
    for axq, (x0, x1, y0, y1) in zip(axqs, zooms):
        axq.imshow(im, extent=ext, origin="upper", aspect="auto", zorder=0)
        for cod, g, _ in finals:
            gpd.GeoSeries([g]).boundary.plot(ax=axq, color="red", lw=1.5, zorder=5)
        for sgeom in shp["sec"]:
            for pp in (sgeom.geoms if sgeom.geom_type.startswith("Multi") else [sgeom]):
                axq.plot(*pp.exterior.xy, "ko", ms=4, zorder=7)
        axq.set_aspect("auto"); axq.set_xlim(x0, x1); axq.set_ylim(y0, y1); axq.grid(alpha=0.3)
    axqs[0].set_title("QA: processed contacts (red) vs raw digitized vertices (black) — should coincide")
    figq.savefig(os.path.join(OUT, "geo_check.png"), dpi=110, bbox_inches="tight"); plt.close(figq)
    print("saved geo_check.png")

# ---------------- shared drawing ----------------
def draw_frame(ax, axs, title):
    # ERT + START/END + topo
    for nm, dd in ert_dist.items():
        ee = topo_at(dd)
        ax.plot([dd, dd], [BASE_GEO, ee], color="crimson", lw=1.0, alpha=0.45, zorder=4)
        ax.plot(dd, ee, "v", color="crimson", ms=11, mec="black", mew=0.7, zorder=9)
        ax.annotate(f"ERT {nm}", (dd, ee), textcoords="offset points", xytext=(0, 13), ha="center",
                    fontsize=11, fontweight="bold", color="crimson", path_effects=halo, zorder=10)
    for dd, lab, ha_, xo in [(0, "START\n(SE divide)", "left", 8), (dmax, "END\n(NW limit)", "right", 2)]:
        ax.plot(dd, topo_at(dd), "o", color="black", ms=7, zorder=9)
        ax.annotate(lab, (dd, topo_at(dd)), textcoords="offset points", xytext=(xo, -14),
                    va="top", ha=ha_, fontsize=10, fontweight="bold", path_effects=halo, zorder=10)
    ax.plot(dists, zz, color="black", lw=1.9, zorder=6)
    ax.set_aspect("auto"); ax.set_ylim(BASE_GEO, YTOP); ax.set_xlim(0, dmax)
    ax.set_ylabel("Elevation (m a.s.l.)"); ax.grid(axis="y", alpha=0.25)
    ax.set_title(title, fontsize=12.5, fontweight="bold")
    ax.annotate("SE", (0.006, 0.985), xycoords="axes fraction", fontsize=13, fontweight="bold", va="top")
    ax.annotate("NW", (0.99, 0.985), xycoords="axes fraction", fontsize=13, fontweight="bold", ha="right", va="top")
    # surface-geology strip
    axs.set_ylim(0, 1)
    for d0, d1, cc in segs:
        axs.fill_betweenx([0, 1], d0, d1, color=UNIT_COLORS.get(cc, "#ddd"))
        if d1 - d0 >= 110:
            axs.text(0.5 * (d0 + d1), 0.5, cc, ha="center", va="center", fontsize=9, fontweight="bold")
    axs.set_yticks([]); axs.set_ylabel("surface\ngeology", fontsize=7.5, rotation=0, ha="right", va="center")
    axs.set_xlabel("Distance along transect (m)   —   SE → NW", labelpad=8)

def plot_units(ax, pairs):
    gdf = gpd.GeoDataFrame({"k": [p[0] for p in pairs]}, geometry=[p[1] for p in pairs])
    for k in dict.fromkeys(gdf["k"]):
        gdf[gdf["k"] == k].plot(ax=ax, color=pairs_color[k], ec="0.3", lw=0.6, zorder=2, aspect=None)

# ============================ FIGURE: geology ============================
pairs_color = UNIT_COLORS
fig = plt.figure(figsize=(16.5, 9.0))
gs = fig.add_gridspec(2, 1, height_ratios=[6, 0.6], hspace=0.10)
ax = fig.add_subplot(gs[0]); axs = fig.add_subplot(gs[1], sharex=ax)
plot_units(ax, [("MAT", mat)] + [[c, g] for c, g, _ in finals])
draw_frame(ax, axs, "CdL / DRYAD — Geological cross-section (SE→NW) — digitized units (codigo), MAT = basal")
present = [c for c in AGE_ORDER if c == "MAT" or any(f[0] == c for f in finals)]
hl = [Patch(fc=UNIT_COLORS[c], ec="0.4", label=f"{c} — {DESC.get(c, '')}") for c in present]
ax.legend(handles=hl, loc="upper right", bbox_to_anchor=(0.992, 0.93), fontsize=9, framealpha=0.93,
          title="Geology (codigo, recent→old)", title_fontsize=9)
fig.savefig(os.path.join(OUT, "cross_section_geology.png"), dpi=170, bbox_inches="tight")
print("saved cross_section_geology.png")

# ============================ FIGURE: hydrostratigraphy (merged) ============================
alluv = [f for f in finals if f[0] == "a"]
nw_a = max(alluv, key=lambda f: f[2]) if alluv else None          # rightmost (NW) 'a' = alluvium (unit 1)
u1 = nw_a[1] if nw_a else Polygon()
u2 = unary_union([g for c, g, cx in finals if c in ("Qdae", "Qt", "Qi'")] +
                 [g for c, g, cx in alluv if nw_a is None or cx != nw_a[2]])   # + the 2 small 'a' polygons
u3 = unary_union([g for c, g, _ in finals if c in ("PSM", "PSA")] + [mat])     # + MAT (basal background)
HCOL = {"u1": "#d9d9d9", "u2": "#6b6b6b", "u3": "#fdf6a8"}   # U1 light gray / U2 dark gray / U3 pale yellow
HK = {"u1": "K = 25 m/d", "u2": "K ≈ 16–77 m/d", "u3": "K ≈ 1–2 m/d"}

fig = plt.figure(figsize=(16.5, 9.0))
gs = fig.add_gridspec(2, 1, height_ratios=[6, 0.6], hspace=0.10)
ax = fig.add_subplot(gs[0]); axs = fig.add_subplot(gs[1], sharex=ax)
pairs_color = {"u3": HCOL["u3"], "u2": HCOL["u2"], "u1": HCOL["u1"]}
plot_units(ax, [("u3", u3), ("u2", u2), ("u1", u1)])
draw_frame(ax, axs, "CdL / DRYAD — Hydrostratigraphic cross-section (SE→NW) — merged units")
def _lab(geom, txt, fs=10):
    if geom.is_empty: return
    p = geom.representative_point()
    ax.text(p.x, p.y, txt, ha="center", va="center", fontsize=fs, fontweight="bold",
            path_effects=halo, zorder=8)
_lab(u3, "Unit 3 — Plio-Miocene\n(MAT + PSA + PSM)\n" + HK["u3"])
_lab(u2, "Unit 2 — terraces\n(Qdae + Qt + Qi')\n" + HK["u2"])
if not u1.is_empty:                                   # alluvium is a small NW wedge -> label with a leader
    p1 = u1.representative_point()
    ax.annotate("Unit 1 — alluvium\n" + HK["u1"], (p1.x, p1.y), xytext=(p1.x - 430, p1.y + 13),
                ha="center", va="center", fontsize=8.5, fontweight="bold", path_effects=halo,
                zorder=10, arrowprops=dict(arrowstyle="-", color="0.3", lw=0.9))
hl = [Patch(fc=HCOL["u1"], ec="0.4", label="Unit 1 — alluvium (a) · K = 25 m/d"),
      Patch(fc=HCOL["u2"], ec="0.4", label="Unit 2 — terraces (Qdae+Qt+Qi') · K ≈ 16–77 m/d"),
      Patch(fc=HCOL["u3"], ec="0.4", label="Unit 3 — Plio-Miocene (MAT+PSA+PSM) · K ≈ 1–2 m/d")]
ax.legend(handles=hl, loc="upper right", bbox_to_anchor=(0.992, 0.93), fontsize=9, framealpha=0.93,
          title="Hydrostratigraphy (recent→old)", title_fontsize=9)
fig.savefig(os.path.join(OUT, "cross_section_hydrostrat.png"), dpi=170, bbox_inches="tight")
print("saved cross_section_hydrostrat.png")
print("DONE")
