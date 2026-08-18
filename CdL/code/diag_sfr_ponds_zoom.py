r"""
Zoom panels of ALL ponds: grid cells, model SFR cells, LAK cells, footprint, stream.

Reads the WRITTEN model (cdl_gwf.nam) and is STALE-PROOF:
  - LAK cells are read from cdl_gwf.lak ONLY if LAK is an active package in cdl_gwf.nam
    (an ss-mode build has no lakes; an orphaned .lak from an older grid must be ignored —
    that stale-file mix produced the 2026-07-04 "LAK cells outside the ponds" artifact).
  - If the current model has no LAK, the PLANNED lake cells (cell centroid inside the
    footprint, as §5b selects them) are shown instead, hatched.
  - Any package node id >= ncpl aborts: the package files belong to a DIFFERENT grid.

Output: diag\sfr_ponds_zoom_p<N>.png (PONDS_PER_PAGE panels per page, all ponds).
"""
import pickle
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, Patch
from matplotlib.collections import PatchCollection
import geopandas as gpd
from shapely.geometry import MultiPolygon, Point
from shapely.ops import unary_union, linemerge
from shapely.validation import make_valid
from pyproj import CRS
from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect

WS      = Path(r"E:\00code_ws\DRYAD\CdL_model")
GPKG    = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\dryad_modelo_NbS.gpkg"
TARGET  = CRS.from_epsg(3763)
MARGIN  = 25.0                  # m around each pond bbox
PONDS_PER_PAGE = 6              # 2 x 3 panels per page

# ---- grid --------------------------------------------------------------------
with open(WS / "voronoi_grid.pkl", "rb") as f:
    gridprops_vg, _ = pickle.load(f)
vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
xc = np.array([vgrid.xcellcenters[i] for i in range(ncpl)])
yc = np.array([vgrid.ycellcenters[i] for i in range(ncpl)])
try:
    ix = GridIntersect(vgrid, method="vertex")
except TypeError:
    ix = GridIntersect(vgrid)

# ---- which packages does the CURRENT model actually have? ---------------------
nam = (WS / "cdl_gwf.nam").read_text().lower()
has_lak = ".lak" in nam
has_sfr = ".sfr" in nam

# ---- generation check: were the written packages built ON THIS grid? ----------
# (node-id bounds alone cannot catch a SMALLER older grid: e.g. 5922-grid ids all fit
# in a 6609 grid but land in the wrong cells.)  The written DISV records its NCPL.
_disv_ncpl = None
for line in (WS / "cdl_gwf.disv").read_text().splitlines():
    if line.strip().upper().startswith("NCPL"):
        _disv_ncpl = int(line.split()[1]); break
stale = _disv_ncpl != ncpl
if stale:
    print(f"!! written model is STALE (cdl_gwf.disv ncpl={_disv_ncpl} vs grid {ncpl}) — "
          f"showing GRID + PLANNED lake cells only; re-run the model to see SFR/LAK.")
    has_lak = has_sfr = False

def _read_pkg_nodes(path, block, col, guard_label):
    """Node ids (0-based) from a package block; abort if any id exceeds this grid."""
    out = []
    blk = False
    for line in path.read_text().splitlines():
        s = line.strip(); low = s.lower()
        if low.startswith(f"begin {block}"): blk = True; continue
        if low.startswith(f"end {block}"): break
        if blk and s and not s.startswith("#"):
            out.append(int(s.split()[col]) - 1)
    bad = [n for n in out if n >= ncpl]
    if bad:
        raise SystemExit(
            f"!! STALE FILE: {path.name} {guard_label} has node ids up to {max(bad)+1} but this grid has "
            f"ncpl={ncpl}. The package file was written for a DIFFERENT grid — re-run the model first.")
    return out

model_sfr = set(_read_pkg_nodes(WS / "cdl_gwf.sfr", "packagedata", 2, "reach")) if has_sfr else set()
lak_cells = {}
if has_lak:
    blk = False
    for line in (WS / "cdl_gwf.lak").read_text().splitlines():
        s = line.strip(); low = s.lower()
        if low.startswith("begin connectiondata"): blk = True; continue
        if low.startswith("end connectiondata"): break
        if blk and s and not s.startswith("#"):
            p = s.split()
            lak_cells[int(p[3]) - 1] = int(p[0]) - 1
    if lak_cells and max(lak_cells) >= ncpl:
        raise SystemExit(f"!! STALE FILE: cdl_gwf.lak node ids exceed ncpl={ncpl} — different grid; re-run the model.")
print(f"current model: SFR={'yes' if has_sfr else 'NO'} ({len(model_sfr)} cells), "
      f"LAK={'yes' if has_lak else 'NO (ss build — showing PLANNED lake cells)'} ({len(lak_cells)} cells)")

# ---- GIS ----------------------------------------------------------------------
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
ws_poly = unary_union(ws.geometry.apply(lambda g: make_valid(g)).buffer(0).values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
streams = gpd.read_file(GPKG, layer="streams_cdl").to_crs(TARGET)
merged = linemerge(unary_union(streams.geometry.apply(make_valid).values))
stream_union = (merged if merged.geom_type == "LineString"
                else unary_union(list(merged.geoms))).intersection(ws_poly)
ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
ponds["geometry"] = ponds.geometry.apply(make_valid)

def _largest(g):
    return max(g.geoms, key=lambda p: p.area) if isinstance(g, MultiPolygon) else g

pond_list = []                                    # (fid, clipped footprint)
for fid in range(len(ponds)):
    g = _largest(ponds.geometry.iloc[fid]).intersection(ws_poly)
    if not g.is_empty:
        pond_list.append((fid, _largest(g)))

# ---- pages ----------------------------------------------------------------------
outdir = WS / "diag"
outdir.mkdir(exist_ok=True)
npages = int(np.ceil(len(pond_list) / PONDS_PER_PAGE))
for page in range(npages):
    chunk = pond_list[page * PONDS_PER_PAGE:(page + 1) * PONDS_PER_PAGE]
    fig, axes = plt.subplots(2, 3, figsize=(16, 11))
    for ax in axes.ravel()[len(chunk):]:
        ax.axis("off")
    for ax, (fid, g) in zip(axes.ravel(), chunk):
        x0, y0, x1, y1 = g.bounds
        x0 -= MARGIN; y0 -= MARGIN; x1 += MARGIN; y1 += MARGIN
        win = [i for i in range(ncpl) if x0 <= xc[i] <= x1 and y0 <= yc[i] <= y1]
        ax.add_collection(PatchCollection(
            [MplPoly(vgrid.get_cell_vertices(i)) for i in win],
            facecolor="none", edgecolor="0.75", lw=0.4))
        sfr_win = [i for i in win if i in model_sfr]
        ax.add_collection(PatchCollection(
            [MplPoly(vgrid.get_cell_vertices(i)) for i in sfr_win],
            facecolor="deepskyblue", alpha=0.55, edgecolor="none"))
        # LAK cells: actual (orange) if the model has LAK; else PLANNED (orange hatch) =
        # the §5b selection rule (cell centroid inside the footprint, nearest-cell fallback)
        cand = {int(c) for c in ix.intersect(g, geo_dataframe=False)["cellids"]}
        if has_lak:
            lak_win = [i for i in win if i in lak_cells]
        else:
            lak_win = sorted(nd for nd in cand if g.contains(Point(float(xc[nd]), float(yc[nd]))))
            if not lak_win and cand:
                c = g.representative_point()
                lak_win = [int(min(cand, key=lambda nd: (xc[nd] - c.x) ** 2 + (yc[nd] - c.y) ** 2))]
        if lak_win:
            ax.add_collection(PatchCollection(
                [MplPoly(vgrid.get_cell_vertices(i)) for i in lak_win],
                facecolor="orange", alpha=0.75, edgecolor="crimson", lw=1.2,
                hatch=None if has_lak else "//"))
        xs, ys = g.exterior.xy
        ax.plot(xs, ys, color="navy", lw=1.8)
        sx = stream_union.intersection(g.buffer(MARGIN))
        for ln in ([sx] if sx.geom_type == "LineString" else getattr(sx, "geoms", [])):
            if ln.geom_type == "LineString":
                lx, ly = ln.xy; ax.plot(lx, ly, color="tab:blue", lw=1.2)
        nsfr = len(cand & model_sfr)
        ax.set_title(f"pond {fid}: {len(lak_win)} LAK cell(s), {nsfr} SFR cell(s) in footprint", fontsize=10)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
    fig.legend(handles=[
        Patch(facecolor="deepskyblue", alpha=0.55, label="model SFR cells (cdl_gwf.sfr)"),
        Patch(facecolor="orange", alpha=0.75, edgecolor="crimson",
              hatch=None if has_lak else "//",
              label="LAK cells (cdl_gwf.lak)" if has_lak else "PLANNED LAK cells (no LAK in current build)"),
        Patch(facecolor="none", edgecolor="navy", label="pond footprint"),
    ], loc="upper center", ncol=3, fontsize=10)
    fig.suptitle(f"Ponds — LAK vs SFR cells (page {page + 1}/{npages})", y=0.955)
    if stale:
        fig.text(0.5, 0.935,
                 f"⚠ written model files are from a DIFFERENT grid (ncpl {_disv_ncpl} ≠ {ncpl}) — "
                 "SFR/LAK overlays HIDDEN; re-run the model (ss then transient) and regenerate this figure",
                 ha="center", color="crimson", fontsize=11, fontweight="bold")
    elif not has_lak:
        fig.text(0.5, 0.935,
                 "current build is the SS run (no LAK package) — orange = PLANNED lake cells; "
                 "run the transient to see the actual LAK + excised SFR",
                 ha="center", color="darkorange", fontsize=11, fontweight="bold")
    out = outdir / f"sfr_ponds_zoom_p{page + 1}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
