r"""
diag_mvr_arrows.py — visualise every MVR mover as a provider->receiver arrow.

Reads the WRITTEN MODFLOW 6 packages (cdl_gwf.mvr / .uzf / .sfr / .lak / .drn) and
maps each mover's package-internal feature id back to its grid cell, so the moves
can be checked on the grid.  The MVR file references features by their PACKAGE
index (UZF -> iuzno, SFR -> reach, LAK -> lake, DRN -> boundary index), NOT by
cell/node number — this tool resolves that mapping and tests face-adjacency.

Outputs (WORKSPACE/diag/):
  mvr_arrows.shp   provider-centroid -> receiver-centroid LineStrings (EPSG:3763)
                   attrs: pkg_from,id_from,node_from, pkg_to,id_to,node_to,
                          factor, length_m, adjacent, kind
  mvr_arrows.csv   the same table (no geometry)
  mvr_arrows.png   grid + arrows (CRR green / SFR-LAK blue; any NON-adjacent red)

Load mvr_arrows.shp over the grid in QGIS: every CRR arrow should join two cells
that share a face.  A red arrow (non-adjacent) would flag a real routing bug.
"""
import pickle
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import geopandas as gpd
from shapely.geometry import LineString, Point, MultiPolygon
from shapely.ops import unary_union
from shapely.validation import make_valid
from pyproj import CRS
import flopy
from flopy.discretization import VertexGrid
from flopy.utils.gridintersect import GridIntersect

WS   = Path(r"E:\00code_ws\DRYAD\CdL_model")
MODEL = "cdl_gwf"
GPKG = r"E:/zzCloud/OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia/DRYAD/GIS/dryad_modelo_NbS.gpkg"
TARGET = CRS.from_epsg(3763)
OUT  = WS / "diag"; OUT.mkdir(exist_ok=True)

# ---- grid: cell centroids + face-neighbour sets (shared Voronoi edge) ---------
with open(WS / "voronoi_grid.pkl", "rb") as f:
    gp_vg, _ = pickle.load(f)
cell2d = gp_vg["cell2d"]; ncpl = gp_vg["ncpl"]
xc = np.array([c[1] for c in cell2d]); yc = np.array([c[2] for c in cell2d])
_e2c = defaultdict(list)
for c2 in cell2d:
    icc = int(c2[0]); nv = int(c2[3]); ivl = [int(v) for v in c2[4:4 + nv]]
    for a, b in zip(ivl, ivl[1:] + ivl[:1]):
        _e2c[(a, b) if a < b else (b, a)].append(icc)
nbrs = defaultdict(set)
for cl in _e2c.values():
    if len(cl) == 2:
        nbrs[cl[0]].add(cl[1]); nbrs[cl[1]].add(cl[0])

# ---- id -> node maps from the written packages (file ids are 1-based) ---------
sim = flopy.mf6.MFSimulation.load(sim_ws=str(WS), verbosity_level=0,
                                  load_only=["disv", "uzf", "sfr", "lak", "drn"])
gwf = sim.get_model(MODEL)

def _node_map_from_pkgdata(pkg, idfield="ifno"):
    return {int(r[idfield]) + 1: int(r["cellid"][1]) for r in pkg.packagedata.array}

maps = {}          # single-cell packages: file id (1-based) -> node (0-based)
info = {}          # file id (1-based) -> (lay1, node1, landflag)  [landflag=-1 where N/A]
uzf = gwf.get_package("UZF")
if uzf is not None:
    maps["UZF"] = _node_map_from_pkgdata(uzf)
    info["UZF"] = {int(r["ifno"]) + 1: (int(r["cellid"][0]) + 1, int(r["cellid"][1]) + 1, int(r["landflag"]))
                   for r in uzf.packagedata.array}
sfr = gwf.get_package("SFR")
if sfr is not None:
    maps["SFR"] = _node_map_from_pkgdata(sfr)
    info["SFR"] = {int(r["ifno"]) + 1: (int(r["cellid"][0]) + 1, int(r["cellid"][1]) + 1, -1)
                   for r in sfr.packagedata.array}
drnseep = gwf.get_package("DRN-SEEP")
if drnseep is not None:
    rec = drnseep.stress_period_data.get_data(0)     # order == boundary index
    maps["DRN-SEEP"] = {i + 1: int(r["cellid"][1]) for i, r in enumerate(rec)}
    info["DRN-SEEP"] = {i + 1: (int(r["cellid"][0]) + 1, int(r["cellid"][1]) + 1, -1) for i, r in enumerate(rec)}

# CROSS-REFERENCE CSV: every UZF object id (1-based) <-> its (layer, node, landflag, ivertcon), all
# 1-based — so any 'UZF <id>' in the .mvr can be looked up unambiguously (id != node != layer!).
if uzf is not None:
    import csv as _csv
    with open(OUT / "uzf_id_cell_crossref.csv", "w", newline="") as _f:
        _w = _csv.writer(_f)
        _w.writerow(["uzf_id_1based", "layer_1based", "node_1based", "landflag", "ivertcon_1based"])
        for r in sorted(uzf.packagedata.array, key=lambda r: int(r["ifno"])):
            _iv = int(r["ivertcon"]); _w.writerow([int(r["ifno"]) + 1, int(r["cellid"][0]) + 1,
                     int(r["cellid"][1]) + 1, int(r["landflag"]), (_iv + 1 if _iv >= 0 else 0)])
    print(f"   wrote uzf_id_cell_crossref.csv ({len(uzf.packagedata.array)} UZF objects)")

# ---- LAK: a lake is a MULTI-CELL footprint (only its 1 embedded connection cell
#      is in the written file). Reconstruct the full footprint per lake so a mover
#      to/from a lake is tested against ALL its cells, not just the connection one.
lake_cells = {}    # LAK file id (1-based) -> set of footprint nodes
lak_lay = {}       # LAK file id (1-based) -> its embedded-connection layer (1-based)
lak = gwf.get_package("LAK")
if lak is not None:
    for r in lak.connectiondata.array:
        lak_lay.setdefault(int(r["ifno"]) + 1, int(r["cellid"][0]) + 1)
    lake_fid = {}
    for r in lak.packagedata.array:
        bn = str(r["boundname"]).strip()
        lake_fid[int(r["ifno"]) + 1] = int("".join(ch for ch in bn if ch.isdigit()))
    ws_g = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(TARGET)
    ws_poly = unary_union(ws_g.geometry.apply(lambda g: make_valid(g)).buffer(0).values)
    if isinstance(ws_poly, MultiPolygon):
        ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
    ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(TARGET)
    ponds["geometry"] = ponds.geometry.apply(make_valid)
    vgrid = VertexGrid(**gp_vg, nlay=1)
    try:
        ix = GridIntersect(vgrid, method="vertex")
    except TypeError:
        ix = GridIntersect(vgrid)
    def _footprint(pg):                              # = the model's pond_cells_of (centroid-inside)
        try:
            res = ix.intersect(pg, geo_dataframe=False)
        except TypeError:
            res = ix.intersect(pg)
        cand = sorted({int(c) for c in res["cellids"]})
        nds = [nd for nd in cand if pg.contains(Point(float(xc[nd]), float(yc[nd])))]
        if not nds and cand:
            c = pg.representative_point()
            nds = [min(cand, key=lambda nd: (xc[nd] - c.x) ** 2 + (yc[nd] - c.y) ** 2)]
        return set(nds)
    for lid, fid in lake_fid.items():
        g = ponds.geometry.iloc[fid].intersection(ws_poly)
        if isinstance(g, MultiPolygon):
            g = max(g.geoms, key=lambda p: p.area)
        lake_cells[lid] = _footprint(g) or {maps.get("LAK-conn", {}).get(lid, -1)}

print("resolved id->cell maps:", {k: len(v) for k, v in maps.items()},
      "| LAK footprints:", {k: len(v) for k, v in lake_cells.items()})

def resolve(pkg, i, other_node):
    """Return (arrow_node, candidate_cells) for a mover endpoint. For LAK, pick the
    footprint cell nearest the OTHER endpoint (= the cell the CRR/coupling logic used)."""
    if pkg == "LAK":
        cand = lake_cells.get(i, set())
        if not cand:
            return None, set()
        node = (min(cand, key=lambda c: (xc[c] - xc[other_node]) ** 2 + (yc[c] - yc[other_node]) ** 2)
                if other_node is not None else next(iter(cand)))
        return node, cand
    node = maps.get(pkg, {}).get(i)
    return node, ({node} if node is not None else set())

def endpoint_info(pkg, fid, arrow_node0):
    """(layer_1based, node_1based, landflag) for a mover endpoint. UZF/SFR/DRN read their own
    object's cell; LAK reports the embedded-connection layer + the resolved footprint node."""
    if pkg in info and fid in info[pkg]:
        return info[pkg][fid]
    if pkg == "LAK":
        return (lak_lay.get(fid, 0), (arrow_node0 + 1) if arrow_node0 is not None else -1, -1)
    return (-1, -1, -1)

# ---- parse the MVR period-1 movers, resolving the non-LAK side first ----------
# columns (all cell ids 1-based, matching the .uzf/.mvr files):
#   pkg_from, uzfid_from (the file id!), lay_from, node_from, lflag_from,
#   pkg_to,   uzfid_to,                  lay_to,   node_to,   lflag_to,
#   factor, length_m, adjacent, kind      (node_from/to are 1-based cell nodes; lflag -1 = N/A)
cols = ["pkg_from", "uzfid_from", "lay_from", "node_from", "lflag_from",
        "pkg_to", "uzfid_to", "lay_to", "node_to", "lflag_to",
        "factor", "length_m", "adjacent", "kind"]
rows = []; geom_nodes = []
blk = False
for ln in (WS / f"{MODEL}.mvr").read_text().splitlines():
    s = ln.strip(); low = s.lower()
    if low.startswith("begin period"):
        blk = True; continue
    if low.startswith("end period"):
        break
    if not blk or not s or s.startswith("#"):
        continue
    p = s.split()
    if len(p) < 6:
        continue
    p1, i1, p2, i2, mtype, val = p[0], int(p[1]), p[2], int(p[3]), p[4], float(p[5])
    kind = "CRR" if p1 in ("UZF", "DRN-SEEP") else "SFR-LAK"
    # resolve the non-LAK side first so the LAK side can pick its nearest footprint cell
    if p1 == "LAK":
        n2, c2set = resolve(p2, i2, None); n1, c1set = resolve(p1, i1, n2)
    else:
        n1, c1set = resolve(p1, i1, None); n2, c2set = resolve(p2, i2, n1)
    la1, nd1, lf1 = endpoint_info(p1, i1, n1)
    la2, nd2, lf2 = endpoint_info(p2, i2, n2)
    if n1 is None or n2 is None:
        rows.append((p1, i1, la1, nd1, lf1, p2, i2, la2, nd2, lf2, val, np.nan, "UNRESOLVED", kind))
        geom_nodes.append((None, None)); continue
    adj = "yes" if (c1set & c2set) or any(b in nbrs[a] for a in c1set for b in c2set) else "NO"
    rows.append((p1, i1, la1, nd1, lf1, p2, i2, la2, nd2, lf2, val,
                 float(np.hypot(xc[n1] - xc[n2], yc[n1] - yc[n2])), adj, kind))
    geom_nodes.append((n1, n2))

# ---- geometry + write ---------------------------------------------------------
geoms, attrs = [], []
for r, (n1, n2) in zip(rows, geom_nodes):
    if n1 is None or n2 is None or n1 < 0 or n2 < 0:
        continue
    if n1 == n2:                                       # self (e.g. reach in provider cell) -> tiny stub
        geoms.append(LineString([(xc[n1], yc[n1]), (xc[n1] + 0.5, yc[n1] + 0.5)]))
    else:
        geoms.append(LineString([(xc[n1], yc[n1]), (xc[n2], yc[n2])]))
    attrs.append(r)
gdf = gpd.GeoDataFrame(attrs, columns=cols, geometry=geoms, crs=TARGET)
gdf.to_file(OUT / "mvr_arrows.shp")
gdf.drop(columns="geometry").to_csv(OUT / "mvr_arrows.csv", index=False)

# ---- summary ------------------------------------------------------------------
import pandas as pd
df = pd.DataFrame(rows, columns=cols)
n_crr = int((df.kind == "CRR").sum()); n_cpl = int((df.kind == "SFR-LAK").sum())
n_bad_crr = int(((df.kind == "CRR") & (df.adjacent == "NO")).sum())
n_unres = int((df.adjacent == "UNRESOLVED").sum())
print(f"\n{len(df)} movers: {n_crr} CRR (slope cascade), {n_cpl} SFR<->LAK coupling, {n_unres} unresolved")
print(f"CRR moves to a NON-face-adjacent cell: {n_bad_crr}   <-- MUST be 0 (CRR routes to face-neighbours only)")
_g = df[df.adjacent.isin(["yes", "NO"])]
if len(_g):
    print(f"arrow length (m): min {_g.length_m.min():.1f}  mean {_g.length_m.mean():.1f}  "
          f"max {_g.length_m.max():.1f}")
    print("by provider->receiver:")
    print(df.groupby(["pkg_from", "pkg_to", "adjacent"]).size().to_string())

# ---- quick-look PNG -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(13, 12))
edges = []
for c2 in cell2d:
    nv = int(c2[3]); vids = [int(v) for v in c2[4:4 + nv]]
    vxy = {v[0]: (v[1], v[2]) for v in gp_vg["vertices"]}
    pts = [vxy[v] for v in vids]
    edges.append(pts + [pts[0]])
ax.add_collection(LineCollection([[(p[0], p[1]) for p in e] for e in edges],
                                 colors="0.85", linewidths=0.2, zorder=1))
_ADJ, _KIND = cols.index("adjacent"), cols.index("kind")   # tuple positions (schema-proof)
def _seg(pred):
    out = []
    for r, (n1, n2) in zip(rows, geom_nodes):
        if n1 is None or n2 is None or n1 < 0 or n2 < 0 or not pred(r):
            continue
        out.append([(xc[n1], yc[n1]), (xc[n2], yc[n2])])
    return out
crr = _seg(lambda r: r[_KIND] == "CRR" and r[_ADJ] == "yes")
cpl = _seg(lambda r: r[_KIND] == "SFR-LAK" and r[_ADJ] == "yes")
bad = _seg(lambda r: r[_ADJ] == "NO")
ax.add_collection(LineCollection(crr, colors="seagreen", linewidths=0.3, alpha=0.5, zorder=2))
ax.add_collection(LineCollection(cpl, colors="royalblue", linewidths=0.8, alpha=0.8, zorder=3))
if bad:
    ax.add_collection(LineCollection(bad, colors="red", linewidths=1.5, zorder=4))
ax.autoscale(); ax.set_aspect("equal")
ax.set_title(f"CdL MVR movers ({len(df)}): CRR green, SFR<->LAK blue, non-adjacent RED "
             f"({n_bad_crr} bad CRR)")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0], [0], color="seagreen", label=f"CRR cascade ({n_crr})"),
                   Line2D([0], [0], color="royalblue", label=f"SFR<->LAK ({n_cpl})"),
                   Line2D([0], [0], color="red", label=f"non-adjacent ({n_bad_crr+ (n_cpl and 0)})")],
          loc="upper right", fontsize=9)
fig.savefig(OUT / "mvr_arrows.png", dpi=150, bbox_inches="tight")
print(f"\nwrote {OUT / 'mvr_arrows.shp'} (+ .csv, .png)")
