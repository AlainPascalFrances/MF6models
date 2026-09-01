"""diag_divergence.py — FORENSICS for an MF6 convergence/divergence failure.

Parses the listing file for the GW cells and advanced-package features (SFR reach, LAK lake, ...) that fail to
converge, resolves them on the Voronoi grid, then writes:
  - <WS>\\diag\\divergence_map.png    : plan-view map (watershed + SFR network + ponds) with every failing feature
                                        marked, plus a zoom on the dominant one (grid coloured by idomain of its layer)
  - <WS>\\diag\\divergence_report.csv : one row per failing feature, with x/y (EPSG:3763), layer, top/botm/thk,
                                        idomain, K, and flags (lake FID / SFR reach(es) on the cell, reach geometry) —
                                        load the CSV as points in QGIS to inspect alongside the GIS layers.

Run it AFTER a failed run (re-run reproduces the same listing):
    conda run -p C:/sw/miniconda3/envs/mf6models python diag_divergence.py [path\\to\\mfsim.lst]
MF6 listing numbering is 1-BASED (layer, node / reach); flopy arrays are 0-based -> the report shows both.
"""
import config
import re, sys, csv, numpy as np, geopandas as gpd, flopy
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

WS = Path(str(config.MODEL))
GPKG = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
LST = Path(sys.argv[1]) if len(sys.argv) > 1 else WS / "mfsim.lst"
# Save the forensics alongside the run's other outputs in _output\<stamp>\ (stamp from
# last_run_stamp.txt, written by the model). Falls back to the _output root if it's missing.
_stamp = (WS / "last_run_stamp.txt").read_text().strip() if (WS / "last_run_stamp.txt").exists() else ""
OUT = (WS / "_output" / _stamp) if _stamp else (WS / "_output")
OUT.mkdir(parents=True, exist_ok=True)
CRS = 3763

# ---------------- 1. parse the listing for failure signatures ----------------
lines = LST.read_text(errors="ignore").splitlines()
gw, pkg, last_pt, fails = {}, {}, None, []
for ln in lines:
    m = re.search(r"Stress period:\s*(\d+)\s+Time step:\s*(\d+)", ln)
    if m:
        last_pt = (int(m.group(1)), int(m.group(2)))
    for mm in re.finditer(r"-\((\d+),\s*(\d+)\)", ln):                       # GW cell: model-(layer,node)
        kkey = (int(mm.group(1)), int(mm.group(2))); gw[kkey] = gw.get(kkey, 0) + 1
    for mm in re.finditer(r"([A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*)-\((\d+)\)-([a-z]+)", ln):  # <model>-PKG-(id)-term
        pp = mm.group(1).split("-")                                  # strip the model prefix (..GWF-) -> bare pname
        pname = "-".join(pp[1:]) if (len(pp) > 1 and pp[0].endswith("GWF")) else mm.group(1)
        kkey = (pname, int(mm.group(2)), mm.group(3)); pkg[kkey] = pkg.get(kkey, 0) + 1
    if "CONVERGENCE FAILURE" in ln or "failed to converge" in ln.lower():
        fails.append(ln.strip())
if not gw and not pkg:
    print(f"No failure signatures in {LST} — did the run converge?"); sys.exit(0)
top_gw = sorted(gw.items(), key=lambda kv: -kv[1])
top_pkg = sorted(pkg.items(), key=lambda kv: -kv[1])
print(f"Listing            : {LST}")
print(f"Last solved SP/TS  : {last_pt}")
print(f"Failing GW cells   : " + ", ".join(f"(L{l},node{n})x{c}" for (l, n), c in top_gw[:5]))
print(f"Failing pkg feats  : " + ", ".join(f"{p}({i})-{t} x{c}" for (p, i, t), c in top_pkg[:5]))

# ---------------- 2. load the written model ----------------
sim = flopy.mf6.MFSimulation.load(sim_ws=str(WS), verbosity_level=0)
gwf = sim.get_model()
mg = gwf.modelgrid; xc, yc = np.array(mg.xcellcenters), np.array(mg.ycellcenters)
disv = gwf.get_package("DISV"); top = disv.top.array; botm = disv.botm.array; idom = disv.idomain.array
k = gwf.get_package("NPF").k.array; nlay, ncpl = botm.shape

# lake cells -> FID (from LAK connectiondata + boundnames)
lake_cell2fid = {}
try:
    lak = gwf.get_package("LAK"); cd = lak.connectiondata.array; pdd = lak.packagedata.array
    LN = "ifno" if "ifno" in pdd.dtype.names else ("lakeno" if "lakeno" in pdd.dtype.names else pdd.dtype.names[0])
    bn = {int(r[LN]): (str(r["boundname"]) if "boundname" in pdd.dtype.names else f"lake{int(r[LN])}") for r in pdd}
    CLN = "ifno" if "ifno" in cd.dtype.names else cd.dtype.names[0]
    for r in cd:
        lake_cell2fid[int(r["cellid"][-1])] = bn.get(int(r[CLN]), "lake")
except Exception as e:
    print("  (LAK not parsed:", e, ")")

# SFR reach -> cell  and  cell -> reaches
sfr_pd = gwf.get_package("SFR").packagedata.array
RN = "ifno" if "ifno" in sfr_pd.dtype.names else sfr_pd.dtype.names[0]
cell2reach = {}
for r in sfr_pd:
    cell2reach.setdefault(int(r["cellid"][-1]), []).append(int(r[RN]))
def find_reach(rid):
    for cand in (rid, rid - 1, rid + 1):                # listing is 1-based; flopy may store 0- or 1-based
        hit = sfr_pd[sfr_pd[RN] == cand]
        if len(hit):
            return cand, int(hit[0]["cellid"][-1]), hit[0]
    return rid, None, None

# ---------------- 3. build the report rows + collect marker points ----------------
rows, markers = [], []          # markers: (x, y, label, color, marker)
def cellrow(cell0, lay0, kind, listing_id, count, note=""):
    r = dict(kind=kind, listing_id=listing_id, count=count, layer=("" if lay0 is None else lay0 + 1),
             cell0=cell0, node1=cell0 + 1, x=round(float(xc[cell0]), 1), y=round(float(yc[cell0]), 1),
             top=round(float(top[cell0]), 2),
             botm=("" if lay0 is None else round(float(botm[lay0, cell0]), 2)),
             thk=("" if lay0 is None else round(float((top[cell0] if lay0 == 0 else botm[lay0 - 1, cell0]) - botm[lay0, cell0]), 2)),
             idomain=("" if lay0 is None else int(idom[lay0, cell0])),
             K=("" if lay0 is None else round(float(k[lay0, cell0]), 3)),
             lake_fid=lake_cell2fid.get(cell0, ""), sfr_reaches=";".join(map(str, cell2reach.get(cell0, []))), note=note)
    rows.append(r); return r

for (lay1, node1), cnt in top_gw[:6]:                    # GW DVMAX cells (1-based -> 0-based)
    cell0, lay0 = node1 - 1, lay1 - 1
    if not (0 <= cell0 < ncpl and 0 <= lay0 < nlay):
        continue
    r = cellrow(cell0, lay0, "GW_cell(DVMAX)", f"(L{lay1},{node1})", cnt,
                note=("lake " + r0 if (r0 := lake_cell2fid.get(cell0, "")) else "") +
                     (" idomain=-1 PASS-THROUGH" if idom[lay0, cell0] == -1 else "") +
                     (" idomain=0 LAKE-host" if idom[lay0, cell0] == 0 else ""))
    markers.append((r["x"], r["y"], f"GW (L{lay1},{node1})", "red", "*"))

for (p, i, t), cnt in top_pkg[:6]:                       # advanced-package features
    if p == "SFR":
        rid, cell0, rr = find_reach(i)
        if cell0 is not None:
            r = cellrow(cell0, 0, f"{p}-{t}", i, cnt,
                        note=f"reach {rid}: rlen={float(rr['rlen']):.1f} rwid={float(rr['rwid']):.2f} "
                             f"rgrd={float(rr['rgrd']):.1e} rhk={float(rr['rhk']):.3f}")
            markers.append((r["x"], r["y"], f"SFR reach {i}", "darkorange", "D"))
    elif p == "LAK":
        cells = [c for c, f in lake_cell2fid.items()]    # all lake cells; tag the named lake
        fid = next((f for f in lake_cell2fid.values()), "")
        rows.append(dict(kind=f"{p}-{t}", listing_id=i, count=cnt, layer="", cell0="", node1="",
                         x="", y="", top="", botm="", thk="", idomain="", K="", lake_fid=f"lake {i}",
                         sfr_reaches="", note="see lake cells on map"))

with open(OUT / "divergence_report.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT / 'divergence_report.csv'} ({len(rows)} feature row(s))")
for r in rows:
    print(f"  {r['kind']:16s} id={r['listing_id']} L{r['layer']} cell0={r['cell0']} "
          f"x={r['x']} y={r['y']} idom={r['idomain']} K={r['K']} {r['note']}")

# ---------------- 4. map: overview + zoom ----------------
ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(CRS)
try:
    ponds = gpd.read_file(GPKG, layer="ponds_cdl").to_crs(CRS)
except Exception:
    ponds = None
sfr_cells = np.array(sorted(cell2reach.keys()))
lake_cells = np.array(sorted(lake_cell2fid.keys())) if lake_cell2fid else np.array([], int)

fig, (axo, axz) = plt.subplots(1, 2, figsize=(20, 10))
for ax in (axo, axz):
    ws.boundary.plot(ax=ax, color="k", lw=0.8)
    ax.scatter(xc[sfr_cells], yc[sfr_cells], s=3, c="steelblue", label="SFR reaches")
    if lake_cells.size:
        ax.scatter(xc[lake_cells], yc[lake_cells], s=6, c="deepskyblue", label="lake cells")
    if ponds is not None:
        ponds.boundary.plot(ax=ax, color="navy", lw=0.8)
    ax.set_aspect("equal")
for (mx, my, lab, col, mk) in markers:
    for ax in (axo, axz):
        ax.scatter([mx], [my], s=260, c=col, marker=mk, edgecolors="k", linewidths=1.2, zorder=9)
        ax.annotate(lab, (mx, my), color=col, fontsize=9, fontweight="bold",
                    xytext=(6, 6), textcoords="offset points", zorder=10)
axo.set_title(f"DIVERGENCE — overview  (SP/TS {last_pt}; {LST.name})", fontsize=11)
_mk_handles = [
    Line2D([0], [0], marker="*", color="w", markerfacecolor="red", markeredgecolor="k", markersize=16,
           label="GW DVMAX cell (diverging)"),
    Line2D([0], [0], marker="D", color="w", markerfacecolor="darkorange", markeredgecolor="k", markersize=11,
           label="failing advanced-pkg feature (SFR/LAK)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=6, label="SFR reaches"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="deepskyblue", markersize=7, label="lake cells"),
    Line2D([0], [0], color="navy", lw=1.2, label="ponds"),
    Line2D([0], [0], color="k", lw=0.8, label="watershed"),
]
axo.legend(handles=_mk_handles, loc="upper right", fontsize=8, framealpha=0.92)

# zoom on the dominant marker (first), colour the grid by idomain of its layer
if markers:
    cx, cy = markers[0][0], markers[0][1]
    lay_for_zoom = (top_gw[0][0][0] - 1) if top_gw else 0
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=axz, layer=lay_for_zoom)
    bc = ListedColormap(["#9ecae1", "crimson", "0.85"]); bn_ = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], bc.N)
    pmv.plot_array(idom[lay_for_zoom], cmap=bc, norm=bn_, alpha=0.55)
    pmv.plot_grid(lw=0.2, color="0.6")
    W = 900
    axz.set_xlim(cx - W, cx + W); axz.set_ylim(cy - W, cy + W)
    axz.set_title(f"ZOOM on {markers[0][2]}  —  grid colour = idomain (layer {lay_for_zoom+1})", fontsize=10)
    axz.legend(handles=[
        Patch(facecolor="#9ecae1", edgecolor="0.6", label="idomain −1 (pass-through)"),
        Patch(facecolor="crimson", edgecolor="0.6", label="idomain 0 (lake host)"),
        Patch(facecolor="0.85", edgecolor="0.6", label="idomain 1 (active)"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor="red", markeredgecolor="k",
               markersize=14, label="diverging cell"),
    ], loc="upper right", fontsize=8, framealpha=0.92)
fig.savefig(OUT / "divergence_map.png", dpi=140, bbox_inches="tight")
print(f"wrote {OUT / 'divergence_map.png'}")
