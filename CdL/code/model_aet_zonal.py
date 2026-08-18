"""
Model ACTUAL-ET aggregated to land-cover ZONES, monthly (mm) — the model side of
the ET observation group for the PEST++ calibration.

Total model AET per cell = UZET (unsaturated-zone ET, from cdl_gwf.uzf.bud)
                         + UZF-GWET (groundwater ET, from cdl_gwf.cbc).
Both are sinks (q<0); returned as POSITIVE depths.  Zone value = area-weighted
mean depth = (Σ ET volume over the zone) / (Σ cell area) × 1000  [mm/period].

Importable:  compute_zonal_aet(workspace, pest_dir, spinup_nper=12, sim_start=...)
Standalone:  runs on the current WORKSPACE outputs, writes model_aet_zonal.csv + .png
"""
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import flopy
from shapely.geometry import Polygon
from flopy.discretization import VertexGrid

ZONE_NAMES = {1: "oak_broadleaf", 2: "pine", 3: "matos", 4: "grass_crops", 5: "bare_built_water"}


def _uzf_iuzno_to_cell(uzf_file, ncpl):
    """iuzno (1-based, row order) -> cell icpl (0-based) from the UZF packagedata.
    Handles both inline packagedata and the externalised OPEN/CLOSE form (PEST run)."""
    uzf_file = Path(uzf_file)
    rows, blk = [], None
    for ln in uzf_file.read_text().splitlines():
        s = ln.strip(); low = s.lower()
        if low.startswith("begin "): blk = low.split()[1]; continue
        if low.startswith("end "): blk = None; continue
        if blk == "packagedata" and s and not s.startswith("#"):
            if low.startswith("open/close"):
                ext = uzf_file.parent / s.split()[1].strip("'\"")
                rows += [l.strip() for l in ext.read_text().splitlines()
                         if l.strip() and not l.strip().startswith("#")]
            else:
                rows.append(s)
    m = {}
    for s in rows:
        p = s.split()
        m[int(p[0])] = (int(p[2]) - 1) % ncpl          # ifno -> cell icpl (0-based); p[2]=cell 1-based
    return m


def _per_period_by_node(cbf, text, kk):
    """list of dict{node(1-based): q} per kstpkper for an advanced-pkg budget text."""
    out = []
    for k in kk:
        acc = {}
        try:
            data = cbf.get_data(text=text, kstpkper=k)
        except Exception:
            data = None
        if data:
            d = data[0]
            if hasattr(d, "dtype") and d.dtype.names and "node" in d.dtype.names:
                for node, q in zip(d["node"], d["q"]):
                    acc[int(node)] = acc.get(int(node), 0.0) + float(q)
        out.append(acc)
    return out


def compute_zonal_aet(workspace, pest_dir, spinup_nper=12,
                      sim_start=None, model_name="cdl_gwf"):
    workspace, pest_dir = Path(workspace), Path(pest_dir)
    if sim_start is None:            # single source of truth = the model's sidecar (fallback 1981)
        _ssf = Path(r"E:\00code_ws\DRYAD\CdL_model\last_sim_start.txt")
        sim_start = pd.Timestamp(_ssf.read_text().strip()) if _ssf.exists() else pd.Timestamp(1981, 1, 1)
    gp = pickle.load(open(workspace / "voronoi_grid.pkl", "rb"))[0]
    vgrid = VertexGrid(**gp, nlay=1); ncpl = vgrid.ncpl
    area = np.array([Polygon(vgrid.get_cell_vertices(i)).area for i in range(ncpl)])
    zone_id = np.load(pest_dir / "cos_zones.npz")["zone_id"]

    iuz2cell = _uzf_iuzno_to_cell(workspace / f"{model_name}.uzf", ncpl)

    # real (post-spinup) periods
    uzb = flopy.utils.CellBudgetFile(str(workspace / f"{model_name}.uzf.bud"))
    kk = [k for k in uzb.get_kstpkper() if k[1] >= spinup_nper] or uzb.get_kstpkper()
    dates = pd.date_range(sim_start, periods=len(kk), freq="MS")
    days = dates.days_in_month.to_numpy(dtype=float)

    uzet = _per_period_by_node(uzb, "UZET", kk); uzb.close()
    cbc = flopy.utils.CellBudgetFile(str(workspace / f"{model_name}.cbc"), precision="double")
    gwet = _per_period_by_node(cbc, "UZF-GWET", kk); cbc.close()

    # ET VOLUME per cell per period (m3), positive
    vol = np.zeros((len(kk), ncpl))
    for t in range(len(kk)):
        for iuz, q in uzet[t].items():
            c = iuz2cell.get(iuz)
            if c is not None:
                vol[t, c] += -q * days[t]                 # m3/d -> m3 over the period
        for node, q in gwet[t].items():
            c = (node - 1) % ncpl
            vol[t, c] += -q * days[t]

    # zone area-weighted mean depth (mm/period)
    zones = sorted(z for z in np.unique(zone_id) if z != 0)
    out = pd.DataFrame(index=dates)
    for z in zones:
        m = zone_id == z
        A = area[m].sum()
        out[ZONE_NAMES.get(int(z), f"zone{z}")] = (vol[:, m].sum(axis=1) / A) * 1000.0
    out.index.name = "date"
    return out


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    WS = Path(r"E:\00code_ws\DRYAD\CdL_model"); PD = Path(r"E:\00code_ws\DRYAD\CdL_pest")
    df = compute_zonal_aet(WS, PD)
    df.to_csv(PD / "model_aet_zonal.csv")
    print("Zonal AET (mm/month) — long-term means:")
    print(df.mean().round(1).to_string())
    print(f"\n{len(df)} months, annual AET (mm/yr) per zone:")
    print((df.resample("YS").sum().mean().round(0)).to_string())
    fig, ax = plt.subplots(figsize=(11, 4))
    df.plot(ax=ax, lw=0.8); ax.set_ylabel("zone-mean AET (mm/month)")
    ax.set_title("CdL — model actual ET by land-cover zone (monthly)"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(PD / "model_aet_zonal.png", dpi=130, bbox_inches="tight")
    print("wrote model_aet_zonal.csv / .png ->", PD)
