"""
================================================================================
Post-processing for the CdL groundwater-flow model (cdl_gwf_model_opusv1.py)
================================================================================
Reads the MODFLOW 6 output already on disk (no need to re-run the model) and
produces, in WORKSPACE/postproc/:

  1. Head TIME SERIES at the obs_points (cdl_gwf.obs.head.csv)   -> figure + tidy CSV
  2. Mean-head MAPS per layer                                    -> figure + CSV
  2b. Mean-DEPTH MAPS per layer (land surface − head)            -> figure + CSV
  3. WATER BUDGET by compartment (surface / unsaturated / aquifer):
       - overall mean rate (m3/d) from the listing budget        -> figure + CSV
       - yearly volumes (m3/yr)                                   -> figure + CSV
       - unsaturated-zone internal budget (UZF .bud)             -> figure + CSV
       - surface-water internal budget (SFR .bud)                -> figure + CSV
       - per-layer storage change (GWF .cbc)                     -> figure + CSV

The spin-up (steady-state period 0 + the 11 transient months that replay year 1)
is EXCLUDED from every average/map so the statistics reflect the real period.

Run in the activated env, e.g.:
  & "C:\\miniconda3\\Scripts\\conda.exe" run -p C:\\miniconda3\\envs\\flopy ^
      --no-capture-output python -u postprocess_cdl.py
================================================================================
"""

import re
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects

import flopy

# --- CONFIG (keep in sync with cdl_gwf_model_opusv1.py) ----------------------
WORKSPACE   = Path(r"E:\00code_ws\DRYAD\CdL_model")
MODEL_NAME  = "cdl_gwf"
# SIM_START MUST match the model run. The model writes its start date to
# last_sim_start.txt; read it so postproc auto-syncs to whatever window the model
# used (1981-2026 full, or 2010-2026 calibration). Fallback = 1981 (full period).
_ssf = WORKSPACE / "last_sim_start.txt"
SIM_START = pd.Timestamp(_ssf.read_text().strip()) if _ssf.exists() else pd.Timestamp(1981, 1, 1)
print(f">> SIM_START = {SIM_START.date()} ({'from last_sim_start.txt' if _ssf.exists() else 'fallback 1981'})")
SPINUP_NPER = 12          # SS period 0 + 11 transient spin-up months -> first 12 SPs
GPKG        = r"E:/zzCloud/OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia/DRYAD/GIS/dryad_modelo_NbS.gpkg"
OBS_LAYER   = "obs_points_cdl"
P_CSV       = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\p_month_198101_202605.csv"
ET_CSV      = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\et0_month_198101_202605.csv"

PIEZO_XLSX  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\piezos_qualitative_month_198101_202605.xlsx"
# SYNTHETIC piezometer heads (regionalised from SNIRH analogs) — the real observation
# overlay for obs_head_timeseries.png, REPLACING the qualitative placeholders above.
# Long format: cdl (P0..P6), snirh, date, wl_elev_m, depth_m. Falls back to the xlsx if absent.
SYN_PIEZO_CSV = r"E:\00code_ws\DRYAD\CdL_pest\snirh_data_availability\cdl_synthetic_piezo.csv"

# Match the model run's start stamp (written by the main script) so preproc\<stamp>\
# and postproc\<stamp>\ correspond to the same run; fall back to now if absent.
_stampf = WORKSPACE / "last_run_stamp.txt"
RUN_STAMP = _stampf.read_text().strip() if _stampf.exists() else datetime.now().strftime("%Y%m%d%H%M")
OUT = WORKSPACE / "_output" / RUN_STAMP
OUT.mkdir(exist_ok=True, parents=True)

HDS  = WORKSPACE / f"{MODEL_NAME}.hds"
CBC  = WORKSPACE / f"{MODEL_NAME}.cbc"
LST  = WORKSPACE / f"{MODEL_NAME}.lst"
UZFB = WORKSPACE / f"{MODEL_NAME}.uzf.bud"
SFRB = WORKSPACE / f"{MODEL_NAME}.sfr.bud"
OBSCSV = WORKSPACE / f"{MODEL_NAME}.obs.head.csv"
OBS_SFR = WORKSPACE / f"{MODEL_NAME}.obs.sfr.csv"        # SFR inlet/outlet ext-flow obs (per timestep)
SFR_INLET_SERIES = WORKSPACE / "sfr_inlet_series.csv"   # regionalised synthetic inlet/outlet targets
MFSIM_LST = WORKSPACE / "mfsim.lst"                     # simulation listing (per-SP iterations + elapsed time)
GRB  = WORKSPACE / f"{MODEL_NAME}.disv.grb"
ROUTED_SHP = WORKSPACE / "streams_cdl_routed.shp"


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
def load_grid():
    """Return (modelgrid, top[ncpl], nlay) from the binary grid file (fallback: pkl)."""
    try:
        grb = flopy.mf6.utils.MfGrdFile(str(GRB))
        mg = grb.modelgrid
        return mg, np.asarray(mg.top, dtype=float).ravel(), mg.nlay
    except Exception as e:
        print(f"   (grb load failed: {e!r}; falling back to voronoi_grid.pkl, no top)")
        import pickle
        with open(WORKSPACE / "voronoi_grid.pkl", "rb") as f:
            gridprops_vg, _ = pickle.load(f)
        # nlay inferred later from heads; build a 1-layer grid just for geometry
        from flopy.discretization import VertexGrid
        mg = VertexGrid(**gridprops_vg, nlay=1)
        return mg, None, None


def real_kstpkper(binfile):
    """(kstp,kper) saved AFTER the spin-up (kper >= SPINUP_NPER)."""
    return [kk for kk in binfile.get_kstpkper() if kk[1] >= SPINUP_NPER]


def real_dates(n):
    """Monthly calendar dates for the n real stress periods (post spin-up)."""
    return pd.date_range(SIM_START, periods=n, freq="MS")


def read_forcing(csv, value_col="MONTHLY_TOTAL"):
    """Monthly totals (mm) as a Series indexed by month-start date (same parse as the model)."""
    df = pd.read_csv(csv, sep=r"\s+|,|;", engine="python")
    datec = [c for c in df.columns if "date" in c.lower()][0]
    if value_col not in df.columns:                       # fall back to last column
        value_col = df.columns[-1]
    df[datec] = pd.to_datetime(df[datec], format="%m/%d/%Y", errors="coerce")
    s = (df[[datec, value_col]].rename(columns={datec: "date", value_col: "val"})
         .dropna().set_index("date").sort_index())
    return s["val"]


def compartment_of(term):
    t = term.upper()
    if "STO" in t:
        return "aquifer (storage)"
    if "UZF" in t:
        return "unsaturated (UZF->GW)"
    if "SFR" in t:
        return "surface (stream)"
    if "DRN" in t:
        return "surface (drain)"
    return "other"


COMP_COLORS = {
    "aquifer (storage)": "tab:brown",
    "unsaturated (UZF->GW)": "tab:green",
    "surface (stream)": "tab:blue",
    "surface (drain)": "tab:cyan",
    "other": "0.6",
}

# Self-explanatory label + colour per RAW listing-budget term.  The three DRN
# packages (name-file order: DRN-OUT, DRN-SEEP, DRN-WEST) are enumerated by MF6 as
# 2026-08-05: DRN-OUT retired (secondary outlet -> explicit DRN-SECONDARY). DRN packages
# now build in order DRN-SEEP, DRN-WEST, DRN-SECONDARY; only DRN-SEEP is a mover, so MF6
# lists the family (verified from cdl_gwf.lst) as:
#   DRN          = DRN-SEEP -> surface seepage leaving directly (topographic sinks)
#   DRN-TO-MVR   = DRN-SEEP -> MVR -> SFR (the seepage baseflow to the streams; the big term)
#   DRN2         = DRN-WEST -> western alluvial underflow out (U1+U3)
#   DRN3         = DRN-SECONDARY -> secondary catchment outlet
# DRN family = BLUE TONES; UZF-GWET = orange.
TERM_STYLE = {
    "STO-SS":      ("storage — specific storage",       "tab:brown"),
    "STO-SY":      ("storage — specific yield",         "#c49a6c"),
    "UZF-GWRCH":   ("UZF recharge (in)",                "tab:green"),
    "UZF-GWET":    ("UZF evapotranspiration (out)",     "tab:orange"),
    "SFR":         ("SFR — stream↔aquifer",             "tab:blue"),
    "LAK":         ("LAK — pond↔aquifer",               "navy"),
    "GHB":         ("GHB east inflow (U1+U3)",          "teal"),
    "DRN":         ("DRN-seep surface seepage (direct)", "#6baed6"),
    "DRN-TO-MVR":  ("DRN-seep→SFR baseflow",            "#4292c6"),
    "DRN2":        ("DRN west — alluvial underflow out", "#08519c"),
    "DRN3":        ("DRN secondary outlet",             "#9ecae1"),
}


def term_label(term):
    return TERM_STYLE.get(term, (term, None))[0]


def term_color(term):
    c = TERM_STYLE.get(term, (None, None))[1]
    return c if c is not None else COMP_COLORS.get(compartment_of(term), "0.6")

# Brown palette for the head curves — ties to the aquifer/storage color (tab:brown)
# and keeps the curves distinct from the blue rainfall / orange ET0 bars.
HEAD_COLORS = ["#5c3211", "#8c564b", "#b3743a", "#cf9d62"]  # dark brown -> tan (#8c564b = tab:brown)


def read_piezo_obs():
    """Observed heads per piezometer: the SYNTHETIC water-table elevations regionalised
    from SNIRH analogs (cdl_synthetic_piezo.csv, wl_elev_m), replacing the qualitative
    placeholders. Returns {POINTNAME: DataFrame[date, head]} keyed to the obs points
    (P0..P6). Falls back to the qualitative Excel if the synthetic csv is unavailable."""
    out = {}
    # --- primary: synthetic SNIRH-regionalised heads (long format) ---
    if Path(SYN_PIEZO_CSV).exists():
        try:
            syn = pd.read_csv(SYN_PIEZO_CSV, parse_dates=["date"])
            for pt, g in syn.groupby("cdl"):
                s = (g[["date", "wl_elev_m"]].rename(columns={"wl_elev_m": "head"})
                     .dropna().sort_values("date"))
                out[str(pt).upper()] = s
            print(f"   piezo overlay: SYNTHETIC heads for {sorted(out)} "
                  f"({sum(len(v) for v in out.values())} obs)")
            return out
        except Exception as e:
            print(f"   (synthetic piezo csv not read: {e!r}; falling back to xlsx)")
    # --- fallback: qualitative Excel (one sheet per point, p0..p6) ---
    try:
        xl = pd.ExcelFile(PIEZO_XLSX)
    except Exception as e:
        print(f"   (piezo Excel not read: {e!r})")
        return out
    for sh in xl.sheet_names:
        m = re.match(r"\s*(p\d+)", sh, re.IGNORECASE)
        name = m.group(1).upper() if m else sh
        d = pd.read_excel(PIEZO_XLSX, sheet_name=sh)
        datec = [c for c in d.columns if "date" in c.lower()]
        headc = [c for c in d.columns if "head" in c.lower()]
        if not datec or not headc:
            continue
        s = (d[[datec[0], headc[0]]].rename(columns={datec[0]: "date", headc[0]: "head"}))
        s["date"] = pd.to_datetime(s["date"], errors="coerce")
        out[name] = s.dropna().sort_values("date")
    return out


# =============================================================================
# 1. HEAD TIME SERIES AT OBS POINTS
# =============================================================================
def obs_timeseries():
    print(">> [1] Obs-point head time series …")
    if not OBSCSV.exists():
        print(f"   !! {OBSCSV.name} not found — was the OBS package written/run? Skipping.")
        return
    df = pd.read_csv(OBSCSV)
    tcol = df.columns[0]
    t = df[tcol].to_numpy(dtype=float)

    # drop the spin-up portion of the record (model time < end of spin-up)
    spin = pd.date_range(SIM_START, periods=SPINUP_NPER + 1, freq="MS")
    spinup_days = (spin[-1] - spin[0]).days
    keep = t >= spinup_days - 1e-6
    dates = SIM_START + pd.to_timedelta(t[keep] - spinup_days, unit="D")

    # group columns by point name (strip _L<k>)
    pts = {}
    for c in df.columns[1:]:
        m = re.match(r"(.+)_L(\d+)$", c, re.IGNORECASE)
        name, lay = (m.group(1), int(m.group(2))) if m else (c, None)
        pts.setdefault(name, []).append((lay, c))

    # tidy long CSV
    tidy = []
    for name, cols in pts.items():
        for lay, c in cols:
            for d, v in zip(dates, df.loc[keep, c].to_numpy()):
                tidy.append((name, lay, d, v))
    pd.DataFrame(tidy, columns=["point", "layer", "date", "head_m"]) \
        .to_csv(OUT / "obs_head_timeseries.csv", index=False)

    # figure: one panel per point, a line per layer + land-surface reference
    mg, top, _ = load_grid()
    try:
        xc = np.array(mg.xcellcenters).ravel()
        yc = np.array(mg.ycellcenters).ravel()
    except Exception:
        xc = yc = None
    pt_node = {}
    if xc is not None:
        import geopandas as gpd
        from scipy.spatial import cKDTree
        gobs = gpd.read_file(GPKG, layer=OBS_LAYER).to_crs(3763)
        tree = cKDTree(np.column_stack([xc, yc]))
        ncol = "Name" if "Name" in gobs.columns else gobs.columns[0]
        for _, r in gobs.iterrows():
            if r.geometry is not None and not r.geometry.is_empty:
                _, idx = tree.query([r.geometry.x, r.geometry.y])
                pt_node[str(r[ncol]).strip().replace(" ", "_")] = int(idx)

    # monthly rainfall & ET0 over the plotted span -> bars hanging from a top axis
    months = pd.date_range(dates.min().to_period("M").to_timestamp(),
                           dates.max().to_period("M").to_timestamp(), freq="MS")
    rain = et = None
    try:
        _r = read_forcing(P_CSV).reindex(months).to_numpy(dtype=float)
        _e = read_forcing(ET_CSV).reindex(months).to_numpy(dtype=float)
        pe_max = float(np.nanmax(np.concatenate([_r, _e])))
        if not np.isfinite(pe_max):
            raise ValueError("no finite rainfall/ET in the plotted span")
        rain, et = _r, _e
    except Exception as e:
        print(f"   (rainfall/ET bars skipped: {e!r})")

    obs_meas = read_piezo_obs()    # observed heads per point (overlay + later PEST targets)

    # --- OUTLET virtual piezometers (outw/outs): they ARE columns in obs.head.csv (so they get
    #     panels), but they're not GPKG obs points, so (a) register their model cell in pt_node so
    #     the DEPTH plot has a land-surface reference, and (b) synthesise their "observed" series =
    #     the SAME virtual-P4 target make_obs calibrates to (outlet_top - P4_depth(t)). (2026-08-15) ---
    OUTLET_META = WORKSPACE / "outlet_cells.csv"
    if OUTLET_META.exists():
        _om = pd.read_csv(OUTLET_META)
        _p4d = None
        if Path(SYN_PIEZO_CSV).exists():
            _sy = pd.read_csv(SYN_PIEZO_CSV, parse_dates=["date"])
            _g4 = _sy[_sy["cdl"].astype(str).str.upper() == "P4"]
            if len(_g4):
                _p4d = _g4[["date", "depth_m"]].dropna().sort_values("date")
        for _, _r in _om.iterrows():
            _nm = str(_r["name"]).strip().upper()                      # OUTW / OUTS (MF6 uppercases the obs-csv column -> pts key)
            _nd = int(_r["node"]) - 1                                  # model wrote 1-based node
            pt_node[_nm] = _nd
            _tp = float(top[_nd]) if (top is not None and np.isfinite(top[_nd])) else float(_r["top"])
            if _p4d is not None:
                obs_meas[_nm] = pd.DataFrame({"date": _p4d["date"].to_numpy(),
                                              "head": _tp - _p4d["depth_m"].to_numpy()})

    names = list(pts.keys())
    n = len(names)

    # common HEAD-axis span across panels so amplitudes are directly comparable
    def _hrange(name):
        vals = []
        for lay, c in pts[name]:
            v = df.loc[keep, c].to_numpy(dtype=float); vals.append(v[np.isfinite(v)])
        if top is not None and name in pt_node and np.isfinite(top[pt_node[name]]):
            vals.append(np.array([top[pt_node[name]]]))
        if name in obs_meas and len(obs_meas[name]):
            om = obs_meas[name]
            _m = (om["date"] >= dates.min()) & (om["date"] <= dates.max())
            ov = om.loc[_m, "head"].to_numpy(dtype=float); vals.append(ov[np.isfinite(ov)])
        allv = np.concatenate([a for a in vals if a.size]) if any(a.size for a in vals) else np.array([])
        return (float(allv.min()), float(allv.max())) if allv.size else None
    _hr = {nm: _hrange(nm) for nm in names}
    _hspans = [hi - lo for r in _hr.values() if r for lo, hi in [r]]
    head_span = (max(_hspans) * 1.10) if _hspans else None

    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))

    # ---- shared-legend + x-axis helper (user 2026-08-15): ONE figure legend listing the layers
    #      present (+ land surface / observed / Rainfall,ET0), hosted in the empty trailing panel;
    #      mm-yyyy date labels (vertical) on the bottom-most VISIBLE panel of EACH column. ----
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    _lays = sorted({lay for cols in pts.values() for lay, _c in cols if lay is not None})
    def _finalize(fig2, axes2d, with_bars):
        hs = [Line2D([0], [0], color=HEAD_COLORS[(l - 1) % len(HEAD_COLORS)], lw=1.6, label=f"L{l}")
              for l in _lays]
        hs += [Line2D([0], [0], color="0.5", ls="--", lw=0.9, label="land surface"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="k", markeredgecolor="w",
                      markeredgewidth=0.4, markersize=6, label="observed")]
        if with_bars:
            hs += [Patch(facecolor="tab:blue", alpha=0.75, label="Rainfall"),
                   Patch(facecolor="tab:orange", alpha=0.75, label="ET$_0$")]
        flat = axes2d.ravel()
        host = flat[n] if n < nrows * ncols else None        # empty trailing panel -> legend host
        for ax in flat[(n + 1 if host is not None else n):]:
            ax.set_visible(False)
        if host is not None:
            host.set_visible(True); host.axis("off")
            host.legend(handles=hs, loc="center", fontsize=9, frameon=True,
                        title="Legend", title_fontsize=9)
        else:
            fig2.legend(handles=hs, loc="lower center", ncol=len(hs), fontsize=8,
                        bbox_to_anchor=(0.5, -0.01))
        for c in range(ncols):                                # bottom-most VISIBLE panel of each column
            rws = [r for r in range(nrows) if r * ncols + c < n]
            if not rws:
                continue
            axb = axes2d[max(rws), c]
            axb.xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))  # Oct 1 = hydro-year start
            axb.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
            axb.tick_params(axis="x", labelbottom=True, labelrotation=90, labelsize=7)

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 2.9 * nrows),
                             squeeze=False, sharex=True)
    for ax, name in zip(axes.ravel(), names):
        for j, (lay, c) in enumerate(sorted(pts[name])):
            col = HEAD_COLORS[((lay - 1) if lay else j) % len(HEAD_COLORS)]
            ax.plot(dates, df.loc[keep, c].to_numpy(), lw=1.3, color=col, zorder=5,
                    label=f"L{lay}" if lay else c)
        if top is not None and name in pt_node:
            ax.axhline(top[pt_node[name]], color="0.5", ls="--", lw=0.8,
                       label="land surface", zorder=4)
        if name in obs_meas and len(obs_meas[name]):       # observed heads (PEST targets)
            om = obs_meas[name]
            _m = (om["date"] >= dates.min()) & (om["date"] <= dates.max())
            ax.scatter(om.loc[_m, "date"], om.loc[_m, "head"], s=26, marker="o",
                       facecolor="k", edgecolor="w", linewidths=0.4, zorder=8,
                       label="observed")
        r = _hr.get(name)
        if head_span and r:                    # identical head-axis span, centred per panel
            mid = 0.5 * (r[0] + r[1]); ax.set_ylim(mid - head_span / 2, mid + head_span / 2)
        ax.set_ylabel("head (m)")
        ax2 = None
        if rain is not None:
            ax2 = ax.twinx()
            w = 12.0                                  # bar width (days); two per month
            mid = months + pd.Timedelta(days=14)
            ax2.bar(mid - pd.Timedelta(days=w / 2), rain, width=w,
                    color="tab:blue", alpha=0.75, label="Rainfall")
            ax2.bar(mid + pd.Timedelta(days=w / 2), et, width=w,
                    color="tab:orange", alpha=0.75, label="ET$_0$")
            ax2.set_ylim(pe_max * 3.2, 0.0)           # 0 at top -> bars descend
            ax2.set_ylabel("P, ET$_0$ (mm/mo)", fontsize=8)
            ax2.tick_params(labelsize=7)
            ax.set_zorder(ax2.get_zorder() + 1)       # head lines in front of bars
            ax.patch.set_visible(False)
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3)
    _finalize(fig, axes, with_bars=(rain is not None))
    fig.suptitle("CdL — simulated head at observation points  (rainfall & ET$_0$ on top; post spin-up; "
                 + (f"identical y-span = {head_span:.2f} m)" if head_span else "post spin-up)"), y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "obs_head_timeseries.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote obs_head_timeseries.png / .csv  ({n} points)")

    # ---- companion: DEPTH-to-water-table (user 2026-08-14) — head expressed as depth below
    #      the model land surface (depth = top - head), y-axis INVERTED (0 = ground at the top,
    #      deeper downward), so the WT depth can be sanity-checked spatially. Shallowest layer =
    #      the water table. Observed depth = top - observed head (same land-surface reference).
    if pt_node:
        tidyd = []
        figd, axesd = plt.subplots(nrows, ncols, figsize=(7 * ncols, 2.9 * nrows),
                                   squeeze=False, sharex=True)
        for ax, name in zip(axesd.ravel(), names):
            g = float(top[pt_node[name]]) if name in pt_node else np.nan
            _dv = []
            for j, (lay, c) in enumerate(sorted(pts[name])):
                col = HEAD_COLORS[((lay - 1) if lay else j) % len(HEAD_COLORS)]
                depth = g - df.loc[keep, c].to_numpy()
                _dv.append(depth[np.isfinite(depth)])
                ax.plot(dates, depth, lw=1.3, color=col, zorder=5, label=f"L{lay}" if lay else c)
                for d, v in zip(dates, depth):
                    tidyd.append((name, lay, d, v))
            ax.axhline(0.0, color="0.5", ls="--", lw=0.8, label="land surface", zorder=4)
            if name in obs_meas and len(obs_meas[name]) and np.isfinite(g):
                om = obs_meas[name]
                _m = (om["date"] >= dates.min()) & (om["date"] <= dates.max())
                _od = g - om.loc[_m, "head"].to_numpy()
                _dv.append(_od[np.isfinite(_od)])
                ax.scatter(om.loc[_m, "date"], _od, s=26, marker="o",
                           facecolor="k", edgecolor="w", linewidths=0.4, zorder=8, label="observed")
            # identical y-SPAN across panels (= the head plot's head_span; depth amplitude ==
            # head amplitude) so the fluctuation amplitude is directly comparable panel-to-panel.
            _allv = np.concatenate([a for a in _dv if a.size]) if any(a.size for a in _dv) else None
            if head_span and _allv is not None and _allv.size:
                _mid = 0.5 * (float(_allv.min()) + float(_allv.max()))
                ax.set_ylim(_mid - head_span / 2, _mid + head_span / 2)
            if rain is not None:                       # Rainfall & ET0 bars (same as the head plot)
                ax2 = ax.twinx()
                _w = 12.0; _mid = months + pd.Timedelta(days=14)
                ax2.bar(_mid - pd.Timedelta(days=_w / 2), rain, width=_w, color="tab:blue", alpha=0.75)
                ax2.bar(_mid + pd.Timedelta(days=_w / 2), et, width=_w, color="tab:orange", alpha=0.75)
                ax2.set_ylim(pe_max * 3.2, 0.0)        # 0 at top -> bars descend
                ax2.set_ylabel("P, ET$_0$ (mm/mo)", fontsize=8); ax2.tick_params(labelsize=7)
                ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
            ax.set_ylabel("depth to WT (m)"); ax.set_title(name, fontsize=10); ax.grid(alpha=0.3)
        for ax in axesd.ravel()[:n]:
            ax.invert_yaxis()                          # depth increases downward; 0 (ground) at the top
        _finalize(figd, axesd, with_bars=(rain is not None))
        figd.suptitle("CdL — simulated DEPTH to water table at obs points  (depth = land surface − head; "
                      "0 = ground; " + (f"identical y-span = {head_span:.2f} m)" if head_span else "post spin-up)"),
                      y=1.0)
        figd.tight_layout()
        figd.savefig(OUT / "obs_depth_timeseries.png", dpi=150, bbox_inches="tight")
        plt.close(figd)
        pd.DataFrame(tidyd, columns=["point", "layer", "date", "depth_m"]).to_csv(
            OUT / "obs_depth_timeseries.csv", index=False)
        print(f"   wrote obs_depth_timeseries.png / .csv  ({n} points)")


# =============================================================================
# 2. MEAN-HEAD & MEAN-DEPTH MAPS PER LAYER
# =============================================================================
def compute_mean_head():
    """Mean head per layer (nlay, ncpl) over the post-spin-up periods, + kstpkper used."""
    hds = flopy.utils.HeadFile(str(HDS))
    kk = real_kstpkper(hds) or hds.get_kstpkper()
    arr = np.array([np.squeeze(hds.get_data(kstpkper=k)) for k in kk])  # (nt,nlay,ncpl)
    if arr.ndim == 2:                       # single layer -> (nt,ncpl)
        arr = arr[:, None, :]
    arr = np.where(np.abs(arr) > 1e29, np.nan, arr)
    return np.nanmean(arr, axis=0), kk      # (nlay,ncpl), list


def load_overlays():
    """(obs_points gdf, routed-streams gdf) for map overlays; None where unavailable."""
    gobs = routed = None
    try:
        import geopandas as gpd
        gobs = gpd.read_file(GPKG, layer=OBS_LAYER).to_crs(3763)
    except Exception:
        pass
    try:
        import geopandas as gpd
        if ROUTED_SHP.exists():
            routed = gpd.read_file(ROUTED_SHP)
    except Exception:
        pass
    return gobs, routed


def _plot_layer_maps(arrays, titles, cbar_label, cmap, suptitle, outfile,
                     mg, gobs, routed):
    """Shared per-layer map figure (one panel per array) used by head & depth maps."""
    n = len(arrays)
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), squeeze=False)
    axf = axes.ravel()
    for i, (a, t) in enumerate(zip(arrays, titles)):
        ax = axf[i]
        pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax, layer=0)
        ca = pmv.plot_array(a, cmap=cmap)
        pmv.plot_grid(lw=0.1, color="0.85")
        if routed is not None:
            routed.plot(ax=ax, color="white", lw=0.8, zorder=3)
        if gobs is not None:
            _ncol = "Name" if "Name" in gobs.columns else gobs.columns[0]
            ax.scatter(gobs.geometry.x, gobs.geometry.y, marker="s", s=40,
                       c="magenta", edgecolors="k", linewidths=0.6, zorder=5)
            for _, _r in gobs.iterrows():
                if _r.geometry is None or _r.geometry.is_empty:
                    continue
                ax.annotate(str(_r[_ncol]), (_r.geometry.x, _r.geometry.y),
                            textcoords="offset points", xytext=(5, 4),
                            fontsize=8, fontweight="bold", color="k", zorder=6,
                            path_effects=[patheffects.withStroke(linewidth=2,
                                                                 foreground="white")])
        fig.colorbar(ca, ax=ax, shrink=0.7, label=cbar_label)
        ax.set_title(t)
        ax.set_aspect("equal")
    for ax in axf[n:]:
        ax.set_visible(False)
    fig.suptitle(suptitle, y=1.0)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def head_maps():
    print(">> [2] Mean-head maps per layer …")
    if not HDS.exists():
        print(f"   !! {HDS.name} not found. Skipping.")
        return
    mean_head, kk = compute_mean_head()
    nlay = mean_head.shape[0]
    mg, top, _ = load_grid()
    gobs, routed = load_overlays()
    pd.DataFrame({"layer": np.arange(1, nlay + 1),
                  "mean_head_m": np.nanmean(mean_head, axis=1)}) \
        .to_csv(OUT / "mean_head_per_layer.csv", index=False)
    _plot_layer_maps(
        [mean_head[l] for l in range(nlay)],
        [f"Mean head — layer {l + 1}" for l in range(nlay)],
        "head (m)", "viridis",
        f"CdL — mean head over {len(kk)} post-spin-up periods",
        OUT / "mean_head_maps.png", mg, gobs, routed)
    print(f"   wrote mean_head_maps.png / mean_head_per_layer.csv  ({nlay} layers)")


def depth_maps():
    """Per-layer mean DEPTH = land surface (topography) − head, mirroring head_maps."""
    print(">> [2b] Mean-depth maps per layer (land surface − head) …")
    if not HDS.exists():
        print(f"   !! {HDS.name} not found. Skipping.")
        return
    mg, top, _ = load_grid()
    if top is None:
        print(f"   !! no top elevation ({GRB.name} missing) — cannot compute depth. Skipping.")
        return
    mean_head, kk = compute_mean_head()
    nlay = mean_head.shape[0]
    gobs, routed = load_overlays()
    depth = np.array([top - mean_head[l] for l in range(nlay)])   # + = head below ground
    pd.DataFrame({"layer": np.arange(1, nlay + 1),
                  "mean_depth_m": np.nanmean(depth, axis=1)}) \
        .to_csv(OUT / "mean_depth_per_layer.csv", index=False)
    _plot_layer_maps(
        [depth[l] for l in range(nlay)],
        [f"Mean depth (surface − head) — layer {l + 1}" for l in range(nlay)],
        "depth below ground (m)", "RdYlBu_r",
        f"CdL — mean depth to head over {len(kk)} post-spin-up periods\n"
        f"(+ = head below land surface,  − = above)",
        OUT / "mean_depth_maps.png", mg, gobs, routed)
    print(f"   wrote mean_depth_maps.png / mean_depth_per_layer.csv  ({nlay} layers)")


# =============================================================================
# 3a. COMPARTMENT BUDGET FROM THE LISTING (overall mean + yearly)
# =============================================================================
def list_budget():
    print(">> [3a] Compartment water budget (listing) …")
    if not LST.exists():
        print(f"   !! {LST.name} not found. Skipping.")
        return
    mflist = flopy.utils.Mf6ListBudget(str(LST))
    df_flux, _ = mflist.get_dataframes(start_datetime=str(SIM_START.date()), diff=False)
    if df_flux is None or len(df_flux) == 0:
        print("   !! no budget parsed. Skipping.")
        return
    flux = df_flux.iloc[SPINUP_NPER:].copy()           # drop SS + spin-up
    n_real = len(flux)
    flux.index = real_dates(n_real)

    drop = {"TOTAL", "PERCENT_DISCREPANCY", "IN-OUT"}
    bases = sorted({c[:-3] for c in flux.columns
                    if c.endswith("_IN") and c[:-3] not in drop})
    net = pd.DataFrame(
        {b: flux.get(b + "_IN", 0.0) - flux.get(b + "_OUT", 0.0) for b in bases},
        index=flux.index)                              # + = source INTO aquifer

    # ---- overall mean rate (m3/d) ----
    summ = pd.DataFrame({
        "term": bases,
        "label": [term_label(b) for b in bases],
        "compartment": [compartment_of(b) for b in bases],
        "mean_rate_m3d": [net[b].mean() for b in bases],
    }).sort_values(["compartment", "term"])
    summ.to_csv(OUT / "budget_compartment_mean.csv", index=False)

    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    labels = [term_label(b) for b in summ["term"]]
    colors = [term_color(b) for b in summ["term"]]
    bars = ax.barh(labels, summ["mean_rate_m3d"], color=colors, edgecolor="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.8)
    for bar, v in zip(bars, summ["mean_rate_m3d"]):     # value labels (small terms readable)
        ax.annotate(f"{v:,.0f}", (v, bar.get_y() + bar.get_height() / 2),
                    xytext=(4 if v >= 0 else -4, 0), textcoords="offset points",
                    va="center", ha="left" if v >= 0 else "right", fontsize=7)
    ax.set_xlabel("mean rate (m³/d)   + = source into aquifer,  − = sink")
    ax.set_title("CdL — mean water budget by term (post spin-up)\n"
                 "DRN family in blue tones (seep / seep→SFR baseflow / west underflow / secondary outlet)")
    ax.margins(x=0.18)
    fig.tight_layout()
    fig.savefig(OUT / "budget_compartment_mean.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- yearly volumes (m3/yr) ----
    days = net.index.days_in_month.to_numpy(dtype=float)
    vol = net.multiply(days, axis=0)                   # m3 per period
    yearly = vol.groupby(vol.index.year).sum()         # year-agnostic across pandas versions
    yearly.to_csv(OUT / "budget_yearly_volume_m3.csv")

    if len(yearly) >= 1:
        fig, ax = plt.subplots(figsize=(10, 5))
        bottom_pos = np.zeros(len(yearly)); bottom_neg = np.zeros(len(yearly))
        x = np.arange(len(yearly))
        for b in bases:
            vals = yearly[b].to_numpy()
            base = np.where(vals >= 0, bottom_pos, bottom_neg)
            ax.bar(x, vals, bottom=base, label=term_label(b),
                   color=term_color(b), edgecolor="k", lw=0.3)
            bottom_pos += np.where(vals >= 0, vals, 0)
            bottom_neg += np.where(vals < 0, vals, 0)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(yearly.index, rotation=90, va="top")
        ax.set_ylabel("volume (m³/yr)   + = source into aquifer")
        ax.set_title("CdL — yearly water budget by term  (DRN family in blue tones)")
        ax.legend(fontsize=6.5, ncol=2, loc="best")
        fig.tight_layout()
        fig.savefig(OUT / "budget_yearly_volume.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"   wrote budget_compartment_mean.* and budget_yearly_volume.*  "
          f"({n_real} periods, {len(yearly)} yr)")


# =============================================================================
# 3b. PACKAGE-INTERNAL BUDGETS (UZF = unsaturated, SFR = surface water)
# =============================================================================
def package_budget(bud_path, title, fname):
    if not bud_path.exists():
        print(f"   !! {bud_path.name} not found. Skipping {title}.")
        return
    print(f">> [3b] {title} internal budget ({bud_path.name}) …")
    cbf = flopy.utils.CellBudgetFile(str(bud_path))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    raw = cbf.get_unique_record_names()
    means = {}
    for rname in raw:
        label = (rname.decode() if isinstance(rname, bytes) else rname).strip()
        if label.upper().startswith("FLOW-JA-FACE"):
            continue
        tot = []
        for k in kk:
            try:
                data = cbf.get_data(text=rname, kstpkper=k)
            except Exception:
                continue
            if not data:
                continue
            d = data[0]
            if hasattr(d, "dtype") and d.dtype.names and "q" in d.dtype.names:
                tot.append(float(np.sum(d["q"])))
            else:
                tot.append(float(np.nansum(np.asarray(d))))
        if tot:
            means[label] = np.mean(tot)
    cbf.close()
    if not means:
        print(f"   !! no records summed for {title}.")
        return
    # NET mover flux: FROM-MVR (+, received) and *-TO-MVR (-, sent) are two sides of
    # the SAME routed water — in a CRR cascade the same parcel is provider for one
    # cell and receiver for the next, so summing them gross double-counts and can read
    # as "fake water" (e.g. FROM-MVR ~= 40% of infiltration here). Report their signed
    # net so the figure is not misread as inflated recharge. See mvr_accounting() for
    # the cross-package split (internal recirculation vs export to SFR/LAK).
    recv = sum(v for n, v in means.items() if "FROM-MVR" in n.upper())
    sent = sum(v for n, v in means.items() if n.upper().endswith("-TO-MVR"))
    if recv or sent:
        means["MVR-NET (recv+sent)"] = recv + sent
    s = pd.Series(means).sort_values()
    s.to_frame("mean_rate_m3d").to_csv(OUT / f"{fname}.csv")
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(s) + 1.5))
    ax.barh(s.index, s.values,
            color=["tab:red" if v < 0 else "tab:blue" for v in s.values],
            edgecolor="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean rate (m³/d)   (MF6 sign: + into the package element)")
    ax.set_title(f"CdL — {title} budget (post spin-up)")
    fig.tight_layout()
    fig.savefig(OUT / f"{fname}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote {fname}.png / .csv  ({len(s)} terms)")


# =============================================================================
# 3b-bis. SFR SURFACE-WATER BUDGET — relabelled overview + inlet/outlet split
# =============================================================================
SFR_LABELS = {
    "EXT-INFLOW":  "spec. inflow (inlet, node 3506)",
    "EXT-OUTFLOW": "outflow at the outlets (total)",
    "FROM-MVR":    "baseflow in (from DRN-seep)",
    "TO-MVR":      "to ponds (LAK, via MVR)",
    "GWF":         "stream↔aquifer leakage",
    "RAINFALL":    "rainfall on channel",
    "EVAPORATION": "evaporation",
    "RUNOFF":      "runoff",
    "STORAGE":     "channel storage",
}


def _reach_means(cbf, text, kk):
    """Mean per-reach q for an SFR budget text -> dict{reach(1-based): mean_q}."""
    acc, cnt = {}, 0
    for k in kk:
        try:
            data = cbf.get_data(text=text, kstpkper=k)
        except Exception:
            continue
        if not data:
            continue
        d = data[0]
        if not (hasattr(d, "dtype") and d.dtype.names
                and "node" in d.dtype.names and "q" in d.dtype.names):
            continue
        for node, q in zip(d["node"], d["q"]):
            acc[int(node)] = acc.get(int(node), 0.0) + float(q)
        cnt += 1
    if cnt:
        for n in list(acc):
            acc[n] /= cnt
    return acc


def sfr_budget():
    if not SFRB.exists():
        print(f"   !! {SFRB.name} not found. Skipping SFR budget.")
        return
    print(">> [3b] surface water (SFR) budget + inlet/outlet split …")
    cbf = flopy.utils.CellBudgetFile(str(SFRB))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    # ---- (a) relabelled internal overview (replaces the generic sfr_budget_mean) ----
    raw = [(r.decode() if isinstance(r, bytes) else r).strip()
           for r in cbf.get_unique_record_names()]
    means = {}
    for label in raw:
        if label.upper().startswith("FLOW-JA-FACE"):
            continue
        m = _reach_means(cbf, label, kk)
        if m:
            means[label] = sum(m.values())
    s = pd.Series(means).sort_values()
    (s.rename(index=lambda t: SFR_LABELS.get(t, t))
       .to_frame("mean_rate_m3d").to_csv(OUT / "sfr_budget_mean.csv"))
    fig, ax = plt.subplots(figsize=(8.5, 0.5 * len(s) + 1.6))
    ax.barh([SFR_LABELS.get(t, t) for t in s.index], s.values,
            color=["tab:red" if v < 0 else "tab:blue" for v in s.values], edgecolor="k", lw=0.4)
    ax.set_xscale("symlog", linthresh=10)
    ax.axvline(0, color="k", lw=0.8)
    for i, v in enumerate(s.values):
        ax.annotate(f"{v:,.0f}", (v, i), xytext=(4, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=7)
    ax.set_xlabel("mean rate (m³/d, symlog)   blue = into stream, red = out of stream")
    ax.set_title("CdL — SFR surface-water budget (post spin-up)")
    fig.tight_layout(); fig.savefig(OUT / "sfr_budget_mean.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- (b) inlet / main-outlet / 2nd-outlet split ----
    inflow = _reach_means(cbf, "EXT-INFLOW", kk)
    outflow = _reach_means(cbf, "EXT-OUTFLOW", kk)
    frommvr = sum(_reach_means(cbf, "FROM-MVR", kk).values())
    tomvr   = sum(_reach_means(cbf, "TO-MVR", kk).values())
    gwf     = sum(_reach_means(cbf, "GWF", kk).values())
    cbf.close()
    inl = [(n, q) for n, q in inflow.items() if abs(q) > 1e-9]
    Q_in = sum(q for _, q in inl)
    r_in = max(inl, key=lambda t: abs(t[1]))[0] if inl else None
    outs = sorted([(n, q) for n, q in outflow.items() if abs(q) > 1e-9],
                  key=lambda t: abs(t[1]), reverse=True)
    r_main, q_main = (outs[0] if outs else (None, 0.0))     # (reach, flow)
    r_sec,  q_sec  = (outs[1] if len(outs) > 1 else (None, 0.0))
    q_out_tot = sum(q for _, q in outs)
    net = Q_in + q_out_tot + frommvr + tomvr + gwf
    rows = [
        (f"inlet — spec. inflow (reach {r_in})", Q_in),
        (f"main outlet (reach {r_main})", q_main),
        (f"2nd outlet — node 2 (reach {r_sec})", q_sec),
        ("outlets total (main + 2nd)", q_out_tot),
        ("baseflow in (from DRN-seep)", frommvr),
        ("to ponds (LAK)", tomvr),
        ("stream↔aquifer (GWF)", gwf),
        ("net balance (Σ)", net),
    ]
    pd.DataFrame(rows, columns=["term", "mean_rate_m3d"]).to_csv(
        OUT / "sfr_inlet_outlet.csv", index=False)
    labels = [r[0] for r in rows][::-1]
    vals = [r[1] for r in rows][::-1]
    cols = ["0.6" if "net" in l else ("tab:blue" if v >= 0 else "tab:red")
            for l, v in zip(labels, vals)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(labels, vals, color=cols, edgecolor="k", lw=0.4)
    ax.set_xscale("symlog", linthresh=10)
    ax.axvline(0, color="k", lw=0.8)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:,.0f}", (v, i), xytext=(4, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=8)
    ax.set_xlabel("mean rate (m³/d, symlog)   blue = into stream, red = out")
    ax.set_title("CdL — SFR inlet / outlets split (post spin-up)")
    fig.tight_layout(); fig.savefig(OUT / "sfr_inlet_outlet.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote sfr_budget_mean.* and sfr_inlet_outlet.*  "
          f"(inlet reach {r_in}, main {r_main}, 2nd {r_sec})")


# =============================================================================
# 3c. PER-LAYER STORAGE CHANGE (aquifer layers) FROM THE GWF .cbc
# =============================================================================
def layer_storage():
    print(">> [3c] Per-layer storage change (GWF .cbc) …")
    if not CBC.exists():
        print(f"   !! {CBC.name} not found. Skipping.")
        return
    cbf = flopy.utils.CellBudgetFile(str(CBC))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    names = [(r.decode() if isinstance(r, bytes) else r).strip()
             for r in cbf.get_unique_record_names()]
    sto_terms = [n for n in names if n.upper().startswith("STO")]
    if not sto_terms:
        print("   !! no STO records (steady-state only run?). Skipping.")
        cbf.close()
        return
    per_layer = None
    for k in kk:
        tot = None
        for term in sto_terms:
            d = np.squeeze(cbf.get_data(text=term, kstpkper=k)[0])  # (nlay,ncpl)
            if d.ndim == 1:
                d = d[None, :]
            tot = d if tot is None else tot + d
        lay_sum = np.nansum(tot, axis=1)                            # (nlay,)
        per_layer = lay_sum[None, :] if per_layer is None else np.vstack([per_layer, lay_sum])
    cbf.close()
    mean_lay = np.nanmean(per_layer, axis=0)
    nlay = mean_lay.size
    pd.DataFrame({"layer": np.arange(1, nlay + 1), "mean_storage_rate_m3d": mean_lay}) \
        .to_csv(OUT / "layer_storage_mean.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(np.arange(1, nlay + 1), mean_lay, color="tab:brown", edgecolor="k")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("layer"); ax.set_ylabel("mean storage rate (m³/d)")
    ax.set_title("CdL — per-layer net storage change (+ = released to flow)")
    ax.set_xticks(np.arange(1, nlay + 1))
    fig.tight_layout()
    fig.savefig(OUT / "layer_storage_mean.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote layer_storage_mean.png / .csv  ({nlay} layers)")


# =============================================================================
# 4. LAKE (POND) STAGE + SEEPAGE
# =============================================================================
def lake_outputs():
    print(">> [4] Lake (pond) stage + seepage …")
    STG = WORKSPACE / f"{MODEL_NAME}.lak.stage"
    LKB = WORKSPACE / f"{MODEL_NAME}.lak.bud"
    if not STG.exists() and not LKB.exists():
        print("   !! no LAK output (no ponds / not run). Skipping.")
        return
    # --- stage time series per lake ---
    if STG.exists():
        try:
            sf = flopy.utils.HeadFile(str(STG), text="STAGE")
            kk = [k for k in sf.get_kstpkper() if k[1] >= SPINUP_NPER] or sf.get_kstpkper()
            stages = np.array([np.ravel(sf.get_data(kstpkper=k)) for k in kk])    # (nt, nlakes)
            stages = np.where(np.abs(stages) > 1e29, np.nan, stages)   # mask MF6 dry-lake sentinel (pond dry)
            d = real_dates(len(kk))
            nlk = stages.shape[1]
            out = pd.DataFrame(stages, columns=[f"lake{L}" for L in range(nlk)])
            out.insert(0, "date", d)
            out.to_csv(OUT / "lake_stage.csv", index=False)
            fig, ax = plt.subplots(figsize=(10, 4))
            for L in range(nlk):
                ax.plot(d, stages[:, L], lw=1.3, label=f"lake {L}")
            ax.set_ylabel("lake stage (m)"); ax.set_title("CdL — pond (LAK) stage")
            ax.legend(fontsize=8, ncol=min(nlk, 5)); ax.grid(alpha=0.3)
            fig.tight_layout(); fig.savefig(OUT / "lake_stage.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"   wrote lake_stage.png / .csv  ({nlk} lakes)")
        except Exception as e:
            print(f"   (lake stage skipped: {e!r})")
        # --- per-pond stage panels with bottom + spill-invert datums (2026-07-04) -----------
        # metadata parsed from the WRITTEN cdl_gwf.lak (boundnames pondN, outlet inverts) and
        # the .lakN.tab files (bottom = min table stage) -> always consistent with the run.
        try:
            names, inverts, bottoms = {}, {}, {}
            blk = None
            for ln in (WORKSPACE / f"{MODEL_NAME}.lak").read_text().splitlines():
                s = ln.strip(); low = s.lower()
                if low.startswith("begin "):
                    blk = low.split()[1]; continue
                if low.startswith("end "):
                    blk = None; continue
                if not s or s.startswith("#") or blk is None:
                    continue
                p = s.split()
                if blk == "packagedata":
                    names[int(p[0]) - 1] = p[-1]
                elif blk == "outlets":
                    inverts[int(p[1]) - 1] = float(p[4])
            tab_stage, tab_vol = {}, {}      # per-lake stage-volume table (cols 0,1 of .lakN.tab)
            for L in range(nlk):
                tab = WORKSPACE / f"{MODEL_NAME}.lak{L + 1}.tab"
                if tab.exists():
                    stg, vol = [], []
                    for tl in tab.read_text().splitlines():
                        p = tl.split()
                        try:
                            stg.append(float(p[0])); vol.append(float(p[1]))
                        except (ValueError, IndexError):
                            pass
                    if stg:
                        bottoms[L] = min(stg)
                        order = np.argsort(stg)
                        tab_stage[L] = np.asarray(stg)[order]
                        tab_vol[L] = np.asarray(vol)[order]

            # common y-SPAN across panels so stage AMPLITUDES are directly comparable
            def _content(L):
                ys = [np.nanmin(stages[:, L]), np.nanmax(stages[:, L])]
                if L in bottoms:
                    ys.append(bottoms[L])
                if L in inverts:
                    ys.append(inverts[L])
                return [y for y in ys if np.isfinite(y)]
            spans = [max(c) - min(c) for c in (_content(L) for L in range(nlk)) if c]
            common = (max(spans) * 1.08) if spans else 1.0
            if common <= 0:
                common = 1.0
            ncols = 4
            nrows = int(np.ceil(nlk / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.1 * nrows),
                                     sharex=True, squeeze=False)
            for L in range(nlk):
                ax = axes.ravel()[L]
                ax.plot(d, stages[:, L], lw=0.7, color="tab:blue")
                if L in bottoms:
                    ax.axhline(bottoms[L], color="saddlebrown", lw=0.8, ls="--")
                if L in inverts:
                    ax.axhline(inverts[L], color="crimson", lw=0.8, ls=":")
                c = _content(L)
                if c:                                    # identical span, centred on this pond
                    mid = 0.5 * (max(c) + min(c))
                    ax.set_ylim(mid - common / 2, mid + common / 2)
                ttl = names.get(L, f"lake {L}")
                if L in bottoms and L in inverts:
                    ttl += f"  (bottom {bottoms[L]:.1f}, spill {inverts[L]:.1f})"
                ax.set_title(ttl, fontsize=9)
                ax.tick_params(labelsize=7)
            for ax in axes.ravel()[nlk:]:
                ax.axis("off")
            fig.suptitle("CdL — pond (LAK) stages  (brown dashes = pond bottom, red dots = spill invert; "
                         f"identical y-span = {common:.2f} m for amplitude comparison)")
            fig.tight_layout()
            fig.savefig(OUT / "lake_stages_panels.png", dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"   wrote lake_stages_panels.png ({nlk} ponds, y-span {common:.2f} m)")

            # --- total pond STORAGE (volume) over time, from the stage-volume tables ---
            if tab_vol:
                per_lake, volt = {}, np.zeros(len(d))
                for L in range(nlk):
                    if L in tab_vol:
                        vL = np.interp(stages[:, L], tab_stage[L], tab_vol[L],
                                       left=tab_vol[L][0], right=tab_vol[L][-1])
                        vL = np.where(np.isfinite(stages[:, L]), vL, 0.0)   # dry pond -> 0 m3
                        per_lake[names.get(L, f"lake{L}")] = vL
                        volt = volt + vL
                dfv = pd.DataFrame(per_lake); dfv.insert(0, "date", d)
                dfv["total_m3"] = volt
                dfv.to_csv(OUT / "lake_volume.csv", index=False)
                fig, ax = plt.subplots(figsize=(10, 4))
                ax.fill_between(d, 0, volt, color="tab:blue", alpha=0.30)
                ax.plot(d, volt, lw=1.5, color="tab:blue")
                ax.set_ylim(bottom=0)
                ax.set_ylabel("total pond storage (m³)")
                ax.set_title("CdL — total pond (LAK) storage volume over time   "
                             f"(mean {volt.mean():,.0f} m³, range {volt.min():,.0f}–{volt.max():,.0f} m³, "
                             f"ΔV {volt.max() - volt.min():,.0f} m³)")
                ax.grid(alpha=0.3)
                fig.tight_layout(); fig.savefig(OUT / "lake_volume.png", dpi=150, bbox_inches="tight")
                plt.close(fig)
                print(f"   wrote lake_volume.png / .csv  (mean {volt.mean():,.0f} m³, "
                      f"ΔV {volt.max() - volt.min():,.0f} m³)")
        except Exception as e:
            print(f"   (lake stage panels / volume skipped: {e!r})")
    # --- lake budget terms ('GWF' = lakebed seepage = the SW-GW flux) ---
    if LKB.exists():
        try:
            cbf = flopy.utils.CellBudgetFile(str(LKB))
            kk = [k for k in cbf.get_kstpkper() if k[1] >= SPINUP_NPER] or cbf.get_kstpkper()
            recs = [(r.decode() if isinstance(r, bytes) else r).strip()
                    for r in cbf.get_unique_record_names()]
            summary = {}
            for r in recs:
                if r.upper().startswith("FLOW-JA-FACE"):
                    continue
                tot = []
                for k in kk:
                    try:
                        data = cbf.get_data(text=r, kstpkper=k)
                    except Exception:
                        continue
                    if data and hasattr(data[0], "dtype") and data[0].dtype.names \
                            and "q" in data[0].dtype.names:
                        tot.append(float(np.sum(data[0]["q"])))
                if tot:
                    summary[r] = np.mean(tot)
            cbf.close()
            if summary:
                s = pd.Series(summary).sort_values()
                s.to_frame("mean_rate_m3d").to_csv(OUT / "lake_budget_mean.csv")
                fig, ax = plt.subplots(figsize=(8, 0.5 * len(s) + 1.5))
                ax.barh(s.index, s.values,
                        color=["tab:red" if v < 0 else "tab:blue" for v in s.values], edgecolor="k")
                ax.axvline(0, color="k", lw=0.8)
                ax.set_xlabel("mean rate (m³/d)   (MF6 sign: + into the lake)")
                ax.set_title("CdL — pond (LAK) budget  ('GWF' = bed seepage = SW–GW flux)")
                fig.tight_layout(); fig.savefig(OUT / "lake_budget_mean.png", dpi=150, bbox_inches="tight")
                plt.close(fig)
                print(f"   wrote lake_budget_mean.png / .csv  ({len(s)} terms)")
        except Exception as e:
            print(f"   (lake budget skipped: {e!r})")


# =============================================================================
# 3c. MVR / CRR MOVER ACCOUNTING — net flux + provider->receiver topology
# =============================================================================
def _mvr_type(pname):
    """map an MVR package name to a compartment label."""
    u = pname.upper()
    if u.startswith("UZF"):      return "UZF"
    if u.startswith("SFR"):      return "SFR"
    if u.startswith("LAK"):      return "LAK"
    if "SEEP" in u:              return "DRN-SEEP"
    if u.startswith("DRN"):      return "DRN-OUT"
    return pname


def mvr_accounting():
    """Cross-package mover (MVR / CRR) flow accounting.

    FROM-MVR (water received from the mover, +) and *-TO-MVR (water given to the
    mover, -) are opposite sides of the SAME routed water.  Summed gross across a
    cascade they double-count (the classic 'fake water' that can exceed rainfall).
    This reports, per package, water GIVEN to and RECEIVED from the mover, the NET,
    and the provider->receiver topology from the static .mvr file so that internal
    UZF->UZF recirculation is distinguished from real export to SFR / LAK.
    """
    print(">> [3c] MVR / CRR mover accounting …")

    def mvr_terms(bud_path):
        res = {}
        if not bud_path.exists():
            return res
        cbf = flopy.utils.CellBudgetFile(str(bud_path))
        kk = real_kstpkper(cbf) or cbf.get_kstpkper()
        for rn in cbf.get_unique_record_names():
            lbl = (rn.decode() if isinstance(rn, bytes) else rn).strip()
            if "MVR" not in lbl.upper():
                continue
            tot = []
            for k in kk:
                try:
                    d = cbf.get_data(text=rn, kstpkper=k)
                except Exception:
                    continue
                if not d:
                    continue
                a = d[0]
                if hasattr(a, "dtype") and a.dtype.names and "q" in a.dtype.names:
                    tot.append(float(np.sum(a["q"])))
                else:
                    tot.append(float(np.nansum(np.asarray(a))))
            if tot:
                res[lbl] = np.mean(tot)
        cbf.close()
        return res

    # ---- (a) per-package given/received from the advanced-package budgets ----
    rows = []
    for pkg, bp in [("UZF", UZFB), ("SFR", SFRB), ("LAK", WORKSPACE / f"{MODEL_NAME}.lak.bud")]:
        for term, val in mvr_terms(bp).items():
            direction = "received (from mover)" if "FROM-MVR" in term.upper() else "given (to mover)"
            rows.append([pkg, term, round(val, 1), direction])
    # DRN mover flux from the listing budget (DRN is a GWF-model package)
    try:
        mflist = flopy.utils.Mf6ListBudget(str(LST))
        dff, _ = mflist.get_dataframes(start_datetime=str(SIM_START.date()), diff=False)
        real = dff.iloc[SPINUP_NPER:] if len(dff) > SPINUP_NPER else dff
        for col in dff.columns:
            cu = col.upper()
            # listing splits each term into _IN/_OUT; DRN->mover is an OUT flux
            if "MVR" in cu and "DRN" in cu and cu.endswith("_OUT"):
                val = -abs(float(real[col].mean()))
                if abs(val) < 0.05:                      # skip inactive DRN-OUT mover
                    continue
                rows.append(["DRN", col[:-4], round(val, 1), "given (to mover)"])
    except Exception as e:
        print(f"   (listing DRN-MVR skipped: {e!r})")

    if not rows:
        print("   !! no MVR terms found — CRR not active? Skipping.")
        return
    df = pd.DataFrame(rows, columns=["package", "term", "mean_rate_m3d", "direction"])
    df.to_csv(OUT / "mvr_accounting.csv", index=False)

    given_tot = -sum(v for *_, v, d in rows if d.startswith("given"))     # positive magnitude
    recv_tot  =  sum(v for *_, v, d in rows if d.startswith("received"))
    print(f"   total GIVEN to mover  = {given_tot:>10.1f} m3/d")
    print(f"   total RECEIVED        = {recv_tot:>10.1f} m3/d  (should match -> mover conserves mass)")
    print(f"   imbalance             = {recv_tot - given_tot:>10.1f} m3/d")

    # ---- (b) provider->receiver topology from the static .mvr file ----
    import collections
    mvr_file = WORKSPACE / f"{MODEL_NAME}.mvr"
    topo = collections.Counter()
    if mvr_file.exists():
        inp = False
        for line in mvr_file.read_text().splitlines():
            s = line.strip()
            if s.upper().startswith("BEGIN PERIOD"):
                inp = True; continue
            if s.upper().startswith("END PERIOD"):
                break
            if inp and s and not s.startswith("#"):
                p = s.split()
                if len(p) >= 4:
                    topo[(_mvr_type(p[0]), _mvr_type(p[2]))] += 1
    topo_rows = [[f"{a} -> {b}", n] for (a, b), n in sorted(topo.items(), key=lambda x: -x[1])]
    if topo_rows:
        pd.DataFrame(topo_rows, columns=["route", "n_movers"]).to_csv(
            OUT / "mvr_topology.csv", index=False)
        internal = sum(n for (a, b), n in topo.items() if a == "UZF" and b == "UZF")
        total = sum(topo.values())
        print(f"   movers: {total} total; UZF->UZF internal recirculation = {internal} "
              f"({100*internal/total:.0f}%); rest = export/other routes")

    # ---- (c) figure: given vs received per package + net ----
    try:
        piv = (df.assign(signed=df["mean_rate_m3d"])
                 .groupby("package")["signed"].sum().sort_values())
        given = df[df.direction.str.startswith("given")].groupby("package")["mean_rate_m3d"].sum()
        recv  = df[df.direction.str.startswith("received")].groupby("package")["mean_rate_m3d"].sum()
        pkgs = sorted(set(given.index) | set(recv.index))
        import numpy as _np
        y = _np.arange(len(pkgs))
        fig, ax = plt.subplots(figsize=(8, 0.6 * len(pkgs) + 2))
        ax.barh(y - 0.2, [recv.get(p, 0.0) for p in pkgs], height=0.38,
                color="tab:blue", label="received (from mover, +)")
        ax.barh(y + 0.2, [given.get(p, 0.0) for p in pkgs], height=0.38,
                color="tab:red", label="given (to mover, -)")
        for i, p in enumerate(pkgs):
            net = recv.get(p, 0.0) + given.get(p, 0.0)
            ax.text(0, i + 0.42, f"net {net:+.0f}", fontsize=7, color="k")
        ax.set_yticks(y); ax.set_yticklabels(pkgs)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("mean rate (m³/d)   (MF6 sign: + into the package)")
        ax.set_title("CdL — MVR/CRR mover accounting per package\n"
                     "(gross FROM-MVR/TO-MVR double-count the cascade; net = actual exchange)")
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout(); fig.savefig(OUT / "mvr_accounting.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"   (mvr figure skipped: {e!r})")
    print(f"   wrote mvr_accounting.csv / .png (+ mvr_topology.csv)")


# =============================================================================
# 3d. DETAILED POND (LAK) WATER BALANCE — streams / GW / ET / rain, gains vs losses
# =============================================================================
def lake_water_balance():
    """Detailed pond (LAK) water balance as one dedicated figure: rainfall,
    evaporation, exchange with STREAMS (MVR in = inflow, MVR out = spill) and with
    GROUNDWATER (lakebed seepage, SPLIT into GW->pond gains and pond->GW losses),
    and storage change. The GWF connection term is split by sign so 'water from GW'
    and 'water lost to GW' are shown separately (the raw GWF term is only the net)."""
    LKB = WORKSPACE / f"{MODEL_NAME}.lak.bud"
    if not LKB.exists():
        print("   (lake water balance: no LAK budget, skipping)")
        return
    print(">> [3d] pond (LAK) detailed water balance …")
    cbf = flopy.utils.CellBudgetFile(str(LKB))
    kk = [k for k in cbf.get_kstpkper() if k[1] >= SPINUP_NPER] or cbf.get_kstpkper()
    recs = {(r.decode() if isinstance(r, bytes) else r).strip().upper()
            for r in cbf.get_unique_record_names()}

    def term(name, sign=None):
        tot = []
        for k in kk:
            try:
                d = cbf.get_data(text=name, kstpkper=k)
            except Exception:
                continue
            if not d:
                continue
            a = d[0]
            if hasattr(a, "dtype") and a.dtype.names and "q" in a.dtype.names:
                q = a["q"]
                if sign == "+":
                    q = q[q > 0]
                elif sign == "-":
                    q = q[q < 0]
                tot.append(float(np.sum(q)) if len(q) else 0.0)
        return np.mean(tot) if tot else 0.0

    bal = {}
    if "RAINFALL" in recs:    bal["Rainfall (in)"]           = term("RAINFALL")
    if "FROM-MVR" in recs:    bal["From streams+runoff (MVR, in)"] = term("FROM-MVR")
    if "GWF" in recs:
        bal["From GW — bed seepage in"] = term("GWF", "+")
        bal["To GW — bed seepage out"]  = term("GWF", "-")
    if "EVAPORATION" in recs: bal["Evaporation (out)"]       = term("EVAPORATION")
    if "TO-MVR" in recs:      bal["To streams — spill (MVR, out)"] = term("TO-MVR")
    for extra in ("RUNOFF", "EXT-INFLOW", "EXT-OUTFLOW", "WITHDRAWAL"):
        if extra in recs and abs(term(extra)) > 1e-6:
            bal[extra.title()] = term(extra)
    if "STORAGE" in recs:     bal["Storage change"]          = term("STORAGE")
    cbf.close()

    s = pd.Series(bal)
    s.to_frame("mean_rate_m3d").to_csv(OUT / "lake_water_balance.csv")
    flux = s.drop(index=[i for i in ["Storage change"] if i in s.index])
    resid = float(s.sum())    # all terms incl storage should ~ 0 (mass balance)
    order = flux.reindex(flux.abs().sort_values().index).index
    colors = ["tab:blue" if flux[i] >= 0 else "tab:red" for i in order]
    fig, ax = plt.subplots(figsize=(8.5, 0.55 * len(order) + 1.8))
    ax.barh(range(len(order)), [flux[i] for i in order], color=colors, edgecolor="k", lw=0.4)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)
    for j, i in enumerate(order):
        ax.text(flux[i], j, f" {flux[i]:+,.0f}", va="center",
                ha="left" if flux[i] >= 0 else "right", fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean rate (m³/d)   (+ into ponds / − out of ponds)")
    ax.set_title("CdL — pond (LAK) detailed water balance, all 19 ponds\n"
                 f"blue = sources, red = sinks · Σ(all incl. storage) = {resid:+.1f} m³/d (≈0 = closed)")
    fig.tight_layout(); fig.savefig(OUT / "lake_water_balance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote lake_water_balance.png / .csv  (Σ residual {resid:+.1f} m³/d)")
    print("   NOTE: 'From streams+runoff (MVR)' is the aggregate mover inflow; run the model "
          "with the new MVR budget_filerecord + mvr_flux_by_route() to split it by source (SFR vs UZF).")


# =============================================================================
# 3e. ROUTE-LEVEL MVR/CRR FLUX — needs the MVR budget file (budget_filerecord)
# =============================================================================
def mvr_flux_by_route():
    """Per-mover (provider->receiver) MVR fluxes from the MVR budget file.
    Requires the MVR package to be built with budget_filerecord=<name>.mvr.bud
    (queued 2026-08-04). Until the model is re-run with that option the file is
    absent and this is skipped. Produces route-level CRR fluxes (m3/d) so the
    aggregate FROM-MVR/TO-MVR can be split by source/target."""
    mvrb = WORKSPACE / f"{MODEL_NAME}.mvr.bud"
    if not mvrb.exists():
        print(f"   (mvr_flux_by_route: {mvrb.name} not found — re-run the model with the "
              f"MVR budget_filerecord to enable this. Skipping.)")
        return
    print(">> [3e] route-level MVR/CRR flux (from mvr.bud) …")
    cbf = flopy.utils.CellBudgetFile(str(mvrb))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    names = [(r.decode() if isinstance(r, bytes) else r).strip() for r in cbf.get_unique_record_names()]
    print(f"   mvr.bud records: {names}")
    # defensive: report the field layout of the first record so labeling can be refined
    try:
        first = cbf.get_data(text=names[0], kstpkper=kk[0])[0]
        if hasattr(first, "dtype") and first.dtype.names:
            print(f"   record fields: {first.dtype.names}")
    except Exception:
        pass
    means = {}
    for nm in names:
        tot = []
        for k in kk:
            try:
                d = cbf.get_data(text=nm, kstpkper=k)
            except Exception:
                continue
            if not d:
                continue
            a = d[0]
            if hasattr(a, "dtype") and a.dtype.names and "q" in a.dtype.names:
                tot.append(float(np.sum(a["q"])))
            else:
                tot.append(float(np.nansum(np.asarray(a))))
        if tot:
            means[nm] = np.mean(tot)
    cbf.close()
    if not means:
        print("   !! no MVR budget records summed.")
        return
    s = pd.Series(means).sort_values()
    s.to_frame("mean_rate_m3d").to_csv(OUT / "mvr_flux_by_route.csv")
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(s) + 1.5))
    ax.barh(s.index, s.values, color="tab:purple", edgecolor="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean MVR flux (m³/d)")
    ax.set_title("CdL — route-level MVR/CRR flux (per receiver package)")
    fig.tight_layout(); fig.savefig(OUT / "mvr_flux_by_route.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote mvr_flux_by_route.png / .csv  ({len(s)} records)")


def sfr_outlet_obs_vs_sim():
    """Outlet streamflow: the regionalised SYNTHETIC observation vs the MODFLOW-computed
    value.  Sim = SFR ext-outflow at the outlet reach (cdl_gwf.obs.sfr.csv, per timestep).
    Obs = sfr_inlet_series.csv Q_outlet = q_donor x (A_upstream + A_CdL); REAL donor months
    (1981-1990) are the true observations, the climatology gap is drawn faint for context.
    -> outlet_obs_vs_sim.png / .csv (+ RMSE / bias / r on the overlapping obs months)."""
    print(">> SFR outlet: synthetic-observed vs MODFLOW …")
    if not OBS_SFR.exists():
        print(f"   !! {OBS_SFR.name} not found — add the SFR obs to the model + rerun. Skipping.")
        return
    sim = pd.read_csv(OBS_SFR)
    ocands = [c for c in sim.columns if c.lower() == "outlet"]
    if not ocands:
        print("   !! no 'outlet' column in the SFR obs csv. Skipping.")
        return
    tcol, ocol = sim.columns[0], ocands[0]
    t = sim[tcol].to_numpy(dtype=float)
    spin = pd.date_range(SIM_START, periods=SPINUP_NPER + 1, freq="MS")
    spinup_days = (spin[-1] - spin[0]).days
    keep = t >= spinup_days - 1e-6
    sdates = SIM_START + pd.to_timedelta(t[keep] - spinup_days, unit="D")
    q_sim = np.abs(sim.loc[keep, ocol].to_numpy(dtype=float))          # ext-outflow magnitude, m3/d
    sim_s = pd.Series(q_sim, index=sdates).resample("MS").mean()       # month-start monthly means

    obs = (pd.read_csv(SFR_INLET_SERIES, parse_dates=["date"]).set_index("date")
           if SFR_INLET_SERIES.exists() else None)
    if obs is None:
        print(f"   (no {SFR_INLET_SERIES.name} — plotting the MODFLOW outlet only)")

    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(sim_s.index, sim_s.values, color="tab:red", lw=1.2, zorder=6,
            label="MODFLOW outlet (SFR ext-outflow)")
    stat = ""
    if obs is not None:
        _real = obs["is_real"].astype(bool)
        gap, real = obs[~_real], obs[_real]
        ax.plot(gap.index, gap["Q_outlet_m3d"], color="0.6", lw=0.6, ls=":", zorder=3,
                label="synthetic outlet — climatology gap (not observed)")
        ax.plot(real.index, real["Q_outlet_m3d"], color="navy", lw=0.7, alpha=0.45, zorder=4)
        ax.scatter(real.index, real["Q_outlet_m3d"], s=15, color="navy", zorder=7,
                   label="synthetic OBSERVED outlet (real donor months)")
        j = real[["Q_outlet_m3d"]].join(sim_s.rename("sim"), how="inner").dropna()
        if len(j):
            o, s = j["Q_outlet_m3d"].to_numpy(), j["sim"].to_numpy()
            rmse, bias = float(np.sqrt(np.mean((s - o) ** 2))), float(np.mean(s - o))
            r = float(np.corrcoef(o, s)[0, 1]) if len(j) > 2 else float("nan")
            stat = f"  |  {len(j)} obs months: RMSE {rmse:,.0f}, bias {bias:+,.0f} m³/d, r={r:.2f}"
            j.rename(columns={"Q_outlet_m3d": "obs_synthetic_m3d", "sim": "modflow_m3d"}).to_csv(
                OUT / "outlet_obs_vs_sim.csv")
    ax.set_yscale("symlog", linthresh=1000)
    ax.set_ylabel("outlet discharge (m³/d, symlog)")
    ax.set_xlim(sim_s.index.min(), sim_s.index.max())
    ax.set_title("CdL outlet — synthetic observed vs MODFLOW-computed" + stat)
    ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(OUT / "outlet_obs_vs_sim.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote outlet_obs_vs_sim.png / .csv{stat}")


def run_timing():
    """Solver effort / estimated wall-clock time per stress period, from mfsim.lst
    ('CALLS TO NUMERICAL SOLUTION' outer + 'TOTAL ITERATIONS' inner per step; time
    apportioned from the reported 'Elapsed run time' — iterations are hardware-independent,
    so which periods are expensive is robust even if wall time was throttled/paused). Overlays
    the inlet inflow to show whether the time-varying boundary drives the cost.
    -> sp_run_timing.png / .csv."""
    print(">> Run timing per stress period …")
    if not MFSIM_LST.exists():
        print(f"   !! {MFSIM_LST.name} not found. Skipping.")
        return
    txt = MFSIM_LST.read_text(errors="ignore")
    rec = re.findall(r"(\d+)\s+CALLS TO NUMERICAL SOLUTION IN TIME STEP\s+(\d+)\s+"
                     r"STRESS PERIOD\s+(\d+)\s+(\d+)\s+TOTAL ITERATIONS", txt)
    if not rec:
        print("   !! no per-step iteration records in mfsim.lst. Skipping.")
        return
    tim = (pd.DataFrame([(int(sp), int(o), int(it)) for (o, _ts, sp, it) in rec],
                        columns=["sp", "outer", "inner"])
           .groupby("sp", as_index=False).agg(outer=("outer", "sum"), inner=("inner", "sum")))
    nper = len(tim)
    m = re.search(r"Elapsed run time:\s*(?:(\d+)\s*Hours?,\s*)?(?:(\d+)\s*Minutes?,\s*)?"
                  r"([\d.]+)\s*Seconds?", txt)
    total_s = (int(m.group(1) or 0) * 3600 + int(m.group(2) or 0) * 60 + float(m.group(3))) if m else None
    if total_s:
        tim["est_s"] = total_s * tim["inner"] / tim["inner"].sum()
    # SP -> date: spin-up = first SPINUP_NPER months replay year 1, then the main series
    main_n = nper - SPINUP_NPER
    if main_n > 0:
        main_dates = pd.date_range(SIM_START, periods=main_n, freq="MS")
        all_dates = list(main_dates[:SPINUP_NPER]) + list(main_dates)
    else:
        all_dates = list(pd.date_range(SIM_START, periods=nper, freq="MS"))
    tim["date"] = all_dates[:nper]
    inflow = None
    if SFR_INLET_SERIES.exists():
        sc = pd.read_csv(SFR_INLET_SERIES, parse_dates=["date"]).set_index("date")["Q_inlet_m3d"]
        inflow = [float(sc.get(pd.Timestamp(d), np.nan)) for d in tim["date"]]
    tim.to_csv(OUT / "sp_run_timing.csv", index=False)

    ycol = "est_s" if total_s else "inner"
    ylab = "est. solve time / period (s)" if total_s else "inner iterations / period"
    ttl = ("CdL run — " + ("estimated time" if total_s else "solver iterations") + " per stress period"
           + (f"  (total {int(total_s // 3600)}h{int(total_s % 3600 // 60)}m over {nper} SPs; SP1=steady-state)"
              if total_s else f"  ({nper} SPs; SP1=steady-state)"))
    if inflow is not None:
        fig, (ax, ax2) = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True,
                                      gridspec_kw={"height_ratios": [2, 1]})
    else:
        fig, ax = plt.subplots(figsize=(13, 5)); ax2 = None
    ax.bar(tim["date"], tim[ycol], width=20, color="tab:red", alpha=0.75)
    ax.set_ylabel(ylab); ax.set_title(ttl)
    axb = ax.twinx(); axb.plot(tim["date"], tim["inner"], color="k", lw=0.5, alpha=0.4)
    axb.set_ylabel("inner iterations", fontsize=8)
    for _, r0 in tim.sort_values("inner", ascending=False).head(3).iterrows():
        ax.annotate(f"SP{int(r0.sp)}", (r0["date"], r0[ycol]), fontsize=7, ha="center", va="bottom")
    if ax2 is not None:
        ax2.plot(tim["date"], inflow, color="tab:blue", lw=0.7)
        ax2.set_ylabel("inlet inflow (m³/d)"); ax2.set_yscale("symlog", linthresh=1000)
        ax2.grid(alpha=0.3)
    (ax2 or ax).set_xlabel("date")
    fig.tight_layout(); fig.savefig(OUT / "sp_run_timing.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    top = tim.sort_values("inner", ascending=False).iloc[0]
    corr = (float(np.corrcoef(tim["inner"], np.nan_to_num(inflow))[0, 1])
            if inflow is not None and nper > 2 else float("nan"))
    print(f"   wrote sp_run_timing.png / .csv  (nper {nper}, total iters {tim['inner'].sum():,}"
          + (f", ~{total_s / 60:.0f} min" if total_s else "")
          + f"; hottest SP{int(top.sp)} = {int(top.inner)} iters"
          + (f"; corr(iters,inflow)={corr:.2f}" if inflow is not None else "") + ")")


# =============================================================================
if __name__ == "__main__":
    print(f">> Post-processing {MODEL_NAME} -> {OUT}")
    for fn in (obs_timeseries, head_maps, depth_maps, list_budget, layer_storage, lake_outputs):
        try:
            fn()
        except Exception as e:
            print(f"   !! {fn.__name__} failed: {e!r}")
    try:
        package_budget(UZFB, "unsaturated zone (UZF)", "uzf_budget_mean")
    except Exception as e:
        print(f"   !! uzf budget failed: {e!r}")
    try:
        sfr_budget()
    except Exception as e:
        print(f"   !! sfr budget failed: {e!r}")
    try:
        sfr_outlet_obs_vs_sim()
    except Exception as e:
        print(f"   !! sfr outlet obs-vs-sim failed: {e!r}")
    try:
        run_timing()
    except Exception as e:
        print(f"   !! run timing failed: {e!r}")
    try:
        mvr_accounting()
    except Exception as e:
        print(f"   !! mvr accounting failed: {e!r}")
    try:
        lake_water_balance()
    except Exception as e:
        print(f"   !! lake water balance failed: {e!r}")
    try:
        mvr_flux_by_route()
    except Exception as e:
        print(f"   !! mvr flux by route failed: {e!r}")
    print(">> Done. Outputs in:", OUT.resolve())
