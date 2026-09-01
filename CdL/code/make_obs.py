"""
Populate observation TARGETS and WEIGHTS in cdl.pst.
  AET    -> MODIS MOD16A2GF zonal monthly mm (mod16_aet_zonal.csv), 2000-2025 only.
  HEAD   -> SYNTHETIC P0-P6 water-table elevations, regionalised from SNIRH analogs
            (cdl_synthetic_piezo.csv, replaces the old qualitative placeholders),
            post spin-up only; confidence-graded weights (P5/P6 > P0 > P3/P4 > P1/P2).
  SFR    -> inlet ext-inflow + outlet ext-outflow (cdl_gwf.obs.sfr.csv). Target = the
            regionalised MONTHLY donor series (sfr_inlet_series.csv: inlet = q_donor x
            A_upstream, outlet = inlet + q_donor x A_CdL), weighted on REAL donor months
            (1981-1990) only, with RELATIVE (1/Q) weighting so flood peaks don't dominate and
            the fit targets baseflow (the groundwater signal). The OUTLET is what the DRN
            conductances (drnseep / drnwest / drnsec) calibrate against.  Fallback: constant
            long-term mean if no sidecar.
Weighting: RELATIVE weights within each group (head confidence; sfr inlet<outlet), then
balance the three groups to comparable phi using the base-run residuals.  Zero weight
where no data.  Writes cdl.pst in the template (backing up to cdl.interface.pst).
"""
import config
import re, shutil
from pathlib import Path
import numpy as np, pandas as pd, pyemu

PEST = Path(str(config.PEST))
TEMPLATE = PEST / "template"
BASE_OBS = PEST / "master" / "cdl.base.obs.csv"         # base (initial-parameter) simulated obs
SNIRH = PEST / "snirh_data_availability"
SYN_PIEZO = SNIRH / "cdl_synthetic_piezo.csv"           # regionalised P0-P6 water-table elevation
INLET_OUTLET = SNIRH / "inlet_outlet_series.csv"        # regionalised inlet/outlet streamflow
# SIM_START is read from the model's sidecar (single source of truth: set the date
# ONLY in cdl_gwf_model_fable_v2.py). Fallback 1981 (full period) if not yet written.
_SSF = Path((str(config.MODEL) + r"\last_sim_start.txt"))
SIM_START = pd.Timestamp(_SSF.read_text().strip()) if _SSF.exists() else pd.Timestamp(1981, 1, 1)
SPINUP_NPER = 12

# Confidence-graded RELATIVE head weights (SNIRH import methodology §3/§5): P5/P6 excellent
# (354 obs, full period) > P0 good (294) > P3/P4 moderate (109/94) > P1/P2 LOW (8 obs).
HEAD_CONF_W = {"P0": 0.7, "P1": 0.15, "P2": 0.15, "P3": 0.4, "P4": 0.4, "P5": 1.0, "P6": 1.0}
# RELATIVE streamflow weights: the outlet (DRN calibration) counts more than the inlet
# (which just pins the sfr_inflow multiplier to the regionalised estimate).
SFR_REL_W = {"inlet": 0.5, "outlet": 1.0}
# MF6 SFR-obs SIGN convention: ext-inflow is POSITIVE (into the reach), ext-outflow is
# NEGATIVE (leaving the model). The regionalised targets are positive discharges, so the
# outlet target MUST be negated to match the simulated sign — otherwise obs(+) vs sim(-)
# gives a spurious residual ~2x the flow and PEST would drive the outlet toward zero.
SFR_SIGN = {"inlet": 1.0, "outlet": -1.0}
# RELATIVE (1/Q) streamflow weighting — heteroscedastic error model sigma = REL*|Q| + FLOOR,
# weight = SFR_REL_W / sigma. Flood peaks (large Q) get large sigma -> LOW weight; the fit
# focuses on the LOW/BASEFLOW months, i.e. the GROUNDWATER-discharge signal (user 2026-08-13,
# "our focus is mainly groundwater"). FLOOR keeps dry-month (Q~0) weights finite and sets the
# baseflow scale below which months weight ~equally. Tune REL up / FLOOR up to de-emphasise
# streamflow further; the 3-group phi balance still normalises the overall sfr influence.
SFR_REL_ERR = 0.20        # relative flow error (20%)
SFR_MIN_ERR = 3000.0      # m3/d floor (~ baseflow scale; prevents 1/Q blow-up at Q->0)
# Each group is scaled to this fraction of the AET group phi at the base run (tune to
# dial a group's overall influence — e.g. raise SFR_PHI_FRAC to lean harder on the outlet).
HEAD_PHI_FRAC = 1.0
SFR_PHI_FRAC = 1.0
# OUTLET-head constraint (user 2026-08-14): the west + secondary DRN-outlet GW level is tied to
# ~P4's mean level. FIXED weight (NOT in the piezo head-group balance). BC inequality weight:
# ghb_lt_drn (mean ghb + mean drn, <0 when GHB in < DRN out) in a 'less_than' group -> penalty
# (weight*residual)^2 ONLY when violated; ~1000 m3/d violation -> phi ~(0.1*1000)^2 = 1e4.
OUTLET_HEAD_W = 0.3       # base weight on each west/secondary outlet-WT obs (same base as a nominal piezo)
VIRTUAL_WEIGHT_FRAC = 0.5 # outlets are VIRTUAL (no real piezo) -> half the real head-data weight (user 2026-08-15)
BC_CONSTRAINT_W = 0.1    # weight on the GHB<DRN inequality residual (m3/d^-1)

pst = pyemu.Pst(str(TEMPLATE / "cdl.pst"))
od = pst.observation_data
od["weight"] = 0.0                                       # start all zero; turn on where we have data

# NOTE: identify obs by OBSNME (index), not obgnme — pyEMU puts the first use-col of each
# group under the bare group name, so an obgnme filter silently drops it. obsnme is reliable:
#   AET :  oname:aet_otype:lst_usecol:<zone>_date:<YYYY-MM-DD>
#   HEAD:  oname:head_otype:lst_usecol:<pt>_l<k>_time:<days>
#   SFR :  oname:sfr_otype:lst_usecol:<inlet|outlet>_time:<days>
is_aet = od.index.str.contains("aet_otype")
is_head = od.index.str.contains("head_otype")
is_sfr = od.index.str.contains("sfr_otype")
is_outlet = is_head & od.index.str.contains(r"usecol:out[ws]_")   # west/secondary outlet-head obs
is_piezo = is_head & ~is_outlet                                   # the P0-P6 piezometers (head balance uses these)
is_bc = od.index.str.contains("bccon_otype")                      # GHB<DRN inequality (less_than group)

# ---------------- AET targets (MODIS) ----------------
mod = pd.read_csv(PEST / "mod16_aet_zonal.csv", parse_dates=["date"]).set_index("date")
n_aet = 0
for i in od.index[is_aet]:
    m = re.search(r"usecol:(.+?)_date:(.+)$", i)
    z, d = m.group(1), pd.Timestamp(m.group(2))
    if z in mod.columns and d in mod.index and np.isfinite(mod.at[d, z]):
        od.at[i, "obsval"] = float(mod.at[d, z]); od.at[i, "weight"] = 1.0; n_aet += 1

# ---------------- obs-time -> monthly DATE map (robust to ATS substeps) ----------------
# The MF6 obs 'time' is cumulative days. ATS can split a stress period into several timesteps
# (e.g. the first period here -> obs at days 1,4,9,19,31), adding EXTRA obs times, so slicing
# times[SPINUP_NPER:] by INDEX is wrong (it mislabels every month by the number of substeps).
# Instead reconstruct each stress period's END day from calendar month lengths — the 12-month
# spin-up replays year 1, then the main series from SIM_START — and match obs times to those
# period ends, mapping each to the period's START date. ATS substeps and spin-up obs don't
# match and are dropped.
head_info = {i: re.search(r"usecol:(p\d+)_l\d+_time:(.+)$", i) for i in od.index[is_piezo]}
_obstimes = np.array(sorted({float(re.search(r"_time:(.+)$", i).group(1))
                             for i in od.index[is_head | is_sfr]}))
_spin_days = float(sum(d.days_in_month for d in pd.date_range(SIM_START, periods=SPINUP_NPER, freq="MS")))
_main_dates = pd.date_range(SIM_START, periods=len(_obstimes), freq="MS")   # >= n main periods (generous)
_end_days = _spin_days + np.cumsum([d.days_in_month for d in _main_dates]).astype(float)
time2date = {}
for _ed, _dt in zip(_end_days, _main_dates):
    _j = int(np.argmin(np.abs(_obstimes - _ed)))
    if abs(_obstimes[_j] - _ed) < 0.75:            # this main period has an obs at its end-day
        time2date[_obstimes[_j]] = _dt
real_times = sorted(time2date)                       # kept for the SFR fallback spin-up cutoff

# synthetic monthly water-table ELEVATION per point (resample the irregular series to MS)
syn = pd.read_csv(SYN_PIEZO, parse_dates=["date"])
piezo = {}
for pt, g in syn.groupby("cdl"):
    s = g.set_index("date")["wl_elev_m"].sort_index()
    s.index = s.index.to_period("M").to_timestamp()
    piezo[str(pt).upper()] = s[~s.index.duplicated()]

n_head = 0
for i, m in head_info.items():
    pt = m.group(1).upper()                                      # p0 -> P0
    d = time2date.get(float(m.group(2)))
    if d is None or pt not in piezo:
        continue
    if d in piezo[pt].index and np.isfinite(piezo[pt].at[d]):
        od.at[i, "obsval"] = float(piezo[pt].at[d])
        od.at[i, "weight"] = HEAD_CONF_W.get(pt, 0.3)            # relative; group-scaled below
        n_head += 1

# ---------------- OUTLET constraint: a "virtual P4" at each outlet, in DEPTH not elevation ----------------
# (user 2026-08-15) The outlet WT must sit at the SAME DEPTH below its own ground as P4, NOT at P4's
# absolute elevation (which would force the high secondary outlet ~10 m under water). So the target
# head at each outlet = outlet_top - P4_depth(t): a virtual piezometer that mimics P4's depth-to-water-
# table behaviour, corrected to the outlet's own model top. The model observes the outlet water table
# (shallowest layer) and writes each outlet's top to outlet_cells.csv.
_p4_depth = None
if "P4" in {str(p).upper() for p in syn["cdl"]}:
    _g4 = syn[syn["cdl"].str.upper() == "P4"].set_index("date")["depth_m"].sort_index()
    _g4.index = _g4.index.to_period("M").to_timestamp()
    _p4_depth = _g4[~_g4.index.duplicated()]
_otop = {}
OUTLET_META = Path(str(config.MODEL)) / "outlet_cells.csv"
if OUTLET_META.exists():
    _om = pd.read_csv(OUTLET_META)
    _otop = {str(k).lower(): float(v) for k, v in zip(_om["name"], _om["top"])}
n_outlet = 0
for i in od.index[is_outlet]:
    _mt = re.search(r"usecol:(out[ws])_l\d+_time:(.+)$", i)
    if _mt is None or _p4_depth is None:
        continue
    _which = _mt.group(1)                                        # 'outw' or 'outs'
    d = time2date.get(float(_mt.group(2)))
    if d is None or _which not in _otop or d not in _p4_depth.index:   # only where the virtual P4 has data
        continue
    od.at[i, "obsval"] = _otop[_which] - float(_p4_depth.at[d])  # outlet_top - P4_depth = same depth as P4
    od.at[i, "weight"] = OUTLET_HEAD_W
    n_outlet += 1

# ---------------- SFR targets (inlet ext-inflow, outlet ext-outflow) ----------------
# TIME-VARYING regionalised series: the model writes sfr_inlet_series.csv (per-month inlet =
# q_donor x A_upstream, outlet = inlet + q_donor x A_CdL, is_real flag). REAL donor months
# (1981-1990) carry weight; the climatology-filled gap months do NOT (exactly like the
# piezometers, which are weighted only where the SNIRH analog has data). Falls back to the
# constant long-term MEAN if the sidecar is absent (i.e. the model ran constant-inflow mode).
SFR_SIDECAR = Path((str(config.MODEL) + r"\sfr_inlet_series.csv"))
sfr_info = {i: re.search(r"usecol:(inlet|outlet)_time:(.+)$", i) for i in od.index[is_sfr]}
n_sfr = 0
if SFR_SIDECAR.exists():
    sc = pd.read_csv(SFR_SIDECAR, parse_dates=["date"]).set_index("date")
    _col = {"inlet": "Q_inlet_m3d", "outlet": "Q_outlet_m3d"}
    for i, m in sfr_info.items():
        if m is None:
            continue
        name, d = m.group(1), time2date.get(float(m.group(2)))
        if d is None or d not in sc.index:                     # skip spin-up / gap outside the series
            continue
        q = float(sc.at[d, _col[name]])                        # positive discharge magnitude
        od.at[i, "obsval"] = SFR_SIGN[name] * q
        od.at[i, "weight"] = (SFR_REL_W[name] / (SFR_REL_ERR * abs(q) + SFR_MIN_ERR)
                              if bool(sc.at[d, "is_real"]) else 0.0)    # relative (1/Q) on real months
        n_sfr += int(od.at[i, "weight"] > 0)
    _sfr_desc = f"time-varying donor series (sidecar); {n_sfr} real-month targets weighted"
else:
    _io = pd.read_csv(INLET_OUTLET)
    SFR_TARGET = {"inlet": float(_io["Q_inlet_m3d"].mean()), "outlet": float(_io["Q_outlet_m3d"].mean())}
    _spin_cut = real_times[0] if real_times else -np.inf       # first post-spin-up model time
    for i, m in sfr_info.items():
        if m is None:
            continue
        name, t = m.group(1), float(m.group(2))
        if t < _spin_cut:                                      # skip spin-up
            continue
        q = SFR_TARGET[name]
        od.at[i, "obsval"] = SFR_SIGN[name] * q
        od.at[i, "weight"] = SFR_REL_W[name] / (SFR_REL_ERR * abs(q) + SFR_MIN_ERR); n_sfr += 1
    _sfr_desc = f"constant mean (inlet {SFR_TARGET['inlet']:.0f}, outlet {SFR_TARGET['outlet']:.0f} m3/d)"

# ---------------- BC inequality: GHB inflow < DRN-west outflow (less_than group, obsval 0) ----
n_bc = 0
for i in od.index[is_bc]:
    od.at[i, "obsval"] = 0.0                                   # derived (mean ghb + mean drn) must be < 0
    od.at[i, "weight"] = BC_CONSTRAINT_W
    n_bc += 1
od.loc[is_bc, "obgnme"] = "less_than_ghbdrn"                   # ensure the pestpp inequality group name

# ---------------- balance the three DATA groups (base residuals) — piezos, not outlet/constraint ----
if not BASE_OBS.exists():
    raise SystemExit(
        f"\n>> {BASE_OBS} not found.\n"
        f">> make_obs.py needs the NOPTMAX=0 base-run residuals to balance the group weights.\n"
        f">> Don't run make_obs.py directly: run run_ies.py with NOPTMAX=0 first — it does the\n"
        f">> base run, writes cdl.base.obs.csv, and then auto-runs make_obs.py.\n")
base = pd.read_csv(BASE_OBS, index_col=0).iloc[0]              # base simulated obs (one row)
res = (od["obsval"] - base.reindex(od.index).astype(float)).fillna(0.0)  # obs - sim
_phi = lambda mask: float((od.weight[mask] ** 2 * res[mask] ** 2).sum())
phi_aet, phi_head, phi_sfr = _phi(is_aet), _phi(is_piezo), _phi(is_sfr)
if phi_head > 0:
    _hscale = np.sqrt(phi_aet * HEAD_PHI_FRAC / phi_head)
    od.loc[is_piezo, "weight"] *= _hscale                       # scale real piezos into the head group
    od.loc[is_outlet, "weight"] *= _hscale * VIRTUAL_WEIGHT_FRAC # virtual P4s: same scaling, HALF the weight
if phi_sfr > 0:
    od.loc[is_sfr, "weight"] *= np.sqrt(phi_aet * SFR_PHI_FRAC / phi_sfr)

if not (TEMPLATE / "cdl.interface.pst").exists():             # preserve the clean (pre-weights) backup
    shutil.copy2(TEMPLATE / "cdl.pst", TEMPLATE / "cdl.interface.pst")
pst.control_data.noptmax = -1                                 # -1 = prior Monte Carlo first; set >0 for IES
pst.write(str(TEMPLATE / "cdl.pst"), version=2)

nz = lambda mask: int(((od.weight > 0) & mask).sum())
print(f">> AET targets set: {n_aet} (nz {nz(is_aet)}); HEAD: {n_head} (nz {nz(is_piezo)}); "
      f"SFR: {n_sfr} (nz {nz(is_sfr)}); OUTLET-WT: {n_outlet} @ {VIRTUAL_WEIGHT_FRAC:.0%} real-weight "
      f"(virtual-P4 depth {(_p4_depth.min() if _p4_depth is not None else float('nan')):.2f}-"
      f"{(_p4_depth.max() if _p4_depth is not None else float('nan')):.2f} m; tops "
      + ", ".join(f"{k}={v:.1f}" for k, v in sorted(_otop.items())) + "); "
      f"BC-ineq: {n_bc}")
print(f">> SFR targets: {_sfr_desc}")
print(f">> group phi (base): aet={phi_aet:.1f}, head={phi_head:.1f}, sfr={phi_sfr:.1f} "
      f"(scaled: head->{phi_aet*HEAD_PHI_FRAC:.1f}, sfr->{phi_aet*SFR_PHI_FRAC:.1f})")
print(f">> total nonzero-weight obs: {int((od.weight>0).sum())} / {pst.nobs}; "
      f"noptmax={pst.control_data.noptmax}")
