"""
Build the PEST++ interface (pyEMU PstFrom) for the CdL model.

PARAMETERS
  - Kh pilot points x 3 layers (multiplier, geostatistical prior)   [Kv=0.1*Kh in forward_run]
  - GHB conductance   (constant multiplier)
  - GHB external head (constant additive offset, m)
  - DRN-OUT / DRN-SEEP / DRN-WEST conductance (constant multipliers)
OBSERVATIONS
  - P0-P6 heads, monthly (cdl_gwf.obs.head.csv)
  - actual ET aggregated to land-cover zones, monthly (model_aet_zonal.csv)
  - SFR inlet ext-inflow + outlet ext-outflow, monthly (cdl_gwf.obs.sfr.csv)

Runs on the external 'org' model built by pest_prep.py; writes the template dir.
NOTE: all work is under `if __name__ == "__main__"` — pyEMU's kriging spawns worker
processes on Windows, which re-import this module; unguarded top-level code would
re-run the whole build in every child and deadlock.
"""
import shutil, pickle, multiprocessing
from pathlib import Path
import numpy as np, pandas as pd, pyemu
from scipy.spatial import cKDTree
from flopy.discretization import VertexGrid

PEST = Path(r"E:\00code_ws\DRYAD\CdL_pest")
ORG, TEMPLATE = PEST / "org", PEST / "template"
CODE = Path(r"E:\00code\flopy\dryad_cdl")


def _flat(p, dt=float):                              # wrapped MF6 array -> flat 1-D
    return np.array(Path(p).read_text().split(), dtype=dt)


def main():
    gp = pickle.load(open(ORG / "voronoi_grid.pkl", "rb"))[0]
    vgrid = VertexGrid(**gp, nlay=1)
    # BINARY active mask (1 active / 0 else). idomain also has -1 (pass-through); pyEMU
    # treats any distinct value as a pilot-point zone, so -1 must NOT leak in.
    idomain = [(_flat(ORG / f"cdl_gwf.disv_idomain_layer{L}.txt", int) == 1).astype(int)
               for L in (1, 2, 3)]

    pf = pyemu.utils.PstFrom(original_d=str(ORG), new_d=str(TEMPLATE), remove_existing=True,
                             longnames=True, zero_based=False, spatial_reference=vgrid)

    # pyEMU reads arrays with np.loadtxt, which needs a RECTANGULAR file. MF6 writes the
    # DISV K arrays WRAPPED at 20 values/line, and ncpl=6609 is not a multiple of 20, so the
    # last line has 9 cols -> "number of columns changed from 20 to 9". Rewrite the K arrays
    # ONE VALUE PER LINE (cell order preserved) in the template so PstFrom reads them; MF6 and
    # forward_run.py read free-format, so the model side is unaffected. (_flat() splits on all
    # whitespace, so it reads the wrapped file regardless of the ragged last line.)
    for L in (1, 2, 3):
        _kf = TEMPLATE / f"cdl_gwf.npf_k_layer{L}.txt"
        np.savetxt(_kf, _flat(_kf).reshape(-1, 1), fmt="%20.10E")

    # Geostatistical structure PER LAYER — correlation range ~4x the pilot-point spacing so the PRIOR
    # is spatially SMOOTH (neighbouring pilot points co-vary). This variogram now ALSO builds the prior
    # COVARIANCE the ensemble is drawn from (pf.build_prior / pf.draw after build_pst) — the piece that
    # was missing before, which let the 3-iter IES overfit to a checkerboard K field (prior Moran's I
    # ~0, posterior negative). Deeper layers (coarser 400 m spacing) get a longer range. (2026-08-16)
    GS_RANGE = {1: 800.0, 2: 1600.0, 3: 1600.0}          # L1 150 m spacing; L2/L3 400 m spacing
    def _make_gs(rng):
        return pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(contribution=1.0, a=rng),
                                        name=f"k_gs_{int(rng)}")

    # ---- explicit pilot-point grid per layer (deterministic count) ----
    xc = np.array(vgrid.xcellcenters).ravel(); yc = np.array(vgrid.ycellcenters).ravel()
    tree = cKDTree(np.c_[xc, yc])

    def pp_df(idom, spacing, base):
        xs = np.arange(xc.min() + spacing / 2, xc.max(), spacing)
        ys = np.arange(yc.min() + spacing / 2, yc.max(), spacing)
        rows, n = [], 0
        for x in xs:
            for y in ys:
                d, ic = tree.query([x, y])
                if d < spacing and idom[ic] == 1:            # keep pp only over active cells
                    rows.append([f"{base}_{n:03d}", float(x), float(y), 1, 1.0]); n += 1
        return pd.DataFrame(rows, columns=["name", "x", "y", "zone", "parval1"])

    # PARAMETER BOUNDS — widened then MODERATED 2026-08-13 (user): the obs are now genuinely
    # constraining (real SNIRH-regionalised piezo heads + relative-weighted streamflow), so let
    # the DATA limit the parameters, not artificially tight bounds (which risk params railing at
    # the limit and reporting a false-low phi at the edge of an over-restricted space). BUT the
    # first widening (kh 0.1-10, conds 0.01-100) drove 60% ensemble non-convergence — and when
    # 2/3 of draws fail, the survivors are BIASED to the tame region, so ultra-wide bounds buy
    # failures, not exploration. SWEET SPOT = wide enough to explore, tame enough that most reals
    # converge (<~30% fail): kh 0.2-5 (x25, was 0.5-1.7), conds 0.05-20. ghb_head additive (±5 m).
    # sfr_inflow is no longer a parameter (inlet = fixed observed forcing).
    PP_SPACING = {1: 150.0, 2: 400.0, 3: 400.0}          # L1 (U1) narrow -> finer
    for L in (1, 2, 3):
        ppdf = pp_df(idomain[L - 1], PP_SPACING[L], f"kh{L}")
        print(f"   layer {L}: {len(ppdf)} pilot points (spacing {PP_SPACING[L]:.0f} m)")
        pf.add_parameters(f"cdl_gwf.npf_k_layer{L}.txt", par_type="pilotpoints",
                          par_name_base=f"kh{L}", pargp=f"kh{L}",
                          zone_array=idomain[L - 1], geostruct=_make_gs(GS_RANGE[L]),
                          pp_options={"pp_space": ppdf, "use_pp_zones": True, "num_threads": 1},
                          lower_bound=0.2, upper_bound=5.0, par_style="multiplier")    # 0.2-5 (x25): MODERATED 2026-08-13 (0.1-10 -> 60% ensemble failure)

    # ---- BC parameters (one constant per package/column) ----
    pf.add_parameters("cdl_gwf.ghb_stress_period_data_1.txt", par_type="constant",
                      index_cols=[0, 1], use_cols=[3], par_name_base="ghb_cond", pargp="ghb_cond",
                      lower_bound=0.05, upper_bound=20.0, par_style="multiplier")      # 0.05-20: MODERATED 2026-08-13 (0.01-100)
    pf.add_parameters("cdl_gwf.ghb_stress_period_data_1.txt", par_type="constant",
                      index_cols=[0, 1], use_cols=[2], par_name_base="ghb_head", pargp="ghb_head",
                      lower_bound=-5.0, upper_bound=5.0, par_style="add")               # was -3..+3 m
    # DRN conductance multipliers. 2026-08-05: DRN-OUT retired (secondary outlet is now
    # the explicit line-based DRN-SECONDARY), so the externalised DRN files shift by one:
    #   cdl_gwf.drn = DRN-SEEP, cdl_gwf_0.drn = DRN-WEST, cdl_gwf_1.drn = DRN-SECONDARY.
    # flopy names DRN files by BUILD ORDER — a re-order would silently point drnseep_cond
    # at a band drain and ruin the outlet calibration. ASSERT the mapping every build:
    #   DRN-SEEP  = surface seepage on ~every active cell   -> thousands of records
    #   DRN-WEST  = ~26-cell band at the WEST outlet         -> tens of records, low x
    #   DRN-SECONDARY = ~26-cell band at the 2nd outlet      -> tens of records, higher x
    drn_files = {"drnseep_cond": "cdl_gwf.drn_stress_period_data_1.txt",
                 "drnwest_cond": "cdl_gwf_0.drn_stress_period_data_1.txt",
                 "drnsec_cond": "cdl_gwf_1.drn_stress_period_data_1.txt"}

    def _drn_nodes(fn):                                   # -> 0-based node ids (col 1 = node, 1-based)
        a = np.loadtxt(ORG / fn, ndmin=2, usecols=(1,))   # ONLY the node col — DRN-WEST now has a string boundname col
        return a[:, 0].astype(int) - 1

    seep, west, sec = (_drn_nodes(drn_files[t]) for t in ("drnseep_cond", "drnwest_cond", "drnsec_cond"))
    assert len(seep) > 500 and len(seep) > 10 * max(len(west), len(sec)), (
        f"DRN mapping BROKEN: {drn_files['drnseep_cond']} has {len(seep)} records but should be the "
        f"whole-domain DRN-SEEP (thousands); west={len(west)}, sec={len(sec)}. Check the DRN build order.")
    xw, xs = float(xc[west].mean()), float(xc[sec].mean())
    assert xw < xs, (f"DRN mapping BROKEN: cdl_gwf_0.drn mean-x {xw:.0f} is not WEST of cdl_gwf_1.drn "
                     f"{xs:.0f} — DRN-WEST / DRN-SECONDARY are swapped. Check the DRN build order.")
    print(f"   DRN mapping OK: seep {len(seep)} cells >> west {len(west)} (x̄={xw:.0f}) / "
          f"sec {len(sec)} (x̄={xs:.0f})")

    # DRN-SEEP (whole-domain surface seepage / baseflow) — one multiplier, [0.05, 20].
    pf.add_parameters(drn_files["drnseep_cond"], par_type="constant", index_cols=[0, 1], use_cols=[3],
                      par_name_base="drnseep_cond", pargp="drnseep_cond",
                      lower_bound=0.05, upper_bound=20.0, par_style="multiplier")
    # DRN-WEST / DRN-SECONDARY outlet drains — SPLIT BY LAYER (user 2026-08-15) so U1 and the deep
    # units tune INDEPENDENTLY: U1 (alluvium, base cond 1) x[0.05,20] -> cond up to 20; deep U2/U3
    # (base cond 0.1) x[0.05,10] -> cond up to 1. use_rows selects the layer rows (col 0 = 1-based).
    for _tag, _fn in (("drnwest", drn_files["drnwest_cond"]), ("drnsec", drn_files["drnsec_cond"])):
        _lay = np.loadtxt(ORG / _fn, ndmin=2, usecols=(0,))[:, 0].astype(int)
        _u1 = [int(r) for r in np.where(_lay == 1)[0]]
        _dp = [int(r) for r in np.where(_lay > 1)[0]]
        if _u1:
            pf.add_parameters(_fn, par_type="constant", index_cols=[0, 1], use_cols=[3], use_rows=_u1,
                              par_name_base=f"{_tag}_u1", pargp=f"{_tag}_u1",
                              lower_bound=0.05, upper_bound=20.0, par_style="multiplier")   # cond up to 20
        if _dp:
            pf.add_parameters(_fn, par_type="constant", index_cols=[0, 1], use_cols=[3], use_rows=_dp,
                              par_name_base=f"{_tag}_dp", pargp=f"{_tag}_dp",
                              lower_bound=0.05, upper_bound=10.0, par_style="multiplier")   # base 0.1 -> cond up to 1

    # ---- SFR external inflow at the inlet reach (293) — NOT CALIBRATED (user 2026-08-13) ----
    # The inlet is now treated as an OBSERVED / KNOWN forcing — the regionalised monthly series
    # (donor gauge 21F/01H x A_upstream), exactly like the piezometer data: a state variable we
    # OBSERVE, not a parameter we optimise. The former `sfr_inflow` multiplier is therefore
    # REMOVED, so the inlet perioddata stay at their physical (mult=1) values in every forward
    # run. The inlet is still carried as an OBSERVATION (sfr-inlet, ext-inflow) in make_obs for
    # tracking/plotting; being a prescribed input its residual is ~0, so the OUTLET is the
    # streamflow that actually informs the parameters. (To re-enable inlet calibration, restore a
    # pf.add_parameters over `cdl_gwf.sfr_perioddata_*.txt` with par_name_base="sfr_inflow".)

    # ---- observations ----
    hcols = list(pd.read_csv(ORG / "cdl_gwf.obs.head.csv", nrows=0).columns)[1:]
    pf.add_observations("cdl_gwf.obs.head.csv", index_cols="time", use_cols=hcols,
                        prefix="head", obsgp="head")
    zcols = list(pd.read_csv(ORG / "model_aet_zonal.csv", nrows=0).columns)[1:]
    pf.add_observations("model_aet_zonal.csv", index_cols="date", use_cols=zcols,
                        prefix="aet", obsgp="aet")
    # SFR streamflow: inlet (ext-inflow, reach 293) + outlet (ext-outflow) — MF6 continuous
    # obs csv (first col 'time', then the obs names). Targets/weights set in make_obs.py;
    # the outlet is what the DRN conductances calibrate against.
    scols = list(pd.read_csv(ORG / "cdl_gwf.obs.sfr.csv", nrows=0).columns)[1:]
    pf.add_observations("cdl_gwf.obs.sfr.csv", index_cols="time", use_cols=scols,
                        prefix="sfr", obsgp="sfr")
    # BC INEQUALITY: GHB inflow < DRN-west outflow (user 2026-08-14). forward_run writes the
    # derived value (mean ghb + mean drn) to bc_constraint.csv; it is < 0 when satisfied. A pestpp
    # group whose name starts with 'less_than' penalises ONLY when the simulated value exceeds its
    # obsval (set to 0 in make_obs) -> a hard one-sided constraint that never distorts when honoured.
    if (ORG / "bc_constraint.csv").exists():
        pf.add_observations("bc_constraint.csv", index_cols="cname", use_cols=["value"],
                            prefix="bccon", obsgp="less_than_ghbdrn")   # NOT 'name' (pyEMU alias of obsnme)

    pst = pf.build_pst(str(TEMPLATE / "cdl.pst"))

    # ---- GEOSTATISTICAL PRIOR COVARIANCE + ENSEMBLE (2026-08-16 anti-checkerboard fix) ----
    # Build the prior covariance from the pilot-point geostructs and DRAW a spatially-correlated
    # prior parameter ensemble. pestpp-ies otherwise draws its prior from a DIAGONAL (bounds-only)
    # covariance -> spatially WHITE pilot points -> the checkerboard/overfit K field diagnosed in the
    # 3-iter IES (prior Moran's I ~0, posterior NEGATIVE). Supplying prior_pe.jcb (via
    # ies_parameter_ensemble in run_ies.py) makes neighbouring pilot points co-vary. Draw EXTRA reals:
    # ~30% fail on MF6 non-convergence, so survivors must still number >~100 for a rank-healthy update.
    # BC constants (ghb/drn) have no geostruct -> drawn from their bounds (uncorrelated), as intended.
    N_PRIOR_REALS = 200
    pe = pf.draw(num_reals=N_PRIOR_REALS, use_specsim=False)   # draws from the geostructs -> correlated
    pe.enforce()                                               # clip any draw outside the bounds
    pe.to_binary(str(TEMPLATE / "prior_pe.jcb"))
    print(f"   geostatistical prior drawn: {N_PRIOR_REALS}-real prior_pe.jcb "
          f"(supply via ies_parameter_ensemble)")

    # ---- FIRST-ORDER TIKHONOV (preferred-difference / smoothness) on the Kh pilot points ----
    # The geostatistical PRIOR (prior_pe.jcb) makes neighbouring pilot points co-vary, but the IES
    # UPDATE still roughens the field (an un-regularised 3-iteration run drove ~59% of pilot points
    # to their bounds -> posterior negative Moran's I / checkerboard). Add explicit Tikhonov
    # regularisation - prior-information equations penalising the LOG-K DIFFERENCE between spatially-
    # correlated pilot points (preferred difference = 0 -> smooth field), built from the geostruct
    # prior COVARIANCE (also saved as prior_cov.jcb). These enter the IES COMPOSITE phi only when
    # ies_reg_factor>0 (set in run_ies.py) - that factor is the smoothness STRENGTH knob. abs_drop_tol
    # keeps only pilot-point pairs whose prior correlation exceeds it (true neighbours), so the
    # equation count stays modest - important: pestpp-ies 5.2.27 heap-crashes when reg is active with
    # a large ensemble x many PI eqs, so ~1-1.5k equations (tol 0.60) keeps the 150-real load safe.
    TIK_DROP_TOL = 0.60
    cov = pf.build_prior(fmt="binary", filename=str(TEMPLATE / "prior_cov.jcb"))
    pyemu.helpers.first_order_pearson_tikhonov(pst, cov, reset=True, abs_drop_tol=TIK_DROP_TOL)
    print(f"   Tikhonov: {pst.prior_information.shape[0]} first-order (smoothness) PI equations "
          f"added; prior_cov.jcb saved. Activate with ies_reg_factor>0 in run_ies.py.")

    # ---- forward run: our forward_run.py (applies mults -> Kv -> mf6 -> AET) ----
    # ⚠ MUST copy AFTER pf.build_pst(): build_pst() writes its OWN apply-ONLY
    # forward_run.py stub (no mf6, no AET). Copying ours afterward overwrites that
    # stub. (Copying before -> build_pst clobbers ours -> the base run does apply
    # only, mf6 never runs, obs files missing -> "model_aet_zonal.csv not found".)
    for f in ("forward_run.py", "model_aet_zonal.py", "obs_collapse.py"):
        shutil.copy2(CODE / f, TEMPLATE / f)
    shutil.copy2(PEST / "cos_zones.npz", TEMPLATE / "cos_zones.npz")

    # full path to the flopy-env python — PEST++ spawns the model command from the
    # system shell, where bare "python" would not resolve to the conda env.
    pst.model_command = r"C:\miniconda3\envs\flopy\python.exe forward_run.py"
    pst.control_data.noptmax = 0                          # 0 = one run to check the interface
    pst.write(str(TEMPLATE / "cdl.pst"), version=2)

    print("\n== PST summary ==")
    print("  parameters:", pst.npar, " groups:", list(pst.par_groups))
    print("  observations:", pst.nobs, " groups:", list(pst.obs_groups))
    print("  template:", TEMPLATE)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
