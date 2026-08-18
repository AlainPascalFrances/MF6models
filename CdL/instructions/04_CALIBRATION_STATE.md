# 04 — Calibration state (as of 2026-08-17)

## Where the calibration stands

A first pestpp-ies calibration (60 realisations, `NOPTMAX=3`) **converged but over-fit**:

- phi fell 30.6M → 2.26M; heads improved strongly (P5 −83 %, P0 −73 %); the GHB<DRN inequality
  and virtual-P4 outlet constraints were satisfied; AET is a ~1.2M parameter-insensitive **floor**
  (~50 % of the final phi).
- **But the parameter field is not trustworthy:** 58 % of parameters (148/257) railed to a bound,
  the ensemble collapsed from 33 to 17 realisations (rank ≤ 16 vs 257 parameters), and the posterior
  log-K field is a **checkerboard** (negative Moran's I in all layers). Running the model with those
  parameters **diverges** (a high-K spike, K≈95–99, at the pilot points pinned to the upper bound).

**Therefore the bundled `pest_optimised.npz` is an over-fit set — kept only to exercise the import
machinery. A `USE_PEST_PARAMS=True` run with it will diverge. Re-calibrate first (below).**

## Root cause

The geostruct in `build_pst.py` was only used to *krige* the pilot points onto the grid, **not** to
build the prior *ensemble* — pestpp-ies drew its prior from a diagonal (bounds-only) covariance, so
the pilot points were spatially white (prior Moran's I ≈ 0) and IES was free to over-fit a checkerboard.

## The fix (already applied in this code — the "3-fix package")

1. **Geostatistical prior ensemble.** `build_pst.py` now draws a spatially-correlated prior with
   `pf.draw(num_reals=200)` → `template/prior_pe.jcb`, using per-layer variogram ranges
   (`GS_RANGE = {1:800, 2:1600, 3:1600}` m). Verified: prior Moran's I now **+0.2 … +0.8** (was ≈ 0).
2. **Localization.** `run_ies.py` sets `ies_parameter_ensemble = "prior_pe.jcb"` and
   `ies_autoadaloc = True` (automatic adaptive localization) so a rank-limited ensemble cannot push
   distant pilot points to their bounds.
3. **Larger ensemble.** `ies_num_reals = 150` (the prior file holds 200 to survive ~30 % attrition).

Also fixed along the way: the DRN outlet multipliers are **split by layer** (`drnwest_u1/dp`,
`drnsec_dp`) end-to-end (`build_pst.py` → `extract_pest_optimised.py` → the model loader/application);
`postproc_ies.py` head-key regex broadened for the OUTW/OUTS virtual-piezometer observations.

## To re-calibrate (recommended, on a capable machine)

150 realisations × ~1–1.5 h each is a server/workstation job, not a laptop one.

1. Rebuild the chain: `cdl_gwf_model_fable_v2.py` (base) → `pest_prep.py` → `build_pst.py`
   (writes `prior_pe.jcb`) → `run_ies.py` `NOPTMAX=0` (base + auto `make_obs.py`).
2. `run_ies.py` `NOPTMAX=3` for the IES.
3. `extract_pest_optimised.py` → fresh `pest_optimised.npz`.
4. `cdl_gwf_model_fable_v2.py` `USE_PEST_PARAMS=True` — should now complete all 557 stress periods.

### Success criteria (how to confirm the checkerboard is gone)

- Prior Moran's I **positive**, posterior **stays positive**.
- Bound-hitting **well below 58 %**.
- The `USE_PEST_PARAMS=True` forward run **converges** through the full transient.

## Observation groups & weighting (context)

257 adjustable parameters (kh1/kh2/kh3 pilot points + ghb_cond/ghb_head + drnseep + split DRN outlets),
4358 non-zero observations: piezometer heads (P0–P6, synthetic SNIRH-regionalised), MODIS AET by
land-cover zone, SFR inlet/outlet streamflow (1/Q baseflow weighting), a GHB<DRN inflow-inequality,
and two "virtual-P4" outlet water-table constraints (targeted in *depth*, at half a real piezometer's weight).
