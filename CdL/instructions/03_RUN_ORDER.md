# 03 — Run order

All scripts live in `CdL/code/`. Run them from an environment where `conda activate mf6models`
is active (the author runs them in **Spyder**; reload a script after editing its flags).

## A. Forward model (no calibration)

The canonical model is **`cdl_gwf_model_fable_v2.py`**. It builds the grid (or reuses the cached
`voronoi_grid.pkl`), assembles the MF6 packages (DIS/NPF/UZF/SFR/LAK/MVR/GHB/DRN), runs MODFLOW 6,
and writes outputs to `<DATA_ROOT>/CdL_model/_output/<stamp>/`.

Key flag near the top:

- `USE_PEST_PARAMS = False` → **base / uncalibrated** run (also the mode used to build the PEST interface).
- `USE_PEST_PARAMS = True`  → apply the calibrated parameter set from `pest_optimised.npz`.
  With `AUTO_EXTRACT_PEST = True` it first re-runs `extract_pest_optimised.py` to refresh the npz
  from the current `master/` dir.

```bash
python cdl_gwf_model_fable_v2.py     # writes MF6 inputs, runs mf6, quick post-process
python postprocess_cdl.py            # figures: head/depth time series, budgets, maps, SFR/LAK
```

On divergence, the model auto-runs `diag_divergence.py`, writing
`divergence_map.png` + `divergence_report.csv` to the run's `_output/<stamp>/` folder.

## B. PEST++ calibration chain

Run in this exact order (the model must run once in base mode first to emit the obs files):

| Step | Script | What it does |
|---|---|---|
| 1 | `cdl_gwf_model_fable_v2.py` (`USE_PEST_PARAMS=False`) | Build + run the base model → obs csv files. |
| 2 | `pest_prep.py` | Build the `org/` dir with external MF6 arrays/lists. |
| 3 | `build_pst.py` | pyEMU PstFrom → `template/cdl.pst` **and draws the geostatistical prior `prior_pe.jcb`**. |
| 4 | `run_ies.py` with `NOPTMAX = 0` | Single base run → auto-runs `make_obs.py` (sets targets + weights from base residuals). |
| 5 | `run_ies.py` with `NOPTMAX = -1` | *(optional)* prior Monte-Carlo — checks the ensemble machinery. |
| 6 | `run_ies.py` with `NOPTMAX = 3` | The IES calibration (uses `prior_pe.jcb` + automatic localization). |
| 7 | `postproc_ies.py` | phi progress, phi-by-group, parameter & head-fit ensembles. |
| 8 | `extract_pest_optimised.py` | Writes `pest_optimised.npz` (min-phi realisation of the final iteration). |
| 9 | `cdl_gwf_model_fable_v2.py` (`USE_PEST_PARAMS=True`) | Run the model with the calibrated parameters (plausibility gate: it should complete all 557 stress periods). |

### Notes

- **`make_obs.py` is not run standalone** — it needs `master/cdl.base.obs.csv`, which the `NOPTMAX=0`
  base run produces; `run_ies.py` auto-runs it after that base run.
- `run_ies.py` header has a **SERVER-MIGRATION CHECKLIST** — set `NUM_WORKERS` to your core count and
  the executable/base paths before the big 150-realisation IES.
- `ies_num_reals = 150` in `run_ies.py`; the prior ensemble file holds 200 (buffer for ~30 % run
  attrition). See `04_CALIBRATION_STATE.md`.

## C. One-time preprocessing (only if regenerating the bundled caches)

You do **not** need these if you use the bundled `input_data` — they are how those caches were built:

| Produces | Script | Needs |
|---|---|---|
| `voronoi_layers.npz` | `build_layers.py` / `make_k_layers.py` | DEM, soils rasters, geology (GPKG) |
| `cos_zones.npz` | `make_cos_zones.py` | the COS raster |
| `cdl_synthetic_piezo.csv` | `snirh_fetch_piezo_synthetic.py` → `correct_synthetic_piezo.py` | SNIRH; GRB (model top) |
| `streamflow_21F_01H.csv` | `snirh_fetch_streamflow.py` | SNIRH donor gauge |
| MODIS AET zonal targets | `download_mod16.py` → `make_mod16_zonal.py` | MOD16A2GF |
