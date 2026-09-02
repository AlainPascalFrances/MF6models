========================================================================
 CdL groundwater-flow model  —  code/  (script index)
========================================================================
MODFLOW 6 / FloPy model of the Casa de Lobos (CdL) catchment, DRYAD / NbS
project.  Voronoi DISV grid (ncpl = 6609, 3 layers), monthly transient
1981-2026, calibrated with PEST++ (pestpp-ies) via pyEMU.

 * ALL machine-specific paths live in config.py — edit that ONE file (or
   set the CDL_* environment variables) when moving to a new machine.
   Every other script imports its paths from it.
 * Most heavy scripts (model run, figures, geopandas/rasterio) are meant
   to be run inside Spyder / the conda env, not head-less.
 * Typical order: config -> data prep -> grid/layers -> main model ->
   post-processing;  calibration chain is its own track (section 5).
 * See ../instructions/ for the setup + run-order cookbook.

------------------------------------------------------------------------
 0. CONFIGURATION
------------------------------------------------------------------------
config.py
    Single place for the machine-specific disk paths (BASE / MODFLOW_DIR /
    PYTHON_EXE and the paths derived from them). Imported by every script.
    Replaces the old per-script find-and-replace path workflow.

------------------------------------------------------------------------
 1. DATA PREPARATION — observations & climatic forcing
    (SNIRH scripts must run where the SNIRH API is reachable, e.g. Spyder)
------------------------------------------------------------------------
snirh_fetch_piezo_synthetic.py
    Fetch SNIRH piezometer series for the regional analog wells and build
    the SYNTHETIC piezometer targets for the catchment (head calibration).
snirh_fetch_streamflow.py
    Fetch the SNIRH streamflow gauge series (donor 21F/01H) used to
    regionalise the SFR inlet inflow and the outlet streamflow target.
snirh_task1_map.py
    Map SNIRH piezometers: reproject to EPSG:3763, filter to the T3/T7
    aquifer within the search buffer (data-availability figure).
snirh_task2_geomorph.py
    SNIRH task-2 geomorphological screening of the piezometer sets.
correct_synthetic_piezo.py
    Post-correct the raw synthetic piezometer series (two fixes) and
    overwrite cdl_synthetic_piezo.csv (consumed by make_obs / postproc).
download_mod16.py
    Download MODIS MOD16A2GF actual-ET granules (gap-filled 8-day, 500 m)
    over the CdL watershed — the RS source for the AET calibration group.
make_mod16_zonal.py
    Aggregate MOD16A2GF actual-ET to land-cover ZONES, monthly (mm) —
    the OBSERVED side of the PEST 'aet' observation group.
model_aet_zonal.py
    Aggregate the model's UZF actual-ET to the same land-cover zones,
    monthly (mm) — the SIMULATED side of the 'aet' group.
make_cos_zones.py
    Build the land-cover zone map (COS-2025 majority class per Voronoi
    cell) used by both AET zonal scripts and the calibration.

------------------------------------------------------------------------
 2. GRID, LAYERS & AQUIFER PROPERTIES (model-input preprocessing)
------------------------------------------------------------------------
preview_transition_grid.py
    Preview the smooth-transition Triangle+Voronoi grid (Daoud 2022 style)
    in a scratch dir WITHOUT touching the model workspace or grid cache.
prototype_layers_per_unit.py
    Prototype of the chosen layering (host + U1 + U2 + U3) — point-in-
    polygon K assignment per unit; a design sandbox for build_layers.
build_layers.py
    Build the top/bottom ELEVATION surfaces + thickness of the 3 hydro-
    stratigraphic units from the DTM + the hydrostrat map -> voronoi_layers.npz.
make_k_layers.py
    Produce the per-layer hydraulic-conductivity (NPF) maps from
    voronoi_layers.npz, exactly as the model assembles them.
make_layer_diagnostics.py
    Layer diagnostics (Kh per layer on a shared ramp, top/bottom surfaces,
    thickness, idomain) to sanity-check the built layers before running.

------------------------------------------------------------------------
 3. MAIN MODEL  (build + run)
------------------------------------------------------------------------
cdl_gwf_model_fable_v2.py
    THE main script — builds the full MODFLOW 6 GWF model (DISV grid,
    NPF/STO/IC, UZF, SFR, LAK, MVR, GHB, DRN, CRR), writes the inputs, runs
    mf6, and writes the pre-processing figures. USE_PEST_PARAMS loads the
    calibrated parameter set. Canonical model builder (v2 "fable").

------------------------------------------------------------------------
 4. MODEL POST-PROCESSING
------------------------------------------------------------------------
postprocess_cdl.py
    Post-process a model run (no re-run): water-balance graphs and stream-
    flow hydrographs in mm/yr, head & AET fits, lake stage/volume, MVR
    accounting, per-pond panels -> _output/<stamp>/.

------------------------------------------------------------------------
 5. PEST++ CALIBRATION  (pestpp-ies via pyEMU)
------------------------------------------------------------------------
pest_prep.py
    PEST prep step 1 — build the 'org' model dir with EXTERNAL input
    arrays/lists so pyEMU PstFrom can template them; collapse the reference
    obs to period-end rows (ATS-safe).
build_pst.py
    Build the PEST++ interface (PstFrom): pilot-point K + BC parameters,
    the geostatistical prior, the observations, and the first-order
    Tikhonov (smoothness) prior-information equations -> cdl.pst.
obs_collapse.py
    Helper — make the MF6 time-series OBS csv output deterministic under
    ATS (collapse sub-steps to stress-period-end rows). Used by pest_prep
    and forward_run so the positional .ins files always align.
make_obs.py
    Populate the observation TARGETS and WEIGHTS in cdl.pst (piezo heads,
    MODIS AET, SFR streamflow, GHB<DRN inequality, virtual-P4 outlet depth).
forward_run.py
    The per-worker forward run: apply the parameter multipliers to the MF6
    inputs, run mf6, collapse the obs, extract the model outputs for PEST.
run_ies.py
    Launch pestpp-ies with local PANTHER workers (multi-hour); sets the
    Tikhonov REG_FACTOR and worker count.
postproc_ies.py
    Post-process the pestpp-ies master outputs (no re-run): phi progress,
    parameter/obs ensembles, and the zonal AET fit -> ies_output/.
extract_pest_optimised.py
    Extract the optimised (min-phi) parameter realisation into pest_optimised.npz
    so the main model can run with the calibrated parameters.

------------------------------------------------------------------------
 6. VALIDATION & DIAGNOSTICS
------------------------------------------------------------------------
plot_initial_conditions.py
    Map the transient run's initial conditions = the no-lake steady-state
    spin-up heads (water-table elevation, depth to water, saturation).
diag_wt_vs_obs.py
    Tabulate the no-lake steady-state water table vs the observed piezo-
    meters at each obs point (model vs observed head & depth-to-water).
ss_obs_vs_computed.py
    Scatter of OBSERVED vs COMPUTED steady-state heads at the 7 piezometers
    (1:1 line + fit statistics).
pond_recharge_analysis.py
    Offline pond-recharge analysis (no-lake deliverable): focused recharge
    from each pond's contributing catchment.
diag_divergence.py
    Forensics for an MF6 convergence/divergence failure — parses the
    listing for the GW cells and advanced-package features that fail.
diag_mvr_arrows.py
    Visualise every MVR mover as a provider->receiver arrow, read from the
    WRITTEN MF6 packages (mvr/uzf/sfr/lak/drn).
diag_lake_budget.py
    Dump the LAK budget terms for the last converged period (evaporation
    vs leakage vs outflow per lake).
diag_lake_stage.py
    Dump the LAK stage trajectory per lake per stress period (lake water-
    balance diagnosis).
diag_sfr828.py
    Trace the source of a large SFR throughflow (lake budget + reaches with
    the largest MVR/JA inflows).
diag_cell2422.py
    Pin a specific SFR runaway: what cell 2422 / reach 828 are and which
    MVR providers feed that reach.
diag_sfr_ponds.py
    Check that on-channel pond cells get SFR reaches where LAK is inactive
    (per-pond comparison on the current grid).
diag_sfr_ponds_zoom.py
    Zoom panels of ALL ponds (grid cells, SFR cells, LAK cells, footprint,
    stream), read stale-proof from the written model.
diag_pond6.py
    Inspect why pond 6 fails at the first time step (cell layer stack +
    lake stage vs GW initial head).
diag_pond_perched.py
    Quantify how perched each charca is: pond bottom vs the steady-state
    water table.
diag_pond_catchments.py
    Derive each pond's contributing catchment area from the DEM (for the
    focused deep-U3 recharge fallback).
diag_buffer_vs_sfr.py
    Clarify the stream-refinement buffer corridor vs the actual SFR channel
    cells (refinement QC).
diag_u1.py
    U1-unit geometry check: dive coordinate measured from the U1 polygon
    edge (U1 overlaps U2).
diag_plan.py
    Plan view of U1/U2/U3 + transect + distance fields, to see how the
    transect crosses U1.

------------------------------------------------------------------------
 7. FIGURES — conceptual model, geology & cross-sections
------------------------------------------------------------------------
build_profile.py
    Step 1/2 of the conceptual cross-section: build the topographic transect
    (SE divide -> ERT P1..P4 -> outlet).
cross_section.py
    Final topographic cross-section with ERT P1-P6 and geological-unit
    boundaries (de-overlapped labels).
geo_section.py
    Build the geological & hydrostratigraphic cross-sections from the
    digitized section shapefile.
hydrostrat_plan.py
    Plan-view map of the 3 hydrostratigraphic units draped over the LiDAR
    hillshade (ERT, streams, ponds overlays).
plot_dem_hillshade.py
    DEM figure: hypsometric colours + hillshade relief, with DEM envelope /
    raster bbox / watershed overlays.
merge_fig1_fig5_map.py
    Combined SNIRH map (geological setting + stations + upstream drainage
    area), geology streamed live from the LNEG WMS.
export_geology_qml.py
    Export a QGIS categorized .qml style for the geology layer (Codigo),
    using the official 100k geological-map symbology.
========================================================================
