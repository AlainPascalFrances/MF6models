# 00 — CdL model overview

## Purpose

A MODFLOW 6 groundwater-flow model of the **Casa de Lobos (CdL) catchment** (~20 km², central
Portugal), built for the DRYAD / Nature-based-Solutions project to quantify the water balance and
the effect of recharge/retention ponds. It is calibrated with PEST++ (pestpp-ies) against
piezometer heads, actual evapotranspiration and streamflow.

- **CRS:** EPSG:3763 (PT-TM06).
- **Period:** monthly transient, 1981–2026 (557 stress periods incl. a 12-month spin-up).
- **Canonical script:** `code/cdl_gwf_model_fable_v2.py` (builds everything and runs MF6).

## Discretization

- **DISV** unstructured Voronoi grid, `ncpl ≈ 6609` cells, **3 layers**:
  - **L1 (U1)** — alluvium / weathered top (higher K),
  - **L2 (U2)** and **L3 (U3)** — deeper, lower-K units.
- The grid is generated with Triangle (FloPy wrapper) and cached as `voronoi_grid.pkl`; conceptual
  layer tops/bottoms come from `voronoi_layers.npz` (built from the DEM, soils and geology).

## Processes & packages

| Package | Role |
|---|---|
| **NPF** | Horizontal K per layer via pilot-point multipliers (calibrated); Kv = 0.1·Kh. |
| **UZF** | Unsaturated-zone recharge from monthly P − ET0; soil hydraulic props from the soils rasters; per-cell ET extinction (rooting) depth from the COS land-cover map. |
| **SFR** | Stream network (DEM-routed); inlet = a regionalised donor-gauge series (SNIRH 21F/01H × upstream area), outlet informs the calibration. |
| **LAK** | Recharge/retention ponds (embedded connection), coupled to SFR via **MVR**. |
| **GHB** | Eastern alluvial inflow boundary (calibrated conductance + head). |
| **DRN** | Western + secondary catchment outlets (alluvial underflow exit; per-layer conductance) and distributed surface seepage/baseflow. |

## Boundary-condition concept

Eastern **GHB** inflow (alluvial underflow in) balances a western **DRN** outflow (alluvial
underflow out), constrained so GHB-inflow < DRN-outflow. Outlet water tables are held near a
"virtual P4" **depth** (not elevation), and the drain datum sits at each layer's base so the full
saturated column can drain.

## Calibration (PEST++ / pyEMU)

- **Parameters (257):** kh1/kh2/kh3 pilot points, GHB conductance & head, DRN conductances
  (seepage + split west/secondary outlet U1/deep).
- **Observations (4358 non-zero):** synthetic SNIRH-regionalised piezometer heads (P0–P6),
  MODIS MOD16 actual-ET by land-cover zone, SFR inlet/outlet streamflow (1/Q baseflow weighting),
  a GHB<DRN inequality, and two virtual-P4 outlet water-table constraints.
- **Method:** pestpp-ies with a geostatistical prior ensemble and automatic localization.
  **See `04_CALIBRATION_STATE.md` for the current status and the important caveat about the
  bundled (over-fit) parameter set.**

## Outputs

Heads/budgets go to `<DATA_ROOT>/CdL_model/`; post-processing figures (head & depth time series,
water balance, SFR/LAK, maps) to `_output/<stamp>/` via `postprocess_cdl.py` (forward model) and
`postproc_ies.py` (calibration). On non-convergence the model auto-writes
`divergence_map.png` + `divergence_report.csv` to the same `_output/<stamp>/` folder.
