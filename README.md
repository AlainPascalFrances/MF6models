# MF6models

MODFLOW 6 (FloPy) groundwater-flow models and their PEST++ (pestpp-ies) calibration
workflows, by Alain Pascal Frances (LNEG, Portugal).

## Models

| Model | Description |
|-------|-------------|
| [`CdL/`](CdL/) | Catchment of **Casa de Lobos** — ~20 km² Voronoi-DISV MF6 model (UZF / SFR / LAK / MVR, monthly transient 1981–2026) with a pyEMU/pestpp-ies calibration against synthetic piezometers, MODIS AET and SNIRH-regionalised streamflow. |

## Start here

- **What is this model?** → [`CdL/instructions/00_OVERVIEW.md`](CdL/instructions/00_OVERVIEW.md).
- **Just run the model, no Python** → [`CdL/input_model/`](CdL/input_model/) holds the complete set of
  MODFLOW 6 input files for the calibrated model; with `mf6.exe` alone you can reproduce the run
  (see its [`readme.txt`](CdL/input_model/readme.txt) — one file must be decompressed first).
- **Set up a new machine (install cook-book)** → **[`CdL/instructions/01_SETUP.md`](CdL/instructions/01_SETUP.md)** —
  Miniconda, the pinned conda env, MODFLOW 6, Triangle, PEST++ and a proposed disk layout. All
  machine-specific paths live in a single file, [`CdL/code/config.py`](CdL/code/config.py): edit its
  **three settings** (`BASE`, `MODFLOW_DIR`, `PYTHON_EXE`) or set the matching `CDL_*` environment
  variables — there is no find-and-replace across the scripts
  (detail in [`02_PATHS_TO_CHANGE.md`](CdL/instructions/02_PATHS_TO_CHANGE.md)).
- **What does each script do?** → [`CdL/code/readme.txt`](CdL/code/readme.txt) — a one-line index of
  every script, grouped by workflow stage.
- **Rebuild / run the whole chain** → [`CdL/instructions/03_RUN_ORDER.md`](CdL/instructions/03_RUN_ORDER.md);
  calibration status in [`04_CALIBRATION_STATE.md`](CdL/instructions/04_CALIBRATION_STATE.md).

## Layout

```
MF6models/
├── README.md                  (this file)
├── environment.yml            (conda environment spec, versions pinned)
├── .gitignore
└── CdL/
    ├── code/                  (all Python scripts — the canonical model is cdl_gwf_model_fable_v2.py)
    │   ├── config.py          (THE only machine-specific paths — edit this one file)
    │   └── readme.txt         (one-line description of every script, grouped by workflow)
    ├── input_model/           (ready-to-run MODFLOW 6 input files + readme.txt — needs only mf6.exe)
    ├── input_data/            (bundled inputs; two large GIS rasters are documented, not committed)
    │   ├── gis/               (master GeoPackage + clipped soils rasters)
    │   ├── forcing/           (monthly precipitation + ET0)
    │   ├── grid_cache/        (Voronoi grid, DISV grb, conceptual layers, spin-up heads)
    │   ├── snirh/             (streamflow donor + synthetic piezometer series)
    │   ├── pest/              (COS zone map + last calibrated parameter set)
    │   └── LARGE_DATA_MANIFEST.md   (the DEM + COS rasters: where to get, where to place)
    ├── docs/                  (model/PEST reports: strategy, initial parameters, crash report)
    └── instructions/
        ├── 00_OVERVIEW.md
        ├── 01_SETUP.md
        ├── 02_PATHS_TO_CHANGE.md   (the three config.py settings)
        ├── 03_RUN_ORDER.md
        └── 04_CALIBRATION_STATE.md
```

## Requirements (summary)

- To **only re-run the calibrated model** from `CdL/input_model/`: **MODFLOW 6** ≥ 6.7.0
  (`mf6.exe`) and nothing else — no Python, no conda.
- To **rebuild / post-process / calibrate**:
  - Python 3.12 conda env — see `environment.yml` (FloPy, pyEMU, geopandas, rasterio, …).
  - **MODFLOW 6** ≥ 6.7.0 and **Triangle** executables.
  - **pestpp-ies** (PEST++ suite) for the calibration.

## License / data

Code is provided for research reproducibility. The GIS/forcing data are LNEG/DRYAD project
data — see `LARGE_DATA_MANIFEST.md` for the two large rasters that are not redistributed here.
