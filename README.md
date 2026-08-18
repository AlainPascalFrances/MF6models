# MF6models

MODFLOW 6 (FloPy) groundwater-flow models and their PEST++ (pestpp-ies) calibration
workflows, by Alain Pascal Frances (LNEG, Portugal).

## Models

| Model | Description |
|-------|-------------|
| [`CdL/`](CdL/) | Catchment of **Casa de Lobos** — ~20 km² Voronoi-DISV MF6 model (UZF / SFR / LAK / MVR, monthly transient 1981–2026) with a pyEMU/pestpp-ies calibration against synthetic piezometers, MODIS AET and SNIRH-regionalised streamflow. |

## Start here

- **New machine, from scratch?** → **[`INSTALL.md`](INSTALL.md)** — a full cook-book (Miniconda,
  the conda env with pinned versions, MODFLOW 6, Triangle, PEST++), with a proposed disk layout.
- **What is this model?** → [`CdL/instructions/00_OVERVIEW.md`](CdL/instructions/00_OVERVIEW.md).
- **Already have the tools?** → [`CdL/instructions/01_SETUP.md`](CdL/instructions/01_SETUP.md) and the
  **six path prefixes** to change in [`02_PATHS_TO_CHANGE.md`](CdL/instructions/02_PATHS_TO_CHANGE.md).

## Layout

```
MF6models/
├── README.md                  (this file)
├── INSTALL.md                 (from-scratch install cook-book + disk layout)
├── environment.yml            (conda environment spec, versions pinned)
├── .gitignore
└── CdL/
    ├── code/                  (all Python scripts — the canonical model is cdl_gwf_model_fable_v2.py)
    ├── input_data/            (bundled inputs; two large GIS rasters are documented, not committed)
    │   ├── gis/               (master GeoPackage + clipped soils rasters)
    │   ├── forcing/           (monthly precipitation + ET0)
    │   ├── grid_cache/        (Voronoi grid, DISV grb, conceptual layers, spin-up heads)
    │   ├── snirh/             (streamflow donor + synthetic piezometer series)
    │   ├── pest/              (COS zone map + last calibrated parameter set)
    │   └── LARGE_DATA_MANIFEST.md   (the DEM + COS rasters: where to get, where to place)
    └── instructions/
        ├── 00_OVERVIEW.md
        ├── 01_SETUP.md
        ├── 02_PATHS_TO_CHANGE.md
        ├── 03_RUN_ORDER.md
        └── 04_CALIBRATION_STATE.md
```

## Requirements (summary)

- Python 3.12 conda env — see `environment.yml` (FloPy, pyEMU, geopandas, rasterio, …).
- **MODFLOW 6** ≥ 6.7.0 and **Triangle** executables.
- **pestpp-ies** (PEST++ suite) for the calibration.

## License / data

Code is provided for research reproducibility. The GIS/forcing data are LNEG/DRYAD project
data — see `LARGE_DATA_MANIFEST.md` for the two large rasters that are not redistributed here.
