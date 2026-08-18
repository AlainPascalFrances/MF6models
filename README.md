# MF6models

MODFLOW 6 (FloPy) groundwater-flow models and their PEST++ (pestpp-ies) calibration
workflows, by Alain Pascal Frances (LNEG, Portugal).

## Models

| Model | Description |
|-------|-------------|
| [`CdL/`](CdL/) | Catchment of **Casa de Lobos** — ~20 km² Voronoi-DISV MF6 model (UZF / SFR / LAK / MVR, monthly transient 1981–2026) with a pyEMU/pestpp-ies calibration against synthetic piezometers, MODIS AET and SNIRH-regionalised streamflow. |

## Start here

Open **[`CdL/instructions/01_SETUP.md`](CdL/instructions/01_SETUP.md)** — it walks through the
conda environment, the external executables (MODFLOW 6, Triangle, pestpp-ies), where to put
the input data, and the **six path prefixes** you must change to run on a new machine
(see [`02_PATHS_TO_CHANGE.md`](CdL/instructions/02_PATHS_TO_CHANGE.md)).

## Layout

```
MF6models/
├── README.md                  (this file)
├── environment.yml            (conda environment spec)
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
