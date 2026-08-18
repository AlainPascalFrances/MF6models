# 01 — Setup on a new machine

## 1. Python environment

From the repo root:

```bash
conda env create -f environment.yml
conda activate mf6models
```

This installs FloPy, pyEMU, geopandas, rasterio, scipy, matplotlib, pyflwdir, etc. (Python 3.12).

> The scripts are run interactively (the author uses **Spyder**). After editing any flag in a
> script, **reload the file in Spyder** before re-running so the change takes effect.

## 2. External executables (not in the conda env)

Download and install, then set their paths in the scripts (see `02_PATHS_TO_CHANGE.md`, prefix #5):

| Tool | Version | Get it from |
|---|---|---|
| **MODFLOW 6** | ≥ 6.7.0 | https://github.com/MODFLOW-USGS/modflow6/releases |
| **Triangle** | any | https://www.cs.cmu.edu/~quake/triangle.html (used to build the Voronoi grid; **skippable if you use the bundled `voronoi_grid.pkl`**) |
| **PEST++ (pestpp-ies)** | recent | https://github.com/usgs/pestpp/releases |

## 3. Working directories

The scripts read/write two working dirs (kept **outside** this repo — they hold multi-GB outputs):

- `<DATA_ROOT>/CdL_model` — model workspace (MF6 inputs, heads, budgets, `_output/`, `_input/`).
- `<DATA_ROOT>/CdL_pest`  — PEST dirs (`org/`, `template/`, `master/`, `workers/`).

Create both, then seed them with the bundled inputs:

| Copy this bundled file… | …to here |
|---|---|
| `input_data/grid_cache/voronoi_grid.pkl` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/cdl_gwf.disv.grb` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/spinup_heads.npy` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/voronoi_layers.npz` | `<DATA_ROOT>/CdL_model/conceptual/layers/` |
| `input_data/pest/cos_zones.npz` | `<DATA_ROOT>/CdL_pest/` |
| `input_data/pest/pest_optimised.npz` | `<DATA_ROOT>/CdL_pest/` *(see caveat in `04_CALIBRATION_STATE.md`)* |
| `input_data/snirh/*.csv` | `<DATA_ROOT>/CdL_pest/snirh_data_availability/` |

The `gis/` and `forcing/` bundled files can stay in `input_data/` — just point the model's
`WATERSHED_GPKG`, `SOLOS_DIR`, `P_CSV`, `ET0_CSV` variables at them (prefix #3).

## 4. Large GIS rasters (not bundled)

The DEM (`dem_cdl.tif`, 619 MB) and COS land-cover raster (187 MB) exceed GitHub's file limit.
See **`input_data/LARGE_DATA_MANIFEST.md`** for what they are, where to obtain them, and which
variables to point at them.

## 5. Verify

A quick sanity check that the env is healthy:

```bash
python -c "import flopy, pyemu, geopandas, rasterio, scipy, matplotlib; print('env OK')"
```

Then follow **`03_RUN_ORDER.md`**.
