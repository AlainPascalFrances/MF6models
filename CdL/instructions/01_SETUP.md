# 01 — Setup on a new machine (install cook-book)

Step-by-step to go from a bare machine to running the CdL model + PEST++ calibration.
Every version below is the one actually used (captured 2026-08-17); newer point releases usually
work, but these are the known-good set.

> Windows is assumed (the executables are the `_win64` builds). On Linux/macOS the conda side is
> identical; download the matching MODFLOW/PEST++ builds instead of the Windows ones.
>
> The scripts are run interactively (the author uses **Spyder**). After editing any flag in a
> script, **reload the file in Spyder** before re-running so the change takes effect.

---

## 0. Proposed folder architecture

Keep **tools**, the **repo**, and the **working data** separate. Suggested layout:

```
C:\sw\                                  ← software / tools
├── miniconda3\                         ← Miniconda (Python 3.12 base + the 'mf6models' env)
└── MODFLOW\                            ← all binaries in one place (this becomes prefix #5 = C:\sw\MODFLOW)
    ├── mf6.7.0\bin\mf6.exe
    ├── triangle\triangle.exe
    └── pestpp-5.2.16\pestpp-ies.exe (+ the other pestpp-*.exe)

D:\DRYAD\                               ← project root (any drive with ≥ ~60 GB free)
├── MF6models\                          ← this git repository (code + input_data + instructions)
│   └── CdL\ ...
├── work\                               ← WORKING dirs, kept OUT of git (multi-GB outputs)
│   ├── CdL_model\                      ← model workspace  (this is <DATA_ROOT>\CdL_model)
│   └── CdL_pest\                       ← PEST org/template/master/workers
└── gis\                               ← large GIS rasters (DEM, COS) — see ../input_data/LARGE_DATA_MANIFEST.md
```

With this layout the **six path prefixes** (full detail in [`02_PATHS_TO_CHANGE.md`](02_PATHS_TO_CHANGE.md))
become:

| # | Old prefix | New value |
|---|---|---|
| 1 | `E:\00code_ws\DRYAD` | `D:\DRYAD\work` (holds `CdL_model` + `CdL_pest`) |
| 2 | `E:\00code\flopy\dryad_cdl` | `D:\DRYAD\MF6models\CdL\code` |
| 3 | `E:\zzCloud\OneDrive …\DRYAD` | `D:\DRYAD\gis` (or point vars at `input_data\`) |
| 4 | `E:\ArcGis_Data\…COS…tif` | `D:\DRYAD\gis\COSc_2025_N3_v0_TM06.tif` |
| 5 | `C:\00MODFLOW` | `C:\sw\MODFLOW` |
| 6 | `C:\miniconda3\envs\flopy` | `C:\sw\miniconda3\envs\mf6models` |

`<DATA_ROOT>` used below = prefix #1 = `D:\DRYAD\work`.

---

## 1. Miniconda (Python 3.12)

1. Download the **Miniconda3 Windows 64-bit** installer: https://docs.conda.io/en/latest/miniconda.html
   (base Python 3.12; the repo env pins Python **3.12.13**).
2. Install to `C:\sw\miniconda3` (or use the *Anaconda Prompt* it provides).
3. Open an **Anaconda Prompt** and verify: `conda --version`.

## 2. Clone the repository

```bash
cd /d D:\DRYAD
git clone https://github.com/AlainPascalFrances/MF6models.git
```
(or download the ZIP from GitHub and extract to `D:\DRYAD\MF6models`).

## 3. Create the conda environment

From the repo root:

```bash
cd /d D:\DRYAD\MF6models
conda env create -f environment.yml
conda activate mf6models
```

This installs the exact stack: FloPy 3.10.0, pyEMU 1.4.0, geopandas 1.1.3, rasterio 1.5.0,
scipy 1.17.1, numpy 2.4.6, pandas 2.3.3, matplotlib 3.10.9, pyflwdir 0.5.12, shapely 2.1.2,
fiona 1.10.1, pyproj 3.7.2, netCDF4 1.7.4, xarray 2026.4.0, openpyxl 3.1.5, requests 2.34.2, Spyder.

Verify:
```bash
python -c "import flopy, pyemu, geopandas, rasterio, scipy, matplotlib; print('env OK', flopy.__version__, pyemu.__version__)"
```

> If the solver is slow, install **mamba** (`conda install -n base -c conda-forge mamba`) and use
> `mamba env create -f environment.yml`. If an exact pin can't be found on your platform, relax that
> line in `environment.yml` to `>=`.

## 4. MODFLOW 6  (v6.7.0)

1. Download the MODFLOW 6 **6.7.0** Windows release: https://github.com/MODFLOW-USGS/modflow6/releases
2. Extract so that `mf6.exe` sits at `C:\sw\MODFLOW\mf6.7.0\bin\mf6.exe`.
3. Test: `C:\sw\MODFLOW\mf6.7.0\bin\mf6.exe -v`  → `mf6.exe: 6.7.0`.

## 5. Triangle

Grid generation (Voronoi via FloPy's Triangle wrapper).

1. Download Triangle: https://www.cs.cmu.edu/~quake/triangle.html (Windows build, or compile).
2. Place `triangle.exe` at `C:\sw\MODFLOW\triangle\triangle.exe`.

> Optional: if you reuse the bundled `input_data/grid_cache/voronoi_grid.pkl`, the model skips
> Triangle entirely — you only need it to re-mesh from scratch.

## 6. PEST++  (pestpp-ies 5.2.16)

1. Download the PEST++ suite release: https://github.com/usgs/pestpp/releases (v5.2.16 used).
2. Extract the `pestpp-*.exe` into `C:\sw\MODFLOW\pestpp-5.2.16\`.
3. Test: `C:\sw\MODFLOW\pestpp-5.2.16\pestpp-ies.exe`  → `version: 5.2.16`.

## 7. Point the scripts at your machine

Do the **six find-and-replace** prefixes from §0 across `CdL/code/*.py`
(VS Code: Ctrl+Shift+H, scope = the `code` folder). The per-variable executable list and the data
variables are in [`02_PATHS_TO_CHANGE.md`](02_PATHS_TO_CHANGE.md).

## 8. Seed the working directories

Create `<DATA_ROOT>\CdL_model` and `<DATA_ROOT>\CdL_pest`, then copy the bundled inputs into them:

| Copy this bundled file… | …to here |
|---|---|
| `input_data/grid_cache/voronoi_grid.pkl` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/cdl_gwf.disv.grb` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/spinup_heads.npy` | `<DATA_ROOT>/CdL_model/` |
| `input_data/grid_cache/voronoi_layers.npz` | `<DATA_ROOT>/CdL_model/conceptual/layers/` |
| `input_data/pest/cos_zones.npz` | `<DATA_ROOT>/CdL_pest/` |
| `input_data/pest/pest_optimised.npz` | `<DATA_ROOT>/CdL_pest/` *(over-fit — see [`04_CALIBRATION_STATE.md`](04_CALIBRATION_STATE.md))* |
| `input_data/snirh/*.csv` | `<DATA_ROOT>/CdL_pest/snirh_data_availability/` |

The `gis/` and `forcing/` bundled files can stay in `input_data/` — just point the model's
`WATERSHED_GPKG`, `SOLOS_DIR`, `P_CSV`, `ET0_CSV` variables at them (prefix #3). Obtain the two large
rasters (DEM, COS) and set `DEM_TIF` / `COS_TIF` — see [`../input_data/LARGE_DATA_MANIFEST.md`](../input_data/LARGE_DATA_MANIFEST.md).

## 9. Smoke test

```bash
conda activate mf6models
cd /d D:\DRYAD\MF6models\CdL\code
python cdl_gwf_model_fable_v2.py          # USE_PEST_PARAMS = False  → base run
```

A successful base run writes MF6 outputs under `<DATA_ROOT>\CdL_model\_output\<stamp>\`.
Then follow **[`03_RUN_ORDER.md`](03_RUN_ORDER.md)** for the PEST++ chain.

---

## Version summary (known-good)

| Software | Version | Software | Version |
|---|---|---|---|
| Python | 3.12.13 | MODFLOW 6 | 6.7.0 |
| FloPy | 3.10.0 | Triangle | any |
| pyEMU | 1.4.0 | PEST++ (pestpp-ies) | 5.2.16 |
| geopandas | 1.1.3 | rasterio | 1.5.0 |
| numpy | 2.4.6 | pandas | 2.3.3 |
| scipy | 1.17.1 | matplotlib | 3.10.9 |
| shapely | 2.1.2 | fiona | 1.10.1 |
| pyproj | 3.7.2 | pyflwdir | 0.5.12 |
| netCDF4 | 1.7.4 | xarray | 2026.4.0 |
