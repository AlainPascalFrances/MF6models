# Paths to change on a new machine

All machine-specific paths are centralised in **one file — [`CdL/code/config.py`](../code/config.py)**.
Every script does `import config` and reads its paths from there, so moving to a new machine means
editing **three settings** (or setting three environment variables). There is **no find-and-replace**
to do across the scripts, and nothing else to edit.

## The three settings

Open [`../code/config.py`](../code/config.py) and edit the three roots at the top (or leave them and
set the matching environment variables, which take precedence):

| Setting | Env var | What it is | Example |
|---|---|---|---|
| `BASE` | `CDL_BASE` | Project **data root** — the folder that holds the two working dirs `CdL_model` (model workspace/outputs, plus `gis/` and `forcing/` under it) and `CdL_pest` (PEST org/template/master/workers). | `D:\DRYAD\work` |
| `MODFLOW_DIR` | `CDL_MODFLOW_DIR` | Folder holding the **MODFLOW 6 / Triangle / PEST++** executables. | `C:\sw\MODFLOW` |
| `PYTHON_EXE` | `CDL_PYTHON_EXE` | The **conda-env `python.exe`** that PEST++ workers call to run `forward_run.py`. | `C:\sw\miniconda3\envs\mf6models\python.exe` |

Everything else is **derived** from these (you normally don't touch it):

| Derived (in `config.py`) | Value |
|---|---|
| `CODE` | this repo's `CdL/code` folder (self-locating — never edit) |
| `MODEL` | `BASE / "CdL_model"` |
| `PEST` | `BASE / "CdL_pest"` |
| `MODIS_DIR` | `MODEL / "RS" / "MODIS4CDL"` |
| `MF6_EXE` | `MODFLOW_DIR / "mf6.7.0_win64" / "bin" / "mf6.exe"` |
| `TRIANGLE_EXE` | `MODFLOW_DIR / "win64" / "triangle.exe"` |
| `PESTPP_IES` | `MODFLOW_DIR / "pestpp-5.2.27-win" / "bin" / "pestpp-ies.exe"` |

## Two things to check inside `config.py`

1. **Executable sub-folder names.** `MF6_EXE`, `TRIANGLE_EXE` and `PESTPP_IES` assume the exact
   sub-folders shown above under `MODFLOW_DIR` (`mf6.7.0_win64\bin`, `win64`, `pestpp-5.2.27-win\bin`).
   Either lay your `MODFLOW_DIR` out the same way, or edit those three lines to match where your
   binaries actually sit.
2. **Environment variables win.** If `CDL_BASE` / `CDL_MODFLOW_DIR` / `CDL_PYTHON_EXE` are set in your
   shell, they override the hardcoded defaults — handy for switching machines without editing the file.

## Where the bundled `input_data/` maps

The scripts expect the data under `BASE\CdL_model` and `BASE\CdL_pest`. Copy the bundled inputs there
(the full seeding table is in [`01_SETUP.md`](01_SETUP.md) §8):

| Bundled file | Goes to (relative to `config.MODEL` / `config.PEST`) |
|---|---|
| `input_data/gis/dryad_modelo_NbS.gpkg`, soils `*_cdl.tif` | `CdL_model/gis/` |
| `input_data/forcing/p_month…csv`, `et0_month…csv` | `CdL_model/forcing/` |
| `input_data/grid_cache/*` | `CdL_model/` (and `voronoi_layers.npz` → `CdL_model/conceptual/layers/`) |
| `input_data/pest/*.npz` | `CdL_pest/` |
| `input_data/snirh/*` | `CdL_pest/snirh_data_availability/` |
| DEM (`dem_cdl.tif`), COS raster | **not bundled** — see [`../input_data/LARGE_DATA_MANIFEST.md`](../input_data/LARGE_DATA_MANIFEST.md) |

Individual data-file paths are defined at the top of each script **relative to `config.MODEL` /
`config.PEST`** — so once the three roots are right and the files are seeded, no per-script edits are
needed. (Legacy note: the old manual six-prefix / `port_paths.py` workflow has been retired in favour
of `config.py`.)
