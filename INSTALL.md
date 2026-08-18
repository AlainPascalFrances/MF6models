# INSTALL — cook-book for a fresh machine

Step-by-step to go from a bare Windows machine to running the CdL model + PEST++ calibration.
Every version below is the one actually used (captured 2026-08-17); newer point releases usually
work, but these are the known-good set.

> Windows is assumed (the executables are the `_win64` builds). On Linux/macOS the conda side is
> identical; download the matching MODFLOW/PEST++ builds instead of the Windows ones.

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
│   ├── CdL_model\                      ← model workspace  (prefix #1 → D:\DRYAD\work)
│   └── CdL_pest\                       ← PEST org/template/master/workers
└── gis\                               ← large GIS rasters (DEM, COS) — see LARGE_DATA_MANIFEST.md
```

With this layout the **six path prefixes** (see `CdL/instructions/02_PATHS_TO_CHANGE.md`) become:

| # | Old prefix | New value |
|---|---|---|
| 1 | `E:\00code_ws\DRYAD` | `D:\DRYAD\work` |
| 2 | `E:\00code\flopy\dryad_cdl` | `D:\DRYAD\MF6models\CdL\code` |
| 3 | `E:\zzCloud\OneDrive …\DRYAD` | `D:\DRYAD\gis` (or point vars at `input_data\`) |
| 4 | `E:\ArcGis_Data\…COS…tif` | `D:\DRYAD\gis\COSc_2025_N3_v0_TM06.tif` |
| 5 | `C:\00MODFLOW` | `C:\sw\MODFLOW` |
| 6 | `C:\miniconda3\envs\flopy` | `C:\sw\miniconda3\envs\mf6models` |

---

## 1. Miniconda (Python 3.12)

1. Download the **Miniconda3 Windows 64-bit** installer: https://docs.conda.io/en/latest/miniconda.html
   (base Python 3.12; the repo env pins Python **3.12.13**).
2. Install to `C:\sw\miniconda3` (tick *"Add to PATH"* only if you know you want it; otherwise use the
   *Anaconda Prompt*).
3. Open an **Anaconda Prompt** and verify:
   ```bash
   conda --version
   ```

---

## 2. Clone the repository

```bash
cd /d D:\DRYAD
git clone https://github.com/AlainPascalFrances/MF6models.git
```
(or download the ZIP from GitHub and extract to `D:\DRYAD\MF6models`).

---

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

> If the solver is slow, install **mamba** first (`conda install -n base -c conda-forge mamba`) and use
> `mamba env create -f environment.yml`. If an exact pin can't be found on your platform, relax that
> line to `>=`.

---

## 4. MODFLOW 6  (v6.7.0)

1. Download the MODFLOW 6 **6.7.0** Windows release: https://github.com/MODFLOW-USGS/modflow6/releases
2. Extract so that `mf6.exe` sits at `C:\sw\MODFLOW\mf6.7.0\bin\mf6.exe`.
3. Test:
   ```bash
   C:\sw\MODFLOW\mf6.7.0\bin\mf6.exe -v      # -> mf6.exe: 6.7.0
   ```

---

## 5. Triangle

Grid generation (Voronoi via FloPy's Triangle wrapper).

1. Download Triangle: https://www.cs.cmu.edu/~quake/triangle.html (Windows build, or compile).
2. Place `triangle.exe` at `C:\sw\MODFLOW\triangle\triangle.exe`.

> Optional: if you reuse the bundled `input_data/grid_cache/voronoi_grid.pkl`, the model skips
> Triangle entirely — you only need it to re-mesh from scratch.

---

## 6. PEST++  (pestpp-ies 5.2.16)

1. Download the PEST++ suite release: https://github.com/usgs/pestpp/releases (v5.2.16 used).
2. Extract the `pestpp-*.exe` into `C:\sw\MODFLOW\pestpp-5.2.16\`.
3. Test:
   ```bash
   C:\sw\MODFLOW\pestpp-5.2.16\pestpp-ies.exe      # prints: version: 5.2.16
   ```

---

## 7. Point the scripts at your machine

Do the **six find-and-replace** prefixes from the table in §0 across `CdL/code/*.py`
(VS Code: Ctrl+Shift+H, scope = the `code` folder). Full detail + the per-variable executable
list is in `CdL/instructions/02_PATHS_TO_CHANGE.md`.

---

## 8. Seed the working directories

Create `D:\DRYAD\work\CdL_model` and `D:\DRYAD\work\CdL_pest`, then copy the bundled inputs into
them per the table in `CdL/instructions/01_SETUP.md §3`. Obtain the two large rasters
(`LARGE_DATA_MANIFEST.md`) and set `DEM_TIF` / `COS_TIF`.

---

## 9. Smoke test

```bash
conda activate mf6models
cd /d D:\DRYAD\MF6models\CdL\code
python cdl_gwf_model_fable_v2.py          # USE_PEST_PARAMS = False  → base run
```

A successful base run writes MF6 outputs under `D:\DRYAD\work\CdL_model\_output\<stamp>\`.
Then follow **`CdL/instructions/03_RUN_ORDER.md`** for the PEST++ chain.

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
