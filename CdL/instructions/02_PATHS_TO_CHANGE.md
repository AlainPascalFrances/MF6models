# Paths to change on a new machine

All hardcoded absolute paths in the scripts reduce to **six prefixes**. Do a
find-and-replace for each across `CdL/code/*.py`. Nothing else needs editing.

> Tip: most editors do folder-wide replace (VS Code: Ctrl+Shift+H). Replace the **prefix**
> only, leaving the rest of each path intact.

### Automated alternative — `CdL/port_paths.py`

Instead of editing by hand, run the helper script [`../port_paths.py`](../port_paths.py). It has the
six prefixes below pre-loaded as its search strings. Workflow:

1. Run it as-is (`DRY_RUN = True`, all NEW values empty) — it **scans** and reports where each old
   prefix appears across `CdL/code`, so you can see exactly what will change.
2. Fill in the `NEW` value for each prefix in the `REPLACEMENTS` list.
3. Run again to **preview** the edits, then set `DRY_RUN = False` to **apply** them.

It replaces both the backslash and forward-slash forms of each Windows path, only touches text files
(skips binaries and `.git`), applies the most specific prefix first, and keeps a `.bak` of every
modified file. Point `TARGET_DIR` at whichever folder you unpacked the code into.

| # | Find (current prefix) | Replace with | What it is |
|---|---|---|---|
| 1 | `E:\00code_ws\DRYAD` | your `<DATA_ROOT>` (e.g. `D:\DRYAD`) | The two working dirs `CdL_model` (model workspace/outputs) and `CdL_pest` (PEST dirs) live under here. **Also written as forward slashes** in a few spots — check both `E:\00code_ws\DRYAD` and `E:/00code_ws/DRYAD`. |
| 2 | `E:\00code\flopy\dryad_cdl` | your `<CODE_DIR>` = the `CdL/code` folder of this repo | Only in `build_pst.py` (`CODE = …`), used to copy `forward_run.py` into the PEST template. |
| 3 | `E:\zzCloud\OneDrive - LNEG …\DRYAD` | your `<GIS_ROOT>` | Master GeoPackage, DEM, forcing CSVs and soils rasters. If you use the bundled `input_data`, point the individual variables at `input_data/…` instead (see below). Appears with both `\` and `/`. |
| 4 | `E:\ArcGis_Data\WorkSpace\COS\COSc2025\COSc_2025_N3_v0_TM06.tif` | your COS raster path | The national land-cover raster (see `LARGE_DATA_MANIFEST.md`). |
| 5 | `C:\00MODFLOW` | your MODFLOW/PEST++ install root | Executables: `mf6.exe`, `triangle.exe`, `pestpp-ies.exe`. |
| 6 | `C:\miniconda3\envs\flopy` | your conda env prefix | The Python that PEST++ calls per forward run (`build_pst.py` `model_command`). Must be the env with FloPy installed. |

## Executables (prefix #5) — set each explicitly

| Variable | File | Set to |
|---|---|---|
| `MF6_EXE` | `cdl_gwf_model_fable_v2.py`; `MF6` in `pest_prep.py` | `…\mf6.exe` (MODFLOW 6 ≥ 6.7.0) |
| `TRIANGLE_EXE` | `cdl_gwf_model_fable_v2.py` | `…\triangle.exe` |
| `PESTPP_IES` | `run_ies.py` | `…\pestpp-ies.exe` |

## Data variables (prefix #3) — where the bundled files map

If you copy `input_data/` somewhere and want to point straight at it, set these in
`cdl_gwf_model_fable_v2.py` (they are all near the top config block, ~lines 85–120 and ~240, ~300):

| Variable | Bundled file |
|---|---|
| `WATERSHED_GPKG` | `input_data/gis/dryad_modelo_NbS.gpkg` |
| `SOLOS_DIR` (→ `ks/ths/wp/fc_cdl.tif`) | `input_data/gis/` |
| `P_CSV`, `ET0_CSV` | `input_data/forcing/` |
| `VORONOI_LAYERS` | `input_data/grid_cache/voronoi_layers.npz` |
| `DONOR_STREAMFLOW_CSV` | `input_data/snirh/streamflow_21F_01H.csv` |
| `DEM_TIF` | **not bundled** — see `LARGE_DATA_MANIFEST.md` |
| `COS_TIF` (prefix #4) | **not bundled** — see `LARGE_DATA_MANIFEST.md` |

Also point these calibration-side paths (they use `<DATA_ROOT>`, prefix #1, so #1 usually covers them):
`correct_synthetic_piezo.py` (GRB + GPKG), `postprocess_cdl.py` (`SYN_PIEZO_CSV`, GPKG),
`make_obs.py` (`SYN_PIEZO`, `last_sim_start.txt`), `extract_pest_optimised.py`.

## The one path pyEMU writes for you — do NOT hand-edit

`build_pst.py` sets `pst.model_command = r"…python.exe forward_run.py"` (prefix #6). That is the
only place the per-forward-run Python is named. Set it to your env's `python.exe`.
