# CdL MODFLOW 6 model — input files

The complete set of MODFLOW 6 **input** files needed to run the calibrated CdL model,
as written by `CdL/code/cdl_gwf_model_fable_v2.py` (with `USE_PEST_PARAMS = True`, i.e. the
Tikhonov-regularized calibrated parameters — realisation 73).

## Run

```
gunzip cdl_gwf.uzf.gz          # see note below
mf6                            # run in this folder (mf6.exe on PATH)
```

## Contents

- `mfsim.nam`, `cdl_gwf.nam` — simulation and model name files
- `cdl_gwf.tdis`, `cdl_gwf.tdis.ats` — time discretisation (557 monthly periods; ATS)
- `cdl_gwf.ims`, `cdl_gwf.oc` — solver, output control
- `cdl_gwf.disv` — DISV grid (6 609 cells/layer, 3 layers)
- `cdl_gwf.npf`, `cdl_gwf.sto`, `cdl_gwf.ic` — flow, storage, initial conditions
- `cdl_gwf.drn`, `cdl_gwf_0.drn`, `cdl_gwf_1.drn` — seepage / western / secondary drains
- `cdl_gwf.ghb`, `cdl_gwf.sfr`, `cdl_gwf.lak` (+ `*.tab`), `cdl_gwf.mvr` — GHB, SFR, LAK, MVR
- `cdl_gwf.uzf.gz` — UZF (unsaturated zone; see note)
- `cdl_gwf.obs`, `cdl_gwf.ghb.obs`, `cdl_gwf.sfr.obs`, `cdl_gwf_0.drn.obs` — observation configs

## Notes

- **`cdl_gwf.uzf` is gzipped.** Uncompressed it is ~285 MB (inline per-period data for 10 393
  UZF cells over 557 stress periods), which exceeds GitHub's 100 MB limit; it compresses to
  ~13 MB. **Decompress it (`gunzip cdl_gwf.uzf.gz`) before running.**
- **Outputs are not included** — heads (`.hds`), budget (`.cbc`), package budgets (`.bud`),
  listing (`.lst`, `mfsim.lst`) and the binary grid (`.grb`) are written when the model runs.
- To rebuild these inputs from scratch (grid, parameters, forcing), run
  `CdL/code/cdl_gwf_model_fable_v2.py`.
