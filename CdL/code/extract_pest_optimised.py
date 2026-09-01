"""
Extract the PEST++ (pestpp-ies) OPTIMISED parameter set into a small npz that
cdl_gwf_model_fable_v2.py can load (USE_PEST_PARAMS=True) to run the model with
the calibrated parameters.

We take the MINIMUM-phi realisation of the final IES iteration (the single best
fit).  For the K pilot points we re-run the pyEMU apply once (fill templates ->
apply_list_and_array_pars) to get the kriged optimised K field, then express it
as a per-cell multiplier on the model's base K.  The BC parameters are constants,
so their optimised value IS the multiplier (or additive offset for ghb_head).

Output:  E:\\00code_ws\\DRYAD\\CdL_pest\\pest_optimised.npz
  kh_mult (3, ncpl)  ghb_cond_mult  ghb_head_add  drnout_mult  drnseep_mult
  drnwest_u1_mult  drnwest_dp_mult  drnsec_u1_mult  drnsec_dp_mult  sfr_inflow_mult  (+ metadata)

CAVEAT: in this IES 5/9 param groups railed to their bounds (uninformative data)
-> these are data-fitting values, not physical estimates.  Kept for inspection
and to exercise the import machinery for the next (data-driven) calibration.
"""
import config
import shutil, multiprocessing
from pathlib import Path
import numpy as np, pandas as pd, pyemu

PEST     = Path(str(config.PEST))
TEMPLATE = PEST / "template"
MASTER   = PEST / "master"
ORG      = PEST / "org"
OUT_NPZ  = PEST / "pest_optimised.npz"
ITER     = 3                                    # final IES iteration


def _flat(p):
    return np.array(Path(p).read_text().split(), dtype=float)


def main():
    # ---- 1. pick the minimum-phi realisation of the final iteration ----
    phi = pd.read_csv(MASTER / "cdl.phi.actual.csv")
    row = phi[phi["iteration"] == ITER].iloc[0]
    meta = ["iteration", "total_runs", "mean", "standard_deviation", "min", "max"]
    real_phi = row.drop(labels=meta).astype(float)
    best = real_phi.idxmin()
    print(f">> iteration {ITER}: best realisation = '{best}'  (phi {real_phi[best]:.0f} "
          f"of {real_phi.min():.0f}-{real_phi.max():.0f})")

    par = pd.read_csv(MASTER / f"cdl.{ITER}.par.csv", index_col=0)
    par.index = par.index.astype(str)
    pvals = par.loc[str(best)]                  # Series: parname -> value

    # ---- 2. BC scalar parameters (their value = the multiplier / offset) ----
    def one(sub, default=None):
        hit = [c for c in pvals.index if sub in c.lower()]
        if not hit and default is not None:      # e.g. drnsec absent in a pre-2026-08-05 posterior
            return float(default)
        assert len(hit) == 1, f"{sub}: expected 1 par, got {hit}"
        return float(pvals[hit[0]])
    ghb_cond_mult = one("ghb_cond"); ghb_head_add = one("ghb_head")
    drnout_mult = one("drnout", default=1.0)     # DRN-OUT retired; kept for loader back-compat
    drnseep_mult = one("drnseep")
    # DRN-WEST / DRN-SECONDARY are SPLIT BY LAYER since build_pst 2026-08-15 (U1 vs deep tune
    # independently), so one("drnwest") would now match BOTH u1+dp and trip the len==1 assert.
    # Extract all four; the secondary outlet's L1 is inactive so it has NO drnsec_u1 row -> 1.0.
    drnwest_u1_mult = one("drnwest_u1", default=1.0); drnwest_dp_mult = one("drnwest_dp", default=1.0)
    drnsec_u1_mult  = one("drnsec_u1",  default=1.0); drnsec_dp_mult  = one("drnsec_dp",  default=1.0)
    sfr_inflow_mult = one("sfr_inflow", default=1.0)   # inlet no longer calibrated (fixed observed forcing) -> 1.0

    # ---- 3. optimised K field via the pyEMU apply (kriges the pilot points) ----
    work = PEST / "_optimtmp"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(TEMPLATE, work)
    pst = pyemu.Pst(str(work / "cdl.pst"))
    pst.parameter_data.loc[:, "parval1"] = pvals.reindex(pst.parameter_data.index).astype(float)
    pst.write(str(work / "cdl.pst"), version=2)
    cwd = Path.cwd()
    import os
    os.chdir(work)
    try:
        pyemu.Pst(str(work / "cdl.pst")).write_input_files()
        pyemu.helpers.apply_list_and_array_pars(arr_par_file="mult2model_info.csv", chunk_len=50)
    finally:
        os.chdir(cwd)

    kh_mult = []
    for L in (1, 2, 3):
        opt = _flat(work / f"cdl_gwf.npf_k_layer{L}.txt")
        base = _flat(ORG / f"cdl_gwf.npf_k_layer{L}.txt")
        m = np.where(base > 0, opt / base, 1.0)
        kh_mult.append(m)
        print(f"   layer {L}: K mult range {m.min():.3f}–{m.max():.3f} (mean {m.mean():.3f})")
    kh_mult = np.vstack(kh_mult)                 # (3, ncpl)

    np.savez(OUT_NPZ,
             kh_mult=kh_mult,
             ghb_cond_mult=ghb_cond_mult, ghb_head_add=ghb_head_add,
             drnout_mult=drnout_mult, drnseep_mult=drnseep_mult,
             drnwest_u1_mult=drnwest_u1_mult, drnwest_dp_mult=drnwest_dp_mult,
             drnsec_u1_mult=drnsec_u1_mult, drnsec_dp_mult=drnsec_dp_mult,
             sfr_inflow_mult=sfr_inflow_mult,
             realisation=str(best), phi=float(real_phi[best]), iteration=ITER, ncpl=kh_mult.shape[1])
    shutil.rmtree(work, ignore_errors=True)

    print(f"\n>> wrote {OUT_NPZ}")
    print(f"   ghb_cond x{ghb_cond_mult:.3g}  ghb_head {ghb_head_add:+.3g} m  drnseep x{drnseep_mult:.3g}  "
          f"drnwest[u1 x{drnwest_u1_mult:.3g}/dp x{drnwest_dp_mult:.3g}]  "
          f"drnsec[u1 x{drnsec_u1_mult:.3g}/dp x{drnsec_dp_mult:.3g}]  sfr_inflow x{sfr_inflow_mult:.3g}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
