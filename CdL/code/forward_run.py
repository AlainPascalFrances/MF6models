"""
PEST++ forward run (executed in every worker dir).
  1. apply pilot-point (K) + BC multipliers to the external MF6 input files
  2. enforce Kv = KV_RATIO * Kh so the fixed anisotropy is preserved after the Kh mults
  3. run MODFLOW 6
  4. post-process the zonal actual-ET (heads are written directly by the MF6 OBS pkg)
"""
import os
import glob
import numpy as np
import pyemu

MF6 = r"C:\00MODFLOW\mf6.7.0_win64\bin\mf6.exe"
KV_RATIO = 0.1                      # k33 = KV_RATIO * k  (model's fixed Kh/Kv = 10)


def main():
    # 1. multipliers -> model input files (PstFrom writes mult2model_info.csv)
    pyemu.helpers.apply_list_and_array_pars(arr_par_file="mult2model_info.csv", chunk_len=50)

    # 2. Kv follows Kh (re-derive every k33 layer file from the just-updated k file)
    for kf in glob.glob("*.npf_k_layer*.txt") + glob.glob("*.npf_k.txt"):
        k33f = kf.replace("npf_k_layer", "npf_k33_layer").replace("npf_k.txt", "npf_k33.txt")
        if os.path.exists(k33f):
            np.savetxt(k33f, KV_RATIO * np.loadtxt(kf), fmt="%15.6E")

    # 3. run MODFLOW 6.  NOTE: pass the exe UNQUOTED — pyemu.os_utils.run splits on
    # spaces and requires the first token to end with "exe"; surrounding quotes make it
    # end with 'exe"', so pyemu mangles the command to ...mf6.exe".exe and it no-ops.
    # (The path has no spaces, so quotes aren't needed.)
    pyemu.os_utils.run(MF6)

    # 4. zonal AET (cos_zones.npz + voronoi_grid.pkl are copied into each worker)
    from model_aet_zonal import compute_zonal_aet
    compute_zonal_aet(".", ".").to_csv("model_aet_zonal.csv")

    # 5. BC inequality constraint (user 2026-08-14): GHB inflow < DRN-west outflow.
    #    MF6 ghb obs is +into-aquifer, drn obs is -out -> derived = mean(ghb) + mean(drn),
    #    which is < 0 exactly when |DRN-west out| > GHB in. build_pst puts it in a pestpp
    #    'less_than' group (obsval 0) so it is penalised ONLY when the constraint is violated.
    import pandas as pd

    def _mean_flow(fn, key):
        if not os.path.exists(fn):
            return 0.0
        d = pd.read_csv(fn)
        cols = [c for c in d.columns if c.strip().lower() != "time" and key in c.lower()]
        return float(d[cols[0]].mean()) if cols else 0.0

    _ghb = _mean_flow("cdl_gwf.obs.ghb.csv", "ghb")
    _drn = _mean_flow("cdl_gwf.obs.drnw.csv", "drn")
    pd.DataFrame({"cname": ["ghb_lt_drn"], "value": [_ghb + _drn]}).to_csv("bc_constraint.csv", index=False)


if __name__ == "__main__":
    # freeze_support() is REQUIRED here: pestpp calls this as bare `python.exe
    # forward_run.py`, and apply_list_and_array_pars() spawns a multiprocessing
    # Pool.  Without freeze_support the process silently terminates right after the
    # apply step (exit 0, MF6 never runs, no obs files written) -> pestpp reports
    # "output file model_aet_zonal.csv not found" and the whole run fails.
    import multiprocessing
    multiprocessing.freeze_support()
    main()
