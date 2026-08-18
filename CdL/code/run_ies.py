"""
Launch pestpp-ies for the CdL calibration with local PANTHER workers.
Run in Spyder / a terminal (multi-hour): 8 workers x ~6-min MF6 run each.

  NOPTMAX = -1  -> prior Monte-Carlo only (recommended first: checks the ensemble
                   machinery + shows prior predictive spread), then set to 3 for IES.
  NOPTMAX =  0  -> single base run (initial parameters) -> base phi report (quick).
  NOPTMAX >  0  -> that many IES iterations (the actual calibration).
"""
import shutil
from pathlib import Path
import pyemu

PEST     = Path(r"E:\00code_ws\DRYAD\CdL_pest")
TEMPLATE = PEST / "template"
MASTER   = PEST / "master"
WORKERROOT = PEST / "workers"
PESTPP_IES = r"C:\00MODFLOW\pestpp\pestpp-ies.exe"

# ======================================================================================
# SERVER-MIGRATION CHECKLIST (moving off the laptop for the big 150-real IES) — edit these:
#   1. PEST / TEMPLATE / MASTER / WORKERROOT   base paths below
#   2. PESTPP_IES  path to pestpp-ies.exe on the server
#   3. build_pst.py line ~189  pst.model_command python path (C:\miniconda3\envs\flopy\python.exe)
#   4. NUM_WORKERS  = number of physical cores you want to give it (server has more -> raise)
#   Re-run the FULL chain on the server (model -> pest_prep -> build_pst -> base run -> make_obs)
#   so the prior_pe.jcb / prior_cov.jcb and org/template are rebuilt with the server's paths.
# ======================================================================================
NUM_WORKERS = 12       # raise on the server (one worker per physical core)
# NOPTMAX = 0 — single base run (~15 min): confirms the base phi with real targets before committing to the ensemble.
# NOPTMAX = -1 — prior Monte-Carlo: checks the ensemble machinery and shows prior spread.
# NOPTMAX = 3 — the IES calibration (hours-to-days at 150 reals).
NOPTMAX     = 3        # <-- set to 0 (base run), -1 (prior MC), or 3 (IES calibration)
NUM_REALS   = 150      # geostatistical prior ensemble (prior_pe.jcb holds 200; pestpp uses the first 150)
PORT        = 4269

if __name__ == "__main__":
    pst = pyemu.Pst(str(TEMPLATE / "cdl.pst"))
    pst.control_data.noptmax = NOPTMAX
    pst.pestpp_options["ies_num_reals"] = NUM_REALS
    # ---- ANTI-CHECKERBOARD / OVERFIT package (2026-08-16) ----
    # (a) GEOSTATISTICAL PRIOR: draw from the spatially-correlated ensemble built by build_pst.py
    #     (prior_pe.jcb) with its covariance (prior_cov.jcb), instead of pestpp's diagonal bounds prior.
    # (b) LOCALIZATION: automatic adaptive localization zeroes spurious param<->obs correlations, so a
    #     rank-limited ensemble can't push distant pilot points to their bounds (the 58%-at-bounds,
    #     negative-Moran's-I result). Needs no localizer matrix.
    if NOPTMAX != 0:                                        # (base run ignores the ensemble)
        pst.pestpp_options["ies_parameter_ensemble"] = "prior_pe.jcb"
        pst.pestpp_options["ies_autoadaloc"] = True        # automatic adaptive localization
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.5          # reject runaway realizations
    pst.pestpp_options["overdue_giveup_fac"] = 2.0         # kill hung/grinding MF6 runs FAST (was 5.0 = 6.5 h)
    pst.pestpp_options["overdue_giveup_minutes"] = 120.0    # ABSOLUTE cap: kill any run >120 min. CRITICAL —
    #   overdue_giveup_fac needs an avg_run_time baseline that doesn't exist until the 1st run finishes;
    #   without this cap, grinders ran 4-9 h in prior-MC #3 before anything could be killed.
    pst.write(str(TEMPLATE / "cdl.pst"), version=2)

    if MASTER.exists():
        shutil.rmtree(MASTER)
    pyemu.os_utils.start_workers(
        str(TEMPLATE), PESTPP_IES, "cdl.pst",
        num_workers=NUM_WORKERS, worker_root=str(WORKERROOT),
        master_dir=str(MASTER), port=PORT, cleanup=False,
    )
    print(">> pestpp-ies finished. Results in", MASTER)

    # After a SUCCESSFUL base run (NOPTMAX=0), auto-run make_obs.py to populate
    # targets + weights from the fresh base residuals — saves a manual step.
    if NOPTMAX == 0:
        base_obs = MASTER / "cdl.base.obs.csv"
        if base_obs.exists():
            import subprocess, sys
            make_obs = Path(__file__).resolve().parent / "make_obs.py"
            print(f">> base run OK ({base_obs.name} present) -> auto-running make_obs.py …")
            subprocess.run([sys.executable, str(make_obs)], check=True)
            print(">> make_obs.py done — weights set + noptmax=-1. Next: run this with NOPTMAX=-1 for the prior MC.")
        else:
            print(f">> !! base run did NOT produce {base_obs} — make_obs.py NOT run (check the base run).")
