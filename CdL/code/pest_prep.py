"""
PEST prep step 1 — build the 'org' model dir with EXTERNAL input arrays/lists so
pyEMU PstFrom can template them.  Loads the already-written MF6 model, flips all
data to external files, and writes a clean copy to CdL_pest\org (inputs only).
Also drops in example output files (heads csv, zonal-AET csv) that add_observations
will read for their format.
"""
import shutil
from pathlib import Path
import flopy

WS   = Path(r"E:\00code_ws\DRYAD\CdL_model")
PEST = Path(r"E:\00code_ws\DRYAD\CdL_pest")
ORG  = PEST / "org"
MF6  = r"C:\00MODFLOW\mf6.7.0_win64\bin\mf6.exe"

if ORG.exists():
    shutil.rmtree(ORG)
ORG.mkdir(parents=True)

print(">> loading MF6 simulation …")
sim = flopy.mf6.MFSimulation.load(sim_ws=str(WS), exe_name=MF6, verbosity_level=0)
sim.set_sim_path(str(ORG))                 # MUST precede set_all_data_external (else it no-ops)
print(">> set_all_data_external …")
sim.set_all_data_external(check_data=False)
print(">> writing external model to", ORG)
sim.write_simulation()

# copy the voronoi grid pkl (forward-run AET post-proc needs it) + example outputs
shutil.copy2(WS / "voronoi_grid.pkl", ORG / "voronoi_grid.pkl")
if (WS / "cdl_gwf.obs.head.csv").exists():
    shutil.copy2(WS / "cdl_gwf.obs.head.csv", ORG / "cdl_gwf.obs.head.csv")
if (WS / "cdl_gwf.obs.sfr.csv").exists():                # SFR inlet/outlet streamflow obs (build_pst reads its cols)
    shutil.copy2(WS / "cdl_gwf.obs.sfr.csv", ORG / "cdl_gwf.obs.sfr.csv")

# Collapse the reference head/sfr obs to STRESS-PERIOD-END rows (drop ATS sub-steps),
# so the positional .ins that build_pst generates expects a DETERMINISTIC NPER rows —
# exactly what forward_run.py collapses every run to.  Without this, realisations whose
# ATS sub-stepping differs from the reference misalign the .ins and fail. (obs_collapse.py)
import sys as _sys
_sys.path.insert(0, str(Path(__file__).parent))
from obs_collapse import period_end_times, collapse_obs_csv
_ends = period_end_times(str(WS / "cdl_gwf.tdis"))
for _o in ("cdl_gwf.obs.head.csv", "cdl_gwf.obs.sfr.csv"):
    if (ORG / _o).exists():
        _nb, _na = collapse_obs_csv(str(ORG / _o), _ends)
        print(f">> collapsed reference {_o}: {_nb} -> {_na} rows (period ends; {_nb - _na} ATS sub-steps dropped)")
# GHB/DRN-west total-flow obs (for the BC inequality) + an example bc_constraint.csv (the derived
# value forward_run.py writes; build_pst reads its format to register the 'less_than' constraint obs)
import pandas as _pd
for _f in ("cdl_gwf.obs.ghb.csv", "cdl_gwf.obs.drnw.csv"):
    if (WS / _f).exists():
        shutil.copy2(WS / _f, ORG / _f)


def _mean_flow(fn, key):
    p = WS / fn
    if not p.exists():
        return 0.0
    d = _pd.read_csv(p)
    cols = [c for c in d.columns if c.strip().lower() != "time" and key in c.lower()]
    return float(d[cols[0]].mean()) if cols else 0.0


_pd.DataFrame({"cname": ["ghb_lt_drn"],
               "value": [_mean_flow("cdl_gwf.obs.ghb.csv", "ghb") + _mean_flow("cdl_gwf.obs.drnw.csv", "drn")]}
              ).to_csv(ORG / "bc_constraint.csv", index=False)
print(f">> wrote example bc_constraint.csv (ghb_lt_drn = "
      f"{_mean_flow('cdl_gwf.obs.ghb.csv', 'ghb') + _mean_flow('cdl_gwf.obs.drnw.csv', 'drn'):.1f} m3/d)")

# example zonal-AET csv (same routine the forward run will call)
import sys
sys.path.insert(0, str(Path(__file__).parent))
from model_aet_zonal import compute_zonal_aet
aet = compute_zonal_aet(WS, PEST)
aet.to_csv(ORG / "model_aet_zonal.csv")
print(f">> wrote example model_aet_zonal.csv ({aet.shape[0]} months x {aet.shape[1]} zones)")

# report the external files that build_pst.py will parameterize
print("\n== external files in org (npf / drn / ghb) ==")
for f in sorted(ORG.glob("*.txt")):
    n = f.name.lower()
    if any(k in n for k in ("npf", "drn", "ghb")):
        print("  ", f.name)
