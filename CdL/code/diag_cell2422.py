"""Nail the SFR-828 runaway: what are cell 2422 and reach 828, and what MVR providers feed reach 828?"""
import config
import numpy as np, flopy
from pathlib import Path
WS = Path(str(config.MODEL))
C = 2422   # 0-based cell id from the listing's (4,2422)
R = 828    # reach number from 1_GWF-SFR-(828)-inflow

sim = flopy.mf6.MFSimulation.load(sim_ws=str(WS), verbosity_level=0)
gwf = sim.get_model()
disv = gwf.get_package("DISV"); sfr = gwf.get_package("SFR")
top = disv.top.array; botm = disv.botm.array; idom = disv.idomain.array
npf = gwf.get_package("NPF"); k = npf.k.array
xc = np.array(gwf.modelgrid.xcellcenters); yc = np.array(gwf.modelgrid.ycellcenters)
print(f"--- cell {C}: x={xc[C]:.0f} y={yc[C]:.0f} top={top[C]:.2f}")
for L in range(botm.shape[0]):
    print(f"   L{L+1}: botm={botm[L,C]:7.2f}  idomain={idom[L,C]:+d}  K={k[L,C]:.2f}  thk={(top[C] if L==0 else botm[L-1,C])-botm[L,C]:.2f}")

# SFR reach R: cell + geometry + connections
pk = sfr.packagedata.array
RN = "ifno" if "ifno" in pk.dtype.names else ("rno" if "rno" in pk.dtype.names else pk.dtype.names[0])
def cellof(rc): return rc[1] if (np.ndim(rc) and len(rc) > 1) else (rc[-1] if np.ndim(rc) else rc)
row = pk[pk[RN] == R]
if len(row):
    row = row[0]
    rc_cell = cellof(row["cellid"])
    print(f"\n--- reach {R}: cellid={row['cellid']} -> cell {rc_cell} (x={xc[rc_cell]:.0f} y={yc[rc_cell]:.0f}); "
          f"{'SAME as 2422' if rc_cell==C else 'different from 2422'}")
    for f in ["rlen", "rwid", "rgrd", "rbth", "rhk", "ncon", "ndv"]:
        if f in pk.dtype.names: print(f"    {f} = {row[f]}")
cd = sfr.connectiondata.array
CN = "ifno" if "ifno" in cd.dtype.names else cd.dtype.names[0]
crow = [r for r in cd if r[CN] == R]
if crow:
    print(f"    reach {R} connections (ic): {crow[0]['ic']}")
on_C = [int(r[RN]) for r in pk if cellof(r["cellid"]) == C]
print(f"\n--- reaches sitting on cell {C}: {on_C}")

# MVR records feeding reach R (receiver = SFR, id R) and any provider/receiver touching cell C's reaches
mvr = gwf.get_package("MVR") or sim.get_package("MVR")
try:
    recs = mvr.perioddata.get_data()[0]
    print(f"\n--- MVR records with receiver SFR reach {R} (1-based id {R+1}):")
    n = 0
    for rec in recs:
        p1, i1, p2, i2 = rec[0], int(rec[1]), rec[2], int(rec[3])
        if str(p2).upper() == "SFR" and i2 == R + 1:
            print(f"    {p1}[{i1}] -> SFR[{i2}]  {rec[4]} {rec[5]}"); n += 1
    print(f"    ({n} providers feed reach {R})")
except Exception as e:
    print("MVR read failed:", e)
