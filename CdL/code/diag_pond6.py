"""Why does pond 6 fail at TS=1? Inspect cell 1840's layer stack + the lake stage vs the GW initial head there."""
import config
import numpy as np, flopy
from pathlib import Path
WS = Path(str(config.MODEL))
sim = flopy.mf6.MFSimulation.load(sim_ws=str(WS), verbosity_level=0)
gwf = sim.get_model()
disv = gwf.get_package("DISV"); top = disv.top.array; botm = disv.botm.array; idom = disv.idomain.array
npf = gwf.get_package("NPF"); k = npf.k.array; ict = npf.icelltype.array
ic = gwf.get_package("IC").strt.array
xc, yc = np.array(gwf.modelgrid.xcellcenters), np.array(gwf.modelgrid.ycellcenters)

for C in (1840, 2094, 2109):                                   # the failing cluster (DVMAX cells, 0-based)
    print(f"\ncell {C}: x={xc[C]:.0f} y={yc[C]:.0f}  top={top[C]:.2f}")
    for L in range(botm.shape[0]):
        thk = (top[C] if L == 0 else botm[L-1, C]) - botm[L, C]
        print(f"   L{L+1}: botm={botm[L,C]:7.2f} thk={thk:5.2f} idom={idom[L,C]:+d} "
              f"K={k[L,C]:5.1f} icelltype={ict[L,C]} | IC head={ic[L,C]:.2f}")

# lake 0 = FID 6: initial stage + which cells/layers it connects to
lak = gwf.get_package("LAK"); pdd = lak.packagedata.array; cd = lak.connectiondata.array
print("\nLAK packagedata (lakeno, strt, nconn, boundname):")
for r in pdd:
    print("  ", tuple(r))
CN = "ifno" if "ifno" in cd.dtype.names else cd.dtype.names[0]
l0 = cd[cd[CN] == cd[CN].min()]                                 # first lake's connections
print(f"\nlake {cd[CN].min()} connections: {len(l0)} ; connected (layer,cell):")
for r in l0[:6]:
    cell = int(r['cellid'][-1]); lay = int(r['cellid'][0])
    print(f"   -> L{lay+1} cell {cell}  bedleak={r['bedleak']}  (IC head {ic[lay,cell]:.2f}, lake strt {pdd[pdd[CN]==cd[CN].min()]['strt'][0]:.2f})")
print(f"\n=> head GAP at the lake bottom = lake_stage - GW_IC_head  (perched charca => lake sits ABOVE the deep WT)")
