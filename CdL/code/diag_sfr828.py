"""Find the source of the ~126000 m3/d SFR throughflow: lake budget + reaches with the largest MVR/JA inflows."""
import numpy as np, flopy
from pathlib import Path
WS = Path(r"E:\00code_ws\DRYAD\CdL_model")

st = flopy.utils.HeadFile(WS / "cdl_gwf.lak.stage", text="STAGE")
print("lake stages (t=%.0f d): %s" % (st.get_times()[-1], st.get_data(totim=st.get_times()[-1]).ravel()))

# --- lake budget: bed seepage (GWF), outlet (TO-MVR), inflow (FROM-MVR), storage ---
lb = flopy.utils.CellBudgetFile(WS / "cdl_gwf.lak.bud", precision="double")
t = lb.get_times()[-1]
print("\nLAK budget records:", [n.strip().decode() for n in lb.get_unique_record_names()])
for nm in ["GWF", "TO-MVR", "FROM-MVR", "RAINFALL", "EVAPORATION", "STORAGE", "EXT-OUTFLOW"]:
    try:
        rec = lb.get_data(totim=t, text=nm)
        if rec:
            a = rec[0]
            q = a["q"] if hasattr(a, "dtype") and a.dtype.names and "q" in a.dtype.names else np.ravel(a)
            # aggregate per lake (GWF has many connections; sum by lake via 'node'? just show totals + per-lake)
            print(f"  {nm:12s}: sum={np.sum(q):12.1f}  vals(head)={np.array2string(np.atleast_1d(q)[:6], precision=1)}")
    except Exception as e:
        print(f"  {nm}: n/a ({e})")

# --- SFR: rank reaches by FROM-MVR and by |FLOW-JA-FACE| inflow ---
sb = flopy.utils.CellBudgetFile(WS / "cdl_gwf.sfr.bud", precision="double")
t = sb.get_times()[-1]
mvr = sb.get_data(totim=t, text="FROM-MVR")[0]
mvr_q = mvr["q"]; mvr_node = mvr["node"]
order = np.argsort(-np.abs(mvr_q))[:8]
print("\nTop SFR reaches by FROM-MVR (reach=node-1):")
for i in order:
    print(f"   reach {mvr_node[i]-1:5d} : FROM-MVR = {mvr_q[i]:12.2f} m3/d")

gwf = sb.get_data(totim=t, text="GWF")[0]
g_q = gwf["q"]; g_node = gwf["node"]
order = np.argsort(-np.abs(g_q))[:8]
print("\nTop SFR reaches by GWF exchange (negative = reach loses to GW):")
for i in order:
    print(f"   reach {g_node[i]-1:5d} : GWF = {g_q[i]:12.2f} m3/d")
