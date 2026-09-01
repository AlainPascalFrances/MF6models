"""Dump the LAK budget terms for the last converged period — to see what drains
the perched lakes (evaporation vs leakage vs outflow)."""
import config
import numpy as np
import flopy

BUD = (str(config.MODEL) + r"\cdl_gwf.lak.bud")
cbb = flopy.utils.CellBudgetFile(BUD, precision="double")
recs = [r.strip() for r in cbb.get_unique_record_names(decode=True)]
kk = cbb.get_kstpkper()
last = kk[-1]
print(f"records: {recs}")
print(f"last converged kstpkper: {last}  (stress period {last[1]+1})\n")
print(f"{'term':<22}{'lake1(pond6)':>14}{'lake2(pond10)':>15}")
print("-" * 51)
for rn in recs:
    try:
        d = cbb.get_data(kstpkper=last, text=rn)[0]
        if hasattr(d, "dtype") and d.dtype.names and "q" in d.dtype.names:
            ids = d["node"] if "node" in d.dtype.names else np.arange(len(d)) + 1
            q = d["q"]
            l1 = float(q[ids == 1].sum()) if np.any(ids == 1) else float("nan")
            l2 = float(q[ids == 2].sum()) if np.any(ids == 2) else float("nan")
            print(f"{rn:<22}{l1:>14.3f}{l2:>15.3f}")
        else:
            print(f"{rn:<22}{str(np.ravel(d)[:3]):>29}")
    except Exception as e:
        print(f"{rn:<22}  ERR {e!r}")
print("\n(+ = into lake, - = out of lake)")
