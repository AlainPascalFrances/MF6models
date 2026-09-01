"""Dump the LAK stage trajectory per lake per stress period (to diagnose a LAK
water-balance failure). Tries the binary STAGE file; falls back to parsing the
PRINT_STAGE tables in the GWF listing."""
import config
import re
import numpy as np
import flopy

WS = str(config.MODEL)
STG = WS + r"\cdl_gwf.lak.stage"
LST = WS + r"\cdl_gwf.lst"

# pond bottoms / rims for context (from diag_pond_perched / the model)
INFO = {1: ("pond6", 16.5, 20.0), 2: ("pond10", 30.7, 34.2)}   # lakeno -> (name, bottom, ~rim)

print("=== LAK stage trajectory ===")
ok = False
try:
    s = flopy.utils.HeadFile(STG, text="STAGE")
    times = s.get_times()
    arr = np.array([s.get_data(totim=t).ravel() for t in times])   # (ntime, nlakes)
    nlk = arr.shape[1]
    print(f"{'period':>6} " + " ".join(f"{INFO.get(k+1,('lk'+str(k+1),))[0]:>10}" for k in range(nlk)))
    for i, t in enumerate(times):
        print(f"{i+1:>6} " + " ".join(f"{arr[i,k]:>10.3f}" for k in range(nlk)))
    print("\nbottoms/rims:", {INFO[k][0]: (INFO[k][1], INFO[k][2]) for k in INFO})
    print("min/max per lake:", {INFO.get(k+1,('lk'+str(k+1),))[0]: (round(float(arr[:,k].min()),2),
                                  round(float(arr[:,k].max()),2)) for k in range(nlk)})
    ok = True
except Exception as e:
    print(f"(binary STAGE read failed: {e!r}; falling back to listing)")

if not ok:
    txt = open(LST, errors="ignore").read()
    for m in re.finditer(r"STAGE.*?LAKE.*?\n(.*?)\n\n", txt, re.S):
        print(m.group(0)[:400])
