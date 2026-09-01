"""
TASK 2 (data + plots) — RUN THIS IN SPYDER WHERE SNIRH IS REACHABLE.

The live SNIRH API is blocked from the assistant's environment, so this script
is prepared for you to run.  It:
  1. reads the offline geomorphological matching (built by task2_geomorph) that
     pairs each CdL piezometer P0-P6 with SNIRH T3/T7 analogs (same landform,
     T7=Aluviões do Tejo alluvial => shallow, nearest elevation);
  2. fetches the matched SNIRH piezometers' time series from SNIRH
     (PIEZOMETRIC_LEVEL = elevation; GWL_DEPTH = depth to water);
  3. builds SYNTHETIC hydrographs for P0-P6 by transferring the matched SNIRH
     station's seasonal fluctuation (anomaly about its own mean), anchored at
     each P's expected mean water level (ground elevation minus a nominal mean
     depth by landform);
  4. plots two time series: (a) ELEVATION and (b) DEPTH TO WATER TABLE.

Outputs -> E:\\00code_ws\\DRYAD\\CdL_pest\\snirh_data_availability\\
  snirh_piezo_series_raw.csv, cdl_synthetic_piezo.csv,
  piezo_synthetic_elevation.png, piezo_synthetic_depth.png
"""
import config
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import get_snirh as gs
from get_snirh.constants import Parameters

DIR   = Path((str(config.PEST) + r"\snirh_data_availability"))
START, END = "1981-01-01", "2026-05-01"
# nominal mean depth-to-water by landform (m) — anchors the synthetic mean level.
# Adjust to your field knowledge (qualitative piezos: ~1.5-3 m valley, ~10-15 m uplands).
MEAN_DEPTH = {"flat/plain": 3.0, "valley": 2.0, "slope": 8.0, "hilltop": 12.0}
SHALLOW_MAX_DEPTH = 50.0        # keep SNIRH analogs whose depth-to-water stays < this


def main():
    match = pd.read_csv(DIR / "cdl_snirh_match.csv")
    geom  = pd.read_csv(DIR / "cdl_piezo_geomorph.csv").set_index("id")
    pz    = gs.load_snapshot("piezometria")
    code2uid = dict(zip(pz.code.astype(str), pz.uid.astype(str)))

    # primary matched SNIRH code per CdL piezo (parse "code(aq,elev,landform)")
    match["snirh_code"] = match["match1"].str.split("(").str[0]
    match["snirh_uid"]  = match["snirh_code"].map(code2uid)
    stations = match[["snirh_uid", "snirh_code"]].dropna().drop_duplicates()
    stations.columns = ["uid", "code"]
    print(">> fetching", len(stations), "SNIRH piezometers:", list(stations.code))

    client = gs.SnirhClient()
    lvl = gs.fetch_timeseries(client, stations, Parameters.PIEZOMETRIC_LEVEL, START, END)
    dep = gs.fetch_timeseries(client, stations, Parameters.GWL_DEPTH, START, END)
    lvl.to_csv(DIR / "snirh_piezo_series_raw.csv", index=False)
    if lvl.empty:
        print("!! no PIEZOMETRIC_LEVEL returned — check SNIRH access / parameter."); return

    def series(df, code):
        s = df[df.code == code].copy()
        if s.empty:
            return None
        s["timestamp"] = pd.to_datetime(s["timestamp"])
        return s.set_index("timestamp")["value"].astype(float).sort_index()

    # keep only shallow analogs (depth-to-water stays < SHALLOW_MAX_DEPTH)
    ok = {}
    for code in stations.code:
        d = series(dep, code)
        if d is None or d.median() < SHALLOW_MAX_DEPTH:
            ok[code] = True
    fig1, ax1 = plt.subplots(figsize=(13, 5))     # ELEVATION
    fig2, ax2 = plt.subplots(figsize=(13, 5))     # DEPTH
    syn_rows = []
    cmap = plt.get_cmap("tab10")
    for i, (_, r) in enumerate(match.iterrows()):
        P, code = r["cdl"], r["snirh_code"]
        wl = series(lvl, code)
        if wl is None:
            print(f"   {P}: matched {code} returned no level series — skipped"); continue
        anom = wl - wl.mean()
        lf = geom.loc[P, "landform"]; g = float(geom.loc[P, "elev"])
        md = MEAN_DEPTH.get(lf, 5.0)
        P_elev = (g - md) + anom                    # synthetic water-level ELEVATION
        P_depth = g - P_elev                         # synthetic DEPTH to water
        col = cmap(i % 10)
        ax1.plot(wl.index, wl.values, color=col, lw=0.7, alpha=0.5)               # SNIRH source (elevation)
        ax1.plot(P_elev.index, P_elev.values, color=col, lw=1.8, label=f"{P} (from {code})")
        ax2.plot(P_depth.index, P_depth.values, color=col, lw=1.6, label=f"{P} (from {code})")
        for t, v in P_elev.items():
            syn_rows.append(dict(cdl=P, snirh=code, date=t, wl_elev_m=round(v, 3),
                                 depth_m=round(g - v, 3)))
    pd.DataFrame(syn_rows).to_csv(DIR / "cdl_synthetic_piezo.csv", index=False)

    ax1.set_ylabel("water-table ELEVATION (m)"); ax1.set_title(
        "CdL synthetic piezometers (bold) transferred from matched SNIRH T3/T7 analogs (thin) — ELEVATION")
    ax1.legend(fontsize=8, ncol=4); ax1.grid(alpha=0.3)
    ax2.set_ylabel("DEPTH to water table (m)"); ax2.invert_yaxis()
    ax2.set_title("CdL synthetic piezometers — DEPTH to water table"); ax2.legend(fontsize=8, ncol=4); ax2.grid(alpha=0.3)
    fig1.tight_layout(); fig1.savefig(DIR / "piezo_synthetic_elevation.png", dpi=150, bbox_inches="tight")
    fig2.tight_layout(); fig2.savefig(DIR / "piezo_synthetic_depth.png", dpi=150, bbox_inches="tight")
    print(">> wrote piezo_synthetic_elevation.png / piezo_synthetic_depth.png + cdl_synthetic_piezo.csv")


if __name__ == "__main__":
    main()
