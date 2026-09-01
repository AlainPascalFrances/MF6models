"""
TASK 3 (streamflow) — RUN THIS IN SPYDER WHERE SNIRH IS REACHABLE.

Prepared for you to run (the assistant's environment can't reach SNIRH).
  Part 1: import daily streamflow at 20E/02H PONTE SANTO ESTEVÃO (Rio Almansor,
          979.67 km2) from 1981 until its extinction (1990).
  Part 2: compare Sorraia (PONTE CORUCHE 20F/02H, 5847 km2), PONTE PRECES
          (19C/04H, 18.4 km2) and PONTE CARDOSAS (20C/03H, 61.2 km2): raw flow
          and specific discharge (normalised by drainage area) time series.

Streamflow is NOT in get_snirh's Parameters enum, so the flow parameter uid is
discovered live via fetch_parameters (label ~ 'caudal').

Outputs -> E:\\00code_ws\\DRYAD\\CdL_pest\\snirh_data_availability\\
  streamflow_20E02H_1981_1990.csv/.png, streamflow_comparison.csv/.png
"""
import config
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import get_snirh as gs

DIR = Path((str(config.PEST) + r"\snirh_data_availability"))
# code -> (uid, drainage_area_km2, river label)
STN = {
    "20E/02H": ("1627759236", 979.67, "Ponte Santo Estêvão (Almansor)"),
    "20F/02H": ("1627759070", 5847.16, "Ponte Coruche (Sorraia)"),
    "19C/04H": ("1627759184", 18.43,  "Ponte Preces (Sant. da Carnota)"),
    "20C/03H": ("1627759050", 61.18,  "Ponte Cardosas (Grande da Pipa)"),
}


def discover_flow_parameter(client, station_uid):
    """Find the streamflow ('caudal') parameter uid live (not in the enum)."""
    nets = gs.fetch_networks(client)
    ncol = next((c for c in nets.columns if nets[c].astype(str).str.contains("idrom", case=False).any()), None)
    hyd = nets[nets[ncol].astype(str).str.contains("idrom", case=False)] if ncol else nets
    uid_col = next(c for c in hyd.columns if "uid" in c.lower())
    net_uid = str(hyd.iloc[0][uid_col])
    pars = gs.fetch_parameters(client, net_uid, station_uid)
    lab_col = next((c for c in pars.columns if pars[c].astype(str).str.contains("caudal", case=False).any()), None)
    puid_col = next(c for c in pars.columns if "uid" in c.lower())
    hit = pars[pars[lab_col].astype(str).str.contains("caudal", case=False)] if lab_col else pars
    # prefer a daily mean flow ('médio diário') if present
    daily = hit[hit[lab_col].astype(str).str.contains("diári", case=False)] if lab_col else hit
    row = (daily if len(daily) else hit).iloc[0]
    print(f">> flow parameter: {row.to_dict()}")
    return str(row[puid_col])


def get_series(client, code, param_uid, start, end):
    uid, area, lab = STN[code]
    st = pd.DataFrame([{"uid": uid, "code": code}])
    ts = gs.fetch_timeseries(client, st, param_uid, start, end)
    if ts.empty:
        print(f"   {code}: no data"); return None
    ts["timestamp"] = pd.to_datetime(ts["timestamp"])
    return ts.set_index("timestamp")["value"].astype(float).sort_index()


def main():
    client = gs.SnirhClient()
    pflow = discover_flow_parameter(client, STN["20E/02H"][0])

    # --- Part 1: 20E/02H, 1981-1990 ---
    s = get_series(client, "20E/02H", pflow, "1981-01-01", "1990-12-31")
    if s is not None:
        s.to_frame("Q_m3s").to_csv(DIR / "streamflow_20E02H_1981_1990.csv")
        fig, ax = plt.subplots(figsize=(13, 4.5))
        ax.plot(s.index, s.values, color="tab:blue", lw=0.8)
        ax.set_ylabel("streamflow (m³/s)")
        ax.set_title("20E/02H Ponte Santo Estêvão (Rio Almansor, 979.7 km²) — daily flow 1981–1990")
        ax.grid(alpha=0.3); fig.tight_layout()
        fig.savefig(DIR / "streamflow_20E02H_1981_1990.png", dpi=150, bbox_inches="tight")
        print(">> wrote streamflow_20E02H_1981_1990.*")

    # --- Part 2: Sorraia / Preces / Cardosas — flow + specific discharge ---
    fig, (axq, axs) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    allrows = []
    for code in ["20F/02H", "19C/04H", "20C/03H"]:
        uid, area, lab = STN[code]
        q = get_series(client, code, pflow, "1980-01-01", "2026-05-01")
        if q is None:
            continue
        qspec = q * 86.4 / area                     # m³/s -> mm/day per km²
        axq.plot(q.index, q.values, lw=0.7, label=f"{code} {lab} ({area:.0f} km²)")
        axs.plot(qspec.index, qspec.values, lw=0.7, label=f"{code} ({area:.0f} km²)")
        for t, v in q.items():
            allrows.append(dict(code=code, area_km2=area, date=t, Q_m3s=v, q_mm_day=v * 86.4 / area))
    pd.DataFrame(allrows).to_csv(DIR / "streamflow_comparison.csv", index=False)
    axq.set_ylabel("flow Q (m³/s)"); axq.set_title("Streamflow comparison — raw discharge"); axq.legend(fontsize=8); axq.grid(alpha=0.3)
    axs.set_ylabel("specific discharge (mm/day)"); axs.set_title("Normalised by drainage area (specific discharge)")
    axs.legend(fontsize=8); axs.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(DIR / "streamflow_comparison.png", dpi=150, bbox_inches="tight")
    print(">> wrote streamflow_comparison.*")


if __name__ == "__main__":
    main()
