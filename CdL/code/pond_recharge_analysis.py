"""
================================================================================
CdL / DRYAD — OFFLINE POND-RECHARGE ANALYSIS  (no-lake deliverable, 2026-06-27)
================================================================================
The 19 charcas (cattle ponds) are a Nature-based Solution: they intercept hillslope
runoff that would otherwise leave the basin, pond it, and let it recharge the aquifer.
The LAK package is INTRACTABLE for these small ponds in MF6 (its lake-GWF connection
cell blew up at every configuration tried), so the model is shipped WITHOUT lakes and
the ponds' recharge contribution is quantified here, offline.

Method (focused-recharge / runoff-capture, after Daoud et al. 2022 reinfiltration):
    captured runoff to pond f  =  C_RUNOFF * P * A_catchment(f)
where A_catchment(f) is the pond's contributing area (D8 on the 20 m DEM, from
diag_pond_catchments.py -> pond_catchments.csv) and C_RUNOFF is the runoff coefficient.
The captured volume is the EXTRA recharge the pond delivers (the perched ponds in
particular convert would-be runoff into aquifer recharge).

Outputs (workspace):
    pond_recharge_summary.csv    per-pond table
    pond_recharge_analysis.png   bar chart (per pond) + seasonal climatology + totals
================================================================================
"""
import config
import csv
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- paths / parameters (mirror the main model) ----------------------------------
WS    = Path(str(config.MODEL))
GPKG  = (str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
P_CSV = (str(config.MODEL) + r"\gis\WP3_modeling\3.1.1\modelo_numerico\p_month_198101_202605.csv")
C_RUNOFF = 0.15        # runoff coefficient (fraction of catchment rainfall captured as runoff); calibratable starter value
RECHARGE_FRAC = 0.20   # natural aquifer recharge as a fraction of P (model SS value; for context only)
# Perched ponds (water table BELOW the pond bottom) — these convert captured runoff into NEW recharge.
# In-contact ponds (4,5,11,15,17,18) sit at/above the water table, so their capture is less "additional".
PERCHED = {0, 1, 2, 3, 6, 7, 8, 9, 10, 12, 13, 14, 16}
SUSPECT = {8}          # pond 8: catchment ratio ~1081 is SUSPECT (20 m DEM conflates its 10 m-offset footprint with the stream)

# ---- precipitation: monthly totals (mm) -> annual mean + monthly climatology -----
pdf = pd.read_csv(P_CSV, sep=r"\s+|,|;", engine="python")
_dcol = [c for c in pdf.columns if "date" in c.lower()][0]
pdf[_dcol] = pd.to_datetime(pdf[_dcol], format="%m/%d/%Y", errors="coerce")
pdf = pdf[[_dcol, "MONTHLY_TOTAL"]].dropna().rename(columns={_dcol: "date", "MONTHLY_TOTAL": "P_mm"})
pdf["year"], pdf["month"] = pdf["date"].dt.year, pdf["date"].dt.month
_ann = pdf.groupby("year")["P_mm"].sum()
_ann = _ann[(_ann.index >= 1981) & (_ann.index <= 2025)]   # complete calendar years only
MEAN_ANNUAL_P = float(_ann.mean())                          # mm/yr
MONTHLY_CLIM = pdf.groupby("month")["P_mm"].mean()          # mean mm per calendar month
print(f"precip: mean annual P = {MEAN_ANNUAL_P:.1f} mm/yr ({len(_ann)} complete years)")

# ---- pond catchments (D8) + pond geometry ----------------------------------------
_cat_csv = WS / "pond_catchments.csv"
if not _cat_csv.exists():
    raise FileNotFoundError(f"{_cat_csv} missing — run diag_pond_catchments.py first.")
catch = {int(r["fid"]): (float(r["pond_area_m2"]), float(r["catchment_area_m2"]))
         for r in csv.DictReader(open(_cat_csv))}

ws = gpd.read_file(GPKG, layer="watershed_cdl_fixed").to_crs(3763)
A_WS = float(ws.geometry.union_all().area)                 # watershed area, m2

# ---- per-pond captured runoff -> recharge ----------------------------------------
rows = []
for fid in sorted(catch):
    a_pond, a_catch = catch[fid]
    captured_m3yr = C_RUNOFF * (MEAN_ANNUAL_P / 1000.0) * a_catch          # m3/yr
    rows.append({
        "fid": fid,
        "type": "perched" if fid in PERCHED else "in-contact",
        "pond_area_m2": round(a_pond, 1),
        "catchment_ha": round(a_catch / 1e4, 2),
        "ratio_catch_pond": round(a_catch / a_pond, 1),
        "captured_m3_yr": round(captured_m3yr, 1),
        "depth_over_pond_mm_yr": round(captured_m3yr / a_pond * 1000.0, 1),  # mm/yr spread over the pond footprint
        "suspect": fid in SUSPECT,
    })
df = pd.DataFrame(rows).sort_values("captured_m3_yr", ascending=False)
df.to_csv(WS / "pond_recharge_summary.csv", index=False)

# ---- totals + context ------------------------------------------------------------
tot_all     = df["captured_m3_yr"].sum()
tot_perched = df.loc[df["type"] == "perched", "captured_m3_yr"].sum()
tot_robust  = df.loc[~df["suspect"], "captured_m3_yr"].sum()                # excl. the suspect pond 8
nat_recharge_m3yr = RECHARGE_FRAC * (MEAN_ANNUAL_P / 1000.0) * A_WS         # natural recharge over the whole watershed
catch_sum_ha = df["catchment_ha"].sum()
print(f"TOTAL captured runoff (all 19): {tot_all:,.0f} m3/yr  | perched-only: {tot_perched:,.0f} m3/yr"
      f"  | excl. suspect pond8: {tot_robust:,.0f} m3/yr")
print(f"natural recharge (~{RECHARGE_FRAC:.0%} of P over {A_WS/1e6:.1f} km2): {nat_recharge_m3yr:,.0f} m3/yr")
print(f"pond capture as % of natural recharge: all={100*tot_all/nat_recharge_m3yr:.1f}%  "
      f"robust={100*tot_robust/nat_recharge_m3yr:.1f}%")

# ---- figure ----------------------------------------------------------------------
fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[2, 1], hspace=0.32, wspace=0.22)

# (a) per-pond captured runoff (bar), coloured by type; hatch the suspect pond
axb = fig.add_subplot(gs[0, :])
dfx = df.sort_values("fid")
colors = ["#2c7fb8" if t == "perched" else "#d95f0e" for t in dfx["type"]]
bars = axb.bar([str(f) for f in dfx["fid"]], dfx["captured_m3_yr"], color=colors, edgecolor="k", lw=0.5)
for b, sus in zip(bars, dfx["suspect"]):
    if sus:
        b.set_hatch("////"); b.set_edgecolor("red")
axb.set_xlabel("pond FID"); axb.set_ylabel("captured runoff (m³/yr)")
axb.set_title(f"Per-pond captured runoff → recharge  (C_runoff={C_RUNOFF}, mean annual P={MEAN_ANNUAL_P:.0f} mm)")
from matplotlib.patches import Patch
axb.legend(handles=[Patch(facecolor="#2c7fb8", label="perched (adds new recharge)"),
                    Patch(facecolor="#d95f0e", label="in-contact (near water table)"),
                    Patch(facecolor="0.8", hatch="////", edgecolor="red", label="pond 8 — suspect catchment")],
           fontsize=8, loc="upper right")
axb.grid(axis="y", alpha=0.3)

# (b) seasonal: monthly captured runoff (total) vs monthly P climatology
axs = fig.add_subplot(gs[1, 0])
mon = np.arange(1, 13)
A_catch_total = sum(catch[f][1] for f in catch)
cap_month = C_RUNOFF * (MONTHLY_CLIM.reindex(mon).values / 1000.0) * A_catch_total   # m3/month, all ponds
axs.bar(mon, cap_month, color="#2c7fb8", alpha=0.85, label="captured runoff (all ponds)")
axs.set_xlabel("month"); axs.set_ylabel("captured runoff (m³/mo)"); axs.set_xticks(mon)
axp = axs.twinx(); axp.plot(mon, MONTHLY_CLIM.reindex(mon).values, "k-o", ms=4, lw=1.3, label="mean P")
axp.set_ylabel("precipitation (mm/mo)")
axs.set_title("Seasonality (monthly climatology)")
axs.grid(axis="y", alpha=0.3)

# (c) totals text panel
axt = fig.add_subplot(gs[1, 1]); axt.axis("off")
txt = (f"WATER BALANCE CONTEXT\n"
       f"────────────────────\n"
       f"Watershed area:        {A_WS/1e6:.1f} km²\n"
       f"Mean annual P:         {MEAN_ANNUAL_P:.0f} mm/yr\n"
       f"Natural recharge\n  (~{RECHARGE_FRAC:.0%} of P):        {nat_recharge_m3yr/1e6:.2f} hm³/yr\n\n"
       f"POND CAPTURE (C={C_RUNOFF})\n"
       f"────────────────────\n"
       f"All 19 ponds:          {tot_all/1e3:,.0f} ×10³ m³/yr\n"
       f"Perched only (13):     {tot_perched/1e3:,.0f} ×10³ m³/yr\n"
       f"Excl. suspect pond8:   {tot_robust/1e3:,.0f} ×10³ m³/yr\n\n"
       f"Capture vs natural\n"
       f"  recharge:    all {100*tot_all/nat_recharge_m3yr:.1f}%\n"
       f"               robust {100*tot_robust/nat_recharge_m3yr:.1f}%\n\n"
       f"⚠ C_runoff is a STARTER value —\n  calibrate against streamflow.\n"
       f"⚠ pond8 catchment SUSPECT.")
axt.text(0.0, 1.0, txt, va="top", ha="left", family="monospace", fontsize=9.5,
         bbox=dict(boxstyle="round", fc="#f7f7f7", ec="0.6"))

fig.suptitle("CdL / DRYAD — cattle-pond recharge contribution (offline; LAK retired as intractable in MF6)",
             fontsize=13, fontweight="bold")
fig.savefig(WS / "pond_recharge_analysis.png", dpi=150, bbox_inches="tight")
print(f"wrote {WS/'pond_recharge_summary.csv'}")
print(f"wrote {WS/'pond_recharge_analysis.png'}")
print("DONE")
