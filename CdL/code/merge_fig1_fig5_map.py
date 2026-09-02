"""
Merge Fig.1 (geological setting + SNIRH stations) and Fig.5 (upstream drainage area) of
CdL_SNIRH_data_import_methodology.docx into ONE map, with every legend item combined.  (user 2026-08-24)

Background geology is streamed LIVE from the LNEG geology WMS (fixes the white edges of the clipped
local raster).  The 1:500 000 map is NOT published as a WMS, and the 1:200 000 / 1:100 000 services
return blank at this regional zoom, so the reliable choice here is the national 1:1 000 000 map
(CGP1M).  The geological UNITS actually present in the view are detected via WMS GetFeatureInfo and
added to the legend with their real colours + names.  EU-Hydro streams are the second background
layer.  Requires internet.

Run in Spyder / the flopy env.  Writes the merged PNG next to the methodology document.
"""
import config
from pathlib import Path
import io, urllib.request, urllib.parse
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")                        # Spyder: comment out for an interactive window
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import geopandas as gpd
from shapely.geometry import Point

# ---------------- config ----------------
SNIRH = Path(str(config.PEST) + r"\snirh_data_availability")
MODEL_GPKG = Path(str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
OUT   = Path(str(config.PEST) + r"\CdL_geology_upstream_combined_map.png")

WMS_URL   = "https://inspire.lneg.pt/arcgis/services/CartografiaGeologica/CGP1M/MapServer/WmsServer"
WMS_LAYER = "PRT_LNEG_EN_1M_Lithology"        # renders the geology AND supports GetFeatureInfo
CRS = 3763                                    # WMS supports EPSG:3763 directly (project CRS)

EUHYDRO  = SNIRH / "euhydro_strahle_region.gpkg"
UPSTREAM = SNIRH / "inlet_upstream_basin.gpkg"
PIEZO    = SNIRH / "snirh_piezo_T3T7_50km.gpkg"
HIDRO    = SNIRH / "snirh_hidro_50km.gpkg"
PIEZO_ANALOGS = ["377/94", "378/129", "390/99", "391/243", "418/4"]
FLOW_STATIONS = ["20F/02H", "20E/02H", "21F/01H"]
INLET_XY = (-59956.4, -90732.1)              # SFR reach 293 (node 3506), EPSG:3763
MARGIN_FRAC = 0.04                             # symmetric map margin (the original tight bounding box)
WMS_WIDTH_PX = 2000                           # geology-image resolution (px along X)
GEO_UNIT_MINFRAC = 0.015                      # min area fraction for a unit to enter the legend


def _wms_bytes(params):
    url = WMS_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=90).read()


def fetch_geology(xmin, ymin, xmax, ymax, W, H):
    data = _wms_bytes(dict(service="WMS", version="1.1.1", request="GetMap", layers=WMS_LAYER, styles="",
                           srs=f"EPSG:{CRS}", bbox=f"{xmin},{ymin},{xmax},{ymax}", width=str(W), height=str(H),
                           format="image/png", transparent="false", bgcolor="0xFFFFFF"))
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"WMS did not return a PNG (got {data[:80]!r})")
    return plt.imread(io.BytesIO(data), format="png")


def detect_geology_units(img, xmin, ymin, xmax, ymax, W, H):
    """Return [(rgb01, label, frac), ...] for the geological units present in the view (>= MINFRAC),
    with their real colours and names read from the WMS via GetFeatureInfo."""
    M = (img[..., :3] * 255).round().astype(int)
    flat = M.reshape(-1, 3)
    uniq, cnt = np.unique(flat, axis=0, return_counts=True)
    common = dict(srs=f"EPSG:{CRS}", bbox=f"{xmin},{ymin},{xmax},{ymax}", width=str(W), height=str(H))
    units = []
    for i in np.argsort(-cnt):
        c = uniq[i]; frac = cnt[i] / len(flat)
        if frac < GEO_UNIT_MINFRAC:
            break                                        # sorted descending -> done
        if c.min() > 238:
            continue                                     # near-white background / rivers
        rows, cols = np.where((M[..., 0] == c[0]) & (M[..., 1] == c[1]) & (M[..., 2] == c[2]))
        k = len(rows) // 2
        try:
            d = _wms_bytes(dict(service="WMS", version="1.1.1", request="GetFeatureInfo", layers=WMS_LAYER,
                                query_layers=WMS_LAYER, info_format="application/json",
                                x=str(int(cols[k])), y=str(int(rows[k])), feature_count="1", **common))
            attr = {}
            for el in ET.fromstring(d).iter():
                if el.tag.endswith("FIELDS"):
                    attr = el.attrib; break
            name = attr.get("featureName") or attr.get("representativeLithology_label") or "unit"
            age = attr.get("representativeAge_label")
            label = f"{name} ({age})" if age and age not in ("Null", "None") else name
            units.append((tuple(np.array(c) / 255.0), label, frac, age))
        except Exception as e:
            print(f"   (unit info failed at RGB {tuple(int(v) for v in c)}: {e!r})")
    return units


# chronological order (youngest -> oldest) for the geology legend
AGE_ORDER = {"Holocene": 0, "Pleistocene": 1, "Pliocene": 2, "Miocene": 3, "Oligocene": 4,
             "Eocene": 5, "Paleocene": 6, "Cretaceous": 7, "Jurassic": 8, "Triassic": 9}


# ---------------- 1. vector layers (all -> EPSG:3763) ----------------
euhydro  = gpd.read_file(EUHYDRO).to_crs(CRS)
upstream = gpd.read_file(UPSTREAM).to_crs(CRS)
ws    = gpd.read_file(MODEL_GPKG, layer="watershed_cdl_fixed").to_crs(CRS)
piezo = gpd.read_file(PIEZO, layer="piezo").to_crs(CRS)
hidro = gpd.read_file(HIDRO, layer="hidro").to_crs(CRS)
inlet = gpd.GeoSeries([Point(*INLET_XY)], crs=CRS)
pz = piezo[piezo["code"].astype(str).isin(PIEZO_ANALOGS)].copy()
hy = hidro[hidro["code"].astype(str).isin(FLOW_STATIONS)].copy()
print(f">> piezo analogs {len(pz)}/{len(PIEZO_ANALOGS)}; flow stations {len(hy)}/{len(FLOW_STATIONS)}")

bounds = [ws.total_bounds, upstream.total_bounds, pz.total_bounds, hy.total_bounds]
xmin = min(b[0] for b in bounds); ymin = min(b[1] for b in bounds)
xmax = max(b[2] for b in bounds); ymax = max(b[3] for b in bounds)
rx = xmax - xmin; ry = ymax - ymin
xmin -= rx * MARGIN_FRAC; xmax += rx * MARGIN_FRAC
ymin -= ry * MARGIN_FRAC; ymax += ry * MARGIN_FRAC

# ---------------- 2. geology background + in-view units (LNEG WMS) ----------------
W = WMS_WIDTH_PX; H = max(1, int(round(W * (ymax - ymin) / (xmax - xmin))))
geo_img = fetch_geology(xmin, ymin, xmax, ymax, W, H)
print(f">> geology WMS image {geo_img.shape} over [{xmin:.0f},{ymin:.0f},{xmax:.0f},{ymax:.0f}]")
try:
    geo_units = detect_geology_units(geo_img, xmin, ymin, xmax, ymax, W, H)
    geo_units.sort(key=lambda u: AGE_ORDER.get(u[3], 99))          # Holocene -> Pleistocene -> Pliocene ...
    print(">> geology units in view:", [u[1] for u in geo_units])
except Exception as e:
    geo_units = []; print(f"   (geology-unit legend skipped: {e!r})")

# ---------------- 3. plot ----------------
fig, ax = plt.subplots(figsize=(12, 12 * (ymax - ymin) / (xmax - xmin)))
ax.imshow(geo_img, extent=[xmin, xmax, ymin, ymax], origin="upper", zorder=0)

euhydro.plot(ax=ax, color="tab:blue", lw=0.6, alpha=0.85, zorder=2)
upstream.plot(ax=ax, facecolor="orange", edgecolor="darkorange", alpha=0.35, lw=1.2, zorder=3)
ws.plot(ax=ax, facecolor="red", edgecolor="red", alpha=0.16, lw=2.2, zorder=4)
ws.boundary.plot(ax=ax, color="red", lw=2.2, zorder=5)
ax.scatter(pz.geometry.x, pz.geometry.y, s=95, marker="o", facecolor="lightblue",
           edgecolor="navy", linewidths=1.2, zorder=7)
for _, r in pz.iterrows():
    ax.annotate(str(r["code"]), (r.geometry.x, r.geometry.y), xytext=(6, 4),
                textcoords="offset points", fontsize=8, color="navy", weight="bold", zorder=8)
ax.scatter(hy.geometry.x, hy.geometry.y, s=180, marker="*", facecolor="royalblue",
           edgecolor="navy", linewidths=1.0, zorder=7)
for _, r in hy.iterrows():
    ax.annotate(str(r["code"]), (r.geometry.x, r.geometry.y), xytext=(7, 4),
                textcoords="offset points", fontsize=8, color="navy", weight="bold", zorder=8)
ax.scatter([inlet.iloc[0].x], [inlet.iloc[0].y], s=210, marker="v", facecolor="lime",
           edgecolor="black", linewidths=1.2, zorder=9)
ax.annotate("inlet (SFR reach 293)", (inlet.iloc[0].x, inlet.iloc[0].y), xytext=(8, -12),
            textcoords="offset points", fontsize=9, weight="bold", zorder=9)

ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_aspect("equal")

# ---------------- 4. two legends: monitoring network & features (top-left), geology (top-right) ----------------
feat_handles = [
    Line2D([0], [0], color="red", lw=2.2, label="CdL watershed / catchment (19.9 km\u00b2)"),
    Patch(facecolor="orange", edgecolor="darkorange", alpha=0.5, label="A_upstream \u2248 100 km\u00b2 (feeds inlet)"),
    Line2D([0], [0], color="tab:blue", lw=1.0, label="river network (EU-Hydro STRAHLE, EEA)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="lightblue", markeredgecolor="navy",
           markersize=10, label="SNIRH groundwater-level station"),
    Line2D([0], [0], marker="*", color="none", markerfacecolor="royalblue", markeredgecolor="navy",
           markersize=15, label="SNIRH streamflow station"),
    Line2D([0], [0], marker="v", color="none", markerfacecolor="lime", markeredgecolor="black",
           markersize=12, label="inlet (SFR reach 293)"),
]
# Legends INSIDE the map frame: monitoring & features -> bottom-left; geological units -> top-left.
leg1 = ax.legend(handles=feat_handles, loc="lower left", fontsize=8, framealpha=0.93,
                 title="Monitoring network & features")
leg1.get_title().set_fontweight("bold"); ax.add_artist(leg1)
if geo_units:
    geo_handles = [Patch(facecolor=rgb, edgecolor="0.35", label=label) for rgb, label, frac, age in geo_units]
    leg2 = ax.legend(handles=geo_handles, loc="upper left", fontsize=8, framealpha=0.93,
                     title="Geological units (LNEG CGP 1:1M)")
    leg2.get_title().set_fontweight("bold")

# 10 km scale bar (EPSG:3763 is metric), TOP-RIGHT
x0 = xmax - (xmax - xmin) * 0.06 - 10000; y0 = ymax - (ymax - ymin) * 0.06
ax.plot([x0, x0 + 10000], [y0, y0], color="black", lw=3, zorder=10)
ax.text(x0 + 5000, y0 - (ymax - ymin) * 0.018, "10 km", ha="center", va="top",
        fontsize=10, weight="bold", zorder=10)

ax.set_title("CdL \u2014 geological setting (LNEG CGP 1:1 000 000, WMS) + EU-Hydro river network "
             "+ SNIRH stations + upstream drainage area", fontsize=11)
ax.set_xlabel(f"X (m, EPSG:{CRS})"); ax.set_ylabel("Y (m)")
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f">> wrote {OUT}")
