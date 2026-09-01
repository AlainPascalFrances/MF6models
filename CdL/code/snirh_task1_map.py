import config
import get_snirh as gs, geopandas as gpd, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely.geometry import Point
OUT=Path((str(config.PEST) + r"\snirh_data_availability")); OUT.mkdir(parents=True,exist_ok=True)
GPKG=(str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
CRS=3763
ws=gpd.read_file(GPKG,layer="watershed_cdl_fixed").to_crs(CRS)
ws_poly=ws.geometry.union_all(); buf=ws_poly.buffer(50_000)     # 50 km buffer
obs=gpd.read_file(GPKG,layer="obs_points_cdl").to_crs(CRS)
print("obs_points_cdl cols:", list(obs.columns), "| n", len(obs))
# SNIRH piezos -> gdf in 3763, filter T3/T7 aquifer + within buffer
def to_gdf(net):
    d=gs.load_snapshot(net)
    g=gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d.longitude,d.latitude), crs=4326).to_crs(CRS)
    return g
pz=to_gdf("piezometria"); hy=to_gdf("hidrometrica")
T37=pz["aquifer_system"].astype(str).str.contains("MARGEM ESQUERDA")|pz["aquifer_system"].astype(str).str.contains("ALUVI")
pzf=pz[T37 & pz.within(buf)].copy()
hyf=hy[hy.within(buf)].copy()
print(f"piezos T3/T7 within 50km: {len(pzf)}  | hidrometrica within 50km: {len(hyf)}")
pzf.to_file(OUT/"snirh_piezo_T3T7_50km.gpkg",layer="piezo",driver="GPKG")
hyf.to_file(OUT/"snirh_hidro_50km.gpkg",layer="hidro",driver="GPKG")
# ---- map ----
fig,ax=plt.subplots(figsize=(12,12))
gpd.GeoSeries([buf],crs=CRS).boundary.plot(ax=ax,color="0.6",ls="--",lw=1.2)
gpd.GeoSeries([ws_poly],crs=CRS).plot(ax=ax,facecolor="red",alpha=0.25,edgecolor="red",lw=2,zorder=5)
pzf.plot(ax=ax,color="tab:blue",markersize=26,zorder=4,label=f"piezometers T3/T7 ({len(pzf)})")
_act=hyf["status"].astype(str).str.upper().eq("ATIVA")
hyf[~_act].plot(ax=ax,color="0.6",marker="v",markersize=70,edgecolor="k",zorder=6,label=f"stream gauge inactive ({(~_act).sum()})")
hyf[_act].plot(ax=ax,color="tab:green",marker="v",markersize=90,edgecolor="k",zorder=7,label=f"stream gauge ATIVA ({_act.sum()})")
for _,r in pzf.iterrows(): ax.annotate(str(r["code"]),(r.geometry.x,r.geometry.y),fontsize=5.5,color="navy",xytext=(2,2),textcoords="offset points")
for _,r in hyf.iterrows(): ax.annotate(str(r["code"]),(r.geometry.x,r.geometry.y),fontsize=6.5,color="darkgreen" if str(r["status"]).upper()=="ATIVA" else "0.3",fontweight="bold",xytext=(3,3),textcoords="offset points")
ax.annotate("CdL",(ws_poly.centroid.x,ws_poly.centroid.y),fontsize=11,fontweight="bold",color="red",ha="center")
ax.set_title("SNIRH stations within 50 km of the CdL catchment\npiezometers in T3 (Tejo-Sado/Margem Esquerda) + T7 (Aluviões do Tejo)")
ax.set_xlabel("X (m, EPSG:3763)"); ax.set_ylabel("Y (m, EPSG:3763)"); ax.set_aspect("equal"); ax.legend(loc="upper right",fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT/"task1_snirh_T3T7_50km_map.png",dpi=150,bbox_inches="tight")
print("wrote task1_snirh_T3T7_50km_map.png")
