import config
import geopandas as gpd, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import uniform_filter
from pathlib import Path
OUT=Path((str(config.PEST) + r"\snirh_data_availability"))
DEM=r"OpenFileGDB:E:\ArcGis_Data\WorkSpace\MDT_ESA\DEM_PI.gdb:DEM_EU_DEM_PI"
GPKG=(str(config.MODEL) + r"\gis\GIS\dryad_modelo_NbS.gpkg")
DEMCRS=3035
# piezo sets
pz=gpd.read_file(OUT/"snirh_piezo_T3T7_50km.gpkg").to_crs(DEMCRS)
obs=gpd.read_file(GPKG,layer="obs_points_cdl").to_crs(DEMCRS)
ws=gpd.read_file(GPKG,layer="watershed_cdl_fixed").to_crs(DEMCRS)
buf=ws.geometry.union_all().buffer(50_000)
minx,miny,maxx,maxy=buf.bounds; m=2000
with rasterio.open(DEM) as src:
    win=from_bounds(minx-m,miny-m,maxx+m,maxy+m,src.transform)
    z=src.read(1,window=win).astype("float32"); wt=src.window_transform(win); res=src.res[0]
z[z>1e30]=np.nan
# TPI (500 m neighbourhood) + slope
rad=int(round(500/res)); k=2*rad+1
zf=np.where(np.isnan(z),np.nanmean(z),z)
tpi=z-uniform_filter(zf,size=k,mode="nearest")
gy,gx=np.gradient(zf,res); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
tmu,tsd=np.nanmean(tpi),np.nanstd(tpi)
tpiz=(tpi-tmu)/tsd
def samp(gdf,namecol):
    rows=[]
    for _,r in gdf.iterrows():
        c=int((r.geometry.x-wt.c)/wt.a); rr=int((r.geometry.y-wt.f)/wt.e)
        if not (0<=rr<z.shape[0] and 0<=c<z.shape[1]): 
            rows.append(dict(id=r[namecol],elev=np.nan,tpi_z=np.nan,slope=np.nan,landform="out")); continue
        e,t,s=z[rr,c],tpiz[rr,c],slope[rr,c]
        lf=("hilltop" if t>1 else "valley" if t<-1 else ("flat/plain" if s<1.5 else "slope"))
        d=dict(id=str(r[namecol]),elev=round(float(e),1),tpi_z=round(float(t),2),slope=round(float(s),2),landform=lf)
        if "aquifer_system" in gdf.columns: d["aquifer"]="T7" if "ALUVI" in str(r["aquifer_system"]) else "T3"
        rows.append(d)
    return pd.DataFrame(rows)
sn=samp(pz,"code"); sn.to_csv(OUT/"snirh_piezo_geomorph.csv",index=False)
cd=samp(obs,"Name"); cd.to_csv(OUT/"cdl_piezo_geomorph.csv",index=False)
print("=== CdL P0-P6 geomorphology (EU-DEM, elevation basis) ===")
print(cd.to_string(index=False))
print("\n=== SNIRH T3/T7 landform distribution ===")
print(sn.groupby(["landform","aquifer"]).size().to_string())
print("elev range SNIRH:",round(sn.elev.min()),"-",round(sn.elev.max()),"m")
# ---- matching: each P -> SNIRH same landform, T7(shallow) preferred, nearest elevation ----
mrows=[]
for _,p in cd.iterrows():
    cand=sn[(sn.landform==p.landform)].copy()
    if len(cand)==0: cand=sn.copy()
    cand["elev_diff"]=(cand.elev-p.elev).abs()
    cand["t7"]=(cand.aquifer=="T7").astype(int)
    cand=cand.sort_values(["t7","elev_diff"],ascending=[False,True])
    top=cand.head(3)
    mrows.append(dict(cdl=p.id,cdl_elev=p.elev,cdl_landform=p.landform,
                      match1=f"{top.iloc[0].id}({top.iloc[0].aquifer},{top.iloc[0].elev}m,{top.iloc[0].landform})",
                      match2=f"{top.iloc[1].id}({top.iloc[1].aquifer},{top.iloc[1].elev}m)" if len(top)>1 else "",
                      match3=f"{top.iloc[2].id}({top.iloc[2].aquifer},{top.iloc[2].elev}m)" if len(top)>2 else ""))
mt=pd.DataFrame(mrows); mt.to_csv(OUT/"cdl_snirh_match.csv",index=False)
print("\n=== MATCHING (CdL piezo -> best SNIRH analogs: same landform, T7 preferred, nearest elevation) ===")
print(mt.to_string(index=False))

# --- landform map (task2_landform_map.png) ----------------------------------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
CRS=3763
_pz=pz.to_crs(CRS).merge(sn,left_on="code",right_on="id",how="left")
_obs=obs.to_crs(CRS).merge(cd,left_on="Name",right_on="id",how="left")
_ws=ws.to_crs(CRS)
COL={"hilltop":"firebrick","slope":"orange","flat/plain":"gold","valley":"tab:green"}
fig,ax=plt.subplots(figsize=(12,11))
gpd.GeoSeries([_ws.geometry.union_all()],crs=CRS).plot(ax=ax,facecolor="none",edgecolor="red",lw=2,zorder=8)
for lf,c in COL.items():
    s=_pz[_pz.landform==lf]
    if len(s): ax.scatter(s.geometry.x,s.geometry.y,c=c,s=30,edgecolor="0.3",lw=0.3,label=f"SNIRH {lf} ({len(s)})",zorder=4)
ax.scatter(_obs.geometry.x,_obs.geometry.y,marker="*",s=420,
           c=[COL.get(l,"grey") for l in _obs.landform],edgecolor="k",lw=1.2,zorder=9,label="CdL P0-P6")
for _,r in _obs.iterrows():
    ax.annotate(f"{r['Name']}\n{r.elev:.0f}m {r.landform}",(r.geometry.x,r.geometry.y),
                fontsize=8,fontweight="bold",xytext=(6,6),textcoords="offset points")
ax.set_title("Geomorphological position (EU-DEM 25 m): CdL piezometers vs SNIRH T3/T7\n"
             "landform by TPI(500 m)+slope | analysis in ELEVATION")
ax.set_xlabel("X (m, EPSG:3763)"); ax.set_ylabel("Y (m, EPSG:3763)"); ax.set_aspect("equal")
ax.legend(fontsize=8,loc="upper right"); ax.grid(alpha=0.3)
cx,cy=_ws.geometry.union_all().centroid.coords[0]; ax.set_xlim(cx-30000,cx+30000); ax.set_ylim(cy-30000,cy+30000)
fig.tight_layout(); fig.savefig(OUT/"task2_landform_map.png",dpi=150,bbox_inches="tight")
print("wrote task2_landform_map.png")
