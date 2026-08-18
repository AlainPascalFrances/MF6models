"""Export a QGIS categorized .qml style for layer gc_35a_cdl, keyed on field 'Codigo',
using the official AML 100k geological-map symbology (user-provided RGB)."""
import os, geopandas as gpd

GIS  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS"
GPKG = os.path.join(GIS, "dryad_modelo_NbS.gpkg")
QML  = os.path.join(GIS, "geology_aml100k_symbology.qml")

# official AML 100k symbology colours (user-provided RGB), keyed on the AML Codigo
UNIT_COLORS = {"a": "#e1e1e1", "Qdae": "#d1d1c7", "Qt": "#bbbdb8", "Qi'": "#a6aba7",
               "PSA": "#ffffb7", "PSM": "#fff091", "MAT": "#ffff00"}
RECENT_TO_OLD = ["a", "Qdae", "Qt", "Qi'", "PSA", "PSM", "MAT"]    # legend order recent -> old

def rgba(hexc):
    h = hexc.lstrip("#")
    return f"{int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},255"

g = gpd.read_file(GPKG, layer="gc_35a_cdl")
present = list(dict.fromkeys(g["Codigo"].dropna()))
desc = g[["Codigo", "NomeOrigin"]].drop_duplicates().set_index("Codigo")["NomeOrigin"].to_dict()
codes = [c for c in RECENT_TO_OLD if c in present] + [c for c in present if c not in RECENT_TO_OLD]

cats, syms = [], []
for i, c in enumerate(codes):
    lab = f"{c} — {desc.get(c, '')}".strip(" —")
    cats.append(f'      <category value="{c}" label="{lab}" symbol="{i}" render="true"/>')
    syms.append(f'''      <symbol type="fill" name="{i}" alpha="1" clip_to_extent="1" force_rhr="0">
        <layer class="SimpleFill" enabled="1" locked="0" pass="0">
          <prop k="color" v="{rgba(UNIT_COLORS.get(c, '#dddddd'))}"/>
          <prop k="style" v="solid"/>
          <prop k="outline_color" v="50,50,50,255"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0.26"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_unit" v="MM"/>
        </layer>
      </symbol>''')

qml = f'''<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <renderer-v2 type="categorizedSymbol" attr="Codigo" forceraster="0" symbollevels="0" enableorderby="0" referencescale="-1">
    <categories>
{chr(10).join(cats)}
    </categories>
    <symbols>
{chr(10).join(syms)}
    </symbols>
  </renderer-v2>
  <blendMode>0</blendMode>
  <featureBlendMode>0</featureBlendMode>
  <layerGeometryType>2</layerGeometryType>
</qgis>
'''
with open(QML, "w", encoding="utf-8") as f:
    f.write(qml)
print("wrote", QML)
print("categories (recent->old):", ", ".join(codes))
