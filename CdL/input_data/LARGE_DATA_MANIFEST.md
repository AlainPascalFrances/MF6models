# Large input data (not committed)

Two rasters exceed GitHub's 100 MB file limit and are **not** in this repo. Obtain them
separately, place them anywhere, and point the corresponding variable in
`code/cdl_gwf_model_fable_v2.py` at them.

| File | Size | Model variable | Purpose | Source |
|---|---|---|---|---|
| `dem_cdl.tif` | ~619 MB | `DEM_TIF` | LIDAR DEM (DGT), clipped to CdL — cell-top elevations & SFR routing | LNEG/DRYAD `GIS/Geodatabase_LIDAR_DGT/Geodatabase_CdL/` (DGT LIDAR) |
| `COSc_2025_N3_v0_TM06.tif` | ~187 MB | `COS_TIF` | COS 2025 N3 national land-cover (EPSG:3763) — per-cell rooting/extinction depths | DGT COS 2025 |

## Notes

- Both are **project/agency data** (LNEG, DGT). They are referenced, not redistributed here.
- If you only need to **run the transient** and you use the bundled `grid_cache/` (the Voronoi grid
  and conceptual `voronoi_layers.npz` are already built), the DEM is still read by the model for
  cell tops and SFR routing — supply it. To avoid the 187 MB national COS raster you may instead
  clip it to the CdL watershed once and repoint `COS_TIF` at the small clip; the derived
  `pest/cos_zones.npz` (already bundled) is what the PEST *forward run* uses.
- Placement suggestion (keeps prefix #3 consistent): drop them under your `<GIS_ROOT>` mirroring the
  original layout, or put them in `input_data/gis/` and set `DEM_TIF` / `COS_TIF` to those paths.
