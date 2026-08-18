# ============================================================================
# Prep: Stage data/output → data/hf_upload (optional mirror before upload)
# ============================================================================
# Canonical build artifacts live in data/output/ (from Rscripts/prep/* and
# run_all_prep.R). Copy them here, then upload data/hf_upload/ to HuggingFace:
#
#   huggingface-cli upload boettiger-lab/sf_biodiv_access data/hf_upload/ . \
#     --repo-type dataset
#
# Or upload files from data/output/ directly via the web UI.
#
# Already on HuggingFace (not staged from data/output): greenspaces_osm_nad83.{shp,...}
# Other manual assets: SF_EastBay_NDVI_Sentinel_10.tif, cbg_vect_sf.Rdata, hotspots/coldspots
#
# GTFS: sf_muni_gtfs.zip contains the feed; timetable .rds + headways .csv are
# precomputed (~20–30 s) so the app does not rebuild them every session.
# ============================================================================

library(glue)

out_dir    <- "data/output"
upload_dir <- "data/hf_upload"
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

artifacts <- c(
  "nearest_greenspace_dist.tif",
  "nearest_greenspace_osmid.tif",
  "nearest_rsfprogram_dist.tif",
  "nearest_rsfprogram_id.tif",
  "gbif_census_ndvi_anno.parquet",
  "gtfs_timetable_monday.rds",
  "gtfs_stop_headways.csv",
  "sf_muni_gtfs.zip",
  "calenviro_sf.gpkg",
  "sf_ej_communities_map.gpkg",
  "RSF_Program_Projects_polygons.gpkg",
  "cbg_greenspace_coverage.csv"
)

for (f in artifacts) {
  src <- file.path(out_dir, f)
  dst <- file.path(upload_dir, f)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
  } else {
    warning(glue("Missing in data/output/: {f}"))
  }
}
