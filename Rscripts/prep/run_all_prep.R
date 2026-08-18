# ============================================================================
# Prep: Run full data build → data/output
# ============================================================================
# Working directory must be the project root (folder containing data/ and Rscripts/).
#
# Order (fast steps first; heavy raster build last):
#   1. GBIF DuckDB → sf-gbif.parquet → gbif_census_ndvi_anno.parquet
#   2. CalEnviroScreen + SF EJ → GeoPackages
#   3. GTFS zip + timetable + stop headways (~20–30 s; cached for app startup)
#   4. Greenspace distance rasters (DuckDB; slow)
#   5. RSF Program nearest rasters (DuckDB; same grid as 4)
#   6. Per-CBG greenspace coverage CSV (precomputed so the app's first isochrone
#      doesn't recompute the CBG x greenspace intersection)
#
# Then upload files from data/output/ to the HuggingFace dataset (manually or
# see comments in upload_to_huggingface.R).

source("Rscripts/prep/create_annotated_gbif_parquet.R")
source("Rscripts/prep/build_equity_layers.R")
source("Rscripts/prep/implement_optimizations.R")
source("Rscripts/prep/making-greenspace-raster.R")
source("Rscripts/prep/making-rsfprogram-raster.R")
source("Rscripts/prep/build_cbg_greenspace_coverage.R")
