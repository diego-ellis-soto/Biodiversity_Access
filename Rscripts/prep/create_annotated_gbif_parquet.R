# ============================================================================
# Prep: GBIF (DuckDB) → sf-gbif.parquet → CBG + NDVI → gbif_census_ndvi_anno.parquet
# ============================================================================
# 1) COPY San Francisco iNaturalist occurrences from attached gbif.duckdb → parquet
# 2) Join census block groups, extract NDVI, write annotated parquet for the app
#
# Requires: local gbif.duckdb path below; data/source/cbg_vect_sf.Rdata,
#           data/source/SF_EastBay_NDVI_Sentinel_10.tif

# library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(lubridate)
library(ggplot2)
library(duckdb)
library(sf)
library(terra)
library(arrow)

# ============================================================================
# 1) Export raw GBIF points
# ============================================================================

scon <- dbConnect(duckdb())
scon |> dbExecute("INSTALL spatial; LOAD spatial;")
scon |> dbExecute("ATTACH '~/Data/Occurrences/GBIF/gbif.duckdb' AS gbifmain;")

scon |> dbExecute(
  "COPY
(SELECT DISTINCT year,kingdom,phylum,class,\"order\",family,genus,species,decimallatitude,decimallongitude,coordinateuncertaintyinmeters,institutioncode FROM gbifmain.gbif
NATURAL JOIN gbifmain.join_lookup
WHERE stateprovince='California'
AND county_name='San Francisco'
AND institutioncode='iNaturalist'
AND decimallatitude IS NOT NULL
AND decimallongitude IS NOT NULL)
TO
'data/output/sf-gbif.parquet'
(FORMAT parquet);
"
)

scon |> dbDisconnect(shutdown = TRUE)

# ============================================================================
# 2) Annotate: CBG + NDVI → gbif_census_ndvi_anno.parquet
# ============================================================================

gbif_in <- "data/output/sf-gbif.parquet"
cbg_rdata <- "data/source/cbg_vect_sf.Rdata"
ndvi_tif <- "data/source/SF_EastBay_NDVI_Sentinel_10.tif"
out_dir <- "data/output"
out_file <- file.path(out_dir, "gbif_census_ndvi_anno.parquet")

gbif <- arrow::read_parquet(gbif_in, as_data_frame = TRUE)
load(cbg_rdata)
ndvi <- terra::rast(ndvi_tif)

cbg_acs <- cbg_vect_sf |>
  select(GEOID, medincE, popE, housingE)

pts <- gbif |>
  filter(!is.na(decimallongitude), !is.na(decimallatitude)) |>
  st_as_sf(coords = c("decimallongitude", "decimallatitude"), crs = 4326, remove = FALSE)

cbg_join <- st_join(pts, cbg_acs, join = st_within, left = TRUE)

cbg_join$ndvi_sentinel <- terra::extract(ndvi, cbg_join, ID = FALSE)$NDVI

out <- cbg_join |>
  mutate(
    geom_wkt = st_as_text(geometry),
    longitude = decimallongitude,
    latitude = decimallatitude
  ) |>
  st_drop_geometry() |>
  rename(
    institutionCode = institutioncode,
    coordinateUncertaintyInMeters = coordinateuncertaintyinmeters
  )

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(out, out_file)
