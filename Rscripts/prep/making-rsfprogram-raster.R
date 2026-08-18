# ============================================================================
# Prep: Nearest RSF Program rasters (30 m, same grid as greenspace)
# ============================================================================
# Produces:
#   - nearest_rsfprogram_dist.tif — distance (m) to nearest RSF program polygon
#   - nearest_rsfprogram_id.tif   — polygon_id for lookup in RSF gpkg (prj_name)
#
# Source polygons: data/source/RSF_Program_Projects_polygons.gpkg
# Template grid: same slope.tif + SF county mask as making-greenspace-raster.R
# ============================================================================

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
library(terra)
library(tidyterra)
library(sf)
library(duckdb)
library(glue)
library(tictoc)
library(tidycensus)

# ============================================================================
# DuckDB + spatial
# ============================================================================
tcon <- dbConnect(duckdb::duckdb())

tcon |> dbExecute("
  INSTALL spatial;
  LOAD spatial;
  SET memory_limit = '200GB';
  SET preserve_insertion_order = false;
  SET threads TO 8;
")

rsf_gpkg <- "data/source/RSF_Program_Projects_polygons.gpkg"

# ============================================================================
# SF county boundary + template grid (match greenspace raster script)
# ============================================================================
sf <- get_acs(
  geography = "county",
  state = "CA",
  variables = "B01003_001",
  year = 2016,
  geometry = TRUE
) |>
  filter(NAME == "San Francisco County, California") |>
  st_transform(3310) |>
  st_cast("POLYGON") |>
  mutate(area = st_area(geometry)) |>
  slice_max(area, n = 1) |>
  select(GEOID, NAME, geometry)

empty.sr <- rast("data/source/slope.tif")

template.sr <- empty.sr %>%
  mutate(cell_id = 1:ncell(.)) %>%
  filter(!is.na(prcnt_slope30)) %>%
  crop(sf %>% st_transform(4326), mask = TRUE) %>%
  crop(ext(c(-123, -122, 37.65, 37.85))) %>%
  trim()

template_pts <- template.sr %>%
  as.points(na.rm = TRUE) %>%
  st_as_sf() %>%
  st_transform(4326)

dir.create("data/intermediate", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)

template_pts |>
  write_sf("data/intermediate/template_pts.gpkg")

# Escape single quotes in path for SQL (Windows paths unlikely here)
rsf_gpkg_sql <- gsub("'", "''", normalizePath(rsf_gpkg, winslash = "/", mustWork = TRUE))

tcon |> dbExecute(glue("
CREATE OR REPLACE TABLE rsf_geo AS
SELECT
  polygon_id,
  prj_name,
  ST_TRANSFORM(geom, 'EPSG:4326', 'EPSG:3310', always_xy := true) AS geom3310,
  ST_SimplifyPreserveTopology(
    ST_TRANSFORM(geom, 'EPSG:4326', 'EPSG:3310', always_xy := true),
    10
  ) AS simple_geom3310
FROM ST_Read('{rsf_gpkg_sql}');

CREATE INDEX idx_rsf_simple_geom ON rsf_geo USING RTREE (simple_geom3310);
"))

tcon |> dbExecute("
CREATE OR REPLACE TABLE template_grid_geo AS
SELECT
  cell_id,
  geom,
  ST_TRANSFORM(geom, 'EPSG:4326', 'EPSG:3310', always_xy := true) AS geom3310
FROM ST_Read('data/intermediate/template_pts.gpkg');

CREATE INDEX idx_template_geom ON template_grid_geo USING RTREE (geom3310);
CREATE INDEX idx_template_cell_id ON template_grid_geo (cell_id);
")

nn_query <- glue("
CREATE OR REPLACE TABLE rsf_distance_complete AS
WITH distances AS (
  SELECT
    template.cell_id,
    template.geom AS template_geom,
    ST_AsText(template.geom) AS geom_wkt,
    rsf.polygon_id,
    ST_Distance(template.geom3310, rsf.simple_geom3310) AS distance_meters
  FROM template_grid_geo AS template, rsf_geo AS rsf
)
SELECT
  cell_id,
  template_geom,
  geom_wkt,
  MIN(distance_meters) AS distance_to_rsf_meters,
  arg_min(polygon_id, distance_meters) AS nearest_rsf_polygon_id
FROM distances
GROUP BY cell_id, template_geom, geom_wkt;
")

cat("Calculating nearest RSF program distances...\n")
tic()
tcon |> dbExecute(nn_query)
toc()

rsf_dist_complete <- tcon |>
  tbl("rsf_distance_complete") |>
  select(cell_id, distance_to_rsf_meters, nearest_rsf_polygon_id) |>
  collect()

empty_raster <- empty.sr
values(empty_raster) <- NA

aligned_results <- tibble(cells = 1:ncell(empty.sr)) |>
  left_join(rsf_dist_complete |> rename(cells = cell_id), by = "cells")

rsf_nearest_dist <- empty_raster
values(rsf_nearest_dist) <- aligned_results$distance_to_rsf_meters
names(rsf_nearest_dist) <- "rsfprogram_nearest_dist"

cropped_rsf_dist <- rsf_nearest_dist %>%
  crop(sf %>% st_transform(4326), mask = TRUE) %>%
  crop(ext(c(-123, -122, 37.65, 37.85))) %>%
  trim()

rsf_nearest_id <- empty_raster
values(rsf_nearest_id) <- as.integer(round(as.numeric(aligned_results$nearest_rsf_polygon_id)))
names(rsf_nearest_id) <- "rsfprogram_nearest_id"

cropped_rsf_id <- rsf_nearest_id %>%
  crop(sf %>% st_transform(4326), mask = TRUE) %>%
  crop(ext(c(-123, -122, 37.65, 37.85))) %>%
  trim()

plot(cropped_rsf_dist)
plot(cropped_rsf_id)

writeRaster(cropped_rsf_dist, "data/output/nearest_rsfprogram_dist.tif", overwrite = TRUE)
writeRaster(cropped_rsf_id, "data/output/nearest_rsfprogram_id.tif", overwrite = TRUE, datatype = "INT4S")

dbDisconnect(tcon, shutdown = TRUE)
