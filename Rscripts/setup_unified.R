# ============================================================================
# Setup: load all data for the app (sourced once at startup by app.R)
# ============================================================================
# Checks local files in data/cached/ first, falling back to HuggingFace
# downloads (see hf_or_local()). Loads everything the app needs in one pass.
#
# The two most expensive products -- the CBG x greenspace intersection and the
# transit-routing timetable -- are precomputed in Rscripts/prep/ and just read
# here, so startup stays fast without splitting the load across files.
# ============================================================================

# ----------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------
# require(shinyjs)      # commented: useShinyjs() enabled features the app never calls
library(shiny)
library(shinydashboard)
library(leaflet)
library(mapboxapi)
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
# library(tidycensus)   # commented: unused in app (CBG data is precomputed) -- speeds startup
library(sf)
library(DT)
# library(RColorBrewer) # commented: brewer.pal only appears in dead if(FALSE) blocks
library(terra)
# Commented out: no detectable use in the app, and loading them added ~1-2 s to
# startup (mapview/sjPlot/sjlabelled each pull large dependency trees). Re-enable
# the relevant line if you start using one.
# library(data.table)
# library(mapview)
# library(sjPlot)
# library(sjlabelled)
library(bslib)
library(shinycssloaders)
library(glue)

# ----------------------------------------------------------------------------
# HuggingFace base URL + cache dir + download helper
# ----------------------------------------------------------------------------
HF_BASE <- "https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main"

# Use data/cached/ when running locally (writable), otherwise fall back to
# /tmp/sf_biodiv_cache/ for read-only environments like HuggingFace Spaces.
cache_dir <- if (file.access(".", mode = 2) == 0) "data/cached" else "/tmp/sf_biodiv_cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: if the file already exists in cache_dir, return that path. Otherwise
# attempt to download from HuggingFace into cache_dir. Returns the destination
# path regardless -- caller must check file.exists() if the download may fail.
hf_or_local <- function(filename) {
  dest <- file.path(cache_dir, filename)
  if (!file.exists(dest)) {
    tryCatch(
      download.file(glue::glue("{HF_BASE}/{filename}"), dest, mode = "wb", quiet = TRUE),
      error   = function(e) warning(glue::glue("HuggingFace download failed for {filename}: {e$message}")),
      warning = function(w) warning(glue::glue("HuggingFace download warning for {filename}: {w$message}"))
    )
  }
  dest
}

message("[setup] loading greenspace, CBG, RSF, CalEnviroScreen, SF EJ…")

# ----------------------------------------------------------------------------
# Greenspace (OSM polygons) -- "Greenspace" map layer + coverage calc
# ----------------------------------------------------------------------------
greenspace_shp <- file.path(cache_dir, "greenspaces_osm_nad83.shp")
if (!file.exists(greenspace_shp)) {
  for (ext in c("shp", "dbf", "prj", "shx")) {
    hf_or_local(glue::glue("greenspaces_osm_nad83.{ext}"))
  }
}
osm_greenspace <- st_read(greenspace_shp, quiet = TRUE) |> st_transform(4326)
if (!"name" %in% names(osm_greenspace)) osm_greenspace$name <- "Unnamed Greenspace"

# ----------------------------------------------------------------------------
# GBIF observations (parquet) -- path only; queried via DuckDB in app.R / server
# ----------------------------------------------------------------------------
gbif_parquet <- hf_or_local("gbif_census_ndvi_anno.parquet")

# ----------------------------------------------------------------------------
# Census block groups (CBG) -- Income / Richness / Data map layers
# ----------------------------------------------------------------------------
load(hf_or_local("cbg_vect_sf.Rdata"))

if (!"unique_species" %in% names(cbg_vect_sf))  cbg_vect_sf$unique_species  <- cbg_vect_sf$n_species
if (!"n_observations" %in% names(cbg_vect_sf)) cbg_vect_sf$n_observations <- cbg_vect_sf$n
if (!"median_inc"     %in% names(cbg_vect_sf)) cbg_vect_sf$median_inc     <- cbg_vect_sf$medincE
if (!"ndvi_mean"      %in% names(cbg_vect_sf)) cbg_vect_sf$ndvi_mean      <- cbg_vect_sf$ndvi_sentinel

# ----------------------------------------------------------------------------
# RSF Program Projects -- "RSF Program Projects" map layer
# ----------------------------------------------------------------------------
rsf_projects <- st_read(hf_or_local("RSF_Program_Projects_polygons.gpkg"), quiet = TRUE) |>
  st_transform(4326)

# ----------------------------------------------------------------------------
# CalEnviroScreen 4.0 (pre-filtered to SF) -- map overlay
# ----------------------------------------------------------------------------
cenv_sf <- tryCatch({
  sf::st_read(hf_or_local("calenviro_sf.gpkg"), quiet = TRUE)
}, error = function(e) {
  warning("CalEnviroScreen failed to load: ", e$message); NULL
})

# ----------------------------------------------------------------------------
# SF Environmental Justice Communities -- map overlay
# ----------------------------------------------------------------------------
sf_ej_sf <- tryCatch({
  sf::st_read(hf_or_local("sf_ej_communities_map.gpkg"), quiet = TRUE) |>
    dplyr::mutate(
      symbol_hex = stringr::str_split(symbol_rgb, ",\\s*") |>
        lapply(function(x) sprintf("#%02X%02X%02X",
                                   as.integer(x[1]), as.integer(x[2]), as.integer(x[3]))) |>
        unlist(),
      ej_label = dplyr::case_when(
        is.na(score) ~ "Not EJ",
        score >= 21  ~ "High EJ burden (21-30)",
        score >= 11  ~ "Moderate EJ burden (11-20)",
        score >= 1   ~ "Low EJ burden (1-10)",
        score == 0   ~ "Score 0",
        TRUE         ~ "Unknown"
      )
    )
}, error = function(e) {
  warning("SF EJ layer failed to load: ", e$message); NULL
})

message("[setup] loading distance rasters + NDVI (terra lazy: headers only)…")

# ----------------------------------------------------------------------------
# Greenspace + RSF distance/id rasters + NDVI (terra is lazy: reads headers, not
# cells, so these open instantly; pixels are read during isochrone analysis).
# ----------------------------------------------------------------------------
greenspace_dist_raster  <- terra::rast(hf_or_local("nearest_greenspace_dist.tif"))
greenspace_osmid_raster <- terra::rast(hf_or_local("nearest_greenspace_osmid.tif"))
rsfprogram_dist_raster  <- terra::rast(hf_or_local("nearest_rsfprogram_dist.tif"))
rsfprogram_id_raster    <- terra::rast(hf_or_local("nearest_rsfprogram_id.tif"))
ndvi <- terra::rast(hf_or_local("SF_EastBay_NDVI_Sentinel_10.tif"))

# ----------------------------------------------------------------------------
# Per-CBG greenspace overlap (used in the socio summary + BAI benchmark).
# Precomputed by Rscripts/prep/build_cbg_greenspace_coverage.R and cached as a
# small CSV (downloaded from HuggingFace like the other prep artifacts). If that
# file is unavailable, fall back to the (slow) st_union + st_intersection here.
# ----------------------------------------------------------------------------
cov_path <- hf_or_local("cbg_greenspace_coverage.csv")
cbg_greenspace_coverage <- if (file.exists(cov_path)) {
  # GEOID forced to character so leading zeros survive (FIPS codes; join key).
  readr::read_csv(cov_path, col_types = readr::cols(GEOID = readr::col_character()),
                  show_col_types = FALSE)
} else {
  message("[setup] precomputed greenspace coverage not found; computing inline (slow)…")
  cbg_proj <- st_transform(cbg_vect_sf[, "GEOID"], 3857) |>
    mutate(cbg_area_m2 = as.numeric(st_area(geometry)))
  gs_proj <- st_transform(osm_greenspace, 3857) |> st_make_valid()
  gs_union <- st_union(gs_proj)
  cbg_gs_inter <- st_intersection(cbg_proj, gs_union)
  coverage <- cbg_gs_inter |>
    mutate(greenspace_m2 = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(greenspace_m2 = sum(greenspace_m2), .groups = "drop") |>
    right_join(cbg_proj |> st_drop_geometry() |> dplyr::select(GEOID, cbg_area_m2), by = "GEOID") |>
    mutate(
      greenspace_m2 = tidyr::replace_na(greenspace_m2, 0),
      GEOID = as.character(GEOID)
    )
  readr::write_csv(coverage, file.path(cache_dir, "cbg_greenspace_coverage.csv"))
  coverage
}

message("[setup] loading GTFS (stops, routes, headways, routing timetable)…")

# ----------------------------------------------------------------------------
# GTFS (SF Muni)
# ----------------------------------------------------------------------------
gtfs_zip_path <- hf_or_local("sf_muni_gtfs.zip")

# Unzip for read.csv(stops.txt, …); gtfsrouter/tidytransit read the .zip directly.
gtfs_unzip_dir <- file.path(cache_dir, "muni_gtfs")
dir.create(gtfs_unzip_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(gtfs_unzip_dir) || length(list.files(gtfs_unzip_dir, pattern = "\\.txt$")) == 0L) {
  unzip(gtfs_zip_path, exdir = gtfs_unzip_dir, overwrite = TRUE)
}
gtfs_path <- gtfs_unzip_dir

# --- Transit stops ("Transit Stops" overlay) --------------------------------
gtfs_stops_sf <- tryCatch({
  read.csv(file.path(gtfs_path, "stops.txt")) |>
    st_as_sf(coords = c("stop_lon", "stop_lat"), crs = 4326)
}, error = function(e) { warning("GTFS stops failed to load: ", e$message); NULL })

# --- Route shapes ("Transit Routes" overlay) --------------------------------
gtfs_routes_sf <- tryCatch({
  gtfs_shapes_raw <- read.csv(file.path(gtfs_path, "shapes.txt"))
  gtfs_trips_raw  <- read.csv(file.path(gtfs_path, "trips.txt"))
  gtfs_routes_raw <- read.csv(file.path(gtfs_path, "routes.txt"))

  shape_route_map <- gtfs_trips_raw |> distinct(shape_id, route_id)
  route_meta <- gtfs_routes_raw |>
    select(route_id, route_short_name, route_long_name, route_color) |>
    mutate(route_color_hex = paste0("#", trimws(route_color)))

  shapes_split <- gtfs_shapes_raw |>
    arrange(shape_id, shape_pt_sequence) |>
    group_by(shape_id) |>
    group_split()

  shape_geoms <- lapply(shapes_split, function(s) {
    st_linestring(cbind(s$shape_pt_lon, s$shape_pt_lat))
  })

  st_sf(
    shape_id = sapply(shapes_split, function(s) s$shape_id[1]),
    geometry = st_sfc(shape_geoms, crs = 4326)
  ) |>
    left_join(shape_route_map, by = "shape_id") |>
    left_join(route_meta, by = "route_id")
}, error = function(e) { warning("GTFS route shapes failed to load: ", e$message); NULL })

# --- Stop headways (AM peak 7-9am): cached as CSV (readable / diffable).
# Joined into gtfs_stops_sf below so they show in the stop popups on the map.
hw_csv <- file.path(cache_dir, "gtfs_stop_headways.csv")
hw_rds <- file.path(cache_dir, "gtfs_stop_headways.rds")
if (!file.exists(hw_csv) && file.exists(hw_rds)) {
  readRDS(hw_rds) |> readr::write_csv(hw_csv)
}

gtfs_stop_headways <- tryCatch({
  headways_path <- hf_or_local("gtfs_stop_headways.csv")
  if (file.exists(headways_path)) {
    readr::read_csv(headways_path, show_col_types = FALSE) |>
      mutate(stop_id = as.character(stop_id))
  } else {
    gt <- tidytransit::read_gtfs(gtfs_zip_path)
    hw <- tidytransit::get_stop_frequency(gt, start_time = 7 * 3600, end_time = 9 * 3600) |>
      group_by(stop_id) |>
      summarise(
        mean_headway_min  = mean(mean_headway, na.rm = TRUE) / 60,
        n_departures_peak = sum(n_departures, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(stop_id = as.character(stop_id))
    readr::write_csv(hw, hw_csv)
    hw
  }
}, error = function(e) { warning("tidytransit headway computation failed: ", e$message); NULL })

if (!is.null(gtfs_stop_headways) && !is.null(gtfs_stops_sf)) {
  gtfs_stops_sf <- gtfs_stops_sf |>
    mutate(stop_id = as.character(stop_id)) |>
    left_join(gtfs_stop_headways, by = "stop_id")
}

# --- gtfsrouter timetable: routes transit isochrones (precomputed .rds) ------
gtfs_router <- tryCatch({
  timetable_path <- hf_or_local("gtfs_timetable_monday.rds")
  if (file.exists(timetable_path)) {
    readRDS(timetable_path)
  } else {
    gr <- gtfsrouter::extract_gtfs(gtfs_zip_path)
    result <- gtfsrouter::gtfs_timetable(gr, day = "Monday")
    saveRDS(result, file.path(cache_dir, "gtfs_timetable_monday.rds"))
    result
  }
}, error = function(e) { warning("gtfsrouter failed to initialise: ", e$message); NULL })

# --- Pre-computed transit isochrone cache (keyed by nearest stop + time) -----
transit_iso_cache <- tryCatch({
  p <- file.path(cache_dir, "transit_iso_cache.rds")
  if (file.exists(p)) readRDS(p) else NULL
}, error = function(e) { NULL })

message("[setup] data load complete.")
