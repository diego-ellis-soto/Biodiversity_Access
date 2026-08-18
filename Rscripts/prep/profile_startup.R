# ============================================================================
# Startup Performance Profiling
# ============================================================================
# This script benchmarks each data loading step in setup.R to identify
# bottlenecks and guide optimization decisions (precomputation, caching, etc.)

# library(tidyverse)
library(profvis)

# Create benchmark dataframe to track timings
benchmarks <- tibble(
  step = character(),
  description = character(),
  time_sec = numeric(),
  time_pct = numeric()
)

start_overall <- Sys.time()

# ============================================================================
# 1. Library Loading
# ============================================================================
t1 <- Sys.time()

library(shiny)
library(shinydashboard)
library(leaflet)
library(mapboxapi)
# library(tidyverse)
library(tidycensus)
library(sf)
library(DT)
library(RColorBrewer)
library(terra)
library(data.table)
library(mapview)
# library(sjPlot)
# library(sjlabelled)
library(bslib)
library(shinycssloaders)
library(DBI)
library(duckdb)
library(dbplyr)
library(gtfsrouter)
library(tidytransit)
library(fmsb)
library(scales)

time_libraries <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "00_libraries",
  description = "Load all 23 packages",
  time_sec = time_libraries,
  time_pct = NA
))

# ============================================================================
# 2. Greenspace (OSM)
# ============================================================================
t1 <- Sys.time()
osm_greenspace <- st_read("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/greenspaces_osm_nad83.shp", quiet = TRUE) |>
  st_transform(4326)
if (!"name" %in% names(osm_greenspace)) {
  osm_greenspace$name <- "Unnamed Greenspace"
}
time_greenspace <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "01_greenspace",
  description = "Load OSM greenspace polygons",
  time_sec = time_greenspace,
  time_pct = NA
))

# ============================================================================
# 3. Greenspace Distance Rasters
# ============================================================================
t1 <- Sys.time()
greenspace_dist_raster <- terra::rast("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/nearest_greenspace_dist.tif")
greenspace_osmid_raster <- terra::rast("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/nearest_greenspace_osmid.tif")
time_gs_rasters <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "02_gs_rasters",
  description = "Load greenspace distance rasters",
  time_sec = time_gs_rasters,
  time_pct = NA
))

# ============================================================================
# 4. NDVI Raster
# ============================================================================
t1 <- Sys.time()
ndvi <- terra::rast("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/SF_EastBay_NDVI_Sentinel_10.tif")
time_ndvi <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "03_ndvi",
  description = "Load NDVI raster",
  time_sec = time_ndvi,
  time_pct = NA
))

# ============================================================================
# 5. CBG Vector Data
# ============================================================================
t1 <- Sys.time()
download.file(
  'https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/cbg_vect_sf.Rdata',
  '/tmp/cbg_vect_sf.Rdata',
  mode = 'wb',
  quiet = TRUE
)
load('/tmp/cbg_vect_sf.Rdata')

if (!"unique_species" %in% names(cbg_vect_sf)) {
  cbg_vect_sf$unique_species <- cbg_vect_sf$n_species
}
if (!"n_observations" %in% names(cbg_vect_sf)) {
  cbg_vect_sf$n_observations <- cbg_vect_sf$n
}
if (!"median_inc" %in% names(cbg_vect_sf)) {
  cbg_vect_sf$median_inc <- cbg_vect_sf$medincE
}
if (!"ndvi_mean" %in% names(cbg_vect_sf)) {
  cbg_vect_sf$ndvi_mean <- cbg_vect_sf$ndvi_sentinel
}

time_cbg <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "04_cbg",
  description = "Download + load census block groups",
  time_sec = time_cbg,
  time_pct = NA
))

# ============================================================================
# 6. Hotspots/Coldspots
# ============================================================================
t1 <- Sys.time()
biodiv_hotspots <- st_read("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/hotspots.shp", quiet = TRUE) |>
  st_transform(4326)
biodiv_coldspots <- st_read("/vsicurl/https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access/resolve/main/coldspots.shp", quiet = TRUE) |>
  st_transform(4326)
time_hotcold <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "05_hotcold",
  description = "Load hotspots and coldspots",
  time_sec = time_hotcold,
  time_pct = NA
))

# ============================================================================
# 7. RSF Projects
# ============================================================================
t1 <- Sys.time()
rsf_projects <- st_read("data/source/RSF_Program_Projects_polygons.gpkg", quiet = TRUE) |>
  st_transform(4326)
time_rsf <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "06_rsf",
  description = "Load RSF program projects",
  time_sec = time_rsf,
  time_pct = NA
))

# ============================================================================
# 8. GBIF Parquet Reference (just get path, don't load yet)
# ============================================================================
time_gbif_setup <- 0.001  # negligible
benchmarks <- bind_rows(benchmarks, tibble(
  step = "07_gbif_setup",
  description = "GBIF parquet path (lazy load in server)",
  time_sec = time_gbif_setup,
  time_pct = NA
))

# ============================================================================
# 9. GTFS Data
# ============================================================================
t1 <- Sys.time()

gtfs_path <- '/Users/diegoellis/Desktop/RSF_next_steps/GPFS_OSM_Transit/sf_muni_gtfs-current/'

# Stops
gtfs_stops_sf <- tryCatch({
  read.csv(file.path(gtfs_path, 'stops.txt')) |>
    st_as_sf(coords = c("stop_lon", "stop_lat"), crs = 4326)
}, error = function(e) {
  warning("GTFS stops failed: ", e$message)
  NULL
})

# Route shapes
gtfs_shapes_raw <- read.csv(file.path(gtfs_path, 'shapes.txt'))
gtfs_trips_raw  <- read.csv(file.path(gtfs_path, 'trips.txt'))
gtfs_routes_raw <- read.csv(file.path(gtfs_path, 'routes.txt'))

shape_route_map <- gtfs_trips_raw |>
  distinct(shape_id, route_id)

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

gtfs_routes_sf <- st_sf(
  shape_id = sapply(shapes_split, function(s) s$shape_id[1]),
  geometry = st_sfc(shape_geoms, crs = 4326)
) |>
  left_join(shape_route_map, by = "shape_id") |>
  left_join(route_meta, by = "route_id")

time_gtfs_basic <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "08_gtfs_basic",
  description = "Load GTFS stops and route shapes",
  time_sec = time_gtfs_basic,
  time_pct = NA
))

# ============================================================================
# 10. gtfsrouter Initialization
# ============================================================================
t1 <- Sys.time()

gtfs_router <- tryCatch({
  gtfs_zip_path <- tempfile(fileext = ".zip")
  old_wd <- getwd()
  setwd(gtfs_path)
  utils::zip(gtfs_zip_path, files = list.files('.', pattern = "\\.txt$"))
  setwd(old_wd)

  gr <- gtfsrouter::extract_gtfs(gtfs_zip_path)
  gtfsrouter::gtfs_timetable(gr, day = "Monday")
}, error = function(e) {
  warning("gtfsrouter failed: ", e$message)
  NULL
})

time_gtfs_router <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "09_gtfs_router",
  description = "Initialize gtfsrouter + timetable",
  time_sec = time_gtfs_router,
  time_pct = NA
))

# ============================================================================
# 11. tidytransit Headways
# ============================================================================
t1 <- Sys.time()

gtfs_stop_headways <- tryCatch({
  gt <- tidytransit::read_gtfs(gtfs_path)
  tidytransit::get_stop_frequency(gt, start_time = 7 * 3600, end_time = 9 * 3600) |>
    group_by(stop_id) |>
    summarise(
      mean_headway_min  = mean(mean_headway, na.rm = TRUE) / 60,
      n_departures_peak = sum(n_departures, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(stop_id = as.character(stop_id))
}, error = function(e) {
  warning("tidytransit failed: ", e$message)
  NULL
})

if (!is.null(gtfs_stop_headways) && !is.null(gtfs_stops_sf)) {
  gtfs_stops_sf <- gtfs_stops_sf |>
    mutate(stop_id = as.character(stop_id)) |>
    left_join(gtfs_stop_headways, by = "stop_id")
}

time_gtfs_headways <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "10_gtfs_headways",
  description = "Compute transit stop headways (AM peak)",
  time_sec = time_gtfs_headways,
  time_pct = NA
))

# ============================================================================
# 12. CalEnviroScreen
# ============================================================================
t1 <- Sys.time()

calenviro_path <- '/Users/diegoellis/Downloads/calenviroscreen40gdb_F_2021.gdb'
if (!file.exists(calenviro_path)) {
  calenviro_path <- '/Users/diegoellis/Desktop/Projects/Presentations/Data_Schell_Lab_Tutorial/calenviroscreen40gdb_F_2021.gdb'
}

cenv_sf <- tryCatch({
  sf::st_read(calenviro_path, quiet = TRUE) |>
    dplyr::filter(grepl("san francisco", County, ignore.case = TRUE), !is.na(CIscore)) |>
    dplyr::select(
      Tract, CIscore, CIscoreP,
      PM2_5, PM2_5_Pctl, Traffic, Traffic_Pctl,
      Poverty, Poverty_Pctl, HousBurd, HousBurd_Pctl,
      County
    ) |>
    sf::st_transform(4326) |>
    sf::st_make_valid()
}, error = function(e) {
  warning("CalEnviroScreen failed: ", e$message)
  NULL
})

time_calenviro <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "11_calenviro",
  description = "Load CalEnviroScreen layer",
  time_sec = time_calenviro,
  time_pct = NA
))

# ============================================================================
# 13. SF EJ Communities
# ============================================================================
t1 <- Sys.time()

sf_ej_path <- '/Users/diegoellis/Downloads/San Francisco Environmental Justice Communities Map_20251217/geo_export_a21b0a0a-7306-46fd-8381-06581cdbe6e9.shp'

sf_ej_sf <- tryCatch({
  sf::st_read(sf_ej_path, quiet = TRUE) |>
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
    ) |>
    sf::st_transform(4326) |>
    sf::st_make_valid()
}, error = function(e) {
  warning("SF EJ layer failed: ", e$message)
  NULL
})

time_ej <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
benchmarks <- bind_rows(benchmarks, tibble(
  step = "12_ej",
  description = "Load SF EJ communities layer",
  time_sec = time_ej,
  time_pct = NA
))

# ============================================================================
# SUMMARY
# ============================================================================
time_overall <- as.numeric(difftime(Sys.time(), start_overall, units = "secs"))

benchmarks <- benchmarks |>
  mutate(time_pct = round(100 * time_sec / time_overall, 1))

# Print results
cat("\n")
cat("================================================================================\n")
cat("STARTUP PERFORMANCE BENCHMARK\n")
cat("================================================================================\n\n")

print(benchmarks |> select(step, description, time_sec, time_pct))

cat("\n")
cat("TOTAL STARTUP TIME: ", round(time_overall, 2), " seconds\n")
cat("================================================================================\n\n")

# ============================================================================
# OPTIMIZATION RECOMMENDATIONS
# ============================================================================
cat("OPTIMIZATION RECOMMENDATIONS:\n")
cat("================================================================================\n\n")

# Flag slow steps (>5 seconds or >10% of total)
slow_steps <- benchmarks |>
  filter(time_sec > 5 | time_pct > 10) |>
  arrange(desc(time_sec))

if (nrow(slow_steps) > 0) {
  cat("CRITICAL BOTTLENECKS (>5s or >10%):\n\n")
  for (i in 1:nrow(slow_steps)) {
    row <- slow_steps[i, ]
    cat("  •", row$step, "(", row$time_sec, "s,", row$time_pct, "%)\n")
    cat("    Description:", row$description, "\n")
    
    if (grepl("libraries", row$step)) {
      cat("    Recommendation: Load packages only in server() if possible. Use lazy loading.\n")
    } else if (grepl("download", row$step)) {
      cat("    Recommendation: Cache downloaded files locally or use precomputed versions.\n")
    } else if (grepl("gtfsrouter|tidytransit", row$step)) {
      cat("    Recommendation: Pre-compute and cache GTFS timetable. Consider lazy loading for session.\n")
    } else if (grepl("calenviro|ej", row$step)) {
      cat("    Recommendation: Pre-filter to SF boundary. Store as .gpkg or parquet locally.\n")
    }
    cat("\n")
  }
} else {
  cat("No critical bottlenecks detected (all steps < 5s).\n\n")
}

# Create visualization
p <- ggplot(benchmarks, aes(x = reorder(step, -time_sec), y = time_sec, fill = time_pct)) +
  geom_col() +
  geom_text(aes(label = paste0(round(time_sec, 2), "s")), vjust = -0.3, size = 3) +
  scale_fill_gradient(low = "green", high = "red", name = "% of Total") +
  labs(
    title = "Shiny App Startup Performance Profile",
    subtitle = paste0("Total time: ", round(time_overall, 2), " seconds"),
    x = "Loading Step",
    y = "Time (seconds)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

print(p)

ggsave("/Users/diegoellis/Desktop/Projects/Postdoc/Biodiversity_Access_Indicator/SF_biodiv_access_shiny/startup_benchmark.png",
       p, width = 12, height = 6, dpi = 150)

cat("\nBenchmark plot saved to: startup_benchmark.png\n")

# ============================================================================
# Export benchmark data
# ============================================================================
write_csv(benchmarks, 
          "/Users/diegoellis/Desktop/Projects/Postdoc/Biodiversity_Access_Indicator/SF_biodiv_access_shiny/startup_benchmarks.csv")

cat("Benchmark data saved to: startup_benchmarks.csv\n\n")
