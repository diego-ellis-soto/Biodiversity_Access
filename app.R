# Currently downloaded the geomestry of San Francsico every time, pretty unefficient

###############################################################################
# Shiny App: San Francisco Biodiversity Access Decision Support Tool
# AUDITED citywide 100 m reference — path-safe app + conservation/displacement context UI
# Base app + multimodal transit + EJ + access scores + spider plot
#
# Scoring invariant:
# ONE user isochrone -> ONE raw metric row -> ONE same-mode × same-time
# citywide reference distribution. No pooling/averaging of user isochrones.
###############################################################################

# Next steps:
# Color code spiderploy by access-community partney-transportation and human movement-habitats (NDVI, Greenspace), species richness, sampling density of GBIF records (both are biodiversity),envrionmental justice (CalRnviro AND SF Screening Tool)
# Make spiderplot labels in boly and legend bigger
# Run the backgroun calculation for each isochrone
# WHAT to do about EJ layers in presidio and golden gate park, etc


# =============================================================================
# PACKAGES
# =============================================================================
# require(shinyjs)      # commented: useShinyjs() enabled features the app never calls
library(shiny)
library(shinydashboard)
library(leaflet)
library(leaflet.extras)
library(mapboxapi)
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
library(DBI)
library(duckdb)
library(dbplyr)
# Called only via pkg:: at runtime, so each namespace loads on first use -- no
# startup attach needed (keeps their dependency trees off startup, ~0.5-1 s):
#   gtfsrouter   -- transit routing (gtfsrouter::gtfs_route / gtfs_isochrone)
#   tidytransit  -- headway computation, fallback only (tidytransit::*)
#   fmsb         -- spider/radar plot (fmsb::radarchart)
# library(gtfsrouter)
# library(tidytransit)
# library(fmsb)
library(scales)

# =============================================================================
# PROJECT SETUP
# =============================================================================
# Resolve the project root automatically.
#
# This app can now be launched either as:
#   shiny::runApp("app_with_strict_single_active_corridor_isochrone_AUDITED_READY_UI_POLISHED_RADAR_CLEAR.R")
# or:
#   shiny::runApp("code/app_with_strict_single_active_corridor_isochrone_AUDITED_READY_UI_POLISHED_RADAR_CLEAR.R")
#
# All scientific inputs remain project-root-relative (code/, data/, outputs/).
.app_start_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (file.exists(file.path(.app_start_dir, "code", "setup_unified.R"))) {
  .project_root <- .app_start_dir
} else if (
  basename(.app_start_dir) == "code" &&
  file.exists(file.path(dirname(.app_start_dir), "code", "setup_unified.R"))
) {
  .project_root <- dirname(.app_start_dir)
} else {
  stop(
    "Could not locate the BAI_100m project root. Expected code/setup_unified.R ",
    "either in the current working directory or one directory above it."
  )
}

setwd(.project_root)

# Keep scientific data paths absolute for the lifetime of the Shiny process.
# Shiny may reset getwd() to the directory containing the app file after startup;
# absolute paths prevent lazy metric loaders from suddenly looking under code/.
project_path <- function(...) {
  normalizePath(
    file.path(.project_root, ...),
    winslash = "/",
    mustWork = FALSE
  )
}

.ISO_SPECIES_CELL_RDS <- project_path(
  "outputs", "results", "gbif_cell_species_100m.rds"
)
.ISO_ROUTE_CELL_RDS <- project_path(
  "outputs", "results", "muni_route_cells_100m.rds"
)

Sys.setenv(
  ISO_METRIC_GRID_GPKG = project_path(
    "outputs", "results", "grid_results_100m.gpkg"
  ),
  ISOCHRONE_REFERENCE_RDS = project_path(
    "outputs", "results", "isochrone_basi_100m.rds"
  )
)

# Fail early with an informative message rather than surfacing the same missing
# relative-path error from multiple Shiny outputs.
.required_harmonized_files <- c(
  grid = Sys.getenv("ISO_METRIC_GRID_GPKG"),
  reference = Sys.getenv("ISOCHRONE_REFERENCE_RDS"),
  species_lookup = .ISO_SPECIES_CELL_RDS,
  route_lookup = .ISO_ROUTE_CELL_RDS
)
.missing_harmonized_files <- .required_harmonized_files[
  !file.exists(.required_harmonized_files)
]
if (length(.missing_harmonized_files) > 0) {
  stop(
    "Required harmonized BAI files were not found:\n",
    paste0(" - ", .missing_harmonized_files, collapse = "\n"),
    "\nProject root resolved as: ", .project_root
  )
}

# setup_unified.R checks local files first, falling back to HuggingFace
# downloads, and loads everything the app needs at startup. The two expensive
# products (CBG x greenspace coverage, transit timetable) are precomputed in
# code/prep/ and just read here, so startup stays fast.
# Load the shared metric definition BEFORE the audited metric helper.
source("code/iso_metric_definitions.R")
source("code/setup_unified.R")

# =============================================================================
# GLOBAL CONFIG
# =============================================================================
mapbox_token <- Sys.getenv("MAPBOX_TOKEN")

if (!nzchar(mapbox_token)) {
  stop("MAPBOX_TOKEN environment variable is not set.")
}
theme <- bs_theme(
  bootswatch   = "minty",
  base_font    = font_google("Roboto"),
  heading_font = font_google("Roboto Slab"),
  bg           = "#f3f8f5",
  fg           = "#3d5c4a"
)

# =============================================================================
# OPTIONAL ENVIRONMENT / EQUITY LAYERS + GTFS
# =============================================================================
# cenv_sf, sf_ej_sf, gtfs_stops_sf, gtfs_routes_sf, gtfs_router,
# gtfs_zip_path, gtfs_stop_headways, transit_iso_cache are all loaded
# by setup_local.R (sourced above).

# =============================================================================
# GBIF UI VALUES FROM PARQUET
# =============================================================================
con_temp <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con_temp, "INSTALL spatial; LOAD spatial;")
dbExecute(con_temp, "INSTALL httpfs; LOAD httpfs;")
gbif_tab_temp <- tbl(con_temp, glue("read_parquet('{gbif_parquet}')"))

gbif_classes <- gbif_tab_temp |>
  distinct(class) |>
  collect() |>
  pull(class) |>
  sort()

gbif_families <- gbif_tab_temp |>
  distinct(family) |>
  collect() |>
  pull(family) |>
  sort()

dbDisconnect(con_temp, shutdown = TRUE)
rm(con_temp, gbif_tab_temp)

# =============================================================================
# HELPERS
# =============================================================================
pretty_mode <- function(x) {
  dplyr::case_when(
    x == "driving"         ~ "Driving",
    x == "walking"         ~ "Walking",
    x == "cycling"         ~ "Cycling",
    x == "driving-traffic" ~ "Driving-Traffic",
    x == "transit"         ~ "Transit",
    x == "walk_transit"    ~ "Walk-Transit",
    TRUE                   ~ tools::toTitleCase(x)
  )
}

mode_palette <- c(
  "Driving"         = "#4393C3",
  "Walking"         = "#74C476",
  "Cycling"         = "#FD8D3C",
  "Driving-Traffic" = "#9E9AC8",
  "Transit"         = "#D6604D",
  "Walk-Transit"    = "#E6AB02"
)

# Consistent conceptual colors used across radar labels, metric cards and
# explanatory legends. These colors describe BAI DIMENSION GROUPS, not
# transportation modes or comparison points.
BAI_GROUP_COLORS <- c(
  access = "#2f6fb0",
  biodiversity = "#c26a00",
  green_environment = "#2f855a",
  equity = "#7a58a6"
)

bai_dimension_group_legend_ui <- function() {
  tags$div(
    style = paste0(
      "display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); ",
      "gap:7px; margin:8px 0 10px 0;"
    ),
    tags$div(
      style = paste0(
        "border-left:5px solid ", BAI_GROUP_COLORS[["access"]],
        "; background:#f5f8fc; border-radius:7px; padding:7px 9px;"
      ),
      tags$b(style = paste0("color:", BAI_GROUP_COLORS[["access"]], ";"), "ACCESS"),
      tags$br(),
      tags$small("Transit access · Muni route access")
    ),
    tags$div(
      style = paste0(
        "border-left:5px solid ", BAI_GROUP_COLORS[["biodiversity"]],
        "; background:#fff8ef; border-radius:7px; padding:7px 9px;"
      ),
      tags$b(style = paste0("color:", BAI_GROUP_COLORS[["biodiversity"]], ";"), "BIODIVERSITY EVIDENCE"),
      tags$br(),
      tags$small("Species richness · GBIF sampling density")
    ),
    tags$div(
      style = paste0(
        "border-left:5px solid ", BAI_GROUP_COLORS[["green_environment"]],
        "; background:#f1f8f3; border-radius:7px; padding:7px 9px;"
      ),
      tags$b(style = paste0("color:", BAI_GROUP_COLORS[["green_environment"]], ";"), "GREEN ENVIRONMENT"),
      tags$br(),
      tags$small("Vegetation NDVI · Greenspace cover")
    ),
    tags$div(
      style = paste0(
        "border-left:5px solid ", BAI_GROUP_COLORS[["equity"]],
        "; background:#f7f3fb; border-radius:7px; padding:7px 9px;"
      ),
      tags$b(style = paste0("color:", BAI_GROUP_COLORS[["equity"]], ";"), "EQUITY CONTEXT"),
      tags$br(),
      tags$small("SF EJ context · outward = lower burden")
    )
  )
}


scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if ((rng[2] - rng[1]) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Muffle only sf's common overlay warning about non-geometry attributes being
# treated as spatially constant. Other warnings remain visible.
quiet_sf_overlay <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(
        "attribute variables are assumed to be spatially constant throughout all geometries",
        conditionMessage(w),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

ecdf01 <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  if (length(ref) == 0) return(rep(NA_real_, length(x)))
  f <- ecdf(ref)
  out <- f(x)
  out[is.na(x)] <- NA_real_
  out
}

standardize_iso_sf <- function(x, mode_name, time_min) {
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x <- st_as_sf(x)
  st_sf(
    mode = mode_name,
    time = as.numeric(time_min),
    geometry = st_geometry(x),
    crs = st_crs(x)
  )
}

choose_existing_sf_object <- function(candidates) {
  for (nm in candidates) {
    if (exists(nm, inherits = TRUE)) {
      obj <- get(nm, inherits = TRUE)
      if (inherits(obj, "sf")) return(obj)
    }
  }
  NULL
}

# Convert the FeatureCollection emitted by leaflet.extras into sf geometry.
# Polygon drawings represent habitat patches; line drawings represent proposed
# corridors and are buffered to the width selected in the scenario controls.
# A temporary GeoJSON file lets GDAL handle GeoJSON nesting without adding
# geojsonio as another app dependency.
draw_features_to_sf <- function(feature_collection) {
  if (is.null(feature_collection) ||
      is.null(feature_collection$features) ||
      length(feature_collection$features) == 0) {
    return(NULL)
  }
  
  tmp_geojson <- tempfile(fileext = ".geojson")
  on.exit(unlink(tmp_geojson), add = TRUE)
  
  json_txt <- jsonlite::toJSON(
    feature_collection,
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  )
  writeLines(json_txt, tmp_geojson, useBytes = TRUE)
  
  out <- tryCatch(
    suppressWarnings(st_read(tmp_geojson, quiet = TRUE)),
    error = function(e) NULL
  )
  if (is.null(out) || !inherits(out, "sf") || nrow(out) == 0) return(NULL)
  
  if (is.na(st_crs(out))) {
    st_crs(out) <- 4326
  } else {
    out <- st_transform(out, 4326)
  }
  
  out <- out |>
    st_zm(drop = TRUE, what = "ZM") |>
    st_make_valid()
  
  geom_type <- as.character(st_geometry_type(out))
  out <- out[
    geom_type %in% c("POLYGON", "MULTIPOLYGON", "LINESTRING", "MULTILINESTRING"),
    , drop = FALSE
  ]
  if (nrow(out) == 0) return(NULL)
  
  out
}

safe_vect_gbif_intersection <- function(poly_i) {
  out <- tryCatch({
    # browser()
    if (exists("vect_gbif")) {
      st_as_sf(intersect(vect_gbif, vect(poly_i)))
    } else if (exists("sf_gbif")) {
      st_intersection(sf_gbif, poly_i)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  out
}

safe_biodiv_hotspots <- function() {
  choose_existing_sf_object(c("biodiv_hotspots", "hotspots"))
}

safe_biodiv_coldspots <- function() {
  choose_existing_sf_object(c("biodiv_coldspots", "coldspots"))
}

safe_partner_orgs <- function() {
  choose_existing_sf_object(c(
    "partner_orgs_sf",
    "partner_organizations_sf",
    "community_orgs_sf",
    "community_organizations_sf",
    "partner_orgs",
    "community_orgs",
    "sf_partner_orgs",
    "community_sites_sf"
  ))
}

walk_speed_m_per_min <- 80

build_first_mile_walkshed <- function(location_sf, walk_minutes, mapbox_token) {
  if (is.null(walk_minutes) || walk_minutes <= 0) return(NULL)
  
  tryCatch(
    mb_isochrone(
      location_sf,
      time = walk_minutes,
      profile = "walking",
      access_token = mapbox_token
    ) |>
      st_as_sf() |>
      st_make_valid() |>
      st_transform(4326),
    error = function(e) NULL
  )
}

get_nearest_stops_by_distance <- function(location_sf, gtfs_stops_sf, walk_minutes, max_n = 5) {
  max_dist_m <- walk_minutes * walk_speed_m_per_min
  
  loc_proj   <- st_transform(location_sf, 3857)
  stops_proj <- st_transform(gtfs_stops_sf, 3857)
  
  dists <- as.numeric(st_distance(loc_proj, stops_proj))
  keep  <- which(dists <= max_dist_m)
  
  if (length(keep) == 0) {
    keep <- order(dists)[seq_len(min(max_n, length(dists)))]
  }
  
  gtfs_stops_sf[keep, ] |>
    mutate(
      dist_to_origin_m = round(dists[keep], 1),
      walk_time_to_origin_min = pmax(0, round(dists[keep] / walk_speed_m_per_min, 1))
    ) |>
    arrange(dist_to_origin_m)
}

get_walk_accessible_stops <- function(location_sf, walk_minutes, gtfs_stops_sf, mapbox_token) {
  if (is.null(walk_minutes) || walk_minutes <= 0) return(NULL)
  
  walk_iso <- build_first_mile_walkshed(location_sf, walk_minutes, mapbox_token)
  
  if (!is.null(walk_iso) && nrow(walk_iso) > 0) {
    stops_in_walkshed <- tryCatch(
      st_intersection(gtfs_stops_sf, st_union(walk_iso)),
      error = function(e) NULL
    )
    
    if (!is.null(stops_in_walkshed) && nrow(stops_in_walkshed) > 0) {
      dists <- as.numeric(st_distance(
        st_transform(location_sf, 3857),
        st_transform(stops_in_walkshed, 3857)
      ))
      
      return(
        stops_in_walkshed |>
          mutate(
            access_method = "mapbox_walkshed",
            dist_to_origin_m = round(dists, 1),
            walk_time_to_origin_min = pmax(0, round(dists / walk_speed_m_per_min, 1))
          ) |>
          arrange(walk_time_to_origin_min, dist_to_origin_m)
      )
    }
  }
  
  fallback_stops <- get_nearest_stops_by_distance(
    location_sf   = location_sf,
    gtfs_stops_sf = gtfs_stops_sf,
    walk_minutes  = walk_minutes,
    max_n         = 5
  )
  
  if (!is.null(fallback_stops) && nrow(fallback_stops) > 0) {
    fallback_stops <- fallback_stops |>
      mutate(access_method = "distance_fallback")
  }
  
  fallback_stops
}

extract_transit_minutes <- function(iso_result, dep_secs) {
  if (is.null(iso_result) || nrow(iso_result) == 0) return(numeric(0))
  
  if ("travel_time" %in% names(iso_result)) {
    return(as.numeric(iso_result$travel_time) / 60)
  }
  if ("duration" %in% names(iso_result)) {
    return(as.numeric(iso_result$duration) / 60)
  }
  if ("time" %in% names(iso_result)) {
    return(as.numeric(iso_result$time) / 60)
  }
  if ("arrival_time" %in% names(iso_result)) {
    return((as.numeric(iso_result$arrival_time) - dep_secs) / 60)
  }
  
  rep(NA_real_, nrow(iso_result))
}

build_last_mile_walkshed <- function(
    reachable_sf,
    remaining_walk_col = "remaining_walk_min",
    mapbox_token,
    walk_from_stop_cap_min = 8,
    max_stops = 12
) {
  if (is.null(reachable_sf) || nrow(reachable_sf) == 0) return(NULL)
  if (!(remaining_walk_col %in% names(reachable_sf))) return(NULL)
  
  rs <- reachable_sf |>
    mutate(
      remaining_walk_min = as.numeric(.data[[remaining_walk_col]]),
      remaining_walk_min = pmin(remaining_walk_min, walk_from_stop_cap_min)
    ) |>
    filter(is.finite(remaining_walk_min), remaining_walk_min > 0.5)
  
  if (nrow(rs) == 0) return(NULL)
  
  if ("n_departures_peak" %in% names(rs)) {
    rs <- rs |>
      arrange(desc(n_departures_peak), desc(remaining_walk_min))
  } else if ("mean_headway_min" %in% names(rs)) {
    rs <- rs |>
      arrange(mean_headway_min, desc(remaining_walk_min))
  } else {
    rs <- rs |>
      arrange(desc(remaining_walk_min))
  }
  
  rs <- rs |>
    slice_head(n = max_stops)
  
  walk_polys <- list()
  
  for (i in seq_len(nrow(rs))) {
    stop_i <- rs[i, ]
    walk_t <- floor(as.numeric(stop_i$remaining_walk_min[[1]]))
    if (!is.finite(walk_t) || walk_t <= 0) next
    
    iso_i <- tryCatch(
      mb_isochrone(
        stop_i,
        time = walk_t,
        profile = "walking",
        access_token = mapbox_token
      ),
      error = function(e) NULL
    )
    
    if (!is.null(iso_i) && nrow(iso_i) > 0) {
      walk_polys[[length(walk_polys) + 1]] <- st_as_sf(iso_i)
    }
  }
  
  if (length(walk_polys) == 0) return(NULL)
  
  walk_geom <- dplyr::bind_rows(walk_polys) |>
    st_as_sf() |>
    st_make_valid() |>
    st_union()
  
  st_sf(geometry = walk_geom, crs = 4326)
}

build_walk_transit_isochrone <- function(
    location_sf, total_time_min, dep_secs,
    walk_to_stop_min, walk_from_stop_min,
    gtfs_stops_sf, gtfs_router, mapbox_token,
    departure_window_min = 10,
    departure_step_min = 5,
    max_last_mile_stops = 12,
    include_first_mile_polygon = TRUE
) {
  if (is.null(gtfs_router)) return(NULL)
  if (is.null(total_time_min) || total_time_min <= 0) return(NULL)
  if (is.null(walk_to_stop_min) || walk_to_stop_min <= 0) return(NULL)
  
  first_mile_walkshed <- build_first_mile_walkshed(
    location_sf = location_sf,
    walk_minutes = walk_to_stop_min,
    mapbox_token = mapbox_token
  )
  
  origin_walk_stops <- get_walk_accessible_stops(
    location_sf   = location_sf,
    walk_minutes  = walk_to_stop_min,
    gtfs_stops_sf = gtfs_stops_sf,
    mapbox_token  = mapbox_token
  )
  
  if (is.null(origin_walk_stops) || nrow(origin_walk_stops) == 0) {
    return(NULL)
  }
  
  origin_walk_stops <- origin_walk_stops |>
    mutate(stop_id_chr = as.character(stop_id))
  
  departure_offsets_min <- seq(
    from = 0,
    to   = max(0, departure_window_min),
    by   = max(1, departure_step_min)
  )
  
  reachable_rows <- list()
  
  for (i in seq_len(nrow(origin_walk_stops))) {
    sid <- as.character(origin_walk_stops$stop_id_chr[[i]])
    first_mile_time_i <- as.numeric(origin_walk_stops$walk_time_to_origin_min[[i]])
    first_mile_time_i <- min(first_mile_time_i, walk_to_stop_min, na.rm = TRUE)
    
    if (!is.finite(first_mile_time_i) || first_mile_time_i >= total_time_min) next
    
    for (wait_offset_min in departure_offsets_min) {
      remaining_budget_before_transit <- total_time_min - first_mile_time_i - wait_offset_min
      if (!is.finite(remaining_budget_before_transit) || remaining_budget_before_transit <= 0) next
      
      start_time_i <- dep_secs + wait_offset_min * 60
      end_time_i   <- start_time_i + remaining_budget_before_transit * 60
      
      iso_result <- tryCatch(
        gtfsrouter::gtfs_isochrone(
          gtfs       = gtfs_router,
          from       = sid,
          start_time = start_time_i,
          end_time   = end_time_i,
          from_is_id = TRUE
        ),
        error = function(e) NULL
      )
      
      if (is.null(iso_result) || nrow(iso_result) == 0 || !("stop_id" %in% names(iso_result))) next
      
      transit_min <- extract_transit_minutes(iso_result, start_time_i)
      
      res_i <- iso_result |>
        mutate(
          stop_id_chr = as.character(stop_id),
          origin_stop_id = sid,
          first_mile_walk_min = first_mile_time_i,
          wait_time_min = wait_offset_min,
          transit_time_min = transit_min
        )
      
      if (!("transit_time_min" %in% names(res_i)) || all(is.na(res_i$transit_time_min))) {
        res_i$transit_time_min <- remaining_budget_before_transit
      }
      
      reachable_rows[[length(reachable_rows) + 1]] <- res_i
    }
  }
  
  if (length(reachable_rows) == 0) return(NULL)
  
  reachable_tbl <- dplyr::bind_rows(reachable_rows) |>
    mutate(total_pre_lastmile_min = first_mile_walk_min + wait_time_min + transit_time_min) |>
    group_by(stop_id_chr) |>
    summarise(
      first_mile_walk_min = min(first_mile_walk_min, na.rm = TRUE),
      wait_time_min       = min(wait_time_min, na.rm = TRUE),
      transit_time_min    = min(transit_time_min, na.rm = TRUE),
      best_total_pre_lastmile_min = min(total_pre_lastmile_min, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      remaining_walk_min = total_time_min - best_total_pre_lastmile_min,
      remaining_walk_min = pmin(remaining_walk_min, walk_from_stop_min)
    ) |>
    filter(is.finite(remaining_walk_min), remaining_walk_min > 0.5)
  
  if (nrow(reachable_tbl) == 0) return(NULL)
  
  reachable_sf <- gtfs_stops_sf |>
    mutate(stop_id_chr = as.character(stop_id)) |>
    inner_join(reachable_tbl, by = "stop_id_chr")
  
  if (nrow(reachable_sf) == 0) return(NULL)
  
  last_mile_walkshed <- build_last_mile_walkshed(
    reachable_sf = reachable_sf,
    remaining_walk_col = "remaining_walk_min",
    mapbox_token = mapbox_token,
    walk_from_stop_cap_min = walk_from_stop_min,
    max_stops = max_last_mile_stops
  )
  
  if (is.null(last_mile_walkshed) || nrow(last_mile_walkshed) == 0) return(NULL)
  
  final_geom <- st_union(st_geometry(last_mile_walkshed))
  
  if (include_first_mile_polygon && !is.null(first_mile_walkshed) && nrow(first_mile_walkshed) > 0) {
    final_geom <- st_union(final_geom, st_union(st_geometry(first_mile_walkshed)))
  }
  
  final_sf <- st_sf(geometry = final_geom, crs = 4326) |>
    st_make_valid()
  
  iso_sf <- st_sf(
    mode = "walk_transit",
    time = as.numeric(total_time_min),
    geometry = st_geometry(final_sf),
    crs = 4326
  )
  
  standardize_iso_sf(iso_sf, mode_name = "walk_transit", time_min = total_time_min)
}

# Shared isochrone -> metrics -> BAI -> radar functions. Both the Isochrone
# Explorer tab and the Isochrone Comparer tab call these, so the scoring logic
# lives in one place. (Reads the setup_unified.R globals + helpers above; the
# per-session gbif_tab is passed in as an argument; percentile reference scoring is added below.)
#
# local = TRUE is required: runApp() evaluates this app.R inside its own app
# environment, so mapbox_token and the helper functions above live there, not in
# the global env. Sourcing with the default (local = FALSE) would define these
# functions in the global env, where they could not see mapbox_token/pretty_mode/
# etc. -> "object 'mapbox_token' not found". local = TRUE defines them right here
# in the app env, alongside the helpers they call.
source("code/iso_metrics_AUDITED.R", local = TRUE)

# The audited helper is sourced into this app environment. Replace only its two
# lazy lookup loaders so they keep using the absolute project-root paths above
# even if Shiny later changes the working directory.
load_iso_species_cell_lookup <- function() {
  if (exists("species", envir = .iso_metric_cache, inherits = FALSE)) {
    return(get("species", envir = .iso_metric_cache, inherits = FALSE))
  }
  
  path <- .ISO_SPECIES_CELL_RDS
  if (!file.exists(path)) {
    warning(
      "Missing ", path, ". Run corrected 03_compute_isochrone_basi.R before ",
      "using definition-matched biodiversity percentiles."
    )
    result <- NULL
  } else {
    result <- readRDS(path)
    ver <- attr(result, "iso_metric_definition_version")
    if (is.null(ver) || !identical(ver, ISO_METRIC_DEFINITION_VERSION)) {
      stop("GBIF cell lookup definition version does not match the Shiny metric definition.")
    }
  }
  
  assign("species", result, envir = .iso_metric_cache)
  result
}

load_iso_route_cell_lookup <- function() {
  if (exists("routes", envir = .iso_metric_cache, inherits = FALSE)) {
    return(get("routes", envir = .iso_metric_cache, inherits = FALSE))
  }
  
  path <- .ISO_ROUTE_CELL_RDS
  if (!file.exists(path)) {
    warning(
      "Missing ", path, ". Run corrected 03_compute_isochrone_basi.R before ",
      "using definition-matched route-access percentiles."
    )
    result <- NULL
  } else {
    result <- readRDS(path)
    ver <- attr(result, "iso_metric_definition_version")
    if (is.null(ver) || !identical(ver, ISO_METRIC_DEFINITION_VERSION)) {
      stop("Muni route-cell lookup definition version does not match the Shiny metric definition.")
    }
  }
  
  assign("routes", result, envir = .iso_metric_cache)
  result
}


# =============================================================================
# PRECOMPUTED ISOCHRONE REFERENCE PERCENTILES — AUDITED
# =============================================================================
# ONE user-generated isochrone -> ONE independently calculated raw metric row ->
# ONE reference distribution with exactly the SAME mode and SAME time.
#
# The corrected Step-5 reference and compute_iso_metrics() share the definition
# version declared in code/iso_metric_definitions.R. Old reference files are
# rejected rather than silently mixing incompatible metric definitions.
iso_reference_path <- Sys.getenv(
  "ISOCHRONE_REFERENCE_RDS",
  unset = file.path("outputs", "results", "isochrone_basi_100m.rds")
)

normalize_reference_mode <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("_", "-", x, fixed = TRUE)
  x <- gsub("\\s+", "-", x)
  dplyr::case_when(
    x %in% c("walk", "walking") ~ "walking",
    x %in% c("bike", "bicycle", "cycling") ~ "cycling",
    x %in% c("drive", "driving") ~ "driving",
    x %in% c("driving-traffic", "driving-with-traffic") ~ "driving-traffic",
    x %in% c("walk-transit", "walk+transit", "walk-plus-transit") ~ "walk-transit",
    TRUE ~ x
  )
}

pick_first_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

# External-observation percentile convention.
# -------------------------------------------
# A user-generated isochrone is NOT one of the fixed reference observations.
# Therefore its percentile is the fraction of finite reference values it
# strictly outperforms, with denominator n (not n-1):
#   higher-is-better: 100 * #{ref < x} / n
#   lower-is-better : 100 * #{ref > x} / n
# This has an immediate interpretation: "the user value is better than X% of
# the matching precomputed reference isochrones." Ties are not counted as wins.
percent_rank_against_reference <- function(x, reference_values, direction = "higher") {
  x <- suppressWarnings(as.numeric(x))
  ref <- suppressWarnings(as.numeric(reference_values))
  ref <- ref[is.finite(ref)]
  if (length(x) != 1 || !is.finite(x) || length(ref) == 0) return(NA_real_)
  
  out <- if (identical(direction, "lower")) {
    100 * sum(ref > x) / length(ref)
  } else {
    100 * sum(ref < x) / length(ref)
  }
  pmin(100, pmax(0, out))
}

# Separate helper for QA only: reproduce the stored dplyr::percent_rank() of an
# observation already INCLUDED in the reference sample. Step 5 uses n-1 and
# min-rank tie behavior, and lower-better metrics are stored as 100 - percent_rank.
percent_rank_existing_reference <- function(x, reference_values, direction = "higher") {
  x <- suppressWarnings(as.numeric(x))
  ref <- suppressWarnings(as.numeric(reference_values))
  ref <- ref[is.finite(ref)]
  if (length(x) != 1 || !is.finite(x) || length(ref) < 2) return(NA_real_)
  
  direct <- 100 * sum(ref < x) / (length(ref) - 1)
  direct <- pmin(100, pmax(0, direct))
  if (identical(direction, "lower")) 100 - direct else direct
}

# Strict app raw -> corrected Step-5 raw mapping. There are intentionally no
# fallback aliases to older incompatible reference columns.
iso_reference_metric_spec <- list(
  GBIF_Species = list(
    raw = "n_species_accessed", pct = "pctile_n_species_accessed",
    out = "pctile_GBIF_Species", axis = "Biodiversity_Potential_std",
    label = "Total species richness", direction = "higher"
  ),
  Bird_Species = list(
    raw = "n_birds_accessed", pct = "pctile_n_birds_accessed",
    out = "pctile_Bird_Species", axis = NA_character_,
    label = "Bird species richness", direction = "higher"
  ),
  Mammal_Species = list(
    raw = "n_mammals_accessed", pct = "pctile_n_mammals_accessed",
    out = "pctile_Mammal_Species", axis = NA_character_,
    label = "Mammal species richness", direction = "higher"
  ),
  Plant_Species = list(
    raw = "n_plants_accessed", pct = "pctile_n_plants_accessed",
    out = "pctile_Plant_Species", axis = NA_character_,
    label = "Plant species richness", direction = "higher"
  ),
  GBIF_Records = list(
    raw = "n_records_accessed", pct = "pctile_n_records_accessed",
    out = "pctile_GBIF_Records", axis = NA_character_,
    label = "GBIF records", direction = "higher"
  ),
  SamplingDensity_km2 = list(
    raw = "sampling_density_km2", pct = "pctile_sampling_density_km2",
    out = "pctile_SamplingDensity_km2", axis = "Observation_Intensity_std",
    label = "GBIF sampling density", direction = "higher"
  ),
  MeanNDVI = list(
    raw = "mean_ndvi", pct = "pctile_mean_ndvi",
    out = "pctile_MeanNDVI", axis = "Environmental_Quality_std",
    label = "Mean NDVI", direction = "higher"
  ),
  Greenspace_percent = list(
    raw = "greenspace_pct", pct = "pctile_greenspace_pct",
    out = "pctile_Greenspace_percent", axis = "Greenspace_Cover_std",
    label = "Greenspace cover", direction = "higher"
  ),
  Transit_Stops = list(
    raw = "n_stops_accessed", pct = "pctile_n_stops_accessed",
    out = "pctile_Transit_Stops", axis = NA_character_,
    label = "Transit stops", direction = "higher"
  ),
  Transit_Access_Score = list(
    raw = "transit_access_score", pct = "pctile_transit_access_score",
    out = "pctile_Transit_Access_Score", axis = "Mobility_Access_std",
    label = "Transit stop density", direction = "higher"
  ),
  Unique_Muni_Routes = list(
    raw = "unique_muni_routes", pct = "pctile_unique_muni_routes",
    out = "pctile_Unique_Muni_Routes", axis = "Route_Access_std",
    label = "Muni route access", direction = "higher"
  ),
  Nearest_Stop_m = list(
    raw = "nearest_stop_m", pct = "pctile_nearest_stop_m",
    out = "pctile_Nearest_Stop_m", axis = NA_character_,
    label = "Nearest-stop distance", direction = "lower"
  ),
  CalEnviro_CIscore = list(
    raw = "calenviro_ci_score", pct = "pctile_calenviro_ci_score",
    out = "pctile_CalEnviro_CIscore", axis = NA_character_,
    label = "CalEnviroScreen burden", direction = "lower"
  ),
  CalEnviro_Traffic_Pctl = list(
    raw = "calenviro_traffic_pct", pct = "pctile_calenviro_traffic_pct",
    out = "pctile_CalEnviro_Traffic_Pctl", axis = NA_character_,
    label = "Traffic burden percentile", direction = "lower"
  ),
  SF_EJ_Score = list(
    raw = "sf_ej_score", pct = "pctile_sf_ej_score",
    out = "pctile_SF_EJ_Score", axis = "Equity_Context_std",
    label = "SF EJ context", direction = "lower"
  )
)

build_iso_reference_metric_map <- function(spec = iso_reference_metric_spec) {
  dplyr::bind_rows(lapply(names(spec), function(app_col) {
    s <- spec[[app_col]]
    data.frame(
      app_col = app_col,
      ref_raw_col = s$raw,
      ref_pct_col = s$pct,
      out_col = s$out,
      n_out_col = paste0("nref_", app_col),
      axis_col = s$axis,
      label = s$label,
      direction = s$direction,
      stringsAsFactors = FALSE
    )
  }))
}

load_iso_reference <- function(path = iso_reference_path) {
  if (!file.exists(path)) {
    warning(
      "Precomputed isochrone reference not found: ", path,
      ". Raw metrics can still be calculated, but same-mode × same-time ",
      "percentiles are unavailable."
    )
    return(NULL)
  }
  
  ref <- tryCatch(readRDS(path), error = function(e) {
    warning("Could not read precomputed isochrone reference: ", e$message)
    NULL
  })
  if (is.null(ref)) return(NULL)
  if (inherits(ref, "sf")) ref <- sf::st_drop_geometry(ref)
  ref <- as.data.frame(ref)
  
  required_identity <- c("mode", "time_min", "reference_schema_version")
  if (!all(required_identity %in% names(ref))) {
    warning(
      "Isochrone reference is an OLD/INCOMPATIBLE schema. Missing: ",
      paste(setdiff(required_identity, names(ref)), collapse = ", "),
      ". Rerun 03_compute_isochrone_basi_AUDITED.R (or the validated recovery output); do not use this RDS for scoring."
    )
    return(NULL)
  }
  
  versions <- unique(na.omit(as.character(ref$reference_schema_version)))
  if (length(versions) != 1 || !identical(versions[[1]], ISO_METRIC_DEFINITION_VERSION)) {
    warning(
      "Isochrone reference metric-definition version mismatch. Expected '",
      ISO_METRIC_DEFINITION_VERSION, "' but found '", paste(versions, collapse = ", "),
      "'. Rerun corrected Step 5."
    )
    return(NULL)
  }
  
  ref$mode <- normalize_reference_mode(ref$mode)
  ref$time_min <- suppressWarnings(as.numeric(ref$time_min))
  
  metric_map <- build_iso_reference_metric_map()
  required_raw <- unique(metric_map$ref_raw_col)
  missing_raw <- setdiff(required_raw, names(ref))
  if (length(missing_raw) > 0) {
    warning(
      "Corrected isochrone reference is missing required raw metrics: ",
      paste(missing_raw, collapse = ", "),
      ". Percentile scoring is disabled rather than substituting another metric."
    )
    return(NULL)
  }
  
  counts <- ref |>
    dplyr::filter(is.finite(time_min)) |>
    dplyr::count(mode, time_min, name = "n_reference") |>
    dplyr::arrange(mode, time_min)
  
  # This audited reference was validated against the complete SF 100 m origin
  # grid: 24,182 origins in each supported mode × time group.
  expected_counts <- tidyr::crossing(
    mode = c("walking", "cycling", "driving"),
    time_min = c(5, 10, 15)
  ) |>
    dplyr::left_join(counts, by = c("mode", "time_min"))
  
  bad_counts <- is.na(expected_counts$n_reference) |
    expected_counts$n_reference != 24182L
  
  if (any(bad_counts)) {
    warning(
      "Isochrone reference is incomplete. Expected 24,182 origin rows in every ",
      "walking/cycling/driving × 5/10/15-minute group. Scoring is disabled. ",
      "Rerun the validated Step-5 recovery/output before launching the app."
    )
    return(NULL)
  }
  
  message("Loaded corrected same-mode × same-time isochrone reference: ", path)
  if (nrow(counts) > 0) {
    message(paste0(
      "Reference groups: ",
      paste(paste0(counts$mode, " ", counts$time_min, " min n=", counts$n_reference), collapse = "; ")
    ))
  }
  ref
}

iso_reference_metric_map <- build_iso_reference_metric_map()
iso_reference <- load_iso_reference()

# Cache each exact mode × time reference subset. No other group may enter a row's
# percentile calculation.
.iso_reference_group_cache <- new.env(parent = emptyenv())
get_iso_reference_group <- function(mode, time_min, ref = iso_reference) {
  if (is.null(ref)) return(NULL)
  mode <- normalize_reference_mode(mode)
  time_min <- suppressWarnings(as.numeric(time_min))
  if (length(mode) != 1 || is.na(mode) || length(time_min) != 1 || !is.finite(time_min)) return(NULL)
  
  if (!(mode %in% c("walking", "cycling", "driving")) || !(time_min %in% c(5, 10, 15))) {
    return(NULL)
  }
  
  key <- paste0(mode, "__", time_min)
  if (exists(key, envir = .iso_reference_group_cache, inherits = FALSE)) {
    return(get(key, envir = .iso_reference_group_cache, inherits = FALSE))
  }
  
  g <- ref[
    ref$mode == mode & is.finite(ref$time_min) & abs(ref$time_min - time_min) < 1e-9,
    , drop = FALSE
  ]
  if (nrow(g) > 0) {
    stopifnot(
      dplyr::n_distinct(g$mode) == 1L,
      dplyr::n_distinct(g$time_min) == 1L,
      all(g$mode == mode),
      all(abs(g$time_min - time_min) < 1e-9)
    )
  }
  assign(key, g, envir = .iso_reference_group_cache)
  g
}

# Same seven-axis BAI distribution for the fixed citywide reference group.
#
# IMPORTANT:
# - one mode × one travel time only;
# - no isochrone groups are pooled;
# - only reference rows with all seven stored percentile dimensions contribute;
# - the USER BAI is treated as an external observation and ranked using the
#   same strict-outperformance convention as every other live percentile.
.reference_bai_pct_cols <- c(
  "pctile_transit_access_score",
  "pctile_unique_muni_routes",
  "pctile_n_species_accessed",
  "pctile_sampling_density_km2",
  "pctile_mean_ndvi",
  "pctile_greenspace_pct",
  "pctile_sf_ej_score"
)

.iso_reference_bai_cache <- new.env(parent = emptyenv())

get_reference_bai_distribution <- function(mode, time_min, ref = iso_reference) {
  mode <- normalize_reference_mode(mode)
  time_min <- suppressWarnings(as.numeric(time_min))
  
  if (length(mode) != 1 || is.na(mode) ||
      length(time_min) != 1 || !is.finite(time_min)) {
    return(numeric(0))
  }
  
  if (!(mode %in% c("walking", "cycling", "driving")) ||
      !(time_min %in% c(5, 10, 15))) {
    return(numeric(0))
  }
  
  key <- paste0("bai__", mode, "__", time_min)
  if (exists(key, envir = .iso_reference_bai_cache, inherits = FALSE)) {
    return(get(key, envir = .iso_reference_bai_cache, inherits = FALSE))
  }
  
  g <- get_iso_reference_group(mode, time_min, ref = ref)
  result <- numeric(0)
  
  if (!is.null(g) && nrow(g) > 0 &&
      all(.reference_bai_pct_cols %in% names(g))) {
    mat <- as.matrix(g[, .reference_bai_pct_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    complete <- rowSums(is.finite(mat)) == length(.reference_bai_pct_cols)
    
    if (any(complete)) {
      # Stored reference percentiles are on 0–100; BAI itself is on 0–1.
      result <- rowMeans(mat[complete, , drop = FALSE]) / 100
      result <- result[is.finite(result)]
    }
  }
  
  assign(key, result, envir = .iso_reference_bai_cache)
  result
}


validate_iso_reference_percentiles <- function(
    ref = iso_reference,
    metric_map = iso_reference_metric_map,
    tolerance_pct_points = 0.11
) {
  if (is.null(ref) || nrow(metric_map) == 0) return(invisible(NULL))
  
  expected_groups <- tidyr::crossing(
    mode = c("walking", "cycling", "driving"),
    time_min = c(5, 10, 15)
  )
  group_check <- ref |>
    dplyr::count(mode, time_min, name = "n") |>
    dplyr::right_join(expected_groups, by = c("mode", "time_min"))
  if (any(is.na(group_check$n))) warning("One or more expected mode × time reference groups are missing.")
  
  diagnostics <- list()
  for (j in seq_len(nrow(metric_map))) {
    m <- metric_map[j, , drop = FALSE]
    raw_col <- m$ref_raw_col[[1]]
    pct_col <- m$ref_pct_col[[1]]
    if (!(raw_col %in% names(ref)) || !(pct_col %in% names(ref))) next
    
    # Use the first complete standard group for a stored-percentile reproduction
    # check. This is deliberately separate from external user scoring.
    for (k in seq_len(nrow(expected_groups))) {
      g <- get_iso_reference_group(expected_groups$mode[k], expected_groups$time_min[k], ref)
      if (is.null(g) || nrow(g) < 3) next
      raw <- suppressWarnings(as.numeric(g[[raw_col]]))
      stored <- suppressWarnings(as.numeric(g[[pct_col]]))
      ok <- is.finite(raw) & is.finite(stored)
      raw <- raw[ok]; stored <- stored[ok]
      if (length(raw) < 3) next
      
      ord <- order(raw)
      picks <- unique(c(ord[1], ord[ceiling(length(ord)/2)], ord[length(ord)]))
      pred <- vapply(
        raw[picks], percent_rank_existing_reference, numeric(1),
        reference_values = raw, direction = m$direction[[1]]
      )
      err <- max(abs(pred - stored[picks]), na.rm = TRUE)
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        metric = m$app_col[[1]], mode = expected_groups$mode[k],
        time_min = expected_groups$time_min[k], max_abs_error = err
      )
      break
    }
  }
  
  diagnostics <- dplyr::bind_rows(diagnostics)
  if (nrow(diagnostics) > 0 && any(diagnostics$max_abs_error > tolerance_pct_points)) {
    bad <- diagnostics$max_abs_error > tolerance_pct_points
    warning(
      "Stored Step-5 percentile reproduction exceeded tolerance for: ",
      paste0(diagnostics$metric[bad], " (", round(diagnostics$max_abs_error[bad], 3), " pp)", collapse = "; ")
    )
  }
  invisible(diagnostics)
}

iso_reference_validation <- validate_iso_reference_percentiles()

score_iso_metrics_against_reference <- function(
    metrics,
    ref = iso_reference,
    metric_map = iso_reference_metric_map
) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  out <- metrics
  old_attrs <- attributes(metrics)
  
  mode_col <- pick_first_column(out, c("mode", "Mode"))
  time_col <- pick_first_column(out, c("time_min", "time", "Time"))
  
  for (nm in unique(metric_map$out_col)) out[[nm]] <- NA_real_
  for (nm in unique(metric_map$n_out_col)) out[[nm]] <- NA_integer_
  out$Reference_Mode <- NA_character_
  out$Reference_Time_min <- NA_real_
  out$Reference_N <- NA_integer_
  out$Reference_Status <- NA_character_
  
  metric_version_ok <- "MetricDefinitionVersion" %in% names(out) &&
    all(as.character(out$MetricDefinitionVersion) == ISO_METRIC_DEFINITION_VERSION)
  
  if (!metric_version_ok) {
    out$Reference_Status <- paste0(
      "N/A — raw metric definition is not ", ISO_METRIC_DEFINITION_VERSION
    )
  } else if (is.na(mode_col) || is.na(time_col) || is.null(ref)) {
    out$Reference_Status <- if (is.null(ref)) {
      "No compatible corrected isochrone reference loaded"
    } else {
      "Metrics table lacks mode/time columns"
    }
  } else {
    for (i in seq_len(nrow(out))) {
      mode_i <- normalize_reference_mode(out[[mode_col]][[i]])
      time_i <- suppressWarnings(as.numeric(out[[time_col]][[i]]))
      out$Reference_Mode[[i]] <- mode_i
      out$Reference_Time_min[[i]] <- time_i
      
      ref_i <- get_iso_reference_group(mode_i, time_i, ref = ref)
      if (is.null(ref_i) || nrow(ref_i) == 0) {
        out$Reference_Status[[i]] <- paste0(
          "N/A — no precomputed ", pretty_mode(mode_i), " ", time_i, "-minute reference"
        )
        next
      }
      
      out$Reference_N[[i]] <- nrow(ref_i)
      out$Reference_Status[[i]] <- paste0(
        "Citywide SF reference: ", nrow(ref_i), " precomputed 100 m origins · ",
        pretty_mode(mode_i), " · ", time_i, " min"
      )
      
      for (j in seq_len(nrow(metric_map))) {
        m <- metric_map[j, , drop = FALSE]
        app_col <- m$app_col[[1]]
        ref_col <- m$ref_raw_col[[1]]
        out_col <- m$out_col[[1]]
        n_out_col <- m$n_out_col[[1]]
        if (!(app_col %in% names(out)) || !(ref_col %in% names(ref_i))) next
        
        x <- suppressWarnings(as.numeric(out[[app_col]][[i]]))
        ref_values <- suppressWarnings(as.numeric(ref_i[[ref_col]]))
        out[[n_out_col]][[i]] <- sum(is.finite(ref_values))
        out[[out_col]][[i]] <- round(
          percent_rank_against_reference(
            x, ref_values, direction = m$direction[[1]]
          ),
          1
        )
      }
    }
  }
  
  custom_attr_names <- setdiff(names(old_attrs), c("names", "row.names", "class"))
  for (nm in custom_attr_names) attr(out, nm) <- old_attrs[[nm]]
  attr(out, "iso_reference_scored") <- TRUE
  out
}

# Full seven-axis BAI only. Missing dimensions are explicit; BAI is NOT silently
# redefined as a mean over a smaller denominator.
add_reference_bai <- function(metrics) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  out <- metrics
  
  axis_from_pct <- c(
    Mobility_Access_std = "pctile_Transit_Access_Score",
    Route_Access_std = "pctile_Unique_Muni_Routes",
    Biodiversity_Potential_std = "pctile_GBIF_Species",
    Observation_Intensity_std = "pctile_SamplingDensity_km2",
    Environmental_Quality_std = "pctile_MeanNDVI",
    Greenspace_Cover_std = "pctile_Greenspace_percent",
    Equity_Context_std = "pctile_SF_EJ_Score"
  )
  
  for (axis in names(axis_from_pct)) {
    pct_col <- axis_from_pct[[axis]]
    out[[axis]] <- if (pct_col %in% names(out)) {
      suppressWarnings(as.numeric(out[[pct_col]])) / 100
    } else NA_real_
  }
  
  axis_cols <- names(axis_from_pct)
  axis_mat <- as.matrix(out[, axis_cols, drop = FALSE])
  storage.mode(axis_mat) <- "numeric"
  n_components <- rowSums(is.finite(axis_mat))
  complete <- n_components == length(axis_cols)
  
  out$BAI <- NA_real_
  if (any(complete)) {
    out$BAI[complete] <- rowMeans(axis_mat[complete, , drop = FALSE])
  }
  
  out$BAI_n_components <- as.integer(n_components)
  out$BAI_reference_complete <- complete
  out$BAI_reference_system <- paste0(
    "same mode × same time precomputed 100 m isochrones; metric definition ",
    ISO_METRIC_DEFINITION_VERSION
  )
  
  # A BAI score (e.g. 51/100) is NOT itself a percentile. Rank that complete
  # seven-axis BAI against the distribution of complete seven-axis BAIs from the
  # exact same citywide mode × time reference group.
  out$pctile_BAI <- NA_real_
  out$nref_BAI <- NA_integer_
  
  mode_col <- pick_first_column(out, c("mode", "Mode"))
  time_col <- pick_first_column(out, c("time_min", "time", "Time"))
  
  if (!is.na(mode_col) && !is.na(time_col)) {
    for (i in seq_len(nrow(out))) {
      if (!isTRUE(complete[[i]]) || !is.finite(out$BAI[[i]])) next
      
      mode_i <- normalize_reference_mode(out[[mode_col]][[i]])
      time_i <- suppressWarnings(as.numeric(out[[time_col]][[i]]))
      ref_bai <- get_reference_bai_distribution(mode_i, time_i)
      
      if (length(ref_bai) == 0) next
      
      out$nref_BAI[[i]] <- length(ref_bai)
      out$pctile_BAI[[i]] <- round(
        percent_rank_against_reference(
          out$BAI[[i]],
          ref_bai,
          direction = "higher"
        ),
        1
      )
    }
  }
  
  out
}

format_pct <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x)) return("N/A")
  paste0(round(x, digits), "th percentile")
}

iso_metric_row_label <- function(df, i) {
  mode_col <- pick_first_column(df, c("Mode", "mode"))
  time_col <- pick_first_column(df, c("Time", "time", "time_min"))
  mode_i <- if (!is.na(mode_col)) as.character(df[[mode_col]][[i]]) else "Isochrone"
  time_i <- if (!is.na(time_col)) suppressWarnings(as.numeric(df[[time_col]][[i]])) else NA_real_
  if (is.finite(time_i)) paste0(mode_i, " — ", time_i, " min") else mode_i
}

# Radar plot for independently scored isochrones. Each polygon/line is one
# isochrone row; there is no averaging or pooling across mode/time combinations.
#
# plot_style controls presentation only:
#   "standard"     = Explorer / multiple independent isochrones
#   "point_compare"= Point A vs Point B
#   "before_after" = existing conditions vs planned conservation action
draw_reference_radar <- function(
    bai_rows,
    labels = NULL,
    colors = NULL,
    title = "Citywide Biodiversity Access Profile",
    subtitle = NULL,
    plot_style = c("standard", "point_compare", "before_after")
) {
  plot_style <- match.arg(plot_style)
  
  if (is.null(bai_rows)) {
    plot.new(); title("No isochrone scores available."); return(invisible(NULL))
  }
  
  if (is.data.frame(bai_rows)) {
    rows <- lapply(seq_len(nrow(bai_rows)), function(i) bai_rows[i, , drop = FALSE])
  } else if (is.list(bai_rows)) {
    rows <- bai_rows
  } else {
    plot.new(); title("No isochrone scores available."); return(invisible(NULL))
  }
  
  rows <- rows[vapply(rows, function(x) !is.null(x) && nrow(x) > 0, logical(1))]
  if (length(rows) == 0) {
    plot.new(); title("No isochrone scores available."); return(invisible(NULL))
  }
  
  if (is.null(labels)) {
    labels <- vapply(rows, function(x) iso_metric_row_label(x, 1), character(1))
  }
  if (length(labels) != length(rows)) labels <- paste0("Isochrone ", seq_along(rows))
  
  # Short enough to read around a radar, but explicit about what each spoke means.
  axis_labels <- c(
    Mobility_Access_std = "Transit\naccess ↑",
    Route_Access_std = "Muni route\naccess ↑",
    Biodiversity_Potential_std = "Total species\nrichness ↑",
    Observation_Intensity_std = "GBIF sampling\ndensity ↑",
    Environmental_Quality_std = "Vegetation\nNDVI ↑",
    Greenspace_Cover_std = "Greenspace\ncover ↑",
    Equity_Context_std = "Equity context ↑\n(lower EJ burden)"
  )
  
  # Pair related BAI dimensions visually so stakeholders can immediately see
  # which spokes belong together conceptually.
  axis_group_colors <- c(
    Mobility_Access_std = BAI_GROUP_COLORS[["access"]],
    Route_Access_std = BAI_GROUP_COLORS[["access"]],
    Biodiversity_Potential_std = BAI_GROUP_COLORS[["biodiversity"]],
    Observation_Intensity_std = BAI_GROUP_COLORS[["biodiversity"]],
    Environmental_Quality_std = BAI_GROUP_COLORS[["green_environment"]],
    Greenspace_Cover_std = BAI_GROUP_COLORS[["green_environment"]],
    Equity_Context_std = BAI_GROUP_COLORS[["equity"]]
  )
  
  all_df <- dplyr::bind_rows(lapply(seq_along(rows), function(i) {
    x <- rows[[i]]
    x$.Series <- labels[[i]]
    x
  }))
  
  available_axes <- names(axis_labels)[vapply(names(axis_labels), function(a) {
    a %in% names(all_df) && any(is.finite(suppressWarnings(as.numeric(all_df[[a]]))))
  }, logical(1))]
  
  if (length(available_axes) < 3) {
    plot.new()
    title("Not enough comparable same-mode × same-time reference dimensions for a radar plot.")
    return(invisible(NULL))
  }
  
  long <- all_df |>
    dplyr::select(.Series, dplyr::all_of(available_axes)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_axes),
      names_to = "Axis",
      values_to = "Score"
    ) |>
    dplyr::mutate(
      Score = suppressWarnings(as.numeric(Score)) * 100,
      Axis_key = Axis,
      Axis = factor(
        Axis,
        levels = available_axes,
        labels = unname(axis_labels[available_axes])
      )
    ) |>
    dplyr::filter(is.finite(Score))
  
  if (nrow(long) == 0) {
    plot.new(); title("No comparable percentile dimensions available."); return(invisible(NULL))
  }
  
  axis_label_df <- long |>
    dplyr::distinct(Axis_key, Axis) |>
    dplyr::mutate(
      axis_label = as.character(Axis),
      axis_color = unname(axis_group_colors[Axis_key]),
      label_y = 108.5
    )
  
  axis_label_layers <- lapply(seq_len(nrow(axis_label_df)), function(i) {
    ggplot2::annotate(
      "text",
      x = axis_label_df$Axis[[i]],
      y = axis_label_df$label_y[[i]],
      label = axis_label_df$axis_label[[i]],
      color = axis_label_df$axis_color[[i]],
      size = 4.0,
      fontface = "bold",
      lineheight = 1.0
    )
  })
  
  if (is.null(colors)) {
    pal <- scales::hue_pal()(length(unique(long$.Series)))
    names(pal) <- unique(long$.Series)
  } else {
    pal <- rep(colors, length.out = length(labels))
    names(pal) <- labels
  }
  
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Citywide percentile within the SAME mode × travel-time reference. ",
      "Center = 0th percentile; outer ring = 100th percentile."
    )
  }
  
  base_plot <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = Axis,
      y = Score,
      group = .Series,
      color = .Series,
      fill = .Series
    )
  )
  
  if (identical(plot_style, "before_after") && length(labels) >= 2) {
    baseline_label <- labels[[1]]
    proposal_label <- labels[[2]]
    
    # Keep both full profiles visible, but deliberately mute them. The actual
    # intervention effect is overlaid only on dimensions whose percentile changed.
    base_plot <- base_plot +
      ggplot2::geom_polygon(
        data = long |> dplyr::filter(.Series == baseline_label),
        fill = NA,
        linewidth = 1.15,
        linetype = "22",
        alpha = 0.55,
        na.rm = TRUE
      ) +
      ggplot2::geom_polygon(
        data = long |> dplyr::filter(.Series == proposal_label),
        linewidth = 1.15,
        alpha = 0.035,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = long,
        size = 2.0,
        alpha = 0.35,
        stroke = 0.7,
        na.rm = TRUE
      )
    
    # Determine changed BAI dimensions. Connectivity is intentionally NOT added
    # here because it is a separate structural outcome rather than an eighth BAI axis.
    delta_df <- long |>
      dplyr::select(.Series, Axis_key, Axis, Score) |>
      tidyr::pivot_wider(names_from = .Series, values_from = Score)
    
    if (all(c(baseline_label, proposal_label) %in% names(delta_df))) {
      delta_df <- delta_df |>
        dplyr::mutate(
          baseline_score = .data[[baseline_label]],
          proposal_score = .data[[proposal_label]],
          delta = proposal_score - baseline_score,
          delta_label = dplyr::if_else(
            is.finite(delta) & abs(delta) >= 0.1,
            sprintf("Δ %+.1f pp", delta),
            NA_character_
          )
        ) |>
        dplyr::filter(!is.na(delta_label))
      
      if (nrow(delta_df) > 0) {
        # A thick radial segment marks the exact BAI spoke that changed.
        base_plot <- base_plot +
          ggplot2::geom_segment(
            data = delta_df,
            ggplot2::aes(
              x = Axis,
              xend = Axis,
              y = baseline_score,
              yend = proposal_score
            ),
            inherit.aes = FALSE,
            color = "#238b45",
            linewidth = 3.0,
            lineend = "round",
            alpha = 0.95,
            na.rm = TRUE
          ) +
          ggplot2::geom_point(
            data = delta_df,
            ggplot2::aes(x = Axis, y = baseline_score),
            inherit.aes = FALSE,
            shape = 21,
            size = 3.8,
            stroke = 1.2,
            fill = "white",
            color = "#636363",
            na.rm = TRUE
          ) +
          ggplot2::geom_point(
            data = delta_df,
            ggplot2::aes(x = Axis, y = proposal_score),
            inherit.aes = FALSE,
            shape = 21,
            size = 4.8,
            stroke = 1.4,
            fill = "#74c476",
            color = "#1b7837",
            na.rm = TRUE
          ) +
          ggplot2::geom_label(
            data = delta_df,
            ggplot2::aes(
              x = Axis,
              y = pmin(98, pmax(baseline_score, proposal_score) + 7),
              label = delta_label
            ),
            inherit.aes = FALSE,
            size = 3.6,
            fontface = "bold",
            label.size = 0.25,
            fill = "#f3fbf5",
            color = "#1b7837",
            alpha = 0.98,
            na.rm = TRUE
          )
      }
    }
  } else if (identical(plot_style, "point_compare")) {
    base_plot <- base_plot +
      ggplot2::geom_polygon(alpha = 0.035, linewidth = 1.55, na.rm = TRUE) +
      ggplot2::geom_point(size = 3.2, stroke = 1.0, na.rm = TRUE)
  } else {
    base_plot <- base_plot +
      ggplot2::geom_polygon(alpha = 0.045, linewidth = 1.35, na.rm = TRUE) +
      ggplot2::geom_point(size = 3.0, stroke = 0.9, na.rm = TRUE)
  }
  
  print(
    base_plot +
      axis_label_layers +
      ggplot2::coord_polar(clip = "off") +
      ggplot2::scale_y_continuous(
        limits = c(0, 112),
        breaks = c(0, 25, 50, 75, 100),
        labels = c("0", "25", "50", "75", "100"),
        expand = ggplot2::expansion(mult = c(0, 0.03))
      ) +
      ggplot2::scale_color_manual(values = pal, name = NULL) +
      ggplot2::scale_fill_manual(values = pal, name = NULL) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        caption = paste0(
          "Farther outward = higher / more favourable citywide percentile. ",
          "Species = more unique GBIF species; greenspace = more cover; NDVI = greener vegetation. ",
          "Paired spoke colors: blue = access, orange = biodiversity, green = vegetation/greenspace, purple = equity context. ",
          "Equity is inverted once: farther outward = lower SF EJ burden. ",
          if (identical(plot_style, "before_after")) {
            "Only changed BAI spokes are highlighted in green; connectivity is reported separately."
          } else {
            ""
          }
        ),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major = ggplot2::element_line(color = "#d8e2dc", linewidth = 0.55),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(
          size = 9.5,
          color = "#6b766f",
          face = "bold"
        ),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "top",
        legend.justification = "center",
        legend.text = ggplot2::element_text(size = 11, face = "bold"),
        legend.key.width = grid::unit(1.5, "cm"),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 16,
          color = "#294337",
          margin = ggplot2::margin(b = 5)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 10.5,
          color = "#596861",
          margin = ggplot2::margin(b = 8)
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          size = 9.5,
          color = "#596861",
          margin = ggplot2::margin(t = 10)
        ),
        plot.margin = ggplot2::margin(18, 36, 22, 36)
      )
  )
}


# Project the existing OSM greenspace once per app process, not once per user
# session. Scenario calculations use a spatial index and union only polygons near
# the user's drawing.
osm_greenspace_3857 <- tryCatch(
  osm_greenspace |> st_make_valid() |> st_transform(3857),
  error = function(e) NULL
)

# Stakeholder-facing coverage metadata. This is the complete mapped OSM
# greenspace layer loaded for San Francisco by setup_unified.R. Scenario
# calculations query this citywide layer and then display only the patches near
# a user's proposal.
osm_greenspace_feature_count <- if (!is.null(osm_greenspace_3857)) {
  nrow(osm_greenspace_3857)
} else {
  0L
}

# Transport-mode choices, shared by the Explorer checkboxes and the Comparer
# single-select dropdowns so both tabs offer the same modes.
transport_mode_choices <- list(
  "Driving"               = "driving",
  "Walking"               = "walking",
  "Cycling"               = "cycling",
  "Driving with Traffic"  = "driving-traffic",
  "Transit (GTFS)"        = "transit",
  "Walk + Transit (Muni)" = "walk_transit"
)

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "SF Biodiversity Access Tool"),
  
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Isochrone Explorer", tabName = "isochrone", icon = icon("map-marker-alt")),
                menuItem("Isochrone Comparer", tabName = "comparer", icon = icon("balance-scale")),
                menuItem("Green Corridor Planner", tabName = "corridor", icon = icon("project-diagram")),
                menuItem("Green Investment & Displacement Context", tabName = "displacement", icon = icon("home")),
                menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "app_pastel.css"),
      
      tags$style(HTML("
        /* ================================================================
           UI polish only — no server, scoring, reactive, or metric changes
           ================================================================ */

        :root {
          --sf-green: #5d8d70;
          --sf-green-dark: #35664e;
          --sf-green-soft: #dff1e6;
          --sf-green-border: #b9dfc8;
          --sf-blue-soft: #dceafb;
          --sf-blue-border: #b5d3ee;
          --sf-gold-soft: #fff0d5;
          --sf-gold-border: #efd29a;
          --sf-bg: #f4f9f6;
          --sf-card: #ffffff;
          --sf-text: #263a31;
          --sf-muted: #66766e;
          --sf-sidebar: #1f2c33;
          --sf-sidebar-active: #273941;
          --sf-shadow: 0 4px 16px rgba(39, 72, 56, 0.08);
        }

        html, body, .wrapper {
          background: var(--sf-bg) !important;
        }

        body {
          color: var(--sf-text);
          font-size: 14px;
        }

        /* Header */
        .main-header .logo {
          background: var(--sf-green) !important;
          color: #fff !important;
          font-weight: 700;
          letter-spacing: 0.15px;
          border-right: 1px solid rgba(255,255,255,0.12);
        }

        .main-header .navbar {
          background: linear-gradient(90deg, #5f8f72 0%, #3d6f57 100%) !important;
          box-shadow: 0 1px 8px rgba(22, 51, 38, 0.20);
        }

        .main-header .navbar .sidebar-toggle:hover {
          background: rgba(255,255,255,0.09) !important;
        }

        /* Sidebar */
        .main-sidebar,
        .left-side {
          background: var(--sf-sidebar) !important;
          box-shadow: 2px 0 10px rgba(16, 34, 27, 0.08);
        }

        .sidebar-menu > li > a {
          color: #c9d4d8 !important;
          padding: 13px 15px;
          border-left: 3px solid transparent;
          transition: background-color .15s ease, color .15s ease, border-color .15s ease;
        }

        .sidebar-menu > li > a:hover {
          color: #fff !important;
          background: #25343b !important;
          border-left-color: #86b999;
        }

        .sidebar-menu > li.active > a,
        .sidebar-menu > li.active > a:hover {
          color: #fff !important;
          background: var(--sf-sidebar-active) !important;
          border-left-color: #72b48a !important;
          font-weight: 600;
        }

        .sidebar-menu > li > a > .fa,
        .sidebar-menu > li > a > .fas,
        .sidebar-menu > li > a > .far {
          width: 20px;
        }

        /* Main page */
        .content-wrapper,
        .right-side {
          background: var(--sf-bg) !important;
        }

        .content {
          padding: 16px 18px 28px 18px !important;
        }

        /* Cards / shinydashboard boxes */
        .box {
          background: var(--sf-card);
          border: 1px solid #dbe8e1 !important;
          border-top: 1px solid #dbe8e1 !important;
          border-radius: 12px !important;
          box-shadow: var(--sf-shadow);
          overflow: hidden;
          margin-bottom: 18px;
        }

        .box-header {
          padding: 12px 14px 11px 14px !important;
          border-bottom: 1px solid #e6eee9;
        }

        .box-header .box-title {
          font-size: 17px;
          line-height: 1.3;
          font-weight: 600;
          color: #294337;
        }

        .box-body {
          padding: 14px !important;
        }

        .box.box-success {
          border-color: var(--sf-green-border) !important;
        }

        .box.box-success > .box-header {
          background: var(--sf-green-soft) !important;
          color: #28543b !important;
          border-bottom-color: var(--sf-green-border) !important;
        }

        .box.box-primary,
        .box.box-info {
          border-color: var(--sf-blue-border) !important;
        }

        .box.box-primary > .box-header,
        .box.box-info > .box-header {
          background: var(--sf-blue-soft) !important;
          color: #294f74 !important;
          border-bottom-color: var(--sf-blue-border) !important;
        }

        .box.box-warning {
          border-color: var(--sf-gold-border) !important;
        }

        .box.box-warning > .box-header {
          background: var(--sf-gold-soft) !important;
          color: #745720 !important;
          border-bottom-color: var(--sf-gold-border) !important;
        }

        /* Inputs */
        .form-group {
          margin-bottom: 15px;
        }

        .control-label,
        .shiny-input-container > label {
          font-weight: 600;
          color: #31463b;
          margin-bottom: 6px;
        }

        .form-control,
        .selectize-input,
        .selectize-control.single .selectize-input {
          min-height: 38px;
          border: 1px solid #d4dfd9 !important;
          border-radius: 8px !important;
          box-shadow: none !important;
          background: #fff;
        }

        .form-control:focus,
        .selectize-input.focus {
          border-color: #86b89a !important;
          box-shadow: 0 0 0 3px rgba(93, 141, 112, 0.10) !important;
        }

        .checkbox label,
        .radio label {
          line-height: 1.55;
        }

        .checkbox input[type='checkbox'],
        .radio input[type='radio'] {
          margin-top: 4px;
        }

        .irs--shiny .irs-bar,
        .irs--shiny .irs-single,
        .irs--shiny .irs-from,
        .irs--shiny .irs-to {
          background: var(--sf-green) !important;
          border-color: var(--sf-green) !important;
        }

        .irs--shiny .irs-handle {
          border-color: var(--sf-green) !important;
        }

        /* Buttons */
        .btn {
          border-radius: 8px !important;
          font-weight: 600;
          padding: 8px 13px;
          transition: transform .08s ease, box-shadow .12s ease, background-color .12s ease;
        }

        .btn:hover {
          box-shadow: 0 3px 9px rgba(31, 63, 49, 0.12);
        }

        .btn:active {
          transform: translateY(1px);
        }

        .btn-success {
          background: #2faa67 !important;
          border-color: #2b9e60 !important;
        }

        .btn-success:hover,
        .btn-success:focus {
          background: #28985b !important;
          border-color: #258d55 !important;
        }

        .btn-default {
          background: #fff !important;
          border-color: #d5dfda !important;
          color: #3e5047 !important;
        }

        /* Leaflet / map panels */
        .leaflet {
          border-radius: 8px;
          overflow: hidden;
        }

        .leaflet-container {
          background: #eef3f0;
          font-family: inherit;
        }

        .leaflet-control-layers,
        .leaflet-bar {
          border: 1px solid rgba(50, 76, 63, 0.16) !important;
          border-radius: 7px !important;
          box-shadow: 0 2px 8px rgba(32, 57, 44, 0.10) !important;
        }

        .leaflet-control-layers {
          padding: 8px 10px;
        }

        /* Wells, alerts, and informational blocks */
        .well {
          background: #fbfdfc !important;
          border: 1px solid #dbe7e0 !important;
          border-radius: 10px !important;
          box-shadow: none !important;
        }

        .alert {
          border-radius: 9px;
          border-width: 1px;
        }

        /* Tables */
        .table {
          background: #fff;
        }

        .table > thead > tr > th {
          border-bottom: 2px solid #dce7e1;
          color: #365044;
          font-weight: 600;
        }

        .table-hover > tbody > tr:hover {
          background-color: #f3f8f5;
        }

        table.dataTable {
          border-collapse: separate !important;
          border-spacing: 0;
        }

        .dataTables_wrapper .dataTables_filter input,
        .dataTables_wrapper .dataTables_length select {
          border: 1px solid #d4dfd9;
          border-radius: 7px;
          padding: 5px 7px;
          background: #fff;
        }

        /* Typography and whitespace */
        h1, h2, h3, h4, h5 {
          color: #294337;
        }

        .help-block,
        .text-muted,
        small {
          color: var(--sf-muted);
        }

        hr {
          border-top-color: #e3ece7;
        }

        /* Loading message */
        #loading {
          color: var(--sf-green-dark) !important;
          background: #edf7f1;
          border: 1px solid #cfe5d8;
          border-radius: 8px;
          padding: 8px 12px;
          font-size: 14px !important;
          font-weight: 600;
        }

        /* Slightly more breathing room between major rows */
        .row {
          margin-bottom: 2px;
        }

        /* Mobile/tablet polish */
        @media (max-width: 991px) {
          .content {
            padding: 12px !important;
          }

          /* Responsive fallback for inline scenario summary grids. */
          .box-body [style*='grid-template-columns:repeat(4'] {
            grid-template-columns: repeat(2, minmax(0, 1fr)) !important;
          }

          .box-body [style*='grid-template-columns:repeat(3'] {
            grid-template-columns: 1fr !important;
          }

          .box {
            margin-bottom: 14px;
          }

          .box-header .box-title {
            font-size: 16px;
          }
        }
      "))
    ),
    theme = theme,
    # useShinyjs(),   # removed with library(shinyjs) -- no shinyjs features are used
    div(id = "loading", style = "display:none; font-size:20px; color:red;", "Calculating..."),
    
    tabItems(
      tabItem(
        tabName = "isochrone",
        fluidRow(
          box(
            title = "Controls", status = "success", solidHeader = TRUE, width = 3,
            
            radioButtons(
              "location_choice",
              "Select Location Method:",
              choices = c("Address (Geocode)" = "address", "Click on Map" = "map_click"),
              selected = "map_click"
            ),
            
            conditionalPanel(
              condition = "input.location_choice == 'address'",
              mapboxGeocoderInput(
                inputId = "geocoder",
                placeholder = "Search for an address",
                access_token = mapbox_token
              )
            ),
            
            checkboxGroupInput(
              "transport_modes",
              "Select Transportation Modes:",
              choices = list(
                "Driving"               = "driving",
                "Walking"               = "walking",
                "Cycling"               = "cycling",
                "Driving with Traffic"  = "driving-traffic",
                "Transit (GTFS)"        = "transit",
                "Walk + Transit (Muni)" = "walk_transit"
              ),
              selected = c("driving", "walking")
            ),
            
            conditionalPanel(
              condition = "input.transport_modes.includes('transit') || input.transport_modes.includes('walk_transit')",
              sliderInput(
                "transit_hour",
                "Transit Departure Hour (24h):",
                min = 5, max = 22, value = 9, step = 1, post = ":00"
              ),
              sliderInput(
                "transit_departure_window_min",
                "Transit departure flexibility window (minutes):",
                min = 0, max = 20, value = 10, step = 5
              ),
              helpText("Several departures after the selected time can be evaluated for walk + transit.")
            ),
            
            conditionalPanel(
              condition = "input.transport_modes.includes('walk_transit')",
              sliderInput(
                "walk_to_stop_min",
                "First-mile walking budget (minutes):",
                min = 1, max = 20, value = 5, step = 1
              ),
              sliderInput(
                "walk_from_stop_min",
                "Maximum last-mile walking budget (minutes):",
                min = 0, max = 20, value = 5, step = 1
              )
            ),
            
            checkboxGroupInput(
              "iso_times",
              "Select Isochrone Times (minutes):",
              choices = list("5" = 5, "10" = 10, "15" = 15),
              selected = c(5, 10)
            ),
            
            actionButton("generate_iso", "Generate Isochrones", icon = icon("play")),
            actionButton("clear_map", "Clear", icon = icon("times")),
            tags$hr(),
            actionButton(
              "evaluate_green_intervention",
              "Evaluate a green intervention within these isochrones",
              icon = icon("seedling"),
              class = "btn-success",
              style = "width:100%; white-space:normal;"
            ),
            helpText("Opens the Green Corridor Planner and carries the selected location and any generated isochrones with you."),
            conditionalPanel(
              condition = "output.iso_ready",
              tags$hr(style = "margin:12px 0;"),
              downloadButton(
                "download_explorer_report",
                "Download Explorer HTML report card",
                icon = icon("file-alt"),
                class = "btn-primary btn-sm",
                style = "width:100%; white-space:normal;"
              )
            )
          ),
          
          box(
            title = "Map", status = "success", solidHeader = TRUE, width = 9,
            leafletOutput("isoMap", height = 600) %>% withSpinner(type = 8, color = "#28a745")
          )
        ),
        
        fluidRow(
          box(title = "Biodiversity Access Score", status = "warning", solidHeader = TRUE, width = 4, uiOutput("bioScoreBox")),
          box(title = "Transit Access Score", status = "primary", solidHeader = TRUE, width = 4, uiOutput("transitScoreBox")),
          box(title = "Isochrone-Benchmarked Biodiversity Access Index", width = 4, uiOutput("biodiversityAccessIndexBox"))
        ),
        
        
        fluidRow(
          box(title = "Closest Greenspace", status = "success", solidHeader = TRUE, width = 6, uiOutput("closestGreenspaceUI")),
          box(title = "Closest RSF Program", status = "success", solidHeader = TRUE, width = 6, uiOutput("closestRSFProgramUI"))
        ),
        
        # Spider plot is hidden until isochrones exist (output.iso_ready), so the
        # large box doesn't sit empty before anything is generated.
        conditionalPanel(
          condition = "output.iso_ready",
          fluidRow(
            box(
              title = "Citywide Biodiversity Access Profile — BAI Dimensions",
              status = "warning", solidHeader = TRUE, width = 12, collapsible = TRUE,
              tags$div(
                style = "display:flex; gap:8px; flex-wrap:wrap; margin-bottom:7px;",
                tags$span(
                  style = "background:#eef6f1; border:1px solid #cfe2d6; border-radius:999px; padding:4px 9px; font-weight:600;",
                  "Center = 0th percentile"
                ),
                tags$span(
                  style = "background:#eef6f1; border:1px solid #cfe2d6; border-radius:999px; padding:4px 9px; font-weight:600;",
                  "Outer ring = 100th percentile"
                ),
                tags$span(
                  style = "background:#fff7e8; border:1px solid #ead7ad; border-radius:999px; padding:4px 9px; font-weight:600;",
                  "Farther out = higher / more favourable"
                )
              ),
              p(
                "Each polygon is one independently scored isochrone. The spoke labels state exactly what is being ranked. ",
                "Species richness uses the total number of unique GBIF species. Equity is oriented so farther outward means lower EJ burden. ",
                "Percentiles are always calculated against the matching SF mode × travel-time reference; profiles are never averaged or pooled."
              ),
              plotOutput("radarPlot", height = "700px") %>% withSpinner(type = 8, color = "#f0ad4e"),
              bai_dimension_group_legend_ui()
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.iso_ready",
          fluidRow(
            box(
              title = "Data Confidence / Observation Coverage",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              collapsed = TRUE,
              uiOutput("dataCoverageBox")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.iso_ready",
          fluidRow(
            box(
              title = "BAI Dimensions — Raw Values + Citywide Percentiles",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              uiOutput("baiDimensionDetailsBox")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Summary Data",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            DTOutput("dataTable") %>% withSpinner(type = 8, color = "#28a745")
          )
        ),
        
        #   fluidRow(
        #     box(
        #       title = "Biodiversity & Transit Metrics by Mode",
        #       status = "primary", solidHeader = TRUE, width = 12,
        #       plotOutput("transitMetricsPlot", height = "450px") %>% withSpinner(type = 8, color = "#005B95")
        #     )
        #   )
        # ),
        fluidRow(
          box(
            title = "Biodiversity & Transit Metrics by Mode",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            plotOutput("transitMetricsPlot", height = "450px") %>% withSpinner(type = 8, color = "#005B95")
          )
        )
        
      ),
      
      # =====================================================================
      # Isochrone Comparer tab
      # =====================================================================
      # Two independent points, one isochrone each (single mode + time), shown
      # side by side. Nothing here computes until the user clicks "Compare"
      # (cmp_results is an eventReactive), and the outputs only render once the
      # tab is opened (Shiny suspends hidden outputs) -- so users who never open
      # this tab pay no cost. All scoring reuses the shared iso_metrics_AUDITED.R functions.
      tabItem(
        tabName = "comparer",
        fluidRow(
          box(
            width = 12, status = "warning", solidHeader = FALSE,
            h4("Compare two locations side by side"),
            p("Place two points (click a map or search an address), pick one transport mode and travel time for each, then press Compare. ",
              "Each point is scored independently against its own matching mode × time reference before the spider plot and difference table are shown."),
            conditionalPanel(
              condition = "output.cmp_ready",
              downloadButton(
                "download_comparer_report",
                "Download Comparer HTML report card",
                icon = icon("file-alt"),
                class = "btn-primary btn-sm"
              )
            )
          )
        ),
        
        fluidRow(
          # --- Point A picker ---
          box(
            title = "Point A", status = "success", solidHeader = TRUE, width = 6,
            mapboxGeocoderInput(
              inputId = "cmp_geocoder_a",
              placeholder = "Search an address (Point A)",
              access_token = mapbox_token
            ),
            leafletOutput("cmp_map_a", height = 400) %>% withSpinner(type = 8, color = "#1b9e77"),
            helpText("Click the map or search an address to place Point A."),
            fluidRow(
              column(7, selectInput("cmp_mode_a", "Mode:", choices = transport_mode_choices, selected = "walking", selectize = FALSE)),
              column(5, selectInput("cmp_time_a", "Time (min):", choices = c(5, 10, 15), selected = 10, selectize = FALSE))
            )
          ),
          # --- Point B picker ---
          box(
            title = "Point B", status = "primary", solidHeader = TRUE, width = 6,
            mapboxGeocoderInput(
              inputId = "cmp_geocoder_b",
              placeholder = "Search an address (Point B)",
              access_token = mapbox_token
            ),
            leafletOutput("cmp_map_b", height = 400) %>% withSpinner(type = 8, color = "#2c7fb8"),
            helpText("Click the map or search an address to place Point B."),
            fluidRow(
              column(7, selectInput("cmp_mode_b", "Mode:", choices = transport_mode_choices, selected = "walking", selectize = FALSE)),
              column(5, selectInput("cmp_time_b", "Time (min):", choices = c(5, 10, 15), selected = 10, selectize = FALSE))
            )
          )
        ),
        
        fluidRow(
          column(
            width = 12, align = "center",
            actionButton("cmp_compare", "Compare", icon = icon("balance-scale"),
                         class = "btn-success btn-lg",
                         style = "padding: 12px 45px; font-size: 20px; font-weight: 600;"),
          )
        ),
        
        # Results stay hidden until Compare has produced something (output.cmp_ready),
        # so the spider boxes don't sit empty before the first comparison.
        conditionalPanel(
          condition = "output.cmp_ready",
          
          fluidRow(
            box(
              title = "Biodiversity Access Profile — Point A vs Point B", status = "warning", solidHeader = TRUE,
              width = 12, collapsible = TRUE,
              tags$div(
                style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap; margin-bottom:6px;",
                tags$span(
                  style = "background:#e9f6f1; border:1px solid #bfe0d2; border-radius:999px; padding:4px 9px; color:#146b51; font-weight:700;",
                  "Point A"
                ),
                tags$span("vs"),
                tags$span(
                  style = "background:#edf5fb; border:1px solid #c8dfee; border-radius:999px; padding:4px 9px; color:#27688d; font-weight:700;",
                  "Point B"
                ),
                tags$span(
                  style = "margin-left:6px; color:#66766e;",
                  "Farther outward = higher / more favourable percentile"
                )
              ),
              p(
                "Each point is ranked independently against the SF reference matching its own selected mode and travel time. ",
                "The chart compares relative citywide standing, while the difference table below pairs the raw values with their citywide percentiles."
              ),
              plotOutput("cmp_radar", height = "680px") %>% withSpinner(type = 8, color = "#f0ad4e"),
              bai_dimension_group_legend_ui()
            )
          ),
          
          fluidRow(
            box(
              title = "Difference (Point B − Point A)", status = "warning", solidHeader = TRUE, width = 12,
              p("Point B minus Point A for each metric. Percentile differences are percentile-point differences after each point has been independently benchmarked against its own matching mode × time reference. Raw values are also retained."),
              DTOutput("cmp_diff_table") %>% withSpinner(type = 8, color = "#f0ad4e")
            )
          )
        )
      ),
      
      
      # =====================================================================
      # Green Corridor Planner tab
      # =====================================================================
      tabItem(
        tabName = "corridor",
        fluidRow(
          box(
            width = 12, status = "primary", solidHeader = FALSE,
            h4(icon("project-diagram"), " Green Corridor Planner"),
            p(
              "Draw a proposed street corridor, habitat patch, pocket park, or other green intervention. ",
              "The planner compares it with mapped existing San Francisco greenspace, identifies potential habitat links, ",
              "and estimates nearby population and equity context."
            ),
            tags$div(
              style = "background:#fff3cd; border-left:4px solid #f0ad4e; padding:9px 11px;",
              tags$b("Screening tool: "),
              "Connectivity is based on a user-selected distance rule. It does not yet model species-specific movement, road barriers, vegetation quality, land ownership, or project feasibility."
            )
          )
        ),
        fluidRow(
          box(
            title = "Scenario Controls", status = "success", solidHeader = TRUE, width = 3,
            textInput(
              "corridor_scenario_name",
              "Scenario name:",
              value = "Green corridor scenario",
              placeholder = "e.g., Mission pollinator corridor"
            ),
            uiOutput("corridorLocationStatus"),
            fluidRow(
              column(
                12,
                actionButton(
                  "use_explorer_location",
                  "Use location + isochrones from Explorer",
                  icon = icon("map-marker-alt"),
                  class = "btn-success btn-sm",
                  style = "width:100%; white-space:normal; margin-bottom:5px;"
                ),
                actionButton(
                  "return_to_isochrone",
                  "Return to Isochrone Explorer",
                  icon = icon("arrow-left"),
                  class = "btn-default btn-sm",
                  style = "width:100%; white-space:normal;"
                )
              )
            ),
            uiOutput("corridorIsochroneSelectorUI"),
            tags$hr(),
            tags$div(
              style = "background:#e8f5e9; border-left:4px solid #2e7d32; padding:9px 10px; margin-bottom:10px;",
              tags$b("Uses existing San Francisco greenspace: "),
              paste0(
                "Your proposal is checked against ",
                scales::comma(osm_greenspace_feature_count),
                " mapped OSM greenspace polygons loaded for the city."
              )
            ),
            tags$p(
              style = "margin-bottom:4px;",
              tags$b("1. Draw the intervention"),
              tags$br(),
              tags$small("Use the line tool on the map for a street corridor or the polygon tool for a park, garden, or habitat patch.")
            ),
            selectInput(
              "corridor_width_m",
              "For a line, choose the total planted width:",
              choices = c(
                "2 m — narrow planting strip" = 2,
                "5 m — street trees or pollinator verge" = 5,
                "10 m — planted median or green street" = 10,
                "20 m — broad multi-layer corridor" = 20,
                "40 m — large conceptual corridor" = 40
              ),
              selected = 10,
              selectize = FALSE
            ),
            uiOutput("corridorWidthPreviewUI"),
            actionButton(
              "inspect_corridor_width",
              "Zoom in to inspect width",
              icon = icon("search-plus"),
              class = "btn-info btn-sm",
              style = "width:100%; white-space:normal; margin-top:5px;"
            ),
            helpText("This setting affects line drawings only. Polygon drawings already have an area and ignore this control."),
            tags$p(
              style = "margin:10px 0 4px 0;",
              tags$b("2. Choose the connection rule")
            ),
            selectInput(
              "connectivity_gap_m",
              "When should two habitat patches count as connected?",
              choices = c(
                "Very close — 25 m" = 25,
                "Small stepping-stone gap — 50 m" = 50,
                "Neighborhood planning screen — 100 m" = 100,
                "Broad planning screen — 250 m" = 250,
                "City-scale exploratory screen — 500 m" = 500,
                "Very broad exploratory screen — 1,000 m" = 1000
              ),
              selected = 100,
              selectize = FALSE
            ),
            helpText("This is a structural planning threshold, not a species-specific movement distance. The 500 m and 1,000 m options are broad exploratory screens and may connect many urban patches."),
            tags$p(
              style = "margin:10px 0 4px 0;",
              tags$b("3. View the result")
            ),
            actionButton(
              "show_corridor_result",
              "Show corridor result",
              icon = icon("project-diagram"),
              class = "btn-primary btn-sm"
            ),
            actionButton(
              "clear_proposed_greenspace",
              "Clear drawing",
              icon = icon("eraser"),
              class = "btn-default btn-sm"
            ),
            conditionalPanel(
              condition = "output.proposal_ready",
              tags$div(
                style = "margin-top:10px; padding:9px 10px; background:#eef6fb; border-left:4px solid #1976d2;",
                tags$b("Export this scenario"),
                tags$br(),
                tags$small("The report card is a self-contained HTML file. Open it in a browser and choose Print → Save as PDF if needed."),
                downloadButton(
                  "download_corridor_report",
                  "Download Corridor HTML report card",
                  icon = icon("file-alt"),
                  class = "btn-primary btn-sm",
                  style = "width:100%; white-space:normal; margin-top:7px;"
                )
              )
            ),
            tags$hr(style = "margin:12px 0;"),
            checkboxInput(
              "show_beneficiary_catchment",
              "Show population-estimate area (purple, optional)",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.show_beneficiary_catchment",
              sliderInput(
                "beneficiary_walk_min",
                "Assumed walking time:",
                min = 5, max = 15, value = 10, step = 5, post = " min"
              ),
              helpText("This purple area is used only to estimate nearby population and equity context. It is not habitat and does not affect connectivity.")
            ),
            tags$div(
              style = "background:#f7f9fb; border:1px solid #d7dde3; padding:9px 10px; margin:10px 0; font-size:12px;",
              tags$b("Map colors"),
              tags$div(style="margin-top:5px;", tags$span(style="display:inline-block;width:12px;height:12px;background:#81c784;border:1px solid #2e7d32;margin-right:6px;"), "Green = existing mapped greenspace"),
              tags$div(tags$span(style="display:inline-block;width:12px;height:12px;background:#9fa8da;border:2px solid #3949ab;margin-right:6px;"), "Active analysis isochrone = the only access area used for scenario BAI"),
              tags$div(tags$span(style="display:inline-block;width:12px;height:12px;background:#42a5f5;border:1px solid #0d47a1;margin-right:6px;"), "Blue = your proposed intervention"),
              tags$div(tags$span(style="display:inline-block;width:12px;height:12px;background:#ffb74d;border:1px solid #e65100;margin-right:6px;"), "Orange = existing greenspaces reached"),
              tags$div(tags$span(style="display:inline-block;width:12px;height:12px;background:#ce93d8;border:1px dashed #6a1b9a;margin-right:6px;"), "Purple = optional population-estimate area")
            ),
            uiOutput("greenspaceScenarioSummary")
          ),
          box(
            title = "Corridor Scenario Map", status = "success", solidHeader = TRUE, width = 9,
            leafletOutput("corridorMap", height = 680) %>% withSpinner(type = 8, color = "#28a745")
          )
        ),
        conditionalPanel(
          condition = "output.proposal_ready",
          fluidRow(
            box(
              title = "What This Proposal Connects",
              status = "primary", solidHeader = TRUE, width = 7, collapsible = TRUE,
              uiOutput("corridorConnectivitySummary")
            ),
            box(
              title = "Existing Greenspaces Reached",
              status = "warning", solidHeader = TRUE, width = 5, collapsible = TRUE, collapsed = TRUE,
              p("These are existing mapped greenspaces within the selected connection rule. They are highlighted in orange and numbered on the map."),
              DTOutput("linkedHabitatTable")
            )
          ),
          fluidRow(
            box(
              title = "Potentially Reached Residents",
              status = "success", solidHeader = TRUE, width = 12, collapsible = TRUE,
              uiOutput("greenspaceBeneficiarySummary")
            )
          )
        ),
        conditionalPanel(
          condition = "output.proposal_ready && !output.iso_ready",
          fluidRow(
            box(
              width = 12, status = "warning", solidHeader = FALSE,
              tags$b("Connectivity results are available now. "),
              "To compare the proposal with the Biodiversity Access Index, first generate isochrones in the Isochrone Explorer, then use ",
              tags$em("Evaluate a green intervention within these isochrones"), "."
            )
          )
        ),
        conditionalPanel(
          condition = "output.iso_ready && output.proposal_ready",
          fluidRow(
            box(
              title = "Planned Conservation Action — BAI Before vs After",
              status = "warning", solidHeader = TRUE, width = 12, collapsible = TRUE,
              uiOutput("corridorImpactSummary"),
              uiOutput("corridorBaiComparisonNote"),
              plotOutput("corridorRadarPlot", height = "650px") %>% withSpinner(type = 8, color = "#f0ad4e"),
              bai_dimension_group_legend_ui(),
              tags$small(
                style = "display:block; color:#66766e; margin-top:4px;",
                "Connectivity is displayed alongside the BAI as a separate structural-connectivity outcome; it is not inserted as an eighth BAI axis."
              )
            )
          ),
          fluidRow(
            box(
              title = "Detailed BAI Change — Selected Isochrone",
              status = "warning", solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = TRUE,
              p("Detailed standardized BAI dimensions for the selected Explorer isochrone (0–100 scale)."),
              DTOutput("greenspaceScenarioDiffTable") %>% withSpinner(type = 8, color = "#f0ad4e")
            )
          )
        )
      ),
      
      tabItem(
        tabName = "displacement",
        
        fluidRow(
          box(
            title = "Green Investment & Displacement Context",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            tags$div(
              style = paste0(
                "background:#fff7e8; border-left:5px solid #d6a542; ",
                "border-radius:8px; padding:10px 12px; margin-bottom:10px;"
              ),
              tags$b("Context screen — not a green-gentrification risk score. "),
              "This tab places greenspace access and SF environmental-justice context side by side. ",
              "The current app does not contain the housing-market variables needed to estimate displacement risk."
            ),
            p(
              "Use this view to identify places where a greening opportunity overlaps with greater environmental-justice burden, ",
              "and therefore where community governance and anti-displacement safeguards may deserve explicit consideration. ",
              "Do not interpret the plot as evidence that a proposed greening action will cause displacement."
            ),
            conditionalPanel(
              condition = "output.iso_ready",
              downloadButton(
                "download_displacement_report",
                "Download Context HTML report card",
                icon = icon("file-alt"),
                class = "btn-primary btn-sm"
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Selected Isochrone Context",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            uiOutput("displacementContextSummary")
          )
        ),
        
        fluidRow(
          box(
            title = "Greenspace Access Gap × EJ Burden Context",
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            collapsible = TRUE,
            p(
              "Each point is one independently scored isochrone. ",
              "Farther right means a larger relative greenspace-access gap; farther up means greater EJ burden context. ",
              "The 50/50 guide lines are descriptive visual references, not validated risk thresholds."
            ),
            plotOutput("displacementContextPlot", height = "540px") %>%
              withSpinner(type = 8, color = "#5d8d70")
          ),
          box(
            title = "What Is Still Needed for Displacement-Risk Screening?",
            status = "warning",
            solidHeader = TRUE,
            width = 4,
            tags$p(
              "A defensible displacement or green-gentrification screen should add housing and tenure indicators rather than infer risk from EJ or income alone."
            ),
            tags$ul(
              tags$li("Renter share / tenure"),
              tags$li("Rent burden"),
              tags$li("Eviction notices or displacement events"),
              tags$li("Recent rent and home-value change"),
              tags$li("Subsidized / permanently affordable housing"),
              tags$li("Housing production and loss"),
              tags$li("An established displacement-vulnerability classification")
            ),
            tags$div(
              style = "background:#f4f4f4; border-radius:7px; padding:9px;",
              tags$b("Current interpretation: "),
              "equity-sensitive greening context, not predicted displacement."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "NGO / Grant Reporting Language",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            uiOutput("ngoReportingSummary")
          )
        )
      ),
      
      tabItem(
        tabName = "about",
        
        # ── Logo banner ──────────────────────────────────────────────────────
        fluidRow(
          column(
            width = 12,
            style = "text-align: center; padding: 20px 0 8px 0;",
            imageOutput("combine_logo", height = "auto")
          )
        ),
        
        # ── Hero tagline ─────────────────────────────────────────────────────
        fluidRow(
          column(
            width = 12,
            style = "text-align: center; padding: 4px 40px 16px 40px;",
            tags$h3(
              style = "color: #2e8b57; font-style: italic; font-weight: 400;",
              "Exploring equitable access to urban biodiversity across San Francisco"
            )
          )
        ),
        
        # ── Green corridor conceptual framework ──────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("project-diagram"), " Green Corridor Scenario Planner"),
            status = "primary", solidHeader = TRUE, width = 12, collapsible = TRUE,
            p(
              "The corridor workflow links three questions: where habitat is fragmented, where a proposed ",
              "corridor or stepping-stone intervention could be placed, and how ecological connectivity and ",
              "equitable human access change under that scenario."
            ),
            tags$img(
              src = "green_corridor_concept.png",
              alt = "Conceptual framework showing a fragmented baseline landscape, a proposed green corridor, and ecological and human-access outcomes.",
              style = "display:block; width:100%; height:auto; margin:10px auto 4px auto; border:1px solid #d9e5dc; border-radius:6px;"
            ),
            tags$small(
              style = "color:#666;",
              "Conceptual values shown in the figure are illustrative; the app reports results calculated from the user's drawing and selected thresholds."
            )
          )
        ),
        
        # ── Row 1: Tool overview + team + GitHub ─────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("leaf"), " About This Tool"),
            status = "success", solidHeader = TRUE, width = 8,
            p(
              "The ", strong("SF Biodiversity Access Decision Support Tool"),
              " is an interactive web application developed by the ",
              strong("Reimagining San Francisco (RSF) Data Working Group"),
              " to investigate how equitably San Francisco residents can reach urban biodiversity ",
              "depending on their transportation options and socioeconomic context."
            ),
            p(
              "Users select any location in San Francisco — by clicking the map or geocoding an address — ",
              "choose one or more transport modes and travel-time thresholds, and the app generates ",
              tags$em("isochrones"),
              " (reachable-area polygons). Within each isochrone the app computes biodiversity, ",
              "greenspace, transit, socioeconomic, environmental quality, and equity metrics, and ",
              "synthesises them into the ", strong("Biodiversity Access Index (BAI)"), "."
            ),
            tags$hr(),
            fluidRow(
              column(6,
                     tags$b(icon("users"), " Team"),
                     tags$ul(
                       style = "padding-left: 18px; margin-top: 6px;",
                       tags$li("Diego Ellis Soto"),
                       tags$li("Avery Hill"),
                       tags$li("Lizzy Edson"),
                       tags$li("Álvaro Casanova"),
                       tags$li("Christopher J. Schell"),
                       tags$li("Carl Boettiger"),
                       tags$li("Rebecca Johnson")
                     )
              ),
              column(6,
                     tags$b(icon("university"), " Institutions"),
                     tags$ul(
                       style = "padding-left: 18px; margin-top: 6px;",
                       tags$li("UC Berkeley — ESPM"),
                       tags$li("California Academy of Sciences")
                     ),
                     tags$br(),
                     tags$b(icon("envelope"), " Contact"),
                     tags$p(
                       style = "margin-top: 4px;",
                       tags$a(
                         href = "mailto:diego.ellissoto@berkeley.edu",
                         "diego.ellissoto@berkeley.edu"
                       )
                     )
              )
            )
          ),
          box(
            title = tagList(icon("seedling"), " Reimagining San Francisco"),
            status = "success", solidHeader = TRUE, width = 4,
            p(
              "Reimagining San Francisco is an initiative integrating ecological, social, and ",
              "technological dimensions to shape a sustainable future for the Bay Area. The RSF Data ",
              "Working Group co-develops frameworks that bring together multiple sources of ",
              "socio-ecological biodiversity information for decision-support."
            ),
            tags$a(
              href = "https://www.calacademy.org", target = "_blank",
              icon("external-link-alt"), " California Academy of Sciences"
            ),
            tags$hr(),
            tags$b(icon("code-branch"), " Source Code"),
            tags$p(
              style = "margin-top: 6px;",
              tags$a(
                href   = "https://github.com/diego-ellis-soto/SF_biodiv_access_shiny",
                target = "_blank",
                style  = "font-size: 14px;",
                icon("github"), " diego-ellis-soto/SF_biodiv_access_shiny"
              )
            ),
            tags$div(
              style = "background-color: #d4edda; border-left: 4px solid #28a745; padding: 8px 10px; margin-top: 8px; font-size: 12px;",
              icon("info-circle"), " The full source code, setup scripts, and data pipeline are publicly available on GitHub."
            )
          )
        ),
        
        # ── Row 2: Why biodiversity access matters ───────────────────────────
        fluidRow(
          box(
            title = tagList(icon("globe-americas"), " Why Biodiversity Access Matters"),
            status = "success", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6,
                     p(
                       "Access to urban biodiversity is deeply unequal. Legacies of redlining, disinvestment, ",
                       "and car-centric planning have concentrated green, biodiverse spaces in wealthier ",
                       "neighbourhoods, while lower-income and environmental-justice communities often face ",
                       "longer travel distances and fewer transit options to reach them."
                     ),
                     p(
                       "Areas with higher biodiversity support essential ecosystem services — pollinators, ",
                       "carbon sequestration, urban heat mitigation — and provide documented cultural, ",
                       "recreational, and mental health benefits to local residents."
                     )
              ),
              column(6,
                     p(
                       "Cities are complex socio-ecological systems shaped by ongoing human pressures and ",
                       "historical decisions. The RSF initiative integrates multiple facets of biodiversity ",
                       "with variables used by city planners, public health practitioners, and equity advocates ",
                       "to support a more integrative, justice-oriented lens for urban sustainability."
                     ),
                     p(
                       "This tool is designed to make those inequities visible, quantifiable, and actionable — ",
                       "surfacing where the gaps are largest and where investment could make the greatest difference."
                     )
              )
            )
          )
        ),
        
        # ── Row 3: How to use ────────────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("map-marked-alt"), " How to Use"),
            status = "success", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6,
                     tags$ol(
                       style = "font-size: 14px; line-height: 1.85;",
                       tags$li(
                         strong("Go to the Isochrone Explorer tab"), " in the left sidebar."
                       ),
                       tags$li(
                         strong("Pick a location."),
                         tags$ul(
                           style = "margin-top: 4px;",
                           tags$li(strong("Click on Map:"), " click anywhere in SF; a red marker confirms selection."),
                           tags$li(strong("Address (Geocode):"), " type a street address in the search box.")
                         )
                       ),
                       tags$li(
                         strong("Select transport modes."),
                         " Driving, Walking, Cycling, and Driving with Traffic use Mapbox. ",
                         "Transit and Walk + Transit use the SF Muni GTFS timetable. ",
                         "Extra sliders appear for departure hour and first/last-mile walking budgets."
                       ),
                       tags$li(
                         strong("Choose time budgets"), " — 5, 10, and/or 15 minutes."
                       ),
                       tags$li(
                         strong("Click Generate Isochrones."),
                         " Shaded polygons appear for each mode × time combination. ",
                         "Click a polygon for a summary popup."
                       )
                     )
              ),
              column(6,
                     tags$ol(
                       start = 6,
                       style = "font-size: 14px; line-height: 1.85;",
                       tags$li(
                         strong("Toggle map layers"), " (top-right layer control) to overlay income, ",
                         "species richness, greenspace, CalEnviroScreen, EJ communities, transit routes, and NDVI."
                       ),
                       tags$li(
                         strong("Explore the panels below the map:"),
                         tags$ul(
                           style = "margin-top: 4px;",
                           tags$li(strong("Score boxes:"), " biodiversity access percentile, transit density, and BAI."),
                           tags$li(strong("Closest Greenspace:"), " nearest OSM green area, distance, and % greenspace cover."),
                           tags$li(strong("Spider / Radar Plot:"), " independently benchmarked per-isochrone profiles; mode × time combinations are never averaged or pooled."),
                           tags$li(strong("Summary Table:"), " full per-isochrone metrics."),
                           tags$li(strong("Metric Plots:"), " species, population, GBIF institutions, and transit by mode.")
                         )
                       ),
                       tags$li(
                         strong("Green Corridor Planner tab"), " — draw a corridor or habitat patch, identify existing greenspaces reached, and compare the proposal with Explorer isochrones."
                       ),
                       tags$li(
                         strong("Click Clear"), " to reset the map for a new query."
                       )
                     )
              )
            )
          )
        ),
        
        # ── Row 4: Map layers ────────────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("layer-group"), " Map Layers"),
            status = "primary", solidHeader = TRUE, width = 12,
            p(
              style = "font-size: 13px; color: #555;",
              "All layers are toggled in the map's layer-control panel (top-right). The default base map is CartoDB Positron; you can switch to Street or Satellite."
            ),
            tags$table(
              class = "table table-hover table-sm",
              style = "font-size: 13px;",
              tags$thead(
                tags$tr(
                  tags$th("Layer"), tags$th("Description"), tags$th("Source")
                )
              ),
              tags$tbody(
                tags$tr(tags$td(icon("dollar-sign"), " Income"),
                        tags$td("Median household income per census block group"),
                        tags$td("ACS 5-yr")),
                tags$tr(tags$td(icon("tree"), " Greenspace"),
                        tags$td("OSM parks, gardens, and public green areas"),
                        tags$td("OpenStreetMap")),
                tags$tr(tags$td(icon("ruler"), " Greenspace Distance"),
                        tags$td("Raster showing distance (m) to the nearest greenspace pixel"),
                        tags$td("Derived from OSM")),
                tags$tr(tags$td(icon("map"), " RSF Program Projects"),
                        tags$td("Partner project polygons from the Reimagining SF Initiative"),
                        tags$td("RSF Initiative")),
                tags$tr(tags$td(icon("ruler"), " RSF Program Distance"),
                        tags$td("Raster showing distance (m) to the nearest RSF program polygon"),
                        tags$td("RSF Initiative")),
                tags$tr(tags$td(icon("fire"), " Hotspots (KnowBR)"),
                        tags$td("Block groups with anomalously high species richness relative to sampling effort"),
                        tags$td("KnowBR / GBIF")),
                tags$tr(tags$td(icon("snowflake"), " Coldspots (KnowBR)"),
                        tags$td("Block groups with anomalously low species richness relative to sampling effort"),
                        tags$td("KnowBR / GBIF")),
                tags$tr(tags$td(icon("dove"), " Species Richness"),
                        tags$td("Unique GBIF species per census block group"),
                        tags$td("GBIF")),
                tags$tr(tags$td(icon("database"), " Data Availability"),
                        tags$td("Total GBIF occurrence records per block group"),
                        tags$td("GBIF")),
                tags$tr(tags$td(icon("smog"), " CalEnviroScreen (CI Score)"),
                        tags$td("Cumulative environmental and health burden by census tract"),
                        tags$td("OEHHA")),
                tags$tr(tags$td(icon("balance-scale"), " SF EJ Communities"),
                        tags$td("SF Environment Dept. environmental justice community burden scores"),
                        tags$td("SF Environment")),
                tags$tr(tags$td(icon("route"), " Transit Routes"),
                        tags$td("All SF Muni routes from GTFS shapes, coloured by official SFMTA route colour"),
                        tags$td("SFMTA GTFS")),
                tags$tr(tags$td(icon("bus"), " Transit Stops"),
                        tags$td("All SF Muni stops with AM peak headway and departure frequency"),
                        tags$td("SFMTA GTFS")),
                tags$tr(tags$td(icon("draw-polygon"), " Isochrones"),
                        tags$td("Generated travel-time polygons (one per mode × time combination)"),
                        tags$td("Mapbox / gtfsrouter")),
                tags$tr(tags$td(icon("leaf"), " NDVI Raster"),
                        tags$td("Sentinel-2 NDVI cropped and masked to the isochrone union"),
                        tags$td("Sentinel-2"))
              )
            )
          )
        ),
        
        # ── Row 5: Transport modes ───────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("car"), " Transportation Modes"),
            status = "primary", solidHeader = TRUE, width = 12,
            p(
              style = "font-size: 13px; color: #555;",
              "Six modes are supported across two routing engines. Walk + Transit is an approximation ",
              "combining a Mapbox first-mile walk, GTFS stop-to-stop reachability (SF Muni), and a last-mile ",
              "walk buffer — not a full door-to-door multimodal isochrone."
            ),
            tags$table(
              class = "table table-hover table-sm",
              style = "font-size: 13px;",
              tags$thead(
                tags$tr(
                  tags$th("Mode"), tags$th("Engine"), tags$th("Data Source"), tags$th("Notes")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td(icon("car"), " Driving"), tags$td("Mapbox"),
                  tags$td("OSM road network"), tags$td("Free-flow speed")
                ),
                tags$tr(
                  tags$td(icon("person-walking"), " Walking"), tags$td("Mapbox"),
                  tags$td("Pedestrian network"), tags$td("Pedestrian paths and crossings")
                ),
                tags$tr(
                  tags$td(icon("bicycle"), " Cycling"), tags$td("Mapbox"),
                  tags$td("Bicycle network"), tags$td("Dedicated cycle lanes where available")
                ),
                tags$tr(
                  tags$td(icon("traffic-light"), " Driving with Traffic"), tags$td("Mapbox"),
                  tags$td("Traffic-aware road network"), tags$td("Real-time + historical congestion")
                ),
                tags$tr(
                  tags$td(icon("bus"), " Transit (GTFS)"), tags$td("gtfsrouter"),
                  tags$td("SF Muni GTFS"), tags$td("Timetable-based stop-to-stop reachability from nearest stop")
                ),
                tags$tr(
                  tags$td(tagList(icon("person-walking"), "+", icon("bus"), " Walk + Transit")),
                  tags$td("Mapbox + gtfsrouter"),
                  tags$td("Pedestrian + SF Muni GTFS"),
                  tags$td("First-mile walk → Muni ride → last-mile walk within one total time budget")
                )
              )
            )
          )
        ),
        
        # ── Row 6: BAI explanation ───────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("chart-area"), " Biodiversity Access Index (BAI) & Spider Plot"),
            status = "warning", solidHeader = TRUE, width = 12,
            p(
              "The ", strong("Biodiversity Access Index (BAI)"),
              " is calculated separately for each generated isochrone. Where a comparable Step-5 reference metric exists, ",
              "the raw value is converted to a percentile using only the precomputed 100 m SF origins with the ",
              strong("same transportation mode and same travel-time budget"),
              ". Walking 10-minute results are therefore compared only with Walking 10-minute references; ",
              "different modes and travel times are never pooled. The validated reference contains 24,182 SF 100 m origin cells in each supported mode × time group. ",
              "For each metric, the percentile uses all finite comparable values from that full matching citywide group. ",
              "Raw scored metrics use the same explicit 100 m grid support as Step 5: a cell contributes when its projected centroid falls within that one isochrone."
            ),
            tags$div(
              style = "background-color:#e3f2fd; border-left:4px solid #1976d2; padding:10px 14px; margin:10px 0;",
              tags$b("Reference rule: "),
              "one user isochrone → one matching mode × time reference distribution. ",
              "If no comparable precomputed reference exists (for example transit, walk + transit, or driving with traffic), ",
              "the app reports the raw metric and an N/A percentile rather than borrowing a different benchmark."
            ),
            tags$table(
              class = "table table-bordered table-sm",
              style = "font-size: 13px; margin-top: 12px;",
              tags$thead(
                tags$tr(
                  tags$th("#"), tags$th("BAI Dimension"), tags$th("Variable"),
                  tags$th("What it measures"), tags$th("Direction")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("1"), tags$td(strong("Mobility Access")),
                  tags$td(code("Transit_Access_Score")),
                  tags$td("Muni stops per km² across the isochrone's 100 m grid support"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("2"), tags$td(strong("Route Access")),
                  tags$td(code("Unique_Muni_Routes")),
                  tags$td("Distinct Muni route IDs intersecting the selected 100 m support cells"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("3"), tags$td(strong("Biodiversity Potential")),
                  tags$td(code("GBIF_Species")),
                  tags$td("Unique GBIF species assigned to selected 100 m support cells"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("4"), tags$td(strong("Sampling Density")),
                  tags$td(code("SamplingDensity_km2")),
                  tags$td("GBIF records per km² across the selected 100 m support cells — proxy for observation coverage"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("5"), tags$td(strong("Environmental Quality")),
                  tags$td(code("MeanNDVI")),
                  tags$td("Area-weighted mean of 100 m cell NDVI values selected by the isochrone"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("6"), tags$td(strong("Greenspace Cover")),
                  tags$td(code("Greenspace_percent")),
                  tags$td("OSM greenspace area as % of the selected 100 m grid-support area"),
                  tags$td(icon("arrow-up"), " Higher = better")
                ),
                tags$tr(
                  tags$td("7"), tags$td(strong("Equity Context")),
                  tags$td(code("SF_EJ_Score (inverted)")),
                  tags$td("Lower EJ burden → more favourable access context"),
                  tags$td(icon("arrow-down"), " Lower burden = better")
                )
              )
            ),
            tags$p(
              style = "margin-top: 12px;",
              "The displayed isochrone-benchmarked BAI is calculated only when all seven percentile-standardised dimensions are available for that single isochrone. ",
              "Unsupported dimensions remain missing and are not replaced with census-block-group benchmarks. ",
              "The score box reports how many comparable reference dimensions are available; an incomplete seven-axis profile returns BAI = N/A rather than silently changing the denominator."
            ),
            tags$pre(
              style = "background: #f8f9fa; border-left: 3px solid #f0ad4e; padding: 10px;",
              "For each isochrone i:\n  percentile(metric_i) = rank within SAME mode × SAME time reference\n  BAI_i = mean(all 7 percentile-standardised dimensions for isochrone i)
  if any required axis is unavailable: BAI_i = N/A\n\nNo averaging or union across different isochrones occurs."
            ),
            tags$div(
              style = "background-color: #fff3cd; border-left: 4px solid #f0ad4e; padding: 10px 14px; margin-top: 12px;",
              icon("exclamation-triangle"), tags$b(" Prototype / Work in Progress: "),
              " Variable weighting, spatial units, and benchmark distributions should be refined through ",
              "stakeholder co-development before the BAI is used as a policy-grade score."
            )
          )
        ),
        
        # ── Row 7: Data sources ──────────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("database"), " Data Sources & Cyberinfrastructure"),
            status = "success", solidHeader = TRUE, width = 12,
            fluidRow(
              column(7,
                     tags$table(
                       class = "table table-hover table-sm",
                       style = "font-size: 13px;",
                       tags$thead(
                         tags$tr(
                           tags$th("Dataset"), tags$th("Source"), tags$th("Format"), tags$th("Use in App")
                         )
                       ),
                       tags$tbody(
                         tags$tr(tags$td(strong("GBIF occurrences")), tags$td("Global Biodiversity Information Facility"),
                                 tags$td("Parquet"), tags$td("Species richness, sampling density, taxonomic breakdowns")),
                         tags$tr(tags$td(strong("ACS / Census")), tags$td("US Census Bureau (tidycensus)"),
                                 tags$td(".Rdata"), tags$td("Population and median income per census block group")),
                         tags$tr(tags$td(strong("NDVI raster")), tags$td("Sentinel-2 (pre-processed)"),
                                 tags$td("GeoTIFF"), tags$td("Vegetation quality within isochrones")),
                         tags$tr(tags$td(strong("OSM Greenspace")), tags$td("OpenStreetMap"),
                                 tags$td("GeoPackage"), tags$td("Greenspace cover %, distance raster, map layer")),
                         tags$tr(tags$td(strong("SF Muni GTFS")), tags$td("SFMTA"),
                                 tags$td("CSV / ZIP"), tags$td("Transit isochrones, stop density, route access, headways")),
                         tags$tr(tags$td(strong("CalEnviroScreen 4.0")), tags$td("OEHHA"),
                                 tags$td("File GDB"), tags$td("Cumulative environmental burden scores")),
                         tags$tr(tags$td(strong("SF EJ Communities")), tags$td("SF Environment Dept."),
                                 tags$td("Shapefile"), tags$td("Environmental justice burden and equity context")),
                         tags$tr(tags$td(strong("Hotspots / Coldspots")), tags$td("KnowBR analysis on GBIF"),
                                 tags$td("Shapefile"), tags$td("Under- and over-observed biodiversity areas")),
                         tags$tr(tags$td(strong("RSF Program Projects")), tags$td("RSF Initiative"),
                                 tags$td("GeoPackage"), tags$td("Partner project areas overlay"))
                       )
                     )
              ),
              column(5,
                     tags$div(
                       style = "background-color: #e8f5e9; border-left: 4px solid #2e8b57; padding: 12px 14px; margin-bottom: 12px;",
                       tags$b(icon("cloud"), " Remote Data Hosting"),
                       tags$p(
                         style = "margin-top: 6px; font-size: 13px;",
                         "Greenspace, CBG, hotspots, NDVI, and GBIF data are hosted on ",
                         tags$a(
                           href = "https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access",
                           target = "_blank",
                           "HuggingFace (boettiger-lab/sf_biodiv_access)"
                         ),
                         ". The cloud deployment streams these via GDAL's ", code("/vsicurl/"),
                         " virtual filesystem — no large files need to be bundled in the Docker image."
                       )
                     ),
                     tags$div(
                       style = "background-color: #e3f2fd; border-left: 4px solid #1976D2; padding: 12px 14px; margin-bottom: 12px;",
                       tags$b(icon("bolt"), " GBIF Queries via DuckDB"),
                       tags$p(
                         style = "margin-top: 6px; font-size: 13px;",
                         "GBIF occurrence records (~3M rows for SF) are stored as a local ",
                         code(".parquet"), " file and queried on-the-fly using ",
                         tags$a(href = "https://duckdb.org/", target = "_blank", "DuckDB"),
                         " with the spatial extension. SQL ", code("ST_Intersects"),
                         " filters records to the isochrone without loading the full dataset into memory."
                       )
                     ),
                     tags$div(
                       style = "background-color: #f3e5f5; border-left: 4px solid #7B1FA2; padding: 12px 14px;",
                       tags$b(icon("github"), " Open Source"),
                       tags$p(
                         style = "margin-top: 6px; font-size: 13px;",
                         "Full source code, data pipeline scripts, and setup instructions are available at:",
                         tags$br(),
                         tags$a(
                           href = "https://github.com/diego-ellis-soto/SF_biodiv_access_shiny",
                           target = "_blank",
                           icon("github"), " github.com/diego-ellis-soto/SF_biodiv_access_shiny"
                         )
                       )
                     )
              )
            )
          )
        ),
        
        # ── Row 8: Roadmap ───────────────────────────────────────────────────
        fluidRow(
          box(
            title = tagList(icon("road"), " Status & Roadmap"),
            status = "warning", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6,
                     tags$div(
                       style = "background-color: #fff3cd; border-left: 4px solid #f0ad4e; padding: 12px 14px;",
                       tags$b(icon("flask"), " Current Status"),
                       tags$p(
                         style = "margin-top: 6px; font-size: 13px;",
                         "This tool is a ", strong("decision-support prototype"), " co-developed with the RSF Data Working Group. ",
                         "The BAI should be treated as an exploratory indicator; variable weights and reference ",
                         "distributions are subject to revision through ongoing stakeholder engagement."
                       )
                     )
              ),
              column(6,
                     tags$b(icon("tasks"), " Planned Additions"),
                     tags$ul(
                       style = "font-size: 13px; line-height: 1.85; margin-top: 8px;",
                       tags$li("Impervious surface coverage layer"),
                       tags$li("National Walkability Index integration"),
                       tags$li("CDC Social Vulnerability Index"),
                       tags$li("NatureServe biodiversity and rarity maps"),
                       tags$li("Frequency-weighted multimodal transit accessibility score"),
                       tags$li("Pre-cached transit isochrones for faster querying"),
                       tags$li("Stakeholder-driven BAI dimension weighting"),
                       tags$li("Historical and comparative isochrone analysis")
                     )
              )
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  
  
  # ---------------------------------------------------------------------------
  # Shared HTML report-card helpers
  # ---------------------------------------------------------------------------
  report_card_escape <- function(x) {
    htmltools::htmlEscape(as.character(x), attribute = FALSE)
  }
  
  report_card_table <- function(df) {
    if (is.null(df) || nrow(df) == 0) return('<div class="note">No rows available.</div>')
    if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    head_html <- paste0(
      "<tr>",
      paste0("<th>", vapply(names(df), report_card_escape, character(1)), "</th>", collapse = ""),
      "</tr>"
    )
    
    body_html <- vapply(seq_len(nrow(df)), function(i) {
      row <- df[i, , drop = FALSE]
      vals <- vapply(row, function(v) {
        v <- v[[1]]
        if (length(v) == 0 || is.na(v)) "—" else report_card_escape(v)
      }, character(1))
      paste0("<tr>", paste0("<td>", vals, "</td>", collapse = ""), "</tr>")
    }, character(1))
    
    paste0(
      '<div class="table-wrap"><table><thead>', head_html,
      '</thead><tbody>', paste(body_html, collapse = ""),
      '</tbody></table></div>'
    )
  }
  
  write_report_card <- function(file, title, subtitle = NULL, sections = character(), caution = NULL) {
    html <- paste0(
      '<!doctype html><html><head><meta charset="utf-8">',
      '<meta name="viewport" content="width=device-width,initial-scale=1">',
      '<title>', report_card_escape(title), '</title>',
      '<style>',
      'body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;margin:0;background:#f4f8f5;color:#263a31;}',
      'main{max-width:1100px;margin:0 auto;padding:30px 24px 50px;}',
      'header{background:linear-gradient(90deg,#5d8d70,#38684f);color:white;border-radius:14px;padding:22px 24px;margin-bottom:18px;}',
      'h1{margin:0 0 5px;font-size:28px;}h2{color:#294337;margin-top:0}.subtitle{opacity:.92;font-size:14px;}',
      '.section{background:white;border:1px solid #dbe8e1;border-radius:12px;padding:18px 20px;margin:14px 0;}',
      '.note{background:#eef6f1;border-left:4px solid #5d8d70;padding:10px 12px;border-radius:7px;margin:8px 0;}',
      '.caution{background:#fff6df;border-left:4px solid #d6a542;padding:10px 12px;border-radius:7px;margin:14px 0;}',
      '.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th{background:#eef6f1;text-align:left}',
      'th,td{border-bottom:1px solid #dde7e1;padding:8px 9px;vertical-align:top}footer{color:#66766e;font-size:12px;margin-top:20px}',
      '@media print{body{background:white}.section{break-inside:avoid}main{max-width:none;padding:10px}}',
      '</style></head><body><main><header><h1>', report_card_escape(title), '</h1>',
      if (!is.null(subtitle)) paste0('<div class="subtitle">', report_card_escape(subtitle), '</div>') else "",
      '</header>',
      if (!is.null(caution)) paste0('<div class="caution">', report_card_escape(caution), '</div>') else "",
      paste(sections, collapse = ""),
      '<footer>SF Biodiversity Access Decision Support Tool · Generated ',
      report_card_escape(format(Sys.time(), "%Y-%m-%d %H:%M")),
      ' · Open in a browser and use Print → Save as PDF if desired.</footer>',
      '</main></body></html>'
    )
    writeLines(html, con = file, useBytes = TRUE)
  }
  
  # Store a lightweight static spatial map inside downloadable HTML report cards.
  # The PNG is embedded as a base64 data URI, so the report remains self-contained.
  report_card_map_uri <- function(writer) {
    map_file <- tempfile(fileext = ".png")
    on.exit(unlink(map_file), add = TRUE)
    
    map_ok <- tryCatch(writer(map_file), error = function(e) FALSE)
    if (!isTRUE(map_ok) || !file.exists(map_file) ||
        !is.finite(file.info(map_file)$size) || file.info(map_file)$size <= 0 ||
        !requireNamespace("base64enc", quietly = TRUE)) {
      return(NULL)
    }
    
    paste0("data:image/png;base64,", base64enc::base64encode(map_file))
  }
  
  report_card_map_section <- function(map_uri, alt = "Static spatial context map") {
    map_html <- if (!is.null(map_uri)) {
      paste0(
        '<img src="', map_uri, '" alt="',
        htmltools::htmlEscape(as.character(alt), attribute = TRUE),
        '" style="display:block;width:100%;height:auto;max-height:620px;',
        'object-fit:contain;border:1px solid #dbe8e1;border-radius:8px;background:white;">'
      )
    } else {
      '<div style="min-height:220px;display:flex;align-items:center;justify-content:center;background:#f4f7f5;color:#66766e;border:1px solid #dbe8e1;border-radius:8px;">Static report map unavailable.</div>'
    }
    
    paste0(
      '<section class="section"><h2>Spatial context</h2>',
      '<div class="note">Static map of the isochrone area(s) used for this report card. The map is embedded in this HTML file so the spatial context is retained with the report.</div>',
      map_html,
      '</section>'
    )
  }
  
  write_report_isochrone_map <- function(
    file,
    iso_list,
    labels = NULL,
    colors = NULL,
    map_title = "Isochrone spatial context"
  ) {
    if (inherits(iso_list, "sf")) {
      iso_list <- lapply(seq_len(nrow(iso_list)), function(i) iso_list[i, , drop = FALSE])
    }
    if (!is.list(iso_list)) iso_list <- list(iso_list)
    
    iso_list <- Filter(
      function(x) !is.null(x) && inherits(x, "sf") && nrow(x) > 0,
      iso_list
    )
    if (length(iso_list) == 0) return(FALSE)
    
    iso_list <- lapply(iso_list, function(x) {
      suppressWarnings(sf::st_transform(sf::st_make_valid(x), 3857))
    })
    
    n_iso <- length(iso_list)
    if (is.null(labels) || length(labels) != n_iso) {
      labels <- paste0("Isochrone ", seq_len(n_iso))
    }
    if (is.null(colors) || length(colors) != n_iso) {
      colors <- rep("#2f6fb0", n_iso)
    }
    colors <- as.character(colors)
    colors[is.na(colors) | !nzchar(colors)] <- "#2f6fb0"
    
    focus_geom <- do.call(c, lapply(iso_list, sf::st_geometry))
    bb <- sf::st_bbox(focus_geom)
    dx <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
    dy <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
    pad_x <- max(dx * 0.18, 400)
    pad_y <- max(dy * 0.18, 400)
    map_bb <- sf::st_bbox(c(
      xmin = bb[["xmin"]] - pad_x,
      ymin = bb[["ymin"]] - pad_y,
      xmax = bb[["xmax"]] + pad_x,
      ymax = bb[["ymax"]] + pad_y
    ), crs = sf::st_crs(3857))
    
    local_cbg <- NULL
    if (exists("cbg_vect_sf", inherits = TRUE)) {
      local_cbg <- tryCatch({
        cbg <- get("cbg_vect_sf", inherits = TRUE)
        suppressWarnings(sf::st_crop(sf::st_transform(cbg, 3857), map_bb))
      }, error = function(e) NULL)
    }
    
    grDevices::png(file, width = 1500, height = 920, res = 170, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(0.4, 0.4, 1.4, 0.4), xaxs = "i", yaxs = "i")
    
    graphics::plot(
      sf::st_geometry(iso_list[[1]]),
      col = NA,
      border = NA,
      xlim = c(map_bb[["xmin"]], map_bb[["xmax"]]),
      ylim = c(map_bb[["ymin"]], map_bb[["ymax"]]),
      axes = FALSE,
      asp = 1
    )
    
    if (!is.null(local_cbg) && nrow(local_cbg) > 0) {
      graphics::plot(
        sf::st_geometry(local_cbg),
        add = TRUE,
        col = "#f6f7f6",
        border = "#d8ddda",
        lwd = 0.5
      )
    }
    
    for (i in seq_len(n_iso)) {
      graphics::plot(
        sf::st_geometry(iso_list[[i]]),
        add = TRUE,
        col = grDevices::adjustcolor(colors[[i]], alpha.f = 0.20),
        border = colors[[i]],
        lwd = 2.2
      )
    }
    
    label_points <- lapply(iso_list, function(x) {
      suppressWarnings(sf::st_point_on_surface(sf::st_union(sf::st_geometry(x))))
    })
    label_xy <- do.call(rbind, lapply(label_points, sf::st_coordinates))
    if (!is.null(label_xy) && nrow(label_xy) == n_iso) {
      graphics::points(
        label_xy[, 1], label_xy[, 2],
        pch = 21, bg = "white", col = colors,
        cex = 1.3, lwd = 1.8
      )
      graphics::text(
        label_xy[, 1], label_xy[, 2],
        labels = seq_len(n_iso),
        cex = 0.72, font = 2, col = colors
      )
    }
    
    graphics::legend(
      "bottomleft",
      legend = paste0(seq_len(n_iso), ". ", labels),
      fill = grDevices::adjustcolor(colors, alpha.f = 0.20),
      border = colors,
      bg = grDevices::adjustcolor("white", alpha.f = 0.92),
      cex = 0.72,
      bty = "o"
    )
    graphics::mtext(map_title, side = 3, line = 0.2, cex = 0.82, font = 2, col = "#294337")
    TRUE
  }
  
  # ---------------------------------------------------------------------------
  # DuckDB connection
  # ---------------------------------------------------------------------------
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  dbExecute(con, "INSTALL spatial; LOAD spatial;")
  dbExecute(con, "INSTALL httpfs; LOAD httpfs;")
  gbif_tab <- tbl(con, glue("read_parquet('{gbif_parquet}')"))
  
  onStop(function() {
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  })
  
  chosen_point <- reactiveVal(NULL)
  corridor_reference_point <- reactiveVal(NULL)
  proposed_greenspace <- reactiveVal(NULL)
  greenspace_draw_active <- reactiveVal(FALSE)
  
  # ---------------------------------------------------------------------------
  # Legacy CBG contextual distributions
  # ---------------------------------------------------------------------------
  # Retained because corridor-planning thresholds still use greenspace_cover.
  # These CBG distributions are NOT used for the same-mode × same-time
  # isochrone percentile scores or the new isochrone-benchmarked BAI.
  city_benchmarks <- local({
    bench_cache <- file.path(cache_dir, "city_benchmarks.rds")
    if (file.exists(bench_cache)) {
      message("Loading city_benchmarks from cache...")
      return(readRDS(bench_cache))
    }
    message("Computing legacy citywide contextual distributions for corridor thresholds (first run only)...")
    
    cbg_bench <- cbg_vect_sf |>
      st_transform(3857) |>
      mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
      st_transform(4326)
    
    transit_density_bench <- rep(0, nrow(cbg_bench))
    if (!is.null(gtfs_stops_sf)) {
      stop_join <- tryCatch(
        st_join(gtfs_stops_sf, cbg_bench[, c("GEOID")], left = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(stop_join) && nrow(stop_join) > 0) {
        transit_density_bench <- stop_join |>
          st_drop_geometry() |>
          count(GEOID, name = "n_stops") |>
          right_join(
            cbg_bench |> st_drop_geometry() |> select(GEOID, area_km2),
            by = "GEOID"
          ) |>
          mutate(
            n_stops = replace_na(n_stops, 0),
            transit_density = ifelse(area_km2 > 0, n_stops / area_km2, 0)
          ) |>
          pull(transit_density)
      }
    }
    
    sampling_density_bench <- cbg_bench |>
      st_drop_geometry() |>
      mutate(
        obs = ifelse(is.na(n_observations), 0, n_observations),
        sampling_density = ifelse(area_km2 > 0, obs / area_km2, 0)
      ) |>
      pull(sampling_density)
    
    biodiversity_bench <- cbg_bench |>
      st_drop_geometry() |>
      mutate(unique_species = replace_na(unique_species, 0)) |>
      pull(unique_species)
    
    ndvi_bench <- cbg_bench |>
      st_drop_geometry() |>
      mutate(ndvi_ref = dplyr::coalesce(ndvi_mean, ndvi_sentinel)) |>
      pull(ndvi_ref)
    
    ej_bench <- if (!is.null(sf_ej_sf) && "score" %in% names(sf_ej_sf)) {
      sf_ej_sf |> st_drop_geometry() |> pull(score) |> na.omit()
    } else {
      numeric(0)
    }
    
    # Greenspace coverage (%) per CBG, read from the precomputed coverage table
    # (code/prep/build_cbg_greenspace_coverage.R) -- the same source the socio
    # summary uses -- instead of recomputing the CBG x greenspace intersection.
    # gs_pct = 100 * greenspace_m2 / cbg_area_m2 reproduces the old value exactly.
    greenspace_cover_bench <- tryCatch({
      cov_path <- hf_or_local("cbg_greenspace_coverage.csv")
      readr::read_csv(cov_path, col_types = readr::cols(GEOID = readr::col_character()),
                      show_col_types = FALSE) |>
        # NA (not 0) for zero-area CBGs so the is.finite() filter below drops them,
        # matching the old inline computation exactly.
        mutate(gs_pct = ifelse(cbg_area_m2 > 0, 100 * greenspace_m2 / cbg_area_m2, NA_real_)) |>
        pull(gs_pct)
    }, error = function(e) {
      warning("Greenspace coverage benchmark failed: ", e$message)
      numeric(0)
    })
    
    # Unique Muni routes accessible per CBG (for route-access axis)
    route_access_bench <- if (!is.null(gtfs_routes_sf)) {
      message("Computing per-CBG unique route counts for BAI benchmark...")
      tryCatch({
        route_join <- st_join(
          gtfs_routes_sf[, "route_id"],
          cbg_bench[, "GEOID"],
          left = FALSE
        )
        route_join |>
          st_drop_geometry() |>
          group_by(GEOID) |>
          summarise(n_routes = n_distinct(route_id), .groups = "drop") |>
          right_join(
            cbg_bench |> st_drop_geometry() |> select(GEOID),
            by = "GEOID"
          ) |>
          mutate(n_routes = replace_na(n_routes, 0L)) |>
          pull(n_routes)
      }, error = function(e) {
        warning("Route access benchmark failed: ", e$message)
        numeric(0)
      })
    } else {
      numeric(0)
    }
    
    bench <- list(
      transit_density  = transit_density_bench[is.finite(transit_density_bench)],
      biodiversity     = biodiversity_bench[is.finite(biodiversity_bench)],
      sampling         = sampling_density_bench[is.finite(sampling_density_bench)],
      ndvi             = ndvi_bench[is.finite(ndvi_bench)],
      ej               = ej_bench[is.finite(ej_bench)],
      route_access     = route_access_bench[is.finite(route_access_bench)],
      greenspace_cover = greenspace_cover_bench[is.finite(greenspace_cover_bench)]
    )
    
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    saveRDS(bench, bench_cache)
    message("City benchmarks saved to cache.")
    bench
  })
  
  # ---------------------------------------------------------------------------
  # Logos
  # ---------------------------------------------------------------------------
  output$combine_logo <- renderImage({
    list(
      src    = file.path("www", "Combined_logos.png"),
      width  = "50%",
      height = "auto",
      alt    = "UC Berkeley ESPM · California Academy of Sciences · Reimagining San Francisco"
    )
  }, deleteFile = FALSE)
  
  # ---------------------------------------------------------------------------
  # Base map
  # ---------------------------------------------------------------------------
  output$isoMap <- renderLeaflet({
    pal_cbg  <- colorNumeric("YlOrRd", cbg_vect_sf$medincE)
    pal_rich <- colorNumeric("YlOrRd", domain = cbg_vect_sf$unique_species)
    pal_data <- colorNumeric("Blues", domain = cbg_vect_sf$n_observations)
    
    hotspot_sf <- safe_biodiv_hotspots()
    coldspot_sf <- safe_biodiv_coldspots()
    
    m <- leaflet() |>
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB.Positron") |>
      addTiles(group = "Street Map (Default)") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite (ESRI)") |>
      addPolygons(
        data = cbg_vect_sf, group = "Income",
        fillColor = ~pal_cbg(medincE), fillOpacity = 0.6,
        color = "white", weight = 1, label = ~GEOID
      ) |>
      addPolygons(
        data = osm_greenspace, group = "Existing SF Greenspace",
        fillColor = "#81c784", fillOpacity = 0.26,
        color = "#2e7d32", weight = 1.2, label = ~name
      ) |>
      addPolygons(
        data = cbg_vect_sf, group = "Species Richness",
        fillColor = ~pal_rich(unique_species), fillOpacity = 0.6,
        color = "white", weight = 1, label = ~unique_species
      ) |>
      addPolygons(
        data = cbg_vect_sf, group = "Data Availability",
        fillColor = ~pal_data(n_observations), fillOpacity = 0.6,
        color = "white", weight = 1, label = ~n_observations
      )
    
    if (FALSE && !is.null(hotspot_sf)) {
      m <- m |>
        addPolygons(
          data = hotspot_sf,
          group = "Hotspots (KnowBR)",
          fillColor = "firebrick", fillOpacity = 0.2,
          color = "firebrick", weight = 2
        )
    }
    
    if (FALSE && !is.null(coldspot_sf)) {
      m <- m |>
        addPolygons(
          data = coldspot_sf,
          group = "Coldspots (KnowBR)",
          fillColor = "navy", fillOpacity = 0.2,
          color = "navy", weight = 2
        )
    }
    
    rsf_proj_sf <- choose_existing_sf_object(c("rsf_projects", "rsf_project_areas", "rsf_polygons", "rsf_programs"))
    if (!is.null(rsf_proj_sf)) {
      m <- m |>
        addPolygons(
          data = rsf_proj_sf,
          group = "RSF Program Projects",
          fillColor = "purple",
          fillOpacity = 0.3,
          color = "purple",
          weight = 1,
          label = ~if (exists("prj_name")) ~prj_name else "RSF Project",
          highlightOptions = highlightOptions(
            weight = 5,
            color = "blue",
            fillOpacity = 0.5,
            bringToFront = TRUE
          ),
          labelOptions = labelOptions(
            style = list("font-weight" = "bold", "color" = "blue"),
            textsize = "12px",
            direction = "auto",
            noHide = FALSE
          )
        )
    }
    
    if (FALSE && exists("greenspace_dist_raster")) {
      greenspace_vals_clean <- values(greenspace_dist_raster) |>
        as.vector() |>
        (\(x) x[is.finite(x)])()
      
      if (length(greenspace_vals_clean) > 0) {
        upper_limit <- quantile(greenspace_vals_clean, 0.998, na.rm = TRUE)
        
        pal_greenspace_dist <- colorNumeric(
          palette = rev(brewer.pal(9, "YlGnBu")),
          domain = c(0, upper_limit),
          na.color = "transparent"
        )
        
        m <- m |>
          addRasterImage(
            x = greenspace_dist_raster,
            colors = pal_greenspace_dist,
            opacity = 0.65,
            project = TRUE,
            group = "Greenspace Distance"
          ) |>
          addLegend(
            position = "bottomleft",
            pal = pal_greenspace_dist,
            values = c(0, upper_limit),
            title = "Distance to<br>Greenspace (m)",
            group = "Greenspace Distance"
          )
      }
    }
    
    if (FALSE && exists("rsfprogram_dist_raster")) {
      rsf_vals_clean <- values(rsfprogram_dist_raster) |>
        as.vector() |>
        (\(x) x[is.finite(x)])()
      
      if (length(rsf_vals_clean) > 0) {
        upper_rsf <- quantile(rsf_vals_clean, 0.998, na.rm = TRUE)
        
        pal_rsf_dist <- colorNumeric(
          palette = rev(brewer.pal(9, "Purples")),
          domain = c(0, upper_rsf),
          na.color = "transparent"
        )
        
        m <- m |>
          addRasterImage(
            x = rsfprogram_dist_raster,
            colors = pal_rsf_dist,
            opacity = 0.65,
            project = TRUE,
            group = "RSF Program Distance"
          ) |>
          addLegend(
            position = "bottomright",
            pal = pal_rsf_dist,
            values = c(0, upper_rsf),
            title = "Distance to<br>RSF program (m)",
            group = "RSF Program Distance"
          )
      }
    }
    
    if (!is.null(gtfs_routes_sf)) {
      m <- m |>
        addPolylines(
          data = gtfs_routes_sf,
          group = "Transit Routes",
          color = ~route_color_hex,
          weight = 2,
          opacity = 0.8,
          label = ~paste0(route_short_name, ": ", route_long_name)
        )
    }
    
    # if (!is.null(gtfs_stops_sf)) {
    #   m <- m |>
    #     addCircleMarkers(
    #       data = gtfs_stops_sf,
    #       group = "Transit Stops",
    #       radius = 4,
    #       color = "#005B95", fillColor = "#005B95",
    #       fillOpacity = 0.7, stroke = FALSE,
    #       label = ~stop_name
    #     )
    # }
    
    # Inserted this to make reactive to URL:
    if (!is.null(gtfs_stops_sf)) {
      stops_for_map <- gtfs_stops_sf
      
      if (!"stop_url" %in% names(stops_for_map)) {
        stops_for_map$stop_url <- NA_character_
      }
      if (!"mean_headway_min" %in% names(stops_for_map)) {
        stops_for_map$mean_headway_min <- NA_real_
      }
      if (!"n_departures_peak" %in% names(stops_for_map)) {
        stops_for_map$n_departures_peak <- NA_real_
      }
      
      stops_for_map$popup_html <- paste0(
        "<strong>", stops_for_map$stop_name, "</strong>",
        "<br>Stop ID: ", stops_for_map$stop_id,
        ifelse(
          !is.na(stops_for_map$mean_headway_min),
          paste0("<br>Mean headway (AM): ", round(stops_for_map$mean_headway_min, 1), " min"),
          ""
        ),
        ifelse(
          !is.na(stops_for_map$n_departures_peak),
          paste0("<br>Departures (7–9am): ", stops_for_map$n_departures_peak),
          ""
        ),
        ifelse(
          !is.na(stops_for_map$stop_url) & stops_for_map$stop_url != "",
          paste0("<br><a href=\"", stops_for_map$stop_url, "\" target=\"_blank\">Open stop page</a>"),
          ""
        )
      )
      
      m <- m |>
        addCircleMarkers(
          data = stops_for_map,
          group = "Transit Stops",
          radius = 4,
          color = "#005B95",
          fillColor = "#005B95",
          fillOpacity = 0.7,
          stroke = FALSE,
          label = ~stop_name,
          popup = ~popup_html
        )
    }
    
    if (!is.null(cenv_sf)) {
      pal_cenv <- colorNumeric("YlOrBr", domain = cenv_sf$CIscore, na.color = "transparent")
      m <- m |>
        addPolygons(
          data = cenv_sf,
          group = "CalEnviroScreen (CI Score)",
          fillColor = ~pal_cenv(CIscore),
          fillOpacity = 0.65,
          color = "white",
          weight = 0.5,
          label = ~paste0("CI Score: ", round(CIscore, 1))
        )
    }
    
    if (!is.null(sf_ej_sf)) {
      m <- m |>
        addPolygons(
          data = sf_ej_sf,
          group = "SF EJ Communities",
          fillColor = ~symbol_hex,
          fillOpacity = 0.7,
          color = "white",
          weight = 0.5,
          label = ~paste0(ej_label, ifelse(is.na(score), "", paste0(" (score: ", score, ")")))
        )
    }
    
    m |>
      setView(lng = -122.4194, lat = 37.7749, zoom = 12) |>
      addLayersControl(
        baseGroups = c("CartoDB.Positron", "Street Map (Default)", "Satellite (ESRI)"),
        overlayGroups = c(
          "Income", "Existing SF Greenspace", "RSF Program Projects",
          "Species Richness", "Data Availability",
          "CalEnviroScreen (CI Score)", "SF EJ Communities",
          "Transit Routes", "Transit Stops",
          "Isochrones", "Transit Isochrones", "NDVI Raster"
        ),
        options = layersControlOptions(collapsed = TRUE)
      ) |>
      hideGroup("Income") |>
      hideGroup("RSF Program Projects") |>
      hideGroup("Species Richness") |>
      hideGroup("Data Availability") |>
      hideGroup("CalEnviroScreen (CI Score)") |>
      hideGroup("SF EJ Communities") |>
      hideGroup("Transit Routes") |>
      hideGroup("Transit Stops")
  })
  
  # A separate, simplified map for intervention planning. Keeping drawing and
  # scenario layers off the Isochrone Explorer prevents the two workflows from
  # competing for controls, map clicks, and visual attention.
  output$corridorMap <- renderLeaflet({
    m <- leaflet() |>
      addMapPane("corridorBasePane", zIndex = 360) |>
      addMapPane("corridorIsochronePane", zIndex = 385) |>
      addMapPane("corridorProposalPane", zIndex = 430) |>
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB.Positron") |>
      addTiles(group = "Street Map") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite (ESRI)") |>
      addPolygons(
        data = osm_greenspace,
        group = "Existing SF Greenspace",
        fillColor = "#81c784",
        fillOpacity = 0.28,
        color = "#2e7d32",
        weight = 1.2,
        label = ~ifelse(is.na(name) | name == "", "Unnamed mapped greenspace", name),
        options = pathOptions(pane = "corridorBasePane")
      )
    
    if (!is.null(sf_ej_sf) && inherits(sf_ej_sf, "sf") && nrow(sf_ej_sf) > 0) {
      m <- m |>
        addPolygons(
          data = sf_ej_sf,
          group = "SF EJ Communities",
          fillColor = ~symbol_hex,
          fillOpacity = 0.40,
          color = "#8e44ad",
          weight = 1,
          label = ~paste0(ej_label, ifelse(is.na(score), "", paste0(" (score: ", score, ")"))),
          options = pathOptions(pane = "corridorBasePane")
        )
    }
    
    m |>
      setView(lng = -122.4194, lat = 37.7749, zoom = 12) |>
      addScaleBar(
        position = "bottomright",
        options = scaleBarOptions(metric = TRUE, imperial = FALSE, maxWidth = 140)
      ) |>
      addDrawToolbar(
        targetGroup = "Your Drawn Route or Shape",
        position = "topleft",
        polylineOptions = drawPolylineOptions(
          shapeOptions = drawShapeOptions(
            color = "#1565c0",
            weight = 5,
            opacity = 0.95
          )
        ),
        polygonOptions = drawPolygonOptions(
          showArea = TRUE,
          metric = TRUE,
          shapeOptions = drawShapeOptions(
            color = "#1565c0",
            weight = 3,
            fillColor = "#90caf9",
            fillOpacity = 0.40
          )
        ),
        circleOptions = FALSE,
        rectangleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE,
        editOptions = editToolbarOptions(
          edit = TRUE,
          remove = TRUE,
          selectedPathOptions = selectedPathOptions(maintainColor = TRUE)
        )
      ) |>
      addLayersControl(
        baseGroups = c("CartoDB.Positron", "Street Map", "Satellite (ESRI)"),
        overlayGroups = c(
          "Existing SF Greenspace", "SF EJ Communities", "Explorer Location",
          "Active Analysis Isochrone",
          "Your Drawn Route or Shape", "Your Proposed Habitat Area",
          "Existing Greenspaces Reached", "Measured Connection Gaps",
          "Population Estimate Area (Optional)"
        ),
        options = layersControlOptions(collapsed = TRUE)
      ) |>
      addControl(
        html = paste0(
          '<div style="background:rgba(255,255,255,0.96); padding:9px 11px; border:1px solid #bbb; border-radius:4px; line-height:1.3; max-width:270px;">',
          '<div style="font-weight:700; margin-bottom:6px;">How to read this scenario</div>',
          '<div style="margin:4px 0;"><span style="display:inline-block;width:16px;height:12px;background:#81c784;border:2px solid #2e7d32;margin-right:7px;vertical-align:middle;"></span><b>Existing greenspace</b><br><span style="margin-left:27px;color:#555;">Mapped parks and green areas already present</span></div>',
          '<div style="margin:4px 0;"><span style="display:inline-block;width:16px;height:12px;background:#9fa8da;border:3px solid #3949ab;margin-right:7px;vertical-align:middle;"></span><b>Active analysis isochrone</b><br><span style="margin-left:27px;color:#555;">The only reachable area used for baseline-vs-proposal greenspace and BAI</span></div>',
          '<div style="margin:4px 0;"><span style="display:inline-block;width:16px;height:12px;background:#42a5f5;border:2px solid #0d47a1;margin-right:7px;vertical-align:middle;"></span><b>Your proposed intervention</b><br><span style="margin-left:27px;color:#555;">For a line, the blue area reflects the selected total width; zoom in to inspect it</span></div>',
          '<div style="margin:4px 0;"><span style="display:inline-block;width:16px;height:12px;background:#ffb74d;border:2px solid #e65100;margin-right:7px;vertical-align:middle;"></span><b>Greenspaces reached</b><br><span style="margin-left:27px;color:#555;">Existing patches within the selected rule</span></div>',
          '<div style="margin:4px 0;"><span style="display:inline-block;width:16px;height:12px;background:#ce93d8;border:2px dashed #6a1b9a;margin-right:7px;vertical-align:middle;"></span><b>Population-estimate area</b><br><span style="margin-left:27px;color:#555;">Optional; hidden until enabled</span></div>',
          '<div style="margin-top:6px;color:#555;"><span style="display:inline-block;width:18px;border-top:3px dashed #ef6c00;margin-right:6px;vertical-align:middle;"></span>Measured ecological gap</div>',
          '</div>'
        ),
        position = "bottomleft",
        className = "corridor-map-key"
      ) |>
      hideGroup("SF EJ Communities")
  })
  
  
  # Transfer the selected Explorer location and generated isochrones into the
  # dedicated corridor workflow. All generated isochrones remain available in
  # the selector, but exactly ONE isochrone is active at a time. No averaging
  # across transport modes or travel times is used for corridor BAI comparisons.
  corridor_isochrone_catalog <- function(iso_data) {
    if (is.null(iso_data) || !inherits(iso_data, "sf") || nrow(iso_data) == 0) {
      return(NULL)
    }
    
    out <- iso_data
    mode_chr <- as.character(out$mode)
    mode_label <- vapply(mode_chr, pretty_mode, character(1))
    time_num <- suppressWarnings(as.numeric(out$time))
    time_label <- ifelse(is.finite(time_num), paste0(format(time_num, trim = TRUE), " min"), "time unavailable")
    
    raw_key <- paste0(mode_chr, "__", ifelse(is.finite(time_num), time_num, "NA"))
    out$.corridor_iso_key <- make.unique(raw_key, sep = "__")
    out$.corridor_iso_label <- paste0(mode_label, " — ", time_label)
    out
  }
  
  # Default analysis context: Walking at the shortest available travel time.
  # If Walking was not generated, use the shortest available isochrone overall.
  corridor_default_iso_index <- function(iso_data) {
    if (is.null(iso_data) || nrow(iso_data) == 0) return(NA_integer_)
    
    mode_chr <- tolower(as.character(iso_data$mode))
    time_num <- suppressWarnings(as.numeric(iso_data$time))
    time_rank <- ifelse(is.finite(time_num), time_num, Inf)
    
    walking_idx <- which(mode_chr == "walking")
    if (length(walking_idx) > 0) {
      return(walking_idx[[which.min(time_rank[walking_idx])]])
    }
    
    if (any(is.finite(time_num))) return(which.min(time_rank))
    1L
  }
  
  corridor_active_iso_info <- reactive({
    iso_data <- corridor_isochrone_catalog(isochrones_data())
    if (is.null(iso_data) || nrow(iso_data) == 0) return(NULL)
    
    selected_key <- input$corridor_active_isochrone
    idx <- match(selected_key, iso_data$.corridor_iso_key)
    if (length(idx) == 0 || is.na(idx)) idx <- corridor_default_iso_index(iso_data)
    
    list(
      index = idx,
      key = iso_data$.corridor_iso_key[[idx]],
      label = iso_data$.corridor_iso_label[[idx]],
      data = iso_data[idx, , drop = FALSE],
      all = iso_data
    )
  })
  
  output$corridorIsochroneSelectorUI <- renderUI({
    iso_data <- corridor_isochrone_catalog(isochrones_data())
    if (is.null(iso_data) || nrow(iso_data) == 0) {
      return(tags$div(
        style = "font-size:12px; color:#777; margin:7px 0;",
        icon("info-circle"),
        " No generated Explorer isochrones are currently available. Generate them in the Explorer, then return here."
      ))
    }
    
    choices <- iso_data$.corridor_iso_key
    names(choices) <- iso_data$.corridor_iso_label
    selected <- isolate(input$corridor_active_isochrone)
    if (is.null(selected) || !(selected %in% choices)) {
      selected <- choices[[corridor_default_iso_index(iso_data)]]
    }
    
    tagList(
      selectInput(
        "corridor_active_isochrone",
        "Active analysis isochrone:",
        choices = choices,
        selected = selected,
        selectize = FALSE
      ),
      tags$div(
        style = "background:#e3f2fd; border-left:4px solid #1976d2; padding:8px 9px; margin:0 0 7px 0; font-size:12px;",
        tags$b("Used for: "),
        "greenspace cover and BAI baseline-versus-proposal results. This is the only isochrone used for those calculations.",
        tags$br(),
        tags$b("Not used for: "),
        "habitat connectivity or population estimates; those are calculated independently of the active isochrone."
      ),
      tags$div(
        style = "font-size:12px; color:#2e7d32; margin:-2px 0 7px 0;",
        icon("check-circle"), " ",
        "Default = Walking at the shortest generated travel time. You can select a different single mode × time above."
      )
    )
  })
  
  draw_explorer_isochrones_on_corridor <- function(
    fit = FALSE,
    active_key = NULL,
    iso_data = NULL
  ) {
    proxy <- leafletProxy("corridorMap") |>
      clearGroup("Active Analysis Isochrone")
    
    if (is.null(iso_data) || !inherits(iso_data, "sf") || nrow(iso_data) == 0) {
      return(invisible(NULL))
    }
    
    iso_display <- corridor_isochrone_catalog(iso_data)
    iso_display <- tryCatch(
      st_transform(st_make_valid(iso_display), 4326),
      error = function(e) NULL
    )
    if (is.null(iso_display) || nrow(iso_display) == 0) return(invisible(NULL))
    
    active_idx <- match(active_key, iso_display$.corridor_iso_key)
    if (length(active_idx) == 0 || is.na(active_idx)) {
      active_idx <- corridor_default_iso_index(iso_display)
    }
    
    poly_i <- iso_display[active_idx, , drop = FALSE]
    mode_i <- as.character(poly_i$mode[[1]])
    time_i <- suppressWarnings(as.numeric(poly_i$time[[1]]))
    mode_label <- pretty_mode(mode_i)
    
    is_transit <- identical(mode_i, "transit")
    is_walk_transit <- identical(mode_i, "walk_transit")
    dash_style <- if (is_transit) "6,4" else if (is_walk_transit) "4,4" else NULL
    
    proxy <- proxy |>
      addPolygons(
        data = poly_i,
        group = "Active Analysis Isochrone",
        color = "#3949ab",
        weight = 4,
        opacity = 1,
        dashArray = dash_style,
        fillColor = "#9fa8da",
        fillOpacity = 0.18,
        label = paste0(
          "ACTIVE ANALYSIS ISOCHRONE — ",
          mode_label,
          if (is.finite(time_i)) paste0(" — ", time_i, " min") else ""
        ),
        popup = paste0(
          "<strong>Active analysis isochrone</strong><br>",
          mode_label,
          if (is.finite(time_i)) paste0(" — ", time_i, " min") else "",
          "<br><small>This is the only Explorer isochrone used for baseline-versus-proposal greenspace and BAI. No averaging across isochrones is performed.</small>"
        ),
        options = pathOptions(pane = "corridorIsochronePane")
      )
    
    if (isTRUE(fit)) {
      bb <- tryCatch(
        st_bbox(st_union(st_geometry(poly_i))),
        error = function(e) NULL
      )
      if (!is.null(bb) && all(is.finite(as.numeric(bb)))) {
        leafletProxy("corridorMap") |>
          fitBounds(
            lng1 = unname(bb[["xmin"]]),
            lat1 = unname(bb[["ymin"]]),
            lng2 = unname(bb[["xmax"]]),
            lat2 = unname(bb[["ymax"]])
          )
      }
    }
    
    invisible(NULL)
  }
  
  focus_corridor_on_reference <- function(
    pt,
    zoom = 15,
    active_key = NULL,
    iso_data = NULL
  ) {
    if (is.null(pt) || length(pt) < 2) return(invisible(NULL))
    
    leafletProxy("corridorMap") |>
      clearGroup("Explorer Location") |>
      addCircleMarkers(
        lng = unname(pt[["lon"]]),
        lat = unname(pt[["lat"]]),
        radius = 7,
        color = "#212121",
        fillColor = "#ffffff",
        fillOpacity = 1,
        weight = 2,
        group = "Explorer Location",
        label = "Location transferred from Isochrone Explorer"
      )
    
    has_iso <- !is.null(iso_data) && inherits(iso_data, "sf") && nrow(iso_data) > 0
    
    draw_explorer_isochrones_on_corridor(
      fit = has_iso,
      active_key = active_key,
      iso_data = iso_data
    )
    
    if (!has_iso) {
      leafletProxy("corridorMap") |>
        flyTo(
          lng = unname(pt[["lon"]]),
          lat = unname(pt[["lat"]]),
          zoom = zoom
        )
    }
    invisible(NULL)
  }
  
  transfer_explorer_location <- function() {
    pt <- chosen_point()
    if (is.null(pt)) {
      showNotification(
        "Select a location in the Isochrone Explorer first.",
        type = "warning"
      )
      return(invisible(NULL))
    }
    
    corridor_reference_point(pt)
    updateTabItems(session, "tabs", selected = "corridor")
    
    iso_data <- isochrones_data()
    if (is.null(iso_data) || !inherits(iso_data, "sf") || nrow(iso_data) == 0) {
      showNotification(
        "The location was transferred. Generate isochrones in the Explorer to display reachable areas here.",
        type = "message",
        duration = 6
      )
    }
    
    # Capture reactive values before entering onFlushed(). The callback itself
    # must not read input$... or reactive expressions directly.
    iso_snapshot <- iso_data
    active_snapshot <- input$corridor_active_isochrone
    session$onFlushed(function() {
      focus_corridor_on_reference(
        pt,
        active_key = active_snapshot,
        iso_data = iso_snapshot
      )
    }, once = TRUE)
    invisible(NULL)
  }
  
  observeEvent(input$evaluate_green_intervention, {
    transfer_explorer_location()
  })
  
  observeEvent(input$use_explorer_location, {
    transfer_explorer_location()
  })
  
  observeEvent(input$return_to_isochrone, {
    updateTabItems(session, "tabs", selected = "isochrone")
  })
  
  # Keep the corridor map synchronized with exactly one active isochrone when
  # the selector changes or Explorer isochrones are regenerated/cleared.
  observe({
    active_key <- input$corridor_active_isochrone
    iso_snapshot <- isochrones_data()
    
    session$onFlushed(function() {
      draw_explorer_isochrones_on_corridor(
        fit = FALSE,
        active_key = active_key,
        iso_data = iso_snapshot
      )
    }, once = TRUE)
  })
  
  observeEvent(input$tabs, {
    if (identical(input$tabs, "corridor")) {
      pt <- corridor_reference_point()
      active_key <- input$corridor_active_isochrone
      iso_snapshot <- isochrones_data()
      
      if (!is.null(pt)) {
        session$onFlushed(function() {
          focus_corridor_on_reference(
            pt,
            active_key = active_key,
            iso_data = iso_snapshot
          )
        }, once = TRUE)
      } else {
        session$onFlushed(function() {
          draw_explorer_isochrones_on_corridor(
            fit = FALSE,
            active_key = active_key,
            iso_data = iso_snapshot
          )
        }, once = TRUE)
      }
    }
  }, ignoreInit = TRUE)
  
  corridor_drawing_mode <- reactive({
    proposal <- proposed_greenspace()
    if (is.null(proposal) || nrow(proposal) == 0) return("none")
    
    types <- as.character(st_geometry_type(st_zm(proposal, drop = TRUE, what = "ZM")))
    has_line <- any(types %in% c("LINESTRING", "MULTILINESTRING"))
    has_polygon <- any(types %in% c("POLYGON", "MULTIPOLYGON"))
    
    if (has_line && has_polygon) return("mixed")
    if (has_line) return("line")
    if (has_polygon) return("polygon")
    "other"
  })
  
  output$corridorWidthPreviewUI <- renderUI({
    width_m <- suppressWarnings(as.numeric(input$corridor_width_m))
    if (!is.finite(width_m) || width_m <= 0) width_m <- 10
    half_width <- width_m / 2
    mode <- corridor_drawing_mode()
    
    # A deliberately exaggerated cross-section makes the selected width legible
    # even when the map is zoomed out. It is explanatory, not map-scale.
    preview_px <- round(34 + 166 * (width_m - 2) / 38)
    preview_px <- max(34, min(200, preview_px))
    
    mode_note <- switch(
      mode,
      line = tags$div(
        style = "margin-top:6px; color:#0d47a1;",
        icon("check-circle"),
        " A line is currently drawn, so this width is being used."
      ),
      polygon = tags$div(
        style = "margin-top:6px; color:#8a6d3b;",
        icon("info-circle"),
        " A polygon is currently drawn. The width setting does not change polygons."
      ),
      mixed = tags$div(
        style = "margin-top:6px; color:#8a6d3b;",
        icon("info-circle"),
        " The width is applied only to the line portion; polygon portions keep their drawn shape."
      ),
      tags$div(
        style = "margin-top:6px; color:#666;",
        "Draw a line to apply this width, or draw a polygon to define area directly."
      )
    )
    
    area_note <- NULL
    if (mode %in% c("line", "mixed")) {
      proposal <- proposed_greenspace()
      proposal_3857 <- tryCatch(
        proposal |> st_transform(3857) |> st_zm(drop = TRUE, what = "ZM"),
        error = function(e) NULL
      )
      if (!is.null(proposal_3857) && nrow(proposal_3857) > 0) {
        types <- as.character(st_geometry_type(proposal_3857))
        line_idx <- types %in% c("LINESTRING", "MULTILINESTRING")
        if (any(line_idx)) {
          route_length_m <- sum(as.numeric(st_length(st_geometry(proposal_3857[line_idx, , drop = FALSE]))), na.rm = TRUE)
          estimated_line_area_ha <- route_length_m * width_m / 10000
          if (is.finite(route_length_m) && route_length_m > 0) {
            area_note <- tags$div(
              style = "margin-top:6px; font-size:12px; color:#444;",
              paste0(
                "Current line length: ", scales::comma(round(route_length_m)),
                " m. Approximate planted area from the line: ",
                round(estimated_line_area_ha, 2), " ha before overlap corrections."
              )
            )
          }
        }
      }
    }
    
    tags$div(
      style = "background:#eef6fb; border:1px solid #b8d8ef; border-left:4px solid #1976d2; padding:8px 9px; margin:6px 0 4px 0;",
      tags$div(
        style = "font-size:12px; color:#333; margin-bottom:5px;",
        tags$b(paste0(width_m, " m total planted width")),
        paste0(" = ", half_width, " m on each side of the line")
      ),
      tags$div(
        style = "height:30px; display:flex; justify-content:center; align-items:center; position:relative;",
        tags$div(
          style = paste0(
            "width:", preview_px, "px; height:22px; background:#42a5f5; border:2px solid #0d47a1; position:relative;"
          ),
          tags$div(style = "position:absolute; left:50%; top:-4px; bottom:-4px; border-left:2px dashed #263238;")
        )
      ),
      tags$div(
        style = "display:flex; justify-content:space-between; font-size:11px; color:#555; max-width:220px; margin:0 auto;",
        tags$span(paste0(half_width, " m")),
        tags$span("drawn centerline"),
        tags$span(paste0(half_width, " m"))
      ),
      tags$small(
        style = "display:block; margin-top:5px; color:#666;",
        "The cross-section is enlarged for explanation. Use the zoom button to inspect the actual map footprint."
      ),
      mode_note,
      area_note
    )
  })
  
  output$corridorLocationStatus <- renderUI({
    pt <- corridor_reference_point()
    if (is.null(pt)) {
      return(tags$div(
        style = "background:#f7f9fb; border-left:4px solid #9e9e9e; padding:8px 10px; margin-bottom:10px;",
        tags$b("No Explorer location transferred."),
        tags$br(),
        tags$small("You can still draw anywhere in San Francisco, or return to the Explorer and select a location first.")
      ))
    }
    
    tags$div(
      style = "background:#f5f5f5; border-left:4px solid #212121; padding:8px 10px; margin-bottom:10px;",
      tags$b("Explorer location carried into this tab"),
      tags$br(),
      tags$small(
        paste0(
          "Longitude ", round(unname(pt[["lon"]]), 5),
          "; latitude ", round(unname(pt[["lat"]]), 5),
          ". The black-and-white point is a reference location only."
        )
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # Proposed greenspace drawing scenario
  # ---------------------------------------------------------------------------
  observeEvent(input$corridorMap_draw_start, {
    greenspace_draw_active(TRUE)
  })
  observeEvent(input$corridorMap_draw_stop, {
    greenspace_draw_active(FALSE)
  })
  observeEvent(input$corridorMap_draw_editstart, {
    greenspace_draw_active(TRUE)
  })
  observeEvent(input$corridorMap_draw_editstop, {
    greenspace_draw_active(FALSE)
  })
  
  observeEvent(input$corridorMap_draw_all_features, {
    proposed_greenspace(
      tryCatch(
        draw_features_to_sf(input$corridorMap_draw_all_features),
        error = function(e) {
          showNotification(
            paste("Could not read the proposed habitat drawing:", e$message),
            type = "error"
          )
          NULL
        }
      )
    )
  }, ignoreNULL = FALSE)
  
  observeEvent(input$clear_proposed_greenspace, {
    proposed_greenspace(NULL)
    leafletProxy("corridorMap") |>
      clearGroup("Your Drawn Route or Shape") |>
      clearGroup("Your Proposed Habitat Area") |>
      clearGroup("Existing Greenspaces Reached") |>
      clearGroup("Measured Connection Gaps") |>
      clearGroup("Population Estimate Area (Optional)")
  })
  
  # Convert all user drawings to habitat area in EPSG:3857. Polygon drawings are
  # used directly; line drawings are buffered by half the selected full corridor
  # width on each side. Geometry is unioned without carrying Leaflet attributes.
  proposed_habitat_3857 <- reactive({
    proposal <- proposed_greenspace()
    if (is.null(proposal) || nrow(proposal) == 0) return(NULL)
    
    proposal_3857 <- proposal |>
      st_transform(3857) |>
      st_zm(drop = TRUE, what = "ZM") |>
      st_make_valid()
    
    geom_type <- as.character(st_geometry_type(proposal_3857))
    polygon_idx <- geom_type %in% c("POLYGON", "MULTIPOLYGON")
    line_idx <- geom_type %in% c("LINESTRING", "MULTILINESTRING")
    parts <- list()
    
    if (any(polygon_idx)) {
      parts[[length(parts) + 1]] <- st_geometry(proposal_3857[polygon_idx, , drop = FALSE])
    }
    
    if (any(line_idx)) {
      corridor_width <- suppressWarnings(as.numeric(input$corridor_width_m))
      if (!is.finite(corridor_width) || corridor_width <= 0) corridor_width <- 10
      parts[[length(parts) + 1]] <- st_buffer(
        st_geometry(proposal_3857[line_idx, , drop = FALSE]),
        dist = corridor_width / 2,
        endCapStyle = "ROUND",
        joinStyle = "ROUND"
      )
    }
    
    if (length(parts) == 0) return(NULL)
    combined <- do.call(c, parts)
    proposal_union <- suppressWarnings(st_union(combined))
    out <- st_sf(geometry = proposal_union, crs = 3857) |> st_make_valid()
    out <- out[!st_is_empty(out), , drop = FALSE]
    if (nrow(out) == 0) NULL else out
  })
  
  # Remove portions of a proposal that already overlap mapped OSM greenspace.
  # This prevents double-counting when a user draws over an existing park.
  net_proposed_greenspace_3857 <- reactive({
    proposal_sf <- proposed_habitat_3857()
    if (is.null(proposal_sf) || nrow(proposal_sf) == 0) return(NULL)
    
    if (is.null(osm_greenspace_3857) || nrow(osm_greenspace_3857) == 0) {
      return(proposal_sf)
    }
    
    hits <- st_intersects(proposal_sf, osm_greenspace_3857, sparse = TRUE)[[1]]
    if (length(hits) == 0) return(proposal_sf)
    
    existing_local <- st_union(st_geometry(osm_greenspace_3857[hits, , drop = FALSE]))
    net_geometry <- tryCatch(
      suppressWarnings(st_difference(st_geometry(proposal_sf), existing_local)),
      error = function(e) st_geometry(proposal_sf)
    )
    
    out <- st_sf(geometry = net_geometry, crs = 3857) |> st_make_valid()
    out <- out[!st_is_empty(out), , drop = FALSE]
    if (nrow(out) == 0) NULL else out
  })
  
  # Generic structural connectivity around the proposed intervention. Existing
  # OSM greenspace polygons are treated as habitat patches. The complete SF-wide
  # layer is searched, but only patches near the proposal are returned and drawn.
  # A user-selected gap threshold is a structural planning proxy rather than a
  # species-specific functional-connectivity model.
  corridor_connectivity <- reactive({
    proposal_area <- proposed_habitat_3857()
    proposal_features <- proposed_greenspace()
    if (is.null(proposal_area) || nrow(proposal_area) == 0 ||
        is.null(osm_greenspace_3857) || nrow(osm_greenspace_3857) == 0) {
      return(list(
        available = FALSE,
        message = "Draw a corridor line or habitat polygon to estimate connectivity."
      ))
    }
    
    gap_m <- suppressWarnings(as.numeric(input$connectivity_gap_m))
    if (!is.finite(gap_m) || gap_m < 0) gap_m <- 100
    
    existing_geom <- st_geometry(osm_greenspace_3857) |> st_make_valid()
    raw_names <- if ("name" %in% names(osm_greenspace_3857)) {
      as.character(osm_greenspace_3857$name)
    } else {
      rep(NA_character_, length(existing_geom))
    }
    raw_names[is.na(raw_names) | trimws(raw_names) == ""] <- NA_character_
    
    patches <- st_sf(
      patch_id = seq_along(existing_geom),
      patch_name = raw_names,
      geometry = existing_geom,
      crs = 3857
    )
    patches$patch_area_ha <- as.numeric(st_area(patches)) / 10000
    
    nearby_idx <- which(
      lengths(st_is_within_distance(patches, proposal_area, dist = gap_m)) > 0
    )
    
    feature_types <- if (is.null(proposal_features)) {
      character(0)
    } else {
      as.character(st_geometry_type(proposal_features))
    }
    line_idx <- feature_types %in% c("LINESTRING", "MULTILINESTRING")
    poly_idx <- feature_types %in% c("POLYGON", "MULTIPOLYGON")
    corridor_length_m <- 0
    if (any(line_idx)) {
      corridor_length_m <- sum(
        as.numeric(st_length(st_transform(proposal_features[line_idx, , drop = FALSE], 3857))),
        na.rm = TRUE
      )
    }
    
    novel <- net_proposed_greenspace_3857()
    novel_ha <- if (is.null(novel) || nrow(novel) == 0) {
      0
    } else {
      sum(as.numeric(st_area(novel)), na.rm = TRUE) / 10000
    }
    
    if (length(nearby_idx) == 0) {
      nearest_idx <- tryCatch(st_nearest_feature(proposal_area, patches), error = function(e) integer(0))
      nearest_patch_m <- if (length(nearest_idx) > 0 && is.finite(nearest_idx[[1]])) {
        tryCatch(
          as.numeric(st_distance(proposal_area, patches[nearest_idx[[1]], , drop = FALSE], by_element = TRUE)),
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      }
      
      return(list(
        available = TRUE,
        classification = "Isolated intervention",
        gap_m = gap_m,
        linked_patches = 0,
        overlapping_patches = 0,
        baseline_components = 0,
        scenario_components = 1,
        components_reduced = 0,
        nearest_patch_m = nearest_patch_m,
        maximum_linked_gap_m = NA_real_,
        existing_linked_area_ha = 0,
        novel_area_ha = novel_ha,
        scenario_connected_area_ha = novel_ha,
        corridor_length_m = corridor_length_m,
        drawn_patch_count = sum(poly_idx),
        linked_sf = patches[0, , drop = FALSE],
        links_sf = st_sf(
          display_id = character(0),
          distance_m = numeric(0),
          geometry = st_sfc(crs = 3857)
        )
      ))
    }
    
    linked <- patches[nearby_idx, , drop = FALSE]
    proposal_union <- st_union(st_geometry(proposal_area))
    prop_dist <- vapply(
      seq_len(nrow(linked)),
      function(i) {
        as.numeric(st_distance(st_geometry(linked)[i], proposal_union))[[1]]
      },
      numeric(1)
    )
    nearest_patch_m <- min(prop_dist, na.rm = TRUE)
    overlapping_patches <- sum(prop_dist <= 0.01, na.rm = TRUE)
    
    # Identify baseline network components among the affected existing patches.
    adjacency <- st_is_within_distance(linked, linked, dist = gap_m)
    n <- nrow(linked)
    visited <- rep(FALSE, n)
    component_id <- integer(n)
    baseline_components <- 0L
    for (i in seq_len(n)) {
      if (visited[i]) next
      baseline_components <- baseline_components + 1L
      queue <- i
      visited[i] <- TRUE
      component_id[i] <- baseline_components
      while (length(queue) > 0) {
        current <- queue[[1]]
        queue <- queue[-1]
        neighbors <- setdiff(adjacency[[current]], current)
        new_nodes <- neighbors[!visited[neighbors]]
        if (length(new_nodes) > 0) {
          visited[new_nodes] <- TRUE
          component_id[new_nodes] <- baseline_components
          queue <- c(queue, new_nodes)
        }
      }
    }
    
    linked$baseline_component <- component_id
    linked$distance_to_proposal_m <- prop_dist
    linked$relationship <- ifelse(
      prop_dist <= 0.01,
      "Direct overlap / expansion",
      paste0(round(prop_dist), " m gap")
    )
    linked$display_id <- paste0("P", seq_len(nrow(linked)))
    linked$display_name <- ifelse(
      is.na(linked$patch_name),
      paste0("Unnamed mapped greenspace ", linked$display_id),
      linked$patch_name
    )
    
    # Draw a dashed nearest-distance segment from the proposal to every affected
    # patch. These segments make the selected connection threshold visible.
    link_geoms <- lapply(seq_len(nrow(linked)), function(i) {
      tryCatch(
        st_nearest_points(proposal_union, st_geometry(linked)[i])[[1]],
        error = function(e) st_geometrycollection()
      )
    })
    links_sf <- st_sf(
      display_id = linked$display_id,
      patch_name = linked$display_name,
      distance_m = linked$distance_to_proposal_m,
      geometry = st_sfc(link_geoms, crs = 3857)
    )
    links_sf <- links_sf[!st_is_empty(links_sf), , drop = FALSE]
    
    linked_area_ha <- as.numeric(st_area(st_union(st_geometry(linked)))) / 10000
    maximum_linked_gap_m <- max(prop_dist, na.rm = TRUE)
    
    classification <- dplyr::case_when(
      nrow(linked) >= 2 && baseline_components >= 2 ~ "Bridge",
      nrow(linked) >= 2 && baseline_components == 1 ~ "Network reinforcement",
      nrow(linked) == 1 && overlapping_patches >= 1 ~ "Patch expansion",
      nrow(linked) >= 1 ~ "Stepping-stone",
      TRUE ~ "Isolated intervention"
    )
    
    list(
      available = TRUE,
      classification = classification,
      gap_m = gap_m,
      linked_patches = nrow(linked),
      overlapping_patches = overlapping_patches,
      baseline_components = baseline_components,
      scenario_components = 1L,
      components_reduced = max(0L, baseline_components - 1L),
      nearest_patch_m = nearest_patch_m,
      maximum_linked_gap_m = maximum_linked_gap_m,
      existing_linked_area_ha = linked_area_ha,
      novel_area_ha = novel_ha,
      scenario_connected_area_ha = linked_area_ha + novel_ha,
      corridor_length_m = corridor_length_m,
      drawn_patch_count = sum(poly_idx),
      linked_sf = linked,
      links_sf = links_sf
    )
  })
  
  # Display the full buffered intervention footprint so line width is visible.
  observe({
    footprint <- proposed_habitat_3857()
    proxy <- leafletProxy("corridorMap") |> clearGroup("Your Proposed Habitat Area")
    if (!is.null(footprint) && nrow(footprint) > 0) {
      width_m <- suppressWarnings(as.numeric(input$corridor_width_m))
      if (!is.finite(width_m) || width_m <= 0) width_m <- 10
      mode <- corridor_drawing_mode()
      gross_area_ha <- sum(as.numeric(st_area(footprint)), na.rm = TRUE) / 10000
      width_text <- if (mode %in% c("line", "mixed")) {
        paste0("Assumed total line width: ", width_m, " m")
      } else {
        "Polygon area follows the shape you drew; line width is not applied"
      }
      
      footprint_map <- st_transform(footprint, 4326)
      footprint_map$popup_html <- paste0(
        "<strong>Your proposed habitat area</strong>",
        "<br>", width_text,
        "<br>Modeled footprint area: ", round(gross_area_ha, 2), " ha",
        "<br><em>Zoom in to distinguish narrow width choices.</em>"
      )
      proxy |> addPolygons(
        data = footprint_map,
        group = "Your Proposed Habitat Area",
        color = "#0d47a1", weight = 4,
        fillColor = "#42a5f5", fillOpacity = 0.42,
        label = if (mode %in% c("line", "mixed")) {
          paste0("Proposed habitat: ", width_m, " m total width")
        } else {
          "Your proposed habitat polygon"
        },
        popup = ~popup_html,
        options = pathOptions(pane = "corridorProposalPane")
      )
    }
  })
  
  # Highlight the specific existing greenspaces affected by the proposal and
  # draw the threshold links used by the connectivity calculation.
  observe({
    connectivity <- corridor_connectivity()
    proxy <- leafletProxy("corridorMap") |>
      clearGroup("Existing Greenspaces Reached") |>
      clearGroup("Measured Connection Gaps")
    
    if (is.null(connectivity) || !isTRUE(connectivity$available) ||
        is.null(connectivity$linked_sf) || nrow(connectivity$linked_sf) == 0) {
      return(invisible(NULL))
    }
    
    linked_map <- st_transform(connectivity$linked_sf, 4326)
    linked_map$popup_html <- paste0(
      "<strong>", linked_map$display_id, ": ", linked_map$display_name, "</strong>",
      "<br>Relationship: ", linked_map$relationship,
      "<br>Existing area: ", round(linked_map$patch_area_ha, 2), " ha",
      "<br>Baseline network: ", linked_map$baseline_component
    )
    
    proxy <- proxy |>
      addPolygons(
        data = linked_map,
        group = "Existing Greenspaces Reached",
        color = "#e65100",
        weight = 4,
        fillColor = "#ffb74d",
        fillOpacity = 0.48,
        label = ~paste0(display_id, ": ", display_name),
        popup = ~popup_html,
        highlightOptions = highlightOptions(
          weight = 6,
          color = "#bf360c",
          fillOpacity = 0.58,
          bringToFront = TRUE
        )
      )
    
    label_points <- suppressWarnings(st_point_on_surface(linked_map))
    proxy <- proxy |>
      addLabelOnlyMarkers(
        data = label_points,
        group = "Existing Greenspaces Reached",
        label = ~display_id,
        labelOptions = labelOptions(
          noHide = TRUE,
          direction = "center",
          textOnly = TRUE,
          style = list(
            "font-weight" = "bold",
            "font-size" = "12px",
            "color" = "#7f2704",
            "background" = "rgba(255,255,255,0.88)",
            "border" = "1px solid #e65100",
            "padding" = "1px 4px"
          )
        )
      )
    
    if (!is.null(connectivity$links_sf) && nrow(connectivity$links_sf) > 0) {
      links_map <- st_transform(connectivity$links_sf, 4326)
      proxy <- proxy |>
        addPolylines(
          data = links_map,
          group = "Measured Connection Gaps",
          color = "#ef6c00",
          weight = 3,
          opacity = 0.9,
          dashArray = "6,7",
          label = ~paste0(
            display_id,
            ifelse(distance_m <= 0.01, ": direct overlap", paste0(": ", round(distance_m), " m gap"))
          )
        )
    }
  })
  
  # At city scale, a 2--40 m buffer can be only a few screen pixels wide.
  # This action centers the map on the middle of a drawn line at street scale so
  # the selected total planted width can be inspected directly.
  observeEvent(input$inspect_corridor_width, {
    proposal <- proposed_greenspace()
    if (is.null(proposal) || nrow(proposal) == 0) {
      showNotification("Draw a line first, then use this button to inspect its planted width.", type = "warning")
      return(invisible(NULL))
    }
    
    proposal_3857 <- tryCatch(
      proposal |> st_transform(3857) |> st_zm(drop = TRUE, what = "ZM"),
      error = function(e) NULL
    )
    if (is.null(proposal_3857) || nrow(proposal_3857) == 0) return(invisible(NULL))
    
    types <- as.character(st_geometry_type(proposal_3857))
    line_idx <- types %in% c("LINESTRING", "MULTILINESTRING")
    if (!any(line_idx)) {
      showNotification(
        "The width control applies only to line drawings. Your polygon already defines its own area.",
        type = "message"
      )
      return(invisible(NULL))
    }
    
    line_geom <- suppressWarnings(st_union(st_geometry(proposal_3857[line_idx, , drop = FALSE])))
    inspect_point <- tryCatch(
      st_line_sample(line_geom, sample = 0.5),
      error = function(e) st_centroid(line_geom)
    )
    if (length(inspect_point) == 0 || all(st_is_empty(inspect_point))) {
      inspect_point <- st_centroid(line_geom)
    }
    
    inspect_sf <- st_sf(geometry = inspect_point, crs = 3857) |> st_transform(4326)
    coords <- st_coordinates(inspect_sf)
    if (nrow(coords) == 0) return(invisible(NULL))
    
    width_m <- suppressWarnings(as.numeric(input$corridor_width_m))
    if (!is.finite(width_m) || width_m <= 0) width_m <- 10
    
    leafletProxy("corridorMap") |>
      showGroup("Your Drawn Route or Shape") |>
      showGroup("Your Proposed Habitat Area") |>
      setView(lng = coords[1, "X"], lat = coords[1, "Y"], zoom = 19)
    
    showNotification(
      paste0(
        "Street-scale view: the blue footprint is ", width_m,
        " m wide in total (", width_m / 2, " m on each side)."
      ),
      type = "message",
      duration = 5
    )
    invisible(NULL)
  })
  
  # One stakeholder-facing action produces a clean corridor view and zooms to
  # the proposal plus the existing patches it affects.
  observeEvent(input$show_corridor_result, {
    footprint <- proposed_habitat_3857()
    connectivity <- corridor_connectivity()
    if (is.null(footprint) || nrow(footprint) == 0) {
      showNotification("Draw a corridor line or habitat polygon first.", type = "warning")
      return()
    }
    
    proxy <- leafletProxy("corridorMap") |>
      showGroup("Existing SF Greenspace") |>
      showGroup("Your Drawn Route or Shape") |>
      showGroup("Your Proposed Habitat Area") |>
      showGroup("Existing Greenspaces Reached") |>
      showGroup("Measured Connection Gaps")
    
    focus_parts <- list(st_geometry(footprint))
    if (!is.null(connectivity$linked_sf) && nrow(connectivity$linked_sf) > 0) {
      focus_parts[[length(focus_parts) + 1]] <- st_geometry(connectivity$linked_sf)
    }
    focus_geom <- do.call(c, focus_parts)
    focus_sf <- st_sf(geometry = focus_geom, crs = 3857) |> st_transform(4326)
    bb <- st_bbox(focus_sf)
    pad_x <- max((bb[["xmax"]] - bb[["xmin"]]) * 0.15, 0.002)
    pad_y <- max((bb[["ymax"]] - bb[["ymin"]]) * 0.15, 0.002)
    proxy |> fitBounds(
      lng1 = bb[["xmin"]] - pad_x,
      lat1 = bb[["ymin"]] - pad_y,
      lng2 = bb[["xmax"]] + pad_x,
      lat2 = bb[["ymax"]] + pad_y
    )
  })
  
  # Straight-line walking catchment around the novel proposed greenspace. This is
  # intentionally labelled as an approximation: it uses the app's existing
  # 80 m/min walking assumption rather than a routed pedestrian network.
  scenario_beneficiary_catchment_3857 <- reactive({
    proposal_net <- net_proposed_greenspace_3857()
    if (is.null(proposal_net) || nrow(proposal_net) == 0) return(NULL)
    
    walk_min <- suppressWarnings(as.numeric(input$beneficiary_walk_min))
    if (!is.finite(walk_min) || walk_min <= 0) return(NULL)
    
    catchment_geom <- st_geometry(proposal_net) |>
      st_union() |>
      st_buffer(dist = walk_min * walk_speed_m_per_min) |>
      st_make_valid()
    
    out <- st_sf(
      walk_minutes = walk_min,
      geometry = catchment_geom,
      crs = 3857
    )
    out <- out[!st_is_empty(out), , drop = FALSE]
    if (nrow(out) == 0) NULL else out
  })
  
  # Area-weighted beneficiary estimates. Population is apportioned by the share
  # of each census block group covered by the catchment, which assumes residents
  # are uniformly distributed within each block group. EJ beneficiaries are the
  # portion of that population also falling inside the SF EJ-community polygons.
  scenario_beneficiaries <- reactive({
    catchment <- scenario_beneficiary_catchment_3857()
    if (is.null(catchment) || nrow(catchment) == 0) {
      proposal <- proposed_greenspace()
      msg <- if (!is.null(proposal) && nrow(proposal) > 0) {
        "No novel proposed area remains after overlap with existing OSM greenspace is removed."
      } else {
        "Draw proposed greenspace to estimate beneficiaries."
      }
      return(list(available = FALSE, message = msg))
    }
    
    required_cbg_cols <- c("GEOID", "popE", "medincE")
    if (!inherits(cbg_vect_sf, "sf") || !all(required_cbg_cols %in% names(cbg_vect_sf))) {
      return(list(
        available = FALSE,
        message = "Beneficiary estimates require GEOID, popE, and medincE in cbg_vect_sf."
      ))
    }
    
    cbg_3857 <- cbg_vect_sf |>
      dplyr::select(dplyr::all_of(required_cbg_cols)) |>
      st_transform(3857) |>
      st_make_valid() |>
      mutate(cbg_area_m2 = as.numeric(st_area(geometry)))
    
    catchment_geom <- st_geometry(catchment)
    cbg_inter <- tryCatch(
      suppressWarnings(st_intersection(cbg_3857, catchment_geom)),
      error = function(e) NULL
    )
    
    if (is.null(cbg_inter) || nrow(cbg_inter) == 0) {
      return(list(
        available = TRUE,
        walk_minutes = catchment$walk_minutes[[1]],
        walk_distance_m = catchment$walk_minutes[[1]] * walk_speed_m_per_min,
        catchment_area_km2 = as.numeric(st_area(catchment_geom)) / 1e6,
        total_population = 0,
        ej_population = 0,
        ej_share = NA_real_,
        ej_vulnerability_threshold = NA_real_,
        ej_vulnerability_percentile = 75,
        low_access_population = 0,
        low_access_threshold = median(city_benchmarks$greenspace_cover, na.rm = TRUE),
        weighted_income = NA_real_,
        mean_ej_score = NA_real_,
        n_cbgs = 0
      ))
    }
    
    cbg_inter <- cbg_inter |>
      mutate(
        inter_area_m2 = as.numeric(st_area(geometry)),
        area_fraction = pmin(1, pmax(0, inter_area_m2 / cbg_area_m2)),
        served_population = pmax(0, as.numeric(popE)) * area_fraction
      )
    
    # Attach existing CBG greenspace cover so beneficiaries in relatively
    # low-greenspace areas can be reported. The threshold is the citywide median.
    low_access_threshold <- median(city_benchmarks$greenspace_cover, na.rm = TRUE)
    if (exists("cbg_greenspace_coverage") &&
        !is.null(cbg_greenspace_coverage) &&
        all(c("GEOID", "greenspace_m2", "cbg_area_m2") %in% names(cbg_greenspace_coverage))) {
      gs_lookup <- cbg_greenspace_coverage |>
        transmute(
          GEOID = as.character(GEOID),
          existing_greenspace_pct = ifelse(
            cbg_area_m2 > 0,
            100 * greenspace_m2 / cbg_area_m2,
            NA_real_
          )
        )
      cbg_inter <- cbg_inter |>
        mutate(GEOID = as.character(GEOID)) |>
        left_join(gs_lookup, by = "GEOID")
    } else {
      cbg_inter$existing_greenspace_pct <- NA_real_
    }
    
    total_population <- sum(cbg_inter$served_population, na.rm = TRUE)
    weighted_income <- if (total_population > 0) {
      sum(cbg_inter$medincE * cbg_inter$served_population, na.rm = TRUE) /
        sum(cbg_inter$served_population[is.finite(cbg_inter$medincE)], na.rm = TRUE)
    } else {
      NA_real_
    }
    if (!is.finite(weighted_income)) weighted_income <- NA_real_
    
    low_access_population <- if (is.finite(low_access_threshold)) {
      sum(
        cbg_inter$served_population[
          is.finite(cbg_inter$existing_greenspace_pct) &
            cbg_inter$existing_greenspace_pct < low_access_threshold
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
    
    ej_population <- 0
    mean_ej_score <- NA_real_
    ej_vulnerability_threshold <- NA_real_
    ej_vulnerability_percentile <- 75
    
    if (!is.null(sf_ej_sf) && inherits(sf_ej_sf, "sf") && nrow(sf_ej_sf) > 0) {
      ej_all_3857 <- sf_ej_sf |>
        st_transform(3857) |>
        st_make_valid()
      
      if ("score" %in% names(ej_all_3857)) {
        city_ej_scores <- suppressWarnings(as.numeric(ej_all_3857$score))
        city_ej_scores <- city_ej_scores[is.finite(city_ej_scores)]
        
        if (length(city_ej_scores) > 0) {
          ej_vulnerability_threshold <- as.numeric(
            stats::quantile(
              city_ej_scores,
              probs = 0.75,
              na.rm = TRUE,
              names = FALSE,
              type = 7
            )
          )
        }
        
        ej_vulnerable_3857 <- ej_all_3857 |>
          dplyr::filter(
            is.finite(score),
            score > 0,
            is.finite(ej_vulnerability_threshold),
            score >= ej_vulnerability_threshold
          )
      } else {
        ej_vulnerable_3857 <- ej_all_3857[0, , drop = FALSE]
      }
      
      vulnerable_union <- if (nrow(ej_vulnerable_3857) > 0) {
        st_union(st_geometry(ej_vulnerable_3857))
      } else {
        st_sfc(crs = 3857)
      }
      
      vulnerable_catchment <- tryCatch(
        suppressWarnings(st_intersection(catchment_geom, vulnerable_union)),
        error = function(e) st_sfc(crs = 3857)
      )
      
      if (length(vulnerable_catchment) > 0 &&
          !all(st_is_empty(vulnerable_catchment))) {
        ej_cbg_inter <- tryCatch(
          suppressWarnings(st_intersection(cbg_3857, vulnerable_catchment)),
          error = function(e) NULL
        )
        
        if (!is.null(ej_cbg_inter) && nrow(ej_cbg_inter) > 0) {
          ej_cbg_inter <- ej_cbg_inter |>
            mutate(
              inter_area_m2 = as.numeric(st_area(geometry)),
              area_fraction = pmin(1, pmax(0, inter_area_m2 / cbg_area_m2)),
              served_population = pmax(0, as.numeric(popE)) * area_fraction
            )
          ej_population <- sum(ej_cbg_inter$served_population, na.rm = TRUE)
        }
      }
      
      if ("score" %in% names(ej_all_3857)) {
        ej_score_inter <- tryCatch(
          suppressWarnings(st_intersection(ej_all_3857[, "score"], catchment_geom)),
          error = function(e) NULL
        )
        if (!is.null(ej_score_inter) && nrow(ej_score_inter) > 0) {
          ej_score_inter$area_m2 <- as.numeric(st_area(ej_score_inter))
          valid <- is.finite(ej_score_inter$score) & ej_score_inter$area_m2 > 0
          if (any(valid)) {
            mean_ej_score <- weighted.mean(
              ej_score_inter$score[valid],
              ej_score_inter$area_m2[valid],
              na.rm = TRUE
            )
          }
        }
      }
    }
    
    list(
      available = TRUE,
      walk_minutes = catchment$walk_minutes[[1]],
      walk_distance_m = catchment$walk_minutes[[1]] * walk_speed_m_per_min,
      catchment_area_km2 = as.numeric(st_area(catchment_geom)) / 1e6,
      total_population = total_population,
      ej_population = min(ej_population, total_population),
      ej_share = if (total_population > 0) 100 * min(ej_population, total_population) / total_population else NA_real_,
      ej_vulnerability_threshold = ej_vulnerability_threshold,
      ej_vulnerability_percentile = ej_vulnerability_percentile,
      low_access_population = low_access_population,
      low_access_threshold = low_access_threshold,
      weighted_income = weighted_income,
      mean_ej_score = mean_ej_score,
      n_cbgs = dplyr::n_distinct(cbg_inter$GEOID)
    )
  })
  
  # Draw the beneficiary catchment as a separate, toggleable map layer.
  observe({
    catchment <- scenario_beneficiary_catchment_3857()
    proxy <- leafletProxy("corridorMap") |>
      clearGroup("Population Estimate Area (Optional)")
    
    if (isTRUE(input$show_beneficiary_catchment) &&
        !is.null(catchment) && nrow(catchment) > 0) {
      catchment_4326 <- st_transform(catchment, 4326)
      proxy |>
        addPolygons(
          data = catchment_4326,
          group = "Population Estimate Area (Optional)",
          color = "#6a1b9a",
          weight = 2,
          dashArray = "6, 4",
          fillColor = "#ce93d8",
          fillOpacity = 0.14,
          label = paste0(
            catchment$walk_minutes[[1]],
            "-minute population-estimate area (approximate)"
          )
        )
    }
  })
  
  # ---------------------------------------------------------------------------
  # Location selection
  # ---------------------------------------------------------------------------
  observeEvent(input$isoMap_click, {
    req(input$location_choice == "map_click")
    req(!isTRUE(greenspace_draw_active()))
    click <- input$isoMap_click
    if (!is.null(click)) {
      chosen_point(c(lon = click$lng, lat = click$lat))
      leafletProxy("isoMap") |>
        clearGroup("selected_point") |>
        addCircleMarkers(
          lng = click$lng, lat = click$lat,
          radius = 6, color = "firebrick",
          group = "selected_point", label = "Map Click Location"
        )
    }
  })
  
  observeEvent(input$geocoder, {
    req(input$location_choice == "address")
    geocode_result <- input$geocoder
    if (!is.null(geocode_result)) {
      xy <- geocoder_as_xy(geocode_result)
      chosen_point(c(lon = xy[1], lat = xy[2]))
      leafletProxy("isoMap") |>
        clearGroup("selected_point") |>
        addCircleMarkers(
          lng = xy[1], lat = xy[2],
          radius = 6, color = "navy",
          group = "selected_point", label = "Geocoded Address"
        ) |>
        flyTo(lng = xy[1], lat = xy[2], zoom = 13)
    }
  })
  
  observeEvent(input$clear_map, {
    chosen_point(NULL)
    iso_store$data <- NULL
    leafletProxy("isoMap") |>
      removeControl("ndvi_legend") |>
      clearGroup("selected_point") |>
      clearGroup("Isochrones") |>
      clearGroup("Transit Isochrones") |>
      clearGroup("NDVI Raster")
  })
  
  # ---------------------------------------------------------------------------
  # Isochrone generation
  # ---------------------------------------------------------------------------
  # iso_store holds the computed sf; only written inside observeEvent(input$generate_iso).
  # Nothing can write to it except the button press, so map clicks cannot trigger
  # Mapbox calls regardless of what other reactives read isochrones_data().
  iso_store       <- reactiveValues(data = NULL)
  isochrones_data <- reactive({ iso_store$data })
  
  # Flag for the conditionalPanel that reveals the spider plot once isochrones
  # exist. suspendWhenHidden = FALSE so the condition is evaluated while the box
  # is still hidden.
  output$iso_ready <- reactive({
    d <- isochrones_data()
    !is.null(d) && nrow(d) > 0
  })
  outputOptions(output, "iso_ready", suspendWhenHidden = FALSE)
  
  output$proposal_ready <- reactive({
    d <- proposed_greenspace()
    !is.null(d) && nrow(d) > 0
  })
  outputOptions(output, "proposal_ready", suspendWhenHidden = FALSE)
  
  observeEvent(input$generate_iso, {
    pt <- chosen_point()
    if (is.null(pt)) return()
    if (length(input$transport_modes) == 0) return()
    if (length(input$iso_times) == 0) return()
    
    leafletProxy("isoMap") |>
      removeControl("ndvi_legend") |>
      clearGroup("Isochrones") |>
      clearGroup("Transit Isochrones") |>
      clearGroup("NDVI Raster")
    
    # All isochrone construction lives in build_isochrones() (code/iso_metrics_AUDITED.R)
    # so the comparer tab builds isochrones the exact same way.
    iso_store$data <- build_isochrones(
      point                        = pt,
      modes                        = input$transport_modes,
      times                        = input$iso_times,
      transit_hour                 = input$transit_hour,
      walk_to_stop_min             = input$walk_to_stop_min,
      walk_from_stop_min           = input$walk_from_stop_min,
      transit_departure_window_min = input$transit_departure_window_min
    )
  })
  
  # ---------------------------------------------------------------------------
  # Render isochrones and NDVI
  # ---------------------------------------------------------------------------
  observeEvent(isochrones_data(), {
    iso_data <- isochrones_data()
    req(iso_data)
    
    leafletProxy("isoMap") |>
      removeControl("ndvi_legend") |>
      clearGroup("Isochrones") |>
      clearGroup("Transit Isochrones") |>
      clearGroup("NDVI Raster")
    
    standard_like_modes <- c("driving", "walking", "cycling", "driving-traffic", "walk_transit")
    standard_iso <- iso_data[iso_data$mode %in% standard_like_modes, ]
    
    if (nrow(standard_iso) > 0) {
      for (i in seq_len(nrow(standard_iso))) {
        poly_i <- standard_iso[i, ]
        mode_i <- as.character(poly_i$mode[[1]])
        time_i <- as.numeric(poly_i$time[[1]])
        mode_label <- pretty_mode(mode_i)
        
        is_walk_transit <- identical(mode_i, "walk_transit")
        line_weight <- if (is_walk_transit) 3 else 2
        dash_style  <- if (is_walk_transit) "5, 5" else NULL
        fill_alpha  <- if (is_walk_transit) 0.28 else 0.35
        
        color_i <- unname(mode_palette[mode_label])
        if (is.na(color_i) || length(color_i) == 0) color_i <- "#666666"
        
        popup_i <- if (is_walk_transit) {
          paste0(
            "<strong>Walk + Transit — ", time_i, " min</strong>",
            "<br><b>First-mile walk:</b> up to ", input$walk_to_stop_min, " min",
            "<br><b>Transit:</b> uses remaining budget after first-mile access",
            "<br><b>Departure flexibility window:</b> ", input$transit_departure_window_min, " min",
            "<br><b>Last-mile walk:</b> up to ", input$walk_from_stop_min, " min",
            "<br><small>Built from reachable transit stops plus walking-network access.</small>"
          )
        } else {
          paste0("<strong>", mode_label, " — ", time_i, " min</strong>")
        }
        
        leafletProxy("isoMap") |>
          addPolygons(
            data = poly_i,
            group = "Isochrones",
            color = color_i,
            weight = line_weight,
            dashArray = dash_style,
            fillOpacity = fill_alpha,
            label = paste0(mode_label, " — ", time_i, " mins"),
            popup = popup_i
          )
      }
    }
    
    transit_iso <- iso_data[iso_data$mode == "transit", ]
    if (nrow(transit_iso) > 0 && !is.null(gtfs_stops_sf)) {
      for (i in seq_len(nrow(transit_iso))) {
        poly_i <- transit_iso[i, ]
        time_i <- as.numeric(poly_i$time[[1]])
        
        n_stops_in <- tryCatch(nrow(st_intersection(gtfs_stops_sf, poly_i)), error = function(e) 0)
        area_km2   <- round(as.numeric(st_area(st_transform(poly_i, 3857))) / 1e6, 2)
        t_score    <- if (area_km2 > 0) round(n_stops_in / area_km2, 2) else NA_real_
        
        iso_wkt_t <- st_as_text(st_geometry(poly_i)[[1]])
        gbif_t <- tryCatch({
          gbif_tab |>
            filter(sql(glue("ST_Intersects(ST_GeomFromText(geom_wkt), ST_GeomFromText('{iso_wkt_t}'))"))) |>
            summarise(n_records = n(), n_species = n_distinct(species)) |>
            collect()
        }, error = function(e) NULL)
        n_gbif_rec <- if (!is.null(gbif_t) && nrow(gbif_t) > 0) gbif_t$n_records[[1]] else 0L
        n_gbif_sp  <- if (!is.null(gbif_t) && nrow(gbif_t) > 0) gbif_t$n_species[[1]] else 0L
        
        popup_html <- paste0(
          "<strong>Transit Isochrone — ", time_i, " min</strong>",
          "<br><b>Transit Access Score:</b> ", t_score, " stops/km²",
          "<br><b>Muni stops within:</b> ", n_stops_in,
          "<br><b>Area:</b> ", area_km2, " km²",
          "<br><b>GBIF records within:</b> ", n_gbif_rec,
          "<br><b>Unique species:</b> ", n_gbif_sp
        )
        
        leafletProxy("isoMap") |>
          addPolygons(
            data = poly_i,
            group = "Transit Isochrones",
            color = mode_palette[["Transit"]],
            weight = 2,
            dashArray = "6, 4",
            fillOpacity = 0.22,
            label = paste0("Transit ", time_i, " min"),
            popup = popup_html
          )
      }
    }
    
    if (exists("ndvi")) {
      iso_union <- st_union(iso_data)
      iso_union_vect <- vect(iso_union)
      
      ndvi_crop <- terra::crop(ndvi, iso_union_vect)
      ndvi_mask <- terra::mask(ndvi_crop, iso_union_vect)
      ndvi_vals <- values(ndvi_mask)
      ndvi_vals <- ndvi_vals[!is.na(ndvi_vals)]
      
      if (length(ndvi_vals) > 0) {
        ndvi_pal <- colorNumeric("YlGn", domain = range(ndvi_vals, na.rm = TRUE), na.color = "transparent")
        leafletProxy("isoMap") |>
          addRasterImage(x = ndvi_mask, colors = ndvi_pal, opacity = 0.7, project = TRUE, group = "NDVI Raster") |>
          addLegend(
            position = "bottomright", pal = ndvi_pal, values = ndvi_vals, title = "NDVI",
            layerId = "ndvi_legend"
          )
      }
    }
  })
  
  # ---------------------------------------------------------------------------
  # Per-isochrone summaries
  # ---------------------------------------------------------------------------
  base_socio_data <- reactive({
    iso_data <- isochrones_data()
    if (is.null(iso_data) || nrow(iso_data) == 0) return(data.frame())
    
    # compute_iso_metrics() preserves one row per generated isochrone. Immediately
    # score EACH row independently against the matching precomputed mode × time
    # distribution. No union/mean across generated isochrones enters percentile
    # calculation.
    raw_metrics <- quiet_sf_overlay(
      compute_iso_metrics(iso_data, isolate(chosen_point()), gbif_tab)
    )
    
    score_iso_metrics_against_reference(raw_metrics)
  })
  
  # Apply the planning scenario to each isochrone. The drawn area modifies only
  # Greenspace_percent; all observed biodiversity, transit, NDVI and equity values
  # remain unchanged. Each scenario row is then re-ranked against the SAME mode ×
  # time precomputed reference and its own per-isochrone composite is updated.
  socio_data <- reactive({
    base <- base_socio_data()
    if (is.null(base) || nrow(base) == 0) return(base)
    
    base$Baseline_Greenspace_percent <- base$Greenspace_percent
    base$Added_Greenspace_m2 <- 0
    base$Added_Greenspace_percent <- 0
    
    proposal_net <- net_proposed_greenspace_3857()
    if (is.null(proposal_net) || nrow(proposal_net) == 0) return(base)
    
    # compute_iso_metrics() returns a non-spatial data.frame. Use the original
    # isochrone sf object for geometry and keep the metrics table separate.
    iso_sf <- isochrones_data()
    if (is.null(iso_sf) || !inherits(iso_sf, "sf") || nrow(iso_sf) == 0) {
      return(base)
    }
    
    # compute_iso_metrics() is expected to preserve one row per isochrone in the
    # same order. Fail safely rather than calling an sf method on the metrics table.
    if (nrow(iso_sf) != nrow(base)) {
      warning(
        "Proposed-greenspace scenario skipped: metric rows (", nrow(base),
        ") do not match isochrone rows (", nrow(iso_sf), ")."
      )
      return(base)
    }
    
    # Greenspace is now a 100 m grid-supported metric in BOTH the live app and
    # Step 5. Therefore the scenario must modify the same support: the union of
    # grid cells whose centroids fall within each individual isochrone.
    if (!("MetricSupportArea_km2" %in% names(base))) {
      warning("Scenario skipped: MetricSupportArea_km2 is missing from harmonized metrics.")
      return(base)
    }
    support_area_m2 <- as.numeric(base$MetricSupportArea_km2) * 1e6
    proposal_geometry <- st_geometry(proposal_net)
    
    added_m2 <- vapply(seq_len(nrow(iso_sf)), function(i) {
      support_geom_i <- tryCatch(
        iso_metric_support_geometry(iso_sf[i, , drop = FALSE]),
        error = function(e) st_sfc(crs = 3857)
      )
      if (length(support_geom_i) == 0) return(0)
      proposal_geometry_i <- tryCatch(
        st_transform(st_sf(geometry = proposal_geometry), st_crs(support_geom_i)) |> st_geometry(),
        error = function(e) st_sfc(crs = st_crs(support_geom_i))
      )
      intersection_i <- tryCatch(
        suppressWarnings(st_intersection(support_geom_i, proposal_geometry_i)),
        error = function(e) st_sfc(crs = st_crs(support_geom_i))
      )
      if (length(intersection_i) == 0) return(0)
      sum(as.numeric(st_area(intersection_i)), na.rm = TRUE)
    }, numeric(1))
    
    added_pct <- ifelse(
      is.finite(support_area_m2) & support_area_m2 > 0,
      100 * added_m2 / support_area_m2,
      NA_real_
    )
    
    base$Added_Greenspace_m2 <- added_m2
    base$Added_Greenspace_percent <- added_pct
    base$Greenspace_percent <- ifelse(
      is.na(base$Baseline_Greenspace_percent),
      NA_real_,
      pmin(100, base$Baseline_Greenspace_percent + added_pct)
    )
    # compute_iso_metrics() stores several summary values as custom attributes.
    # Reattach non-percentile attributes, then RE-SCORE each scenario isochrone
    # independently so a changed greenspace value receives a new percentile
    # against its own matching mode × time reference distribution.
    for (nm in c(
      "closest_greenspace", "closest_greenspace_dist_m",
      "closest_rsf_program", "closest_rsf_program_dist_m"
    )) {
      attr(base, nm) <- attr(base_socio_data(), nm)
    }
    
    score_iso_metrics_against_reference(base)
  })
  
  # ---------------------------------------------------------------------------
  # Biodiversity Access Index — same-mode × same-time isochrone benchmarks
  # ---------------------------------------------------------------------------
  baseline_biodiversity_access_index <- reactive({
    base <- base_socio_data()
    if (is.null(base) || nrow(base) == 0) return(base)
    add_reference_bai(base)
  })
  
  biodiversity_access_index <- reactive({
    scenario <- socio_data()
    if (is.null(scenario) || nrow(scenario) == 0) return(scenario)
    add_reference_bai(scenario)
  })
  
  # Corridor access comparisons use exactly one active Explorer isochrone.
  # NEVER average different modes or travel times for the corridor scenario.
  # The selected row alone enters the baseline-versus-proposal BAI calculation.
  corridor_selected_isochrone_metrics <- reactive({
    info <- corridor_active_iso_info()
    if (is.null(info)) return(NULL)
    
    idx <- info$index
    base_raw <- base_socio_data()
    scenario_raw <- socio_data()
    base_bai <- baseline_biodiversity_access_index()
    scenario_bai <- biodiversity_access_index()
    
    objects <- list(base_raw, scenario_raw, base_bai, scenario_bai)
    if (any(vapply(objects, function(x) is.null(x) || nrow(x) < idx, logical(1)))) {
      return(NULL)
    }
    
    list(
      index = idx,
      key = info$key,
      label = info$label,
      baseline_raw = base_raw[idx, , drop = FALSE],
      scenario_raw = scenario_raw[idx, , drop = FALSE],
      baseline_bai = base_bai[idx, , drop = FALSE],
      scenario_bai = scenario_bai[idx, , drop = FALSE]
    )
  })
  
  greenspace_scenario_comparison <- reactive({
    selected <- corridor_selected_isochrone_metrics()
    proposal <- proposed_greenspace()
    
    if (is.null(proposal) || nrow(proposal) == 0 || is.null(selected)) {
      return(NULL)
    }
    
    baseline <- selected$baseline_bai
    scenario <- selected$scenario_bai
    axis_cols <- c(
      "Mobility_Access_std", "Route_Access_std",
      "Biodiversity_Potential_std", "Observation_Intensity_std",
      "Environmental_Quality_std", "Greenspace_Cover_std",
      "Equity_Context_std"
    )
    if (!all(c(axis_cols, "BAI") %in% names(baseline)) ||
        !all(c(axis_cols, "BAI") %in% names(scenario))) {
      return(NULL)
    }
    
    list(
      baseline = baseline,
      scenario = scenario,
      baseline_raw = selected$baseline_raw,
      scenario_raw = selected$scenario_raw,
      axis_cols = axis_cols,
      index = selected$index,
      label = selected$label
    )
  })
  
  output$corridorImpactSummary <- renderUI({
    comparison <- greenspace_scenario_comparison()
    connectivity <- corridor_connectivity()
    
    if (is.null(comparison)) return(NULL)
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else x[[1]]
    }
    
    base_bai <- 100 * safe_num(comparison$baseline$BAI)
    scenario_bai <- 100 * safe_num(comparison$scenario$BAI)
    delta_bai <- scenario_bai - base_bai
    
    base_bai_pctile <- safe_num(comparison$baseline$pctile_BAI)
    scenario_bai_pctile <- safe_num(comparison$scenario$pctile_BAI)
    delta_bai_pctile <- scenario_bai_pctile - base_bai_pctile
    
    base_gs_pctile <- 100 * safe_num(comparison$baseline$Greenspace_Cover_std)
    scenario_gs_pctile <- 100 * safe_num(comparison$scenario$Greenspace_Cover_std)
    delta_gs_pctile <- scenario_gs_pctile - base_gs_pctile
    
    base_gs_raw <- safe_num(comparison$baseline_raw$Greenspace_percent)
    scenario_gs_raw <- safe_num(comparison$scenario_raw$Greenspace_percent)
    delta_gs_raw <- scenario_gs_raw - base_gs_raw
    
    # Identify which validated BAI axes actually changed in this modeled scenario.
    bai_axis_cols <- c(
      "Mobility_Access_std",
      "Route_Access_std",
      "Biodiversity_Potential_std",
      "Observation_Intensity_std",
      "Environmental_Quality_std",
      "Greenspace_Cover_std",
      "Equity_Context_std"
    )
    axis_delta <- vapply(bai_axis_cols, function(nm) {
      if (!(nm %in% names(comparison$baseline)) ||
          !(nm %in% names(comparison$scenario))) return(NA_real_)
      100 * (
        safe_num(comparison$scenario[[nm]]) -
          safe_num(comparison$baseline[[nm]])
      )
    }, numeric(1))
    changed_axes <- names(axis_delta)[is.finite(axis_delta) & abs(axis_delta) >= 0.1]
    
    greenspace_only <- identical(changed_axes, "Greenspace_Cover_std")
    
    impact_banner <- if (greenspace_only) {
      tags$div(
        style = paste0(
          "background:#e8f5e9; border-left:5px solid #2f7d4b; ",
          "border-radius:7px; padding:9px 12px; margin:8px 0 10px 0;"
        ),
        tags$b("Modeled BAI gain is driven by greenspace cover only. "),
        "The other six validated BAI dimensions are held constant in this scenario."
      )
    } else if (length(changed_axes) > 0) {
      tags$div(
        style = paste0(
          "background:#eef6f1; border-left:5px solid #5d8d70; ",
          "border-radius:7px; padding:9px 12px; margin:8px 0 10px 0;"
        ),
        tags$b("Modeled BAI change is concentrated in: "),
        paste(changed_axes, collapse = ", "),
        "."
      )
    } else {
      tags$div(
        style = paste0(
          "background:#f4f4f4; border-left:5px solid #888; ",
          "border-radius:7px; padding:9px 12px; margin:8px 0 10px 0;"
        ),
        tags$b("No validated BAI dimension changed under the current scenario settings.")
      )
    }
    
    fmt_arrow <- function(a, b, suffix = "", digits = 1) {
      if (!is.finite(a) || !is.finite(b)) return("N/A")
      paste0(round(a, digits), suffix, " → ", round(b, digits), suffix)
    }
    fmt_delta <- function(x, suffix = "", digits = 1) {
      if (!is.finite(x)) return("N/A")
      paste0(if (x >= 0) "+" else "", round(x, digits), suffix)
    }
    
    connectivity_value <- "Unavailable"
    connectivity_detail <- "Structural connectivity could not be evaluated."
    connectivity_subdetail <- NULL
    
    if (!is.null(connectivity) && isTRUE(connectivity$available)) {
      connectivity_value <- connectivity$classification
      
      if (identical(connectivity$classification, "Network reinforcement") &&
          is.finite(connectivity$baseline_components) &&
          is.finite(connectivity$scenario_components)) {
        connectivity_detail <- paste0(
          "The habitat network remains ",
          connectivity$baseline_components, " → ", connectivity$scenario_components,
          "; the proposal reinforces habitat around ",
          connectivity$linked_patches, " reached greenspaces."
        )
      } else {
        connectivity_detail <- paste0(
          connectivity$baseline_components, " → ", connectivity$scenario_components,
          " habitat networks; ", connectivity$linked_patches, " greenspaces reached."
        )
      }
      
      connectivity_subdetail <- paste0(
        "Separate structural-connectivity result — not included as a BAI axis."
      )
    }
    
    connectivity_sentence <- if (!is.null(connectivity) && isTRUE(connectivity$available)) {
      dplyr::case_when(
        identical(connectivity$classification, "Network reinforcement") ~
          paste0(
            "It reinforces an existing habitat network rather than joining separate networks under the ",
            connectivity$gap_m, " m connection rule."
          ),
        identical(connectivity$classification, "Bridge") ~
          paste0(
            "It joins previously separate habitat networks under the ",
            connectivity$gap_m, " m connection rule."
          ),
        identical(connectivity$classification, "Patch expansion") ~
          "It expands an existing mapped habitat patch rather than creating a new network connection.",
        identical(connectivity$classification, "Stepping-stone") ~
          "It adds a potential stepping stone but does not join separate habitat networks.",
        TRUE ~
          "It remains structurally isolated from mapped greenspaces under the selected connection rule."
      )
    } else {
      "Structural connectivity is unavailable for the current proposal."
    }
    
    stakeholder_summary <- if (
      is.finite(base_gs_raw) && is.finite(scenario_gs_raw) &&
      is.finite(base_gs_pctile) && is.finite(scenario_gs_pctile) &&
      is.finite(base_bai) && is.finite(scenario_bai)
    ) {
      paste0(
        "This intervention increases accessible greenspace from ",
        round(base_gs_raw, 1), "% to ", round(scenario_gs_raw, 1),
        "%, moving the greenspace benchmark from ",
        round(base_gs_pctile, 1), " to ", round(scenario_gs_pctile, 1),
        " percentile points and the composite BAI from ",
        round(base_bai, 1), " to ", round(scenario_bai, 1), "/100. ",
        connectivity_sentence
      )
    } else {
      paste0("The modeled intervention changes accessible greenspace and the selected-isochrone BAI. ", connectivity_sentence)
    }
    
    tags$div(
      style = "background:#f7fbf8; border:1px solid #cfe2d6; border-radius:10px; padding:12px 12px 10px 12px; margin-bottom:10px;",
      
      tags$div(
        style = "display:flex; justify-content:space-between; align-items:flex-start; gap:10px; flex-wrap:wrap;",
        tags$div(
          tags$div(
            style = "font-size:19px; font-weight:700; color:#28543b;",
            "Conservation action impact dashboard"
          ),
          tags$small(
            "Two different outcomes are shown together: validated BAI change and structural connectivity."
          )
        ),
        tags$span(
          style = "background:#e8f5e9; border:1px solid #c6dfca; border-radius:999px; padding:5px 10px; font-weight:700; color:#2f6b45;",
          comparison$label
        )
      ),
      
      tags$div(
        style = "background:#ffffff; border-left:5px solid #2f855a; border-radius:8px; padding:10px 12px; margin:9px 0;",
        tags$b("Stakeholder summary: "),
        stakeholder_summary
      ),
      
      impact_banner,
      
      tags$div(
        style = "display:grid; grid-template-columns:3fr 2fr; gap:10px; margin-bottom:9px;",
        
        tags$div(
          style = "background:white; border:1px solid #cfe2d6; border-radius:9px; padding:10px;",
          tags$div(
            style = "font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:#557064; font-weight:700;",
            "BAI IMPACT — validated 7-axis accessibility index"
          ),
          tags$div(
            style = "display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:8px; margin-top:7px;",
            tags$div(
              style = "background:#f8fbf9; border:1px solid #e0ebe4; border-radius:8px; padding:9px;",
              tags$small(tags$b("Composite BAI")),
              tags$div(style = "font-size:21px; font-weight:700;", fmt_arrow(base_bai, scenario_bai)),
              tags$small(
                style = if (is.finite(delta_bai) && delta_bai > 0) "color:#2f7d4b; font-weight:700;" else "color:#66766e;",
                paste0("Change: ", fmt_delta(delta_bai))
              )
            ),
            tags$div(
              style = "background:#f3f7fb; border:1px solid #d9e4ef; border-radius:8px; padding:9px;",
              tags$small(tags$b("BAI citywide percentile")),
              tags$div(style = "font-size:21px; font-weight:700;", fmt_arrow(base_bai_pctile, scenario_bai_pctile)),
              tags$small(
                style = if (is.finite(delta_bai_pctile) && delta_bai_pctile > 0) "color:#2f7d4b; font-weight:700;" else "color:#66766e;",
                paste0("Change: ", fmt_delta(delta_bai_pctile, " pp"))
              )
            ),
            tags$div(
              style = "background:#f8fbf9; border:1px solid #e0ebe4; border-radius:8px; padding:9px;",
              tags$small(tags$b("Greenspace percentile")),
              tags$div(style = "font-size:21px; font-weight:700;", fmt_arrow(base_gs_pctile, scenario_gs_pctile)),
              tags$small(
                style = if (is.finite(delta_gs_pctile) && delta_gs_pctile > 0) "color:#2f7d4b; font-weight:700;" else "color:#66766e;",
                paste0("Change: ", fmt_delta(delta_gs_pctile, " pp"))
              )
            ),
            tags$div(
              style = "background:#f8fbf9; border:1px solid #e0ebe4; border-radius:8px; padding:9px;",
              tags$small(tags$b("Greenspace cover (raw)")),
              tags$div(style = "font-size:21px; font-weight:700;", fmt_arrow(base_gs_raw, scenario_gs_raw, "%")),
              tags$small(
                style = if (is.finite(delta_gs_raw) && delta_gs_raw > 0) "color:#2f7d4b; font-weight:700;" else "color:#66766e;",
                paste0("Change: ", fmt_delta(delta_gs_raw, " pp"))
              )
            )
          )
        ),
        
        tags$div(
          style = "background:#fffaf0; border:1px solid #ecd8aa; border-radius:9px; padding:10px;",
          tags$div(
            style = "font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:#7a622b; font-weight:700;",
            "CONNECTIVITY OUTCOME — separate structural metric"
          ),
          tags$div(
            style = "font-size:23px; line-height:1.1; font-weight:800; color:#5f4a19; margin-top:8px;",
            connectivity_value
          ),
          tags$div(
            style = "font-size:13px; line-height:1.35; margin-top:6px; color:#4f493c;",
            connectivity_detail
          ),
          if (!is.null(connectivity_subdetail)) {
            tags$small(
              style = "display:block; margin-top:7px; color:#766b55;",
              connectivity_subdetail
            )
          }
        )
      )
    )
  })
  
  
  output$corridorBaiComparisonNote <- renderUI({
    info <- corridor_active_iso_info()
    if (is.null(info)) {
      return(tags$p("Choose an Explorer isochrone to compare baseline and proposal conditions."))
    }
    tags$p(
      tags$b("How to read this scenario: "), info$label, " is the only active access area. ",
      "The proposal changes greenspace cover inside that same isochrone and re-ranks it against the same citywide mode × time reference. ",
      "Observed species richness, NDVI, transit access, route access, sampling density, and equity context are held constant. ",
      "Therefore any outward movement in the radar identifies the modeled access-score effect of the planned greenspace action rather than a prediction of future biodiversity."
    )
  })
  
  
  # ---------------------------------------------------------------------------
  # HTML report cards
  # ---------------------------------------------------------------------------
  output$download_explorer_report <- downloadHandler(
    filename = function() paste0("sf_biodiversity_access_explorer_", format(Sys.Date(), "%Y%m%d"), ".html"),
    contentType = "text/html",
    content = function(file) {
      df <- baseline_biodiversity_access_index()
      if (is.null(df) || nrow(df) == 0) stop("Generate isochrones before downloading this report card.")
      
      colnum <- function(nm, mult = 1) {
        if (!(nm %in% names(df))) return(rep(NA_real_, nrow(df)))
        suppressWarnings(as.numeric(df[[nm]])) * mult
      }
      
      report_df <- data.frame(
        Isochrone = vapply(seq_len(nrow(df)), function(i) iso_metric_row_label(df, i), character(1)),
        `Transit raw (stops/km²)` = round(colnum("Transit_Access_Score"), 2),
        `Transit percentile` = round(colnum("pctile_Transit_Access_Score"), 1),
        `Muni routes raw` = round(colnum("Unique_Muni_Routes")),
        `Muni route percentile` = round(colnum("pctile_Unique_Muni_Routes"), 1),
        `Total species raw` = round(colnum("GBIF_Species")),
        `Species percentile` = round(colnum("pctile_GBIF_Species"), 1),
        `Sampling density raw` = round(colnum("SamplingDensity_km2"), 1),
        `Sampling-density percentile` = round(colnum("pctile_SamplingDensity_km2"), 1),
        `Mean NDVI raw` = round(colnum("MeanNDVI"), 3),
        `NDVI percentile` = round(colnum("pctile_MeanNDVI"), 1),
        `Greenspace cover raw (%)` = round(colnum("Greenspace_percent"), 1),
        `Greenspace percentile` = round(colnum("pctile_Greenspace_percent"), 1),
        `SF EJ score raw` = round(colnum("SF_EJ_Score"), 1),
        `Favorable equity percentile` = round(colnum("pctile_SF_EJ_Score"), 1),
        `BAI score (/100)` = round(colnum("BAI", 100), 1),
        `BAI citywide percentile` = round(colnum("pctile_BAI"), 1),
        `GBIF records` = round(colnum("GBIF_Records")),
        check.names = FALSE
      )
      
      map_iso <- isochrones_data()
      map_labels <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        vapply(seq_len(nrow(map_iso)), function(i) {
          paste0(pretty_mode(as.character(map_iso$mode[[i]])), " ", as.numeric(map_iso$time[[i]]), " min")
        }, character(1))
      } else character(0)
      map_colors <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        vapply(seq_len(nrow(map_iso)), function(i) {
          clr <- unname(mode_palette[pretty_mode(as.character(map_iso$mode[[i]]))])
          if (length(clr) == 0 || is.na(clr) || !nzchar(clr)) "#666666" else clr
        }, character(1))
      } else character(0)
      map_uri <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        report_card_map_uri(function(path) {
          write_report_isochrone_map(
            path,
            map_iso,
            labels = map_labels,
            colors = map_colors,
            map_title = "Isochrone Explorer — spatial context"
          )
        })
      } else NULL
      
      write_report_card(
        file,
        "Isochrone Explorer — Biodiversity Access Report Card",
        "Citywide same-mode × same-time benchmarking",
        sections = paste0(
          report_card_map_section(map_uri, "Map of Isochrone Explorer analysis areas"),
          '<section class="section"><h2>Results</h2>',
          '<div class="note">Every BAI dimension is shown as a raw value paired with its citywide percentile. BAI score /100 is the seven-axis mean; BAI citywide percentile is a separate rank of that composite score against complete same-mode × same-time reference BAIs.</div>',
          report_card_table(report_df),
          '</section>'
        )
      )
    }
  )
  
  output$download_comparer_report <- downloadHandler(
    filename = function() paste0("sf_biodiversity_access_comparison_", format(Sys.Date(), "%Y%m%d"), ".html"),
    contentType = "text/html",
    content = function(file) {
      res <- cmp_results()
      if (is.null(res)) stop("Run the Point A vs Point B comparison first.")
      
      side_row <- function(side, side_name) {
        if (is.null(side) || is.null(side$metrics) || nrow(side$metrics) == 0) {
          return(data.frame(Point = side_name, Location = "Unavailable", check.names = FALSE))
        }
        m <- side$metrics
        b <- side$bai
        
        gm <- function(nm) if (nm %in% names(m)) suppressWarnings(as.numeric(m[[nm]][[1]])) else NA_real_
        gb <- function(nm) if (!is.null(b) && nm %in% names(b)) suppressWarnings(as.numeric(b[[nm]][[1]])) else NA_real_
        
        data.frame(
          Point = side_name,
          Location = if (!is.null(side$label)) side$label else iso_metric_row_label(m, 1),
          `Transit raw (stops/km²)` = round(gm("Transit_Access_Score"), 2),
          `Transit percentile` = round(gm("pctile_Transit_Access_Score"), 1),
          `Muni routes raw` = round(gm("Unique_Muni_Routes")),
          `Muni route percentile` = round(gm("pctile_Unique_Muni_Routes"), 1),
          `Total species raw` = round(gm("GBIF_Species")),
          `Species percentile` = round(gm("pctile_GBIF_Species"), 1),
          `Sampling density raw` = round(gm("SamplingDensity_km2"), 1),
          `Sampling-density percentile` = round(gm("pctile_SamplingDensity_km2"), 1),
          `Mean NDVI raw` = round(gm("MeanNDVI"), 3),
          `NDVI percentile` = round(gm("pctile_MeanNDVI"), 1),
          `Greenspace raw (%)` = round(gm("Greenspace_percent"), 1),
          `Greenspace percentile` = round(gm("pctile_Greenspace_percent"), 1),
          `SF EJ raw` = round(gm("SF_EJ_Score"), 1),
          `Favorable equity percentile` = round(gm("pctile_SF_EJ_Score"), 1),
          `BAI score (/100)` = round(100 * gb("BAI"), 1),
          `BAI citywide percentile` = round(gb("pctile_BAI"), 1),
          check.names = FALSE
        )
      }
      
      report_df <- dplyr::bind_rows(side_row(res$a, "Point A"), side_row(res$b, "Point B"))
      
      map_iso <- list()
      map_labels <- character()
      map_colors <- character()
      if (!is.null(res$a) && !is.null(res$a$iso) && nrow(res$a$iso) > 0) {
        map_iso[[length(map_iso) + 1L]] <- res$a$iso
        map_labels <- c(map_labels, paste0("Point A — ", res$a$label))
        map_colors <- c(map_colors, cmp_color_a)
      }
      if (!is.null(res$b) && !is.null(res$b$iso) && nrow(res$b$iso) > 0) {
        map_iso[[length(map_iso) + 1L]] <- res$b$iso
        map_labels <- c(map_labels, paste0("Point B — ", res$b$label))
        map_colors <- c(map_colors, cmp_color_b)
      }
      map_uri <- if (length(map_iso) > 0) {
        report_card_map_uri(function(path) {
          write_report_isochrone_map(
            path,
            map_iso,
            labels = map_labels,
            colors = map_colors,
            map_title = "Isochrone Comparer — spatial context"
          )
        })
      } else NULL
      
      write_report_card(
        file,
        "Isochrone Comparer — Report Card",
        "Point A and Point B are benchmarked independently",
        sections = paste0(
          report_card_map_section(map_uri, "Map of Point A and Point B isochrone analysis areas"),
          '<section class="section"><h2>Comparison</h2>',
          '<div class="note">Each point uses its own matching SF transportation-mode × travel-time reference. Raw values and percentiles remain paired, and the two BAI quantities are reported separately.</div>',
          report_card_table(report_df),
          '</section>'
        )
      )
    }
  )
  
  output$download_displacement_report <- downloadHandler(
    filename = function() paste0("sf_green_investment_displacement_context_", format(Sys.Date(), "%Y%m%d"), ".html"),
    contentType = "text/html",
    content = function(file) {
      df <- baseline_biodiversity_access_index()
      if (is.null(df) || nrow(df) == 0) stop("Generate isochrones before downloading this context report card.")
      
      colnum <- function(nm, mult = 1) {
        if (!(nm %in% names(df))) return(rep(NA_real_, nrow(df)))
        suppressWarnings(as.numeric(df[[nm]])) * mult
      }
      gs <- colnum("pctile_Greenspace_percent")
      eq <- colnum("pctile_SF_EJ_Score")
      
      report_df <- data.frame(
        Isochrone = vapply(seq_len(nrow(df)), function(i) iso_metric_row_label(df, i), character(1)),
        `Greenspace cover raw (%)` = round(colnum("Greenspace_percent"), 1),
        `Greenspace percentile` = round(gs, 1),
        `Relative greenspace-access gap` = round(100 - gs, 1),
        `SF EJ score raw` = round(colnum("SF_EJ_Score"), 1),
        `Favorable equity-context percentile` = round(eq, 1),
        `EJ burden context` = round(100 - eq, 1),
        `Median-income context` = round(colnum("MedianIncome")),
        `BAI score (/100)` = round(colnum("BAI", 100), 1),
        `BAI citywide percentile` = round(colnum("pctile_BAI"), 1),
        check.names = FALSE
      )
      
      map_iso <- isochrones_data()
      map_labels <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        vapply(seq_len(nrow(map_iso)), function(i) {
          paste0(pretty_mode(as.character(map_iso$mode[[i]])), " ", as.numeric(map_iso$time[[i]]), " min")
        }, character(1))
      } else character(0)
      map_colors <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        vapply(seq_len(nrow(map_iso)), function(i) {
          clr <- unname(mode_palette[pretty_mode(as.character(map_iso$mode[[i]]))])
          if (length(clr) == 0 || is.na(clr) || !nzchar(clr)) "#666666" else clr
        }, character(1))
      } else character(0)
      map_uri <- if (!is.null(map_iso) && nrow(map_iso) > 0) {
        report_card_map_uri(function(path) {
          write_report_isochrone_map(
            path,
            map_iso,
            labels = map_labels,
            colors = map_colors,
            map_title = "Green Investment & Displacement Context — spatial context"
          )
        })
      } else NULL
      
      write_report_card(
        file,
        "Green Investment & Displacement Context — Report Card",
        "Equity-sensitive greening context; not a displacement-risk score",
        sections = paste0(
          report_card_map_section(map_uri, "Map of isochrone areas used for the green investment and displacement context report"),
          '<section class="section"><h2>Context</h2>',
          report_card_table(report_df),
          '</section>',
          '<section class="section"><h2>Interpretation</h2>',
          '<p>This screen does <strong>not</strong> calculate green-gentrification or displacement risk. A descriptive overlap between greenspace need and EJ burden is a reason to consider community governance and anti-displacement safeguards, not evidence that greening will cause displacement.</p>',
          '<p>A defensible displacement-risk screen would require housing tenure, rent burden, eviction/displacement events, housing-cost trends, affordable-housing protection, and related housing-market indicators.</p>',
          '</section>'
        ),
        caution = "Do not interpret this context screen as evidence that a proposed greening action will cause displacement."
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # Green Investment & Displacement Context
  # ---------------------------------------------------------------------------
  output$displacementContextSummary <- renderUI({
    df <- baseline_biodiversity_access_index()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tags$div(
        style = "color:#66766e;",
        "Generate one or more isochrones in the Explorer to populate this context screen."
      ))
    }
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else x[[1]]
    }
    
    cards <- lapply(seq_len(nrow(df)), function(i) {
      gs_raw <- if ("Greenspace_percent" %in% names(df)) safe_num(df$Greenspace_percent[[i]]) else NA_real_
      gs_pct <- if ("pctile_Greenspace_percent" %in% names(df)) safe_num(df$pctile_Greenspace_percent[[i]]) else NA_real_
      equity_pct <- if ("pctile_SF_EJ_Score" %in% names(df)) safe_num(df$pctile_SF_EJ_Score[[i]]) else NA_real_
      ej_raw <- if ("SF_EJ_Score" %in% names(df)) safe_num(df$SF_EJ_Score[[i]]) else NA_real_
      income <- if ("MedianIncome" %in% names(df)) safe_num(df$MedianIncome[[i]]) else NA_real_
      bai <- if ("BAI" %in% names(df)) 100 * safe_num(df$BAI[[i]]) else NA_real_
      bai_pct <- if ("pctile_BAI" %in% names(df)) safe_num(df$pctile_BAI[[i]]) else NA_real_
      
      greenspace_gap <- if (is.finite(gs_pct)) 100 - gs_pct else NA_real_
      ej_burden_context <- if (is.finite(equity_pct)) 100 - equity_pct else NA_real_
      
      planning_context <- if (!is.finite(greenspace_gap) || !is.finite(ej_burden_context)) {
        "Insufficient comparable percentile information for a descriptive planning context."
      } else if (greenspace_gap >= 50 && ej_burden_context >= 50) {
        "Greening opportunity overlaps with greater EJ burden context. Consider community-led planning and explicit anti-displacement safeguards alongside ecological investment."
      } else if (greenspace_gap >= 50) {
        "Larger relative greenspace-access gap, but lower EJ burden context on this screen. Greening need is visible; displacement vulnerability is not estimated."
      } else if (ej_burden_context >= 50) {
        "Greater EJ burden context despite comparatively stronger greenspace access. Equity considerations remain important, but this is not a displacement-risk classification."
      } else {
        "Neither the relative greenspace-access gap nor EJ burden context falls in the upper half of these descriptive citywide scales."
      }
      
      tags$div(
        style = paste0(
          "border:1px solid #dbe8e1; border-radius:9px; padding:10px 12px; ",
          "margin-bottom:9px; background:white;"
        ),
        tags$div(
          style = "font-size:15px; font-weight:700; color:#294337; margin-bottom:7px;",
          iso_metric_row_label(df, i)
        ),
        tags$div(
          style = "display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:8px;",
          tags$div(
            style = paste0("border-top:4px solid ", BAI_GROUP_COLORS[["green_environment"]], "; padding-top:5px;"),
            tags$small(tags$b("Greenspace cover")),
            tags$div(style = "font-size:19px; font-weight:700;", if (is.finite(gs_raw)) paste0(round(gs_raw, 1), "%") else "N/A"),
            tags$small("Citywide rank: ", format_pct(gs_pct))
          ),
          tags$div(
            style = paste0("border-top:4px solid ", BAI_GROUP_COLORS[["equity"]], "; padding-top:5px;"),
            tags$small(tags$b("SF EJ context")),
            tags$div(style = "font-size:19px; font-weight:700;", if (is.finite(ej_raw)) round(ej_raw, 1) else "N/A"),
            tags$small("Favorable citywide rank: ", format_pct(equity_pct))
          ),
          tags$div(
            tags$small(tags$b("Median-income context")),
            tags$div(
              style = "font-size:19px; font-weight:700;",
              if (is.finite(income)) scales::dollar(income, accuracy = 1) else "N/A"
            ),
            tags$small("Context only; not a vulnerability score")
          ),
          tags$div(
            tags$small(tags$b("BAI score")),
            tags$div(style = "font-size:19px; font-weight:700;", if (is.finite(bai)) paste0(round(bai, 1), "/100") else "N/A"),
            tags$small("Seven-axis mean")
          ),
          tags$div(
            tags$small(tags$b("BAI citywide rank")),
            tags$div(style = "font-size:19px; font-weight:700;", format_pct(bai_pct)),
            tags$small("Same mode × same travel time")
          )
        ),
        tags$div(
          style = paste0(
            "background:#fff7e8; border-left:4px solid #d6a542; ",
            "border-radius:6px; padding:7px 9px; margin-top:8px; font-size:12px;"
          ),
          tags$b("Planning interpretation — descriptive, not a risk classification: "),
          planning_context
        ),
        tags$div(
          style = paste0(
            "background:#f7f7f7; border-left:4px solid #888; border-radius:6px; ",
            "padding:7px 9px; margin-top:6px; font-size:12px;"
          ),
          tags$b("Displacement-risk estimate: not calculated. "),
          "Housing tenure, rent burden, eviction, housing-cost trend, and affordable-housing data are not currently included."
        )
      )
    })
    
    tags$div(cards)
  })
  
  output$displacementContextPlot <- renderPlot({
    df <- baseline_biodiversity_access_index()
    
    if (is.null(df) || nrow(df) == 0 ||
        !("pctile_Greenspace_percent" %in% names(df)) ||
        !("pctile_SF_EJ_Score" %in% names(df))) {
      plot.new()
      title("Generate supported isochrones to view the context plot.")
      return(invisible(NULL))
    }
    
    plot_df <- data.frame(
      label = vapply(
        seq_len(nrow(df)),
        function(i) iso_metric_row_label(df, i),
        character(1)
      ),
      greenspace_gap = 100 - suppressWarnings(as.numeric(df$pctile_Greenspace_percent)),
      ej_burden_context = 100 - suppressWarnings(as.numeric(df$pctile_SF_EJ_Score)),
      stringsAsFactors = FALSE
    )
    
    plot_df <- plot_df[
      is.finite(plot_df$greenspace_gap) &
        is.finite(plot_df$ej_burden_context),
      ,
      drop = FALSE
    ]
    
    if (nrow(plot_df) == 0) {
      plot.new()
      title("No comparable same-mode × same-time context values are available.")
      return(invisible(NULL))
    }
    
    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = greenspace_gap,
        y = ej_burden_context,
        label = label
      )
    ) +
      ggplot2::annotate(
        "rect",
        xmin = 50, xmax = 100,
        ymin = 50, ymax = 100,
        fill = "#fff4d6",
        alpha = 0.45
      ) +
      ggplot2::geom_vline(
        xintercept = 50,
        linetype = "dashed",
        color = "#9aa7a0"
      ) +
      ggplot2::geom_hline(
        yintercept = 50,
        linetype = "dashed",
        color = "#9aa7a0"
      ) +
      ggplot2::geom_point(
        size = 4.2,
        shape = 21,
        fill = "#5d8d70",
        color = "#294337",
        stroke = 1.1
      ) +
      ggplot2::geom_text(
        nudge_y = 3.2,
        size = 4,
        fontface = "bold",
        check_overlap = TRUE
      ) +
      ggplot2::annotate(
        "text",
        x = 75,
        y = 96,
        label = "High greening opportunity\\n+ greater EJ burden context",
        fontface = "bold",
        size = 4.2,
        color = "#795b16"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100)
      ) +
      ggplot2::labs(
        x = "Relative greenspace-access gap\\n100 − greenspace-cover percentile",
        y = "EJ burden context\\n100 − favorable equity-context percentile",
        title = "Equity-Sensitive Greening Context",
        subtitle = paste0(
          "Upper-right highlights overlap between a larger greenspace-access gap ",
          "and greater EJ burden context. It is not a displacement-risk classification."
        ),
        caption = paste0(
          "Each point is one independently scored isochrone. ",
          "The 50/50 guide lines are descriptive only."
        )
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "#e1e9e4"),
        plot.title = ggplot2::element_text(
          face = "bold",
          color = "#294337",
          size = 16
        ),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(color = "#596861"),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          color = "#66766e"
        )
      )
  })
  
  output$ngoReportingSummary <- renderUI({
    df <- baseline_biodiversity_access_index()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tags$p(
        "Generate an isochrone to create reporting-ready citywide benchmark language."
      ))
    }
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else x[[1]]
    }
    
    rows <- lapply(seq_len(nrow(df)), function(i) {
      species_raw <- if ("GBIF_Species" %in% names(df)) safe_num(df$GBIF_Species[[i]]) else NA_real_
      species_pct <- if ("pctile_GBIF_Species" %in% names(df)) safe_num(df$pctile_GBIF_Species[[i]]) else NA_real_
      greenspace_raw <- if ("Greenspace_percent" %in% names(df)) safe_num(df$Greenspace_percent[[i]]) else NA_real_
      greenspace_pct <- if ("pctile_Greenspace_percent" %in% names(df)) safe_num(df$pctile_Greenspace_percent[[i]]) else NA_real_
      equity_raw <- if ("SF_EJ_Score" %in% names(df)) safe_num(df$SF_EJ_Score[[i]]) else NA_real_
      equity_pct <- if ("pctile_SF_EJ_Score" %in% names(df)) safe_num(df$pctile_SF_EJ_Score[[i]]) else NA_real_
      bai <- if ("BAI" %in% names(df)) 100 * safe_num(df$BAI[[i]]) else NA_real_
      bai_pct <- if ("pctile_BAI" %in% names(df)) safe_num(df$pctile_BAI[[i]]) else NA_real_
      
      sentence_parts <- c(
        paste0(iso_metric_row_label(df, i), ":"),
        if (is.finite(species_raw) && is.finite(species_pct)) {
          paste0(
            scales::comma(round(species_raw)),
            " accessible unique GBIF species (", format_pct(species_pct), ")"
          )
        } else NULL,
        if (is.finite(greenspace_raw) && is.finite(greenspace_pct)) {
          paste0(
            round(greenspace_raw, 1),
            "% greenspace cover (", format_pct(greenspace_pct), ")"
          )
        } else NULL,
        if (is.finite(bai)) {
          paste0(
            "BAI score ", round(bai, 1), "/100",
            if (is.finite(bai_pct)) paste0(" with a citywide BAI rank of ", format_pct(bai_pct)) else ""
          )
        } else NULL,
        if (is.finite(equity_raw) && is.finite(equity_pct)) {
          paste0(
            "SF EJ raw score ", round(equity_raw, 1),
            " with a favorable equity-context rank of ", format_pct(equity_pct),
            " (higher rank = lower EJ burden)"
          )
        } else NULL
      )
      
      tags$div(
        style = paste0(
          "background:#f3f8f5; border-left:4px solid #5d8d70; ",
          "border-radius:7px; padding:9px 11px; margin-bottom:8px;"
        ),
        paste0(paste(sentence_parts, collapse = "; "), ".")
      )
    })
    
    tags$div(
      rows,
      tags$small(
        style = "display:block; color:#66766e; margin-top:7px;",
        "Recommended wording attributes each benchmark to the specific isochrone/site and matching transportation scenario—not to the NGO as an organization. BAI score and BAI citywide percentile are reported separately."
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # Summary outputs
  # ---------------------------------------------------------------------------
  output$dataTable <- renderDT({
    df <- base_socio_data()
    if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No isochrones generated yet.")))
    
    out <- df |>
      mutate(
        MedianIncome = ifelse(is.na(MedianIncome), NA, dollar(MedianIncome)),
        MeanNDVI = round(MeanNDVI, 3)
      )
    
    DT::datatable(
      out,
      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$bioScoreBox <- renderUI({
    df <- base_socio_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    fmt_count <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 0 || !is.finite(x)) return("N/A")
      scales::comma(round(x))
    }
    
    # One independently benchmarked result per generated isochrone. Never show a
    # percentile for the union/mean of multiple mode × time polygons.
    cards <- lapply(seq_len(nrow(df)), function(i) {
      total_raw <- if ("GBIF_Species" %in% names(df)) df$GBIF_Species[[i]] else NA_real_
      bird_raw <- if ("Bird_Species" %in% names(df)) df$Bird_Species[[i]] else NA_real_
      mammal_raw <- if ("Mammal_Species" %in% names(df)) df$Mammal_Species[[i]] else NA_real_
      plant_raw <- if ("Plant_Species" %in% names(df)) df$Plant_Species[[i]] else NA_real_
      
      total_pct <- if ("pctile_GBIF_Species" %in% names(df)) df$pctile_GBIF_Species[[i]] else NA_real_
      bird_pct <- if ("pctile_Bird_Species" %in% names(df)) df$pctile_Bird_Species[[i]] else NA_real_
      mammal_pct <- if ("pctile_Mammal_Species" %in% names(df)) df$pctile_Mammal_Species[[i]] else NA_real_
      plant_pct <- if ("pctile_Plant_Species" %in% names(df)) df$pctile_Plant_Species[[i]] else NA_real_
      n_ref <- if ("nref_GBIF_Species" %in% names(df)) df$nref_GBIF_Species[[i]] else NA_integer_
      
      tags$div(
        style = "padding:10px 0 12px 0; border-bottom:1px solid #e6e6e6;",
        tags$b(iso_metric_row_label(df, i)),
        tags$div(
          style = "display:flex; gap:12px; align-items:flex-end; flex-wrap:wrap; margin-top:5px;",
          tags$div(
            tags$div(
              style = "font-size:25px; line-height:1.05; font-weight:700; color:#28543b;",
              paste0(fmt_count(total_raw), " species")
            ),
            tags$small("Total unique GBIF species in this isochrone")
          ),
          tags$div(
            style = "padding-left:12px; border-left:3px solid #b9dfc8;",
            tags$div(
              style = "font-size:20px; line-height:1.05; font-weight:700;",
              format_pct(total_pct)
            ),
            tags$small("Citywide rank for the same mode × time")
          )
        ),
        tags$div(
          style = "display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:6px; margin-top:9px;",
          tags$div(
            style = "background:#f6faf7; border:1px solid #e0ebe4; border-radius:6px; padding:6px 8px;",
            tags$b("Birds"), tags$br(),
            tags$small(paste0(fmt_count(bird_raw), " species · ", format_pct(bird_pct)))
          ),
          tags$div(
            style = "background:#f6faf7; border:1px solid #e0ebe4; border-radius:6px; padding:6px 8px;",
            tags$b("Mammals"), tags$br(),
            tags$small(paste0(fmt_count(mammal_raw), " species · ", format_pct(mammal_pct)))
          ),
          tags$div(
            style = "background:#f6faf7; border:1px solid #e0ebe4; border-radius:6px; padding:6px 8px;",
            tags$b("Plants"), tags$br(),
            tags$small(paste0(fmt_count(plant_raw), " species · ", format_pct(plant_pct)))
          )
        ),
        if (is.finite(total_pct)) {
          tags$div(
            style = paste0(
              "background:#f0f7f3; border-left:4px solid #5d8d70; ",
              "border-radius:6px; padding:7px 9px; margin-top:8px; font-size:12px;"
            ),
            tags$b("Reporting sentence: "),
            paste0(
              iso_metric_row_label(df, i),
              " ranks at the ",
              format_pct(total_pct),
              " for accessible total species richness among comparable SF origins with the same mode and travel time."
            )
          )
        },
        if (is.finite(n_ref)) {
          tags$small(
            style = "display:block; color:#66766e; margin-top:6px;",
            paste0("Percentile reference uses ", scales::comma(n_ref), " finite comparable SF origins.")
          )
        }
      )
    })
    
    wellPanel(
      tags$div(
        style = "max-height:390px; overflow-y:auto;",
        cards
      ),
      tags$small(
        style = "display:block; color:#666; margin-top:7px;",
        "Raw species counts and citywide percentiles are shown separately. Each result is compared only with precomputed SF origins having the same transportation mode and travel time. Isochrones are never pooled."
      )
    )
  })
  
  output$transitScoreBox <- renderUI({
    df <- base_socio_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    cards <- lapply(seq_len(nrow(df)), function(i) {
      raw <- if ("Transit_Access_Score" %in% names(df)) {
        suppressWarnings(as.numeric(df$Transit_Access_Score[[i]]))
      } else NA_real_
      pct <- if ("pctile_Transit_Access_Score" %in% names(df)) {
        suppressWarnings(as.numeric(df$pctile_Transit_Access_Score[[i]]))
      } else NA_real_
      n_ref <- if ("nref_Transit_Access_Score" %in% names(df)) df$nref_Transit_Access_Score[[i]] else NA_integer_
      
      tags$div(
        style = "padding:8px 0; border-bottom:1px solid #e6e6e6;",
        tags$b(iso_metric_row_label(df, i)),
        tags$div(
          style = "font-size:21px; line-height:1.15; margin-top:2px;",
          if (is.finite(raw)) paste0(round(raw, 2), " stops/km²") else "N/A"
        ),
        tags$small(
          "Citywide isochrone percentile: ", format_pct(pct),
          if (is.finite(n_ref)) paste0(" · n=", scales::comma(n_ref)) else ""
        )
      )
    })
    
    wellPanel(
      tags$div(style = "max-height:310px; overflow-y:auto;", cards),
      tags$small(
        style = "display:block; color:#666; margin-top:7px;",
        "Raw transit values are retained. A percentile is shown only when a comparable same-mode × same-time precomputed reference exists."
      )
    )
  })
  
  output$dataCoverageBox <- renderUI({
    df <- baseline_biodiversity_access_index()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else x[[1]]
    }
    
    cards <- lapply(seq_len(nrow(df)), function(i) {
      records <- if ("GBIF_Records" %in% names(df)) safe_num(df$GBIF_Records[[i]]) else NA_real_
      density <- if ("SamplingDensity_km2" %in% names(df)) safe_num(df$SamplingDensity_km2[[i]]) else NA_real_
      pct <- if ("pctile_SamplingDensity_km2" %in% names(df)) safe_num(df$pctile_SamplingDensity_km2[[i]]) else NA_real_
      n_ref <- if ("nref_SamplingDensity_km2" %in% names(df)) safe_num(df$nref_SamplingDensity_km2[[i]]) else NA_real_
      
      interpretation <- if (!is.finite(pct)) {
        "No comparable observation-density benchmark is available."
      } else if (pct >= 75) {
        "Relatively well sampled compared with matched SF isochrones."
      } else if (pct >= 25) {
        "Mid-range sampling intensity compared with matched SF isochrones."
      } else {
        "Relatively lightly sampled compared with matched SF isochrones."
      }
      
      tags$div(
        style = paste0(
          "border-left:5px solid ", BAI_GROUP_COLORS[["biodiversity"]],
          "; background:#fffaf3; border-radius:8px; padding:9px 11px; ",
          "margin-bottom:8px;"
        ),
        tags$b(iso_metric_row_label(df, i)),
        tags$div(
          style = "display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:8px; margin-top:7px;",
          tags$div(
            tags$small(tags$b("GBIF records")),
            tags$div(
              style = "font-size:20px; font-weight:700;",
              if (is.finite(records)) scales::comma(round(records)) else "N/A"
            )
          ),
          tags$div(
            tags$small(tags$b("Sampling density")),
            tags$div(
              style = "font-size:20px; font-weight:700;",
              if (is.finite(density)) paste0(round(density, 1), " records/km²") else "N/A"
            )
          ),
          tags$div(
            tags$small(tags$b("Citywide sampling-density rank")),
            tags$div(
              style = "font-size:20px; font-weight:700;",
              format_pct(pct)
            ),
            if (is.finite(n_ref)) tags$small(paste0("n=", scales::comma(round(n_ref)))) else NULL
          )
        ),
        tags$div(
          style = "margin-top:7px; font-size:12px; color:#5d5245;",
          tags$b("Interpretation: "),
          interpretation
        )
      )
    })
    
    tagList(
      cards,
      tags$small(
        style = "display:block; color:#66766e; margin-top:5px;",
        "This is an observation-coverage indicator, not statistical confidence and not proof that true species richness is completely sampled. The descriptive low/mid/high wording is based on the matched sampling-density percentile."
      )
    )
  })
  
  output$baiDimensionDetailsBox <- renderUI({
    df <- baseline_biodiversity_access_index()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else x[[1]]
    }
    
    spec <- list(
      list(group = "ACCESS", color = BAI_GROUP_COLORS[["access"]],
           label = "Transit access", raw = "Transit_Access_Score",
           pct = "pctile_Transit_Access_Score", nref = "nref_Transit_Access_Score",
           fmt = function(x) if (is.finite(x)) paste0(round(x, 2), " stops/km²") else "N/A"),
      list(group = "ACCESS", color = BAI_GROUP_COLORS[["access"]],
           label = "Muni route access", raw = "Unique_Muni_Routes",
           pct = "pctile_Unique_Muni_Routes", nref = "nref_Unique_Muni_Routes",
           fmt = function(x) if (is.finite(x)) paste0(round(x), " routes") else "N/A"),
      list(group = "BIODIVERSITY EVIDENCE", color = BAI_GROUP_COLORS[["biodiversity"]],
           label = "Total species richness", raw = "GBIF_Species",
           pct = "pctile_GBIF_Species", nref = "nref_GBIF_Species",
           fmt = function(x) if (is.finite(x)) paste0(scales::comma(round(x)), " species") else "N/A"),
      list(group = "BIODIVERSITY EVIDENCE", color = BAI_GROUP_COLORS[["biodiversity"]],
           label = "GBIF sampling density", raw = "SamplingDensity_km2",
           pct = "pctile_SamplingDensity_km2", nref = "nref_SamplingDensity_km2",
           fmt = function(x) if (is.finite(x)) paste0(round(x, 1), " records/km²") else "N/A"),
      list(group = "GREEN ENVIRONMENT", color = BAI_GROUP_COLORS[["green_environment"]],
           label = "Vegetation NDVI", raw = "MeanNDVI",
           pct = "pctile_MeanNDVI", nref = "nref_MeanNDVI",
           fmt = function(x) if (is.finite(x)) format(round(x, 3), nsmall = 3) else "N/A"),
      list(group = "GREEN ENVIRONMENT", color = BAI_GROUP_COLORS[["green_environment"]],
           label = "Greenspace cover", raw = "Greenspace_percent",
           pct = "pctile_Greenspace_percent", nref = "nref_Greenspace_percent",
           fmt = function(x) if (is.finite(x)) paste0(round(x, 1), "%") else "N/A"),
      list(group = "EQUITY CONTEXT", color = BAI_GROUP_COLORS[["equity"]],
           label = "SF EJ context", raw = "SF_EJ_Score",
           pct = "pctile_SF_EJ_Score", nref = "nref_SF_EJ_Score",
           fmt = function(x) if (is.finite(x)) paste0("raw score ", round(x, 1)) else "N/A")
    )
    
    iso_sections <- lapply(seq_len(nrow(df)), function(i) {
      metric_cards <- lapply(spec, function(s) {
        raw <- if (s$raw %in% names(df)) safe_num(df[[s$raw]][[i]]) else NA_real_
        pct <- if (s$pct %in% names(df)) safe_num(df[[s$pct]][[i]]) else NA_real_
        nref <- if (s$nref %in% names(df)) safe_num(df[[s$nref]][[i]]) else NA_real_
        
        tags$div(
          style = paste0(
            "border-top:4px solid ", s$color,
            "; background:white; border-radius:8px; padding:8px 9px; ",
            "border-left:1px solid #e0e6e2; border-right:1px solid #e0e6e2; ",
            "border-bottom:1px solid #e0e6e2;"
          ),
          tags$div(
            style = paste0("font-size:11px; font-weight:800; color:", s$color, "; letter-spacing:.03em;"),
            s$group
          ),
          tags$b(s$label),
          tags$div(
            style = "font-size:18px; font-weight:700; margin-top:3px;",
            s$fmt(raw)
          ),
          tags$small(
            "Citywide rank: ", format_pct(pct),
            if (is.finite(nref)) paste0(" · n=", scales::comma(round(nref))) else ""
          ),
          if (identical(s$raw, "SF_EJ_Score")) {
            tags$small(
              style = "display:block; color:#6b5a78; margin-top:2px;",
              "Percentile is inverted once: higher rank = lower EJ burden."
            )
          } else NULL
        )
      })
      
      tags$div(
        style = "margin-bottom:14px;",
        tags$div(
          style = "font-size:16px; font-weight:700; margin-bottom:7px; color:#294337;",
          iso_metric_row_label(df, i)
        ),
        tags$div(
          style = "display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:8px;",
          metric_cards
        )
      )
    })
    
    tagList(
      iso_sections,
      tags$small(
        style = "display:block; color:#66766e;",
        "Every raw value is paired with its independently calculated citywide percentile from the same transportation mode × travel-time reference."
      )
    )
  })
  
  output$biodiversityAccessIndexBox <- renderUI({
    df <- baseline_biodiversity_access_index()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    cards <- lapply(seq_len(nrow(df)), function(i) {
      bai <- if ("BAI" %in% names(df)) suppressWarnings(as.numeric(df$BAI[[i]])) else NA_real_
      bai_pct <- if ("pctile_BAI" %in% names(df)) suppressWarnings(as.numeric(df$pctile_BAI[[i]])) else NA_real_
      n_bai_ref <- if ("nref_BAI" %in% names(df)) suppressWarnings(as.numeric(df$nref_BAI[[i]])) else NA_real_
      n_comp <- if ("BAI_n_components" %in% names(df)) df$BAI_n_components[[i]] else NA_integer_
      
      tags$div(
        style = "padding:9px 0; border-bottom:1px solid #e6e6e6;",
        tags$b(iso_metric_row_label(df, i)),
        tags$div(
          style = "display:grid; grid-template-columns:1fr 1fr; gap:9px; margin-top:5px;",
          tags$div(
            style = "background:#fff9ee; border-radius:7px; padding:7px 9px;",
            tags$small(tags$b("BAI score")),
            tags$div(
              style = "font-size:24px; line-height:1.15; font-weight:700;",
              if (is.finite(bai)) paste0(round(100 * bai, 1), "/100") else "N/A"
            ),
            tags$small("Mean of the seven percentile-standardized dimensions")
          ),
          tags$div(
            style = "background:#f3f7fb; border-radius:7px; padding:7px 9px;",
            tags$small(tags$b("BAI citywide rank")),
            tags$div(
              style = "font-size:21px; line-height:1.15; font-weight:700;",
              format_pct(bai_pct)
            ),
            tags$small(
              if (is.finite(n_bai_ref)) {
                paste0("Compared with ", scales::comma(round(n_bai_ref)), " complete seven-axis reference isochrones")
              } else {
                "No complete matching BAI reference distribution"
              }
            )
          )
        ),
        tags$small(
          style = "display:block; margin-top:6px;",
          if (is.finite(n_comp)) {
            paste0(n_comp, " of 7 same-mode × same-time percentile dimensions available")
          } else {
            "No comparable precomputed reference dimensions"
          }
        )
      )
    })
    
    wellPanel(
      tags$div(style = "max-height:380px; overflow-y:auto;", cards),
      tags$small(
        style = "display:block; color:#666; margin-top:7px;",
        "BAI score and BAI percentile are different quantities. The score is the seven-axis mean; the citywide rank compares that score with complete seven-axis BAIs from the exact same mode × travel-time reference. Each user isochrone is ranked independently."
      )
    )
  })
  
  output$greenspaceScenarioSummary <- renderUI({
    proposal <- proposed_greenspace()
    if (is.null(proposal) || nrow(proposal) == 0) {
      return(tags$small(
        style = "color:#666;",
        "No proposed habitat or corridor drawn."
      ))
    }
    
    safe_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else x[[1]]
    }
    
    gross_m2 <- tryCatch({
      footprint <- proposed_habitat_3857()
      if (is.null(footprint) || nrow(footprint) == 0) NA_real_ else
        sum(as.numeric(st_area(footprint)), na.rm = TRUE)
    }, error = function(e) NA_real_)
    
    proposal_net <- net_proposed_greenspace_3857()
    net_m2 <- if (is.null(proposal_net) || nrow(proposal_net) == 0) {
      0
    } else {
      sum(as.numeric(st_area(proposal_net)), na.rm = TRUE)
    }
    
    selected_metrics <- corridor_selected_isochrone_metrics()
    cover_line <- NULL
    bai_line <- NULL
    
    if (!is.null(selected_metrics)) {
      base_cover <- safe_num(selected_metrics$baseline_raw$Greenspace_percent)
      scenario_cover <- safe_num(selected_metrics$scenario_raw$Greenspace_percent)
      base_cover_pct <- safe_num(selected_metrics$baseline_bai$pctile_Greenspace_percent)
      scenario_cover_pct <- safe_num(selected_metrics$scenario_bai$pctile_Greenspace_percent)
      
      if (is.finite(base_cover) && is.finite(scenario_cover)) {
        cover_line <- tags$p(
          style = "margin:4px 0 0 0;",
          tags$b(paste0(selected_metrics$label, " greenspace cover: ")),
          paste0(
            round(base_cover, 1), "% (", format_pct(base_cover_pct), ") → ",
            round(scenario_cover, 1), "% (", format_pct(scenario_cover_pct), ")"
          )
        )
      }
      
      base_bai_value <- 100 * safe_num(selected_metrics$baseline_bai$BAI)
      scenario_bai_value <- 100 * safe_num(selected_metrics$scenario_bai$BAI)
      base_bai_pct <- safe_num(selected_metrics$baseline_bai$pctile_BAI)
      scenario_bai_pct <- safe_num(selected_metrics$scenario_bai$pctile_BAI)
      
      if (is.finite(base_bai_value) && is.finite(scenario_bai_value)) {
        bai_line <- tags$p(
          style = "margin:2px 0 0 0;",
          tags$b("Selected-isochrone BAI: "),
          paste0(
            round(base_bai_value, 1), "/100 (", format_pct(base_bai_pct), ") → ",
            round(scenario_bai_value, 1), "/100 (", format_pct(scenario_bai_pct), ")"
          )
        )
      }
    }
    
    tagList(
      tags$p(
        style = "margin:6px 0 0 0;",
        tags$b("Your proposed habitat area: "),
        if (is.finite(gross_m2)) paste0(round(gross_m2 / 10000, 2), " ha") else "N/A"
      ),
      tags$p(
        style = "margin:2px 0 0 0;",
        tags$b("New area not already mapped as greenspace: "),
        paste0(round(net_m2 / 10000, 2), " ha")
      ),
      tags$small(
        style = "display:block; color:#555; margin:4px 0 6px 0;",
        "The blue area is what the app treats as the intervention. A drawn line is widened using the selected corridor planting width; a drawn polygon is used directly."
      ),
      cover_line,
      bai_line,
      tags$small(
        style = "color:#666;",
        "Scenario changes greenspace cover and BAI only; it does not predict future NDVI or biodiversity."
      )
    )
  })
  
  output$corridorConnectivitySummary <- renderUI({
    c <- corridor_connectivity()
    if (is.null(c) || !isTRUE(c$available)) {
      return(tags$small(
        style = "color:#666;",
        if (!is.null(c$message)) c$message else "Connectivity estimate unavailable."
      ))
    }
    
    fmt_num <- function(x, digits = 0) {
      if (is.null(x) || !is.finite(x)) "N/A" else scales::comma(round(x, digits))
    }
    fmt_ha <- function(x) {
      if (is.null(x) || !is.finite(x)) "N/A" else paste0(round(x, 2), " ha")
    }
    
    interpretation <- dplyr::case_when(
      c$classification == "Bridge" ~ paste0(
        "This proposal brings ", c$linked_patches,
        " existing mapped greenspaces from ", c$baseline_components,
        " separate baseline networks into one potential corridor network."
      ),
      c$classification == "Network reinforcement" ~ paste0(
        "This proposal adds habitat around ", c$linked_patches,
        " existing patches that are already structurally connected under the selected threshold."
      ),
      c$classification == "Patch expansion" ~
        "This proposal directly overlaps and enlarges one existing mapped greenspace.",
      c$classification == "Stepping-stone" ~
        "This proposal creates a potential stepping stone near an existing mapped greenspace, but it does not join separate baseline networks.",
      TRUE ~ paste0(
        "No existing mapped greenspace lies within the selected ", c$gap_m,
        " m connection threshold. The intervention is currently classified as isolated."
      )
    )
    
    tagList(
      tags$div(
        style = "background:#e3f2fd; border-left:5px solid #1976d2; padding:12px 14px; margin-bottom:12px;",
        tags$h3(style = "margin:0 0 4px 0;", c$classification),
        tags$p(style = "margin:0; font-size:15px;", interpretation)
      ),
      fluidRow(
        column(
          6,
          wellPanel(
            tags$h3(style = "margin-top:0;", fmt_num(c$linked_patches)),
            tags$b("Existing greenspaces reached"),
            tags$br(),
            tags$small("Orange and numbered P1, P2, … on the map")
          )
        ),
        column(
          6,
          wellPanel(
            tags$h3(style = "margin-top:0;", paste0(c$baseline_components, " → ", c$scenario_components)),
            tags$b("Habitat networks: baseline → scenario"),
            tags$br(),
            tags$small(paste0("Using a ", c$gap_m, " m connection rule"))
          )
        ),
        column(
          6,
          wellPanel(
            tags$h3(style = "margin-top:0;", fmt_ha(c$scenario_connected_area_ha)),
            tags$b("Potential connected habitat area"),
            tags$br(),
            tags$small(paste0(fmt_ha(c$existing_linked_area_ha), " existing + ", fmt_ha(c$novel_area_ha), " novel"))
          )
        ),
        column(
          6,
          wellPanel(
            tags$h3(
              style = "margin-top:0;",
              if (is.finite(c$maximum_linked_gap_m)) paste0(fmt_num(c$maximum_linked_gap_m), " m") else "N/A"
            ),
            tags$b("Largest proposal-to-patch gap"),
            tags$br(),
            tags$small("Dashed orange lines show these measured gaps")
          )
        )
      ),
      if (is.finite(c$corridor_length_m) && c$corridor_length_m > 0) {
        tags$p(
          style = "margin:8px 0 0 0; color:#555;",
          paste0(
            "Drawn corridor length: ", fmt_num(c$corridor_length_m),
            " m; planted width: ", input$corridor_width_m, " m."
          )
        )
      },
      tags$small(
        style = "color:#666; display:block; margin-top:6px;",
        paste0(
          "Coverage note: the analysis searches all ",
          scales::comma(osm_greenspace_feature_count),
          " mapped OSM greenspace polygons loaded for San Francisco, but OSM may not represent every real or ecologically suitable habitat. This result describes structural proximity only; it does not model vegetation quality, roads, mortality risk, or species-specific movement."
        )
      )
    )
  })
  
  output$linkedHabitatTable <- renderDT({
    c <- corridor_connectivity()
    if (is.null(c) || !isTRUE(c$available) ||
        is.null(c$linked_sf) || nrow(c$linked_sf) == 0) {
      return(DT::datatable(
        data.frame(Message = "No existing mapped greenspace is affected under the selected connection rule."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    linked_table <- c$linked_sf |>
      st_drop_geometry() |>
      transmute(
        ID = display_id,
        `Existing greenspace` = display_name,
        Relationship = relationship,
        `Distance to proposal (m)` = round(distance_to_proposal_m),
        `Existing area (ha)` = round(patch_area_ha, 2),
        `Baseline network` = baseline_component
      )
    
    DT::datatable(
      linked_table,
      rownames = FALSE,
      escape = TRUE,
      options = list(
        pageLength = 7,
        lengthChange = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })
  
  
  # ---------------------------------------------------------------------------
  # Corridor scenario exports
  # ---------------------------------------------------------------------------
  corridor_scenario_name <- reactive({
    raw_name <- input$corridor_scenario_name
    nm <- if (is.null(raw_name)) "" else trimws(as.character(raw_name))
    if (!nzchar(nm)) "Green corridor scenario" else nm
  })
  
  safe_export_stem <- function(x) {
    x <- iconv(x, to = "ASCII//TRANSLIT", sub = "")
    x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) "green_corridor_scenario" else substr(x, 1, 80)
  }
  
  # Build a mixed-geometry GeoJSON package containing the analyzed intervention,
  # novel habitat, reached existing greenspaces, measured gap lines, the optional
  # population-estimate area, and the transferred Explorer location when present.
  scenario_export_sf <- reactive({
    footprint <- proposed_habitat_3857()
    if (is.null(footprint) || nrow(footprint) == 0) return(NULL)
    
    scenario_name <- corridor_scenario_name()
    generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    gap_rule <- suppressWarnings(as.numeric(input$connectivity_gap_m))
    width_m <- suppressWarnings(as.numeric(input$corridor_width_m))
    connectivity <- corridor_connectivity()
    beneficiaries <- scenario_beneficiaries()
    active_iso_info <- corridor_active_iso_info()
    active_iso_label <- if (is.null(active_iso_info)) NA_character_ else active_iso_info$label
    connectivity_class <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$classification else NA_character_
    baseline_network_count <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$baseline_components else NA_integer_
    scenario_network_count <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$scenario_components else NA_integer_
    linked_patch_count <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$linked_patches else NA_integer_
    potential_population <- if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) beneficiaries$total_population else NA_real_
    ej_population <- if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) beneficiaries$ej_population else NA_real_
    
    make_layer <- function(x, feature_type, feature_id = NULL,
                           feature_name = NULL, relationship = NULL,
                           distance_m = NULL, baseline_network = NULL,
                           walk_minutes = NULL) {
      if (is.null(x) || !inherits(x, "sf") || nrow(x) == 0) return(NULL)
      x <- tryCatch(
        st_zm(x, drop = TRUE, what = "ZM"),
        error = function(e) x
      )
      x <- suppressWarnings(st_transform(x, 4326))
      n <- nrow(x)
      recycle <- function(value, default = NA) {
        if (is.null(value)) return(rep(default, n))
        rep(value, length.out = n)
      }
      st_sf(
        scenario_name = rep(scenario_name, n),
        feature_type = rep(feature_type, n),
        feature_id = as.character(recycle(feature_id, NA_character_)),
        feature_name = as.character(recycle(feature_name, NA_character_)),
        relationship = as.character(recycle(relationship, NA_character_)),
        distance_m = suppressWarnings(as.numeric(recycle(distance_m, NA_real_))),
        baseline_network = suppressWarnings(as.integer(recycle(baseline_network, NA_integer_))),
        connection_rule_m = rep(gap_rule, n),
        selected_explorer_isochrone = rep(active_iso_label, n),
        assumed_planted_width_m = rep(width_m, n),
        population_walk_minutes = suppressWarnings(as.numeric(recycle(walk_minutes, NA_real_))),
        connectivity_classification = rep(connectivity_class, n),
        baseline_network_count = rep(baseline_network_count, n),
        scenario_network_count = rep(scenario_network_count, n),
        linked_patch_count = rep(linked_patch_count, n),
        potentially_reached_population = rep(potential_population, n),
        ej_population = rep(ej_population, n),
        generated_at = rep(generated_at, n),
        geometry = st_geometry(x),
        crs = 4326
      )
    }
    
    layers <- list(
      make_layer(
        footprint,
        feature_type = "proposed_habitat_area",
        feature_id = "proposal",
        feature_name = scenario_name,
        relationship = "Analyzed footprint"
      )
    )
    
    novel <- net_proposed_greenspace_3857()
    if (!is.null(novel) && nrow(novel) > 0) {
      layers[[length(layers) + 1]] <- make_layer(
        novel,
        feature_type = "novel_habitat_area",
        feature_id = "novel_area",
        feature_name = paste0(scenario_name, " — area not already mapped as greenspace"),
        relationship = "Excludes overlap with existing mapped OSM greenspace"
      )
    }
    
    if (!is.null(connectivity) && isTRUE(connectivity$available)) {
      linked <- connectivity$linked_sf
      if (!is.null(linked) && nrow(linked) > 0) {
        layers[[length(layers) + 1]] <- make_layer(
          linked,
          feature_type = "existing_greenspace_reached",
          feature_id = linked$display_id,
          feature_name = linked$display_name,
          relationship = linked$relationship,
          distance_m = linked$distance_to_proposal_m,
          baseline_network = linked$baseline_component
        )
      }
      
      links <- connectivity$links_sf
      if (!is.null(links) && nrow(links) > 0) {
        layers[[length(layers) + 1]] <- make_layer(
          links,
          feature_type = "measured_connection_gap",
          feature_id = links$display_id,
          feature_name = links$patch_name,
          relationship = ifelse(
            links$distance_m <= 0.01,
            "Direct overlap",
            paste0(round(links$distance_m), " m measured gap")
          ),
          distance_m = links$distance_m
        )
      }
    }
    
    if (isTRUE(input$show_beneficiary_catchment)) {
      catchment <- scenario_beneficiary_catchment_3857()
      if (!is.null(catchment) && nrow(catchment) > 0) {
        layers[[length(layers) + 1]] <- make_layer(
          catchment,
          feature_type = "population_estimate_area",
          feature_id = "population_area",
          feature_name = paste0(catchment$walk_minutes[[1]], "-minute approximate population-estimate area"),
          relationship = "Straight-line proximity screen; not a routed walkshed",
          walk_minutes = catchment$walk_minutes
        )
      }
    }
    
    ref_pt <- corridor_reference_point()
    if (!is.null(ref_pt) && length(ref_pt) >= 2) {
      ref_sf <- st_as_sf(
        data.frame(
          lon = unname(ref_pt[["lon"]]),
          lat = unname(ref_pt[["lat"]])
        ),
        coords = c("lon", "lat"),
        crs = 4326
      )
      layers[[length(layers) + 1]] <- make_layer(
        ref_sf,
        feature_type = "explorer_reference_location",
        feature_id = "explorer_location",
        feature_name = "Location transferred from Isochrone Explorer",
        relationship = "Context point only"
      )
    }
    
    layers <- Filter(function(x) !is.null(x) && nrow(x) > 0, layers)
    if (length(layers) == 0) return(NULL)
    do.call(rbind, layers)
  })
  
  
  html_value <- function(x) {
    htmltools::htmlEscape(as.character(x), attribute = FALSE)
  }
  
  report_metric <- function(label, value, detail = NULL) {
    paste0(
      '<div class="metric"><div class="metric-label">', html_value(label), '</div>',
      '<div class="metric-value">', html_value(value), '</div>',
      if (!is.null(detail) && nzchar(detail)) paste0('<div class="metric-detail">', html_value(detail), '</div>') else "",
      '</div>'
    )
  }
  
  # Draw a lightweight static spatial summary for the downloadable report. This
  # avoids external map tiles and keeps the report self-contained and printable.
  write_corridor_report_map <- function(file) {
    footprint <- proposed_habitat_3857()
    connectivity <- corridor_connectivity()
    if (is.null(footprint) || nrow(footprint) == 0) return(FALSE)
    
    linked <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$linked_sf else NULL
    links <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$links_sf else NULL
    catchment <- if (isTRUE(input$show_beneficiary_catchment)) scenario_beneficiary_catchment_3857() else NULL
    
    focus_parts <- list(st_geometry(footprint))
    if (!is.null(linked) && nrow(linked) > 0) focus_parts[[length(focus_parts) + 1]] <- st_geometry(linked)
    focus_geom <- do.call(c, focus_parts)
    bb <- st_bbox(focus_geom)
    pad_x <- max(as.numeric(bb[["xmax"]] - bb[["xmin"]]) * 0.18, 250)
    pad_y <- max(as.numeric(bb[["ymax"]] - bb[["ymin"]]) * 0.18, 250)
    map_bb <- st_bbox(c(
      xmin = bb[["xmin"]] - pad_x,
      ymin = bb[["ymin"]] - pad_y,
      xmax = bb[["xmax"]] + pad_x,
      ymax = bb[["ymax"]] + pad_y
    ), crs = st_crs(3857))
    
    local_green <- tryCatch(
      suppressWarnings(st_crop(osm_greenspace_3857, map_bb)),
      error = function(e) osm_greenspace_3857[0, , drop = FALSE]
    )
    local_cbg <- tryCatch(
      suppressWarnings(st_crop(st_transform(cbg_vect_sf, 3857), map_bb)),
      error = function(e) NULL
    )
    
    grDevices::png(file, width = 1500, height = 920, res = 170, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(0.3, 0.3, 1.1, 0.3), xaxs = "i", yaxs = "i")
    
    graphics::plot(
      st_geometry(footprint),
      col = NA, border = NA,
      xlim = c(map_bb[["xmin"]], map_bb[["xmax"]]),
      ylim = c(map_bb[["ymin"]], map_bb[["ymax"]]),
      axes = FALSE, asp = 1
    )
    if (!is.null(local_cbg) && nrow(local_cbg) > 0) {
      graphics::plot(st_geometry(local_cbg), add = TRUE, col = "#f5f5f5", border = "#d8d8d8", lwd = 0.5)
    }
    if (!is.null(catchment) && nrow(catchment) > 0) {
      graphics::plot(st_geometry(catchment), add = TRUE, col = grDevices::adjustcolor("#ce93d8", alpha.f = 0.20), border = "#6a1b9a", lty = 2, lwd = 2)
    }
    if (!is.null(local_green) && nrow(local_green) > 0) {
      graphics::plot(st_geometry(local_green), add = TRUE, col = "#b7d9b0", border = "#2e7d32", lwd = 0.8)
    }
    if (!is.null(linked) && nrow(linked) > 0) {
      graphics::plot(st_geometry(linked), add = TRUE, col = "#ffb74d", border = "#e65100", lwd = 2.2)
    }
    graphics::plot(st_geometry(footprint), add = TRUE, col = grDevices::adjustcolor("#42a5f5", alpha.f = 0.55), border = "#0d47a1", lwd = 2.4)
    if (!is.null(links) && nrow(links) > 0) {
      graphics::plot(st_geometry(links), add = TRUE, col = "#ef6c00", lty = 2, lwd = 2)
    }
    if (!is.null(linked) && nrow(linked) > 0) {
      pts <- suppressWarnings(st_point_on_surface(linked))
      xy <- st_coordinates(pts)
      graphics::text(xy[, 1], xy[, 2], labels = linked$display_id, cex = 0.72, font = 2, col = "#7f2704")
    }
    graphics::legend(
      "bottomleft",
      legend = c("Existing greenspace", "Proposal", "Existing greenspaces reached", "Measured gap"),
      fill = c("#b7d9b0", "#42a5f5", "#ffb74d", NA),
      border = c("#2e7d32", "#0d47a1", "#e65100", NA),
      lty = c(NA, NA, NA, 2),
      col = c(NA, NA, NA, "#ef6c00"),
      lwd = c(NA, NA, NA, 2),
      bg = grDevices::adjustcolor("white", alpha.f = 0.90),
      cex = 0.72,
      bty = "o"
    )
    graphics::mtext("Screening map — schematic context without street basemap", side = 3, line = 0.15, cex = 0.72, col = "#555555")
    TRUE
  }
  
  output$download_corridor_report <- downloadHandler(
    filename = function() {
      paste0(safe_export_stem(corridor_scenario_name()), "_screening_report.html")
    },
    contentType = "text/html",
    content = function(file) {
      footprint <- proposed_habitat_3857()
      if (is.null(footprint) || nrow(footprint) == 0) {
        stop("Draw a corridor line or habitat polygon before downloading a report.")
      }
      
      connectivity <- corridor_connectivity()
      beneficiaries <- scenario_beneficiaries()
      gross_ha <- sum(as.numeric(st_area(footprint)), na.rm = TRUE) / 10000
      novel <- net_proposed_greenspace_3857()
      novel_ha <- if (is.null(novel) || nrow(novel) == 0) 0 else sum(as.numeric(st_area(novel)), na.rm = TRUE) / 10000
      
      fmt_num <- function(x, digits = 0) if (is.null(x) || !is.finite(x)) "Not available" else scales::comma(round(x, digits))
      fmt_ha <- function(x) if (is.null(x) || !is.finite(x)) "Not available" else paste0(round(x, 2), " ha")
      fmt_pct <- function(x) if (is.null(x) || !is.finite(x)) "Not available" else paste0(round(x, 1), "%")
      fmt_money <- function(x) if (is.null(x) || !is.finite(x)) "Not available" else scales::dollar(round(x))
      
      class_value <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$classification else "Not available"
      linked_value <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$linked_patches else NA_real_
      network_value <- if (!is.null(connectivity) && isTRUE(connectivity$available)) paste0(connectivity$baseline_components, " → ", connectivity$scenario_components) else "Not available"
      connected_area <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$scenario_connected_area_ha else NA_real_
      largest_gap <- if (!is.null(connectivity) && isTRUE(connectivity$available)) connectivity$maximum_linked_gap_m else NA_real_
      gap_rule <- suppressWarnings(as.numeric(input$connectivity_gap_m))
      width_m <- suppressWarnings(as.numeric(input$corridor_width_m))
      
      selected_access <- corridor_selected_isochrone_metrics()
      active_iso_label <- if (is.null(selected_access)) "No Explorer isochrone selected" else selected_access$label
      bai_base <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$baseline_bai$BAI[[1]])) * 100 else NA_real_
      bai_scenario <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$scenario_bai$BAI[[1]])) * 100 else NA_real_
      bai_pct_base <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$baseline_bai$pctile_BAI[[1]])) else NA_real_
      bai_pct_scenario <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$scenario_bai$pctile_BAI[[1]])) else NA_real_
      cover_base <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$baseline_raw$Greenspace_percent[[1]])) else NA_real_
      cover_scenario <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$scenario_raw$Greenspace_percent[[1]])) else NA_real_
      cover_pct_base <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$baseline_bai$pctile_Greenspace_percent[[1]])) else NA_real_
      cover_pct_scenario <- if (!is.null(selected_access)) suppressWarnings(as.numeric(selected_access$scenario_bai$pctile_Greenspace_percent[[1]])) else NA_real_
      
      map_file <- tempfile(fileext = ".png")
      on.exit(unlink(map_file), add = TRUE)
      map_ok <- tryCatch(write_corridor_report_map(map_file), error = function(e) FALSE)
      map_uri <- if (isTRUE(map_ok) && file.exists(map_file) && requireNamespace("base64enc", quietly = TRUE)) {
        paste0("data:image/png;base64,", base64enc::base64encode(map_file))
      } else {
        NULL
      }
      
      corridor_detail <- if (!is.null(connectivity) && isTRUE(connectivity$available) &&
                             is.finite(connectivity$corridor_length_m) && connectivity$corridor_length_m > 0) {
        paste0("Route length: ", fmt_num(connectivity$corridor_length_m), " m; assumed planted width: ", fmt_num(width_m), " m")
      } else {
        paste0("Assumed line width: ", fmt_num(width_m), " m; polygons retain their drawn area")
      }
      access_value <- if (is.finite(cover_base) && is.finite(cover_scenario) &&
                          is.finite(bai_base) && is.finite(bai_scenario)) {
        paste0("Cover ", round(cover_base, 1), "% → ", round(cover_scenario, 1), "%")
      } else {
        "Generate isochrones to calculate"
      }
      access_detail <- if (is.finite(bai_base) && is.finite(bai_scenario)) {
        paste0(
          active_iso_label,
          "; greenspace percentile: ", round(cover_pct_base, 1), " → ", round(cover_pct_scenario, 1),
          "; BAI score: ", round(bai_base, 1), " → ", round(bai_scenario, 1), "/100",
          "; BAI citywide percentile: ", round(bai_pct_base, 1), " → ", round(bai_pct_scenario, 1),
          "; only greenspace cover changes"
        )
      } else {
        "BAI comparison requires a selected isochrone from the Explorer"
      }
      
      metric_html <- paste0(
        report_metric("Connectivity classification", class_value, paste0("Using a ", fmt_num(gap_rule), " m connection rule")),
        report_metric("Existing greenspaces reached", fmt_num(linked_value), "Mapped OSM greenspaces within the selected rule"),
        report_metric("Habitat networks", network_value, "Separate baseline groups → potential scenario group"),
        report_metric("Potential connected habitat", fmt_ha(connected_area), "Existing reached area plus novel proposed area"),
        report_metric("Proposed habitat area", fmt_ha(gross_ha), corridor_detail),
        report_metric("Novel habitat area", fmt_ha(novel_ha), "Excludes overlap with existing mapped greenspace"),
        report_metric("Largest proposal-to-patch gap", if (is.finite(largest_gap)) paste0(fmt_num(largest_gap), " m") else "Not available", "Measured straight-line separation"),
        report_metric(
          "Potentially reached residents",
          if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) fmt_num(beneficiaries$total_population) else "Not available",
          if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) paste0("Approximate ", beneficiaries$walk_minutes, "-minute area; lower-greenspace residents: ", fmt_num(beneficiaries$low_access_population)) else ""
        ),
        report_metric(
          "Residents in highest-burden SF EJ quartile",
          if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) fmt_num(beneficiaries$ej_population) else "Not available",
          if (!is.null(beneficiaries) && isTRUE(beneficiaries$available)) {
            paste0(
              fmt_pct(beneficiaries$ej_share),
              " of potentially reached residents",
              if (is.finite(beneficiaries$ej_vulnerability_threshold)) {
                paste0(
                  "; SF EJ score ≥ ",
                  round(beneficiaries$ej_vulnerability_threshold, 1),
                  " (75th-percentile cutoff)"
                )
              } else ""
            )
          } else ""
        ),
        report_metric("Selected isochrone access comparison", access_value, access_detail)
      )
      
      linked_rows <- ""
      linked_note <- "No existing mapped greenspaces were reached under the selected rule."
      if (!is.null(connectivity) && isTRUE(connectivity$available) &&
          !is.null(connectivity$linked_sf) && nrow(connectivity$linked_sf) > 0) {
        linked_df <- connectivity$linked_sf |>
          st_drop_geometry() |>
          transmute(
            ID = display_id,
            Name = display_name,
            Relationship = relationship,
            `Area (ha)` = round(patch_area_ha, 2),
            `Baseline network` = baseline_component
          ) |>
          slice_head(n = 10)
        linked_rows <- paste0(apply(linked_df, 1, function(row) {
          paste0(
            "<tr>",
            paste0("<td>", vapply(row, html_value, character(1)), "</td>", collapse = ""),
            "</tr>"
          )
        }), collapse = "")
        linked_note <- paste0(
          "Showing ", nrow(linked_df), " of ", nrow(connectivity$linked_sf),
          " reached greenspaces. The table below summarizes the mapped features included in this report card."
        )
      }
      
      map_html <- if (!is.null(map_uri)) {
        paste0('<img class="map-image" src="', map_uri, '" alt="Static map of the proposed corridor and existing greenspaces reached">')
      } else {
        '<div class="map-unavailable">Static report map unavailable.</div>'
      }
      
      scenario_name <- html_value(corridor_scenario_name())
      generated <- html_value(format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
      source_count <- html_value(scales::comma(osm_greenspace_feature_count))
      
      html <- paste0(
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width,initial-scale=1">',
        '<title>', scenario_name, ' — corridor screening report</title>',
        '<style>',
        '@page{size:letter landscape;margin:0.35in}*{box-sizing:border-box}',
        'body{margin:0;background:#eef3f1;color:#23352d;font-family:Arial,Helvetica,sans-serif;font-size:10px}',
        '.page{width:100%;max-width:11in;min-height:7.6in;margin:0 auto;background:white;padding:18px 22px}',
        'header{display:flex;justify-content:space-between;gap:20px;border-bottom:3px solid #2e7d32;padding-bottom:8px;margin-bottom:10px}',
        'h1{font-size:22px;margin:0;color:#205c35}h2{font-size:13px;margin:0 0 6px;color:#205c35}',
        '.subtitle{margin-top:4px;color:#56645d;font-size:11px}.meta{text-align:right;color:#56645d;white-space:nowrap}',
        '.main{display:grid;grid-template-columns:1.05fr .95fr;gap:12px}.panel{border:1px solid #cfd9d3;border-radius:7px;padding:9px;background:#fff}',
        '.map-image{display:block;width:100%;height:305px;object-fit:contain}.map-unavailable{height:305px;display:flex;align-items:center;justify-content:center;background:#f4f4f4;color:#666}',
        '.metrics{display:grid;grid-template-columns:1fr 1fr;gap:7px}.metric{border-left:4px solid #1976d2;background:#f3f8fb;padding:7px 8px;min-height:58px}',
        '.metric-label{font-size:9px;color:#4d5d55;text-transform:uppercase;letter-spacing:.03em}.metric-value{font-size:16px;font-weight:700;color:#173d2a;margin-top:2px}',
        '.metric-detail{font-size:8.5px;color:#66736c;margin-top:2px;line-height:1.25}',
        '.lower{display:grid;grid-template-columns:1.15fr .85fr;gap:12px;margin-top:10px}',
        'table{width:100%;border-collapse:collapse;font-size:8.5px}th,td{border-bottom:1px solid #dfe6e1;padding:3px 4px;text-align:left;vertical-align:top}th{background:#edf5ef;color:#2f5f3e}',
        '.note{background:#fff8e1;border-left:4px solid #f9a825;padding:8px 10px;line-height:1.35}.limitations{margin:4px 0 0 16px;padding:0}.limitations li{margin:2px 0}',
        '.screening{display:inline-block;background:#e3f2fd;color:#165b8b;border-radius:999px;padding:3px 8px;font-weight:700;margin-bottom:5px}',
        '.print-btn{position:fixed;right:14px;top:14px;background:#205c35;color:white;border:0;border-radius:5px;padding:8px 11px;cursor:pointer}',
        '.footer{margin-top:8px;padding-top:6px;border-top:1px solid #dfe6e1;color:#68756e;font-size:8.5px;display:flex;justify-content:space-between;gap:12px}',
        '@media print{body{background:white}.page{padding:0;max-width:none}.print-btn{display:none}}',
        '@media(max-width:850px){.main,.lower{grid-template-columns:1fr}.page{min-height:auto}.map-image{height:auto}}',
        '</style></head><body>',
        '<button class="print-btn" onclick="window.print()">Print / Save as PDF</button>',
        '<main class="page">',
        '<header><div><div class="screening">Preliminary screening report</div><h1>', scenario_name, '</h1>',
        '<div class="subtitle">Green Corridor Planner — habitat connectivity, greenspace access, and equity context</div></div>',
        '<div class="meta"><b>Generated</b><br>', generated, '<br><b>Connection rule</b><br>', html_value(fmt_num(gap_rule)), ' m</div></header>',
        '<section class="main"><div class="panel"><h2>Scenario map</h2>', map_html,
        '<div style="margin-top:5px;color:#647169">Green = existing mapped greenspace; blue = proposal; orange = existing greenspaces reached; dashed orange = measured gaps.</div></div>',
        '<div class="panel"><h2>Key screening results</h2><div class="metrics">', metric_html, '</div></div></section>',
        '<section class="lower"><div class="panel"><h2>Existing greenspaces reached</h2>',
        '<table><thead><tr><th>ID</th><th>Name</th><th>Relationship</th><th>Area (ha)</th><th>Baseline network</th></tr></thead><tbody>', linked_rows, '</tbody></table>',
        '<div style="margin-top:5px;color:#68756e">', html_value(linked_note), '</div></div>',
        '<div class="note"><h2>Interpretation and limitations</h2><ul class="limitations">',
        '<li>Connectivity is structural proximity under the selected distance rule; it is not evidence of actual wildlife movement.</li>',
        '<li>The analysis uses ', source_count, ' mapped OSM greenspace polygons loaded for San Francisco. OSM may omit or misclassify ecologically useful spaces.</li>',
        '<li>Population estimates use area-weighted census block groups and assume residents are uniformly distributed within each block group.</li>',
        '<li>The population-estimate area is a straight-line proximity buffer, not a routed pedestrian walkshed.</li>',
        '<li>Land ownership, street design, utilities, slopes, vegetation quality, road barriers, maintenance, cost, and community preferences are not assessed.</li>',
        '<li>Field assessment, agency review, ecological expertise, and community consultation are required before implementation.</li>',
        '</ul></div></section>',
        '<div class="footer"><span>Print or save this HTML report card as PDF if needed.</span><span>SF Biodiversity Access Decision Support Tool</span></div>',
        '</main></body></html>'
      )
      writeLines(html, con = file, useBytes = TRUE)
    }
  )
  
  output$greenspaceBeneficiarySummary <- renderUI({
    b <- scenario_beneficiaries()
    if (is.null(b) || !isTRUE(b$available)) {
      return(tags$small(style = "color:#666;", if (!is.null(b$message)) b$message else "Beneficiary estimates unavailable."))
    }
    
    format_people <- function(x) {
      if (!is.finite(x)) "N/A" else scales::comma(round(x))
    }
    format_pct <- function(x) {
      if (!is.finite(x)) "N/A" else paste0(round(x, 1), "%")
    }
    format_income <- function(x) {
      if (!is.finite(x)) "N/A" else scales::dollar(round(x))
    }
    
    fluidRow(
      column(
        3,
        wellPanel(
          tags$h3(format_people(b$total_population)),
          tags$b("Potentially served residents"),
          tags$br(),
          tags$small(paste0("Within the ", b$walk_minutes, "-minute catchment"))
        )
      ),
      column(
        3,
        wellPanel(
          tags$h3(format_people(b$ej_population)),
          tags$b("Residents in highest-burden SF EJ quartile"),
          tags$br(),
          tags$small(
            paste0(
              format_pct(b$ej_share),
              " of potentially served residents",
              if (is.finite(b$ej_vulnerability_threshold)) {
                paste0(
                  " · SF EJ score ≥ ",
                  round(b$ej_vulnerability_threshold, 1),
                  " (75th-percentile cutoff)"
                )
              } else ""
            )
          )
        )
      ),
      column(
        3,
        wellPanel(
          tags$h3(format_people(b$low_access_population)),
          tags$b("Residents in lower-greenspace CBGs"),
          tags$br(),
          tags$small(
            if (is.finite(b$low_access_threshold)) {
              paste0("Existing cover below city median: ", round(b$low_access_threshold, 1), "%")
            } else {
              "Existing greenspace benchmark unavailable"
            }
          )
        )
      ),
      column(
        3,
        wellPanel(
          tags$h3(format_income(b$weighted_income)),
          tags$b("Income context"),
          tags$br(),
          tags$small("Population-weighted average of CBG median incomes")
        )
      ),
      column(
        12,
        tags$p(
          style = "margin:0; color:#555;",
          tags$b("Method: "),
          paste0(
            b$walk_minutes, " minutes × ", walk_speed_m_per_min,
            " m/min = ", round(b$walk_distance_m),
            " m straight-line buffer; ", b$n_cbgs,
            " census block groups intersect; catchment area ",
            round(b$catchment_area_km2, 2), " km²."
          )
        ),
        if (is.finite(b$mean_ej_score)) {
          tags$p(
            style = "margin:2px 0 0 0; color:#555;",
            paste0("Area-weighted SF EJ score within the catchment: ", round(b$mean_ej_score, 1), ".")
          )
        },
        tags$small(
          style = "color:#666;",
          "Population is apportioned by intersected CBG area and therefore assumes residents are uniformly distributed within each block group. The vulnerable-community count uses only SF EJ polygons in the highest-burden quartile (raw score at or above the citywide 75th-percentile cutoff, excluding score 0). This is not a network-routed walkshed."
        )
      )
    )
  })
  
  output$greenspaceScenarioDiffTable <- renderDT({
    comparison <- greenspace_scenario_comparison()
    if (is.null(comparison)) {
      return(DT::datatable(data.frame(Message = "Draw proposed greenspace after generating isochrones.")))
    }
    
    safe_value <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else x[[1]]
    }
    
    base <- comparison$baseline
    scen <- comparison$scenario
    base_raw <- comparison$baseline_raw
    scen_raw <- comparison$scenario_raw
    
    spec <- list(
      list("ACCESS — Transit access (raw stops/km²)", "raw", "Transit_Access_Score", 1),
      list("ACCESS — Transit access percentile", "pct", "pctile_Transit_Access_Score", 1),
      list("ACCESS — Muni routes (raw)", "raw", "Unique_Muni_Routes", 0),
      list("ACCESS — Muni route percentile", "pct", "pctile_Unique_Muni_Routes", 1),
      list("BIODIVERSITY — Total species (raw)", "raw", "GBIF_Species", 0),
      list("BIODIVERSITY — Species richness percentile", "pct", "pctile_GBIF_Species", 1),
      list("BIODIVERSITY — Sampling density (raw records/km²)", "raw", "SamplingDensity_km2", 1),
      list("BIODIVERSITY — Sampling-density percentile", "pct", "pctile_SamplingDensity_km2", 1),
      list("GREEN — Mean NDVI (raw)", "raw", "MeanNDVI", 3),
      list("GREEN — NDVI percentile", "pct", "pctile_MeanNDVI", 1),
      list("GREEN — Greenspace cover (raw %)", "raw", "Greenspace_percent", 1),
      list("GREEN — Greenspace-cover percentile", "pct", "pctile_Greenspace_percent", 1),
      list("EQUITY — SF EJ score (raw)", "raw", "SF_EJ_Score", 1),
      list("EQUITY — Favorable equity-context percentile", "pct", "pctile_SF_EJ_Score", 1),
      list("COMPOSITE — BAI score (/100)", "bai", "BAI", 1),
      list("COMPOSITE — BAI citywide percentile", "pct_bai", "pctile_BAI", 1)
    )
    
    one <- function(s, which_side = c("baseline", "scenario")) {
      which_side <- match.arg(which_side)
      source_raw <- if (which_side == "baseline") base_raw else scen_raw
      source_bai <- if (which_side == "baseline") base else scen
      kind <- s[[2]]
      nm <- s[[3]]
      
      value <- if (kind == "raw") {
        if (nm %in% names(source_raw)) safe_value(source_raw[[nm]]) else NA_real_
      } else if (kind == "bai") {
        if (nm %in% names(source_bai)) 100 * safe_value(source_bai[[nm]]) else NA_real_
      } else {
        if (nm %in% names(source_bai)) safe_value(source_bai[[nm]]) else NA_real_
      }
      value
    }
    
    baseline_values <- vapply(spec, one, numeric(1), which_side = "baseline")
    scenario_values <- vapply(spec, one, numeric(1), which_side = "scenario")
    
    tbl <- data.frame(
      Metric = vapply(spec, function(x) x[[1]], character(1)),
      Baseline = vapply(seq_along(spec), function(i) round(baseline_values[[i]], spec[[i]][[4]]), numeric(1)),
      Scenario = vapply(seq_along(spec), function(i) round(scenario_values[[i]], spec[[i]][[4]]), numeric(1)),
      Difference = vapply(seq_along(spec), function(i) round(scenario_values[[i]] - baseline_values[[i]], spec[[i]][[4]]), numeric(1)),
      check.names = FALSE
    )
    
    DT::datatable(
      tbl,
      rownames = FALSE,
      options = list(dom = "t", ordering = FALSE, paging = FALSE)
    ) |>
      DT::formatStyle("Baseline", color = "#555555", fontWeight = "bold") |>
      DT::formatStyle("Scenario", color = "#238b45", fontWeight = "bold") |>
      DT::formatStyle(
        "Difference",
        fontWeight = "bold",
        backgroundColor = DT::styleInterval(c(-0.05, 0.05), c("#fde0dd", "#f7f7f7", "#e5f5e0"))
      )
  })
  
  output$closestGreenspaceUI <- renderUI({
    df <- base_socio_data()
    if (nrow(df) == 0) return(NULL)
    
    gs_name <- attr(df, "closest_greenspace")
    gs_dist <- attr(df, "closest_greenspace_dist_m")
    
    if (is.null(gs_name) || is.na(gs_name)) gs_name <- "None"
    
    cover_rows <- lapply(seq_len(nrow(df)), function(i) {
      raw_cover <- if ("Greenspace_percent" %in% names(df)) {
        suppressWarnings(as.numeric(df$Greenspace_percent[[i]]))
      } else NA_real_
      pct_cover <- if ("pctile_Greenspace_percent" %in% names(df)) {
        suppressWarnings(as.numeric(df$pctile_Greenspace_percent[[i]]))
      } else NA_real_
      
      tags$div(
        style = "padding:5px 0; border-bottom:1px solid #eee;",
        tags$b(iso_metric_row_label(df, i)),
        tags$br(),
        tags$small(
          "Cover: ",
          if (is.finite(raw_cover)) paste0(round(raw_cover, 1), "%") else "N/A",
          " · percentile: ", format_pct(pct_cover)
        )
      )
    })
    
    tagList(
      strong("Closest Greenspace to selected point:"),
      p(gs_name),
      if (!is.null(gs_dist) && !is.na(gs_dist)) p(paste0("Distance: ", round(gs_dist, 1), " m")),
      hr(),
      strong("Greenspace Cover by Isochrone:"),
      tags$div(style = "max-height:220px; overflow-y:auto;", cover_rows),
      tags$small(
        style = "display:block; color:#666; margin-top:5px;",
        "No averaging across generated isochrones."
      )
    )
  })
  
  output$closestRSFProgramUI <- renderUI({
    df <- base_socio_data()
    if (nrow(df) == 0) return(NULL)
    
    rsf_name <- attr(df, "closest_rsf_program")
    rsf_dist <- attr(df, "closest_rsf_program_dist_m")
    
    if (is.null(rsf_name) || is.na(rsf_name)) rsf_name <- "None"
    
    tagList(
      strong("Nearest RSF program (to your selected point):"),
      p(rsf_name),
      if (!is.null(rsf_dist) && !is.na(rsf_dist)) p(paste0("Distance: ", round(rsf_dist, 1), " m"))
    )
  })
  
  # ---------------------------------------------------------------------------
  # GBIF tab
  # ---------------------------------------------------------------------------
  output$classTable <- renderDT({
    if (is.null(gbif_tab)) {
      return(DT::datatable(data.frame(Message = "GBIF parquet connection not available.")))
    }
    
    q <- gbif_tab
    
    if (input$class_filter != "All") {
      q <- q |> filter(class == input$class_filter)
    }
    if (input$family_filter != "All") {
      q <- q |> filter(family == input$family_filter)
    }
    
    species_counts <- tryCatch({
      q |>
        group_by(species) |>
        summarise(
          n_records = n(),
          .groups = "drop"
        ) |>
        arrange(desc(n_records)) |>
        collect()
    }, error = function(e) NULL)
    
    if (is.null(species_counts) || nrow(species_counts) == 0) {
      return(DT::datatable(data.frame(Message = "No records for that combination.")))
    }
    
    DT::datatable(species_counts, options = list(pageLength = 10), rownames = FALSE)
  })
  
  filtered_data <- reactive({
    data <- cbg_vect_sf
    if ("class" %in% names(data) && input$class_filter != "All") {
      data <- data[data$class == input$class_filter, ]
    }
    if ("family" %in% names(data) && input$family_filter != "All") {
      data <- data[data$family == input$family_filter, ]
    }
    data
  })
  
  output$obsVsSpeciesPlot <- renderPlot({
    data <- filtered_data()
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for selected filters.")
      return(NULL)
    }
    
    ggplot(data, aes(x = log(n_observations + 1), y = log(unique_species + 1))) +
      geom_point(color = "blue", alpha = 0.6) +
      labs(
        x = "Log(Number of Observations + 1)",
        y = "Log(Species Richness + 1)",
        title = "Data Availability vs. Species Richness"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  })
  
  # ---------------------------------------------------------------------------
  # Main plots
  # ---------------------------------------------------------------------------
  output$transitMetricsPlot <- renderPlot({
    df <- base_socio_data()
    if (nrow(df) == 0) return(NULL)
    
    df_long <- df |>
      mutate(IsoLabel = paste0(Mode, "\n", Time, " min")) |>
      select(IsoLabel, Mode, GBIF_Records, GBIF_Species, Transit_Stops, Unique_Muni_Routes) |>
      pivot_longer(
        cols = c(GBIF_Records, GBIF_Species, Transit_Stops, Unique_Muni_Routes),
        names_to = "Metric",
        values_to = "Value"
      ) |>
      mutate(Metric = recode(
        Metric,
        "GBIF_Records"       = "GBIF Records",
        "GBIF_Species"       = "Unique Species",
        "Transit_Stops"      = "Muni Stops (n)",
        "Unique_Muni_Routes" = "Muni Routes (n)"
      ))
    
    ggplot(df_long, aes(x = IsoLabel, y = Value, fill = Mode)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_text(aes(label = Value), vjust = -0.4, size = 3.2, color = "grey30") +
      scale_fill_manual(values = mode_palette, name = "Mode") +
      facet_wrap(~Metric, scales = "free_y", ncol = 2) +
      labs(
        x = NULL,
        y = "Count",
        title = "Biodiversity & Transit Metrics by Transport Mode"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom"
      )
  })
  
  output$radarPlot <- renderPlot({
    df <- baseline_biodiversity_access_index()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    draw_reference_radar(
      df,
      title = "Citywide Biodiversity Access Profile",
      subtitle = "Each line is one independently scored isochrone; farther outward = higher / more favourable percentile.",
      plot_style = "standard"
    )
  })
  
  output$corridorRadarPlot <- renderPlot({
    comparison <- greenspace_scenario_comparison()
    if (is.null(comparison)) return(NULL)
    
    draw_reference_radar(
      bai_rows = list(comparison$baseline, comparison$scenario),
      labels = c("Existing conditions", "Planned conservation action"),
      colors = c("#6b6b6b", "#238b45"),
      title = "Planned Conservation Action — Before vs After",
      subtitle = paste0(
        comparison$label,
        " · dashed = existing conditions · solid green = planned action · Δ labels show changed percentile dimensions"
      ),
      plot_style = "before_after"
    )
  })
  
  # ---------------------------------------------------------------------------
  # Community science tab
  # ---------------------------------------------------------------------------
  community_sf <- safe_partner_orgs()
  
  output$communityMap <- renderLeaflet({
    if (is.null(community_sf) || nrow(community_sf) == 0) {
      return(
        leaflet() |>
          addTiles() |>
          setView(lng = -122.4194, lat = 37.7749, zoom = 12)
      )
    }
    
    name_col <- intersect(c("name", "Name", "organization", "organization_name", "org_name"), names(community_sf))
    popup_col <- if (length(name_col) > 0) name_col[1] else names(community_sf)[1]
    
    leaflet(community_sf) |>
      addTiles() |>
      addCircleMarkers(
        radius = 6,
        color = "#2e8b57",
        fillOpacity = 0.9,
        label = ~as.character(.data[[popup_col]]),
        popup = ~as.character(.data[[popup_col]])
      ) |>
      fitBounds(
        lng1 = min(st_coordinates(st_centroid(st_union(community_sf)))[,1], na.rm = TRUE) - 0.05,
        lat1 = min(st_coordinates(st_centroid(st_union(community_sf)))[,2], na.rm = TRUE) - 0.05,
        lng2 = max(st_coordinates(st_centroid(st_union(community_sf)))[,1], na.rm = TRUE) + 0.05,
        lat2 = max(st_coordinates(st_centroid(st_union(community_sf)))[,2], na.rm = TRUE) + 0.05
      )
  })
  
  output$communityTable <- renderDT({
    if (is.null(community_sf) || nrow(community_sf) == 0) {
      return(DT::datatable(data.frame(Message = "No partner/community organization layer found in setup.R.")))
    }
    
    DT::datatable(
      st_drop_geometry(community_sf),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  # ---------------------------------------------------------------------------
  # Isochrone Comparer tab
  # ---------------------------------------------------------------------------
  # Two point pickers (A and B), one Compare button, two spider plots, and a
  # difference table. Every score reuses the shared iso_metrics_AUDITED.R functions, so
  # a comparison is computed exactly like the main Explorer tab. Nothing here
  # runs until Compare is pressed (cmp_results is an eventReactive).
  
  # One colour per point, reused everywhere (map marker, map isochrone, radar
  # line) so Point A always reads green and Point B always reads blue.
  cmp_color_a <- "#1b9e77"
  cmp_color_b <- "#2c7fb8"
  
  cmp_point_a <- reactiveVal(NULL)
  cmp_point_b <- reactiveVal(NULL)
  
  # Wire one side's map + geocoder to its point reactiveVal. Called once per side
  # so the two pickers stay identical without duplicating the handler code.
  setup_cmp_point <- function(map_id, geocoder_id, point_rv, marker_color) {
    output[[map_id]] <- renderLeaflet({
      leaflet() |>
        addProviderTiles(providers$CartoDB.Positron) |>
        setView(lng = -122.4194, lat = 37.7749, zoom = 12)
    })
    
    observeEvent(input[[paste0(map_id, "_click")]], {
      click <- input[[paste0(map_id, "_click")]]
      if (is.null(click)) return()
      point_rv(c(lon = click$lng, lat = click$lat))
      leafletProxy(map_id) |>
        clearGroup("cmp_pt") |>
        addCircleMarkers(lng = click$lng, lat = click$lat, radius = 7,
                         color = marker_color, group = "cmp_pt")
    })
    
    observeEvent(input[[geocoder_id]], {
      g <- input[[geocoder_id]]
      if (is.null(g)) return()
      xy <- geocoder_as_xy(g)
      point_rv(c(lon = xy[1], lat = xy[2]))
      leafletProxy(map_id) |>
        clearGroup("cmp_pt") |>
        addCircleMarkers(lng = xy[1], lat = xy[2], radius = 7,
                         color = marker_color, group = "cmp_pt") |>
        flyTo(lng = xy[1], lat = xy[2], zoom = 13)
    })
  }
  
  setup_cmp_point("cmp_map_a", "cmp_geocoder_a", cmp_point_a, cmp_color_a)
  setup_cmp_point("cmp_map_b", "cmp_geocoder_b", cmp_point_b, cmp_color_b)
  
  # Build + score one side. Returns NULL if no isochrone could be built (e.g. a
  # transit mode with no reachable stops) so the outputs can show a clear note.
  # Keeps the isochrone sf so it can be drawn back onto the side's map.
  cmp_score_side <- function(point, mode, time) {
    iso <- quiet_sf_overlay(
      build_isochrones(point, modes = mode, times = as.numeric(time))
    )
    if (is.null(iso) || nrow(iso) == 0) return(NULL)
    
    raw_metrics <- quiet_sf_overlay(
      compute_iso_metrics(iso, point, gbif_tab)
    )
    
    # Each comparer side contains one isochrone and is scored only against its
    # own matching mode × time reference distribution. Point A is never ranked
    # directly against Point B.
    metrics <- score_iso_metrics_against_reference(raw_metrics)
    
    list(
      iso     = iso,
      metrics = metrics,
      bai     = add_reference_bai(metrics),
      label   = paste0(pretty_mode(mode), " ", time, " min")
    )
  }
  
  cmp_results <- eventReactive(input$cmp_compare, {
    pa <- cmp_point_a()
    pb <- cmp_point_b()
    
    missing_points <- c(
      if (is.null(pa)) "Point A" else NULL,
      if (is.null(pb)) "Point B" else NULL
    )
    if (length(missing_points) > 0) {
      shiny::showNotification(
        paste0(
          "Place ", paste(missing_points, collapse = " and "),
          " by clicking its map or searching for an address, then press Compare."
        ),
        type = "warning",
        duration = 6
      )
      return(NULL)
    }
    
    list(
      a = cmp_score_side(pa, input$cmp_mode_a, input$cmp_time_a),
      b = cmp_score_side(pb, input$cmp_mode_b, input$cmp_time_b)
    )
  }, ignoreInit = TRUE)
  
  # Flag for the conditionalPanel that reveals the comparison results. FALSE
  # until Compare has produced a value (cmp_results errors before its first
  # event, and on validation failure -- both caught here as "not ready").
  output$cmp_ready <- reactive({
    tryCatch(!is.null(cmp_results()), error = function(e) FALSE)
  })
  outputOptions(output, "cmp_ready", suspendWhenHidden = FALSE)
  
  # Draw each side's isochrone back onto its picker map, and frame it. Fires only
  # when Compare produced results (cmp_results is an eventReactive).
  observeEvent(cmp_results(), {
    res <- cmp_results()
    draw_side_iso <- function(map_id, side, color) {
      proxy <- leafletProxy(map_id) |> clearGroup("cmp_iso")
      if (!is.null(side) && !is.null(side$iso)) {
        bb <- as.numeric(sf::st_bbox(side$iso))
        proxy |>
          addPolygons(data = side$iso, group = "cmp_iso",
                      color = color, weight = 2, fillColor = color, fillOpacity = 0.25) |>
          fitBounds(bb[1], bb[2], bb[3], bb[4])
      }
    }
    draw_side_iso("cmp_map_a", res$a, cmp_color_a)
    draw_side_iso("cmp_map_b", res$b, cmp_color_b)
  })
  
  # One consolidated radar: Point A vs Point B on a single axis, coloured by point.
  output$cmp_radar <- renderPlot({
    res <- cmp_results()
    rows <- list(); labs <- character(); cols <- character()
    if (!is.null(res$a) && !is.null(res$a$bai)) {
      rows <- c(rows, list(res$a$bai)); labs <- c(labs, paste0("Point A: ", res$a$label)); cols <- c(cols, cmp_color_a)
    }
    if (!is.null(res$b) && !is.null(res$b$bai)) {
      rows <- c(rows, list(res$b$bai)); labs <- c(labs, paste0("Point B: ", res$b$label)); cols <- c(cols, cmp_color_b)
    }
    if (length(rows) == 0) {
      plot.new(); title("No isochrones could be built for either point."); return()
    }
    draw_reference_radar(
      rows,
      labels = labs,
      colors = cols,
      title = "Point A vs Point B — Citywide Biodiversity Access",
      subtitle = "Each point is independently benchmarked against its own matching SF mode × travel-time reference.",
      plot_style = "point_compare"
    )
  })
  
  output$cmp_diff_table <- renderDT({
    res <- cmp_results()
    
    side_scores <- function(side) {
      if (is.null(side) || is.null(side$metrics)) return(NULL)
      m <- side$metrics
      b <- side$bai
      
      gm <- function(nm) {
        if (!(nm %in% names(m))) return(NA_real_)
        suppressWarnings(as.numeric(m[[nm]][[1]]))
      }
      gb <- function(nm) {
        if (is.null(b) || !(nm %in% names(b))) return(NA_real_)
        suppressWarnings(as.numeric(b[[nm]][[1]]))
      }
      
      list(
        reference = if ("Reference_Status" %in% names(m)) as.character(m$Reference_Status[[1]]) else "N/A",
        transit_raw = gm("Transit_Access_Score"),
        transit_pct = gm("pctile_Transit_Access_Score"),
        routes_raw = gm("Unique_Muni_Routes"),
        routes_pct = gm("pctile_Unique_Muni_Routes"),
        species_raw = gm("GBIF_Species"),
        species_pct = gm("pctile_GBIF_Species"),
        sampling_raw = gm("SamplingDensity_km2"),
        sampling_pct = gm("pctile_SamplingDensity_km2"),
        ndvi_raw = gm("MeanNDVI"),
        ndvi_pct = gm("pctile_MeanNDVI"),
        greenspace_raw = gm("Greenspace_percent"),
        greenspace_pct = gm("pctile_Greenspace_percent"),
        ej_raw = gm("SF_EJ_Score"),
        ej_pct = gm("pctile_SF_EJ_Score"),
        bai = 100 * gb("BAI"),
        bai_pct = gb("pctile_BAI")
      )
    }
    
    a <- side_scores(res$a)
    b <- side_scores(res$b)
    
    if (is.null(a)) a <- list()
    if (is.null(b)) b <- list()
    
    val <- function(x, digits = 1, suffix = "") {
      if (is.null(x) || length(x) == 0 || !is.finite(x)) return("—")
      paste0(round(x, digits), suffix)
    }
    dint <- function(x) {
      if (is.null(x) || length(x) == 0 || !is.finite(x)) return("—")
      scales::comma(round(x))
    }
    delta <- function(x, y, digits = 1, suffix = "") {
      if (is.null(x) || is.null(y) || length(x) == 0 || length(y) == 0 ||
          !is.finite(x) || !is.finite(y)) return("—")
      d <- y - x
      paste0(if (d >= 0) "+" else "", round(d, digits), suffix)
    }
    
    get <- function(z, nm) if (nm %in% names(z)) z[[nm]] else NA_real_
    
    rows <- list(
      c("Reference group used for percentile scoring", get(a, "reference"), get(b, "reference"), "—"),
      
      c("ACCESS — Transit access (raw stops/km²)", val(get(a,"transit_raw"),2), val(get(b,"transit_raw"),2), delta(get(a,"transit_raw"),get(b,"transit_raw"),2)),
      c("ACCESS — Transit access percentile", val(get(a,"transit_pct")), val(get(b,"transit_pct")), delta(get(a,"transit_pct"),get(b,"transit_pct"),1," pp")),
      c("ACCESS — Muni routes (raw)", dint(get(a,"routes_raw")), dint(get(b,"routes_raw")), delta(get(a,"routes_raw"),get(b,"routes_raw"),0)),
      c("ACCESS — Muni route percentile", val(get(a,"routes_pct")), val(get(b,"routes_pct")), delta(get(a,"routes_pct"),get(b,"routes_pct"),1," pp")),
      
      c("BIODIVERSITY — Total species (raw)", dint(get(a,"species_raw")), dint(get(b,"species_raw")), delta(get(a,"species_raw"),get(b,"species_raw"),0)),
      c("BIODIVERSITY — Species richness percentile", val(get(a,"species_pct")), val(get(b,"species_pct")), delta(get(a,"species_pct"),get(b,"species_pct"),1," pp")),
      c("BIODIVERSITY — Sampling density (raw records/km²)", val(get(a,"sampling_raw"),1), val(get(b,"sampling_raw"),1), delta(get(a,"sampling_raw"),get(b,"sampling_raw"),1)),
      c("BIODIVERSITY — Sampling-density percentile", val(get(a,"sampling_pct")), val(get(b,"sampling_pct")), delta(get(a,"sampling_pct"),get(b,"sampling_pct"),1," pp")),
      
      c("GREEN — Mean NDVI (raw)", val(get(a,"ndvi_raw"),3), val(get(b,"ndvi_raw"),3), delta(get(a,"ndvi_raw"),get(b,"ndvi_raw"),3)),
      c("GREEN — NDVI percentile", val(get(a,"ndvi_pct")), val(get(b,"ndvi_pct")), delta(get(a,"ndvi_pct"),get(b,"ndvi_pct"),1," pp")),
      c("GREEN — Greenspace cover (raw %)", val(get(a,"greenspace_raw"),1,"%"), val(get(b,"greenspace_raw"),1,"%"), delta(get(a,"greenspace_raw"),get(b,"greenspace_raw"),1," pp")),
      c("GREEN — Greenspace-cover percentile", val(get(a,"greenspace_pct")), val(get(b,"greenspace_pct")), delta(get(a,"greenspace_pct"),get(b,"greenspace_pct"),1," pp")),
      
      c("EQUITY — SF EJ score (raw; lower burden is better)", val(get(a,"ej_raw")), val(get(b,"ej_raw")), delta(get(a,"ej_raw"),get(b,"ej_raw"))),
      c("EQUITY — Favorable equity-context percentile", val(get(a,"ej_pct")), val(get(b,"ej_pct")), delta(get(a,"ej_pct"),get(b,"ej_pct"),1," pp")),
      
      c("COMPOSITE — BAI score (/100)", val(get(a,"bai")), val(get(b,"bai")), delta(get(a,"bai"),get(b,"bai"))),
      c("COMPOSITE — BAI citywide percentile", val(get(a,"bai_pct")), val(get(b,"bai_pct")), delta(get(a,"bai_pct"),get(b,"bai_pct"),1," pp"))
    )
    
    tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(tbl) <- c("Metric", "Point A", "Point B", "Difference")
    
    DT::datatable(
      tbl,
      rownames = FALSE,
      options = list(dom = "t", ordering = FALSE, paging = FALSE)
    ) |>
      DT::formatStyle("Point A", color = cmp_color_a, fontWeight = "bold") |>
      DT::formatStyle("Point B", color = cmp_color_b, fontWeight = "bold") |>
      DT::formatStyle("Difference", fontWeight = "bold", backgroundColor = "#fff8e1")
  })
  
} # end server

# =============================================================================
# RUN APP
# =============================================================================
shinyApp(ui, server)
