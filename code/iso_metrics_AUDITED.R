# ============================================================================
# Shared isochrone + metric functions (sourced once by app.R at startup)
# ============================================================================
# These four functions are the single source of truth for turning a point into
# an isochrone, scoring it, and drawing its BAI spider plot. Both the main
# "Isochrone Explorer" tab and the "Isochrone Comparer" tab call them, so the
# scoring logic lives in exactly one place.
#
#   build_isochrones()    point + modes/times      -> isochrone sf (mode,time,geom)
#   compute_iso_metrics() isochrone sf + point     -> per-isochrone metric data.frame
#   draw_radar()          BAI df                   -> spider/radar plot (base graphics)
#
# These functions read the *static* objects loaded by Rscripts/setup_unified.R
# directly as globals (cbg_vect_sf, osm_greenspace, the distance/NDVI rasters,
# rsf_projects, cbg_greenspace_coverage, gtfs_stops_sf, gtfs_routes_sf,
# gtfs_router, transit_iso_cache, cenv_sf, sf_ej_sf) plus the helpers/config
# defined at the top of app.R (mapbox_token, mode_palette, pretty_mode, ecdf01,
# standardize_iso_sf, build_walk_transit_isochrone, safe_biodiv_hotspots/coldspots).
#
# Two objects are created per-session inside server() and so are passed in as
# explicit arguments rather than read as globals:
#   gbif_tab        -- the session's DuckDB handle on the GBIF parquet
#
# compute_iso_metrics() reports progress via withProgress(), so it must be
# called from within a Shiny reactive context (which both tabs satisfy).
# ============================================================================


# Shared scientific definition used by both interactive and Step-5 calculations.
if (!exists("ISO_METRIC_DEFINITION_VERSION", inherits = TRUE)) {
  source("Rscripts/iso_metric_definitions.R", local = TRUE)
}

# ----------------------------------------------------------------------------
# 100 m grid-supported metric cache
# ----------------------------------------------------------------------------
# Every scored raw metric uses the same support as corrected Step 5: the set of
# SF 100 m grid cells whose projected centroids fall within ONE isochrone.
# This makes user and reference metrics mathematically identical at 100 m support.
.iso_metric_cache <- new.env(parent = emptyenv())

load_iso_metric_grid <- function() {
  if (exists("grid", envir = .iso_metric_cache, inherits = FALSE)) {
    return(get("grid", envir = .iso_metric_cache, inherits = FALSE))
  }
  
  path <- Sys.getenv(
    "ISO_METRIC_GRID_GPKG",
    unset = file.path("outputs", "results", "grid_results_100m.gpkg")
  )
  if (!file.exists(path)) {
    stop("Harmonized isochrone metrics require: ", path)
  }
  
  g <- sf::st_read(path, quiet = TRUE) |> sf::st_make_valid()
  required <- c(
    "cell_id", "n_records", "mean_ndvi", "gs_area_m2", "n_stops",
    "nearest_stop_m", "calenviro_ci_score", "calenviro_traffic_pct", "sf_ej_score"
  )
  missing <- setdiff(required, names(g))
  if (length(missing) > 0) {
    stop("grid_results_100m.gpkg lacks required metric columns: ", paste(missing, collapse = ", "))
  }
  
  gp <- sf::st_transform(g, ISO_METRIC_CRS)
  gp$cell_area_m2 <- as.numeric(sf::st_area(gp))
  cp <- sf::st_centroid(gp) |> dplyr::select(cell_id)
  
  out <- list(
    grid_metric = gp,
    centroids_metric = cp,
    attrs = gp |>
      sf::st_drop_geometry() |>
      dplyr::transmute(
        cell_id, cell_area_m2,
        n_records = tidyr::replace_na(as.numeric(n_records), 0),
        mean_ndvi = as.numeric(mean_ndvi),
        gs_area_m2 = tidyr::replace_na(as.numeric(gs_area_m2), 0),
        n_stops = tidyr::replace_na(as.numeric(n_stops), 0),
        nearest_stop_m = as.numeric(nearest_stop_m),
        calenviro_ci_score = as.numeric(calenviro_ci_score),
        calenviro_traffic_pct = as.numeric(calenviro_traffic_pct),
        sf_ej_score = as.numeric(sf_ej_score)
      )
  )
  assign("grid", out, envir = .iso_metric_cache)
  out
}

load_iso_species_cell_lookup <- function() {
  if (exists("species", envir = .iso_metric_cache, inherits = FALSE)) {
    return(get("species", envir = .iso_metric_cache, inherits = FALSE))
  }
  path <- file.path("outputs", "results", "gbif_cell_species_100m.rds")
  if (!file.exists(path)) {
    warning(
      "Missing ", path, ". Run corrected 03_compute_isochrone_basi.R before ",
      "using definition-matched biodiversity percentiles."
    )
    out <- NULL
  } else {
    out <- readRDS(path)
    ver <- attr(out, "iso_metric_definition_version")
    if (is.null(ver) || !identical(ver, ISO_METRIC_DEFINITION_VERSION)) {
      stop("GBIF cell lookup definition version does not match the Shiny metric definition.")
    }
  }
  assign("species", out, envir = .iso_metric_cache)
  out
}

load_iso_route_cell_lookup <- function() {
  if (exists("routes", envir = .iso_metric_cache, inherits = FALSE)) {
    return(get("routes", envir = .iso_metric_cache, inherits = FALSE))
  }
  path <- file.path("outputs", "results", "muni_route_cells_100m.rds")
  if (!file.exists(path)) {
    warning(
      "Missing ", path, ". Run corrected 03_compute_isochrone_basi.R before ",
      "using definition-matched route-access percentiles."
    )
    out <- NULL
  } else {
    out <- readRDS(path)
    ver <- attr(out, "iso_metric_definition_version")
    if (is.null(ver) || !identical(ver, ISO_METRIC_DEFINITION_VERSION)) {
      stop("Muni route-cell lookup definition version does not match the Shiny metric definition.")
    }
  }
  assign("routes", out, envir = .iso_metric_cache)
  out
}

iso_metric_support <- function(poly_i, include_geometry = FALSE) {
  cache <- load_iso_metric_grid()
  poly_metric <- sf::st_transform(sf::st_make_valid(poly_i), ISO_METRIC_CRS)
  hit <- lengths(sf::st_within(cache$centroids_metric, sf::st_geometry(poly_metric))) > 0
  ids <- cache$centroids_metric$cell_id[hit]
  attrs <- cache$attrs[cache$attrs$cell_id %in% ids, , drop = FALSE]
  if (!include_geometry) return(attrs)
  geom <- cache$grid_metric[cache$grid_metric$cell_id %in% ids, , drop = FALSE]
  list(attrs = attrs, geometry = geom)
}

iso_metric_support_geometry <- function(poly_i) {
  support <- iso_metric_support(poly_i, include_geometry = TRUE)
  if (is.null(support$geometry) || nrow(support$geometry) == 0) {
    return(sf::st_sfc(crs = ISO_METRIC_CRS))
  }
  sf::st_union(sf::st_geometry(support$geometry))
}

summarise_iso_metric_support <- function(poly_i) {
  support <- iso_metric_support(poly_i, include_geometry = FALSE)
  if (is.null(support) || nrow(support) == 0) {
    return(list(
      cell_ids = integer(0), n_cells = 0L, area_m2 = NA_real_, area_km2 = NA_real_,
      n_records = NA_real_, n_species = NA_real_, n_birds = NA_real_,
      n_mammals = NA_real_, n_plants = NA_real_, sampling_density = NA_real_,
      mean_ndvi = NA_real_, greenspace_pct = NA_real_, n_stops = NA_real_,
      transit_access = NA_real_, unique_routes = NA_real_, nearest_stop_m = NA_real_,
      calenviro_ci = NA_real_, traffic_pct = NA_real_, sf_ej = NA_real_
    ))
  }
  
  ids <- unique(support$cell_id)
  area_m2 <- sum(support$cell_area_m2[is.finite(support$cell_area_m2)], na.rm = TRUE)
  area_km2 <- area_m2 / 1e6
  n_records <- sum(tidyr::replace_na(support$n_records, 0), na.rm = TRUE)
  n_stops <- sum(tidyr::replace_na(support$n_stops, 0), na.rm = TRUE)
  
  spp <- load_iso_species_cell_lookup()
  if (is.null(spp)) {
    n_species <- n_birds <- n_mammals <- n_plants <- NA_real_
  } else {
    s <- spp[spp$cell_id %in% ids, , drop = FALSE]
    n_species <- dplyr::n_distinct(s$species, na.rm = TRUE)
    n_birds <- dplyr::n_distinct(s$species[s$class == "Aves"], na.rm = TRUE)
    n_mammals <- dplyr::n_distinct(s$species[s$class == "Mammalia"], na.rm = TRUE)
    n_plants <- dplyr::n_distinct(s$species[s$class %in% ISO_PLANT_CLASSES], na.rm = TRUE)
  }
  
  route_lookup <- load_iso_route_cell_lookup()
  unique_routes <- if (is.null(route_lookup)) {
    NA_real_
  } else {
    dplyr::n_distinct(route_lookup$route_id[route_lookup$cell_id %in% ids], na.rm = TRUE)
  }
  
  list(
    cell_ids = ids,
    n_cells = length(ids),
    area_m2 = area_m2,
    area_km2 = area_km2,
    n_records = n_records,
    n_species = n_species,
    n_birds = n_birds,
    n_mammals = n_mammals,
    n_plants = n_plants,
    sampling_density = if (is.finite(area_km2) && area_km2 > 0) n_records / area_km2 else NA_real_,
    mean_ndvi = iso_weighted_mean(support$mean_ndvi, support$cell_area_m2),
    greenspace_pct = if (is.finite(area_m2) && area_m2 > 0) {
      pmin(100, 100 * sum(tidyr::replace_na(support$gs_area_m2, 0), na.rm = TRUE) / area_m2)
    } else NA_real_,
    n_stops = n_stops,
    transit_access = if (is.finite(area_km2) && area_km2 > 0) n_stops / area_km2 else NA_real_,
    unique_routes = unique_routes,
    nearest_stop_m = iso_safe_min(support$nearest_stop_m),
    calenviro_ci = iso_weighted_mean(support$calenviro_ci_score, support$cell_area_m2),
    traffic_pct = iso_weighted_mean(support$calenviro_traffic_pct, support$cell_area_m2),
    sf_ej = iso_weighted_mean(support$sf_ej_score, support$cell_area_m2)
  )
}

# ----------------------------------------------------------------------------
# build_isochrones(): point + chosen modes/times -> combined isochrone sf
# ----------------------------------------------------------------------------
# point: named numeric c(lon = , lat = ) (the shape chosen_point() stores), or NULL.
# modes: character vector from c("driving","walking","cycling","driving-traffic",
#        "transit","walk_transit"). times: numeric minutes. The transit_* / walk_*
# arguments only matter when a transit mode is selected; they default to the
# main tab's defaults so the comparer can pass a single mode/time and ignore them.
# Returns an sf with columns mode, time, geometry (EPSG:4326), or NULL if nothing
# could be built.
build_isochrones <- function(point, modes, times,
                             transit_hour = 9,
                             walk_to_stop_min = 5,
                             walk_from_stop_min = 5,
                             transit_departure_window_min = 10) {
  if (is.null(point) || length(modes) == 0 || length(times) == 0) return(NULL)
  
  location_sf <- st_as_sf(
    data.frame(lon = point["lon"], lat = point["lat"]),
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  iso_list <- list()
  times <- as.numeric(times)
  
  # --- Mapbox modes (driving / walking / cycling / driving-traffic) ----------
  mapbox_modes <- intersect(modes, c("driving", "walking", "cycling", "driving-traffic"))
  for (mode in mapbox_modes) {
    for (t in times) {
      iso <- tryCatch(
        mb_isochrone(location_sf, time = t, profile = mode, access_token = mapbox_token),
        # Surface the Mapbox error to the log instead of silently dropping it --
        # otherwise a failed call just looks like "no isochrone" with no clue why.
        error = function(e) { warning("mb_isochrone failed (", mode, " ", t, " min): ", conditionMessage(e)); NULL }
      )
      if (!is.null(iso)) {
        iso_std <- standardize_iso_sf(iso, mode_name = mode, time_min = t)
        if (!is.null(iso_std)) iso_list <- append(iso_list, list(iso_std))
      }
    }
  }
  
  # --- Transit (GTFS) --------------------------------------------------------
  if ("transit" %in% modes && !is.null(gtfs_router) && !is.null(gtfs_stops_sf)) {
    stop_dists  <- st_distance(location_sf, gtfs_stops_sf)
    nearest_idx <- which.min(stop_dists)
    nearest_id  <- as.character(gtfs_stops_sf$stop_id[nearest_idx])
    dep_secs    <- as.numeric(transit_hour) * 3600
    
    for (t in times) {
      iso_poly <- NULL
      
      if (!is.null(transit_iso_cache) &&
          !is.null(transit_iso_cache[[nearest_id]]) &&
          !is.null(transit_iso_cache[[nearest_id]][[as.character(t)]])) {
        iso_poly <- transit_iso_cache[[nearest_id]][[as.character(t)]]
      }
      
      if (is.null(iso_poly)) {
        iso_result <- tryCatch(
          gtfsrouter::gtfs_isochrone(
            gtfs       = gtfs_router,
            from       = nearest_id,
            start_time = dep_secs,
            end_time   = dep_secs + t * 60,
            from_is_id = TRUE
          ),
          error = function(e) NULL
        )
        
        if (!is.null(iso_result) && nrow(iso_result) > 2) {
          reachable_sf <- gtfs_stops_sf |>
            filter(stop_id %in% as.character(iso_result$stop_id))
          
          if (nrow(reachable_sf) > 2) {
            iso_poly <- st_convex_hull(st_union(reachable_sf))
          } else if (nrow(reachable_sf) > 0) {
            iso_poly <- st_union(st_buffer(st_transform(reachable_sf, 3857), 100)) |>
              st_transform(4326)
          }
        }
      }
      
      if (!is.null(iso_poly)) {
        iso_sf <- st_sf(
          mode = "transit",
          time = as.numeric(t),
          geometry = st_geometry(st_as_sf(iso_poly)),
          crs = 4326
        )
        iso_sf <- standardize_iso_sf(iso_sf, mode_name = "transit", time_min = t)
        iso_list <- append(iso_list, list(iso_sf))
      }
    }
  }
  
  # --- Walk + Transit (Muni) -------------------------------------------------
  if ("walk_transit" %in% modes && !is.null(gtfs_router) && !is.null(gtfs_stops_sf)) {
    dep_secs    <- as.numeric(transit_hour) * 3600
    valid_times <- times[times > walk_to_stop_min]
    
    for (t in valid_times) {
      wt_iso <- tryCatch(
        build_walk_transit_isochrone(
          location_sf            = location_sf,
          total_time_min         = t,
          dep_secs               = dep_secs,
          walk_to_stop_min       = walk_to_stop_min,
          walk_from_stop_min     = walk_from_stop_min,
          gtfs_stops_sf          = gtfs_stops_sf,
          gtfs_router            = gtfs_router,
          mapbox_token           = mapbox_token,
          departure_window_min   = transit_departure_window_min,
          departure_step_min     = 5,
          max_last_mile_stops    = 12,
          include_first_mile_polygon = TRUE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(wt_iso) && nrow(wt_iso) > 0) {
        iso_list <- append(iso_list, list(wt_iso))
      }
    }
  }
  
  if (length(iso_list) == 0) return(NULL)
  
  dplyr::bind_rows(iso_list) |>
    st_as_sf() |>
    st_make_valid() |>
    st_transform(4326)
}


# ----------------------------------------------------------------------------
# compute_iso_metrics(): one independently calculated row per isochrone
# ----------------------------------------------------------------------------
# Scored metrics use the shared 100 m centroid-support definition above. This is
# intentionally identical to corrected Step 5. `gbif_tab` is retained in the
# signature for compatibility but is no longer used for scored richness/effort.
compute_iso_metrics <- function(iso_data, point, gbif_tab) {
  if (is.null(iso_data) || nrow(iso_data) == 0) return(data.frame())
  
  hotspot_union <- safe_biodiv_hotspots()
  coldspot_union <- safe_biodiv_coldspots()
  if (!is.null(hotspot_union)) hotspot_union <- st_union(hotspot_union)
  if (!is.null(coldspot_union)) coldspot_union <- st_union(coldspot_union)
  
  acs_wide <- cbg_vect_sf |>
    mutate(population = popE, med_income = medincE)
  
  user_point_sf <- NULL
  if (!is.null(point)) {
    user_point_sf <- st_as_sf(
      data.frame(lon = point["lon"], lat = point["lat"]),
      coords = c("lon", "lat"), crs = 4326
    )
  }
  
  min_dist_val_global <- NA_real_
  osm_greenspace_name_global <- NA_character_
  if (!is.null(user_point_sf) && exists("greenspace_dist_raster") && exists("greenspace_osmid_raster")) {
    try({
      min_dist_val_global <- (greenspace_dist_raster |> extract(vect(user_point_sf)) |> pull(2))[1]
      user_point_osm_id <- (greenspace_osmid_raster |> extract(vect(user_point_sf)) |> pull(2))[1]
      osm_greenspace_name_global <- osm_greenspace |>
        mutate(osm_id = as.numeric(osm_id)) |>
        filter(osm_id == user_point_osm_id) |>
        pull(name)
      if (length(osm_greenspace_name_global) == 0 || is.na(osm_greenspace_name_global[1])) {
        osm_greenspace_name_global <- "Unnamed Greenspace"
      } else {
        osm_greenspace_name_global <- osm_greenspace_name_global[1]
      }
    }, silent = TRUE)
  }
  
  min_rsf_dist_global <- NA_real_
  rsf_program_name_global <- NA_character_
  if (!is.null(user_point_sf) && exists("rsfprogram_dist_raster") && exists("rsfprogram_id_raster") && exists("rsf_projects")) {
    try({
      min_rsf_dist_global <- (rsfprogram_dist_raster |> extract(vect(user_point_sf)) |> pull(2))[1]
      user_point_rsf_pid <- (rsfprogram_id_raster |> extract(vect(user_point_sf)) |> pull(2))[1]
      rsf_program_name_global <- rsf_projects |>
        dplyr::filter(as.numeric(.data$polygon_id) == as.numeric(user_point_rsf_pid)) |>
        dplyr::pull(.data$prj_name)
      if (length(rsf_program_name_global) == 0 || is.na(rsf_program_name_global[1])) {
        rsf_program_name_global <- "Unknown RSF program"
      } else {
        rsf_program_name_global <- rsf_program_name_global[1]
      }
    }, silent = TRUE)
  }
  
  results <- data.frame()
  n_isos <- nrow(iso_data)
  
  withProgress(message = "Analyzing isochrones...", value = 0, {
    for (i in seq_len(n_isos)) {
      poly_i <- iso_data[i, ]
      incProgress(
        1 / n_isos,
        detail = paste0(pretty_mode(as.character(poly_i$mode[[1]])), " – ", poly_i$time[[1]], " min")
      )
      
      support <- summarise_iso_metric_support(poly_i)
      iso_area_m2 <- as.numeric(st_area(st_transform(poly_i, ISO_METRIC_CRS)))
      iso_area_km2 <- iso_area_m2 / 1e6
      
      dist_hot_km <- if (!is.null(hotspot_union)) {
        round(as.numeric(min(st_distance(poly_i, hotspot_union))) / 1000, 3)
      } else NA_real_
      dist_cold_km <- if (!is.null(coldspot_union)) {
        round(as.numeric(min(st_distance(poly_i, coldspot_union))) / 1000, 3)
      } else NA_real_
      
      # Population/income remain descriptive, non-BAI fields and therefore retain
      # the existing polygon-intersection implementation.
      inter_acs <- tryCatch(intersect(vect(acs_wide), vect(poly_i)) |> st_as_sf(), error = function(e) NULL)
      pop_total <- 0
      w_income <- NA_real_
      if (!is.null(inter_acs) && nrow(inter_acs) > 0) {
        inter_acs <- inter_acs |>
          mutate(
            area_num = as.numeric(st_area(st_transform(geometry, ISO_METRIC_CRS))),
            weighted_pop = population * (area_num / sum(area_num, na.rm = TRUE))
          )
        pop_total <- round(sum(inter_acs$weighted_pop, na.rm = TRUE))
        w_income <- sum(inter_acs$med_income * inter_acs$area_num, na.rm = TRUE) /
          sum(inter_acs$area_num, na.rm = TRUE)
      }
      
      # Frequency/headway are descriptive only; they are not reference-scored.
      freq_weighted_score <- NA_real_
      mean_headway_iso <- NA_real_
      nearest_stop_name <- NA_character_
      if (!is.null(gtfs_stops_sf)) {
        inter_transit <- tryCatch(st_intersection(gtfs_stops_sf, poly_i), error = function(e) NULL)
        if (!is.null(inter_transit) && nrow(inter_transit) > 0 &&
            "mean_headway_min" %in% names(inter_transit)) {
          hw <- inter_transit$mean_headway_min
          hw <- hw[is.finite(hw) & hw > 0]
          if (length(hw) > 0 && is.finite(support$area_km2) && support$area_km2 > 0) {
            freq_weighted_score <- round(sum(60 / hw) / support$area_km2, 2)
            mean_headway_iso <- round(mean(hw), 1)
          }
        }
        if (is.finite(support$nearest_stop_m)) {
          support_geom <- iso_metric_support_geometry(poly_i)
          if (length(support_geom) > 0) {
            support_sf <- st_sf(geometry = support_geom)
            d <- tryCatch(st_distance(support_sf, st_transform(gtfs_stops_sf, st_crs(support_sf))), error = function(e) NULL)
            if (!is.null(d) && length(d) > 0) {
              nearest_stop_name <- gtfs_stops_sf$stop_name[which.min(as.numeric(d))]
            }
          }
        }
      }
      
      row_i <- data.frame(
        Mode = pretty_mode(as.character(poly_i$mode[[1]])),
        Time = as.numeric(poly_i$time[[1]]),
        IsochroneArea_km2 = round(iso_area_km2, 3),
        MetricSupportArea_km2 = round(support$area_km2, 3),
        MetricSupportNCells = as.integer(support$n_cells),
        MetricDefinitionVersion = ISO_METRIC_DEFINITION_VERSION,
        DistToHotspot_km = dist_hot_km,
        DistToColdspot_km = dist_cold_km,
        EstimatedPopulation = pop_total,
        MedianIncome = round(w_income, 2),
        MeanNDVI = as.numeric(support$mean_ndvi),
        GBIF_Records = as.numeric(support$n_records),
        GBIF_Species = as.numeric(support$n_species),
        Bird_Species = as.numeric(support$n_birds),
        Mammal_Species = as.numeric(support$n_mammals),
        Plant_Species = as.numeric(support$n_plants),
        SamplingDensity_km2 = as.numeric(support$sampling_density),
        Greenspace_percent = as.numeric(support$greenspace_pct),
        Transit_Stops = as.numeric(support$n_stops),
        Unique_Muni_Routes = as.numeric(support$unique_routes),
        Transit_Access_Score = as.numeric(support$transit_access),
        Freq_Weighted_Score = freq_weighted_score,
        Mean_Headway_min = mean_headway_iso,
        Nearest_Stop_m = as.numeric(support$nearest_stop_m),
        Nearest_Stop_Name = nearest_stop_name,
        CalEnviro_CIscore = as.numeric(support$calenviro_ci),
        CalEnviro_Traffic_Pctl = as.numeric(support$traffic_pct),
        SF_EJ_Score = as.numeric(support$sf_ej),
        closest_greenspace = osm_greenspace_name_global,
        closest_greenspace_dist_m = min_dist_val_global,
        closest_rsf_program = rsf_program_name_global,
        closest_rsf_program_dist_m = min_rsf_dist_global,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, row_i)
    }
  })
  
  # No union-based or mean-across-isochrone scoring attributes are created here.
  # Each row remains the sole unit of percentile and BAI scoring.
  if (nrow(results) > 0) {
    closest_gs <- results |>
      filter(!is.na(closest_greenspace_dist_m)) |>
      slice_min(closest_greenspace_dist_m, n = 1)
    attr(results, "closest_greenspace") <- if (nrow(closest_gs) > 0) closest_gs$closest_greenspace[1] else "None"
    attr(results, "closest_greenspace_dist_m") <- if (nrow(closest_gs) > 0) closest_gs$closest_greenspace_dist_m[1] else NA_real_
    
    closest_rsf <- results |>
      filter(!is.na(closest_rsf_program_dist_m)) |>
      slice_min(closest_rsf_program_dist_m, n = 1)
    attr(results, "closest_rsf_program") <- if (nrow(closest_rsf) > 0) closest_rsf$closest_rsf_program[1] else "None"
    attr(results, "closest_rsf_program_dist_m") <- if (nrow(closest_rsf) > 0) closest_rsf$closest_rsf_program_dist_m[1] else NA_real_
  }
  
  results
}


# Legacy census-block-group add_bai() scoring was removed. The app now uses
# add_reference_bai() from app.R, which requires same-mode × same-time
# precomputed isochrone references and never falls back to CBG ECDFs.


# Axis order shared by both radar functions, and the column each axis reads from.
# Keep these two vectors aligned -- the group labels below assume this clockwise
# order (1-2 Urban Access, 3-4 Biodiversity, 5-6 Environment, 7 EJ).
RADAR_AXIS_COLS <- c("Mobility_Access_std", "Route_Access_std", "Biodiversity_Potential_std",
                     "Observation_Intensity_std", "Environmental_Quality_std",
                     "Greenspace_Cover_std", "Equity_Context_std")
RADAR_AXIS_LABELS <- c("Stop\nDensity", "Route\nDiversity", "Species\nRichness",
                       "Obs.\nIntensity", "Vegetation\n(NDVI)", "Greenspace\nCover", "EJ\nContext")

# Radial group subheaders placed just beyond the axis labels. Shared by both
# radar functions so the grouping (Urban Access / Biodiversity / Environment /
# EJ) is drawn identically. Call after fmsb::radarchart() has drawn the chart.
draw_radar_group_labels <- function() {
  n_ax  <- 7
  angs  <- pi/2 - (0:(n_ax - 1)) * (2 * pi / n_ax)
  cat_r <- 1.48  # radius just outside axis labels (~1.2-1.3)
  
  text(cat_r * cos(mean(angs[1:2])), cat_r * sin(mean(angs[1:2])),
       "Urban Access",           cex = 0.72, col = "#2166ac", font = 2, xpd = TRUE)
  text(cat_r * cos(mean(angs[3:4])), cat_r * sin(mean(angs[3:4])),
       "Biodiversity",           cex = 0.72, col = "#1b7837", font = 2, xpd = TRUE)
  text(cat_r * cos(mean(angs[5:6])), cat_r * sin(mean(angs[5:6])),
       "Environment",            cex = 0.72, col = "#762a83", font = 2, xpd = TRUE)
  text(cat_r * cos(angs[7]),         cat_r * sin(angs[7]),
       "Environmental\nJustice", cex = 0.72, col = "#b2182b", font = 2, xpd = TRUE)
}

# ----------------------------------------------------------------------------
# draw_radar(): BAI df -> spider/radar plot (base graphics, via fmsb)
# ----------------------------------------------------------------------------
# bai_df: the data.frame returned by same-mode × same-time reference scoring (one row per isochrone). title:
# the chart title. Draws directly to the current graphics device, so call it
# from inside a renderPlot({ }). Lines are coloured by transport mode.
draw_radar <- function(bai_df, title = "Biodiversity Access Index Profile") {
  if (is.null(bai_df) || nrow(bai_df) == 0) return(NULL)
  
  radar_df <- bai_df |>
    mutate(
      ModeTime = paste0(Mode, "_", Time, "m"),
      `Stop\nDensity`     = Mobility_Access_std,
      `Route\nDiversity`  = Route_Access_std,
      `Species\nRichness` = Biodiversity_Potential_std,
      `Obs.\nIntensity`   = Observation_Intensity_std,
      `Vegetation\n(NDVI)`= Environmental_Quality_std,
      `Greenspace\nCover` = Greenspace_Cover_std,
      `EJ\nContext`       = Equity_Context_std
    ) |>
    select(
      ModeTime,
      `Stop\nDensity`,
      `Route\nDiversity`,
      `Species\nRichness`,
      `Obs.\nIntensity`,
      `Vegetation\n(NDVI)`,
      `Greenspace\nCover`,
      `EJ\nContext`
    )
  
  radar_mat <- as.data.frame(radar_df[, -1])
  rownames(radar_mat) <- radar_df$ModeTime
  
  radar_mat <- rbind(
    rep(1, ncol(radar_mat)),
    rep(0, ncol(radar_mat)),
    radar_mat
  )
  
  labels <- rownames(radar_mat)[-(1:2)]
  line_cols <- mode_palette[gsub("_(.*)$", "", labels)]
  line_cols[is.na(line_cols)] <- "#666666"
  
  # Extra margin so two-line axis labels aren't clipped
  par(mar = c(2, 2, 3, 2))
  
  fmsb::radarchart(
    radar_mat,
    axistype = 1,
    pcol = line_cols,
    plwd = 2,
    plty = 1,
    cglcol = "grey80",
    cglty = 1,
    cglwd = 0.8,
    axislabcol = "grey40",
    vlcex = 0.88,
    title = title
  )
  
  draw_radar_group_labels()
  
  legend(
    "topright",
    legend = labels,
    col = line_cols,
    lty = 1,
    lwd = 2,
    cex = 0.75,
    bty = "n"
  )
}


# ----------------------------------------------------------------------------
# draw_compare_radar(): overlay several BAI rows on ONE radar (Comparer tab)
# ----------------------------------------------------------------------------
# bai_rows: a list of single-row reference-scored BAI data.frames (one per point).
# labels:   legend text, one per row (e.g. "Point A: Walking 5 min").
# colors:   one line colour per row -- this is where the A-vs-B colour scheme
#           lives, so keep it in sync with the point markers / isochrones.
# Unlike draw_radar() (coloured by mode), this colours by *point* so two
# isochrones of the same mode are still distinguishable. Draws to the current
# device; call from inside renderPlot({ }).
draw_compare_radar <- function(bai_rows, labels, colors,
                               title = "Biodiversity Access Index — Point A vs Point B") {
  if (length(bai_rows) == 0) return(NULL)
  
  mat <- as.data.frame(do.call(rbind, lapply(bai_rows, function(d) as.numeric(d[1, RADAR_AXIS_COLS]))))
  names(mat) <- RADAR_AXIS_LABELS
  # fmsb wants the axis max (1) and min (0) as the first two rows.
  mat <- rbind(rep(1, ncol(mat)), rep(0, ncol(mat)), mat)
  
  par(mar = c(2, 2, 3, 2))
  fmsb::radarchart(
    mat,
    axistype = 1,
    pcol  = colors,
    pfcol = grDevices::adjustcolor(colors, alpha.f = 0.2),  # translucent fills so overlaps read
    plwd  = 3,
    plty  = 1,
    cglcol = "grey80",
    cglty = 1,
    cglwd = 0.8,
    axislabcol = "grey40",
    vlcex = 0.9,
    title = title
  )
  
  draw_radar_group_labels()
  
  legend(
    "topright",
    legend = labels,
    col = colors,
    lty = 1,
    lwd = 3,
    cex = 0.85,
    bty = "n"
  )
}