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
#   add_bai()             metric df + benchmarks   -> df + 7 BAI axes + composite BAI
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
#   city_benchmarks -- the citywide ECDF reference distributions for the BAI
#
# compute_iso_metrics() reports progress via withProgress(), so it must be
# called from within a Shiny reactive context (which both tabs satisfy).
# ============================================================================


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
# compute_iso_metrics(): isochrone sf -> per-isochrone metric data.frame
# ----------------------------------------------------------------------------
# iso_data: the sf returned by build_isochrones(). point: the c(lon, lat) the
# isochrones were built from (used for the nearest-greenspace / nearest-RSF
# lookups), or NULL. gbif_tab: the session's DuckDB handle on the GBIF parquet.
# Returns a data.frame (one row per isochrone) carrying aggregate summaries as
# attributes: bio_percentile, city/mean transit scores, closest_greenspace(+dist),
# closest_rsf_program(+dist).
compute_iso_metrics <- function(iso_data, point, gbif_tab) {
  if (is.null(iso_data) || nrow(iso_data) == 0) return(data.frame())

  hotspot_union <- safe_biodiv_hotspots()
  coldspot_union <- safe_biodiv_coldspots()

  if (!is.null(hotspot_union)) hotspot_union <- st_union(hotspot_union)
  if (!is.null(coldspot_union)) coldspot_union <- st_union(coldspot_union)

  acs_wide <- cbg_vect_sf |>
    mutate(population = popE, med_income = medincE)

  results <- data.frame()
  n_isos <- nrow(iso_data)

  user_point_sf <- NULL
  if (!is.null(point)) {
    user_point_sf <- st_as_sf(
      data.frame(lon = point["lon"], lat = point["lat"]),
      coords = c("lon", "lat"),
      crs = 4326
    )
  }

  min_dist_val_global <- NA_real_
  osm_greenspace_name_global <- NA_character_
  if (!is.null(user_point_sf) && exists("greenspace_dist_raster") && exists("greenspace_osmid_raster")) {
    try({
      # Distance from the selected point to its nearest greenspace, read straight
      # from the distance raster. terra::extract() returns data.frame(ID, value),
      # so the distance is column 2 -- pulling column 1 (the point ID) was the bug
      # that made this come back as "1" instead of the real distance in metres.
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
      # Distance is column 2 (column 1 is the point ID) -- same fix as greenspace above.
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

  withProgress(message = "Analyzing isochrones...", value = 0, {

  for (i in seq_len(n_isos)) {
    poly_i <- iso_data[i, ]
    vect_poly_i <- vect(poly_i)

    incProgress(
      1 / n_isos,
      detail = paste0(
        pretty_mode(as.character(poly_i$mode[[1]])),
        " – ", poly_i$time[[1]], " min"
      )
    )

    dist_hot_km <- if (!is.null(hotspot_union)) {
      round(as.numeric(min(st_distance(poly_i, hotspot_union))) / 1000, 3)
    } else NA_real_

    dist_cold_km <- if (!is.null(coldspot_union)) {
      round(as.numeric(min(st_distance(poly_i, coldspot_union))) / 1000, 3)
    } else NA_real_

    inter_acs <- tryCatch(intersect(vect(acs_wide), vect_poly_i) |> st_as_sf(), error = function(e) NULL)

    pop_total <- 0
    w_income  <- NA_real_
    if (!is.null(inter_acs) && nrow(inter_acs) > 0) {
      inter_acs <- inter_acs |>
        mutate(
          area_num = as.numeric(st_area(st_transform(geometry, 3857))),
          weighted_pop = population * (area_num / sum(area_num, na.rm = TRUE))
        )
      pop_total <- round(sum(inter_acs$weighted_pop, na.rm = TRUE))
      w_income  <- sum(inter_acs$med_income * inter_acs$area_num, na.rm = TRUE) /
        sum(inter_acs$area_num, na.rm = TRUE)
    }

    iso_area_m2 <- as.numeric(st_area(st_transform(poly_i, 3857)))
    iso_area_km2 <- round(iso_area_m2 / 1e6, 3)

    # Greenspace area within this isochrone: use pre-cached per-CBG coverage
    # (replaces slow per-isochrone raster cellSize() which took ~40-55 s each).
    # Method: intersect isochrone with CBGs, scale each CBG's greenspace_m2 by
    # the fraction of that CBG falling inside the isochrone, then sum.
    gs_area_m2 <- tryCatch({
      if (exists("cbg_greenspace_coverage") && !is.null(cbg_greenspace_coverage)) {
        poly_proj  <- st_transform(poly_i, 3857)
        cbg_proj_i <- st_transform(cbg_vect_sf[, "GEOID"], 3857) |>
          mutate(cbg_area_m2 = as.numeric(st_area(geometry))) |>
          st_make_valid()
        inter_df <- st_intersection(cbg_proj_i, poly_proj) |>
          mutate(inter_area_m2 = as.numeric(st_area(geometry))) |>
          st_drop_geometry() |>
          left_join(
            cbg_greenspace_coverage[, c("GEOID", "greenspace_m2", "cbg_area_m2")],
            by = "GEOID", suffix = c("_iso", "_cbg")
          ) |>
          mutate(
            greenspace_m2  = tidyr::replace_na(greenspace_m2, 0),
            cbg_area_m2    = dplyr::coalesce(cbg_area_m2_cbg, cbg_area_m2_iso),
            contrib        = ifelse(cbg_area_m2 > 0,
                                    inter_area_m2 / cbg_area_m2 * greenspace_m2,
                                    0)
          )
        sum(inter_df$contrib, na.rm = TRUE)
      } else if (exists("greenspace_dist_raster")) {
        # Fallback to raster method if cache is unavailable
        dist_crop <- terra::crop(greenspace_dist_raster, vect_poly_i)
        dist_mask <- terra::mask(dist_crop, vect_poly_i)
        is_greenspace <- dist_mask == 0
        cell_areas <- terra::cellSize(is_greenspace, unit = "m")
        as.numeric(terra::global(cell_areas * is_greenspace, "sum", na.rm = TRUE)[1, 1])
      } else {
        0
      }
    }, error = function(e) 0)
    gs_percent <- ifelse(iso_area_m2 > 0, 100 * gs_area_m2 / iso_area_m2, 0)

    mean_ndvi <- NA_real_
    if (exists("ndvi")) {
      ndvi_vals <- tryCatch(values(terra::mask(terra::crop(ndvi, vect_poly_i), vect_poly_i)), error = function(e) NA)
      ndvi_vals <- ndvi_vals[!is.na(ndvi_vals)]
      mean_ndvi <- ifelse(length(ndvi_vals) > 0, round(mean(ndvi_vals, na.rm = TRUE), 3), NA_real_)
    }

    iso_wkt <- st_as_text(st_geometry(poly_i)[[1]])
    gbif_summary <- tryCatch({
      gbif_tab |>
        filter(sql(glue("ST_Intersects(ST_GeomFromText(geom_wkt), ST_GeomFromText('{iso_wkt}'))"))) |>
        summarise(
          n_records = n(),
          n_species = n_distinct(species),
          n_birds   = n_distinct(case_when(class == "Aves"     ~ species, TRUE ~ NA_character_)),
          n_mammals = n_distinct(case_when(class == "Mammalia" ~ species, TRUE ~ NA_character_)),
          n_plants  = n_distinct(case_when(
            class %in% c("Magnoliopsida", "Liliopsida", "Pinopsida", "Polypodiopsida",
                         "Equisetopsida", "Bryopsida",  "Marchantiopsida") ~ species,
            TRUE ~ NA_character_
          ))
        ) |>
        collect()
    }, error = function(e) NULL)
    n_records <- if (!is.null(gbif_summary) && nrow(gbif_summary) > 0) gbif_summary$n_records[[1]] else 0L
    n_species <- if (!is.null(gbif_summary) && nrow(gbif_summary) > 0) gbif_summary$n_species[[1]] else 0L
    n_birds   <- if (!is.null(gbif_summary) && nrow(gbif_summary) > 0) gbif_summary$n_birds[[1]]   else 0L
    n_mammals <- if (!is.null(gbif_summary) && nrow(gbif_summary) > 0) gbif_summary$n_mammals[[1]] else 0L
    n_plants  <- if (!is.null(gbif_summary) && nrow(gbif_summary) > 0) gbif_summary$n_plants[[1]]  else 0L

    n_transit_stops <- NA_real_
    transit_access_score <- NA_real_
    freq_weighted_score <- NA_real_
    mean_headway_iso <- NA_real_
    nearest_stop_m <- NA_real_
    nearest_stop_name <- NA_character_

    if (!is.null(gtfs_stops_sf)) {
      inter_transit <- tryCatch(st_intersection(gtfs_stops_sf, poly_i), error = function(e) NULL)
      n_transit_stops <- if (!is.null(inter_transit)) nrow(inter_transit) else 0

      dist_transit <- st_distance(poly_i, gtfs_stops_sf)
      nearest_stop_m <- round(as.numeric(min(dist_transit)), 0)
      nearest_stop_name <- gtfs_stops_sf$stop_name[which.min(dist_transit)]

      transit_access_score <- ifelse(iso_area_km2 > 0, round(n_transit_stops / iso_area_km2, 2), NA_real_)

      freq_weighted_score <- if (!is.null(inter_transit) &&
                                 "mean_headway_min" %in% names(inter_transit) &&
                                 nrow(inter_transit) > 0 &&
                                 iso_area_km2 > 0) {
        hw <- inter_transit$mean_headway_min
        hw <- hw[!is.na(hw) & hw > 0]
        if (length(hw) > 0) round(sum(60 / hw) / iso_area_km2, 2) else NA_real_
      } else NA_real_

      mean_headway_iso <- if (!is.null(inter_transit) &&
                              "mean_headway_min" %in% names(inter_transit) &&
                              nrow(inter_transit) > 0) {
        round(mean(inter_transit$mean_headway_min, na.rm = TRUE), 1)
      } else NA_real_
    }

    # Unique Muni route IDs whose shapes intersect this isochrone
    n_unique_routes <- 0L
    if (!is.null(gtfs_routes_sf)) {
      routes_inter <- tryCatch(
        st_intersection(gtfs_routes_sf[, "route_id"], poly_i),
        error = function(e) NULL
      )
      n_unique_routes <- if (!is.null(routes_inter) && "route_id" %in% names(routes_inter)) {
        length(unique(routes_inter$route_id))
      } else 0L
    }

    sampling_density_km2 <- ifelse(iso_area_km2 > 0, round(n_records / iso_area_km2, 2), NA_real_)

    mean_ciscore <- if (!is.null(cenv_sf)) {
      tryCatch({
        ce_inter <- st_intersection(cenv_sf, poly_i)
        if (nrow(ce_inter) > 0) {
          ce_inter$a <- as.numeric(st_area(st_transform(ce_inter, 3857)))
          round(weighted.mean(ce_inter$CIscore, w = ce_inter$a, na.rm = TRUE), 1)
        } else NA_real_
      }, error = function(e) NA_real_)
    } else NA_real_

    mean_traffic_pctl <- if (!is.null(cenv_sf)) {
      tryCatch({
        ce_inter <- st_intersection(cenv_sf, poly_i)
        if (nrow(ce_inter) > 0) {
          ce_inter$a <- as.numeric(st_area(st_transform(ce_inter, 3857)))
          round(weighted.mean(ce_inter$Traffic_Pctl, w = ce_inter$a, na.rm = TRUE), 1)
        } else NA_real_
      }, error = function(e) NA_real_)
    } else NA_real_

    mean_ej_score <- if (!is.null(sf_ej_sf)) {
      tryCatch({
        ej_inter <- st_intersection(sf_ej_sf, poly_i)
        if (nrow(ej_inter) > 0) {
          ej_inter$a <- as.numeric(st_area(st_transform(ej_inter, 3857)))
          valid <- ej_inter[!is.na(ej_inter$score), ]
          if (nrow(valid) > 0) {
            round(weighted.mean(valid$score, w = valid$a, na.rm = TRUE), 1)
          } else NA_real_
        } else NA_real_
      }, error = function(e) NA_real_)
    } else NA_real_

    row_i <- data.frame(
      Mode                   = pretty_mode(as.character(poly_i$mode[[1]])),
      Time                   = as.numeric(poly_i$time[[1]]),
      IsochroneArea_km2      = iso_area_km2,
      DistToHotspot_km       = dist_hot_km,
      DistToColdspot_km      = dist_cold_km,
      EstimatedPopulation    = pop_total,
      MedianIncome           = round(w_income, 2),
      MeanNDVI               = mean_ndvi,
      GBIF_Records           = n_records,
      GBIF_Species           = n_species,
      Bird_Species           = n_birds,
      Mammal_Species         = n_mammals,
      Plant_Species          = n_plants,
      SamplingDensity_km2    = sampling_density_km2,
      Greenspace_percent     = round(gs_percent, 2),
      Transit_Stops          = n_transit_stops,
      Unique_Muni_Routes     = n_unique_routes,
      Transit_Access_Score   = transit_access_score,
      Freq_Weighted_Score    = freq_weighted_score,
      Mean_Headway_min       = mean_headway_iso,
      Nearest_Stop_m         = nearest_stop_m,
      Nearest_Stop_Name      = nearest_stop_name,
      CalEnviro_CIscore      = mean_ciscore,
      CalEnviro_Traffic_Pctl = mean_traffic_pctl,
      SF_EJ_Score            = mean_ej_score,
      closest_greenspace     = osm_greenspace_name_global,
      closest_greenspace_dist_m = min_dist_val_global,
      closest_rsf_program    = rsf_program_name_global,
      closest_rsf_program_dist_m = min_rsf_dist_global,
      stringsAsFactors       = FALSE
    )

    results <- rbind(results, row_i)
  }

  }) # end withProgress

  union_wkt <- st_as_text(st_geometry(st_union(iso_data))[[1]])
  union_n_species <- tryCatch({
    gbif_tab |>
      filter(sql(glue("ST_Intersects(ST_GeomFromText(geom_wkt), ST_GeomFromText('{union_wkt}'))"))) |>
      summarise(n_species = n_distinct(species)) |>
      collect() |>
      pull(n_species)
  }, error = function(e) 0L)

  attr(results, "bio_percentile") <- round(100 * ecdf(cbg_vect_sf$unique_species)(union_n_species), 1)

  if (!is.null(gtfs_stops_sf)) {
    sf_city_area_km2 <- 121.4
    attr(results, "city_transit_score") <- round(nrow(gtfs_stops_sf) / sf_city_area_km2, 2)
    attr(results, "mean_transit_score") <- round(mean(results$Transit_Access_Score, na.rm = TRUE), 2)
    attr(results, "mean_transit_stops") <- round(mean(results$Transit_Stops, na.rm = TRUE), 1)
    attr(results, "mean_muni_routes")   <- round(mean(results$Unique_Muni_Routes, na.rm = TRUE), 1)
  } else {
    attr(results, "city_transit_score") <- NA_real_
    attr(results, "mean_transit_score") <- NA_real_
    attr(results, "mean_transit_stops") <- NA_real_
    attr(results, "mean_muni_routes")   <- NA_real_
  }

  if (nrow(results) > 0) {
    closest_gs <- results |>
      filter(!is.na(closest_greenspace_dist_m)) |>
      slice_min(closest_greenspace_dist_m, n = 1)
    if (nrow(closest_gs) > 0) {
      attr(results, "closest_greenspace") <- closest_gs$closest_greenspace[1]
      attr(results, "closest_greenspace_dist_m") <- closest_gs$closest_greenspace_dist_m[1]
    } else {
      attr(results, "closest_greenspace") <- "None"
      attr(results, "closest_greenspace_dist_m") <- NA_real_
    }

    closest_rsf <- results |>
      filter(!is.na(closest_rsf_program_dist_m)) |>
      slice_min(closest_rsf_program_dist_m, n = 1)
    if (nrow(closest_rsf) > 0) {
      attr(results, "closest_rsf_program") <- closest_rsf$closest_rsf_program[1]
      attr(results, "closest_rsf_program_dist_m") <- closest_rsf$closest_rsf_program_dist_m[1]
    } else {
      attr(results, "closest_rsf_program") <- "None"
      attr(results, "closest_rsf_program_dist_m") <- NA_real_
    }
  }

  results
}


# ----------------------------------------------------------------------------
# add_bai(): metric df + citywide benchmarks -> df + 7 BAI axes + composite BAI
# ----------------------------------------------------------------------------
# df: a data.frame from compute_iso_metrics(). city_benchmarks: the list of
# citywide reference distributions (the ECDFs each axis is scored against).
# Returns df with the seven *_std columns (each 0-1) plus BAI (their row mean),
# or NULL if df is empty.
add_bai <- function(df, city_benchmarks) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  ref <- city_benchmarks

  eq_ref_inverted <- if (length(ref$ej) > 0) max(ref$ej, na.rm = TRUE) - ref$ej else numeric(0)
  eq_obs_inverted <- ifelse(is.na(df$SF_EJ_Score), NA_real_, max(ref$ej, na.rm = TRUE) - df$SF_EJ_Score)

  tmp <- df |>
    mutate(
      Mobility_Access_std        = ecdf01(Transit_Access_Score, ref$transit_density),
      Biodiversity_Potential_std = ecdf01(GBIF_Species, ref$biodiversity),
      Observation_Intensity_std  = ecdf01(SamplingDensity_km2, ref$sampling),
      Environmental_Quality_std  = ecdf01(MeanNDVI, ref$ndvi),
      Greenspace_Cover_std       = ecdf01(Greenspace_percent, ref$greenspace_cover),
      Equity_Context_std         = ecdf01(eq_obs_inverted, eq_ref_inverted),
      Route_Access_std           = ecdf01(Unique_Muni_Routes, ref$route_access)
    )

  tmp$BAI <- rowMeans(
    tmp[, c(
      "Mobility_Access_std",
      "Biodiversity_Potential_std",
      "Observation_Intensity_std",
      "Environmental_Quality_std",
      "Greenspace_Cover_std",
      "Equity_Context_std",
      "Route_Access_std"
    )],
    na.rm = TRUE
  )

  tmp
}


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
# bai_df: the data.frame returned by add_bai() (one row per isochrone). title:
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
# bai_rows: a list of single-row add_bai() data.frames (one per point).
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
