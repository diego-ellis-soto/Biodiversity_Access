# =============================================================================
# SF Find Neighborhoods — coarse-scale summaries of the SF Biodiversity Access
# 100 m products
# =============================================================================
# Usage (from anywhere):
#   Rscript summarise_sf_find_neighborhoods.R <space_dir> <neighborhoods.geojson> \
#           <hf_data_dir> <out_dir>
#
#   space_dir     this repository (or a clone of the HF Space
#                 diegoellissoto/Biodiversity_Access); needs
#                 outputs/results/*_100m.* and code/iso_metric*.R
#   neighborhoods SF Find Neighborhoods polygons (data.sf.gov dataset gfpk-269f,
#                 the table behind the pty2-tcw4 map), GeoJSON:
#                 https://data.sf.gov/resource/gfpk-269f.geojson?$limit=5000
#   hf_data_dir   folder holding sf_ej_communities_map.gpkg and
#                 plant_provenance.csv (data/cached, or the HF dataset
#                 boettiger-lab/sf_biodiv_access)
#   out_dir       where the CSV / GeoPackage outputs are written
#
# From the repository root:
#   Rscript code/prep/summarise_sf_find_neighborhoods.R . \
#     data/cached/sf_find_neighborhoods.geojson data/cached \
#     outputs/tables/sf_find_neighborhoods
#
# Two families of summaries are produced:
#
# A. "What is inside the neighborhood" (neighborhood_area_summary.csv)
#    The neighborhood polygon is passed to the app's own
#    summarise_iso_metric_support() from code/iso_metrics_AUDITED.R, so the
#    support is exactly the app's: 100 m cells whose UTM-10N centroids fall in
#    the polygon. Species/birds/mammals/plants are true unions of species
#    (not sums of cell counts); routes are unique Muni route_ids.
#    Socio-economic columns are population-weighted means of the cell values.
#
# B. "What residents can reach from the neighborhood"
#    (neighborhood_bai_by_mode_time.csv)
#    Every 100 m origin cell already has a precomputed isochrone for each
#    mode x time in isochrone_basi_100m.rds. Its 7-axis BAI is computed exactly
#    like the app's get_reference_bai_distribution(): mean of the seven stored
#    same-mode x same-time percentiles / 100, complete rows only. Those origin
#    BAIs are then summarised over the origins inside each neighborhood.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript summarise_sf_find_neighborhoods.R <space_dir> <neighborhoods.geojson> <hf_data_dir> <out_dir>")
}
space_dir <- normalizePath(args[[1]], mustWork = TRUE)
nbhd_path <- normalizePath(args[[2]], mustWork = TRUE)
hf_dir    <- normalizePath(args[[3]], mustWork = TRUE)
out_dir   <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir   <- normalizePath(out_dir, mustWork = TRUE)

sf_use_s2(FALSE)

# -----------------------------------------------------------------------------
# 0. Load the app's own metric code (definitions + audited helpers)
# -----------------------------------------------------------------------------
# The audited helpers read outputs/results/* through relative paths, exactly
# as they do inside the app, so run them from the Space root.
old_wd <- setwd(space_dir)
on.exit(setwd(old_wd), add = TRUE)
Sys.setenv(ISO_METRIC_GRID_GPKG = file.path(space_dir, "outputs", "results", "grid_results_100m.gpkg"))

app_env <- new.env()
sys.source(file.path("code", "iso_metric_definitions.R"), envir = app_env)
sys.source(file.path("code", "iso_metrics_AUDITED.R"), envir = app_env)

# % native plants uses the same lookup the app builds in setup_unified.R
prov <- read.csv(file.path(hf_dir, "plant_provenance.csv"), colClasses = "character")
app_env$plant_establishment <- setNames(prov$establishment_means, prov$species)

CRS_M <- app_env$ISO_METRIC_CRS          # 26910
stopifnot(identical(app_env$ISO_METRIC_DEFINITION_VERSION, "sf_bai_grid100m_centroid_support_v1"))

# -----------------------------------------------------------------------------
# 1. Neighborhoods and the cell -> neighborhood lookup (centroid-in-polygon)
# -----------------------------------------------------------------------------
nbhd <- st_read(nbhd_path, quiet = TRUE) |>
  st_make_valid() |>
  select(neighborhood = name) |>
  arrange(neighborhood)
stopifnot(!anyDuplicated(nbhd$neighborhood))
nbhd_m <- st_transform(nbhd, CRS_M)
nbhd_m$polygon_area_km2 <- as.numeric(st_area(nbhd_m)) / 1e6

grid_cache <- app_env$load_iso_metric_grid()     # same object the app uses
cent <- grid_cache$centroids_metric
hits <- st_within(cent, nbhd_m)
cell_nbhd <- data.table(
  cell_id = rep(cent$cell_id, lengths(hits)),
  neighborhood = nbhd_m$neighborhood[unlist(hits)]
)
n_multi <- sum(lengths(hits) > 1)
n_none  <- sum(lengths(hits) == 0)
message(sprintf("Cells: %d total | %d in exactly one neighborhood | %d in >1 | %d in none",
                nrow(cent), sum(lengths(hits) == 1), n_multi, n_none))

# socio-economic cell attributes (not part of the audited metric cache)
grid_socio <- st_read(Sys.getenv("ISO_METRIC_GRID_GPKG"), quiet = TRUE) |>
  st_drop_geometry() |>
  select(cell_id, est_population, med_income_k, pct_nonwhite, pct_snap) |>
  as.data.table()

popw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# -----------------------------------------------------------------------------
# 2. SF EJ score from the source polygons, excluding the 999 no-data code
# -----------------------------------------------------------------------------
# The SF EJ Communities layer contains one polygon whose score is 999 (no-data
# code; real scores run 0-30). Before code/prep/fix_sf_ej_nodata.R,
# grid_results_100m.gpkg area-averaged it into cell values, inflating sf_ej_app
# wherever that polygon touches. sf_ej_clean recomputes an area-weighted score
# straight from the source polygons with the 999 polygon dropped; on the fixed
# grid the two should agree closely (they differ only by cell-centroid support
# vs exact polygon overlap).
ej <- st_read(file.path(hf_dir, "sf_ej_communities_map.gpkg"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(CRS_M)
ej_valid <- ej[is.finite(ej$score) & ej$score < 999, "score"]
ej_nodata <- ej[is.finite(ej$score) & ej$score >= 999, ]

ej_int <- suppressWarnings(st_intersection(nbhd_m[, "neighborhood"], ej_valid))
ej_int$a <- as.numeric(st_area(ej_int))
ej_clean <- ej_int |>
  st_drop_geometry() |>
  group_by(neighborhood) |>
  summarise(sf_ej_clean = sum(score * a) / sum(a),
            sf_ej_scored_km2 = sum(a) / 1e6, .groups = "drop")

ej_nd <- suppressWarnings(st_intersection(nbhd_m[, "neighborhood"], st_union(ej_nodata)))
ej_nd <- if (nrow(ej_nd) > 0) {
  data.frame(neighborhood = ej_nd$neighborhood, sf_ej_nodata_km2 = as.numeric(st_area(ej_nd)) / 1e6)
} else data.frame(neighborhood = character(0), sf_ej_nodata_km2 = numeric(0))

# -----------------------------------------------------------------------------
# 3. A. Inside-the-neighborhood summary (app's summarise_iso_metric_support)
# -----------------------------------------------------------------------------
area_rows <- lapply(seq_len(nrow(nbhd_m)), function(i) {
  s <- app_env$summarise_iso_metric_support(nbhd_m[i, ])
  soc <- grid_socio[cell_id %in% s$cell_ids]
  data.frame(
    neighborhood        = nbhd_m$neighborhood[i],
    polygon_area_km2    = nbhd_m$polygon_area_km2[i],
    n_cells             = s$n_cells,
    support_area_km2    = s$area_km2,
    gbif_records        = s$n_records,
    gbif_records_km2    = s$sampling_density,
    species             = s$n_species,
    bird_species        = s$n_birds,
    mammal_species      = s$n_mammals,
    plant_species       = s$n_plants,
    plant_native        = s$n_plants_native,
    plant_introduced    = s$n_plants_introduced,
    pct_plants_native   = ifelse((s$n_plants_native + s$n_plants_introduced) > 0,
                                 100 * s$n_plants_native / (s$n_plants_native + s$n_plants_introduced),
                                 NA_real_),
    mean_ndvi           = s$mean_ndvi,
    greenspace_pct      = s$greenspace_pct,
    transit_stops       = s$n_stops,
    transit_stops_km2   = s$transit_access,
    muni_routes         = s$unique_routes,
    calenviro_ci        = s$calenviro_ci,
    calenviro_traffic_pctl = s$traffic_pct,
    sf_ej_app           = s$sf_ej,
    est_population      = sum(soc$est_population, na.rm = TRUE),
    med_income_k_popw   = popw(soc$med_income_k, soc$est_population),
    pct_nonwhite_popw   = 100 * popw(soc$pct_nonwhite, soc$est_population),
    pct_snap_popw       = 100 * popw(soc$pct_snap, soc$est_population)
  )
})
area_summary <- bind_rows(area_rows) |>
  left_join(ej_clean, by = "neighborhood") |>
  left_join(ej_nd, by = "neighborhood") |>
  mutate(sf_ej_nodata_km2 = replace_na(sf_ej_nodata_km2, 0)) |>
  relocate(sf_ej_clean, sf_ej_scored_km2, sf_ej_nodata_km2, .after = sf_ej_app)

# -----------------------------------------------------------------------------
# 4. B. Reachable-from-the-neighborhood BAI, by mode x travel time
# -----------------------------------------------------------------------------
ref <- as.data.table(readRDS(file.path("outputs", "results", "isochrone_basi_100m.rds")))
# Transit rows were appended later (code/transit/03_append_to_reference.R) and
# carry NA in reference_schema_version; every non-NA value must match.
stopifnot(all(na.omit(ref$reference_schema_version) == app_env$ISO_METRIC_DEFINITION_VERSION))

# Same column -> axis mapping as add_reference_bai() / .reference_bai_pct_cols
axis_map <- c(
  Mobility_Access_std        = "pctile_transit_access_score",
  Route_Access_std           = "pctile_unique_muni_routes",
  Biodiversity_Potential_std = "pctile_n_species_accessed",
  Observation_Intensity_std  = "pctile_sampling_density_km2",
  Environmental_Quality_std  = "pctile_mean_ndvi",
  Greenspace_Cover_std       = "pctile_greenspace_pct",
  Equity_Context_std         = "pctile_sf_ej_score"
)
for (ax in names(axis_map)) ref[, (ax) := get(axis_map[[ax]]) / 100]
axis_mat <- as.matrix(ref[, names(axis_map), with = FALSE])
complete <- rowSums(is.finite(axis_mat)) == length(axis_map)
ref[, BAI := NA_real_]
ref[complete, BAI := rowMeans(axis_mat[complete, , drop = FALSE])]

# Citywide quartile thresholds of origin BAI within each mode x time
ref[, bai_city_p75 := quantile(BAI, 0.75, na.rm = TRUE), by = .(mode, time_min)]

ref_n <- merge(ref, cell_nbhd, by = "cell_id", allow.cartesian = TRUE)
ref_n <- merge(ref_n, grid_socio[, .(cell_id, est_population)], by = "cell_id", all.x = TRUE)

bai_long <- ref_n[, c(
  list(
    n_origins          = .N,
    n_origins_with_bai = sum(is.finite(BAI)),
    bai_mean           = mean(BAI, na.rm = TRUE),
    bai_median         = median(BAI, na.rm = TRUE),
    bai_p25            = as.numeric(quantile(BAI, 0.25, na.rm = TRUE, names = FALSE)),
    bai_p75            = as.numeric(quantile(BAI, 0.75, na.rm = TRUE, names = FALSE)),
    bai_popw           = popw(BAI, est_population),
    pct_origins_city_top_quartile = 100 * mean(BAI[is.finite(BAI)] >= bai_city_p75[1]),
    med_species_reachable     = median(n_species_accessed, na.rm = TRUE),
    med_greenspace_pct_reach  = median(greenspace_pct, na.rm = TRUE),
    med_ndvi_reach            = median(mean_ndvi, na.rm = TRUE),
    med_stops_km2_reach       = median(transit_access_score, na.rm = TRUE),
    med_routes_reach          = median(unique_muni_routes, na.rm = TRUE),
    med_records_km2_reach     = median(sampling_density_km2, na.rm = TRUE)
  ),
  lapply(.SD, function(v) mean(v, na.rm = TRUE))
), by = .(neighborhood, mode, time_min), .SDcols = names(axis_map)]

bai_long[!is.finite(bai_mean), bai_mean := NA_real_]
bai_long[, bai_rank := frank(-bai_mean, ties.method = "min", na.last = "keep"),
         by = .(mode, time_min)]
setorder(bai_long, mode, time_min, bai_rank)
setcolorder(bai_long, c("neighborhood", "mode", "time_min", "bai_rank"))

# -----------------------------------------------------------------------------
# 5. Write outputs
# -----------------------------------------------------------------------------
round_df <- function(d, digits = 4) {
  d <- as.data.frame(d)
  num <- vapply(d, is.double, logical(1))
  d[num] <- lapply(d[num], round, digits)
  d
}

write.csv(round_df(area_summary, 3), file.path(out_dir, "neighborhood_area_summary.csv"),
          row.names = FALSE, na = "")
write.csv(round_df(bai_long, 4), file.path(out_dir, "neighborhood_bai_by_mode_time.csv"),
          row.names = FALSE, na = "")

# GeoPackage: polygons + area summary + mean BAI for each mode x time (wide)
bai_wide <- dcast(bai_long, neighborhood ~ paste0("bai_", mode, "_", time_min, "min"),
                  value.var = "bai_mean")
gpkg <- nbhd |>
  left_join(round_df(area_summary, 3), by = "neighborhood") |>
  left_join(round_df(bai_wide, 4), by = "neighborhood")
gpkg_path <- file.path(out_dir, "sf_find_neighborhoods_summary.gpkg")
if (file.exists(gpkg_path)) file.remove(gpkg_path)
st_write(gpkg, gpkg_path, layer = "sf_find_neighborhoods_summary", quiet = TRUE)

# -----------------------------------------------------------------------------
# 6. Checks printed to the console
# -----------------------------------------------------------------------------
cat("\n== Checks ==\n")
cat(sprintf("Neighborhoods: %d | with >=1 cell: %d\n", nrow(area_summary), sum(area_summary$n_cells > 0)))
cat(sprintf("Cells assigned to >1 neighborhood: %d | to none: %d\n", n_multi, n_none))
tot_rec <- sum(grid_cache$attrs$n_records)
in_rec  <- sum(grid_cache$attrs$n_records[grid_cache$attrs$cell_id %in% cell_nbhd$cell_id])
cat(sprintf("GBIF records in cells covered by a neighborhood: %.1f%% of grid total\n", 100 * in_rec / tot_rec))
cat(sprintf("Origin BAI rows complete: %d of %d\n", sum(complete), nrow(ref)))
cat(sprintf("Neighborhoods touching the SF EJ 999 polygon: %d\n", sum(area_summary$sf_ej_nodata_km2 > 0.001)))
cat("Wrote:\n ", paste(list.files(out_dir, full.names = FALSE), collapse = "\n  "), "\n")
