# =============================================================================
# Shared raw-metric definitions for SF Biodiversity Access isochrone scoring
# =============================================================================
# This file defines the scientific support used by BOTH:
#   * Rscripts/iso_metrics.R (interactive user-generated isochrones)
#   * 03_compute_isochrone_basi.R (citywide precomputed reference isochrones)
#
# NON-NEGOTIABLE SUPPORT DEFINITION
# ---------------------------------
# For each isochrone independently, the metric support is the set of SF 100 m
# grid cells whose projected grid-cell centroids fall within that single
# isochrone polygon. Metrics are then calculated from those selected cells.
# This deliberately makes the interactive and precomputed estimands identical.
#
# It is a 100 m grid-supported access statistic, not a claim of pixel-perfect
# clipping at the isochrone boundary.
# =============================================================================

ISO_METRIC_DEFINITION_VERSION <- "sf_bai_grid100m_centroid_support_v1"
# UTM Zone 10N is appropriate for San Francisco metric area/distance calculations.
# Do not use Web Mercator (EPSG:3857) for area-based density denominators.
ISO_METRIC_CRS <- 26910

ISO_PLANT_CLASSES <- c(
  "Magnoliopsida", "Liliopsida", "Pinopsida", "Polypodiopsida",
  "Equisetopsida", "Bryopsida", "Marchantiopsida"
)

iso_weighted_mean <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

iso_safe_sum <- function(x, na_as_zero = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  if (na_as_zero) x[!is.finite(x)] <- 0
  if (!any(is.finite(x))) return(NA_real_)
  sum(x[is.finite(x)])
}

iso_safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

# Clean, explicit reference schema. `direction` describes whether a larger raw
# value is better for the access score. Step 5 stores raw values plus pctile_*
# columns; the Shiny app uses the RAW columns for external-observation scoring.
ISO_REFERENCE_METRICS <- list(
  n_species_accessed      = list(direction = "higher", label = "Species richness"),
  n_birds_accessed        = list(direction = "higher", label = "Bird richness"),
  n_mammals_accessed      = list(direction = "higher", label = "Mammal richness"),
  n_plants_accessed       = list(direction = "higher", label = "Plant richness"),
  n_records_accessed      = list(direction = "higher", label = "GBIF records"),
  sampling_density_km2    = list(direction = "higher", label = "GBIF records / km2"),
  mean_ndvi               = list(direction = "higher", label = "Mean NDVI"),
  greenspace_pct          = list(direction = "higher", label = "Greenspace cover"),
  n_stops_accessed        = list(direction = "higher", label = "Transit stops"),
  transit_access_score    = list(direction = "higher", label = "Transit stops / km2"),
  unique_muni_routes      = list(direction = "higher", label = "Unique Muni routes"),
  nearest_stop_m          = list(direction = "lower",  label = "Nearest stop distance"),
  calenviro_ci_score      = list(direction = "lower",  label = "CalEnviroScreen burden"),
  calenviro_traffic_pct   = list(direction = "lower",  label = "Traffic burden percentile"),
  sf_ej_score             = list(direction = "lower",  label = "SF EJ burden")
)
