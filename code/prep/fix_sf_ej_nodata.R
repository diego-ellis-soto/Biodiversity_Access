###############################################################################
# fix_sf_ej_nodata.R — remove the SF EJ "999" no-data code from the BAI inputs
#
# PROBLEM
#   sf_ej_communities_map.gpkg has one polygon (12.3 km² of parks and open
#   space, drawn grey 225,225,225 in the source map) with score = 999. That is
#   the map's no-data class; real scores run 0-30. The grid builder
#   (BAI_100m/code/city_grid_analysis.R, section 4b) area-weights `score` per
#   cell and drops only NA, so 999 leaked into grid_results_100m.gpkg:
#   ~1,600 cells are exactly 999 and ~740 more are blends of 999 with a real
#   score. Because a higher EJ score means more burden, Golden Gate Park,
#   McLaren Park, Lincoln Park, etc. became the most "burdened" places in the
#   city, and every isochrone touching them was pushed down the Equity axis
#   of the BAI.
#
# FIX (only sf_ej_score and pctile_sf_ej_score change; nothing else is touched)
#   1. Grid: cells overlapping the 999 polygon get the area-weighted mean of
#      their non-999 pieces, or NA if they have none (the grid builder's stated
#      intent: "parks, water, Presidio have no census-tract-based EJ scores").
#      All other cells keep their stored value byte-for-byte.
#   2. Walking / cycling / driving isochrones: the support cells are rebuilt
#      from the saved Mapbox polygons with the exact Step-5 rule (s2 on for
#      st_make_valid, UTM 10N centroid within polygon). The script first proves
#      it reproduces every stored n_cells_accessed and sf_ej_score, then
#      recomputes sf_ej_score from the corrected cells.
#   3. transit / walk_transit: the support cells are rebuilt with the same
#      reach query as code/transit/02_citywide_isochrones.R, the stored values
#      are reproduced first, then sf_ej_score is recomputed (see note below).
#   4. pctile_sf_ej_score is recomputed within mode x time with the exact
#      convention of code/transit/03_append_to_reference.R.
#
# NOTE ON TRANSIT
#   In 02_citywide_isochrones.R the final SELECT joins G to the species (S)
#   and route (RT) lookups before aggregating, so AVG()/SUM() over G columns
#   count each cell once per species x route row. The corrected transit
#   sf_ej_score here uses each support cell once, area-weighted, the same
#   definition as walking/cycling/driving and the app. The other fanned-out
#   transit columns are left as they are; see the console report.
#
# Run from the repository root:
#   Rscript code/prep/fix_sf_ej_nodata.R            # build + validate into outputs/results/ej_fix/
#   Rscript code/prep/fix_sf_ej_nodata.R --install  # ...then back up originals and copy fixed files in place
###############################################################################

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(data.table)
  library(FNN); library(DBI); library(duckdb)
})

INSTALL <- "--install" %in% commandArgs(trailingOnly = TRUE)
stopifnot("run from the repository root" = file.exists("code/iso_metric_definitions.R"))
source("code/iso_metric_definitions.R")

GRID_FILE <- "outputs/results/grid_results_100m.gpkg"
REF_FILE  <- "outputs/results/isochrone_basi_100m.rds"
EJ_FILE   <- "data/cached/sf_ej_communities_map.gpkg"
ISO_FILE  <- path.expand(Sys.getenv(
  "BAI_ISO_POLYGONS",
  "~/Desktop/Projects/Postdoc/BAI_100m/outputs/isochrones/sf_isochrones_100m.gpkg"
))
TRANSIT_DIR <- "outputs/results/transit"
OUT_DIR   <- "outputs/results/ej_fix"
EJ_NODATA <- 999
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(GRID_FILE, REF_FILE, EJ_FILE, ISO_FILE, file.path(TRANSIT_DIR, "stop_to_stop_15min.rds"))) {
  if (!file.exists(f)) stop("Missing input: ", f)
}
t_start <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")
same_val <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)

# =============================================================================
# 1. GRID: per-cell SF EJ score without the 999 polygon
# =============================================================================
say("[1/5] Grid cells")
grid_wgs84 <- st_read(GRID_FILE, quiet = TRUE)
ej <- st_read(EJ_FILE, quiet = TRUE)
stopifnot(sum(ej$score >= EJ_NODATA, na.rm = TRUE) == 1)

# Same construction as city_grid_analysis.R section 4b (EPSG:3857 areas)
grid_3857 <- st_transform(grid_wgs84, 3857)
ej_3857   <- st_transform(ej, 3857) |> st_make_valid()
ej_inter  <- suppressWarnings(st_intersection(ej_3857[, "score"], grid_3857[, "cell_id"]))
ej_pieces <- data.table(
  cell_id = ej_inter$cell_id,
  score   = ej_inter$score,
  area_w  = as.numeric(st_area(st_geometry(ej_inter)))
)[!is.na(score)]

cells_rebuilt <- ej_pieces[, .(
  ej_with_999 = round(weighted.mean(score, area_w), 2),
  ej_clean    = if (any(score < EJ_NODATA & area_w > 0))
                  round(weighted.mean(score[score < EJ_NODATA], area_w[score < EJ_NODATA]), 2)
                else NA_real_,
  nodata_m2   = sum(area_w[score >= EJ_NODATA])
), by = cell_id]

cell <- data.table(cell_id = grid_wgs84$cell_id, ej_old = grid_wgs84$sf_ej_score)
cell <- merge(cell, cells_rebuilt, by = "cell_id", all.x = TRUE, sort = FALSE)
cell <- cell[match(grid_wgs84$cell_id, cell_id)]

# Validation: the rebuilt-with-999 value reproduces the stored grid value
repro <- cell[is.finite(ej_old)]
err <- abs(repro$ej_with_999 - repro$ej_old)
say("  stored cell scores reproduced (|diff| <= 0.01): %d of %d, max |diff| %.4f",
    sum(err <= 0.01 + 1e-9, na.rm = TRUE), nrow(repro), max(err, na.rm = TRUE))
stopifnot("grid EJ reconstruction does not reproduce the stored values" =
            mean(err <= 0.01 + 1e-9, na.rm = TRUE) > 0.999)

affected <- cell[is.finite(nodata_m2) & nodata_m2 > 0, cell_id]
cell[, ej_new := ej_old]
cell[cell_id %in% affected, ej_new := ej_clean]
say("  cells overlapping the 999 polygon: %d -> %d become NA, %d get a re-weighted real score",
    length(affected), sum(is.na(cell[cell_id %in% affected, ej_new])),
    sum(!is.na(cell[cell_id %in% affected, ej_new])))
say("  max cell sf_ej_score: before %.2f, after %.2f",
    max(cell$ej_old, na.rm = TRUE), max(cell$ej_new, na.rm = TRUE))
stopifnot(max(cell$ej_new, na.rm = TRUE) < EJ_NODATA / 10)

# Write the corrected grid with GDAL (same layer name, column order and
# geometry). An in-place SQLite UPDATE is not possible without SpatiaLite: the
# GeoPackage's R-tree triggers call ST_IsEmpty() on every row update.
grid_out <- file.path(OUT_DIR, basename(GRID_FILE))
layer <- st_layers(GRID_FILE)$name[[1]]
grid_fixed <- grid_wgs84
grid_fixed$sf_ej_score <- cell$ej_new
st_write(grid_fixed, grid_out, layer = layer, delete_dsn = TRUE, quiet = TRUE)

# Read it back: geometry, CRS and every non-EJ column must be unchanged.
chk <- st_read(grid_out, quiet = TRUE)
other_cols <- setdiff(names(st_drop_geometry(grid_wgs84)), "sf_ej_score")
stopifnot(
  identical(names(st_drop_geometry(chk)), names(st_drop_geometry(grid_wgs84))),
  identical(chk$cell_id, grid_wgs84$cell_id),
  isTRUE(all.equal(chk$sf_ej_score, cell$ej_new)),
  identical(st_drop_geometry(chk)[, other_cols], st_drop_geometry(grid_wgs84)[, other_cols]),
  identical(st_crs(chk), st_crs(grid_wgs84)),
  identical(st_coordinates(chk), st_coordinates(grid_wgs84))
)
n_upd <- sum(!same_val(grid_wgs84$sf_ej_score, chk$sf_ej_score))
stopifnot(n_upd == sum(!same_val(cell$ej_old, cell$ej_new)))
say("  wrote %s (%d cell values changed; geometry and other columns identical)", grid_out, n_upd)

# Cell lookups in the Step-5 metric frame (s2 ON for st_make_valid, as Step 5 ran)
sf_use_s2(TRUE)
grid_proj <- st_transform(st_make_valid(grid_wgs84), ISO_METRIC_CRS)
cell[, area_m2 := as.numeric(st_area(grid_proj))]
centroids <- st_centroid(st_geometry(grid_proj))
idx_of <- setNames(seq_len(nrow(cell)), cell$cell_id)

# Weighted means over a support given as (group, cell index) pairs
support_means <- function(grp, ci) {
  d <- data.table(g = grp, a = cell$area_m2[ci], o = cell$ej_old[ci], n = cell$ej_new[ci])
  d[, .(n_cells = .N,
        ej_old = { ok <- is.finite(o) & a > 0; if (any(ok)) sum(o[ok] * a[ok]) / sum(a[ok]) else NA_real_ },
        ej_new = { ok <- is.finite(n) & a > 0; if (any(ok)) sum(n[ok] * a[ok]) / sum(a[ok]) else NA_real_ }),
    by = g]
}

# =============================================================================
# 2. WALKING / CYCLING / DRIVING: rebuild support from the Mapbox polygons
# =============================================================================
say("[2/5] Walking / cycling / driving isochrones")
ref <- as.data.table(readRDS(REF_FILE))
ref_cols <- copy(names(ref))
iso_layer <- st_layers(ISO_FILE)$name[[1]]
groups <- unique(ref[mode %in% c("walking", "cycling", "driving"), .(mode, time_min)])
setorder(groups, mode, time_min)

wcd <- vector("list", nrow(groups))
for (k in seq_len(nrow(groups))) {
  MODE <- groups$mode[k]; TT <- groups$time_min[k]; t0 <- Sys.time()
  iso <- st_read(ISO_FILE, quiet = TRUE, query = sprintf(
    "SELECT cell_id, mode, time_min, geom FROM \"%s\" WHERE mode = '%s' AND time_min = %s",
    iso_layer, MODE, TT))
  stopifnot(nrow(iso) == 24182L, !anyDuplicated(iso$cell_id))
  iso_proj <- st_transform(st_make_valid(iso), ISO_METRIC_CRS)
  chunks <- split(seq_len(nrow(iso_proj)), ceiling(seq_len(nrow(iso_proj)) / 2000))
  res <- rbindlist(lapply(chunks, function(ix) {
    h <- st_contains(st_geometry(iso_proj)[ix], centroids)
    out <- support_means(rep(ix, lengths(h)), unlist(h))
    # zero-support isochrones keep a row with NA metrics, as in Step 5
    empty <- setdiff(ix, out$g)
    if (length(empty)) out <- rbind(out, data.table(g = empty, n_cells = 0L, ej_old = NA_real_, ej_new = NA_real_))
    out
  }))
  res[, `:=`(cell_id = iso$cell_id[g], mode = MODE, time_min = TT)][, g := NULL]
  wcd[[k]] <- res
  say("  %-8s %2d min | %s | %.0f s", MODE, TT,
      format(sum(res$n_cells), big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
wcd <- rbindlist(wcd)

# =============================================================================
# 3. TRANSIT / WALK_TRANSIT: rebuild support with the stage-2 reach query
# =============================================================================
say("[3/5] Transit isochrones (code/transit/02_citywide_isochrones.R logic)")
WALK <- 60
TIMES <- c(5, 10, 15)
MODES <- list(
  transit      = list(charge_wait = FALSE, access_min = 10, egress_cap = Inf, origin_cap = Inf),
  walk_transit = list(charge_wait = TRUE,  access_min = 5,  egress_cap = 5,   origin_cap = 5))
MAX_ACCESS_MIN <- max(vapply(MODES, function(m) m$access_min, numeric(1)))

XY  <- st_coordinates(st_centroid(st_geometry(st_transform(grid_wgs84, ISO_METRIC_CRS))))
cid <- grid_wgs84$cell_id
tt  <- readRDS("data/cached/gtfs_timetable_monday_regional_transfers.rds")
stops_sf <- st_as_sf(tt$stops, coords = c("stop_lon", "stop_lat"), crs = 4326) |> st_transform(ISO_METRIC_CRS)
SXY <- st_coordinates(stops_sf); sid <- as.character(tt$stops$stop_id)
B   <- as.data.table(readRDS(file.path(TRANSIT_DIR, "stop_to_stop_15min.rds")))
hw  <- fread("data/cached/gtfs_stop_headways_regional.csv")[, .(stop_id = as.character(stop_id), mean_headway_min)]
wait_half <- hw$mean_headway_min[match(sid, hw$stop_id)] / 2
wait_half[is.na(wait_half)] <- median(hw$mean_headway_min, na.rm = TRUE) / 2

rad <- function(q, t, k, r, lab) {
  k <- min(k, nrow(t)); nn <- FNN::get.knnx(t, q, k = k)
  keep <- nn$nn.dist <= r
  if (any(keep[, k])) stop(sprintf("k=%d too small for %s", k, lab))
  i <- which(keep, arr.ind = TRUE)
  data.table(q = i[, 1], t = nn$nn.index[keep], d = nn$nn.dist[keep])
}
MAXR <- max(TIMES) * WALK
C <- rad(SXY, XY, 900, MAXR, "stop->cell")[, .(stop = sid[q], cell_id = cid[t], dist_m = d)]
A <- rad(XY, SXY, 200, MAX_ACCESS_MIN * WALK, "origin->stop")[, .(origin = cid[q], stop = sid[t], access_min = d / WALK)]
A[, wait_min := wait_half[match(stop, sid)]]
W <- rad(XY, XY, 900, MAXR, "origin->cell")[, .(origin = cid[q], cell_id = cid[t], dist_m = d)]

# Fan-out weights used by the stage-2 SELECT (rows per cell after LEFT JOIN S, RT)
spp <- unique(as.data.table(readRDS("outputs/results/gbif_cell_species_100m.rds"))[!is.na(species), .(cell_id, species, class)])
RT  <- unique(as.data.table(readRDS("outputs/results/muni_route_cells_100m.rds"))[, .(cell_id, route_id)])
fan <- data.table(cell_id = cid)
fan <- merge(fan, spp[, .(ns = .N), by = cell_id], by = "cell_id", all.x = TRUE)
fan <- merge(fan, RT[, .(nr = .N), by = cell_id], by = "cell_id", all.x = TRUE)
fan[is.na(ns), ns := 0L][is.na(nr), nr := 0L][, w := pmax(1L, ns) * pmax(1L, nr)]
cell[, fan_w := fan$w[match(cell$cell_id, fan$cell_id)]]
cell[, n_records := grid_wgs84$n_records]

con <- dbConnect(duckdb::duckdb()); invisible(dbExecute(con, "PRAGMA threads=6"))
duckdb_register(con, "A", A); duckdb_register(con, "B", B)
duckdb_register(con, "C", C); duckdb_register(con, "W", W)
origins <- sort(unique(W$origin)); ochunks <- split(origins, ceiling(seq_along(origins) / 2000))
BIG <- 1e6
tr <- list()
for (MD in names(MODES)) {
  cfg <- MODES[[MD]]
  wait_expr  <- if (cfg$charge_wait) "A.wait_min" else "0"
  egress_cap <- if (is.finite(cfg$egress_cap)) cfg$egress_cap else BIG
  origin_cap <- if (is.finite(cfg$origin_cap)) cfg$origin_cap else BIG
  for (TT in TIMES) {
    t0 <- Sys.time()
    acc <- lapply(ochunks, function(ch) {
      dbExecute(con, sprintf("CREATE OR REPLACE TEMP TABLE ch AS SELECT UNNEST([%s]) AS origin",
                             paste(ch, collapse = ",")))
      as.data.table(dbGetQuery(con, sprintf("
        WITH reach AS (
          SELECT A.origin, B.to_stop, MIN(A.access_min + %s + B.ride_min) AS cost
          FROM A JOIN ch USING (origin) JOIN B ON A.stop = B.from_stop
          WHERE A.access_min <= %f
          GROUP BY 1,2 HAVING MIN(A.access_min + %s + B.ride_min) <= %f),
        cells AS (
          SELECT DISTINCT r.origin, C.cell_id FROM reach r JOIN C ON C.stop = r.to_stop
            WHERE C.dist_m <= LEAST(%f - r.cost, %f) * %f
          UNION
          SELECT W.origin, W.cell_id FROM W JOIN ch USING (origin) WHERE W.dist_m <= LEAST(%f, %f) * %f)
        SELECT origin, cell_id FROM cells",
        wait_expr, cfg$access_min, wait_expr, TT, TT, egress_cap, WALK, TT, origin_cap, WALK)))
    })
    pairs <- rbindlist(acc)
    ci <- idx_of[as.character(pairs$cell_id)]
    d <- data.table(origin = pairs$origin, a = cell$area_m2[ci], o = cell$ej_old[ci],
                    n = cell$ej_new[ci], w = cell$fan_w[ci], r = cell$n_records[ci])
    res <- d[, .(
      n_cells        = .N,
      ej_old_fanout  = { ok <- is.finite(o); if (any(ok)) sum(o[ok] * w[ok]) / sum(w[ok]) else NA_real_ },
      ej_new         = { ok <- is.finite(n) & a > 0; if (any(ok)) sum(n[ok] * a[ok]) / sum(a[ok]) else NA_real_ },
      records_fanout = sum(r * w), records_distinct = sum(r)
    ), by = origin]
    setnames(res, "origin", "cell_id")
    res[, `:=`(mode = MD, time_min = TT)]
    tr[[paste(MD, TT)]] <- res
    say("  %-12s %2d min | %s pairs | %.0f s", MD, TT, format(nrow(pairs), big.mark = ","),
        as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
}
dbDisconnect(con, shutdown = TRUE)
tr <- rbindlist(tr)

# =============================================================================
# 4. VALIDATE reproduction of the stored reference, then rebuild EJ columns
# =============================================================================
say("[4/5] Validation and percentile rebuild")
key <- c("cell_id", "mode", "time_min")
v_wcd <- merge(ref[, c(key, "n_cells_accessed", "sf_ej_score"), with = FALSE], wcd, by = key)
v_tr  <- merge(ref[, c(key, "n_cells_accessed", "sf_ej_score", "n_records_accessed"), with = FALSE], tr, by = key)
stopifnot(nrow(v_wcd) == nrow(wcd), nrow(v_tr) == nrow(tr),
          nrow(v_wcd) + nrow(v_tr) == nrow(ref))

report <- function(lbl, ok, n) say("  %-58s %s of %s", lbl, format(sum(ok), big.mark = ","), format(n, big.mark = ","))
same_num <- function(a, b, tol) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
ok1 <- v_wcd$n_cells_accessed == v_wcd$n_cells
ok2 <- same_num(v_wcd$sf_ej_score, v_wcd$ej_old, 1e-6)
ok3 <- v_tr$n_cells_accessed == v_tr$n_cells
ok4 <- same_num(v_tr$sf_ej_score, v_tr$ej_old_fanout, 1e-6)
report("walk/cycle/drive: n_cells_accessed reproduced", ok1, nrow(v_wcd))
report("walk/cycle/drive: stored sf_ej_score reproduced (1e-6)", ok2, nrow(v_wcd))
report("transit: n_cells_accessed reproduced", ok3, nrow(v_tr))
report("transit: stored sf_ej_score reproduced (incl. fan-out)", ok4, nrow(v_tr))
stopifnot("support reconstruction is not exact" = mean(ok1) > 0.9999 && mean(ok2) > 0.9999 &&
            mean(ok3) > 0.9999 && mean(ok4) > 0.9999)

fan_ratio <- v_tr[records_distinct > 0, n_records_accessed / records_distinct]
say("  transit n_records_accessed / distinct-cell records: median %.1f, p90 %.1f (1.0 = no fan-out)",
    median(fan_ratio), quantile(fan_ratio, 0.9))

new_ej <- rbind(wcd[, c(key, "ej_new"), with = FALSE], tr[, c(key, "ej_new"), with = FALSE])
fixed <- merge(ref, new_ej, by = key, sort = FALSE)
stopifnot(nrow(fixed) == nrow(ref))
fixed <- fixed[match(paste(ref$cell_id, ref$mode, ref$time_min), paste(cell_id, mode, time_min))]
fixed[, sf_ej_score := ej_new][, ej_new := NULL]

# Percentile convention of code/transit/03_append_to_reference.R
pctile <- function(x) {
  ok <- is.finite(x); o <- rep(NA_real_, length(x)); n <- sum(ok)
  if (n < 2) return(o)
  o[ok] <- 100 * (data.table::frank(x[ok], ties.method = "min") - 1) / (n - 1); o
}
conv <- ref[, .(stored = pctile_sf_ej_score, rebuilt = 100 - pctile(sf_ej_score)), by = .(mode, time_min)]
say("  stored pctile_sf_ej_score reproduced from the old raw values: max |diff| %.2e",
    max(abs(conv$stored - conv$rebuilt), na.rm = TRUE))
stopifnot(max(abs(conv$stored - conv$rebuilt), na.rm = TRUE) < 1e-6)
fixed[, pctile_sf_ej_score := 100 - pctile(sf_ej_score), by = .(mode, time_min)]
setcolorder(fixed, ref_cols)
attributes(fixed)[setdiff(names(attributes(ref)), c("names", "row.names", "class", ".internal.selfref", "index"))] <-
  attributes(ref)[setdiff(names(attributes(ref)), c("names", "row.names", "class", ".internal.selfref", "index"))]

untouched <- setdiff(ref_cols, c("sf_ej_score", "pctile_sf_ej_score"))
stopifnot("a column other than the EJ pair changed" =
            isTRUE(all.equal(as.data.frame(fixed[, ..untouched]), as.data.frame(ref[, ..untouched]), check.attributes = FALSE)))

# BAI before/after summary (7-axis mean of stored percentiles, as in the app)
bai_cols <- c("pctile_transit_access_score", "pctile_unique_muni_routes", "pctile_n_species_accessed",
              "pctile_sampling_density_km2", "pctile_mean_ndvi", "pctile_greenspace_pct", "pctile_sf_ej_score")
bai <- function(d) { m <- as.matrix(d[, ..bai_cols]); ifelse(rowSums(is.finite(m)) == 7, rowMeans(m) / 100, NA_real_) }
delta <- data.table(ref[, ..key], ej_old = ref$sf_ej_score, ej_new = fixed$sf_ej_score,
                    bai_old = bai(ref), bai_new = bai(fixed))
delta[, d_bai := bai_new - bai_old]
summ <- delta[, .(max_raw_ej_before = round(max(ej_old, na.rm = TRUE), 1),
                  max_raw_ej_after  = round(max(ej_new, na.rm = TRUE), 1),
                  pct_rows_bai_moved_gt_0.01 = round(100 * mean(abs(d_bai) > 0.01, na.rm = TRUE), 1),
                  max_bai_gain = round(max(d_bai, na.rm = TRUE), 3),
                  max_bai_drop = round(min(d_bai, na.rm = TRUE), 3),
                  rows_with_bai_before = sum(is.finite(bai_old)),
                  rows_with_bai_after  = sum(is.finite(bai_new))),
              by = .(mode, time_min)]
setorder(summ, mode, time_min)
print(summ)
fwrite(summ, file.path(OUT_DIR, "ej_fix_bai_change_by_mode_time.csv"))
fwrite(delta, file.path(OUT_DIR, "ej_fix_row_changes.csv.gz"))

ref_out <- file.path(OUT_DIR, basename(REF_FILE))
saveRDS(fixed, ref_out)
say("  wrote %s", ref_out)

# =============================================================================
# 5. INSTALL (optional): back up the originals, copy the fixed files in place
# =============================================================================
if (INSTALL) {
  say("[5/5] Installing")
  bk <- file.path(OUT_DIR, "pre_fix_backup"); dir.create(bk, showWarnings = FALSE)
  for (f in c(GRID_FILE, REF_FILE)) {
    b <- file.path(bk, basename(f))
    if (!file.exists(b)) stopifnot(file.copy(f, b))
    stopifnot(file.copy(file.path(OUT_DIR, basename(f)), f, overwrite = TRUE))
    say("  %s replaced (original kept at %s)", f, b)
  }
} else {
  say("[5/5] Not installed. Rerun with --install to copy the fixed files into outputs/results/.")
}
say("Done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
