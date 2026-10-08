# Stage 2: a transit isochrone for every origin x time budget, in the same
# schema as the walking/cycling/driving reference so the four modes can be
# percentiled on equal footing.
#
#   reachable = walkshed(origin, T)
#               U { walkshed(stop s, T - access(s0) - wait(s0) - ride(s0 -> s)) }
#
# Walk speed 60 m/min, CALIBRATED against the Mapbox walking isochrones already
# in the reference (their median equal-area radius at 5/10/15 min implies
# 62.3/60.1/59.4), so transit is comparable to walking rather than inflated.
suppressMessages({library(sf); library(terra); library(data.table); library(FNN)
                  library(duckdb); library(DBI)})
setwd("~/Downloads/Biodiversity_Access_Score")
source("code/iso_metric_definitions.R")
OUT <- "outputs/results/transit"

WALK <- 60; MAX_ACCESS_MIN <- 10; CHARGE_WAIT <- TRUE   # half-headway, per the paper
TIMES <- c(5, 10, 15)

grid <- st_read("outputs/results/grid_results_100m.gpkg", quiet=TRUE)
gp   <- st_transform(grid, ISO_METRIC_CRS)
cid  <- grid$cell_id
XY   <- st_coordinates(st_centroid(st_geometry(gp)))

tt <- readRDS("data/cached/gtfs_timetable_monday_transfers.rds")
stops_sf <- st_as_sf(tt$stops, coords=c("stop_lon","stop_lat"), crs=4326) |>
            st_transform(ISO_METRIC_CRS)
SXY <- st_coordinates(stops_sf); sid <- as.character(tt$stops$stop_id)
B <- as.data.table(readRDS(file.path(OUT,"stop_to_stop_15min.rds")))

hw <- fread("data/cached/gtfs_stop_headways.csv")[, .(stop_id=as.character(stop_id), mean_headway_min)]
wait <- hw$mean_headway_min[match(sid, hw$stop_id)] / 2
wait[is.na(wait)] <- median(hw$mean_headway_min, na.rm=TRUE)/2
if (!CHARGE_WAIT) wait[] <- 0
cat(sprintf("stops %d | grid %d | median wait %.1f min\n", length(sid), nrow(XY), median(wait)))

# radius queries; k must be large enough that no query saturates
rad <- function(q, t, k, r, lab) {
  k <- min(k, nrow(t)); nn <- FNN::get.knnx(t, q, k=k)
  keep <- nn$nn.dist <= r
  if (any(keep[,k])) stop(sprintf("k=%d too small for %s", k, lab))
  i <- which(keep, arr.ind=TRUE)
  data.table(q=i[,1], t=nn$nn.index[keep], d=nn$nn.dist[keep])
}
MAXR <- max(TIMES) * WALK
C <- rad(SXY, XY, 900, MAXR, "stop->cell")[, .(stop=sid[q], cell_id=cid[t], dist_m=d)]
A <- rad(XY, SXY, 200, MAX_ACCESS_MIN*WALK, "origin->stop")[
       , .(origin=cid[q], stop=sid[t], access_min=d/WALK)]
A[, wait_min := wait[match(stop, sid)]]
W <- rad(XY, XY, 900, MAXR, "origin->cell")[, .(origin=cid[q], cell_id=cid[t], dist_m=d)]
cat(sprintf("lookups: stop->cell %s | origin->stop %s | origin->cell %s\n",
    format(nrow(C),big.mark=","), format(nrow(A),big.mark=","), format(nrow(W),big.mark=",")))

spp <- unique(as.data.table(readRDS("outputs/results/gbif_cell_species_100m.rds"))[
              !is.na(species), .(cell_id, species, class)])
g0 <- as.data.table(st_drop_geometry(grid))[, .(cell_id, n_records, mean_ndvi, gs_area_m2,
        n_stops, nearest_stop_m, calenviro_ci_score, calenviro_traffic_pct, sf_ej_score)]
area_m2 <- as.numeric(st_area(gp)); names(area_m2) <- as.character(cid)

con <- dbConnect(duckdb::duckdb()); dbExecute(con,"PRAGMA memory_limit='6GB'"); dbExecute(con,"PRAGMA threads=6")
duckdb_register(con,"A",A); duckdb_register(con,"B",B); duckdb_register(con,"C",C)
duckdb_register(con,"W",W); duckdb_register(con,"S",spp); duckdb_register(con,"G",g0)
# Muni routes touching each cell: needed for the unique_muni_routes axis, which
# is one of the seven IB-BAI radar spokes. Without it the transit radar has holes.
RT <- unique(as.data.table(readRDS("outputs/results/muni_route_cells_100m.rds"))[, .(cell_id, route_id)])
duckdb_register(con,"RT",RT)

origins <- sort(unique(W$origin)); chunks <- split(origins, ceiling(seq_along(origins)/2000))
PLANT <- c("Magnoliopsida","Liliopsida","Pinopsida","Polypodiopsida","Bryopsida","Lycopodiopsida")
out <- list(); t0 <- Sys.time()
for (TT in TIMES) {
  acc <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    dbExecute(con, sprintf("CREATE OR REPLACE TEMP TABLE ch AS SELECT UNNEST([%s]) AS origin",
                           paste(chunks[[i]], collapse=",")))
    acc[[i]] <- dbGetQuery(con, sprintf("
      WITH reach AS (
        SELECT A.origin, B.to_stop, MIN(A.access_min + A.wait_min + B.ride_min) AS cost
        FROM A JOIN ch USING (origin) JOIN B ON A.stop = B.from_stop
        GROUP BY 1,2 HAVING MIN(A.access_min + A.wait_min + B.ride_min) <= %f),
      cells AS (
        SELECT DISTINCT r.origin, C.cell_id FROM reach r JOIN C ON C.stop = r.to_stop
          WHERE C.dist_m <= (%f - r.cost) * %f
        UNION
        SELECT W.origin, W.cell_id FROM W JOIN ch USING (origin) WHERE W.dist_m <= %f * %f)
      SELECT c.origin AS cell_id,
             COUNT(DISTINCT c.cell_id)                         AS n_cells_accessed,
             SUM(G.n_records)                                  AS n_records_accessed,
             AVG(G.mean_ndvi)                                  AS mean_ndvi,
             SUM(G.gs_area_m2)                                 AS greenspace_area_m2,
             SUM(G.n_stops)                                    AS n_stops_accessed,
             MIN(G.nearest_stop_m)                             AS nearest_stop_m,
             AVG(G.calenviro_ci_score)                         AS calenviro_ci_score,
             AVG(G.calenviro_traffic_pct)                      AS calenviro_traffic_pct,
             AVG(G.sf_ej_score)                                AS sf_ej_score,
             COUNT(DISTINCT S.species)                         AS n_species_accessed,
             COUNT(DISTINCT CASE WHEN S.class='Aves' THEN S.species END)      AS n_birds_accessed,
             COUNT(DISTINCT CASE WHEN S.class='Mammalia' THEN S.species END)  AS n_mammals_accessed,
             COUNT(DISTINCT CASE WHEN S.class IN ('%s') THEN S.species END)   AS n_plants_accessed,
             COUNT(DISTINCT RT.route_id)                       AS unique_muni_routes
      FROM cells c JOIN G ON G.cell_id = c.cell_id
                   LEFT JOIN S  ON S.cell_id  = c.cell_id
                   LEFT JOIN RT ON RT.cell_id = c.cell_id
      GROUP BY 1", TT, TT, WALK, TT, WALK, paste(PLANT, collapse="','")))
  }
  d <- rbindlist(acc); d[, time_min := TT]; out[[as.character(TT)]] <- d
  cat(sprintf("T=%2d | %s origins | median cells %s | median species %s | %.1f min\n", TT,
      format(nrow(d),big.mark=","), format(median(d$n_cells_accessed)),
      format(median(d$n_species_accessed)), as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
dbDisconnect(con, shutdown=TRUE)
TR <- rbindlist(out)
TR[, `:=`(mode = "transit",
          metric_support_area_m2 = n_cells_accessed * median(area_m2))]
TR[, `:=`(metric_support_area_km2 = metric_support_area_m2/1e6)]
TR[, `:=`(sampling_density_km2 = n_records_accessed/metric_support_area_km2,
          greenspace_pct = 100*greenspace_area_m2/metric_support_area_m2,
          # verified against the existing reference: transit_access_score is
          # exactly stops per km2 (correlation 1.000)
          transit_access_score = n_stops_accessed/metric_support_area_km2)]
saveRDS(TR, file.path(OUT,"transit_reference_rows.rds"))
cat("WROTE transit_reference_rows.rds\n")
