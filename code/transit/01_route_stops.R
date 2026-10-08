# Stage 1: stop -> stop travel times for the whole Muni network.
#
# THE KEY EFFICIENCY POINT: route once per STOP (3,259 calls), not once per
# origin (24,182 x 3 = 72,546 calls). A transit isochrone depends only on which
# stop you board and how much budget is left, so every origin reuses the same
# per-stop routing. ~4.5 minutes instead of many hours, and no Mapbox cost -
# gtfsrouter runs entirely locally.
suppressMessages({library(gtfsrouter); library(data.table)})
setwd("~/Downloads/Biodiversity_Access_Score")
OUT <- "outputs/results/transit"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

DEP <- 9 * 3600          # 09:00 Monday, matching the app's default
MAX_RIDE_MIN <- 15       # cannot exceed the largest travel budget

# gtfs_traveltimes() returns `duration` as an "HH:MM:SS" STRING. as.numeric()
# on it yields NA, and a later NA-valued logical subset silently returns
# all-NA rows - the failure mode that hid this bug for months.
hms_to_min <- function(x) {
  if (inherits(x, "difftime")) return(as.numeric(x, units="mins"))
  p <- do.call(rbind, strsplit(as.character(x), ":", fixed=TRUE))
  (as.numeric(p[,1])*3600 + as.numeric(p[,2])*60 + as.numeric(p[,3]))/60
}

# gtfs_traveltimes() REFUSES to run without a transfers table, and
# gtfs_transfer_table() must be applied to the RAW feed BEFORE gtfs_timetable().
enriched <- "data/cached/gtfs_timetable_monday_transfers.rds"
g <- if (file.exists(enriched)) readRDS(enriched) else {
  raw <- extract_gtfs("data/cached/sf_muni_gtfs.zip", quiet=TRUE)
  raw <- gtfs_transfer_table(raw, d_limit=200, min_transfer_time=120)
  tt  <- gtfs_timetable(raw, day="Monday", quiet=TRUE)
  saveRDS(tt, enriched); tt
}
stopifnot(!is.null(g$transfers))

sid <- as.character(g$stops$stop_id)
cat("stops in feed:", length(sid), "\n")
t0 <- Sys.time(); res <- vector("list", length(sid))
for (i in seq_along(sid)) {
  res[[i]] <- tryCatch({
    tt <- gtfs_traveltimes(g, from=sid[i], from_is_id=TRUE,
            start_time_limits=c(DEP, DEP+1800), max_traveltime=MAX_RIDE_MIN*60)
    if (is.null(tt) || !nrow(tt)) NULL else
      data.frame(from_stop=sid[i], to_stop=as.character(tt$stop_id),
                 ride_min=hms_to_min(tt$duration), stringsAsFactors=FALSE)
  }, error=function(e) NULL)
  if (i %% 500 == 0) cat(sprintf("  %d/%d (%.1f min)\n", i, length(sid),
      as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
B <- rbindlist(res[!vapply(res, is.null, logical(1))])
B <- B[!is.na(ride_min) & !is.na(to_stop) & ride_min <= MAX_RIDE_MIN]
stopifnot(nrow(B) > 0, !anyNA(B))
cat(sprintf("pairs %s | median reachable stops %d | %.1f min\n",
    format(nrow(B), big.mark=","), median(as.integer(table(B$from_stop))),
    as.numeric(difftime(Sys.time(),t0,units="mins"))))
saveRDS(B, file.path(OUT, "stop_to_stop_15min.rds"))
cat("WROTE stop_to_stop_15min.rds\n")
