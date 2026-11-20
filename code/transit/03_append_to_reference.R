# Stage 3: append transit to the reference and recompute percentiles.
#
# IMPORTANT: BAS is the percentile of accessible richness WITHIN the same mode
# and the same travel time. Adding transit therefore creates a new percentile
# pool of its own and does NOT move walking, cycling or driving by a single
# point. That is the property that makes "relative percentages" comparable
# across modes: each mode is ranked against itself, so a 70th-percentile
# transit origin and a 70th-percentile walking origin both mean "better than
# 70% of the city, for someone travelling that way".
suppressMessages({library(data.table)})
setwd("~/Downloads/Biodiversity_Access_Score")
REF <- "outputs/results/isochrone_basi_100m.rds"

# Must match percent_rank_existing_reference() in app.R EXACTLY, which is
# 100 * sum(ref < x) / (n - 1): the share strictly BELOW, with no tie averaging.
# Tie-averaged ranks look almost identical until a metric has many ties - with
# a tenth of cells at zero mammals the two conventions differ by 6 percentage
# points, which is what the app's own Step-5 self-check caught.
pctile <- function(x) {
  ok <- is.finite(x); o <- rep(NA_real_, length(x)); n <- sum(ok)
  if (n < 2) return(o)
  o[ok] <- 100 * (data.table::frank(x[ok], ties.method = "min") - 1) / (n - 1); o
}

ref <- as.data.table(readRDS(REF))
old_modes <- sort(unique(ref$mode))
TRANSIT_MODES <- c("transit", "walk_transit")
drop <- intersect(TRANSIT_MODES, old_modes)
if (length(drop)) { ref <- ref[!mode %in% drop]
  cat("replacing existing rows for:", paste(drop, collapse=", "), "\n") }
TR  <- as.data.table(readRDS("outputs/results/transit/transit_reference_rows.rds"))

# snapshot walking BAS so we can prove the other modes are untouched
before <- ref[mode=="walking" & time_min==15, .(cell_id, old = pctile(n_species_accessed))]

for (cl in setdiff(names(ref), names(TR))) TR[[cl]] <- NA
TR <- TR[, names(ref), with=FALSE]
full <- rbind(ref, TR)

# recompute every percentile column within mode x time
pcols <- grep("^pctile_", names(full), value=TRUE)
for (pc in pcols) {
  src <- sub("^pctile_", "", pc)
  if (!src %in% names(full)) next
  if (src == "sf_ej_score") {           # inverted: high burden = low percentile
    full[, (pc) := 100 - pctile(get(src)), by=.(mode, time_min)]
  } else if (src %in% c("nearest_stop_m","calenviro_ci_score","calenviro_traffic_pct")) {
    full[, (pc) := 100 - pctile(get(src)), by=.(mode, time_min)]
  } else {
    full[, (pc) := pctile(get(src)), by=.(mode, time_min)]
  }
}

after <- full[mode=="walking" & time_min==15, .(cell_id, new = pctile_n_species_accessed)]
chk <- merge(before, after, by="cell_id")
cat(sprintf("\nwalking-15 BAS unchanged by rebuilding transit: max shift %.6f points\n",
            max(abs(chk$old - chk$new), na.rm=TRUE)))
stopifnot("rebuilding transit perturbed walking" = max(abs(chk$old - chk$new), na.rm=TRUE) < 1e-9)

cat("\n=== reference now holds ===\n")
print(full[, .(origins=.N, median_species=median(n_species_accessed, na.rm=TRUE),
               median_BAS=round(median(pctile_n_species_accessed, na.rm=TRUE),1)),
           by=.(mode, time_min)][order(mode, time_min)])

bk <- sub("\\.rds$", "_pre_transit.rds", REF)
if (!file.exists(bk)) file.copy(REF, bk)
saveRDS(full, REF)
cat(sprintf("\nWROTE %s (%s rows, %d modes). Backup: %s\n",
            REF, format(nrow(full), big.mark=","), length(unique(full$mode)), basename(bk)))
