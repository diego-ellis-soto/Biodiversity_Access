# ============================================================================
# Prep: GTFS timetable, stop headways, and feed zip → data/output
# ============================================================================
# Inputs:  data/source/muni_gtfs-current/*.txt
# Outputs: data/output/gtfs_timetable_monday.rds
#          data/output/gtfs_stop_headways.csv
#          data/output/sf_muni_gtfs.zip
#
# Run after updating the GTFS extract. Then upload data/output/* to HuggingFace.

# library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(lubridate)
library(ggplot2)
library(gtfsrouter)
library(tidytransit)
library(glue)

out_dir <- "data/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir_abs <- normalizePath(out_dir, mustWork = TRUE)

gtfs_dir <- "data/source/muni_gtfs-current"
if (!dir.exists(gtfs_dir)) {
  stop(glue("GTFS folder not found: {gtfs_dir}"))
}
gtfs_txt <- list.files(gtfs_dir, pattern = "\\.txt$", full.names = FALSE)
if (length(gtfs_txt) == 0L) {
  stop(glue("No .txt files under {gtfs_dir}"))
}

# --- Zip feed (same bytes consumers will download from HuggingFace) ----------
# Absolute zip path before setwd(); normalizePath(zip) is NA if the file does not exist yet
zip_abs <- file.path(out_dir_abs, "sf_muni_gtfs.zip")
old_wd <- getwd()
setwd(gtfs_dir)
utils::zip(zip_abs, files = gtfs_txt)
setwd(old_wd)

# --- gtfsrouter Monday timetable ------------------------------------------------
gr <- gtfsrouter::extract_gtfs(zip_abs, quiet = TRUE)
tt <- gtfsrouter::gtfs_timetable(gr, day = "Monday")
saveRDS(tt, file.path(out_dir_abs, "gtfs_timetable_monday.rds"), compress = "gzip")

# --- AM peak headways (inspectable CSV) ---------------------------------------
gt <- tidytransit::read_gtfs(zip_abs)
hw <- tidytransit::get_stop_frequency(gt, start_time = 7 * 3600, end_time = 9 * 3600) |>
  group_by(stop_id) |>
  summarise(
    mean_headway_min    = mean(mean_headway, na.rm = TRUE) / 60,
    n_departures_peak   = sum(n_departures, na.rm = TRUE),
    .groups             = "drop"
  ) |>
  mutate(stop_id = as.character(stop_id))

readr::write_csv(hw, file.path(out_dir_abs, "gtfs_stop_headways.csv"))
