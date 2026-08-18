# ============================================================================
# Prep: CalEnviroScreen + SF EJ → GeoPackages in data/output
# ============================================================================
# Reads local inputs from data/source/, writes derived layers for HuggingFace upload.

library(sf)
library(dplyr)

out_dir <- "data/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- CalEnviroScreen (SF tracts only) -----------------------------------------

sf::st_read("data/source/calenviroscreen40gdb_F_2021.gdb", quiet = TRUE) |>
  filter(grepl("san francisco", County, ignore.case = TRUE), !is.na(CIscore)) |>
  select(
    Tract, CIscore, CIscoreP,
    PM2_5, PM2_5_Pctl, Traffic, Traffic_Pctl,
    Poverty, Poverty_Pctl, HousBurd, HousBurd_Pctl,
    County
  ) |>
  sf::st_transform(4326) |>
  sf::st_make_valid() |>
  sf::st_write(file.path(out_dir, "calenviro_sf.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# --- SF Environmental Justice Communities -------------------------------------

sf::st_read(
  "data/source/San Francisco Environmental Justice Communities Map_20260401/geo_export_e303e420-19f7-4166-9736-c1dbeda5b82e.shp",
  quiet = TRUE
) |>
  sf::st_transform(4326) |>
  sf::st_make_valid() |>
  sf::st_write(file.path(out_dir, "sf_ej_communities_map.gpkg"), delete_dsn = TRUE, quiet = TRUE)
