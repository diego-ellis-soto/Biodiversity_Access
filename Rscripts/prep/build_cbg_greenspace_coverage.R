# ============================================================================
# Prep: Per-CBG greenspace coverage -> data/output
# ============================================================================
# Inputs:  cbg_vect_sf + osm_greenspace (loaded via Rscripts/setup_unified.R,
#          which pulls them from data/cached/ or HuggingFace)
# Output:  data/output/cbg_greenspace_coverage.csv  (GEOID, greenspace_m2, cbg_area_m2)
#
# The CBG x greenspace st_union + st_intersection is the most expensive single
# step at app startup. The inputs are static, so we precompute it here once and
# let setup_unified.R just read this small CSV instead of recomputing every run.
# Re-run when cbg_vect_sf or the greenspace polygons change, then upload
# data/output/cbg_greenspace_coverage.csv to HuggingFace (see upload_to_huggingface.R).
#
# Working directory must be the project root (folder containing data/ and Rscripts/).
# ============================================================================

# Reuse the app's canonical loader so coverage is computed from the exact same
# cbg_vect_sf + osm_greenspace the app uses (also gives us sf/dplyr, etc.).
source("Rscripts/setup_unified.R")

out_dir <- "data/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cbg_proj <- st_transform(cbg_vect_sf[, "GEOID"], 3857) |>
  mutate(cbg_area_m2 = as.numeric(st_area(geometry)))
gs_proj <- st_transform(osm_greenspace, 3857) |> st_make_valid()
gs_union <- st_union(gs_proj)
cbg_gs_inter <- st_intersection(cbg_proj, gs_union)

cbg_greenspace_coverage <- cbg_gs_inter |>
  mutate(greenspace_m2 = as.numeric(st_area(geometry))) |>
  st_drop_geometry() |>
  group_by(GEOID) |>
  summarise(greenspace_m2 = sum(greenspace_m2), .groups = "drop") |>
  right_join(cbg_proj |> st_drop_geometry() |> dplyr::select(GEOID, cbg_area_m2), by = "GEOID") |>
  mutate(
    greenspace_m2 = tidyr::replace_na(greenspace_m2, 0),
    GEOID = as.character(GEOID)
  )

readr::write_csv(cbg_greenspace_coverage, file.path(out_dir, "cbg_greenspace_coverage.csv"))
message(glue::glue(
  "[prep] wrote {nrow(cbg_greenspace_coverage)} CBG rows -> data/output/cbg_greenspace_coverage.csv"
))
