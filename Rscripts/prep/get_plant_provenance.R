# ============================================================================
# Prep: vascular-plant provenance (native / introduced) from the iNaturalist API
# ============================================================================
# 1) Pull every distinct vascular-plant species (phylum Tracheophyta) from the
#    app's biodiversity parquet (data/output/gbif_census_ndvi_anno.parquet).
# 2) For each species name, query /v1/taxa with preferred_place_id = 854
#    (San Francisco County) and record the iNat taxon id plus the
#    establishment means iNat reports for that place. When SF has no listing,
#    cascade to California (14) then United States (1); the place the label
#    came from is recorded in establishment_place.
# 3) Write data/output/plant_provenance.csv keyed on `species`, which is the
#    only species identifier in the parquet (there is no numeric taxon key).
#
# Runtime: ~1 request/s (iNat's recommended rate) → ~20 min for ~1,200 species.
# Run via Rscripts/prep/get_plant_provenance.sbatch (compute nodes on `main`
# have outbound internet).

library(dplyr)
library(duckdb)
library(httr2)
library(purrr)
library(readr)
library(glue)
library(tictoc)

gbif_parquet <- "data/output/gbif_census_ndvi_anno.parquet"
out_csv <- "data/output/plant_provenance.csv"
sf_place_id <- 854 # iNaturalist place: San Francisco County, CA, US

# ============================================================================
# 1) Distinct vascular plant species in the app data
# ============================================================================

con <- dbConnect(duckdb())
plant_species <- con |>
  dbGetQuery(glue("
    SELECT species, count(*) AS n_records
    FROM read_parquet('{gbif_parquet}')
    WHERE phylum = 'Tracheophyta' AND species IS NOT NULL
    GROUP BY species
    ORDER BY species
  "))
con |> dbDisconnect(shutdown = TRUE)

cat(glue("{nrow(plant_species)} distinct vascular plant species\n\n"))

# ============================================================================
# 2) Query iNaturalist /v1/taxa for each species
# ============================================================================

inat_base <- request("https://api.inaturalist.org/v1/taxa") |>
  req_user_agent("SF_biodiv_access_shiny (aph82@cornell.edu)") |>
  req_throttle(capacity = 1, fill_time_s = 1) |>
  req_retry(max_tries = 4, backoff = ~ 10)

# Places to try, most specific first. iNat only reports an establishment
# means for the preferred place when a listing exists there (or it happens to
# fall back), so cascade outward when SF returns nothing.
place_cascade <- c(san_francisco = 854, california = 14, united_states = 1)

# Normalise a scientific name for comparison: lowercase, drop hybrid sign,
# collapse whitespace ("Crocosmia × crocosmiiflora" -> "crocosmia crocosmiiflora")
norm_name <- function(x) {
  x |> tolower() |> gsub("×", "", x = _) |> gsub("\\s+", " ", x = _) |> trimws()
}

# Return the iNat taxon record for a species name, or NULL.
# Tiers: exact species name; species reached via synonym (matched_term);
# infraspecific/hybrid taxon whose name or matched_term equals the binomial
# (or the nominate trinomial, e.g. "Raphanus sativus sativus").
match_taxon <- function(sp, results) {
  q <- norm_name(sp)
  q_nominate <- paste(q, sub(".* ", "", q))
  name_of <- \(t) norm_name(t$name)
  term_of <- \(t) norm_name(t$matched_term %||% "")

  species <- results |> keep(\(t) identical(t$rank, "species"))
  tiers <- list(
    exact = species |> keep(\(t) name_of(t) == q),
    synonym = species |> keep(\(t) term_of(t) == q),
    infraspecific = results |> keep(\(t)
      name_of(t) == q | term_of(t) %in% c(q, q_nominate))
  )
  tiers <- tiers |> discard(\(x) length(x) == 0)
  if (length(tiers) == 0) return(NULL)

  cand <- tiers[[1]]
  active <- cand |> keep(\(t) isTRUE(t$is_active))
  t <- if (length(active) > 0) active[[1]] else cand[[1]]
  t$match_type <- names(tiers)[1]
  t
}

lookup_species <- function(sp) {
  res <- inat_base |>
    req_url_query(q = sp, preferred_place_id = place_cascade[[1]], per_page = 10) |>
    req_perform() |>
    resp_body_json()

  t <- match_taxon(sp, res$results)
  if (is.null(t)) {
    return(tibble(
      species = sp, inat_taxon_id = NA_integer_, inat_name = NA_character_,
      inat_rank = NA_character_, match_type = "none", is_active = NA,
      establishment_means = NA_character_, establishment_place = NA_character_,
      establishment_place_id = NA_integer_
    ))
  }

  # Cascade outward through places until iNat reports an establishment means
  for (pid in place_cascade[-1]) {
    if (!is.null(t$preferred_establishment_means)) break
    detail <- inat_base |>
      req_url_path_append(t$id) |>
      req_url_query(preferred_place_id = pid) |>
      req_perform() |>
      resp_body_json()
    t$establishment_means <- detail$results[[1]]$establishment_means
    t$preferred_establishment_means <- detail$results[[1]]$preferred_establishment_means
  }
  place <- t$establishment_means$place

  tibble(
    species = sp,
    inat_taxon_id = as.integer(t$id),
    inat_name = t$name,
    inat_rank = t$rank,
    match_type = t$match_type,
    is_active = isTRUE(t$is_active),
    establishment_means = t$preferred_establishment_means %||% NA_character_,
    establishment_place = place$name %||% NA_character_,
    establishment_place_id = as.integer(place$id %||% NA)
  )
}

tic("iNaturalist lookups")
provenance <- plant_species$species |>
  map(lookup_species, .progress = TRUE) |>
  list_rbind() |>
  left_join(plant_species, by = "species") |>
  relocate(n_records, .after = species)
toc()

# ============================================================================
# 3) Write CSV + summary
# ============================================================================

provenance |> write_csv(out_csv, na = "")
cat(glue("\nWrote {out_csv}\n\n"))

provenance |> count(match_type) |> print()
provenance |> count(establishment_means) |> print()
provenance |>
  count(establishment_place, establishment_place_id, sort = TRUE) |>
  print()
