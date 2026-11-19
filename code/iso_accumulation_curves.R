# =============================================================================
# Species accumulation curves inside isochrones
# =============================================================================
# Incidence-based rarefaction / extrapolation of species richness (Hill number
# q = 0). The sampling unit is a 100 m grid cell: a species counts as detected
# in a cell when the GBIF cell lookup holds at least one record of it there.
# This is the same approach as the urban-biodiversity app, which ran iNEXT on
# species x H3 hexagon incidence.
#
# Grid cells with at least one GBIF record as the sampling units follows
# Example 4 of Chao et al. (2020), "Quantifying sample completeness and
# comparing diversities among assemblages", Ecological Research 35:292-314,
# which treats GBIF/OBIS coral records the same way (0.01 degree sub-grids).
# The same paper is the reason the expected total is worded "at least": the
# Chao2 estimate is a lower bound while the curve is still rising.
#
# Cells, not observations, are the sampling unit because the GBIF parquet is
# deduplicated (SELECT DISTINCT on year / species / coordinates), so record
# counts are not true abundances.
#
# The estimators are the ones iNEXT uses for incidence data at q = 0 (Colwell
# et al. 2012, J. Plant Ecol. 5:3-21; Chao et al. 2014, Ecol. Monogr. 84:45-67),
# written out here instead of calling iNEXT. iNEXT::iNEXT() and
# iNEXT::estimateD() also compute asymptotic diversity profiles that this plot
# does not use, which costs about 20 s per curve for a city-wide isochrone
# (and that again for every bootstrap replicate) -- too slow for the app. The
# functions below give the same numbers as iNEXT::estimateD() in milliseconds.
#
# Everything here works from objects the app already holds:
#   species_cells -- gbif_cell_species_100m.rds lookup (cell_id, species, class)
#   cell_ids      -- the 100 m cells whose centroid falls inside one isochrone
#
# Requires ISO_PLANT_CLASSES from code/iso_metric_definitions.R.

# Taxonomic groups drawn as the rows of the plot, top to bottom. Each entry
# holds the GBIF classes of the group; NULL keeps every species.
ACCUMULATION_GROUPS <- list(
  "All species" = NULL,
  "Plants"      = ISO_PLANT_CLASSES,
  "Birds"       = "Aves",
  "Mammals"     = "Mammalia"
)

# A curve is drawn only with at least this many sampled cells and species.
ACCUMULATION_MIN_CELLS   <- 5
ACCUMULATION_MIN_SPECIES <- 2

# Line colour of each travel time (the app offers 5, 10 and 15 minutes).
ACCUMULATION_TIME_COLOURS <- c(
  "5 min"  = "#1B9E77",
  "10 min" = "#D95F02",
  "15 min" = "#7570B3"
)

# Chao2 estimate of the number of species present in the sampled area but not
# detected in any cell, from the species detected in exactly one cell (Q1) and
# in exactly two cells (Q2). Observed richness plus this number is the expected
# total richness, the level the accumulation curve climbs towards.
accumulation_n_undetected <- function(cells_detected, n_sampled_cells) {
  n_uniques    <- sum(cells_detected == 1)
  n_duplicates <- sum(cells_detected == 2)
  if (n_duplicates > 0) {
    (n_sampled_cells - 1) / n_sampled_cells * n_uniques^2 / (2 * n_duplicates)
  } else {
    (n_sampled_cells - 1) / n_sampled_cells * n_uniques * (n_uniques - 1) / 2
  }
}

# Expected species richness when `n_cells_levels` cells are sampled.
#   cells_detected  -- for each observed species, the number of cells it was
#                      detected in
#   n_sampled_cells -- the number of cells actually sampled
# Levels up to n_sampled_cells are rarefaction; levels beyond it are
# extrapolation towards the Chao2 richness estimate.
accumulation_richness <- function(cells_detected, n_sampled_cells, n_cells_levels) {
  n_species_observed <- length(cells_detected)

  # Species detected in exactly one cell (Q1) and in exactly two cells (Q2)
  n_uniques    <- sum(cells_detected == 1)
  n_duplicates <- sum(cells_detected == 2)

  # Chao2 estimate of the number of species present but not detected
  n_undetected <- accumulation_n_undetected(cells_detected, n_sampled_cells)

  # Share of the undetected species that each extra cell is expected to reveal
  if (n_uniques > 0) {
    new_species_rate <- n_uniques / (n_sampled_cells * n_undetected + n_uniques)
  } else {
    new_species_rate <- 0
  }

  richness <- numeric(length(n_cells_levels))
  for (level in seq_along(n_cells_levels)) {
    n_cells <- n_cells_levels[level]
    if (n_cells <= n_sampled_cells) {
      # Rarefaction: each species drops out with the probability that none of
      # the n_cells drawn cells is one it was detected in
      prob_missed <- exp(
        lchoose(n_sampled_cells - cells_detected, n_cells) -
          lchoose(n_sampled_cells, n_cells)
      )
      richness[level] <- n_species_observed - sum(prob_missed)
    } else {
      # Extrapolation: observed richness plus the undetected species expected
      # to turn up in the extra cells
      n_extra_cells <- n_cells - n_sampled_cells
      richness[level] <- n_species_observed +
        n_undetected * (1 - (1 - new_species_rate)^n_extra_cells)
    }
  }
  richness
}

# Per-cell detection probabilities of the assemblage that the bootstrap
# resamples: every observed species (its detection rate, shrunk a little to
# leave room for unseen species) followed by the Chao2 number of undetected
# species sharing the remaining probability. Same construction as iNEXT.
accumulation_bootstrap_probabilities <- function(cells_detected, n_sampled_cells) {
  n_uniques    <- sum(cells_detected == 1)
  n_duplicates <- sum(cells_detected == 2)

  # Chao2 estimate of the number of species present but not detected
  n_undetected <- accumulation_n_undetected(cells_detected, n_sampled_cells)

  detection_rate <- cells_detected / n_sampled_cells
  if (n_undetected == 0) return(detection_rate)

  # Total detection probability set aside for the undetected species
  prob_undetected_total <- n_uniques / n_sampled_cells *
    n_sampled_cells * n_undetected / (n_sampled_cells * n_undetected + n_uniques)

  # Shrink the observed species' rates by that amount, mostly the rare ones
  shrink_weight <- prob_undetected_total /
    sum(detection_rate * (1 - detection_rate)^n_sampled_cells)
  prob_observed <- detection_rate *
    (1 - shrink_weight * (1 - detection_rate)^n_sampled_cells)

  n_undetected_whole <- ceiling(n_undetected)
  prob_undetected <- rep(prob_undetected_total / n_undetected_whole, n_undetected_whole)

  c(prob_observed, prob_undetected)
}

# One accumulation curve: species richness against the number of 100 m cells
# sampled, for one isochrone and one taxonomic group. Returns one row per point
# on the curve, or NULL when there are too few cells or species to draw one.
iso_accumulation_curve <- function(species_cells, cell_ids, group_classes, nboot = 50) {
  # Species detections in the cells of this isochrone
  detections <- species_cells |>
    filter(cell_id %in% cell_ids, !is.na(species))

  # Keep only the classes of this taxonomic group
  if (!is.null(group_classes)) {
    detections <- detections |>
      filter(class %in% group_classes)
  }

  # One row per species per cell
  detections <- detections |>
    distinct(cell_id, species)

  # Number of cells each species was detected in
  cells_detected <- count(detections, species)$n

  # Sampled cells = cells with at least one record of this group
  n_sampled_cells <- n_distinct(detections$cell_id)

  if (n_sampled_cells < ACCUMULATION_MIN_CELLS) return(NULL)
  if (length(cells_detected) < ACCUMULATION_MIN_SPECIES) return(NULL)

  # Points along the curve: rarefaction from 1 cell up to the sampled cells,
  # then extrapolation to twice the sampled cells (iNEXT's default endpoint)
  n_cells_levels <- c(
    seq(1, n_sampled_cells, length.out = 30),
    seq(n_sampled_cells, 2 * n_sampled_cells, length.out = 30)
  ) |>
    round() |>
    unique()

  richness <- accumulation_richness(cells_detected, n_sampled_cells, n_cells_levels)

  # Bootstrap: redraw how many cells each species is detected in, recompute
  # the curve, and take the spread across replicates as the standard error
  boot_probabilities <- accumulation_bootstrap_probabilities(cells_detected, n_sampled_cells)
  boot_richness <- matrix(NA_real_, nrow = nboot, ncol = length(n_cells_levels))
  for (replicate in seq_len(nboot)) {
    boot_detected <- rbinom(length(boot_probabilities), n_sampled_cells, boot_probabilities)
    boot_detected <- boot_detected[boot_detected > 0]
    boot_richness[replicate, ] <- accumulation_richness(boot_detected, n_sampled_cells, n_cells_levels)
  }
  richness_se <- apply(boot_richness, 2, sd)

  tibble(
    n_cells      = n_cells_levels,
    method       = case_when(
      n_cells_levels < n_sampled_cells ~ "Rarefaction",
      n_cells_levels == n_sampled_cells ~ "Observed",
      TRUE ~ "Extrapolation"
    ),
    richness     = richness,
    richness_lcl = richness - qnorm(0.975) * richness_se,
    richness_ucl = richness + qnorm(0.975) * richness_se,
    # Expected total: observed species plus the Chao2 number of undetected ones
    # (the same value on every row of the curve)
    richness_expected = length(cells_detected) +
      accumulation_n_undetected(cells_detected, n_sampled_cells)
  )
}

# The curves of every taxonomic group for one isochrone, stacked, with a
# `group` column. Groups with too little data are left out.
iso_accumulation_curves_by_group <- function(species_cells, cell_ids, nboot = 50) {
  curves <- list()
  for (group_name in names(ACCUMULATION_GROUPS)) {
    curve <- iso_accumulation_curve(
      species_cells,
      cell_ids,
      group_classes = ACCUMULATION_GROUPS[[group_name]],
      nboot = nboot
    )
    if (is.null(curve)) next
    curves[[group_name]] <- curve |>
      mutate(group = group_name)
  }
  bind_rows(curves)
}

# Species recorded and cells sampled in one isochrone, for every taxonomic
# group: one row per group with n_cells (100 m cells with at least one record
# of the group) and richness (species recorded). Cheap, so it is computed for
# every isochrone, including those too sparse for a curve.
iso_accumulation_observed_by_group <- function(species_cells, cell_ids) {
  # Species detections in the cells of this isochrone
  detections <- species_cells |>
    filter(cell_id %in% cell_ids, !is.na(species))

  observed <- list()
  for (group_name in names(ACCUMULATION_GROUPS)) {
    group_classes <- ACCUMULATION_GROUPS[[group_name]]
    group_detections <- detections
    if (!is.null(group_classes)) {
      group_detections <- detections |>
        filter(class %in% group_classes)
    }
    observed[[group_name]] <- tibble(
      group    = group_name,
      n_cells  = n_distinct(group_detections$cell_id),
      richness = n_distinct(group_detections$species)
    )
  }
  bind_rows(observed)
}

# One column of the grid: one transport mode, one panel per taxonomic group.
# Each panel has one accumulation curve per travel time, each with a point at
# the species actually recorded. Only the longest travel time is labelled with
# its expected total: it has the most data, and separate estimates for nested
# isochrones can contradict each other (a 5 min estimate above the 10 min one).
#   curves   -- curve points of this mode, every travel time (see
#               plot_iso_accumulation_curves)
#   observed -- one row per travel time x group of this mode: time, group,
#               n_cells, richness
#   time_levels -- travel-time labels of the whole figure, so that every column
#                  carries the same legend and they merge into one
#   n_columns   -- number of mode columns in the whole figure
plot_accumulation_mode_column <- function(mode_name, curves, observed, time_levels, n_columns) {
  group_levels <- names(ACCUMULATION_GROUPS)

  # With many modes side by side the columns are narrow: smaller type, and
  # long mode names wrapped onto two lines
  narrow <- n_columns > 4
  heading_size <- if (narrow) 14 else 18
  label_size   <- if (narrow) 3.1 else 3.5

  # Factor order sets the panel order and the legend
  curves <- curves |>
    mutate(
      group      = factor(group, levels = group_levels),
      time_label = factor(glue::glue("{time} min"), levels = time_levels)
    )
  observed <- observed |>
    filter(n_cells > 0) |>
    mutate(
      group      = factor(group, levels = group_levels),
      time_label = factor(glue::glue("{time} min"), levels = time_levels)
    )

  # The observed point belongs to both line segments so they join up
  rarefaction   <- curves |> filter(method %in% c("Rarefaction", "Observed"))
  extrapolation <- curves |> filter(method %in% c("Observed", "Extrapolation"))

  # One label per panel, for the longest travel time: the expected total
  # richness and the sample completeness (species recorded / species expected).
  # Worded "at least" and "at most" because the Chao2 estimate is a lower bound
  # while the curve is still rising (Chao et al. 2020)
  curve_labels <- curves |>
    filter(method == "Observed", time == max(time)) |>
    mutate(
      completeness = richness / richness_expected,
      label = glue::glue(
        "Within {time} min:\n",
        "at least {scales::comma(round(richness_expected))} species expected\n",
        "at most {scales::percent(completeness, accuracy = 1)} recorded so far"
      )
    )

  # Panels of groups with too little data for a curve
  no_curve <- tibble(group = setdiff(group_levels, as.character(curves$group))) |>
    mutate(group = factor(group, levels = group_levels))

  column_plot <- ggplot(mapping = aes(x = n_cells, y = richness, colour = time_label, fill = time_label)) +
    # show.legend = TRUE keeps a legend key for travel times this mode lacks
    geom_ribbon(
      data = curves,
      aes(ymin = richness_lcl, ymax = richness_ucl),
      alpha = 0.18, colour = NA, show.legend = TRUE
    ) +
    geom_line(data = rarefaction, linewidth = 0.9, show.legend = TRUE) +
    geom_line(data = extrapolation, linewidth = 0.9, linetype = "dashed", show.legend = FALSE) +
    # Points come from `observed`, so a travel time too sparse for a curve
    # still shows what was recorded there
    geom_point(data = observed, size = 2.6, show.legend = TRUE) +
    geom_text(
      data = curve_labels,
      aes(x = Inf, y = -Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.04, vjust = -0.3, size = label_size, lineheight = 0.95, colour = "grey25"
    ) +
    scale_colour_manual(values = ACCUMULATION_TIME_COLOURS[time_levels], name = "Travel time", drop = FALSE) +
    scale_fill_manual(values = ACCUMULATION_TIME_COLOURS[time_levels], name = "Travel time", drop = FALSE) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(vars(group), ncol = 1, scales = "free", drop = FALSE) +
    labs(
      x = "100 m cells sampled",
      y = "Species richness",
      title = stringr::str_wrap(mode_name, width = 18)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, size = heading_size, face = "bold", lineheight = 0.9),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )

  # Every group has its curve: nothing more to label. (ggplot fails to draw
  # the panels if these layers are added with no rows.)
  if (nrow(no_curve) == 0) return(column_plot)

  column_plot +
    # A panel with no records at all has no axis range; two invisible points
    # give it a 0-1 range so that its label is drawn
    geom_blank(data = no_curve, aes(x = 0, y = 0), inherit.aes = FALSE) +
    geom_blank(data = no_curve, aes(x = 1, y = 1), inherit.aes = FALSE) +
    geom_text(
      data = no_curve,
      aes(x = -Inf, y = Inf),
      label = "Too few records for a curve",
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.6, size = 3.4, colour = "grey40"
    )
}

# Grid of accumulation curves: one column per transport mode (headed by the
# mode), one row per taxonomic group. Every panel has its own axes, because
# the number of cells and species differs by orders of magnitude between, say,
# all species by car and mammals on foot.
#   curves   -- one row per curve point: mode, time, group, n_cells, method,
#               richness, richness_lcl, richness_ucl, richness_expected
#   observed -- one row per isochrone x group: mode, time, group, n_cells,
#               richness. Lays out the columns and supplies the points.
plot_iso_accumulation_curves <- function(curves, observed) {
  mode_levels <- unique(observed$mode)

  # Travel times on the map, in legend order
  time_levels <- intersect(
    names(ACCUMULATION_TIME_COLOURS),
    glue::glue("{unique(observed$time)} min")
  )

  # One column plot per mode
  column_plots <- list()
  for (mode_name in mode_levels) {
    column_plots[[mode_name]] <- plot_accumulation_mode_column(
      mode_name,
      curves      = curves |> filter(mode == mode_name),
      observed    = observed |> filter(mode == mode_name),
      time_levels = time_levels,
      n_columns   = length(mode_levels)
    )
  }

  # Columns side by side, with one shared legend and one set of axis titles
  patchwork::wrap_plots(column_plots, nrow = 1) +
    patchwork::plot_layout(guides = "collect", axis_titles = "collect") +
    patchwork::plot_annotation(
      caption = paste(
        "Solid line: rarefaction. Point: species recorded. Dashed line: extrapolation to twice the sampled cells. Band: 95% bootstrap confidence interval.",
        "Species expected: Chao2 estimate of all species present, recorded or not; a lower bound, shown for the longest travel time of each mode.",
        "Cells sampled: 100 m cells with at least one record of the group.",
        "Estimators as in iNEXT (incidence data, q = 0).",
        sep = "\n"
      ),
      theme = theme(plot.caption = element_text(hjust = 0, colour = "grey35", size = 10))
    ) &
    theme(legend.position = "bottom")
}
