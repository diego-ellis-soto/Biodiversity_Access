# =============================================================================
# SF Find Neighborhoods viewer - a small Shiny + Leaflet app for the
# neighborhood roll-ups built by code/prep/summarise_sf_find_neighborhoods.R
#
# Run from the repository root:
#   Rscript -e "shiny::runApp('code/neighborhood_viewer')"
# Reads outputs/tables/sf_find_neighborhoods/ (override with NBHD_DIR).
# =============================================================================

library(shiny)
library(bslib)
library(leaflet)
library(sf)
library(DT)

# shiny::runApp() sets the working directory to this folder
nbhd_dir <- Sys.getenv("NBHD_DIR", file.path("..", "..", "outputs", "tables", "sf_find_neighborhoods"))
gpkg <- file.path(nbhd_dir, "sf_find_neighborhoods_summary.gpkg")
if (!file.exists(gpkg)) stop("Run code/prep/summarise_sf_find_neighborhoods.R first: ", gpkg, " not found")

nb <- st_read(gpkg, quiet = TRUE) |> st_transform(4326)

# Sequential blue, light -> dark (same ramp as the standalone HTML map)
RAMP <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
NODATA <- "#e8e6e1"

inside <- data.frame(
  id = c("species", "bird_species", "mammal_species", "plant_species", "pct_plants_native",
         "gbif_records_km2", "mean_ndvi", "greenspace_pct", "transit_stops_km2", "muni_routes",
         "calenviro_ci", "sf_ej_app", "est_population", "med_income_k_popw",
         "pct_nonwhite_popw", "pct_snap_popw"),
  label = c("Species richness", "Bird species", "Mammal species", "Plant species",
            "Native plants (% of plant species)", "GBIF records per sq km",
            "Vegetation (mean NDVI)", "Greenspace cover (%)", "Transit stops per sq km",
            "Muni routes", "CalEnviroScreen burden (higher = more)",
            "SF EJ burden (0-30, higher = more)", "Population",
            "Median household income ($k)", "Non-white residents (%)", "Households on SNAP (%)"),
  digits = c(0, 0, 0, 0, 1, 0, 3, 1, 1, 0, 1, 1, 0, 0, 1, 1),
  stringsAsFactors = FALSE
)
modes <- c("Walking" = "walking", "Cycling" = "cycling", "Driving" = "driving",
           "Transit (no wait)" = "transit", "Walk + transit" = "walk_transit")

metric_choices <- list(
  "Biodiversity Access Index" = c("Mean BAI of origins (0-1)" = "bai"),
  "Inside the neighborhood" = setNames(inside$id, inside$label)
)

ui <- page_sidebar(
  title = "SF neighborhood biodiversity access",
  theme = bs_theme(version = 5, primary = "#256abf"),
  sidebar = sidebar(
    width = 300,
    selectInput("metric", "Colour neighborhoods by", choices = metric_choices, selected = "bai"),
    conditionalPanel(
      "input.metric == 'bai'",
      selectInput("mode", "Travel mode", choices = modes, selected = "walking"),
      radioButtons("time", "Travel time", choices = c("5 min" = 5, "10 min" = 10, "15 min" = 15),
                   selected = 10, inline = TRUE)
    ),
    checkboxInput("fill", "Show colours (untick for outlines only)", TRUE),
    helpText("Colours are quintiles across the 117 SF Find Neighborhoods. Click a",
             "neighborhood on the map, or a row in the table, for details.")
  ),
  card(full_screen = TRUE, card_header(textOutput("map_title")), leafletOutput("map", height = 560)),
  card(card_header("Ranked neighborhoods"), DTOutput("table"))
)

server <- function(input, output, session) {

  sel <- reactive({
    if (identical(input$metric, "bai")) {
      col <- paste0("bai_", input$mode, "_", input$time, "min")
      list(col = col, digits = 2,
           label = sprintf("Mean BAI - %s, %s min", names(modes)[modes == input$mode], input$time))
    } else {
      r <- inside[inside$id == input$metric, ]
      list(col = r$id, digits = r$digits, label = r$label)
    }
  })

  pal_for <- function(x) {
    br <- unique(quantile(x, seq(0, 1, 0.2), na.rm = TRUE, names = FALSE))
    if (length(br) < 3) return(colorBin(RAMP, x, bins = 5, na.color = NODATA))
    colorBin(RAMP, x, bins = br, na.color = NODATA)
  }

  fmt <- function(v, d) ifelse(is.finite(v), formatC(v, format = "f", digits = d, big.mark = ","), "no data")

  popup_html <- function(i) {
    p <- st_drop_geometry(nb[i, ])
    rows <- paste0("<tr><td>", inside$label, "</td><td style='text-align:right'>",
                   vapply(seq_len(nrow(inside)), function(k) fmt(p[[inside$id[k]]], inside$digits[k]), ""),
                   "</td></tr>", collapse = "")
    bai <- paste0("<tr><td>", names(modes), "</td>",
                  vapply(modes, function(m) paste0("<td style='text-align:right;padding-left:10px'>",
                         fmt(unlist(p[paste0("bai_", m, "_", c(5, 10, 15), "min")]), 2), "</td>", collapse = ""), ""),
                  "</tr>", collapse = "")
    hdr <- paste0("<tr><td></td>", paste0("<td style='text-align:right;padding-left:10px;color:#55544f'>",
                  c("5 min", "10 min", "15 min"), "</td>", collapse = ""), "</tr>")
    paste0("<b>", p$neighborhood, "</b><table style='font-size:12px'>", rows,
           "</table><br><b>Mean BAI</b><table style='font-size:12px'>", hdr, bai, "</table>")
  }
  popups <- vapply(seq_len(nrow(nb)), popup_html, "")

  output$map <- renderLeaflet({
    leaflet(nb) |>
      addProviderTiles(providers$Esri.WorldGrayCanvas) |>
      fitBounds(-122.515, 37.708, -122.357, 37.832)
  })

  observe({
    s <- sel()
    x <- nb[[s$col]]
    pal <- pal_for(x)
    show <- isTRUE(input$fill)
    leafletProxy("map", data = nb) |>
      clearShapes() |> clearControls() |>
      addPolygons(
        layerId = ~neighborhood,
        fillColor = if (show) pal(x) else "transparent",
        fillOpacity = if (show) 0.82 else 0,
        color = if (show) "#ffffff" else "#55544f", weight = if (show) 1.2 else 1,
        label = sprintf("%s: %s", nb$neighborhood, fmt(x, s$digits)),
        popup = popups,
        popupOptions = popupOptions(maxHeight = 320, maxWidth = 340),
        highlightOptions = highlightOptions(weight = 3, color = "#0d366b", bringToFront = TRUE)
      ) |>
      (\(m) if (show) addLegend(m, "bottomleft", pal = pal, values = x, title = s$label, opacity = 0.9,
                                labFormat = labelFormat(digits = s$digits)) else m)()
  })

  output$map_title <- renderText(if (isTRUE(input$fill)) sel()$label else "Outlines only")

  tbl <- reactive({
    s <- sel()
    d <- data.frame(neighborhood = nb$neighborhood, value = round(nb[[s$col]], s$digits))
    d <- d[order(-d$value, na.last = TRUE), ]
    d$rank <- seq_len(nrow(d)); d$rank[is.na(d$value)] <- NA
    names(d) <- c("Neighborhood", s$label, "Rank")
    d[, c(3, 1, 2)]
  })
  output$table <- renderDT(datatable(tbl(), rownames = FALSE, selection = "single",
                                     options = list(pageLength = 10, dom = "ftip")))

  observeEvent(input$table_rows_selected, {
    name <- tbl()$Neighborhood[input$table_rows_selected]
    bb <- st_bbox(nb[nb$neighborhood == name, ])
    leafletProxy("map") |>
      flyToBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]],
                  options = list(maxZoom = 15, padding = c(40, 40))) |>
      clearPopups() |>
      addPopups(mean(bb[c("xmin", "xmax")]), mean(bb[c("ymin", "ymax")]),
                popups[nb$neighborhood == name],
                options = popupOptions(maxHeight = 320, maxWidth = 340))
  })
}

shinyApp(ui, server)
