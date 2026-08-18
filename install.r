options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c(
  
  # Shiny / UI
  "shiny",
  "shinydashboard",
  "shinycssloaders",
  "bslib",
  "DT",
  "leaflet",
  "leaflet.extras",
  "htmltools",
  "base64enc",
  
  # Core data wrangling / plotting
  # Install components directly rather than tidyverse meta-package
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "tibble",
  "stringr",
  "forcats",
  "lubridate",
  "ggplot2",
  "scales",
  
  # Spatial
  "sf",
  "terra",
  "mapboxapi",
  
  # Database
  "DBI",
  "dbplyr",
  "duckdb",
  
  # Transit
  "gtfsrouter",
  "tidytransit",
  "tidycensus",
  
  # Radar plots
  "fmsb",
  
  # Helpers
  "glue",
  "jsonlite",
  "data.table"
))