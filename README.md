---
title: Sf Biodiv Access Shiny
emoji: 📚
colorFrom: blue
colorTo: yellow
sdk: docker
pinned: false
---

# SF Biodiversity Access Decision Support Tool

<img src="www/Combined_logos.png" width="60%">

An interactive Shiny app built for the **Reimagining San Francisco (RSF) Initiative** to explore how equitably San Francisco residents can access urban biodiversity across different transportation modes and socioeconomic contexts.

---

## What the App Does

Users select a location anywhere in San Francisco — either by clicking the map or geocoding an address — and choose one or more transportation modes and travel-time thresholds (5, 10, or 15 minutes). The app then:

1. **Generates travel-time isochrones** — polygons representing the area reachable within the chosen time budget — using Mapbox (driving, walking, cycling, traffic-aware driving) and SF Muni GTFS data (transit, walk + transit).
2. **Summarizes biodiversity** within each isochrone: GBIF occurrence records, unique species richness, and breakdowns by taxonomic group (birds, mammals, plants).
3. **Summarizes greenspace access**: OSM greenspace coverage (%), distance to nearest greenspace via raster, and NDVI from Sentinel-2.
4. **Summarizes socioeconomic and environmental justice context**: area-weighted median income and population from ACS Census block groups, CalEnviroScreen Cumulative Impact scores, and SF Environmental Justice community burden scores.
5. **Computes the Biodiversity Access Index (BAI)** — a composite score benchmarked against citywide empirical distributions — and displays it as a spider/radar plot across seven dimensions: mobility access, route access, biodiversity potential, observation intensity, environmental quality, greenspace cover, and equity context.
6. **Displays partner RSF Program Projects** as a toggleable map layer.

---

## App Tabs

| Tab | Description |
|-----|-------------|
| **Isochrone Explorer** | Interactive map, isochrone generation, BAI spider plot, summary table, and metric plots |
| **Isochrone Comparer** | Pick two locations (one mode + travel time each) and compare them side by side: two BAI spider plots plus a difference table (transit score, biodiversity score, BAI, closest greenspace). Loads only when opened. |
| **GBIF Summaries** | Filter GBIF records by taxonomic class and family within the isochrone; species richness vs. sampling effort plot |
| **Community Science** | Map and table of partner community organizations |
| **About** | Full methodology, data sources, transport mode descriptions, and BAI explanation |

---

## Transportation Modes

| Mode | Engine | Notes |
|------|--------|-------|
| Driving | Mapbox Navigation API | Free-flow road network |
| Walking | Mapbox Navigation API | Pedestrian paths and crossings |
| Cycling | Mapbox Navigation API | Dedicated cycle lanes where available |
| Driving with Traffic | Mapbox Navigation API | Real-time + historical congestion |
| Transit (GTFS) | gtfsrouter | SF Muni timetable-based stop-to-stop routing |
| Walk + Transit (Muni) | Mapbox + gtfsrouter | First-mile walk + Muni ride + last-mile walk, all within total time budget |

---

## Map Layers

All layers are toggleable in the map's layer control panel:

- **Income** — Median household income per census block group (ACS 5-yr)
- **Greenspace** — OSM parks and green areas
- **Greenspace Distance** — Raster showing distance (m) to nearest greenspace
- **RSF Program Distance** — Raster showing distance (m) to nearest RSF program polygon
- **RSF Program Projects** — Partner project areas from the RSF Initiative
- **Hotspots / Coldspots (KnowBR)** — Block groups with anomalously high/low species richness relative to sampling effort
- **Species Richness** — Unique GBIF species per census block group
- **Data Availability** — GBIF occurrence records per block group
- **CalEnviroScreen (CI Score)** — Cumulative environmental burden by census tract
- **SF EJ Communities** — SF Environmental Justice community burden scores
- **Transit Routes** — All SF Muni routes from GTFS shapes, colored by official SFMTA route color
- **Transit Stops** — All SF Muni stops with AM peak headway info
- **Isochrones** — Generated travel-time polygons
- **NDVI Raster** — Sentinel-2 NDVI cropped to the isochrone union

---

## Biodiversity Access Index (BAI)

The BAI is a composite indicator benchmarked against **citywide empirical distributions** (ECDFs across all SF census block groups), not just the current session. It comprises seven equally-weighted dimensions:

| Dimension | Variable | Direction |
|-----------|----------|-----------|
| Mobility Access | Transit stops / km² | Higher = better |
| Route Access | Unique Muni routes crossing isochrone | Higher = better |
| Biodiversity Potential | Unique GBIF species | Higher = better |
| Observation Intensity | GBIF records / km² | Higher = better |
| Environmental Quality | Mean NDVI | Higher = better |
| Greenspace Cover | % OSM greenspace area | Higher = better |
| Equity Context | SF EJ burden score (inverted) | Lower burden = better |

All axes scaled 0–1. BAI = unweighted mean of all seven standardized components.

---

## Repository Structure

```
SF_biodiv_access_shiny/
├── app.R                 # Main app (sources Rscripts/setup_unified.R at startup)
├── Dockerfile            # HF Spaces: install.r + shiny::runApp('app.R', …)
├── install.r             # R package list for Docker
├── www/                  # Static assets (e.g. app_pastel.css, logos)
├── Rscripts/
│   ├── setup_unified.R   # Loads all app data: local data/cached + HuggingFace fallback
│   ├── iso_metrics.R     # Shared point→isochrone→metrics→BAI→radar functions (used by both isochrone tabs)
│   └── prep/             # One-off builds → data/output/ (see run_all_prep.R)
└── data/
    ├── cached/           # Downloaded / runtime cache (often gitignored)
    ├── output/           # Prep outputs (often gitignored)
    └── source/           # Raw inputs for prep (e.g. RSF polygons, GTFS extract)
```

> **Note:** The `data/` directory is gitignored. See setup instructions below.

---

## Cyberinfrastructure & Data Sources

| Dataset | Source | Format | How Loaded |
|---------|--------|--------|-----------|
| GBIF occurrences (SF) | Global Biodiversity Information Facility | Parquet | DuckDB spatial queries in server |
| Census block groups + ACS | US Census / tidycensus | `.Rdata` | Downloaded from HuggingFace at startup |
| NDVI raster | Sentinel-2 (pre-processed) | GeoTIFF | Downloaded from HuggingFace at startup |
| OSM greenspace polygons | OpenStreetMap | Shapefile bundle | Downloaded from HuggingFace at startup |
| Greenspace distance rasters | Derived from OSM (see `making-greenspace-raster.R`) | GeoTIFF | Downloaded from HuggingFace at startup |
| RSF Program distance rasters | Derived from RSF polygons (see `making-rsfprogram-raster.R`) | GeoTIFF | Downloaded from HuggingFace at startup |
| Biodiversity hotspots/coldspots | KnowBR analysis on GBIF | Shapefile | Downloaded from HuggingFace at startup |
| SF Muni GTFS | SFMTA official GTFS feed | zip + rds + csv | Downloaded from HuggingFace at startup |
| CalEnviroScreen 4.0 | OEHHA | GeoPackage | Downloaded from HuggingFace at startup |
| SF EJ Communities | SF Environment | GeoPackage | Downloaded from HuggingFace at startup |
| RSF Program Projects | RSF Initiative | GeoPackage | Downloaded from HuggingFace at startup |

**Remote data** is hosted on HuggingFace at  
[`boettiger-lab/sf_biodiv_access`](https://huggingface.co/datasets/boettiger-lab/sf_biodiv_access) and cached locally in `data/cached/` by `setup_unified.R`. The two most expensive products — the CBG × greenspace intersection and the transit-routing timetable — are precomputed in `Rscripts/prep/` and read at startup rather than recomputed, so the app starts quickly.

**GBIF queries** use [DuckDB](https://duckdb.org/) with the spatial extension, querying a local `.parquet` file directly via SQL `ST_Intersects` — avoiding loading the full dataset into memory.

---

## Local Development Setup

### Prerequisites

Install required R packages:

```r
install.packages(c(
  "shiny", "shinydashboard", "leaflet", "mapboxapi", "tidyverse",
  "tidycensus", "sf", "DT", "RColorBrewer", "terra", "data.table",
  "mapview", "sjPlot", "sjlabelled", "bslib", "shinycssloaders",
  "DBI", "duckdb", "dbplyr", "gtfsrouter", "tidytransit", "fmsb", "scales"
))
```

### First-Time Setup

**Step 1: Build GBIF parquets** (one-time; needs local `gbif.duckdb` path in script)
```r
source("Rscripts/prep/create_annotated_gbif_parquet.R")
# Creates: data/output/sf-gbif.parquet, data/output/gbif_census_ndvi_anno.parquet
```

**Step 2: Pre-compute GTFS timetable and greenspace / RSF rasters** (one-time, slow)
```r
source("Rscripts/prep/run_all_prep.R")
```

**Step 3: Run the app**
```r
shiny::runApp("app.R")
```

### Startup Performance

Typical startup time: **~6–12 seconds** depending on whether caches exist.

---

## Cloud Deployment (HuggingFace Spaces)

The app is deployed via Docker on HuggingFace Spaces using `app.R` + `Rscripts/setup_unified.R`. Data is downloaded from the HuggingFace dataset repository into `data/cached/` at startup — so no large files need to be bundled in the image.

See `Dockerfile` for the deployment configuration.

---

## Authors

**Diego Ellis Soto**, Avery Hill, Christopher J. Schell, Carl Boettiger, Rebecca Johnson  
University of California Berkeley (ESPM) · California Academy of Sciences  
Contact: [diego.ellissoto@berkeley.edu](mailto:diego.ellissoto@berkeley.edu)

---

## Status & Roadmap

This tool is a **decision-support prototype** co-developed with the RSF Data Working Group. The BAI should be treated as an exploratory indicator; variable weights and reference distributions are subject to revision through stakeholder co-development.

**Planned additions:**
- Impervious surface coverage
- National Walkability Index
- CDC Social Vulnerability Index
- NatureServe biodiversity / rarity maps
- Frequency-weighted multimodal transit accessibility
- Pre-cached transit isochrones for faster queries

---

<img src="www/hexbin_RSF_logo.png" width="80">
