###############################################################################
# build_neighborhood_map_html.R — standalone interactive map of the SF Find
# Neighborhood summaries (one HTML file, data inline, Leaflet from cdnjs)
#
# Run from the repository root after summarise_sf_find_neighborhoods.R:
#   Rscript code/prep/build_neighborhood_map_html.R \
#     [outputs/tables/sf_find_neighborhoods]
#
# Writes <dir>/sf_find_neighborhoods_map.html. Open it in any browser (needs
# internet for the basemap tiles and the Leaflet library).
###############################################################################

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(jsonlite); library(rmapshaper)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[[1]] else "outputs/tables/sf_find_neighborhoods"
gpkg <- file.path(in_dir, "sf_find_neighborhoods_summary.gpkg")
out_html <- file.path(in_dir, "sf_find_neighborhoods_map.html")
stopifnot(file.exists(gpkg))

nb <- st_read(gpkg, quiet = TRUE)

# Topology-aware simplification keeps shared borders shared (no slivers/gaps)
nb_s <- ms_simplify(nb, keep = 0.35, keep_shapes = TRUE) |> st_transform(4326)

# ---------------------------------------------------------------------------
# Metric catalogue: what can be switched on, how to label and format it
# ---------------------------------------------------------------------------
inside <- list(
  list(id = "species",          label = "Species richness",              unit = "unique species", digits = 0),
  list(id = "bird_species",     label = "Bird species",                  unit = "unique species", digits = 0),
  list(id = "mammal_species",   label = "Mammal species",                unit = "unique species", digits = 0),
  list(id = "plant_species",    label = "Plant species",                 unit = "unique species", digits = 0),
  list(id = "pct_plants_native",label = "Native plants",                 unit = "% of plant species", digits = 1),
  list(id = "gbif_records_km2", label = "GBIF records per km&sup2;",     unit = "records / km&sup2;", digits = 0),
  list(id = "mean_ndvi",        label = "Vegetation (mean NDVI)",        unit = "NDVI", digits = 3),
  list(id = "greenspace_pct",   label = "Greenspace cover",              unit = "% of area", digits = 1),
  list(id = "transit_stops_km2",label = "Transit stops per km&sup2;",    unit = "stops / km&sup2;", digits = 1),
  list(id = "muni_routes",      label = "Muni routes",                   unit = "unique routes", digits = 0),
  list(id = "calenviro_ci",     label = "CalEnviroScreen burden",        unit = "CI score, higher = more burden", digits = 1),
  list(id = "sf_ej_app",        label = "SF EJ burden",                  unit = "score 0&ndash;30, higher = more burden", digits = 1),
  list(id = "est_population",   label = "Population",                    unit = "residents", digits = 0),
  list(id = "med_income_k_popw",label = "Median household income",       unit = "$k, population-weighted", digits = 0),
  list(id = "pct_nonwhite_popw",label = "Non-white residents",           unit = "%", digits = 1),
  list(id = "pct_snap_popw",    label = "Households on SNAP",            unit = "%", digits = 1)
)
modes <- c(walking = "Walking", cycling = "Cycling", driving = "Driving",
           transit = "Transit (no wait)", walk_transit = "Walk + transit")
times <- c(5, 10, 15)
bai_cols <- as.vector(outer(names(modes), times, function(m, t) paste0("bai_", m, "_", t, "min")))
stopifnot(all(bai_cols %in% names(nb_s)), all(vapply(inside, function(m) m$id %in% names(nb_s), TRUE)))

keep <- c("neighborhood", "polygon_area_km2", "n_cells", "gbif_records",
          vapply(inside, `[[`, "", "id"), bai_cols)
nb_out <- nb_s[, keep]

# Round attributes for a compact file
num <- vapply(st_drop_geometry(nb_out), is.numeric, TRUE)
for (cl in names(num)[num]) nb_out[[cl]] <- round(nb_out[[cl]], 4)

tmp <- tempfile(fileext = ".geojson")
st_write(nb_out, tmp, driver = "GeoJSON", quiet = TRUE,
         layer_options = c("COORDINATE_PRECISION=5", "RFC7946=YES"))
geojson <- paste(readLines(tmp, warn = FALSE), collapse = "")
unlink(tmp)

meta <- list(
  inside = inside,
  modes = lapply(names(modes), function(k) list(id = k, label = modes[[k]])),
  times = times
)
meta_json <- toJSON(meta, auto_unbox = TRUE)

# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------
html <- '<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SF Neighborhood Biodiversity Access</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css">
<style>
  :root {
    --surface: #ffffff;
    --surface-2: #f6f5f2;
    --ink: #1f1f1d;
    --ink-2: #55544f;
    --ink-3: #8a8983;
    --rule: #e2e0da;
    --accent: #256abf;
    --nodata: #e8e6e1;
    --hover: #0d366b;
  }
  * { box-sizing: border-box; }
  html, body { height: 100%; margin: 0; background: var(--surface); color: var(--ink);
    font: 14px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
  #app { display: flex; height: 100%; }
  #panel { width: 300px; flex: 0 0 300px; overflow-y: auto; border-right: 1px solid var(--rule);
    padding: 16px; background: var(--surface); }
  #map { flex: 1 1 auto; }
  h1 { font-size: 16px; margin: 0 0 4px; }
  .sub { color: var(--ink-2); font-size: 12px; margin: 0 0 14px; }
  .group { margin: 14px 0 6px; font-size: 11px; font-weight: 600; letter-spacing: .04em;
    text-transform: uppercase; color: var(--ink-3); }
  .opt { display: flex; align-items: center; gap: 8px; padding: 5px 6px; border-radius: 6px; cursor: pointer; }
  .opt:hover { background: var(--surface-2); }
  .opt input { margin: 0; accent-color: var(--accent); }
  .opt.on { background: #eaf2fd; }
  .row { display: flex; gap: 6px; margin: 6px 0 0 26px; }
  select { font: inherit; font-size: 12px; padding: 3px 4px; border: 1px solid var(--rule);
    border-radius: 6px; background: var(--surface); color: var(--ink); }
  .note { color: var(--ink-3); font-size: 11px; margin-top: 16px; border-top: 1px solid var(--rule); padding-top: 10px; }
  .legend { background: var(--surface); padding: 10px 12px; border-radius: 8px;
    box-shadow: 0 1px 4px rgba(0,0,0,.18); min-width: 180px; }
  .legend .t { font-weight: 600; font-size: 13px; }
  .legend .u { color: var(--ink-2); font-size: 11px; margin-bottom: 6px; }
  .legend .k { display: flex; align-items: center; gap: 8px; font-size: 12px; color: var(--ink); margin: 2px 0; }
  .legend .sw { width: 22px; height: 12px; border-radius: 3px; border: 1px solid rgba(0,0,0,.08); }
  .leaflet-tooltip.nb-tip { font-size: 12px; padding: 4px 8px; }
  .pop h3 { margin: 0 0 6px; font-size: 14px; }
  .pop table { border-collapse: collapse; font-size: 12px; width: 100%; }
  .pop td { padding: 2px 6px 2px 0; vertical-align: top; }
  .pop td.v { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
  .pop .sec { font-weight: 600; color: var(--ink-3); font-size: 10px; text-transform: uppercase;
    letter-spacing: .04em; padding-top: 8px; }
  .pop tr.cur td { font-weight: 700; }
  .pop .bai td, .pop .bai th { text-align: right; padding: 1px 4px; font-variant-numeric: tabular-nums; }
  .pop .bai th { font-weight: 600; color: var(--ink-2); font-size: 11px; }
  .pop .bai td:first-child, .pop .bai th:first-child { text-align: left; }
  .leaflet-popup-content { max-height: 360px; overflow-y: auto; min-width: 260px; }
  @media (max-width: 720px) {
    #app { flex-direction: column; }
    #panel { width: 100%; flex: 0 0 auto; max-height: 42%; border-right: 0; border-bottom: 1px solid var(--rule); }
  }
</style>
</head>
<body>
<div id="app">
  <aside id="panel">
    <h1>SF neighborhood biodiversity access</h1>
    <p class="sub">117 SF Find Neighborhoods. Choose a value to colour the map; click a neighborhood for all its values.</p>
    <label class="opt" id="opt-none"><input type="radio" name="metric" value="__none"> Outlines only (colouring off)</label>
    <div class="group">Biodiversity Access Index</div>
    <label class="opt" id="opt-bai"><input type="radio" name="metric" value="__bai"> Mean BAI of origins (0–1)</label>
    <div class="row"><select id="mode"></select><select id="time"></select></div>
    <div class="group">Inside the neighborhood</div>
    <div id="inside"></div>
    <p class="note">Data: SF Biodiversity Access 100 m grid and precomputed isochrones, rolled up to SF Find Neighborhoods (data.sf.gov). Colours are quintiles of the selected value across the 117 neighborhoods. Built by code/prep/build_neighborhood_map_html.R.</p>
  </aside>
  <div id="map"></div>
</div>
<script type="application/json" id="nb-data">__GEOJSON__</script>
<script type="application/json" id="nb-meta">__META__</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<script>
(function () {
  const data = JSON.parse(document.getElementById("nb-data").textContent);
  const meta = JSON.parse(document.getElementById("nb-meta").textContent);
  // Sequential blue, light -> dark (validated reference ramp: steps 100, 250, 400, 550, 700)
  const RAMP = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"];
  const NODATA = getComputedStyle(document.documentElement).getPropertyValue("--nodata").trim();

  const byId = {};
  meta.inside.forEach(m => byId[m.id] = m);
  const modeLabel = {}; meta.modes.forEach(m => modeLabel[m.id] = m.label);

  // ---- panel ----
  const insideBox = document.getElementById("inside");
  meta.inside.forEach(m => {
    const l = document.createElement("label");
    l.className = "opt"; l.id = "opt-" + m.id;
    l.innerHTML = `<input type="radio" name="metric" value="${m.id}"> ${m.label}`;
    insideBox.appendChild(l);
  });
  const modeSel = document.getElementById("mode"), timeSel = document.getElementById("time");
  meta.modes.forEach(m => modeSel.add(new Option(m.label, m.id)));
  meta.times.forEach(t => timeSel.add(new Option(t + " min", t)));
  timeSel.value = "10";

  // ---- map ----
  const map = L.map("map", { zoomSnap: 0.25 }).setView([37.758, -122.44], 12.25);
  L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}", {
    attribution: "Tiles &copy; Esri &mdash; Esri, HERE, Garmin, &copy; OpenStreetMap contributors", maxZoom: 16
  }).addTo(map);

  let current = null; // {key, label, unit, digits, breaks}
  const fmt = (v, d) => (v === null || v === undefined || !isFinite(v)) ? "no data"
    : Number(v).toLocaleString(undefined, { minimumFractionDigits: d, maximumFractionDigits: d });

  function quantileBreaks(key) {
    const vals = data.features.map(f => f.properties[key]).filter(v => v !== null && isFinite(v)).sort((a, b) => a - b);
    const q = p => { const i = (vals.length - 1) * p, lo = Math.floor(i), hi = Math.ceil(i);
      return vals[lo] + (vals[hi] - vals[lo]) * (i - lo); };
    return [vals[0], q(0.2), q(0.4), q(0.6), q(0.8), vals[vals.length - 1]];
  }
  function colourFor(v) {
    if (!current || v === null || v === undefined || !isFinite(v)) return NODATA;
    const b = current.breaks;
    for (let i = 1; i < b.length; i++) if (v <= b[i]) return RAMP[i - 1];
    return RAMP[RAMP.length - 1];
  }
  function style(f) {
    const v = current ? f.properties[current.key] : null;
    return { color: "#ffffff", weight: 1.2, fillColor: current ? colourFor(v) : "transparent",
             fillOpacity: current ? 0.82 : 0, opacity: 1 };
  }

  const layer = L.geoJSON(data, {
    style,
    onEachFeature: (f, lyr) => {
      lyr.on("mouseover", e => { e.target.setStyle({ weight: 3, color: getComputedStyle(document.documentElement).getPropertyValue("--hover").trim() }); e.target.bringToFront(); });
      lyr.on("mouseout", e => layer.resetStyle(e.target));
      lyr.bindTooltip(() => tipHtml(f), { sticky: true, className: "nb-tip" });
      lyr.bindPopup(() => popupHtml(f), { maxWidth: 340 });
    }
  }).addTo(map);
  // outline-only mode still needs a visible edge
  const outline = L.geoJSON(data, { style: { color: "#55544f", weight: 1, fill: false }, interactive: false });

  function tipHtml(f) {
    const p = f.properties;
    if (!current) return `<b>${p.neighborhood}</b>`;
    return `<b>${p.neighborhood}</b><br>${current.label}: ${fmt(p[current.key], current.digits)}`;
  }
  function popupHtml(f) {
    const p = f.properties;
    let rows = `<tr><td class="sec" colspan="2">Inside the neighborhood</td></tr>`;
    meta.inside.forEach(m => {
      const cur = current && current.key === m.id ? " class=\\"cur\\"" : "";
      rows += `<tr${cur}><td>${m.label}</td><td class="v">${fmt(p[m.id], m.digits)}</td></tr>`;
    });
    rows += `<tr><td>Area (km²)</td><td class="v">${fmt(p.polygon_area_km2, 2)}</td></tr>`;
    rows += `<tr><td>100 m cells</td><td class="v">${fmt(p.n_cells, 0)}</td></tr>`;
    rows += `<tr><td>GBIF records</td><td class="v">${fmt(p.gbif_records, 0)}</td></tr>`;
    let bai = `<table class="bai"><tr><th>Mean BAI</th>${meta.times.map(t => `<th>${t} min</th>`).join("")}</tr>`;
    meta.modes.forEach(m => {
      bai += `<tr><td>${m.label}</td>` + meta.times.map(t => {
        const k = `bai_${m.id}_${t}min`;
        const b = current && current.key === k ? ` style="font-weight:700"` : "";
        return `<td${b}>${fmt(p[k], 2)}</td>`;
      }).join("") + `</tr>`;
    });
    bai += `</table>`;
    return `<div class="pop"><h3>${p.neighborhood}</h3><table>${rows}</table>
      <div class="sec" style="font-weight:600;color:#8a8983;font-size:10px;text-transform:uppercase;letter-spacing:.04em;padding-top:10px">Access from the neighborhood</div>${bai}</div>`;
  }

  // ---- legend ----
  const legend = L.control({ position: "bottomleft" });
  legend.onAdd = () => { const d = L.DomUtil.create("div", "legend"); d.id = "legend"; return d; };
  legend.addTo(map);
  function drawLegend() {
    const d = document.getElementById("legend");
    if (!current) { d.style.display = "none"; return; }
    d.style.display = "";
    const b = current.breaks;
    let h = `<div class="t">${current.label}</div><div class="u">${current.unit}</div>`;
    for (let i = RAMP.length - 1; i >= 0; i--) {
      h += `<div class="k"><span class="sw" style="background:${RAMP[i]}"></span>${fmt(b[i], current.digits)} – ${fmt(b[i + 1], current.digits)}</div>`;
    }
    const nNA = data.features.filter(f => { const v = f.properties[current.key]; return v === null || !isFinite(v); }).length;
    if (nNA) h += `<div class="k"><span class="sw" style="background:${NODATA}"></span>No data (${nNA})</div>`;
    d.innerHTML = h;
  }

  // ---- selection ----
  function select(value) {
    map.closePopup();
    document.querySelectorAll(".opt").forEach(o => o.classList.remove("on"));
    const radio = document.querySelector(`input[name=metric][value="${value}"]`);
    if (radio) { radio.checked = true; radio.closest(".opt").classList.add("on"); }
    if (value === "__none") {
      current = null;
      outline.addTo(map);
    } else {
      map.removeLayer(outline);
      if (value === "__bai") {
        const key = `bai_${modeSel.value}_${timeSel.value}min`;
        current = { key, label: `Mean BAI · ${modeLabel[modeSel.value]}, ${timeSel.value} min`,
                    unit: "0–1, higher = better access", digits: 2 };
      } else {
        const m = byId[value];
        current = { key: m.id, label: m.label, unit: m.unit, digits: m.digits };
      }
      current.breaks = quantileBreaks(current.key);
    }
    layer.setStyle(style);
    drawLegend();
    try { localStorage.setItem("nb-metric", value); localStorage.setItem("nb-mode", modeSel.value); localStorage.setItem("nb-time", timeSel.value); } catch (e) {}
  }
  document.querySelectorAll("input[name=metric]").forEach(r => r.addEventListener("change", e => select(e.target.value)));
  [modeSel, timeSel].forEach(s => s.addEventListener("change", () => select("__bai")));

  let start = "__bai";
  try {
    const m = localStorage.getItem("nb-metric"), md = localStorage.getItem("nb-mode"), t = localStorage.getItem("nb-time");
    if (md && modeLabel[md]) modeSel.value = md;
    if (t) timeSel.value = t;
    if (m && document.querySelector(`input[name=metric][value="${m}"]`)) start = m;
  } catch (e) {}
  select(start);

  // Fit to the city once the map has a real size (a page opened in a hidden or
  // zero-height frame would otherwise fit to the maximum zoom).
  const bounds = layer.getBounds();
  let fitted = false;
  function fitCity() {
    map.invalidateSize();
    const s = map.getSize();
    if (s.x > 0 && s.y > 0) { map.fitBounds(bounds, { padding: [10, 10], maxZoom: 14 }); fitted = true; }
  }
  fitCity();
  if (window.ResizeObserver) {
    new ResizeObserver(() => { if (!fitted) fitCity(); else map.invalidateSize(); })
      .observe(document.getElementById("map"));
  }
  window.addEventListener("load", () => { if (!fitted) fitCity(); });
})();
</script>
</body>
</html>
'

html <- sub("__GEOJSON__", geojson, html, fixed = TRUE)
html <- sub("__META__", meta_json, html, fixed = TRUE)
stopifnot(!grepl("</script", geojson, ignore.case = TRUE))
writeLines(html, out_html, useBytes = TRUE)
cat(sprintf("Wrote %s (%.0f KB)\n", out_html, file.size(out_html) / 1024))
