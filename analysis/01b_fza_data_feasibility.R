# Stage 1 — Feasibility of L3 satellite products for the Foz do Amazonas case studies
# (decision D2).
#
# For every candidate variable this script records the product, server, processing level,
# period and how it can be accessed. For products that can be read anonymously from ERDDAP
# it downloads the monthly fields over the study rectangle on a 0.25-degree grid and
# measures, month by month, the fraction of cells of the ocean case study (marine FZA
# basin) that hold a valid value. Cloud cover under the ITCZ is the main expected loss.
#
# Requires results/fza/study_area.rds (analysis/01a_fza_study_area.R).
# Raw downloads go to data/raw/fza/erddap/ and are never modified.
# Outputs:
#   - results/fza/candidates.csv          catalogue of candidate products
#   - results/fza/coverage_monthly.csv    coverage per product and month
#   - results/fza/coverage_summary.csv    coverage summary per product
#   - results/figs/fza_coverage.png       coverage by calendar month

suppressMessages({
  library(sf)
  library(ggplot2)
})
sf_use_s2(FALSE)

raw_dir <- file.path("data", "raw", "fza", "erddap")
out_dir <- file.path("results", "fza")
fig_dir <- file.path("results", "figs")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

study <- readRDS(file.path(out_dir, "study_area.rds"))
rect <- study$rectangle
params <- list(
  period = c(start = "2013-01-01", end = "2024-12-31"),  # window inspected; not the final period
  grid_deg = 0.25                                         # common grid of the coverage check
)

# --- Candidate catalogue -------------------------------------------------------------
# erddap_id is filled only for products readable without login. Products in NASA archives
# (PO.DAAC, LP DAAC, NSIDC, GES DISC) need a free NASA Earthdata login, set up by the user.
pfeg <- "https://coastwatch.pfeg.noaa.gov/erddap"
candidates <- read.csv(text = '
case,variable,product,level,native_resolution,access,server,erddap_id,erddap_var,stride,note
ocean,chlorophyll-a,MODIS-Aqua R2022 SQ monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMH1chlamday_R2022SQ,chlor_a,6,standard algorithm biased in turbid plume waters
ocean,Kd490,MODIS-Aqua R2022 SQ monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMH1kd490mday_R2022SQ,Kd_490,6,
ocean,PAR,MODIS-Aqua R2022 SQ monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMH1par0mday_R2022SQ,par,6,
ocean,fluorescence line height,MODIS-Aqua monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMH1cflhmday,cflh,6,ends 2022
ocean,particulate organic carbon,MODIS-Aqua monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMPOCmday,poc,6,ends 2022; derived from the same reflectances as chlorophyll
ocean,particulate inorganic carbon,MODIS-Aqua monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMPICmday,pic,6,ends 2022
ocean,sea surface temperature (day),MODIS-Aqua R2019 SQ monthly,L3 mapped,4 km,anonymous ERDDAP,PFEG,erdMH1sstdmdayR20190SQ,sstMasked,6,ends 2021 on this server
ocean,wind (x and y components),Metop ASCAT monthly,L3 gridded,0.25 deg,anonymous ERDDAP,PFEG,erdQMwindmday_LonPM180,x_wind,1,starts 2013-09; ends 2022 on this server
ocean,sea level anomaly,NOAA CoastWatch multi-mission altimetry daily,L3 declared by provider,0.25 deg,anonymous ERDDAP,PFEG,,,,excluded (L3 rule kept 2026-10-20): gridded multi-mission altimetry is interpolated (L4 in Copernicus terms)
ocean,sea surface salinity,SMAP RSS or JPL monthly,L3 mapped,0.25 deg (about 40 km footprint),NASA Earthdata login,PO.DAAC,,,,shows the Amazon plume; starts 2015-04
ocean,precipitation,GPM IMERG monthly,L3,0.1 deg,NASA Earthdata login,GES DISC,,,,
ocean,CDOM absorption (adg 443),MODIS-Aqua GIOP monthly,L3 mapped,4 km,NASA Earthdata login,OB.DAAC,,,,more reliable than chlorophyll in turbid waters
land,NDVI and EVI,MODIS MOD13A3 monthly,L3,1 km,NASA Earthdata login,LP DAAC,,,,large cloud gaps expected
land,land surface temperature (day and night),MODIS MOD11C3 monthly,L3,0.05 deg,NASA Earthdata login,LP DAAC,,,,
land,albedo,MODIS MCD43C3 daily,L3,0.05 deg,NASA Earthdata login,LP DAAC,,,,
land,soil moisture,SMAP SPL3SMP daily,L3,36 km,NASA Earthdata login,NSIDC,,,,starts 2015-03
land,precipitation,GPM IMERG monthly,L3,0.1 deg,NASA Earthdata login,GES DISC,,,,
land,precipitation (alternative),CHIRPS v2.0 monthly,blended satellite and gauges,0.05 deg,anonymous ERDDAP,PFEG,,,,not L3: merges stations and satellite
land,evapotranspiration,MODIS MOD16A2,L4,500 m,NASA Earthdata login,LP DAAC,,,,L4 (model output); excluded by the L3 rule unless D2 is revisited
land,GPP,MODIS MOD17A2H,L4,500 m,NASA Earthdata login,LP DAAC,,,,L4 (model output); excluded by the L3 rule unless D2 is revisited
', stringsAsFactors = FALSE, na.strings = "")
candidates$stride <- as.integer(candidates$stride)

# --- Ocean mask on the 0.25-degree grid ----------------------------------------------
grid <- expand.grid(lon = seq(rect[["lon_min"]] + params$grid_deg / 2, rect[["lon_max"]], params$grid_deg),
                    lat = seq(rect[["lat_min"]] + params$grid_deg / 2, rect[["lat_max"]], params$grid_deg))
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
grid$ocean_case <- lengths(st_intersects(grid_sf, study$ocean_polygon)) > 0
ocean_cells <- grid[grid$ocean_case, c("lon", "lat")]
ocean_key <- paste(floor(ocean_cells$lon / params$grid_deg), floor(ocean_cells$lat / params$grid_deg))

# --- Download and coverage -----------------------------------------------------------
readable <- candidates[!is.na(candidates$erddap_id), ]
coverage <- list()
for (i in seq_len(nrow(readable))) {
  cand <- readable[i, ]
  info <- read.csv(sprintf("%s/info/%s/index.csv", pfeg, cand$erddap_id), stringsAsFactors = FALSE)
  dims <- info$Variable.Name[info$Row.Type == "dimension"]
  t_range <- info$Value[info$Variable.Name == "time" & info$Attribute.Name == "actual_range"]
  t_range <- as.POSIXct(as.numeric(strsplit(t_range, ",")[[1]]), origin = "1970-01-01", tz = "UTC")
  t0 <- max(as.POSIXct(params$period[["start"]], tz = "UTC"), t_range[1])
  t1 <- min(as.POSIXct(params$period[["end"]], tz = "UTC"), t_range[2])
  file <- file.path(raw_dir, sprintf("%s_%s_%s.csv", cand$erddap_id, format(t0, "%Y%m"), format(t1, "%Y%m")))
  if (!file.exists(file)) {
    # latitude may be stored in descending order: try north-to-south first, then the reverse
    for (lat_order in list(c(rect[["lat_max"]], rect[["lat_min"]]), c(rect[["lat_min"]], rect[["lat_max"]]))) {
      sel <- vapply(dims, function(d) switch(d,
        time = sprintf("[(%s):1:(%s)]", format(t0, "%Y-%m-%dT%H:%M:%SZ"), format(t1, "%Y-%m-%dT%H:%M:%SZ")),
        latitude = sprintf("[(%g):%d:(%g)]", lat_order[1], cand$stride, lat_order[2]),
        longitude = sprintf("[(%g):%d:(%g)]", rect[["lon_min"]], cand$stride, rect[["lon_max"]]),
        "[0]"), "")
      url <- sprintf("%s/griddap/%s.csv?%s%s", pfeg, cand$erddap_id, cand$erddap_var,
                     utils::URLencode(paste(sel, collapse = ""), reserved = TRUE))
      ok <- tryCatch({ download.file(url, file, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (ok) break
      unlink(file)
    }
    if (!file.exists(file)) stop("could not download ", cand$erddap_id)
  }
  d <- read.csv(file, skip = 2, header = FALSE, stringsAsFactors = FALSE)
  names(d)[c(1, ncol(d) - 2, ncol(d) - 1, ncol(d))] <- c("time", "lat", "lon", "value")
  # ERDDAP returns the centre of each native cell picked by the stride; assign it to the
  # 0.25-degree cell that contains it
  key <- paste(floor(d$lon / params$grid_deg), floor(d$lat / params$grid_deg))
  d <- d[key %in% ocean_key, ]
  d$month <- substr(d$time, 1, 7)
  # ERDDAP snaps the time bounds to the nearest composite, which can fall outside the window
  d <- d[d$month >= substr(params$period[["start"]], 1, 7) & d$month <= substr(params$period[["end"]], 1, 7), ]
  cov <- aggregate(value ~ month, data = transform(d, value = !is.na(value)), FUN = mean)
  names(cov)[2] <- "coverage"
  cov$n_cells <- as.vector(table(d$month)[cov$month])
  cov$variable <- cand$variable
  cov$erddap_id <- cand$erddap_id
  coverage[[i]] <- cov
  message(sprintf("%-32s %s to %s, median coverage %.2f", cand$variable, min(cov$month), max(cov$month),
                  median(cov$coverage)))
}
coverage <- do.call(rbind, coverage)
coverage$calendar_month <- as.integer(substr(coverage$month, 6, 7))

summary_tab <- do.call(rbind, lapply(split(coverage, coverage$variable), function(x) data.frame(
  variable = x$variable[1], erddap_id = x$erddap_id[1],
  first_month = min(x$month), last_month = max(x$month), n_months = nrow(x),
  cells_per_month = max(x$n_cells),
  coverage_median = median(x$coverage), coverage_p10 = unname(quantile(x$coverage, 0.1)),
  coverage_min = min(x$coverage),
  wettest_quarter_median = median(x$coverage[x$calendar_month %in% 3:5]),
  driest_quarter_median = median(x$coverage[x$calendar_month %in% 9:11]))))

write.csv(candidates, file.path(out_dir, "candidates.csv"), row.names = FALSE)
write.csv(coverage, file.path(out_dir, "coverage_monthly.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "coverage_summary.csv"), row.names = FALSE)

# --- Figure: coverage by calendar month ----------------------------------------------
fig <- aggregate(coverage ~ variable + calendar_month, data = coverage, FUN = median)
p <- ggplot(fig, aes(factor(calendar_month), reorder(variable, coverage), fill = coverage)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f", 100 * coverage)), size = 3) +
  scale_fill_gradient(low = "#f8ecd9", high = "#0d6b72", limits = c(0, 1), labels = scales::percent) +
  labs(x = "Calendar month", y = NULL, fill = "Valid cells",
       title = "Ocean case study: valid L3 cells in the marine FZA basin",
       subtitle = sprintf("Median over %s to %s; 0.25-degree grid", min(coverage$month), max(coverage$month))) +
  theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "fza_coverage.png"), p, width = 8, height = 4, dpi = 150, bg = "white")

print(summary_tab, digits = 2)
