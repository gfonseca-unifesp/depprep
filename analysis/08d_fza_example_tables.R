# Stage 8 — two example tables of the Foz do Amazonas ocean case, shipped with the package
# (inst/extdata), so that anyone can try the app and the functions on real data.
#
# Why two small gap-free tables. The monthly L3 fields have cloud gaps, and the package
# does not support missing values yet; the space-time table of the real case (stage 9)
# needs decisions on the period, the variables and the gaps. These two tables need none:
# * foz_ocean_map.csv (spatial): one row per 0.25-degree cell of the marine FZA basin, the
#   mean of each variable over its valid months;
# * foz_ocean_series.csv (temporal): one row per month, the basin mean of each variable
#   over its valid cells.
# They are examples, not the stage 9 data: a mean over the valid months leans towards the
# drier, less cloudy months, and a basin mean mixes the plume with the open shelf.
#
# Variables: the anonymous ERDDAP products of analysis/01b_fza_data_feasibility.R with a
# median monthly coverage of at least 0.7 (results/fza/coverage_summary.csv); fluorescence
# line height and particulate inorganic carbon (about 0.6) are left out, and so is Kd490,
# whose monthly composites of 2018-03, 2019-03 and 2020-08 are missing on the server
# (chlorophyll and particulate organic carbon carry close information).
# Period: 2015-04 (the start of the real-case window, set by SMAP salinity) to 2021-10
# (the last month of MODIS SST on this server), 79 months.
# Gap rules: a cell enters the map when every variable is valid in at least half of the
# months; every month must have each variable valid in at least 10% of those cells (the
# script stops otherwise, so the series has no gap). The lowest is particulate organic
# carbon in 2016-02, with 14% of the cells; a first rule of 25% failed on that month.
# Coordinates: the cell centres twice, projected to SIRGAS 2000 / Brazil Mercator
# (EPSG:5641, x_km and y_km) and in degrees (lon and lat); either pair can be used, since
# near the equator (latitudes up to 9 degrees here) a degree of longitude is at most about
# 1% shorter than a degree of latitude.
# No random numbers are used.
#
# Requires results/fza/study_area.rds (analysis/01a_fza_study_area.R) and the raw files of
# analysis/01b_fza_data_feasibility.R in data/raw/fza/erddap/ (read only).
# Outputs (tracked, part of the package): inst/extdata/foz_ocean_map.csv,
#   inst/extdata/foz_ocean_series.csv, inst/extdata/foz_ocean_sources.csv
# Output (not tracked): results/fza/example_tables_summary.csv

suppressMessages(library(sf))
sf_use_s2(FALSE)

raw_dir <- file.path("data", "raw", "fza", "erddap")
out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
grid_deg <- 0.25
window <- c("2015-04", "2021-10")
min_months <- 0.5
min_cells <- 0.10

sources <- read.csv(text = '
variable,name,unit,erddap_id,erddap_var,file
chl,chlorophyll-a,mg m-3,erdMH1chlamday_R2022SQ,chlor_a,erdMH1chlamday_R2022SQ_201301_202412.csv
par,photosynthetically available radiation,einstein m-2 day-1,erdMH1par0mday_R2022SQ,par,erdMH1par0mday_R2022SQ_201301_202309.csv
poc,particulate organic carbon,mg m-3,erdMPOCmday,poc,erdMPOCmday_201301_202205.csv
sst,sea surface temperature (day),degree C,erdMH1sstdmdayR20190SQ,sstMasked,erdMH1sstdmdayR20190SQ_201301_202110.csv
wind_x,eastward wind,m s-1,erdQMwindmday_LonPM180,x_wind,erdQMwindmday_LonPM180_201309_202212.csv
', stringsAsFactors = FALSE)

# --- Cells of the ocean case, as in analysis/01b ----------------------------------------
study <- readRDS(file.path("results", "fza", "study_area.rds"))
rect <- study$rectangle
grid <- expand.grid(lon = seq(rect[["lon_min"]] + grid_deg / 2, rect[["lon_max"]], grid_deg),
                    lat = seq(rect[["lat_min"]] + grid_deg / 2, rect[["lat_max"]], grid_deg))
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
cells <- grid[lengths(st_intersects(grid_sf, study$ocean_polygon)) > 0, ]
cells$key <- paste(floor(cells$lon / grid_deg), floor(cells$lat / grid_deg))
months <- format(seq(as.Date(paste0(window[1], "-01")), as.Date(paste0(window[2], "-01")), by = "month"), "%Y-%m")

# --- One cell x month matrix per variable ------------------------------------------------
values <- lapply(seq_len(nrow(sources)), function(i) {
  d <- read.csv(file.path(raw_dir, sources$file[i]), skip = 2, header = FALSE, stringsAsFactors = FALSE)
  names(d)[c(1, ncol(d) - 2, ncol(d) - 1, ncol(d))] <- c("time", "lat", "lon", "value")
  d$key <- paste(floor(d$lon / grid_deg), floor(d$lat / grid_deg))
  d$month <- substr(d$time, 1, 7)
  d <- d[d$key %in% cells$key & d$month %in% months, ]
  m <- matrix(NA_real_, nrow(cells), length(months))
  m[cbind(match(d$key, cells$key), match(d$month, months))] <- d$value
  m[!is.finite(m)] <- NA
  m
})
names(values) <- sources$variable

# --- The map: cells with every variable valid in at least half of the months -------------
valid_months <- sapply(values, function(m) rowMeans(!is.na(m)))
keep <- apply(valid_months >= min_months, 1, all)
xy <- st_coordinates(st_transform(st_as_sf(cells[keep, ], coords = c("lon", "lat"), crs = 4326), 5641)) / 1000
map <- data.frame(x_km = round(xy[, 1], 2), y_km = round(xy[, 2], 2), lon = cells$lon[keep], lat = cells$lat[keep],
                  sapply(values, function(m) signif(rowMeans(m[keep, ], na.rm = TRUE), 5)))

# --- The series: basin means over the cells of the map ------------------------------------
valid_cells <- sapply(values, function(m) colMeans(!is.na(m[keep, ])))
if (any(valid_cells < min_cells)) {
  bad <- which(valid_cells < min_cells, arr.ind = TRUE)
  stop("months with too few valid cells: ", paste(months[bad[, 1]], colnames(valid_cells)[bad[, 2]], collapse = "; "))
}
series <- data.frame(month = months, time = seq_along(months),
                     sapply(values, function(m) signif(colMeans(m[keep, ], na.rm = TRUE), 5)))
stopifnot(!anyNA(map), !anyNA(series))

write.csv(map, file.path(out_dir, "foz_ocean_map.csv"), row.names = FALSE)
write.csv(series, file.path(out_dir, "foz_ocean_series.csv"), row.names = FALSE)
src <- sources[, c("variable", "name", "unit", "erddap_id")]
src$server <- "https://coastwatch.pfeg.noaa.gov/erddap"
src$provider <- ifelse(src$variable == "wind_x", "NOAA NMFS SWFSC ERD (Metop ASCAT)", "NASA GSFC OBPG (Aqua MODIS)")
src$license <- ifelse(src$variable %in% c("sst", "wind_x"),
                      "free use and redistribution (NOAA ERD); no warranty",
                      "NASA Earth science data and information policy (free use, with attribution)")
write.csv(src, file.path(out_dir, "foz_ocean_sources.csv"), row.names = FALSE)

summary_tab <- data.frame(
  table = c("map", "series"), rows = c(nrow(map), nrow(series)),
  of = c(nrow(cells), length(months)),
  min_share_valid = c(min(valid_months[keep, ]), min(valid_cells)),
  rule = c(sprintf("every variable valid in >= %.2f of the months", min_months),
           sprintf("every variable valid in >= %.2f of the map cells", min_cells)))
dir.create(file.path("results", "fza"), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_tab, file.path("results", "fza", "example_tables_summary.csv"), row.names = FALSE)
print(summary_tab, row.names = FALSE)
print(colMeans(valid_months[keep, ]))
