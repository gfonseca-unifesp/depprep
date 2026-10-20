# Stage 1 — Study area for the Foz do Amazonas (FZA) case studies (decision D2).
#
# Two case studies share one rectangle:
#   - ocean: the marine part of the FZA sedimentary basin (official ANP limit), which
#     covers the shelf, slope and deep water of the continental margin;
#   - land: land cells within 200 km of the stretch of coast that borders that basin.
# The rectangle is the smallest box (rounded outwards to 0.25 degrees) that contains both.
#
# Inputs (downloaded once to data/raw/fza/, never modified):
#   - ANP sedimentary basins (WFS, geomaps.anp.gov.br)
#   - ETOPO1 topography/bathymetry (NOAA CoastWatch ERDDAP, dataset etopo180)
# Outputs:
#   - results/fza/study_area.rds   polygons, rectangle, parameters and sources
#   - results/fza/study_area_bbox.csv
#   - results/figs/fza_study_area.png

suppressMessages({
  library(sf)
  library(ggplot2)
})
sf_use_s2(FALSE)

raw_dir <- file.path("data", "raw", "fza")
out_dir <- file.path("results", "fza")
fig_dir <- file.path("results", "figs")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  inland_distance_km = 200,     # D2: land case study reaches 200 km inland
  coast_touch_km = 10,          # coast cells this close to the marine basin border it
  rectangle_step_deg = 0.25,    # rectangle rounded outwards to this step
  search_box = c(lon_min = -58, lon_max = -40, lat_min = -6, lat_max = 12),
  etopo_stride = 3,             # ETOPO1 is 1 arc-minute; stride 3 = 3 arc-minutes
  crs_metric = 5880             # SIRGAS 2000 / Brazil Polyconic, for distances
)

# --- ANP basin limits ------------------------------------------------------------
anp_file <- file.path(raw_dir, "anp_bacias.json")
anp_url <- paste0("https://geomaps.anp.gov.br/geoserver/wfs?service=wfs&version=1.0.0",
                  "&request=GetFeature&typeName=sisroc:bacias_gishub_db",
                  "&outputFormat=application/json")
if (!file.exists(anp_file)) download.file(anp_url, anp_file, mode = "wb", quiet = TRUE)
basins <- st_read(anp_file, quiet = TRUE)
fza_sea <- st_union(basins[basins$name == "Bacia da Foz do Amazonas_Mar", ])
fza_land_basin <- st_union(basins[basins$name == "Bacia da Foz do Amazonas_Terra", ])
stopifnot(length(fza_sea) == 1, length(fza_land_basin) == 1)

# --- ETOPO1 land/sea mask ----------------------------------------------------------
b <- params$search_box
etopo_file <- file.path(raw_dir, "etopo1_search_box.csv")
etopo_url <- sprintf(paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.csv?",
                            "altitude%%5B(%g):%d:(%g)%%5D%%5B(%g):%d:(%g)%%5D"),
                     b[["lat_min"]], params$etopo_stride, b[["lat_max"]],
                     b[["lon_min"]], params$etopo_stride, b[["lon_max"]])
if (!file.exists(etopo_file)) download.file(etopo_url, etopo_file, mode = "wb", quiet = TRUE)
etopo <- read.csv(etopo_file, skip = 2, header = FALSE,
                  col.names = c("lat", "lon", "altitude"))
etopo_pts <- st_transform(st_as_sf(etopo, coords = c("lon", "lat"), crs = 4326),
                          params$crs_metric)
land <- etopo_pts[etopo_pts$altitude > 0, ]

# --- Coast that borders the marine basin, and land within 200 km of it -------------
fza_sea_m <- st_transform(fza_sea, params$crs_metric)
touch <- st_is_within_distance(land, fza_sea_m, dist = params$coast_touch_km * 1000)
coast <- land[lengths(touch) > 0, ]
near <- st_is_within_distance(land, st_union(coast), dist = params$inland_distance_km * 1000)
land_study <- land[lengths(near) > 0, ]
land_study_ll <- st_transform(land_study, 4326)

# --- Rectangle ---------------------------------------------------------------------
bb <- st_bbox(st_union(st_union(land_study_ll), st_transform(fza_sea, 4326)))
step <- params$rectangle_step_deg
rect <- c(lon_min = floor(bb[["xmin"]] / step) * step, lon_max = ceiling(bb[["xmax"]] / step) * step,
          lat_min = floor(bb[["ymin"]] / step) * step, lat_max = ceiling(bb[["ymax"]] / step) * step)
rect_poly <- st_as_sfc(st_bbox(c(xmin = rect[["lon_min"]], ymin = rect[["lat_min"]],
                                 xmax = rect[["lon_max"]], ymax = rect[["lat_max"]]), crs = 4326))
inside <- function(p) lengths(st_intersects(p, rect_poly)) > 0
stopifnot(all(inside(land_study_ll)))

area_km2 <- c(
  ocean_basin = as.numeric(st_area(fza_sea_m)) / 1e6,
  land_strip = nrow(land_study) * prod(c(111.32 * cos(mean(st_coordinates(land_study_ll)[, 2]) * pi / 180),
                                          111.32) * params$etopo_stride / 60),
  rectangle = as.numeric(st_area(st_transform(rect_poly, params$crs_metric))) / 1e6
)

# --- Save --------------------------------------------------------------------------
study_area <- list(
  rectangle = rect,
  ocean_polygon = st_transform(fza_sea, 4326),
  land_points = land_study_ll,
  coast_points = st_transform(coast, 4326),
  area_km2 = area_km2,
  params = params,
  sources = list(anp = anp_url, etopo = etopo_url),
  downloaded = c(anp = format(file.mtime(anp_file), "%Y-%m-%d"),
                 etopo = format(file.mtime(etopo_file), "%Y-%m-%d"))
)
saveRDS(study_area, file.path(out_dir, "study_area.rds"))
write.csv(data.frame(t(rect)), file.path(out_dir, "study_area_bbox.csv"), row.names = FALSE)

# --- Map ---------------------------------------------------------------------------
crop <- etopo[etopo$lon >= rect[["lon_min"]] - 1 & etopo$lon <= rect[["lon_max"]] + 1 &
              etopo$lat >= rect[["lat_min"]] - 1 & etopo$lat <= rect[["lat_max"]] + 1, ]
crop$surface <- ifelse(crop$altitude > 0, "land", "sea")
p <- ggplot() +
  geom_raster(data = crop, aes(lon, lat, fill = surface)) +
  scale_fill_manual(values = c(land = "grey85", sea = "#dbe9f1"), guide = "none") +
  geom_sf(data = st_transform(land_study, 4326), colour = "#2c7340", size = 0.3, alpha = 0.5) +
  geom_sf(data = st_transform(fza_sea, 4326), fill = NA, colour = "#0d6b72", linewidth = 0.7) +
  geom_sf(data = rect_poly, fill = NA, colour = "#9a5806", linewidth = 0.8, linetype = "dashed") +
  coord_sf(xlim = c(rect[["lon_min"]] - 1, rect[["lon_max"]] + 1),
           ylim = c(rect[["lat_min"]] - 1, rect[["lat_max"]] + 1), expand = FALSE) +
  labs(x = NULL, y = NULL,
       title = "Foz do Amazonas case studies",
       subtitle = sprintf("Blue: marine FZA basin (ANP). Green: land within %d km of its coast. Dashed: study rectangle.",
                          params$inland_distance_km)) +
  theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "fza_study_area.png"), p, width = 7, height = 7, dpi = 150, bg = "white")

print(rect)
print(round(area_km2))
