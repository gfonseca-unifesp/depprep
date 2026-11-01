# Stage 9 — the ocean cube of the Foz do Amazonas case study: every variable of the
# protocol (section 12) on the 0.25-degree cells of the marine FZA basin, month by month,
# April 2015 to December 2022 (93 months), with the gaps left as gaps (L3 rule).
#
# Variables and sources (decided 2026-10-20; access checked in analysis/01b and 01c):
# * MODIS-Aqua L3 mapped monthly 4 km, OB.DAAC (OPeNDAP, Earthdata login): night sea
#   surface temperature (NSST, quality levels 0-1), chlorophyll-a, Kd490, PAR and CDOM
#   absorption (adg 443, GIOP). One archive and one reader for all ocean colour; the
#   ERDDAP copies of chlorophyll, Kd490 and PAR (analysis/08d) miss some months.
# * SMAP sea surface salinity, RSS v6 monthly, PO.DAAC: the 40 km product (sss_smap_40km)
#   rather than the 70 km one, which is smoothed over about three cells of this grid.
# * GPM IMERG V07 monthly precipitation, 0.1 degree, GES DISC (the gpm2 host, which serves
#   the same files as gpm1).
# * Metop ASCAT monthly wind (x and y components), 0.25 degree, NOAA ERDDAP (anonymous).
# * Temporal drivers: Amazon discharge at Obidos (ANA station 17050001, HidroWeb web
#   service, monthly mean, consisted data preferred), ONI (Nino 3.4) and the tropical
#   North and South Atlantic indices (TNA, TSA) from NOAA PSL.
#
# Aggregation. Each native pixel goes to the 0.25-degree cell that contains its centre
# (the grid of analysis/01b; a cell belongs to the basin when its centre is inside the
# ocean polygon). The value of a cell in a month is the mean of its valid pixels. Its
# valid fraction is the number of valid pixels over the number of water pixels of the
# cell (pixels valid in at least one month of the period), so that land does not count
# as a gap. The minimum valid fraction is a stage 9 decision (protocol section 12,
# PENDING): the cube keeps the fraction, and the census of the gaps (analysis/09b)
# reports them at several thresholds.
# No random numbers are used.
#
# Access: the Earthdata login is read by curl from ~/_netrc; this script never prints it.
# Requires results/fza/study_area.rds (analysis/01a_fza_study_area.R).
# Raw subsets go to data/raw/fza/cube/ (one file per product and month, kept as
# downloaded) and data/raw/fza/drivers/; the script resumes from the files present.
# Outputs (not tracked):
#   - results/fza/ocean_cube.csv       cell x month x variable: value, valid fraction
#   - results/fza/ocean_cells.csv      the cells: id, lon, lat, water pixels per product
#   - results/fza/ocean_coverage.csv   per variable and month, share of cells with a value
#                                      at several minimum valid fractions
#   - results/fza/drivers.csv          monthly drivers, with their source and quality

suppressMessages({
  library(sf)
  library(curl)
})
sf_use_s2(FALSE)

raw_dir <- file.path("data", "raw", "fza", "cube")
drv_dir <- file.path("data", "raw", "fza", "drivers")
out_dir <- file.path("results", "fza")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(drv_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  months = seq(as.Date("2015-04-01"), as.Date("2022-12-01"), by = "month"),
  grid_deg = 0.25,
  sst_quality_max = 1,                    # NSST quality levels kept: 0 (best) and 1 (good)
  min_fraction = c(0, 0.25, 0.5, 0.75)    # thresholds reported in the coverage table
)
stopifnot(length(params$months) == 93)

netrc <- file.path(Sys.getenv("USERPROFILE", Sys.getenv("HOME")), "_netrc")
if (!file.exists(netrc)) netrc <- file.path(Sys.getenv("HOME"), ".netrc")
stopifnot(file.exists(netrc))
cookies <- tempfile("earthdata_cookies")
earthdata <- function() new_handle(netrc = 1L, netrc_file = netrc, followlocation = TRUE,
                                   cookiefile = cookies, cookiejar = cookies, timeout = 600)

# --- Cells of the basin --------------------------------------------------------------
study <- readRDS(file.path(out_dir, "study_area.rds"))
rect <- study$rectangle
g <- params$grid_deg
grid <- expand.grid(lon = seq(rect[["lon_min"]] + g / 2, rect[["lon_max"]], g),
                    lat = seq(rect[["lat_min"]] + g / 2, rect[["lat_max"]], g))
inside <- lengths(st_intersects(st_as_sf(grid, coords = c("lon", "lat"), crs = 4326), study$ocean_polygon)) > 0
cells <- grid[inside, ]
cells$key <- paste(floor(cells$lon / g), floor(cells$lat / g))
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 628)
# window of the downloads: the cells of the basin, on whole cells
box <- c(lon_min = min(cells$lon) - g / 2, lon_max = max(cells$lon) + g / 2,
         lat_min = min(cells$lat) - g / 2, lat_max = max(cells$lat) + g / 2)

# --- Readers -------------------------------------------------------------------------
# Hyrax "dap.csv": one line per 1-D variable ("/lat, v1, ...") and one line per row of a
# 2-D variable ("/var[k], v1, ..."), as in analysis/01c.
read_dap_csv <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[grepl("^/", lines)]
  name <- sub("(\\[[0-9]+\\])?,.*$", "", lines)
  vals <- lapply(strsplit(sub("^[^,]*,", "", lines), ","), function(x) as.numeric(trimws(x)))
  out <- list()
  for (n in unique(name)) {
    rows <- vals[name == n]
    out[[n]] <- if (length(rows) == 1 && !grepl("\\[", lines[name == n][1])) rows[[1]] else do.call(rbind, rows)
  }
  out
}
# scale_factor, add_offset and _FillValue of one variable, from the .dmr of a granule
dmr_attributes <- function(url, var) {
  r <- curl_fetch_memory(paste0(url, ".dmr"), handle = earthdata())
  if (r$status_code != 200) stop("dmr of ", basename(url), ": HTTP ", r$status_code)
  txt <- iconv(rawToChar(r$content), "latin1", "UTF-8")
  block <- regmatches(txt, regexpr(sprintf('(?s)<(Int16|Byte|Float32|UInt8)[^>]* name="%s">.*?</(Int16|Byte|Float32|UInt8)>', var), txt, perl = TRUE))
  att <- function(a, default) {
    m <- regmatches(block, regexpr(sprintf('<Attribute name="%s"[^>]*>\\s*<Value>[^<]*', a), block, perl = TRUE))
    if (length(m)) as.numeric(sub(".*<Value>", "", m)) else default
  }
  c(scale = att("scale_factor", 1), offset = att("add_offset", 0), fill = att("_FillValue", NA))
}
# TRUE: downloaded; FALSE: the server has no such file (404); NA: no answer after the
# retries, which must not be taken for a gap in the data
fetch <- function(url, file, handle = earthdata()) {
  if (file.exists(file)) return(TRUE)
  tmp <- paste0(file, ".part")
  for (attempt in 1:4) {                   # archives answer 5xx now and then
    res <- tryCatch(curl_fetch_disk(url, tmp, handle = handle), error = function(e) list(status_code = NA))
    if (identical(res$status_code, 200L) || identical(res$status_code, 200)) { file.rename(tmp, file); return(TRUE) }
    unlink(tmp)
    if (!is.na(res$status_code) && res$status_code == 404) return(FALSE)
    Sys.sleep(10 * attempt)
  }
  NA
}
enc <- function(ce) gsub("\\]", "%5D", gsub("\\[", "%5B", ce))
win <- function(a, b) sprintf("[%d:1:%d]", a, b)

# --- Products ------------------------------------------------------------------------
# idx: rows and columns of the basin window on each native grid. lat_desc: rows run from
# north to south. lon_major: the field is stored [lon][lat].
oc <- list(r0 = floor((90 - box[["lat_max"]]) * 24), r1 = ceiling((90 - box[["lat_min"]]) * 24) - 1,
           c0 = floor((box[["lon_min"]] + 180) * 24), c1 = ceiling((box[["lon_max"]] + 180) * 24) - 1)
smap <- list(r0 = floor((box[["lat_min"]] + 90) / 0.25), r1 = ceiling((box[["lat_max"]] + 90) / 0.25) - 1,
             c0 = floor((box[["lon_min"]] + 360) / 0.25), c1 = ceiling((box[["lon_max"]] + 360) / 0.25) - 1)
imerg <- list(r0 = floor((box[["lat_min"]] + 90) / 0.1), r1 = ceiling((box[["lat_max"]] + 90) / 0.1) - 1,
              c0 = floor((box[["lon_min"]] + 180) / 0.1), c1 = ceiling((box[["lon_max"]] + 180) / 0.1) - 1)
obdaac <- function(suite, var) function(m) {
  last <- seq(m, by = "month", length.out = 2)[2] - 1
  sprintf("https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/%s/%s/AQUA_MODIS.%s_%s.L3m.MO.%s.%s.4km.nc",
          format(m, "%Y"), format(m, "%m%d"), format(m, "%Y%m%d"), format(last, "%Y%m%d"), suite, var)
}
products <- list(
  list(key = "sst_night", var = "sst", url = obdaac("NSST", "sst"), idx = oc, quality = "qual_sst"),
  list(key = "chlorophyll", var = "chlor_a", url = obdaac("CHL", "chlor_a"), idx = oc),
  list(key = "kd490", var = "Kd_490", url = obdaac("KD", "Kd_490"), idx = oc),
  list(key = "par", var = "par", url = obdaac("PAR", "par"), idx = oc),
  list(key = "adg443", var = "adg_443", url = obdaac("IOP", "adg_443"), idx = oc),
  list(key = "salinity", var = "sss_smap_40km", idx = smap, lat_asc = TRUE, lon360 = TRUE,
       url = function(m) sprintf("https://opendap.earthdata.nasa.gov/collections/C2832226365-POCLOUD/granules/RSS_smap_SSS_L3_monthly_%s_FNL_v06.0",
                                 format(m, "%Y_%m"))),
  list(key = "precipitation", var = "precipitation", idx = imerg, lat_asc = TRUE, lon_major = TRUE,
       url = function(m) sprintf("https://gpm2.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGM.07/%s/3B-MO.MS.MRG.3IMERG.%s-S000000-E235959.%s.V07B.HDF5",
                                 format(m, "%Y"), format(m, "%Y%m%d"), format(m, "%m")))
)

# one product and month: download the window, return the pixels (lon, lat, value)
pixels_dap <- function(p, m, att) {
  url <- p$url(m)
  file <- file.path(raw_dir, p$key, paste0(format(m, "%Y-%m"), ".csv"))
  i <- p$idx
  field <- if (isTRUE(p$lon_major)) sprintf("/%s[0]%s%s", p$var, win(i$c0, i$c1), win(i$r0, i$r1)) else
    if (p$key == "salinity") sprintf("/%s%s%s", p$var, win(i$r0, i$r1), win(i$c0, i$c1)) else
      sprintf("/%s%s%s", p$var, win(i$r0, i$r1), win(i$c0, i$c1))
  ce <- paste(c(sprintf("/lat%s", win(i$r0, i$r1)), sprintf("/lon%s", win(i$c0, i$c1)), field,
                if (!is.null(p$quality)) sprintf("/%s%s%s", p$quality, win(i$r0, i$r1), win(i$c0, i$c1))), collapse = ";")
  got <- fetch(sprintf("%s.dap.csv?dap4.ce=%s", url, enc(ce)), file)
  if (is.na(got)) stop(sprintf("%s %s: the server did not answer; run the script again later", p$key, format(m, "%Y-%m")), call. = FALSE)
  if (!got) return(NULL)
  d <- read_dap_csv(file)
  val <- d[[sub("\\[.*$", "", sub("^/", "/", field))]]
  if (is.null(val)) val <- d[[grep(p$var, names(d), value = TRUE)[1]]]
  # cells must follow expand.grid(lon, lat): lon varies fastest
  x <- if (isTRUE(p$lon_major)) as.vector(val) else as.vector(t(val))
  bad <- is.na(x) | (!is.na(att[["fill"]]) & x == att[["fill"]])
  x <- x * att[["scale"]] + att[["offset"]]
  if (!is.null(p$quality)) bad <- bad | is.na(as.vector(t(d[[paste0("/", p$quality)]]))) |
    as.vector(t(d[[paste0("/", p$quality)]])) > params$sst_quality_max
  if (p$key == "salinity") bad <- bad | x < 0 | x > 45
  if (p$key == "precipitation") bad <- bad | x < 0
  x[bad] <- NA
  lon <- d[["/lon"]]; lat <- d[["/lat"]]
  if (isTRUE(p$lon360)) lon <- ifelse(lon > 180, lon - 360, lon)
  g2 <- expand.grid(lon = lon, lat = lat)
  data.frame(lon = g2$lon, lat = g2$lat, value = x)
}

# per cell and month: sum and number of valid pixels; per pixel: valid in some month
aggregate_month <- function(px) {
  key <- paste(floor(px$lon / g), floor(px$lat / g))
  cell <- cells$cell[match(key, cells$key)]
  keep <- !is.na(cell)
  px <- px[keep, ]; cell <- cell[keep]
  ok <- !is.na(px$value)
  s <- as.vector(tapply(ifelse(ok, px$value, 0), factor(cell, levels = cells$cell), sum))
  n <- as.vector(tapply(ok, factor(cell, levels = cells$cell), sum))
  s[is.na(s)] <- 0; n[is.na(n)] <- 0                      # cells with no pixel of this grid
  # a pixel is named by its centre on a 1/480-degree lattice: servers print the same
  # centre with different digits from one month to the next
  list(pixel = paste(round(px$lon * 480), round(px$lat * 480)), cell = cell, valid = ok, sum = s, n = n)
}

cube <- list(); water <- list()
for (p in products) {
  dir.create(file.path(raw_dir, p$key), showWarnings = FALSE)
  att <- NULL; months_ok <- list(); ever <- NULL; cell_of <- NULL
  for (k in seq_along(params$months)) {
    m <- params$months[k]
    if (is.null(att)) att <- if (p$key %in% c("salinity", "precipitation")) c(scale = 1, offset = 0, fill = -9999) else
      dmr_attributes(p$url(m), p$var)
    px <- pixels_dap(p, m, att)
    if (is.null(px)) { message(sprintf("%-13s %s missing on the server", p$key, format(m, "%Y-%m"))); next }
    a <- aggregate_month(px)
    if (is.null(ever)) { ever <- setNames(a$valid, a$pixel); cell_of <- setNames(a$cell, a$pixel) } else {
      new <- setdiff(a$pixel, names(ever))
      if (length(new)) { ever[new] <- FALSE; cell_of[new] <- a$cell[match(new, a$pixel)] }
      ever[a$pixel] <- ever[a$pixel] | a$valid
    }
    months_ok[[format(m, "%Y-%m")]] <- a
    if (k %% 12 == 0) message(sprintf("%-13s %s done", p$key, format(m, "%Y-%m")))
  }
  n_water <- as.vector(tapply(ever, factor(cell_of, levels = cells$cell), sum)); n_water[is.na(n_water)] <- 0
  water[[p$key]] <- as.vector(n_water)
  cube[[p$key]] <- do.call(rbind, lapply(names(months_ok), function(mm) {
    a <- months_ok[[mm]]
    data.frame(cell = cells$cell, month = mm, variable = p$key,
               value = ifelse(a$n > 0, a$sum / pmax(a$n, 1), NA), valid_pixels = as.vector(a$n),
               valid_fraction = ifelse(n_water > 0, as.vector(a$n) / as.vector(n_water), 0))
  }))
  message(sprintf("%-13s %d of %d months", p$key, length(months_ok), length(params$months)))
}

# --- Wind (ERDDAP, 0.25 degree, one request for the period) --------------------------
pfeg <- "https://coastwatch.pfeg.noaa.gov/erddap"
wind_file <- file.path(raw_dir, "wind", "erdQMwindmday_LonPM180_201504_202212.csv")
dir.create(dirname(wind_file), showWarnings = FALSE)
if (!file.exists(wind_file)) {
  sel <- function(lat_order) sprintf("[(2015-04-01T00:00:00Z):1:(2022-12-31T00:00:00Z)][0][(%g):1:(%g)][(%g):1:(%g)]",
                                     lat_order[1], lat_order[2], box[["lon_min"]], box[["lon_max"]])
  for (lat_order in list(c(box[["lat_min"]], box[["lat_max"]]), c(box[["lat_max"]], box[["lat_min"]]))) {
    s <- utils::URLencode(sel(lat_order), reserved = TRUE)
    if (isTRUE(fetch(sprintf("%s/griddap/erdQMwindmday_LonPM180.csv?x_wind%s,y_wind%s", pfeg, s, s), wind_file, new_handle(timeout = 600)))) break
  }
  stopifnot(file.exists(wind_file))
}
w <- read.csv(wind_file, skip = 2, header = FALSE, stringsAsFactors = FALSE)
names(w) <- c("time", "altitude", "lat", "lon", "x_wind", "y_wind")
w$month <- substr(w$time, 1, 7)
w <- w[w$month %in% format(params$months, "%Y-%m"), ]
for (v in c("x_wind", "y_wind")) {
  key <- c(x_wind = "wind_u", y_wind = "wind_v")[[v]]
  by_month <- lapply(split(w, w$month), function(d) aggregate_month(data.frame(lon = d$lon, lat = d$lat, value = d[[v]])))
  ever <- Reduce(`|`, lapply(by_month, `[[`, "valid"))
  n_water <- as.vector(tapply(ever, factor(by_month[[1]]$cell, levels = cells$cell), sum)); n_water[is.na(n_water)] <- 0
  water[[key]] <- as.vector(n_water)
  cube[[key]] <- do.call(rbind, lapply(names(by_month), function(mm) { a <- by_month[[mm]]
    data.frame(cell = cells$cell, month = mm, variable = key, value = ifelse(a$n > 0, a$sum / pmax(a$n, 1), NA),
               valid_pixels = as.vector(a$n), valid_fraction = ifelse(n_water > 0, as.vector(a$n) / as.vector(n_water), 0)) }))
}
unlink(cookies)

cube <- do.call(rbind, cube)
# every variable and month present, so that a missing composite shows as a gap
full <- expand.grid(cell = cells$cell, month = format(params$months, "%Y-%m"), variable = unique(cube$variable),
                    stringsAsFactors = FALSE)
cube <- merge(full, cube, all.x = TRUE)
cube$valid_pixels[is.na(cube$valid_pixels)] <- 0; cube$valid_fraction[is.na(cube$valid_fraction)] <- 0
cube <- cube[order(cube$variable, cube$month, cube$cell), ]
cells_out <- cbind(cells[, c("cell", "lon", "lat")], setNames(as.data.frame(water), paste0("water_pixels_", names(water))))
write.csv(cube, file.path(out_dir, "ocean_cube.csv"), row.names = FALSE)
write.csv(cells_out, file.path(out_dir, "ocean_cells.csv"), row.names = FALSE)

coverage <- do.call(rbind, lapply(params$min_fraction, function(f) {
  ok <- !is.na(cube$value) & cube$valid_fraction >= f & (f > 0 | cube$valid_pixels > 0)
  d <- aggregate(ok ~ variable + month, data = data.frame(variable = cube$variable, month = cube$month, ok = ok), FUN = mean)
  names(d)[3] <- "share_cells"; d$min_fraction <- f; d
}))
write.csv(coverage, file.path(out_dir, "ocean_coverage.csv"), row.names = FALSE)

# --- Drivers ---------------------------------------------------------------------------
# Obidos: the web service returns one record per month and consistency level (1 raw,
# 2 consisted); consisted data are kept where they exist.
ana <- do.call(rbind, lapply(2015:2022, function(y) {
  file <- file.path(drv_dir, sprintf("ana_17050001_%d.xml", y))
  url <- sprintf("http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroSerieHistorica?codEstacao=17050001&dataInicio=01/01/%d&dataFim=31/12/%d&tipoDados=3&nivelConsistencia=", y, y)
  stopifnot(fetch(url, file, new_handle(timeout = 300)))
  s <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "")
  recs <- regmatches(s, gregexpr("<SerieHistorica[^>]*>.*?</SerieHistorica>", s, perl = TRUE))[[1]]
  tag <- function(r, t) { m <- regmatches(r, regexpr(sprintf("<%s>[^<]*</%s>", t, t), r)); if (length(m)) sub("<[^>]*>([^<]*)<.*", "\\1", m) else NA }
  if (!length(recs)) return(NULL)
  data.frame(month = substr(vapply(recs, tag, "", "DataHora"), 1, 7), level = as.integer(vapply(recs, tag, "", "NivelConsistencia")),
             discharge = as.numeric(vapply(recs, tag, "", "Media")), row.names = NULL)
}))
ana <- ana[order(ana$month, -ana$level), ]; ana <- ana[!duplicated(ana$month) & !is.na(ana$discharge), ]
psl <- function(name) {
  file <- file.path(drv_dir, paste0(name, ".data"))
  stopifnot(fetch(sprintf("https://psl.noaa.gov/data/correlation/%s.data", name), file, new_handle(timeout = 120)))
  l <- strsplit(trimws(readLines(file, warn = FALSE)), "[[:space:]]+")
  l <- Filter(function(x) length(x) == 13 && grepl("^(19|20)[0-9]{2}$", x[1]), l)
  d <- data.frame(month = as.vector(t(outer(vapply(l, `[`, "", 1), sprintf("%02d", 1:12), paste, sep = "-"))),
                  value = as.numeric(unlist(lapply(l, `[`, -1))))
  d$value[d$value <= -99] <- NA
  d
}
mm <- format(params$months, "%Y-%m")
drivers <- data.frame(month = mm,
  discharge_obidos = ana$discharge[match(mm, ana$month)], discharge_level = ana$level[match(mm, ana$month)],
  oni = psl("oni")$value[match(mm, psl("oni")$month)],
  tna = psl("tna")$value[match(mm, psl("tna")$month)],
  tsa = psl("tsa")$value[match(mm, psl("tsa")$month)])
write.csv(drivers, file.path(out_dir, "drivers.csv"), row.names = FALSE)

# --- Summary ---------------------------------------------------------------------------
cs <- coverage[coverage$min_fraction == 0, ]
print(do.call(rbind, lapply(split(cs, cs$variable), function(x) data.frame(variable = x$variable[1],
  months = sum(x$share_cells > 0), median = median(x$share_cells), p10 = unname(quantile(x$share_cells, 0.1)), min = min(x$share_cells)))),
  digits = 2, row.names = FALSE)
cat("drivers: missing months", colSums(is.na(drivers[, -1])), "\n")
