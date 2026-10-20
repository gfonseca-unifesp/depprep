# Stage 1 — Feasibility of NASA Earthdata L3 products for the Foz do Amazonas case studies
# (decision D2). Companion of 01b_fza_data_feasibility.R, which covers the products that
# ERDDAP serves without login.
#
# Products are subset to the study rectangle with OPeNDAP (no global downloads). Coverage
# is measured on a sample of months (params$sample_months) to see the seasonal cycle of
# cloud and swath gaps; the final period is chosen later.
#
# Access: needs a NASA Earthdata login stored by the maintainer in ~/_netrc
# ("machine urs.earthdata.nasa.gov login ... password ..."). This script never prints it.
#
# Requires results/fza/study_area.rds (analysis/01a_fza_study_area.R).
# Raw subsets go to data/raw/fza/nasa/ and are never modified.
# Outputs:
#   - results/fza/nasa_access.csv             HTTP status per archive
#   - results/fza/nasa_coverage_monthly.csv   valid-cell fraction per product and month
#   - results/fza/nasa_coverage_summary.csv
#   - results/fza/mod13c2_reliability.csv     share of each MOD13C2 reliability rank
#   - results/figs/fza_nasa_coverage.png

suppressMessages({
  library(sf)
  library(curl)
  library(ggplot2)
})
sf_use_s2(FALSE)

raw_dir <- file.path("data", "raw", "fza", "nasa")
out_dir <- file.path("results", "fza")
fig_dir <- file.path("results", "figs")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

netrc <- file.path(Sys.getenv("USERPROFILE", Sys.getenv("HOME")), "_netrc")
if (!file.exists(netrc)) netrc <- file.path(Sys.getenv("HOME"), ".netrc")
stopifnot(file.exists(netrc))
cookies <- tempfile("earthdata_cookies")
handle <- new_handle(netrc = 1L, netrc_file = netrc, followlocation = TRUE,
                     cookiefile = cookies, cookiejar = cookies)

study <- readRDS(file.path(out_dir, "study_area.rds"))
rect <- study$rectangle
params <- list(
  sample_months = format(seq(as.Date("2018-01-01"), as.Date("2019-12-01"), by = "month"), "%Y-%m"),
  smap_sm_months = c("2018-03", "2018-09"),   # daily product: two months, wet and dry season
  land_match_km = 5                           # cell centre this close to a land-strip point
)
cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json"

# --- Product table -------------------------------------------------------------------
# idx: index windows of the study rectangle on each native grid (checked against the
# coordinates returned by the server below).
idx_cmg <- c(r0 = floor((90 - rect[["lat_max"]]) / 0.05), r1 = ceiling((90 - rect[["lat_min"]]) / 0.05) - 1,
             c0 = floor((rect[["lon_min"]] + 180) / 0.05), c1 = ceiling((rect[["lon_max"]] + 180) / 0.05) - 1)
idx_smap_sss <- c(r0 = floor((rect[["lat_min"]] + 90) / 0.25), r1 = ceiling((rect[["lat_max"]] + 90) / 0.25) - 1,
                  c0 = floor((rect[["lon_min"]] + 360) / 0.25), c1 = ceiling((rect[["lon_max"]] + 360) / 0.25) - 1)
idx_oc4 <- c(r0 = floor((90 - rect[["lat_max"]]) * 24), r1 = ceiling((90 - rect[["lat_min"]]) * 24) - 1,
             c0 = floor((rect[["lon_min"]] + 180) * 24), c1 = ceiling((rect[["lon_max"]] + 180) * 24) - 1)
# GPM IMERG 0.1 degree: precipitation is stored [time][lon][lat], lat ascending from -90
idx_imerg <- c(r0 = floor((rect[["lat_min"]] + 90) / 0.1), r1 = ceiling((rect[["lat_max"]] + 90) / 0.1) - 1,
               c0 = floor((rect[["lon_min"]] + 180) / 0.1), c1 = ceiling((rect[["lon_max"]] + 180) / 0.1) - 1)
# EASE-Grid 2.0 global, 36 km (EPSG:6933): corner and cell size from the NSIDC grid definition
ease <- list(x0 = -17367530.445161499, y0 = 7314540.830638, cell = 36032.220840584)
corners <- unname(st_coordinates(st_transform(st_as_sf(data.frame(
  lon = c(rect[["lon_min"]], rect[["lon_max"]]), lat = c(rect[["lat_max"]], rect[["lat_min"]])),
  coords = c("lon", "lat"), crs = 4326), 6933)))
idx_ease <- c(r0 = floor((ease$y0 - corners[1, 2]) / ease$cell), r1 = floor((ease$y0 - corners[2, 2]) / ease$cell),
              c0 = floor((corners[1, 1] - ease$x0) / ease$cell), c1 = floor((corners[2, 1] - ease$x0) / ease$cell))
w <- function(i, s = 1) sprintf("[%d:%d:%d]", i[["r0"]], s, i[["r1"]])  # row window
v <- function(i, s = 1) sprintf("[%d:%d:%d]", i[["c0"]], s, i[["c1"]])  # column window

# Products with a fixed file-naming scheme are addressed directly (fewer catalogue calls)
imerg_url <- function(first, last) sprintf(
  "https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGM.07/%s/3B-MO.MS.MRG.3IMERG.%s-S000000-E235959.%s.V07B.HDF5",
  format(first, "%Y"), format(first, "%Y%m%d"), format(first, "%m"))

products <- list(
  list(key = "smap_sss", case = "ocean", variable = "sea surface salinity", archive = "PO.DAAC",
       collection = "C2832226365-POCLOUD", pattern = "monthly_%Y_%m_",
       ce = sprintf("/lat%s;/lon%s;/sss_smap%s%s", w(idx_smap_sss), v(idx_smap_sss), w(idx_smap_sss), v(idx_smap_sss)),
       field = "/sss_smap", valid = function(x) !is.na(x) & x > -9999),
  list(key = "modis_adg443", case = "ocean", variable = "CDOM absorption (adg 443)", archive = "OB.DAAC",
       collection = "C3380709177-OB_CLOUD", pattern = "%Y%m01_.*\\.MO\\.IOP\\.adg_443\\.4km",
       url = function(first, last) sprintf("https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/%s/%s/AQUA_MODIS.%s_%s.L3m.MO.IOP.adg_443.4km.nc",
                                           format(first, "%Y"), format(first, "%m%d"), format(first, "%Y%m%d"), format(last, "%Y%m%d")),
       ce = sprintf("/lat%s;/lon%s;/adg_443%s%s", w(idx_oc4, 6), v(idx_oc4, 6), w(idx_oc4, 6), v(idx_oc4, 6)),
       field = "/adg_443", valid = function(x) !is.na(x) & x != -32767),
  list(key = "mod13c2", case = "land", variable = "NDVI and EVI (reliability 0-1)", archive = "LP DAAC",
       collection = "C2565788914-LPCLOUD", pattern = "MOD13C2\\.A%Y%j\\.",
       ce = sprintf("/MOD_Grid_monthly_CMG_VI/Latitude%s;/MOD_Grid_monthly_CMG_VI/Longitude%s;/MOD_Grid_monthly_CMG_VI/Data_Fields/CMG_0_05_Deg_Monthly_pixel_reliability%s%s",
                    w(idx_cmg), v(idx_cmg), w(idx_cmg), v(idx_cmg)),
       field = "/MOD_Grid_monthly_CMG_VI/Data_Fields/CMG_0_05_Deg_Monthly_pixel_reliability",
       valid = function(x) !is.na(x) & x %in% c(0, 1)),
  list(key = "mod11c3_day", case = "land", variable = "land surface temperature (day)", archive = "LP DAAC",
       collection = "C2565788897-LPCLOUD", pattern = "MOD11C3\\.A%Y%j\\.",
       ce = sprintf("/MODIS_MONTHLY_0_05DEG_CMG_LST/Latitude%s;/MODIS_MONTHLY_0_05DEG_CMG_LST/Longitude%s;/MODIS_MONTHLY_0_05DEG_CMG_LST/Data_Fields/LST_Day_CMG%s%s;/MODIS_MONTHLY_0_05DEG_CMG_LST/Data_Fields/LST_Night_CMG%s%s",
                    w(idx_cmg), v(idx_cmg), w(idx_cmg), v(idx_cmg), w(idx_cmg), v(idx_cmg)),
       field = "/MODIS_MONTHLY_0_05DEG_CMG_LST/Data_Fields/LST_Day_CMG",
       field2 = c(key = "mod11c3_night", variable = "land surface temperature (night)",
                  field = "/MODIS_MONTHLY_0_05DEG_CMG_LST/Data_Fields/LST_Night_CMG"),
       valid = function(x) !is.na(x) & x > 0),
  list(key = "imerg_ocean", case = "ocean", variable = "precipitation (GPM IMERG)", archive = "GES DISC",
       collection = "C2723754851-GES_DISC", pattern = "IMERG\\.%Y%m01-", lon_major = TRUE, url = imerg_url,
       ce = sprintf("/lat%s;/lon%s;/precipitation[0]%s%s",
                    w(idx_imerg), v(idx_imerg), v(idx_imerg), w(idx_imerg)),
       field = "/precipitation[0]", valid = function(x) !is.na(x) & x >= 0),
  list(key = "imerg_land", case = "land", variable = "precipitation (GPM IMERG)", archive = "GES DISC",
       collection = "C2723754851-GES_DISC", pattern = "IMERG\\.%Y%m01-", lon_major = TRUE, url = imerg_url,
       ce = sprintf("/lat%s;/lon%s;/precipitation[0]%s%s",
                    w(idx_imerg), v(idx_imerg), v(idx_imerg), w(idx_imerg)),
       field = "/precipitation[0]", valid = function(x) !is.na(x) & x >= 0),
  list(key = "smap_sm", case = "land", variable = "soil moisture (SMAP, AM)", archive = "NSIDC",
       collection = "C2938664585-NSIDC_CPRD", pattern = "_%Y%m%d_", daily = TRUE,
       ce = sprintf("/Soil_Moisture_Retrieval_Data_AM/latitude%s%s;/Soil_Moisture_Retrieval_Data_AM/soil_moisture%s%s",
                    w(idx_ease), v(idx_ease), w(idx_ease), v(idx_ease)),
       field = "/Soil_Moisture_Retrieval_Data_AM/soil_moisture", valid = function(x) !is.na(x) & x > -9999)
)

# --- Masks ---------------------------------------------------------------------------
land_pts <- st_transform(study$land_points, 5880)
ocean_poly <- study$ocean_polygon
in_case <- function(lon, lat, case) {
  p <- st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  if (case == "ocean") return(lengths(st_intersects(p, ocean_poly)) > 0)
  lengths(st_is_within_distance(st_transform(p, 5880), land_pts, dist = params$land_match_km * 1000)) > 0
}

# --- DAP4 CSV reader -------------------------------------------------------------------
# Hyrax "dap.csv": one line per 1-D variable ("/lat, v1, v2, ...") and one line per row of
# a 2-D variable ("/var[k], v1, v2, ...").
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

# --- Granule lists from CMR (public) -------------------------------------------------
cmr_granules <- function(collection, start, end, name = NULL) {
  res <- list(); page <- 1
  repeat {
    url <- sprintf("%s?collection_concept_id=%s&temporal=%s,%s&page_size=2000&page_num=%d",
                   cmr, collection, start, end, page)
    if (!is.null(name)) url <- paste0(url, "&readable_granule_name=", URLencode(name, reserved = TRUE),
                                      "&options%5Breadable_granule_name%5D%5Bpattern%5D=true")
    r <- NULL
    for (attempt in 1:4) {  # CMR answers 5xx now and then
      r <- curl_fetch_memory(url)
      if (r$status_code == 200) break
      Sys.sleep(5 * attempt)
    }
    if (r$status_code != 200) stop("CMR failed with HTTP ", r$status_code)
    entries <- jsonlite::fromJSON(rawToChar(r$content), simplifyVector = FALSE)$feed$entry
    if (!length(entries)) break
    res <- c(res, entries); if (length(entries) < 2000) break; page <- page + 1
  }
  data.frame(title = vapply(res, `[[`, "", "title"),
             opendap = vapply(res, function(e) {
               h <- vapply(e$links, `[[`, "", "href")
               o <- h[grepl("opendap", h)]
               if (length(o)) o[1] else NA_character_
             }, ""), stringsAsFactors = FALSE)
}

# --- Access check per archive (GES DISC included) -------------------------------------
access_urls <- c(
  "PO.DAAC" = "https://opendap.earthdata.nasa.gov/collections/C2832226365-POCLOUD/granules/RSS_smap_SSS_L3_monthly_2018_06_FNL_v06.0.dmr",
  "LP DAAC" = "https://opendap.earthdata.nasa.gov/collections/C2565788914-LPCLOUD/granules/MOD13C2.A2018152.061.2021349073614.dmr",
  "NSIDC" = "https://opendap.earthdata.nasa.gov/collections/C2938664585-NSIDC_CPRD/granules/SMAP_L3_SM_P_20180601_R19240_001.h5.dmr",
  "OB.DAAC" = "https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2018/0601/AQUA_MODIS.20180601_20180630.L3m.MO.IOP.adg_443.4km.nc.dmr",
  "GES DISC" = "https://data.gesdisc.earthdata.nasa.gov/data/GPM_L3/GPM_3IMERGM.07/2018/3B-MO.MS.MRG.3IMERG.20180601-S000000-E235959.06.V07B.HDF5"
)
access <- data.frame(archive = names(access_urls), product = c("SMAP SSS", "MOD13C2 / MOD11C3", "SMAP SPL3SMP",
                     "MODIS-Aqua IOP", "GPM IMERG"), http_status = NA_integer_, stringsAsFactors = FALSE)
for (i in seq_along(access_urls)) {
  # .dmr documents are small and generated on request (no byte ranges); data files are
  # probed with the first kilobyte only
  probe <- new_handle(netrc = 1L, netrc_file = netrc, followlocation = TRUE,
                      cookiefile = cookies, cookiejar = cookies)
  if (!grepl("\\.dmr$", access_urls[[i]])) handle_setopt(probe, range = "0-1023")
  for (attempt in 1:4) {  # 5xx answers are transient server errors, not refusals
    access$http_status[i] <- curl_fetch_memory(access_urls[[i]], handle = probe)$status_code
    if (access$http_status[i] < 500) break
    Sys.sleep(10 * attempt)
  }
}
access$ok <- access$http_status %in% c(200, 206)
access$note <- ifelse(access$ok, "access granted",
                      ifelse(access$http_status %in% c(401, 403),
                             "refused: authorise this archive's application in the Earthdata profile",
                             "server error at run time (not an access refusal); cached subsets used if present"))
access$checked <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
write.csv(access, file.path(out_dir, "nasa_access.csv"), row.names = FALSE)
print(access)

# --- Download subsets and measure coverage ---------------------------------------------
coverage <- list(); reliability <- list()
for (p in products) {
  can_download <- access$ok[match(p$archive, access$archive)]
  if (!can_download && !dir.exists(file.path(raw_dir, p$key))) { message("skip ", p$key, ": no access"); next }
  if (!can_download) message(p$key, ": no access now; using subsets already downloaded")
  months <- if (isTRUE(p$daily)) params$smap_sm_months else params$sample_months
  gran <- tryCatch(do.call(rbind, lapply(months, function(m) {
    first <- as.Date(paste0(m, "-01"))
    last <- seq(first, by = "month", length.out = 2)[2] - 1
    if (!is.null(p$url)) {
      u <- p$url(first, last)
      return(data.frame(title = basename(u), opendap = u, month = m, stringsAsFactors = FALSE))
    }
    g <- cmr_granules(p$collection, paste0(first, "T00:00:00Z"), paste0(last, "T23:59:59Z"), p$cmr_name)
    days <- if (isTRUE(p$daily)) seq(first, last, by = "day") else first
    keep <- unlist(lapply(days, function(d) which(grepl(format(d, p$pattern), g$title))))
    g <- g[unique(keep), , drop = FALSE]
    if (nrow(g)) g$month <- m
    g
  })), error = function(e) {
    # catalogue unavailable: rebuild the list from the subsets already downloaded, reading
    # the date in each file name (AYYYYDDD, YYYY_MM or YYYYMMDD)
    files <- list.files(file.path(raw_dir, p$key), pattern = "\\.csv$")
    if (!length(files)) stop(conditionMessage(e))
    message(p$key, ": catalogue unavailable (", conditionMessage(e), "); using local file names")
    date <- as.Date(rep(NA_character_, length(files)))
    i <- grepl("A[0-9]{7}", files)
    if (any(i)) {
      doy <- regmatches(files[i], regexpr("A[0-9]{7}", files[i]))
      date[i] <- as.Date(paste0(substr(doy, 2, 5), "-01-01")) + as.integer(substr(doy, 6, 8)) - 1
    }
    i <- is.na(date) & grepl("_(19|20)[0-9]{6}_", files)
    if (any(i)) date[i] <- as.Date(substr(regmatches(files[i], regexpr("_(19|20)[0-9]{6}_", files[i])), 2, 9),
                                   format = "%Y%m%d")
    i <- is.na(date) & grepl("_(19|20)[0-9]{2}_[0-9]{2}_", files)
    if (any(i)) {
      ym <- regmatches(files[i], regexpr("(19|20)[0-9]{2}_[0-9]{2}", files[i]))
      date[i] <- as.Date(paste0(sub("_", "-", ym), "-01"))
    }
    out <- data.frame(title = sub("\\.csv$", "", files), opendap = NA_character_,
                      month = format(date, "%Y-%m"), stringsAsFactors = FALSE)
    out[out$month %in% months, ]
  })
  dir.create(file.path(raw_dir, p$key), showWarnings = FALSE)
  mask <- NULL
  for (j in seq_len(nrow(gran))) {
    base <- if (!is.na(gran$opendap[j])) gran$opendap[j] else
      sprintf("https://opendap.earthdata.nasa.gov/collections/%s/granules/%s", p$collection, gran$title[j])
    base <- sub("\\.html$", "", base)
    file <- file.path(raw_dir, p$key, paste0(gsub("[^A-Za-z0-9._-]", "_", gran$title[j]), ".csv"))
    if (!file.exists(file) && !can_download) next
    if (!file.exists(file)) {
      ce <- gsub("\\]", "%5D", gsub("\\[", "%5B", p$ce))
      res <- curl_fetch_disk(sprintf("%s.dap.csv?dap4.ce=%s", base, ce), file, handle = handle)
      if (res$status_code != 200) { unlink(file); message("failed ", gran$title[j], " ", res$status_code); next }
    }
    d <- read_dap_csv(file)
    if (is.null(mask)) {
      # coordinates of the window, checked against what the server returned
      if (p$key == "smap_sm") {
        rr <- idx_ease[["r0"]]:idx_ease[["r1"]]; cc <- idx_ease[["c0"]]:idx_ease[["c1"]]
        xy <- expand.grid(c = cc, r = rr)
        pts <- st_transform(st_as_sf(data.frame(x = ease$x0 + (xy$c + 0.5) * ease$cell,
                                                y = ease$y0 - (xy$r + 0.5) * ease$cell),
                                     coords = c("x", "y"), crs = 6933), 4326)
        ll <- st_coordinates(pts)
        lat_server <- as.vector(t(d[["/Soil_Moisture_Retrieval_Data_AM/latitude"]]))
        ok <- lat_server > -9999
        stopifnot(!any(ok) || max(abs(lat_server[ok] - ll[ok, 2])) < 0.05)
        lon <- ll[, 1]; lat <- ll[, 2]
      } else {
        lat_name <- grep("lat", names(d), ignore.case = TRUE, value = TRUE)[1]
        lon_name <- grep("lon", names(d), ignore.case = TRUE, value = TRUE)[1]
        g <- expand.grid(lon = d[[lon_name]], lat = d[[lat_name]])
        lon <- ifelse(g$lon > 180, g$lon - 360, g$lon); lat <- g$lat
      }
      mask <- in_case(lon, lat, p$case)
    }
    fields <- list(c(key = p$key, variable = p$variable, field = p$field))
    if (!is.null(p$field2)) fields <- c(fields, list(p$field2))
    for (f in fields) {
      m <- d[[f[["field"]]]]
      # cells must follow expand.grid(lon, lat): lon varies fastest
      x <- if (isTRUE(p$lon_major)) as.vector(m) else as.vector(t(m))
      coverage[[length(coverage) + 1]] <- data.frame(
        key = f[["key"]], case = p$case, variable = f[["variable"]], archive = p$archive,
        month = gran$month[j], granule = gran$title[j], n_cells = sum(mask),
        coverage = mean(p$valid(x[mask])))
    }
    if (p$key == "mod13c2") {
      x <- as.vector(t(d[[p$field]]))[mask]
      reliability[[length(reliability) + 1]] <- data.frame(month = gran$month[j],
        rank = c(-1, 0, 1, 2, 3, 4), share = vapply(c(-1, 0, 1, 2, 3, 4), function(k) mean(x == k), 0))
    }
  }
  message(sprintf("%-14s %d granules", p$key, nrow(gran)))
}
unlink(cookies)
coverage <- do.call(rbind, coverage)

# daily product: a month counts a cell as covered if it has at least one valid day
cov_month <- aggregate(coverage ~ key + case + variable + archive + month + n_cells,
                       data = coverage[coverage$key != "smap_sm", ], FUN = mean)
if (any(coverage$key == "smap_sm")) {
  sm <- coverage[coverage$key == "smap_sm", ]
  cov_month <- rbind(cov_month, aggregate(coverage ~ key + case + variable + archive + month + n_cells,
                                          data = sm, FUN = mean))
}
cov_month$calendar_month <- as.integer(substr(cov_month$month, 6, 7))
summary_tab <- do.call(rbind, lapply(split(cov_month, cov_month$key), function(x) data.frame(
  key = x$key[1], case = x$case[1], variable = x$variable[1], archive = x$archive[1],
  months = paste(min(x$month), max(x$month), sep = " to "), n_months = nrow(x), cells = max(x$n_cells),
  coverage_median = median(x$coverage), coverage_min = min(x$coverage),
  wettest_quarter_median = median(x$coverage[x$calendar_month %in% 3:5]),
  driest_quarter_median = median(x$coverage[x$calendar_month %in% 9:11]))))
summary_tab$note <- ifelse(summary_tab$key == "smap_sm",
  "daily product: mean fraction of cells with a valid retrieval per day (sampled months only)", "")

write.csv(cov_month, file.path(out_dir, "nasa_coverage_monthly.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "nasa_coverage_summary.csv"), row.names = FALSE)
if (length(reliability)) {
  rel <- aggregate(share ~ rank, data = do.call(rbind, reliability), FUN = mean)
  rel$meaning <- c("fill", "good", "marginal", "snow/ice", "cloudy", "estimated from MODIS history")[match(rel$rank, -1:4)]
  write.csv(rel, file.path(out_dir, "mod13c2_reliability.csv"), row.names = FALSE)
  print(rel)
}

fig <- aggregate(coverage ~ case + variable + calendar_month, data = cov_month, FUN = median)
p <- ggplot(fig, aes(factor(calendar_month), variable, fill = coverage)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f", 100 * coverage)), size = 3) +
  facet_grid(case ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#f8ecd9", high = "#0d6b72", limits = c(0, 1), labels = scales::percent) +
  labs(x = "Calendar month", y = NULL, fill = "Valid cells",
       title = "NASA Earthdata L3 products: valid cells in the FZA case studies",
       subtitle = sprintf("Monthly products: %s to %s. SMAP soil moisture: daily mean, %s.",
                          min(params$sample_months), max(params$sample_months),
                          paste(params$smap_sm_months, collapse = " and "))) +
  theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "fza_nasa_coverage.png"), p, width = 8, height = 4.5, dpi = 150, bg = "white")

print(summary_tab, digits = 2)
