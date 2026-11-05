# Manifest of the raw data: for every file under data/raw/, its size, its SHA-256, the
# date it was written (the files are kept as downloaded, so this is the access date), its
# source and the product with its version. The raw data are not tracked (some need an
# Earthdata login, and the L3 products are reprocessed from time to time), so the manifest
# is how a later run knows whether it has the same inputs.
#
# Usage, from the root of the repository:
#   Rscript analysis/00_raw_manifest.R          writes data/raw_manifest.csv (tracked)
#   Rscript analysis/00_raw_manifest.R check    compares data/raw/ with the manifest and
#                                               lists the files missing, new or changed
# No random numbers. data/raw/ is only read.

raw <- file.path("data", "raw")
manifest_file <- file.path("data", "raw_manifest.csv")

# Source and product of each folder (first match wins), from the scripts that wrote them.
sources <- data.frame(stringsAsFactors = FALSE, matrix(byrow = TRUE, ncol = 4, dimnames = list(NULL, c("pattern", "source", "product", "script")), c(
  "^fza/anp_bacias\\.json$", "ANP WFS (geomaps.anp.gov.br)", "sedimentary basins of Brazil", "01a",
  "^fza/etopo1_search_box\\.csv$", "NOAA CoastWatch ERDDAP", "ETOPO1 (etopo180)", "01a",
  "^fza/cube/(sst_night|chlorophyll|kd490|par|adg443)/", "NASA OB.DAAC OPeNDAP (Earthdata login)", "MODIS-Aqua L3 mapped monthly 4 km", "09a",
  "^fza/cube/salinity/", "NASA PO.DAAC OPeNDAP (Earthdata login)", "RSS SMAP L3 monthly sea surface salinity", "09a",
  "^fza/cube/precipitation/", "NASA GES DISC OPeNDAP (Earthdata login)", "GPM IMERG monthly (GPM_3IMERGM)", "09a",
  "^fza/cube/wind/", "NOAA CoastWatch ERDDAP", "Metop ASCAT monthly wind", "09a",
  "^fza/drivers/ana_", "ANA HidroWeb web service", "discharge at Obidos, station 17050001", "09a",
  "^fza/drivers/(oni|tna|tsa)\\.data$", "NOAA PSL", "climate index (ONI, TNA or TSA)", "09a",
  "^fza/erddap/", "NOAA CoastWatch ERDDAP", "ERDDAP dataset (id in the file name)", "01b",
  "^fza/nasa/imerg_", "NASA GES DISC OPeNDAP (Earthdata login)", "GPM IMERG monthly", "01c",
  "^fza/nasa/mod11c3", "NASA LP DAAC OPeNDAP (Earthdata login)", "MOD11C3 land surface temperature", "01c",
  "^fza/nasa/mod13c2", "NASA LP DAAC OPeNDAP (Earthdata login)", "MOD13C2 vegetation indices", "01c",
  "^fza/nasa/modis_adg443", "NASA OB.DAAC OPeNDAP (Earthdata login)", "MODIS-Aqua L3 mapped monthly IOP adg 443", "01c",
  "^fza/nasa/smap_sm", "NASA NSIDC OPeNDAP (Earthdata login)", "SMAP L3 soil moisture (SPL3SMP)", "01c",
  "^fza/nasa/smap_sss", "NASA PO.DAAC OPeNDAP (Earthdata login)", "RSS SMAP L3 monthly sea surface salinity", "01c")))

# The granule or dataset of a file: the "Dataset:" line that OPeNDAP writes at the top of a
# CSV subset, else the file name; the version is the version token of that name.
granule_of <- function(path) {
  first <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  if (length(first) && grepl("^Dataset: ", first)) sub("^Dataset: ", "", first) else basename(path)
}
version_of <- function(g) {
  v <- regmatches(g, regexpr("V0[0-9][A-Z]?|v[0-9]+\\.[0-9]+|R20[0-9]{2}[A-Z]*|\\.06[0-9]\\.|_R[0-9]{5}_", g))
  if (length(v)) gsub("^[._]|[._]$", "", v) else ""
}

files <- sort(list.files(raw, recursive = TRUE, all.files = FALSE))
files <- files[files != ".gitkeep"]
rel <- files
full <- file.path(raw, files)
src <- t(vapply(rel, function(f) {
  i <- which(vapply(sources$pattern, grepl, TRUE, x = f))[1]
  if (is.na(i)) c("unknown", "unknown", "") else unlist(sources[i, c("source", "product", "script")])
}, character(3)))
granule <- vapply(full, granule_of, "")
manifest <- data.frame(file = rel, bytes = file.size(full), sha256 = vapply(full, function(f) digest::digest(file = f, algo = "sha256"), ""),
                       written = format(file.mtime(full), "%Y-%m-%d"), source = src[, 1], product = src[, 2],
                       granule = granule, version = vapply(granule, version_of, ""), script = src[, 3],
                       row.names = NULL, stringsAsFactors = FALSE)

if (identical(commandArgs(trailingOnly = TRUE)[1], "check")) {
  old <- utils::read.csv(manifest_file, stringsAsFactors = FALSE, colClasses = "character")
  missing <- setdiff(old$file, manifest$file)
  new <- setdiff(manifest$file, old$file)
  both <- intersect(old$file, manifest$file)
  changed <- both[old$sha256[match(both, old$file)] != manifest$sha256[match(both, manifest$file)]]
  cat(sprintf("%d files in the manifest, %d in data/raw: %d missing, %d new, %d changed\n",
              nrow(old), nrow(manifest), length(missing), length(new), length(changed)))
  for (x in list(missing = missing, new = new, changed = changed)) if (length(x)) print(utils::head(x, 20))
  quit(status = as.integer(length(missing) + length(changed) > 0))
}
utils::write.csv(manifest, manifest_file, row.names = FALSE)
cat(sprintf("%d files, %.0f MB; unknown source: %d\n", nrow(manifest), sum(manifest$bytes) / 1e6, sum(manifest$source == "unknown")))
print(table(manifest$script))
