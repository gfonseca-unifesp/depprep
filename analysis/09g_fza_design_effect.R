# Stage 9 — with the same process and the same data, how much do the conclusions change with
# the extent of the grid? (protocol-v1.3.1, section 17.4; decision D31; descriptive.)
#
# Designs, with cells of 0.25 degrees and valid fraction 0.25: the whole basin (the main and
# conditional runs of analysis/09d, reused), the shelf (depth above 200 m) and the slope and
# deep water. Each part goes through the flow of 09d again, with a new diagnosis and grouping,
# for R with whitening, with the conditional mode and with conditional + AR, the default
# chosen in D29 (run analysis/09d first). A pair depends on the design when
# |r(design) - r(whole basin)| >= 0.2 under the same mode.
#
# Seed 17500, as in 09d, which records it in its session files.
# Outputs (not tracked): results/fza/ocean_shelf_* and ocean_slope_* (from 09d),
# results/fza/ocean_design_pairs.csv, ocean_design_summary.csv

base_seed <- 17500L
set.seed(base_seed)
out_dir <- file.path("results", "fza")
regions <- c(shelf = "shelf", slope = "slope and deep")
rscript <- file.path(R.home("bin"), "Rscript")
for (r in regions) {
  Sys.setenv(DEPPREP_OCEAN_REGION = r)
  status <- system2(rscript, "analysis/09d_fza_ocean_analysis.R")
  if (status != 0) stop("09d failed for the region ", r, call. = FALSE)
}
Sys.unsetenv("DEPPREP_OCEAN_REGION")

designs <- c(whole = "ocean_", shelf = "ocean_shelf_", slope = "ocean_slope_")
read_all <- function(name) do.call(rbind, lapply(names(designs), function(d) {
  x <- utils::read.csv(file.path(out_dir, paste0(designs[[d]], name, ".csv")), stringsAsFactors = FALSE)
  x <- x[x$run %in% c("main", "conditional", "conditional_ar"), ]
  cbind(design = d, x)
}))
pairs <- read_all("pairs")
whole <- pairs[pairs$design == "whole", ]
key <- function(x) paste(x$run, x$a, x$b)
m <- match(key(pairs), key(whole))
pairs$r_whole <- whole$r_corrected[m]
pairs$r_raw_whole <- whole$r_raw[m]
pairs$depends_on_design <- abs(pairs$r_corrected - pairs$r_whole) >= 0.2
pairs$raw_depends_on_design <- abs(pairs$r_raw - pairs$r_raw_whole) >= 0.2
utils::write.csv(pairs, file.path(out_dir, "ocean_design_pairs.csv"), row.names = FALSE)

run_rows <- read_all("runs")
pca <- read_all("pca")
first_axis <- pca[pca$axis == 1, ]
classes <- read_all("classes")
kmeans_rows <- classes[classes$method == "kmeans", ]
at <- function(x, d, r, col) { v <- x[[col]][x$design == d & x$run == r]; if (length(v)) v[1] else NA }
summary_rows <- do.call(rbind, lapply(names(designs), function(d) do.call(rbind, lapply(c("main", "conditional", "conditional_ar"), function(r) {
  p <- pairs[pairs$design == d & pairs$run == r, ]
  data.frame(design = d, run = r, cells = at(run_rows, d, r, "rows") / 93, q3 = at(run_rows, d, "main", "q3_statistic"),
             pairs_change_raw_corrected = sum(p$changes), pairs_depend_on_design = sum(p$depends_on_design),
             raw_pairs_depend_on_design = sum(p$raw_depends_on_design), predicted_hold = sum(p$prediction_holds, na.rm = TRUE),
             pc1_share_raw = at(first_axis, d, r, "raw_share"), pc1_congruence = at(first_axis, d, r, "congruence"),
             pc1_verdict = at(first_axis, d, r, "verdict"), reading_1_corrected = at(first_axis, d, r, "reading_1_corrected"),
             kmeans_k_raw = at(kmeans_rows, d, r, "k_raw"), kmeans_k_corrected = at(kmeans_rows, d, r, "k_corrected"),
             kmeans_ari = at(kmeans_rows, d, r, "ari_raw_corrected"), base_seed = base_seed)
}))))
utils::write.csv(summary_rows, file.path(out_dir, "ocean_design_summary.csv"), row.names = FALSE)
print(summary_rows, digits = 3, row.names = FALSE)
print(pairs[!is.na(pairs$predicted_sign), c("design", "run", "a", "b", "r_raw", "r_corrected", "prediction_holds")], digits = 2, row.names = FALSE)
