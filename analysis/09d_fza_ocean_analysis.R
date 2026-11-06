# Stage 9 — the ocean case of the Foz do Amazonas, as registered in protocol-v1.3
# (interno/docs/protocol.md, section 16; decisions D26 and D27).
#
# Data: the cube of analysis/09a (628 cells of 0.25 degrees x 93 months, April 2015 to
# December 2022, 9 variables). A cell-month of a variable is observed when at least
# `min_fraction` of the water pixels of the cell are valid (products on the 0.25-degree
# grid itself have fraction 1 when observed); otherwise it is missing (D26: the package
# takes missing values, no imputation). Chlorophyll, Kd490 and CDOM absorption enter as
# base-10 logarithms.
#
# Runs (section 16.6), each compared with the raw data at the same valid fraction:
#   main          fraction 0.25, no drivers, R by D7 (whitening with the group bases)
#   drivers       fraction 0.25, drivers of D27, R by D7, with the variation partition
#   fraction_0    fraction 0,    no drivers, R by D7
#   fraction_0.5  fraction 0.5,  no drivers, R by D7
#   conditional   fraction 0.25, no drivers, R in the conditional mode on the group bases
#                 of the main run (only its R is used; it is the Q3 fit of the main run)
# The data matrix X is always corrected with one basis per variable (D17).
#
# Output methods (sections 16.7 and 16.8):
# * PCA on the raw pairwise correlations and on R: axes above the broken-stick, matched
#   by the largest absolute Tucker congruence of the loadings (same >= 0.95, similar
#   0.85-0.95, changed < 0.85); reading 1 (salinity opposite to CDOM, Kd490 and
#   chlorophyll, each |loading| >= 0.3) for each axis;
# * R by pair: changes when |r_corrected - r_raw| >= 0.2, reverses when the signs differ
#   with both |r| >= 0.2; three predicted pairs keep their sign with |r| >= 0.2;
# * k-means (k by the gap statistic, k <= 6), the SOM (5 x 5 hexagonal; units cut into 3
#   classes, as in section 6 without a truth), HDBSCAN and DBSCAN (D18), and classical MDS, on the
#   cell-months where every variable has a value, fitted on a seeded subsample of 3,000
#   such rows (section 6), the same rows for raw and corrected data. The ARI between the
#   raw and corrected partitions is computed on those rows (changes below 0.65); every
#   complete row is then assigned (nearest centre, nearest map unit, or the dbscan
#   prediction) to give the modal class of each cell, whose contiguity (share of queen
#   neighbour pairs with the same modal class) is compared with 999 permutations of the
#   modal classes over the cells: regions when above the 95th percentile;
# * verify() (main and drivers runs): local Moran's I of the residuals by region (shelf:
#   depth above 200 m, as in analysis/09b), and the Q3 warning with the space-time
#   threshold of protocol-v1.2, computed as calibrated (analysis/07d): on a conditional
#   correction with the group bases of the run (with the drivers in the drivers run) and
#   the data matrix on their union, verify(nperm = 0);
# * the variation partition of each variable (drivers run).
#
# Seed: 17500 (diagnosis permutations, grouping, bases, verification, subsample and
# techniques through the package's streams; the contiguity permutations use
# set.seed(17500 + 1)). Recorded in results/fza/ocean_session.csv with the package commit.
#
# Timing: DEPPREP_OCEAN_SITES = n runs on the n cells nearest the centre of the basin, and
# writes to results/fza/ocean_subset_*.
#
# Resuming: each correction is saved, with its verification and the conditional fit of
# its Q3 statistic, in results/fza/ocean_correction_<run>.rds; a run that finds that file
# reads it instead of correcting again (same seed, same result). Delete the files to run
# again from the start.
#
# Requires results/fza/ocean_cube.csv, ocean_cells.csv, drivers.csv (analysis/09a) and the
# cached files of data/raw/fza/drivers and data/raw/fza/etopo1_search_box.csv.
# Outputs (not tracked), in results/fza/: ocean_session.csv, ocean_groups.csv,
# ocean_runs.csv, ocean_pca.csv, ocean_loadings.csv, ocean_pairs.csv, ocean_classes.csv,
# ocean_mds.csv, ocean_modal_classes.csv, ocean_partition.csv, ocean_lisa.csv,
# ocean_verify_variables.csv and the figures ocean_classes_<run>.png and
# ocean_pca_scores_<run>.png.

suppressMessages(devtools::load_all(quiet = TRUE))

params <- list(
  base_seed = 17500L,
  min_fraction = 0.25,
  log_variables = c("chlorophyll", "kd490", "adg443"),
  variables = c("sst_night", "salinity", "precipitation", "wind_u", "wind_v", "par", "kd490", "adg443", "chlorophyll"),
  predicted_pairs = list(c("salinity", "adg443", -1), c("chlorophyll", "kd490", 1), c("par", "precipitation", -1)),
  congruence = c(same = 0.95, similar = 0.85),
  delta_r = 0.2, ari_changes = 0.65, n_perm = 999, max_n = 3000, k_max = 6, B_gap = 20,
  cores = max(1L, min(16L, parallel::detectCores() - 4L)))
seed <- params$base_seed
out_dir <- file.path("results", "fza")
subset_sites <- as.integer(Sys.getenv("DEPPREP_OCEAN_SITES", "0"))
region <- Sys.getenv("DEPPREP_OCEAN_REGION", "")           # "shelf" or "slope and deep": a design of protocol-v1.3.1, section 17.4
prefix <- if (subset_sites > 0) "ocean_subset_" else if (nzchar(region)) paste0("ocean_", sub(" .*", "", region), "_") else "ocean_"
out_file <- function(name) file.path(out_dir, paste0(prefix, name))
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")      # socket workers of correct(cores =)

runs <- data.frame(run = c("main", "drivers", "fraction_0", "fraction_0.5", "conditional"),
                   min_fraction = c(0.25, 0.25, 0, 0.5, 0.25), drivers = c(FALSE, TRUE, FALSE, FALSE, FALSE),
                   r_mode = c("whiten", "whiten", "whiten", "whiten", "conditional"), stringsAsFactors = FALSE)

# --- Data ------------------------------------------------------------------------------
cube <- read.csv(file.path(out_dir, "ocean_cube.csv"), stringsAsFactors = FALSE)
cells <- read.csv(file.path(out_dir, "ocean_cells.csv"))
topo <- read.csv(file.path("data", "raw", "fza", "etopo1_search_box.csv"), skip = 2, header = FALSE,
                 col.names = c("lat", "lon", "altitude"))
depth <- tapply(topo$altitude, paste(floor(topo$lon / 0.25), floor(topo$lat / 0.25)), stats::median)
cells$region <- ifelse(as.vector(depth[paste(floor(cells$lon / 0.25), floor(cells$lat / 0.25))]) > -200, "shelf", "slope and deep")
if (nzchar(region)) {
  cells <- cells[cells$region == region, ]
  runs <- runs[runs$run %in% c("main", "conditional"), ]
}
if (subset_sites > 0) {
  d <- sqrt((cells$lon - mean(cells$lon))^2 + (cells$lat - mean(cells$lat))^2)
  cells <- cells[sort(order(d)[seq_len(subset_sites)]), ]
}
months <- sort(unique(cube$month))
cube <- cube[cube$cell %in% cells$cell, ]

# Rows: every cell at every month, month by month; NA where a variable is not observed.
ocean_data <- function(f) {
  ok <- !is.na(cube$value) & cube$valid_pixels > 0 & cube$valid_fraction >= f
  v <- ifelse(ok, cube$value, NA_real_)
  lg <- cube$variable %in% params$log_variables
  v[lg] <- log10(v[lg])
  site <- match(cube$cell, cells$cell); time <- match(cube$month, months)
  X <- matrix(NA_real_, nrow(cells) * length(months), length(params$variables),
              dimnames = list(NULL, params$variables))
  X[cbind((time - 1) * nrow(cells) + site, match(cube$variable, params$variables))] <- v
  X
}
index <- data.frame(site = rep(seq_len(nrow(cells)), length(months)), time = rep(seq_along(months), each = nrow(cells)))
coords <- as.matrix(cells[, c("lon", "lat")])
half <- ifelse(as.integer(substr(months, 6, 7)) <= 6, "January-June", "July-December")[index$time]

# Drivers (D27): discharge at lags 0-2 (the months before April 2015 from the cached ANA
# file of 2015, consisted data preferred, as in analysis/09a) and ONI, TNA, TSA at lag 0.
drv <- read.csv(file.path(out_dir, "drivers.csv"), stringsAsFactors = FALSE)
s <- paste(readLines(file.path("data", "raw", "fza", "drivers", "ana_17050001_2015.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
recs <- regmatches(s, gregexpr("<SerieHistorica[^>]*>.*?</SerieHistorica>", s, perl = TRUE))[[1]]
tag <- function(r, t) { m <- regmatches(r, regexpr(sprintf("<%s>[^<]*</%s>", t, t), r)); if (length(m)) sub("<[^>]*>([^<]*)<.*", "\\1", m) else NA }
ana <- data.frame(month = substr(vapply(recs, tag, "", "DataHora"), 1, 7), level = as.integer(vapply(recs, tag, "", "NivelConsistencia")),
                  discharge = as.numeric(vapply(recs, tag, "", "Media")), row.names = NULL)
ana <- ana[order(ana$month, -ana$level), ]; ana <- ana[!duplicated(ana$month) & !is.na(ana$discharge), ]
before <- c("2015-02", "2015-03")
q <- c(ana$discharge[match(before, ana$month)], drv$discharge_obidos[match(months, drv$month)])
stopifnot(!anyNA(q), !anyNA(drv[match(months, drv$month), c("oni", "tna", "tsa")]))
n_t <- length(months)
covariates <- data.frame(time = seq_len(n_t), discharge = q[3:(n_t + 2)], discharge_lag1 = q[2:(n_t + 1)],
                         discharge_lag2 = q[1:n_t], drv[match(months, drv$month), c("oni", "tna", "tsa")], row.names = NULL)

# --- Helpers ---------------------------------------------------------------------------
broken_stick <- function(values) {
  p <- length(values); bs <- rev(cumsum(1 / rev(seq_len(p))))
  above <- values > bs
  if (above[1]) which.min(c(above, FALSE)) - 1 else 0
}
pca_of <- function(R) {
  e <- eigen(R, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
  dimnames(L) <- list(colnames(R), paste0("PC", seq_len(ncol(R))))
  list(values = e$values, loadings = L, n_bs = broken_stick(e$values))
}
tucker <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
reading_1 <- function(l) {
  o <- l[c("adg443", "kd490", "chlorophyll")]
  abs(l[["salinity"]]) >= 0.3 && all(abs(o) >= 0.3) && all(sign(o) == -sign(l[["salinity"]]))
}
mode_of <- function(x) { t <- table(x); as.integer(names(t)[which.max(t)]) }
# queen neighbours on the 0.25-degree grid
nb_pairs <- local({
  D <- pmax(abs(outer(coords[, 1], coords[, 1], "-")), abs(outer(coords[, 2], coords[, 2], "-")))
  w <- which(D > 0 & D < 0.25 * 1.01 & upper.tri(D), arr.ind = TRUE)
  unname(w)
})
contiguity <- function(cl, perm_seed) {          # cl: modal class of each cell (NA: none)
  ok <- !is.na(cl[nb_pairs[, 1]]) & !is.na(cl[nb_pairs[, 2]])
  pr <- nb_pairs[ok, , drop = FALSE]
  share <- function(z) mean(z[pr[, 1]] == z[pr[, 2]])
  obs <- share(cl)
  set.seed(perm_seed)
  have <- which(!is.na(cl))
  null <- replicate(params$n_perm, { z <- cl; z[have] <- z[have][sample.int(length(have))]; share(z) })
  c(same_share = obs, p95 = unname(stats::quantile(null, 0.95)), regions = as.numeric(obs > stats::quantile(null, 0.95)))
}
nearest <- function(Z, centres) {
  d <- as.matrix(stats::dist(rbind(centres, Z)))[-seq_len(nrow(centres)), seq_len(nrow(centres)), drop = FALSE]
  max.col(-d, ties.method = "first")
}
nearest_chunked <- function(Z, centres, chunk = 5000) {
  unlist(lapply(split(seq_len(nrow(Z)), ceiling(seq_len(nrow(Z)) / chunk)), function(i) nearest(Z[i, , drop = FALSE], centres)))
}

# Partitions of one data matrix: fitted on the subsample `sub` of the complete rows Z,
# then every complete row assigned. Settings as in run_downstream().
partitions <- function(Z, sub) {
  Zs <- Z[sub, , drop = FALSE]
  .with_stream(seed, .streams[["downstream"]], {
    km <- function(x, k) stats::kmeans(x, k, nstart = 10, iter.max = 50)
    gap <- cluster::clusGap(Zs, km, K.max = params$k_max, B = params$B_gap, verbose = FALSE)
    k <- cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
    fit <- if (k > 1) km(Zs, k) else NULL
    kmeans_all <- if (k > 1) nearest_chunked(Z, fit$centers) else rep(1L, nrow(Z))
    grid <- kohonen::somgrid(5, 5, "hexagonal")
    sm <- kohonen::som(Zs, grid = grid, rlen = 100)
    codes <- sm$codes[[1]]
    unit_class <- stats::cutree(stats::hclust(stats::dist(codes), method = "ward.D2"), 3)
    som_all <- unit_class[nearest_chunked(Z, codes)]
    min_pts <- max(10L, as.integer(ceiling(0.02 * nrow(Zs))))
    hd <- dbscan::hdbscan(Zs, minPts = min_pts)
    hd_all <- if (any(hd$cluster > 0)) stats::predict(hd, newdata = Z, data = Zs) else rep(0L, nrow(Z))   # all noise
    kd <- dbscan::kNNdist(Zs, k = min_pts - 1)
    eps_grid <- unique(stats::quantile(kd, seq(0.1, 0.9, by = 0.1), names = FALSE))
    Dm <- as.matrix(stats::dist(Zs))
    sil <- vapply(eps_grid, function(e) {
      cl <- dbscan::dbscan(Zs, eps = e, minPts = min_pts)$cluster; ok <- cl > 0
      if (length(unique(cl[ok])) < 2) return(NA_real_)
      mean(cluster::silhouette(cl[ok], stats::as.dist(Dm[ok, ok, drop = FALSE]))[, "sil_width"])
    }, 0)
    eps <- if (any(is.finite(sil))) eps_grid[which.max(sil)] else stats::median(kd)
    db <- dbscan::dbscan(Zs, eps = eps, minPts = min_pts)
    db_all <- if (any(db$cluster > 0)) stats::predict(db, newdata = Z, data = Zs) else rep(0L, nrow(Z))
    conf <- stats::cmdscale(stats::dist(Zs), k = 2)
    list(k = k, labels = list(kmeans = kmeans_all, som = som_all, hdbscan = hd_all, dbscan = db_all),
         eps = eps, silhouette = if (any(is.finite(sil))) max(sil, na.rm = TRUE) else NA_real_,
         mds = conf, mds_stress = sqrt(sum((stats::dist(Zs) - stats::dist(conf))^2) / sum(stats::dist(Zs)^2)))
  })
}

# Strip the bases (large in space-time) before saving a correction.
slim <- function(cr) {
  cr$kn_variables <- if (is.null(cr$variable_bases)) NA_real_ else
    vapply(cr$variable_bases, function(b) if (is.null(b)) 0 else ncol(b) / nrow(cr$X), 0)
  cr[c("bases", "variable_bases", "covariates")] <- NULL
  cr
}

# --- Diagnosis and grouping per valid fraction ----------------------------------------
t_start <- Sys.time()
package_commit <- system("git rev-parse HEAD", intern = TRUE)   # at the start: commits made during the run do not count
fractions <- unique(runs$min_fraction)
data_of <- setNames(lapply(fractions, ocean_data), fractions)
groupings <- list(); signatures <- list()
for (f in fractions) {
  key <- as.character(f)
  sig <- diagnose_variables(data_of[[key]], coords = coords, index = index, nsim = 499, perm_seed = seed)
  signatures[[key]] <- sig
  groupings[[key]] <- group_regimes(sig, seed = seed)
}
groups_table <- do.call(rbind, lapply(fractions, function(f) {
  key <- as.character(f); gr <- groupings[[key]]
  data.frame(min_fraction = f, variable = params$variables, missing = colMeans(is.na(data_of[[key]])),
             regime = gr$type[match(params$variables, gr$variable)], group = gr$group[match(params$variables, gr$variable)])
}))
cat("diagnosis and grouping:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")
print(groups_table, digits = 3, row.names = FALSE)

# --- Corrections and verification -----------------------------------------------------
results <- list(); q3_fits <- list()
for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]; key <- as.character(r$min_fraction)
  X <- data_of[[key]]; gr <- groupings[[key]]
  saved_file <- out_file(paste0("correction_", r$run, ".rds"))
  if (r$r_mode != "conditional" && file.exists(saved_file)) {     # resume an interrupted run
    sv <- readRDS(saved_file)
    if (!is.null(sv$verification) && is.null(sv$q3_fit)) {        # saved before the Q3 fit was kept:
      fit <- correct(X, gr, coords = coords, index = index, mode = "conditional", scope = "groups",   # rebuilt
                     data_basis = "union", seed = seed, covariates = if (r$drivers) covariates else NULL)
      stopifnot(isTRUE(all.equal(fit$kn, sv$correction$kn)))     # with the same group bases
      sv$q3_fit <- slim(fit); rm(fit)
      sv$q3_fit$q3 <- sv$verification$q3
      saveRDS(sv, saved_file)
    }
    if (!is.null(sv$q3_fit)) q3_fits[[r$run]] <- sv$q3_fit
    results[[r$run]] <- structure(sv$correction, verification = sv$verification)
    cat(sprintf("%s: resumed from %s\n", r$run, saved_file))
    next
  }
  t0 <- Sys.time()
  cr <- if (r$r_mode == "conditional") q3_fits$main else
    correct(X, gr, coords = coords, index = index, seed = seed, cores = params$cores,
            covariates = if (r$drivers) covariates else NULL)
  time_correct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  vf <- NULL
  if (r$run %in% c("main", "drivers")) {
    t1 <- Sys.time()
    vf <- verify(cr, X, coords = coords, index = index, seed = seed)
    # the Q3 statistic as calibrated: conditional mode on the group bases (whitening has no
    # projected part); the bases of the drivers run already hold the drivers
    fit <- correct(X, gr, coords = coords, index = index, bases = cr$bases, mode = "conditional",
                   scope = "groups", data_basis = "union", seed = seed)
    q3 <- verify(fit, X, coords = coords, index = index, nperm = 0, seed = seed)$q3
    q3_fits[[r$run]] <- slim(fit); rm(fit)
    q3_fits[[r$run]]$q3 <- q3
    vf$q3 <- q3
    vf$warnings <- c(vf$warnings, if (q3$warning) sprintf("Q3: %.3f > %.2f (conditional fit on the group bases)", q3$statistic, q3$threshold))
    vf$time <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  }
  if (r$r_mode != "conditional") cr <- slim(cr)
  attr(cr, "time") <- time_correct
  results[[r$run]] <- structure(cr, verification = vf)
  saveRDS(list(correction = cr, verification = if (!is.null(vf)) vf[setdiff(names(vf), c("R_large", "cov_large", "cov_residual"))],
               q3_fit = q3_fits[[r$run]]), saved_file)
  cat(sprintf("%s: correct() %.0f s%s\n", r$run, time_correct, if (is.null(vf)) "" else sprintf(", verify() %.0f s", vf$time)))
}

# --- Output methods --------------------------------------------------------------------
pairs <- t(utils::combn(params$variables, 2))
pca_rows <- list(); loading_rows <- list(); pair_rows <- list(); class_rows <- list(); mds_rows <- list(); modal_rows <- list()
for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]; key <- as.character(r$min_fraction)
  X <- data_of[[key]]; cr <- results[[r$run]]
  R_raw <- stats::cor(X, use = "pairwise.complete.obs")
  R_cor <- cr$R[params$variables, params$variables]
  # PCA
  pr <- pca_of(R_raw); pc <- pca_of(R_cor)
  q <- max(1, pr$n_bs, pc$n_bs)
  C <- abs(outer(seq_len(q), seq_len(q), Vectorize(function(a, b) tucker(pr$loadings[, a], pc$loadings[, b]))))
  match_of <- integer(q); free <- seq_len(q)
  for (a in order(-apply(C, 1, max))) { b <- free[which.max(C[a, free])]; match_of[a] <- b; free <- setdiff(free, b) }
  pca_rows[[r$run]] <- data.frame(run = r$run, axis = seq_len(q), raw_share = pr$values[seq_len(q)] / length(pr$values),
    raw_above_broken_stick = seq_len(q) <= pr$n_bs, corrected_axis = match_of,
    corrected_share = pc$values[match_of] / length(pc$values), corrected_above_broken_stick = match_of <= pc$n_bs,
    congruence = C[cbind(seq_len(q), match_of)],
    verdict = cut(C[cbind(seq_len(q), match_of)], c(-Inf, params$congruence[["similar"]], params$congruence[["same"]], Inf),
                  c("changed", "similar", "same"), right = FALSE),
    reading_1_raw = vapply(seq_len(q), function(a) reading_1(pr$loadings[, a]), TRUE),
    reading_1_corrected = vapply(match_of, function(b) reading_1(pc$loadings[, b]), TRUE))
  loading_rows[[r$run]] <- rbind(
    data.frame(run = r$run, data = "raw", variable = params$variables, pr$loadings[, seq_len(q), drop = FALSE], row.names = NULL),
    data.frame(run = r$run, data = "corrected", variable = params$variables,
               stats::setNames(as.data.frame(pc$loadings[, match_of, drop = FALSE]), paste0("PC", seq_len(q))), row.names = NULL))
  # R by pair
  rr <- R_raw[pairs]; rc <- R_cor[pairs]
  pair_rows[[r$run]] <- data.frame(run = r$run, a = pairs[, 1], b = pairs[, 2], r_raw = rr, r_corrected = rc, delta = rc - rr,
    changes = abs(rc - rr) >= params$delta_r, reverses = sign(rc) != sign(rr) & abs(rc) >= params$delta_r & abs(rr) >= params$delta_r,
    predicted_sign = NA_real_, prediction_holds = NA)
  for (pp in params$predicted_pairs) {
    j <- which((pairs[, 1] == pp[1] & pairs[, 2] == pp[2]) | (pairs[, 1] == pp[2] & pairs[, 2] == pp[1]))
    s0 <- as.numeric(pp[3])
    pair_rows[[r$run]]$predicted_sign[j] <- s0
    pair_rows[[r$run]]$prediction_holds[j] <- sign(rc[j]) == s0 && abs(rc[j]) >= params$delta_r
  }
  if (r$run == "conditional") next                         # only its R is used
  # techniques on observations
  Xc <- cr$X[, params$variables]
  whole <- which(stats::complete.cases(X) & stats::complete.cases(Xc))
  sub <- sort(.with_stream(seed, .streams[["subsample"]], sample(length(whole), min(params$max_n, length(whole)))))
  zs <- function(M) { Z <- scale(M[whole, , drop = FALSE]); Z[!is.finite(Z)] <- 0; Z }
  p_raw <- partitions(zs(X), sub); p_cor <- partitions(zs(Xc), sub)
  mds_rows[[r$run]] <- data.frame(run = r$run, rows = length(sub), stress_raw = p_raw$mds_stress, stress_corrected = p_cor$mds_stress,
                                  procrustes_raw_corrected = .procrustes_r(p_raw$mds, p_cor$mds))
  for (m in names(p_raw$labels)) {
    la <- p_raw$labels[[m]]; lb <- p_cor$labels[[m]]
    modal <- function(l, rows = seq_along(whole)) {
      cl <- rep(NA_integer_, nrow(cells))
      md <- tapply(l[rows], index$site[whole][rows], mode_of)
      cl[as.integer(names(md))] <- md
      cl
    }
    cg_raw <- contiguity(modal(la), seed + 1L); cg_cor <- contiguity(modal(lb), seed + 1L)
    class_rows[[paste(r$run, m)]] <- data.frame(run = r$run, method = m, complete_rows = length(whole), subsample = length(sub),
      k_raw = length(setdiff(unique(la[sub]), 0)), k_corrected = length(setdiff(unique(lb[sub]), 0)),
      noise_raw = mean(la[sub] == 0), noise_corrected = mean(lb[sub] == 0),
      ari_raw_corrected = mclust::adjustedRandIndex(la[sub], lb[sub]),
      changes = mclust::adjustedRandIndex(la[sub], lb[sub]) < params$ari_changes,
      same_share_raw = cg_raw[["same_share"]], p95_raw = cg_raw[["p95"]], regions_raw = cg_raw[["regions"]] == 1,
      same_share_corrected = cg_cor[["same_share"]], p95_corrected = cg_cor[["p95"]], regions_corrected = cg_cor[["regions"]] == 1,
      shelf_share_of_main_raw_class = {
        cl <- modal(la); tab <- table(cl, cells$region); top <- which.max(tab[, "shelf"]); tab[top, "shelf"] / sum(tab[top, ])
      })
    for (h in c("all", "January-June", "July-December")) {
      rows <- if (h == "all") seq_along(whole) else which(half[whole] == h)
      modal_rows[[paste(r$run, m, h)]] <- data.frame(run = r$run, method = m, season = h, cell = cells$cell,
                                                    raw = modal(la, rows), corrected = modal(lb, rows))
    }
  }
}

# --- Partition, verification and LISA ------------------------------------------------
partition <- if (!is.null(results$drivers$partition)) {
  p <- results$drivers$partition
  p$rank_drivers_only <- rank(-p$drivers_only, ties.method = "min")
  p
} else NULL
lisa_rows <- list(); verify_rows <- list()
for (run in intersect(c("main", "drivers"), names(results))) {
  vf <- attr(results[[run]], "verification")
  verify_rows[[run]] <- data.frame(run = run, vf$variables, row.names = NULL)
  for (when in c("before", "after")) {
    L <- if (when == "after") vf$lisa else vf$lisa_before
    if (is.null(L)) next
    L$region <- cells$region[L$site]
    a <- stats::aggregate(list(share_structure = L$cluster != "none", share_p05 = L$p < 0.05), L[c("variable", "region")], mean)
    lisa_rows[[paste(run, when)]] <- data.frame(run = run, data = if (when == "after") "residuals" else "raw site means", a)
  }
}
run_table <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]; cr <- results[[r$run]]; vf <- attr(cr, "verification")
  data.frame(r, rows = nrow(cr$X), missing = mean(is.na(data_of[[as.character(r$min_fraction)]])),
             complete_rows = sum(stats::complete.cases(data_of[[as.character(r$min_fraction)]])),
             time_correct = attr(cr, "time"), adjustment = cr$adjustment, kn_groups = paste(sprintf("%s %.3f", names(cr$kn), cr$kn), collapse = "; "),
             kn_variables_mean = mean(cr$kn_variables), kn_variables_max = max(cr$kn_variables),
             q3_statistic = if (is.null(vf)) NA else vf$q3$statistic, q3_threshold = if (is.null(vf)) NA else vf$q3$threshold,
             q3_warning = if (is.null(vf)) NA else vf$q3$warning,
             warnings = paste(unique(c(cr$warnings, if (!is.null(vf)) vf$warnings)), collapse = " | "))
}))

w <- function(d, name) if (!is.null(d) && NROW(d)) utils::write.csv(d, out_file(name), row.names = FALSE)
w(data.frame(seed = seed, perm_seed_contiguity = seed + 1L, package_commit = package_commit,
             date = format(Sys.time(), "%Y-%m-%d %H:%M"), cores = params$cores, sites = nrow(cells), months = length(months),
             minutes = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)), "session.csv")
w(groups_table, "groups.csv"); w(run_table, "runs.csv")
w(do.call(rbind, pca_rows), "pca.csv"); w(do.call(rbind, loading_rows), "loadings.csv")
w(do.call(rbind, pair_rows), "pairs.csv"); w(do.call(rbind, class_rows), "classes.csv")
w(do.call(rbind, mds_rows), "mds.csv"); w(do.call(rbind, modal_rows), "modal_classes.csv")
w(partition, "partition.csv"); w(do.call(rbind, lisa_rows), "lisa.csv"); w(do.call(rbind, verify_rows), "verify_variables.csv")

# --- Figures ---------------------------------------------------------------------------
modal_all <- do.call(rbind, modal_rows)
for (run in setdiff(runs$run, "conditional")) {
  grDevices::png(out_file(paste0("classes_", run, ".png")), width = 2400, height = 1800, res = 200)
  graphics::par(mfrow = c(2, 4), mar = c(2, 2, 2.5, 0.5))
  for (m in c("kmeans", "som")) for (h in c("January-June", "July-December")) for (d in c("raw", "corrected")) {
    x <- modal_all[modal_all$run == run & modal_all$method == m & modal_all$season == h, ]
    cl <- x[[d]][match(cells$cell, x$cell)]
    graphics::plot(coords, col = .class_colours(cl), pch = 15, cex = 1.1, asp = 1, xlab = "", ylab = "",
                   main = sprintf("%s, %s, %s", m, d, h), cex.main = 0.9)
  }
  grDevices::dev.off()
  # scores of the first raw axis and of its matched corrected axis, cell means by season
  key <- as.character(runs$min_fraction[runs$run == run])
  pr <- pca_of(stats::cor(data_of[[key]], use = "pairwise.complete.obs")); pc <- pca_of(results[[run]]$R[params$variables, params$variables])
  b <- do.call(rbind, pca_rows)[do.call(rbind, pca_rows)$run == run, "corrected_axis"][1]
  score <- function(M, v) { Z <- scale(M); Z[is.na(Z)] <- 0; as.vector(Z %*% (v / sqrt(sum(v^2)))) }
  sc <- list(raw = score(data_of[[key]], pr$loadings[, 1]), corrected = score(results[[run]]$X[, params$variables], pc$loadings[, b]))
  grDevices::png(out_file(paste0("pca_scores_", run, ".png")), width = 2000, height = 1800, res = 200)
  graphics::par(mfrow = c(2, 2), mar = c(2, 2, 2.5, 0.5))
  for (d in names(sc)) for (h in c("January-June", "July-December")) {
    m <- tapply(sc[[d]][half == h], index$site[half == h], mean)
    z <- m[as.character(seq_len(nrow(cells)))]
    col <- grDevices::hcl.colors(11, "Blue-Red 3")[cut(z, stats::quantile(z, seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE)]
    graphics::plot(coords, col = col, pch = 15, cex = 1.1, asp = 1, xlab = "", ylab = "",
                   main = sprintf("PC%d %s, %s", if (d == "raw") 1 else b, d, h), cex.main = 0.9)
  }
  grDevices::dev.off()
}

# --- Summary ---------------------------------------------------------------------------
print(run_table[, c("run", "time_correct", "kn_variables_mean", "q3_statistic", "q3_warning")], digits = 3, row.names = FALSE)
print(do.call(rbind, pca_rows), digits = 3, row.names = FALSE)
print(do.call(rbind, class_rows)[, c("run", "method", "k_raw", "k_corrected", "ari_raw_corrected", "regions_raw", "regions_corrected")],
      digits = 3, row.names = FALSE)
pt <- do.call(rbind, pair_rows); print(pt[!is.na(pt$predicted_sign) | pt$changes, ], digits = 2, row.names = FALSE)
if (!is.null(partition)) print(partition, digits = 3, row.names = FALSE)
cat("total:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")
