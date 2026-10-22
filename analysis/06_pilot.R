# Stage 6 — Pilot: does the correction change the results of PCA, MDS, k-means and SOM?
#
# Scenarios (analysis/d4_scenarios.R, pilot_scenarios): S n=400 high rough and smooth,
# T T=480 high and ST 50x60, each in three truths — null (no association, no classes:
# any structure found is spurious, H1), structured (R_true in blocks + 3 random classes:
# recovery, H3/H4) and contiguous (3 classes contiguous in space or time: Q3).
# 20 replicates per scenario. D14 adds two truths to each family: blocks that cross the
# regime groups ("crossed") and correlation 0.3 between groups (ids 13 on).
# Contiguous truths also get contiguity-constrained clustering (D14).
#
# Approaches: none; independent (one basis per variable); default (decision D7: global
# basis in space and time, whitening in space-time); groups (one basis per group);
# variables (D14: each variable selects its basis, the union is shared);
# whiten (space and time only; it is the default in space-time); broad (projection on
# the broad band of scales, contiguous truth only: what the scale mode keeps); sbss and
# maf (spatial only). Each approach gives a correlation matrix (PCA) and a data matrix
# (MDS with Euclidean and Gower distances, k-means, SOM), scored by run_downstream().
#
# Seeds: data 6600 + 1000 * scenario + replicate; the same number seeds the bases,
# permutations, k-means and SOM.
#
# Outputs (not tracked by git): results/pilot/scenario_NN.rds, results/pilot_replicates.rds,
#   results/pilot_summary.csv, results/pilot_contrasts.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 6600L
n_rep <- 20L
out_dir <- file.path("results", "pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

one_pilot <- function(sc, r) {
  seed <- base_seed + 1000L * sc$id + r
  t_rep <- Sys.time()
  s <- simulate(sc, seed, sc$truth)
  truth_tab <- truth_table(s$groups)
  n <- nrow(s$X)
  large_part <- sweep(s$components$large, 2, sqrt(truth_tab$share), "*")
  position <- switch(s$family, spatial = s$coords, temporal = matrix(s$time),
                     spatiotemporal = cbind(s$coords[s$index$site, ], s$index$time / max(s$index$time)))
  truth_all <- list(R_true = s$R_true, classes = if (sc$truth == "null") NULL else s$classes,
                    signal = s$X - large_part, position = position)
  rows_tail <- function(m) { keep <- utils::tail(seq_len(n), m)                 # whitening drops the first times
    list(R_true = truth_all$R_true, classes = truth_all$classes[keep], signal = truth_all$signal[keep, , drop = FALSE],
         position = truth_all$position[keep, , drop = FALSE]) }
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))

  sig <- diagnose_variables(s)
  type <- unname(classify_regime(sig))
  gr <- group_regimes(sig)
  approaches <- list()
  timed <- function(expr) { t0 <- Sys.time(); v <- expr; attr(v, "seconds") <- as.numeric(difftime(Sys.time(), t0, units = "secs")); v }
  approaches$none <- timed(list(X = s$X, R = stats::cor(s$X)))
  approaches$independent <- timed({
    E <- vapply(seq_len(ncol(s$X)), function(j) if (type[j] == "none") s$X[, j] - mean(s$X[, j]) else
      .resid_on(s$X[, j, drop = FALSE], do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j],
        period = if (is.finite(sig$period[j])) round(sig$period[j]) else 12, seed = seed + j), args))), numeric(n))
    list(X = E, R = stats::cor(E))
  })
  approaches$default <- timed(correct(s, gr, seed = seed)[c("X", "R")])
  approaches$groups <- timed(correct(s, gr, mode = "conditional", scope = "groups", seed = seed)[c("X", "R")])
  # D14: each variable selects its basis, the union is shared (Peres-Neto & Legendre 2010)
  approaches$variables <- timed(correct(s, gr, mode = "conditional", scope = "variables", seed = seed)[c("X", "R")])
  if (s$family != "spatiotemporal") approaches$whiten <- timed(correct(s, gr, mode = "whiten", seed = seed)[c("X", "R")])
  if (sc$truth == "contiguous" && s$family != "spatiotemporal") approaches$broad <- timed({
    cand <- if (s$family == "spatial") .spatial_candidates(s$coords) else .temporal_candidates(s$time, "ar")
    band <- cand[, seq_len(ceiling(ncol(cand) / 3)), drop = FALSE]
    P <- band %*% crossprod(band, scale(s$X, scale = FALSE))
    list(X = P, R = stats::cor(P))
  })
  if (s$family == "spatial") {
    approaches$sbss <- timed(decorrelate_latent(s$X, s$coords, "sbss", seed = seed)[c("X", "R")])
    approaches$maf <- timed(decorrelate_latent(s$X, s$coords, "maf", seed = seed)[c("X", "R")])
  }
  res <- do.call(rbind, lapply(names(approaches), function(a) {
    A <- approaches[[a]]
    tr <- if (nrow(A$X) < n) rows_tail(nrow(A$X)) else truth_all
    m <- rbind(run_downstream(A$X, "pca", truth = tr, R = A$R),
               run_downstream(A$X, "mds", truth = tr, distance = "euclidean", seed = seed),
               run_downstream(A$X, "mds", truth = tr, distance = "gower", seed = seed),
               run_downstream(A$X, "kmeans", truth = tr, seed = seed),
               run_downstream(A$X, "som", truth = tr, seed = seed),
               # D14: contiguity-constrained clustering, the reference for contiguous classes
               if (sc$truth == "contiguous") run_downstream(A$X, "constrained", truth = tr, seed = seed))
    cbind(approach = a, m, seconds_correction = attr(A, "seconds"))
  }))
  cbind(scenario = sc$id, replicate = r, seed = seed, res,
        seconds_replicate = as.numeric(difftime(Sys.time(), t_rep, units = "secs")))
}

# --- Run -------------------------------------------------------------------------------------
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
clusterExport(cl, c("base_seed", "one_pilot"))
cat(format(Sys.time()), "start\n", file = log_file, append = TRUE)
for (k in pilot_scenarios$id) {
  sc <- pilot_scenarios[pilot_scenarios$id == k, ]
  f <- file.path(out_dir, sprintf("scenario_%02d.rds", k))
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- parLapply(cl, seq_len(n_rep), function(r, sc) one_pilot(sc, r), sc = sc)
  saveRDS(do.call(rbind, res), f)
  cat(format(Sys.time()), sc$label, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n",
      file = log_file, append = TRUE)
}
stopCluster(cl)

# --- Summaries ---------------------------------------------------------------------------------
all <- do.call(rbind, lapply(pilot_scenarios$id, function(k) readRDS(file.path(out_dir, sprintf("scenario_%02d.rds", k)))))
all <- merge(all, pilot_scenarios[, c("id", "label", "family", "smoothness", "truth", "between")], by.x = "scenario", by.y = "id")
saveRDS(all, file.path("results", "pilot_replicates.rds"))
summ <- do.call(rbind, lapply(split(all, list(all$scenario, all$approach, all$method, all$metric), drop = TRUE), function(d) {
  v <- d$value[is.finite(d$value)]
  data.frame(d[1, c("scenario", "label", "family", "smoothness", "truth", "between", "approach", "method", "metric")],
             n = length(v), mean = if (length(v)) mean(v) else NA, median = if (length(v)) stats::median(v) else NA,
             q05 = if (length(v)) stats::quantile(v, 0.05) else NA, q95 = if (length(v)) stats::quantile(v, 0.95) else NA,
             sd = if (length(v) > 1) stats::sd(v) else NA, row.names = NULL)
}))
write.csv(summ, file.path("results", "pilot_summary.csv"), row.names = FALSE)

# paired contrasts against "none" (same data in every replicate): effect, its SD, and the
# replicates needed to detect it (two-sided 5%, power 80%) at stage 7
contr <- do.call(rbind, lapply(split(all, list(all$scenario, all$method, all$metric), drop = TRUE), function(d) {
  base <- d[d$approach == "none", c("replicate", "value")]
  do.call(rbind, lapply(setdiff(unique(d$approach), "none"), function(a) {
    x <- merge(d[d$approach == a, c("replicate", "value")], base, by = "replicate")
    diff <- x$value.x - x$value.y
    diff <- diff[is.finite(diff)]
    if (length(diff) < 3) return(NULL)
    m <- mean(diff); sdv <- stats::sd(diff)
    data.frame(d[1, c("scenario", "label", "family", "truth", "method", "metric")], approach = a, n = length(diff),
               effect = m, sd = sdv, replicates_needed = if (abs(m) > 1e-12 && sdv > 0) ceiling(((1.96 + 0.84) * sdv / abs(m))^2) else NA,
               row.names = NULL)
  }))
}))
write.csv(contr, file.path("results", "pilot_contrasts.csv"), row.names = FALSE)
times <- aggregate(seconds_replicate ~ scenario + label, data = unique(all[, c("scenario", "label", "replicate", "seconds_replicate")]), FUN = mean)
write.csv(times, file.path("results", "pilot_times.csv"), row.names = FALSE)
cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
print(times)
