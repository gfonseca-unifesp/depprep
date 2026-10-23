# Stage 6, before the protocol is frozen — two quick checks on the pilot data:
# 1. HDBSCAN and DBSCAN (decision D18): does the density-based clustering change with
#    the correction, and with how much spread? HDBSCAN's contrasts enter the replicate
#    count of D8; DBSCAN (radius by the silhouette) is compared with it.
# 2. The PCA metric against the true signal (freeze plan, item 4): does the corrected PCA
#    stop moving away from the true axes once the axes include the class effects?
#
# The 12 scenarios of the pilot (analysis/d4_scenarios.R, pilot_scenarios ids 1-12), the
# same seeds (6600 + 1000 * scenario + replicate, 20 replicates) and the same corrections
# as analysis/06_pilot.R, with the data matrix of the shared bases on the union, as the
# pilot ran (before D17). Only HDBSCAN and PCA are scored.
#
# Outputs (not tracked): results/pilot_hdbscan/scenario_NN.rds,
#   results/pilot_hdbscan_summary.csv, results/pilot_hdbscan_contrasts.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 6600L
n_rep <- 20L
out_dir <- file.path("results", "pilot_hdbscan")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

one_check <- function(sc, r) {
  seed <- base_seed + 1000L * sc$id + r
  s <- simulate(sc, seed, sc$truth)
  truth_tab <- truth_table(s$groups)
  n <- nrow(s$X)
  large_part <- sweep(s$components$large, 2, sqrt(truth_tab$share), "*")
  tr <- list(R_true = s$R_true, classes = if (sc$truth == "null") NULL else s$classes, signal = s$X - large_part)
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))
  sig <- diagnose_variables(s)
  type <- unname(classify_regime(sig))
  gr <- group_regimes(sig)
  ap <- list(none = list(X = s$X, R = stats::cor(s$X)))
  E <- vapply(seq_len(ncol(s$X)), function(j) if (type[j] == "none") s$X[, j] - mean(s$X[, j]) else
    .resid_on(s$X[, j, drop = FALSE], do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j],
      period = if (is.finite(sig$period[j])) round(sig$period[j]) else 12, seed = seed + j), args))), numeric(n))
  ap$independent <- list(X = E, R = stats::cor(E))
  ap$default <- correct(s, gr, seed = seed, data_basis = "union")[c("X", "R")]
  ap$groups <- correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union")[c("X", "R")]
  ap$variables <- correct(s, gr, mode = "conditional", scope = "variables", seed = seed, data_basis = "union")[c("X", "R")]
  if (s$family != "spatiotemporal") ap$whiten <- correct(s, gr, mode = "whiten", seed = seed, data_basis = "union")[c("X", "R")]
  res <- do.call(rbind, lapply(names(ap), function(a) {
    A <- ap[[a]]
    keep <- utils::tail(seq_len(n), nrow(A$X))             # whitening drops the first times
    t_a <- list(R_true = tr$R_true, classes = tr$classes[keep], signal = tr$signal[keep, , drop = FALSE])
    rbind(cbind(approach = a, run_downstream(A$X, "pca", truth = t_a, R = A$R)),
          cbind(approach = a, run_downstream(A$X, "hdbscan", truth = t_a, seed = seed)),
          cbind(approach = a, run_downstream(A$X, "dbscan", truth = t_a, seed = seed)))
  }))
  cbind(scenario = sc$id, replicate = r, seed = seed, res)
}

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
clusterExport(cl, c("base_seed", "one_check"))
cat(format(Sys.time()), "start\n", file = log_file, append = TRUE)
todo <- pilot_scenarios[pilot_scenarios$id <= 12, ]
for (k in todo$id) {
  sc <- todo[todo$id == k, ]
  f <- file.path(out_dir, sprintf("scenario_%02d.rds", k))
  if (file.exists(f)) next
  t0 <- Sys.time()
  saveRDS(do.call(rbind, parLapply(cl, seq_len(n_rep), function(r, sc) one_check(sc, r), sc = sc)), f)
  cat(format(Sys.time()), sc$label, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n",
      file = log_file, append = TRUE)
}
stopCluster(cl)

all <- do.call(rbind, lapply(todo$id, function(k) readRDS(file.path(out_dir, sprintf("scenario_%02d.rds", k)))))
all <- merge(all, todo[, c("id", "label", "family", "truth")], by.x = "scenario", by.y = "id")
summ <- do.call(rbind, lapply(split(all, list(all$scenario, all$approach, all$method, all$metric), drop = TRUE), function(d) {
  v <- d$value[is.finite(d$value)]
  data.frame(d[1, c("scenario", "label", "family", "truth", "approach", "method", "metric")], n = length(v),
             mean = if (length(v)) mean(v) else NA, sd = if (length(v) > 1) stats::sd(v) else NA,
             q05 = if (length(v)) stats::quantile(v, 0.05) else NA, q95 = if (length(v)) stats::quantile(v, 0.95) else NA,
             row.names = NULL)
}))
write.csv(summ, file.path("results", "pilot_hdbscan_summary.csv"), row.names = FALSE)
# paired contrasts against no correction, as in analysis/06_pilot.R
contr <- do.call(rbind, lapply(split(all, list(all$scenario, all$method, all$metric), drop = TRUE), function(d) {
  base <- d[d$approach == "none", c("replicate", "value")]
  do.call(rbind, lapply(setdiff(unique(d$approach), "none"), function(a) {
    x <- merge(d[d$approach == a, c("replicate", "value")], base, by = "replicate")
    diff <- (x$value.x - x$value.y)[is.finite(x$value.x - x$value.y)]
    if (length(diff) < 3) return(NULL)
    data.frame(d[1, c("scenario", "label", "family", "truth", "method", "metric")], approach = a, n = length(diff),
               effect = mean(diff), sd = stats::sd(diff), row.names = NULL)
  }))
}))
write.csv(contr, file.path("results", "pilot_hdbscan_contrasts.csv"), row.names = FALSE)
cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
