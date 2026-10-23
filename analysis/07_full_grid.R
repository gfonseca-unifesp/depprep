# Stage 7 — the full simulation grid of protocol-v1 (interno/docs/protocol.md, frozen on
# 2026-10-23 with the tag protocol-v1). Nothing here may differ from the protocol.
#
# Grid (protocol section 4): the 21 base scenarios (analysis/d4_scenarios.R: D4 ids 1-10,
# stage 5 ids 11-18, D14 ids 19-21) x five truths (null; structured; contiguous; crossed;
# structured with 0.3 between groups). Space-time 100 x 120 (ids 10 and 16): null and
# structured only. Perturbations, structured truth: counts (one third of the variables) in
# ids 2, 6 and 9; non-stationary range (0.5 to 1.5 times) in id 2.
# Replicates: 100 (spatial, temporal), 50 (space-time 50 x 60), 10 (space-time
# 100 x 120; 5 when the first two replicates show the run would exceed about 24 hours).
# Seeds: 8800 + 1000 * combination + replicate (the combination number is fixed below).
#
# Approaches (section 5): none; independent (one basis per variable: the proposal for the
# data matrix, D17, and comparator 2 for R); default (D7: global basis in space and time,
# whitening in space-time; its R is the proposal for R, its data matrix on the union of the
# bases is the shared-basis comparator); groups (D6); variables (union of the variable
# bases, D14); whiten (space and time); sbss and maf (spatial); broad (contiguous truth,
# space and time).
# Output methods (section 6): PCA on R; MDS (Euclidean, Gower), k-means, SOM, HDBSCAN and
# DBSCAN on the data matrix (constrained clustering too in the contiguous truth), on a
# seeded subsample above 3,000 observations (2,000 in space-time 100 x 120).
# Also per replicate: the error of each R against the true residual correlations and
# against R_true; Delta R of the proposal; Delta R between cor(X) of the proposal and its R;
# the Q3 statistic and warning of verify(); k/n of the unions and whether the union of the
# variable bases was capped.
#
# Runs in parallel and resumes: one file per combination in results/stage7/.
# Outputs (not tracked): results/stage7/combo_NNN.rds, results/stage7/grid.csv,
#   results/stage7/progress.log, results/stage7/timing_st100.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 8800L
out_dir <- file.path("results", "stage7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

# --- Grid ------------------------------------------------------------------------------------
bases <- rbind(scenarios, extra_scenarios, d14_scenarios)
bases <- bases[, setdiff(names(bases), "between")]        # the truth sets it
truths <- data.frame(truth = c("null", "structured", "contiguous", "crossed", "structured"),
                     between = c(0, 0, 0, 0, 0.3), truth_label = c("null", "structured", "contiguous", "crossed", "between03"),
                     stringsAsFactors = FALSE)
big <- bases$id[bases$family == "spatiotemporal" & bases$n == 100]
grid <- merge(bases, truths)
grid <- grid[!(grid$id %in% big & !grid$truth_label %in% c("null", "structured")), ]
grid$perturbation <- "none"
pert <- function(id, what) cbind(bases[bases$id == id, ], truth = "structured", between = 0, truth_label = "structured",
                                 perturbation = what, stringsAsFactors = FALSE)
grid <- rbind(grid, pert(2, "counts"), pert(6, "counts"), pert(9, "counts"), pert(2, "nonstationary"))
fam_order <- match(ifelse(grid$family == "fields", "spatial", grid$family), c("spatial", "temporal", "spatiotemporal"))
grid <- grid[order(grid$id %in% big, grid$perturbation != "none", fam_order, grid$id,
                   match(grid$truth_label, truths$truth_label)), ]
grid$base_id <- grid$id
grid$combo <- seq_len(nrow(grid))
grid$n_rep <- ifelse(grid$id %in% big, 10L, ifelse(grid$family == "spatiotemporal", 50L, 100L))
grid$max_n <- ifelse(grid$id %in% big, 2000L, 3000L)
write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)

# --- One replicate -----------------------------------------------------------------------------
one_rep <- function(g, r) {
  seed <- base_seed + 1000L * g$combo + r
  t_rep <- Sys.time()
  s <- simulate(g, seed, g$truth, between = g$between, counts = if (g$perturbation == "counts") 1 / 3 else 0,
                nonstationary = if (g$perturbation == "nonstationary") c(0.5, 1.5) else NULL)
  fam <- s$family
  truth_tab <- truth_table(s$groups)
  n <- nrow(s$X)
  latent <- if (!is.null(s$X_latent)) s$X_latent else s$X
  large_part <- sweep(s$components$large, 2, sqrt(truth_tab$share), "*")
  position <- switch(fam, spatial = s$coords, temporal = matrix(s$time),
                     spatiotemporal = cbind(s$coords[s$index$site, ], s$index$time / max(s$index$time)))
  truth_all <- list(R_true = s$R_true, classes = if (g$truth == "null") NULL else s$classes,
                    signal = latent - large_part, position = position)
  rows_tail <- function(m) { keep <- utils::tail(seq_len(n), m)
    list(R_true = truth_all$R_true, classes = truth_all$classes[keep], signal = truth_all$signal[keep, , drop = FALSE],
         position = truth_all$position[keep, , drop = FALSE]) }
  args <- switch(fam, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))
  sig <- diagnose_variables(s)
  type <- unname(classify_regime(sig))
  gr <- group_regimes(sig)
  ap <- list(none = list(X = s$X, R = stats::cor(s$X)))
  E <- vapply(seq_len(ncol(s$X)), function(j) if (type[j] == "none") s$X[, j] - mean(s$X[, j]) else
    .resid_on(s$X[, j, drop = FALSE], do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j],
      period = if (is.finite(sig$period[j])) round(sig$period[j]) else 12, seed = seed + j), args))), numeric(n))
  ap$independent <- list(X = E, R = stats::cor(E))
  cr_def <- correct(s, gr, seed = seed, data_basis = "union")
  ap$default <- cr_def[c("X", "R")]
  cr_grp <- correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union")
  ap$groups <- cr_grp[c("X", "R")]
  cr_var <- correct(s, gr, mode = "conditional", scope = "variables", seed = seed, data_basis = "union")
  ap$variables <- cr_var[c("X", "R")]
  if (fam != "spatiotemporal") ap$whiten <- correct(s, gr, mode = "whiten", seed = seed, data_basis = "union")[c("X", "R")]
  if (fam == "spatial") {
    ap$sbss <- decorrelate_latent(s$X, s$coords, "sbss", seed = seed)[c("X", "R")]
    ap$maf <- decorrelate_latent(s$X, s$coords, "maf", seed = seed)[c("X", "R")]
  }
  if (g$truth_label == "contiguous" && fam != "spatiotemporal") {
    cand <- if (fam == "spatial") .spatial_candidates(s$coords) else .temporal_candidates(s$time, "ar")
    band <- cand[, seq_len(ceiling(ncol(cand) / 3)), drop = FALSE]
    P <- band %*% crossprod(band, scale(s$X, scale = FALSE))
    ap$broad <- list(X = P, R = stats::cor(P))
  }
  methods <- function(A, tr) rbind(
    run_downstream(A$X, "pca", truth = tr, R = A$R),
    run_downstream(A$X, "mds", truth = tr, distance = "euclidean", seed = seed, max_n = g$max_n),
    run_downstream(A$X, "mds", truth = tr, distance = "gower", seed = seed, max_n = g$max_n),
    run_downstream(A$X, "kmeans", truth = tr, seed = seed, max_n = g$max_n),
    run_downstream(A$X, "som", truth = tr, seed = seed, max_n = g$max_n),
    run_downstream(A$X, "hdbscan", truth = tr, seed = seed, max_n = g$max_n),
    run_downstream(A$X, "dbscan", truth = tr, seed = seed, max_n = g$max_n),
    if (g$truth_label == "contiguous") run_downstream(A$X, "constrained", truth = tr, seed = seed, max_n = g$max_n))
  R_res <- stats::cor(s$components$residual)
  ut <- upper.tri(R_res)
  res <- do.call(rbind, lapply(names(ap), function(a) {
    A <- ap[[a]]
    tr <- if (nrow(A$X) < n) rows_tail(nrow(A$X)) else truth_all
    err <- data.frame(method = "correlation", metric = c("rms", "rms_vs_R_true"),
                      value = c(sqrt(mean((A$R - R_res)[ut]^2)), sqrt(mean((A$R - s$R_true)[ut]^2))))
    cbind(approach = a, rbind(methods(A, tr), err))
  }))
  # Delta R of the proposal for R (section 6) and between cor(X) of the proposal for the
  # data matrix and that R
  d_raw <- abs(stats::cor(s$X) - s$R_true)[ut]; d_cor <- abs(cr_def$R - s$R_true)[ut]
  extra <- data.frame(approach = "default", method = "delta_R",
                      metric = c("correctly_removed", "wrongly_changed", "mean_abs_cor_X_minus_R"),
                      value = c(mean(d_raw - d_cor >= 0.05), mean(d_cor - d_raw >= 0.05),
                                mean(abs(stats::cor(E) - cr_def$R)[ut])))
  # Q3 warning (verify(); space-time: the group bases, as calibrated) and k/n of the unions
  q3_fit <- if (fam == "spatiotemporal") cr_grp else cr_def
  v <- verify(q3_fit, s, nperm = 0, seed = seed)
  kn_def <- if (fam == "spatiotemporal") NA_real_ else attr(cr_def$bases[[1]], "kn")   # the global basis
  extra <- rbind(extra, data.frame(approach = c("default", "default", "default", "variables", "variables"),
    method = c("verify", "verify", "kn", "kn", "kn"),
    metric = c("q3_statistic", "q3_warning", "union_kn", "union_kn", "union_capped"),
    value = c(v$q3$statistic, as.numeric(v$q3$warning), if (is.null(kn_def)) NA_real_ else kn_def,
              attr(cr_var$bases$union, "kn"), as.numeric(any(grepl("capped", cr_var$warnings))))))
  out <- rbind(res, extra)
  cbind(combo = g$combo, replicate = r, seed = seed, out,
        seconds_replicate = as.numeric(difftime(Sys.time(), t_rep, units = "secs")))
}

# --- Run -----------------------------------------------------------------------------------------
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
n_workers <- max(1, min(24, detectCores() - 4))
cl <- makeCluster(n_workers, setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
clusterExport(cl, c("base_seed", "one_rep"))
cat(format(Sys.time()), "start;", nrow(grid), "combinations\n", file = log_file, append = TRUE)
run_reps <- function(g, reps) do.call(rbind, parLapply(cl, reps, function(r, g) one_rep(g, r), g = g))
timing_file <- file.path(out_dir, "timing_st100.csv")
for (k in grid$combo) {
  g <- grid[grid$combo == k, ]
  f <- file.path(out_dir, sprintf("combo_%03d.rds", k))
  if (file.exists(f)) next
  t0 <- Sys.time()
  if (g$base_id %in% big) {
    # space-time 100 x 120: time two replicates first (protocol section 4)
    if (!file.exists(timing_file)) {
      first <- run_reps(g, 1:2)
      sec <- mean(unique(first[, c("replicate", "seconds_replicate")])$seconds_replicate)
      tasks <- sum(grid$base_id %in% big) * 10L
      hours <- ceiling(tasks / n_workers) * sec / 3600
      reps_big <- if (hours > 24) 5L else 10L
      write.csv(data.frame(seconds_per_replicate = sec, workers = n_workers, estimated_hours = hours, replicates = reps_big),
                timing_file, row.names = FALSE)
      cat(format(Sys.time()), sprintf("space-time 100 x 120 timed: %.0f s per replicate, about %.1f h -> %d replicates\n",
                                      sec, hours, reps_big), file = log_file, append = TRUE)
      saveRDS(first, file.path(out_dir, sprintf("combo_%03d_first.rds", k)))
    }
    reps_big <- utils::read.csv(timing_file)$replicates
    ff <- file.path(out_dir, sprintf("combo_%03d_first.rds", k))
    done <- if (file.exists(ff)) readRDS(ff) else NULL
    todo <- setdiff(seq_len(reps_big), if (is.null(done)) integer(0) else unique(done$replicate))
    res <- rbind(done, if (length(todo)) run_reps(g, todo))
  } else {
    res <- run_reps(g, seq_len(g$n_rep))
  }
  saveRDS(res, f)
  cat(format(Sys.time()), sprintf("combo %d (%s %s, %s) done in %.1f min\n", k, g$label, g$truth_label, g$perturbation,
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_file, append = TRUE)
}
stopCluster(cl)
cat(format(Sys.time()), "all combinations done\n", file = log_file, append = TRUE)
