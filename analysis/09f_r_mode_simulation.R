# Stage 9 — which mode best estimates R in space-time when the association between variables
# is itself autocorrelated and each variable has noise of its own (protocol-v1.3.1, section
# 17.3; decision D29). The rule of 17.3 gives a recommendation; the decision is the maintainer's.
#
# Grid: the groups of scenario 9 at 100 sites x 93 months; association white (the conditions of
# D7), spatial or spatio-temporal; noise share 0 or 0.5; truth null or structured (blocks, 0.5
# within groups and 0.3 across them, no classes). 12 cells, 20 replicates.
# Modes, with the true and with the estimated grouping: whitening with the groups (D7);
# conditional with the groups and the union basis; conditional + AR (the residuals of the
# conditional mode, then an AR filter of each site with the coefficients of the group).
# Truth: the sample correlation of the residual components.
# The diagnosis uses 499 permutations, as in 09d: in space-time the spatial p-value is the
# smallest over 6 time slices times 6, so with 99 it could never fall below 0.05.
#
# Usage: Rscript analysis/09f_r_mode_simulation.R          every replicate, resumable
#        Rscript analysis/09f_r_mode_simulation.R pilot    one replicate of two cells, timed
#
# Seed: 18000 + 1000 * cell + replicate, recorded in the output.
# Outputs (not tracked): results/r_mode_simulation/ (one file per replicate),
# results/r_mode_simulation.csv, results/r_mode_summary.csv, results/r_mode_decision.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))
base_seed <- 18000L
set.seed(base_seed)
n_sites <- 100L; n_times <- 93L; n_rep <- 20L
pilot <- identical(commandArgs(trailingOnly = TRUE)[1], "pilot")
grid <- expand.grid(association = c("white", "spatial", "spatio-temporal"), noise_share = c(0, 0.5),
                    truth = c("null", "structured"), stringsAsFactors = FALSE)
grid$cell <- seq_len(nrow(grid))
range_of <- c(white = 0, spatial = 0.2, "spatio-temporal" = 0.2)
phi_of <- c(white = 0, spatial = 0, "spatio-temporal" = 0.6)
regime_of <- c(sep_long = "separable_ar1", sep_short = "separable_seasonal", moving = "moving",
               common = "temporal_common", none = "none")
out_dir <- file.path("results", "r_mode_simulation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The residuals of the conditional mode, then the AR filter of each site, pair of groups by
# pair of groups as in correct().
conditional_ar <- function(s, cc) {
  X <- s$X; grp <- cc$group; ug <- unique(grp)
  n_s <- nrow(s$coords); n_t <- max(s$index$time)
  phi <- stats::setNames(lapply(ug, function(g) {
    if (g == "none") return(NULL)
    environment(.whitener(X[, grp == g, drop = FALSE], g, "", NULL, "spatiotemporal", s$coords, NULL, s$index))$phi
  }), ug)
  filter_time <- function(E, g) {
    if (is.null(phi[[g]])) return(E)
    p <- length(phi[[g]])
    do.call(cbind, lapply(seq_len(ncol(E)), function(j) {
      M <- t(matrix(E[, j], n_s, n_t))
      out <- M[-(1:p), , drop = FALSE]
      for (i in seq_len(p)) out <- out - phi[[g]][i] * M[(p + 1 - i):(n_t - i), , drop = FALSE]
      as.vector(t(out))
    }))
  }
  basis_of <- function(g) if (g == "none") NULL else cc$bases[[g]]
  R <- diag(ncol(X)); dimnames(R) <- list(colnames(X), colnames(X))
  for (a in seq_along(ug)) for (b in a:length(ug)) {
    ga <- ug[a]; gb <- ug[b]; ia <- which(grp == ga); ib <- which(grp == gb)
    U <- if (ga == gb) basis_of(ga) else if (ga == "none") basis_of(gb) else if (gb == "none") basis_of(ga) else
      suppressWarnings(union_basis(basis_of(ga), basis_of(gb), kn_max = 1))
    res <- function(i) if (is.null(U)) scale(X[, i, drop = FALSE], scale = FALSE) else .resid_on(X[, i, drop = FALSE], U)
    Ea <- filter_time(res(ia), ga); Eb <- filter_time(res(ib), gb)
    m <- min(nrow(Ea), nrow(Eb))
    R[ia, ib] <- stats::cor(utils::tail(Ea, m), utils::tail(Eb, m))
    R[ib, ia] <- t(R[ia, ib])
  }
  diag(R) <- 1
  near <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  dimnames(near) <- dimnames(R)
  near
}

# The measures of 17.3 for one estimate of R, over all pairs and by class of pair.
measures <- function(R, truth, group) {
  up <- upper.tri(truth); p <- nrow(truth)
  pair_class <- ifelse(outer(group, group, "=="), paste("within", regime_of[group])[row(truth)], "across")
  classes <- c("all", unique(pair_class[up]))
  do.call(rbind, lapply(classes, function(k) {
    w <- up & (k == "all" | pair_class == k)
    strong <- w & abs(truth) >= 0.2
    data.frame(class = k, pairs = sum(w), rms = sqrt(mean((R - truth)[w]^2)), strong_pairs = sum(strong),
               mean_abs_true = if (any(strong)) mean(abs(truth[strong])) else NA_real_,
               mean_abs_r = if (any(strong)) mean(abs(R[strong])) else NA_real_,
               sign_kept = if (any(strong)) mean(sign(R[strong]) == sign(truth[strong]) & abs(R[strong]) >= 0.2) else NA_real_)
  }))
}

one_replicate <- function(k, r) {
  cell <- grid[grid$cell == k, ]
  seed <- base_seed + 1000L * k + r
  g <- st_groups()
  R_true <- if (cell$truth == "structured") make_R_true(g, "blocks", within = 0.5, between = 0.3) else NULL
  s <- sim_spatiotemporal(n_sites, n_times, g, R_true = R_true, seed = seed, residual_range = range_of[[cell$association]],
                          residual_phi = phi_of[[cell$association]], noise_share = cell$noise_share)
  truth <- stats::cor(s$components$residual)
  tt <- truth_table(s$groups)
  true_group <- attr(g, "var_group")
  oracle <- data.frame(variable = colnames(s$X), group = true_group, group_type = ifelse(true_group == "none", "none", tt$true_type))
  class(oracle) <- c("depprep_grouping", "data.frame")
  attr(oracle, "period") <- tapply(tt$true_period, true_group, function(x) x[1])
  t0 <- Sys.time()
  estimated <- group_regimes(diagnose_variables(s, nsim = 499, perm_seed = seed), seed = seed)
  time_grouping <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  types_match <- mean(estimated$type == oracle$group_type)
  rows <- list()
  for (grouping in c("true", "estimated")) {
    gr <- if (grouping == "true") oracle else estimated
    t1 <- Sys.time()
    cc <- correct(s, gr, mode = "conditional", scope = "groups", data_basis = "union", seed = seed)
    t2 <- Sys.time()
    cw <- correct(s, gr, bases = cc$bases, mode = "whiten", scope = "groups", data_basis = "union", seed = seed)
    t3 <- Sys.time()
    ca <- conditional_ar(s, cc)
    t4 <- Sys.time()
    est <- list(whitening = cw$R, conditional = cc$R, conditional_ar = ca)
    seconds <- c(whitening = as.numeric(difftime(t3, t2, units = "secs")), conditional = as.numeric(difftime(t2, t1, units = "secs")),
                 conditional_ar = as.numeric(difftime(t4, t3, units = "secs")))
    for (m in names(est)) {
      first <- eigen(est[[m]], symmetric = TRUE, only.values = TRUE)$values[1] / ncol(truth)
      rows[[paste(grouping, m)]] <- cbind(cell[rep(1, 1), ], replicate = r, seed = seed, grouping = grouping, mode = m,
                                          measures(est[[m]], truth, true_group), first_axis_share = first,
                                          first_above_broken_stick = first > sum(1 / seq_len(ncol(truth))) / ncol(truth),
                                          seconds = seconds[[m]], seconds_grouping = time_grouping, types_match = types_match)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

# One replicate, resumable; the warnings of R are kept with its rows.
one <- function(k, r) {
  file <- file.path(out_dir, sprintf("cell_%02d_rep_%02d.rds", k, r))
  if (file.exists(file)) return(readRDS(file))
  warned <- character(0)
  out <- withCallingHandlers(one_replicate(k, r), warning = function(w) {
    warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
  out$warnings <- paste(unique(warned), collapse = " | ")
  saveRDS(out, file)
  out
}

tasks <- if (pilot) data.frame(k = c(1L, 12L), r = 1L) else expand.grid(r = seq_len(n_rep), k = grid$cell)
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
t_start <- Sys.time()
if (pilot) {
  res <- lapply(seq_len(nrow(tasks)), function(i) one(tasks$k[i], tasks$r[i]))
  print(do.call(rbind, res)[do.call(rbind, res)$class == "all", c("cell", "grouping", "mode", "rms", "seconds", "seconds_grouping")], digits = 3)
  cat("minutes:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "\n")
  unlink(file.path(out_dir, sprintf("cell_%02d_rep_%02d.rds", tasks$k, tasks$r)))
  quit(save = "no")
}
cl <- makeCluster(max(1L, min(12L, detectCores() - 4L)), setup_timeout = 900)
invisible(clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd()))
clusterExport(cl, c("one", "one_replicate", "grid", "base_seed", "n_sites", "n_times", "range_of", "phi_of", "regime_of", "out_dir", "conditional_ar", "measures"))
res <- clusterApplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(t) one(t$k, t$r))
stopCluster(cl)
all_rows <- do.call(rbind, res)
utils::write.csv(all_rows, file.path("results", "r_mode_simulation.csv"), row.names = FALSE)

summ <- aggregate(cbind(rms, mean_abs_true, mean_abs_r, sign_kept, first_above_broken_stick, seconds, types_match) ~
                    association + noise_share + truth + grouping + mode + class,
                  data = all_rows, FUN = mean, na.rm = TRUE, na.action = stats::na.pass)
utils::write.csv(summ, file.path("results", "r_mode_summary.csv"), row.names = FALSE)

# The rule of 17.3: test cells, estimated grouping, over all pairs.
a <- all_rows[all_rows$class == "all" & all_rows$grouping == "estimated", ]
test <- a[a$truth == "structured" & a$association != "white" & a$noise_share == 0.5, ]
err <- tapply(test$rms, test$mode, mean)
nulls <- a[a$truth == "null", ]
spurious <- tapply(nulls$first_above_broken_stick, nulls$mode, mean)
best_other <- min(err[c("conditional", "conditional_ar")])
recommended <- if (err[["whitening"]] <= 1.1 * best_other) "whitening" else
  if (err[["conditional"]] <= 1.1 * err[["conditional_ar"]]) "conditional" else "conditional_ar"
guard <- spurious[[recommended]] <= spurious[["whitening"]] + 0.10
regime <- all_rows[all_rows$class %in% paste("within", regime_of[1:4]) & all_rows$grouping == "estimated" &
                     all_rows$truth == "structured" & all_rows$association != "white" & all_rows$noise_share == 0.5, ]
by_regime <- tapply(regime$rms, list(regime$class, regime$mode), mean)
best <- apply(by_regime, 1, function(x) names(x)[which.min(x)])
decision <- data.frame(error_whitening = err[["whitening"]], error_conditional = err[["conditional"]], error_conditional_ar = err[["conditional_ar"]],
                       spurious_whitening = spurious[["whitening"]], spurious_conditional = spurious[["conditional"]],
                       spurious_conditional_ar = spurious[["conditional_ar"]], recommended = if (guard) recommended else "none (guard)",
                       guard_met = guard, best_by_regime = paste(names(best), best, sep = ": ", collapse = "; "), base_seed = base_seed)
utils::write.csv(decision, file.path("results", "r_mode_decision.csv"), row.names = FALSE)
print(decision)
print(by_regime, digits = 3)
cat("minutes:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "\n")
