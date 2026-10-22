# Stage 5 — The correction on the simulation grid.
#
# Scenarios: the D4 grid (ids 1-10, the very data of stages 3 and 4) and the stage 5
# additions (ids 11-18): smooth large-scale fields (Matern 2.5, decision D12) and the
# "same regime, different fields" case (decision D11). Replicates: 50 per scenario
# (spatial, temporal, fields), 20 for ST 50x60 and 10 for ST 100x120 (cost).
#
# Approaches, for every replicate:
#   none          raw correlations
#   independent   one basis per variable, selected on that variable alone
#   global        one basis for all variables (the union of the bases selected on all
#                 dependent variables for each regime type present)
#   groups        the proposal: Gaussian-mixture groups (D6) and one basis per group
#   groups_split  the same plus the second level by subspace overlap (D11)
#   type_split    regime type as the first level plus the subspace level (D11)
#   oracle_groups the true groups, with bases selected on them
#   whiten        the D6 groups in whitening mode
# Error = RMS difference between the estimated and the true (sample) residual
# correlations over all pairs, by kind of pair (within a true group, across dependent
# groups, involving a variable without dependence); also against R_true. Capture: share
# of each variable's true large-scale component removed by the basis applied to it.
#
# Seeds: data 3300 + 1000 * scenario + replicate (the stage 3 seeds for ids 1-10);
# basis selection 5500 + 1000 * (10 + scenario) + replicate. The D11 threshold comes
# from results/d11_threshold.csv (analysis/05a_h2_mechanism.R).
#
# Outputs (not tracked by git): results/correction/scenario_NN.rds,
#   results/correction_replicates.rds, results/correction_summary.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

data_seed <- 3300L
stage_seed <- 5500L
out_dir <- file.path("results", "correction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")
thr_tab <- read.csv(file.path("results", "d11_threshold.csv"))
split_threshold <- thr_tab$threshold[which.max(thr_tab$youden)]
all_scenarios <- rbind(scenarios, extra_scenarios)
n_rep <- function(sc) if (sc$family != "spatiotemporal") 50L else if (sc$n == 50) 20L else 10L

# --- One replicate -----------------------------------------------------------------------
one_replicate <- function(sc, r, sig_d4) {
  seed <- data_seed + 1000L * sc$id + r
  bseed <- stage_seed + 1000L * (10L + sc$id) + r
  s <- simulate(sc, seed)
  sig <- if (!is.null(sig_d4)) sig_d4 else diagnose_variables(s)
  truth <- truth_table(s$groups)
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))
  p <- ncol(s$X)
  R_res <- stats::cor(s$components$residual)
  true_grp <- ifelse(truth$true_type == "none", "none", truth$group)
  kind <- outer(true_grp, true_grp, function(a, b) ifelse(a == "none" | b == "none", "involving_none",
                                                          ifelse(a == b, "within_group", "across_groups")))
  ut <- upper.tri(R_res)
  score <- function(R) {
    e <- (R - R_res)[ut]
    c(rms = sqrt(mean(e^2)), rms_vs_R_true = sqrt(mean((R - s$R_true)[ut]^2)),
      vapply(c("within_group", "across_groups", "involving_none"),
             function(k) if (any(kind[ut] == k)) sqrt(mean(e[kind[ut] == k]^2)) else NA_real_, 0))
  }
  capture_of <- function(basis_by_var) mean(vapply(which(truth$true_type != "none"), function(j)
    .capture_absorption(basis_by_var[[j]], s$components$large[, j, drop = FALSE], s$components$residual[, j, drop = FALSE])[["capture"]], 0))
  type <- unname(classify_regime(sig))
  per <- function(x) if (is.finite(x)) round(x) else 12
  n <- nrow(s$X)
  empty <- structure(matrix(0, n, 0), class = "depprep_basis", type = "none", k = 0, kn = 0, terms = character(0))
  timed <- function(expr) { t0 <- Sys.time(); v <- expr; list(value = v, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
  rows <- list()
  add <- function(name, R, seconds, capture = NA_real_, kn = NA_real_, autocor = NA_real_, n_groups = NA_real_) {
    rows[[name]] <<- data.frame(approach = name, t(score(R)), capture = capture, kn_mean = kn, autocor_after = autocor,
                                n_groups = n_groups, seconds = seconds)
  }
  add("none", stats::cor(s$X), 0)

  # independent: one basis per variable
  ind <- timed(lapply(seq_len(p), function(j) if (type[j] == "none") empty else
    do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j], period = per(sig$period[j]),
                                      seed = bseed + j), args))))
  E <- vapply(seq_len(p), function(j) .resid_on(s$X[, j, drop = FALSE], ind$value[[j]]), numeric(n))
  add("independent", stats::cor(E), ind$seconds, capture_of(ind$value), mean(vapply(ind$value, attr, 0, "kn")))

  # global: one basis for all variables
  glob <- timed({
    dep <- type != "none"
    bs <- lapply(unique(type[dep]), function(t) do.call(build_group_basis, c(list(X_group = s$X[, dep, drop = FALSE], type = t,
      period = per(stats::median(sig$period[dep], na.rm = TRUE)), seed = bseed), args)))
    if (length(bs)) suppressWarnings(do.call(union_basis, bs)) else empty
  })
  add("global", stats::cor(.resid_on(s$X, glob$value)), glob$seconds, capture_of(rep(list(glob$value), p)), attr(glob$value, "kn"))

  # grouped corrections
  grouped <- function(name, gr, mode = "conditional") {
    cr <- timed(correct(s, gr, mode = mode, scope = "groups", seed = bseed))
    v <- if (mode == "conditional") verify(cr$value, s) else NULL
    cap <- if (mode == "conditional") {
      g <- if (inherits(gr, "data.frame")) gr$group else gr
      capture_of(lapply(g, function(x) if (x == "none") empty else cr$value$bases[[x]]))
    } else NA_real_
    add(name, cr$value$R, cr$seconds, cap, if (length(cr$value$kn)) mean(cr$value$kn) else 0,
        if (!is.null(v)) stats::median(v$variables$autocor_after[type != "none"]) else NA_real_,
        length(setdiff(unique(if (inherits(gr, "data.frame")) gr$group else gr), "none")))
  }
  g6 <- timed(group_regimes(sig))
  grouped("groups", g6$value)
  grouped("groups_split", group_regimes(sig, data = s, split = "subspace", split_threshold = split_threshold, seed = bseed))
  grouped("type_split", group_regimes(sig, data = s, method = "type", split = "subspace", split_threshold = split_threshold, seed = bseed))
  oracle <- data.frame(variable = colnames(s$X), group = true_grp, group_type = ifelse(true_grp == "none", "none", truth$true_type))
  class(oracle) <- c("depprep_grouping", "data.frame")
  attr(oracle, "period") <- tapply(truth$true_period, true_grp, function(x) x[1])
  grouped("oracle_groups", oracle)
  grouped("whiten", g6$value, mode = "whiten")
  cbind(scenario = sc$id, replicate = r, seed = seed, do.call(rbind, rows), row.names = NULL)
}

# --- Run -------------------------------------------------------------------------------------
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE)); NULL }, getwd())
clusterExport(cl, c("data_seed", "stage_seed", "split_threshold", "spatial_field", "spatial_groups", "fields_groups",
                    "temporal_groups", "st_groups", "truth_table", "simulate", "one_replicate"))
cat(format(Sys.time()), "start; D11 threshold", split_threshold, "\n", file = log_file, append = TRUE)
for (k in seq_len(nrow(all_scenarios))) {
  sc <- all_scenarios[k, ]
  f <- file.path(out_dir, sprintf("scenario_%02d.rds", sc$id))
  if (file.exists(f)) next
  t0 <- Sys.time()
  d4 <- if (sc$id <= 10) readRDS(file.path("results", "diagnostics", sprintf("scenario_%02d.rds", sc$id))) else NULL
  tasks <- lapply(seq_len(n_rep(sc)), function(r) list(r = r, sig = if (!is.null(d4)) d4[d4$replicate == r, ] else NULL))
  res <- parLapply(cl, tasks, function(t, sc) one_replicate(sc, t$r, t$sig), sc = sc)
  saveRDS(do.call(rbind, res), f)
  cat(format(Sys.time()), sc$label, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n",
      file = log_file, append = TRUE)
}
stopCluster(cl)

# --- Summary -----------------------------------------------------------------------------------
all <- do.call(rbind, lapply(all_scenarios$id, function(k) readRDS(file.path(out_dir, sprintf("scenario_%02d.rds", k)))))
all$label <- all_scenarios$label[match(all$scenario, all_scenarios$id)]
all$family <- all_scenarios$family[match(all$scenario, all_scenarios$id)]
all$smoothness <- all_scenarios$smoothness[match(all$scenario, all_scenarios$id)]
saveRDS(all, file.path("results", "correction_replicates.rds"))
num_cols <- c("rms", "rms_vs_R_true", "within_group", "across_groups", "involving_none", "capture", "kn_mean",
              "autocor_after", "n_groups", "seconds")
summ <- aggregate(all[num_cols], all[c("scenario", "label", "family", "smoothness", "approach")],
                  function(x) mean(x, na.rm = TRUE))
summ <- summ[order(summ$scenario, match(summ$approach, c("none", "independent", "global", "groups", "groups_split",
                                                         "type_split", "oracle_groups", "whiten"))), ]
write.csv(summ, file.path("results", "correction_summary.csv"), row.names = FALSE)
cat(format(Sys.time()), "summary written\n", file = log_file, append = TRUE)
print(summ[, c("label", "approach", "rms", "within_group", "across_groups", "capture", "kn_mean")], digits = 2)
