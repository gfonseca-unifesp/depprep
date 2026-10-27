# Stage 7 — amendment protocol-v1.1 (registered on 2026-10-27, after the stage 7 results):
# H3-R and H2b recomputed against the correlation of the true signal.
#
# Why. Protocol-v1 scores the correlations against the sample correlation of the true
# residual components. In the truths with classes (structured, crossed, between 0.3),
# the class effects add to every variable a correlation that is not in the residuals, so
# every correction, and the raw data, look wrong by about the same amount and H3-R and H2b
# cannot be judged. The same flaw had been found and fixed for the PCA metric before the
# freeze (subspace_distance_signal); it was missed for the correlation error. The fix is
# the same: the truth becomes the correlation of the true signal, residual plus class
# effects (in the null truth the two coincide). The pre-registered results stay the main
# ones; these are reported beside them, labelled as an amendment.
#
# The data, groupings and corrections are those of analysis/07_full_grid.R (the same
# seeds: 8800 + 1000 * combination + replicate; results/stage7/grid.csv), for the
# combinations without perturbation and with the structured, crossed and between-0.3
# truths. Only one basis per variable (comparator 2) and the default of D7 (the proposal
# for R) are rebuilt; no output method is run.
#
# Outputs (not tracked): results/stage7_amendment/combo_NNN.rds,
#   results/stage7_amendment_tests.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 8800L
grid <- read.csv(file.path("results", "stage7", "grid.csv"), stringsAsFactors = FALSE)
todo <- grid[grid$perturbation == "none" & grid$truth_label %in% c("structured", "crossed", "between03"), ]
out_dir <- file.path("results", "stage7_amendment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

one_amend <- function(g, r) {
  seed <- base_seed + 1000L * g$combo + r                 # as in analysis/07_full_grid.R
  s <- simulate(g, seed, g$truth, between = g$between)
  truth_tab <- truth_table(s$groups)
  n <- nrow(s$X)
  latent <- if (!is.null(s$X_latent)) s$X_latent else s$X
  signal <- latent - sweep(s$components$large, 2, sqrt(truth_tab$share), "*")
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))
  sig <- diagnose_variables(s)
  type <- unname(classify_regime(sig))
  gr <- group_regimes(sig)
  E <- vapply(seq_len(ncol(s$X)), function(j) if (type[j] == "none") s$X[, j] - mean(s$X[, j]) else
    .resid_on(s$X[, j, drop = FALSE], do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j],
      period = if (is.finite(sig$period[j])) round(sig$period[j]) else 12, seed = seed + j), args))), numeric(n))
  R_def <- correct(s, gr, seed = seed, data_basis = "union")$R
  R_sig <- stats::cor(signal)
  R_res <- stats::cor(s$components$residual)
  ut <- upper.tri(R_sig)
  rms <- function(R, T) sqrt(mean((R - T)[ut]^2))
  data.frame(combo = g$combo, replicate = r, seed = seed,
             approach = c("none", "independent", "default"),
             rms_signal = c(rms(stats::cor(s$X), R_sig), rms(stats::cor(E), R_sig), rms(R_def, R_sig)),
             rms = c(rms(stats::cor(s$X), R_res), rms(stats::cor(E), R_res), rms(R_def, R_res)))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "one_amend"))
  cat(format(Sys.time()), "start;", nrow(todo), "combinations\n", file = log_file, append = TRUE)
  for (k in todo$combo) {
    g <- todo[todo$combo == k, ]
    f <- file.path(out_dir, sprintf("combo_%03d.rds", k))
    if (file.exists(f)) next
    t0 <- Sys.time()
    reps <- if (g$base_id %in% c(10, 16)) seq_len(utils::read.csv(file.path("results", "stage7", "timing_st100.csv"))$replicates) else
      seq_len(g$n_rep)
    saveRDS(do.call(rbind, parLapply(cl, reps, function(r, g) one_amend(g, r), g = g)), f)
    cat(format(Sys.time()), sprintf("combo %d (%s %s) done in %.1f min\n", k, g$label, g$truth_label,
                                    as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_file, append = TRUE)
  }
  stopCluster(cl)

  # tests as registered in protocol-v1 (section 2.3), with the signal-based truth
  am <- do.call(rbind, lapply(todo$combo, function(k) readRDS(file.path(out_dir, sprintf("combo_%03d.rds", k)))))
  am <- merge(am, transform(todo[, c("combo", "base_id", "family")], data_family = ifelse(family == "fields", "spatial", family)), by = "combo")
  strat <- function(d, sign = 1) { by <- split(sign * d$diff, d$combo)
    c(estimate = mean(vapply(by, mean, 0)), se = sqrt(sum(vapply(by, function(x) stats::var(x) / length(x), 0))) / length(by)) }
  rows <- list()
  for (f in c("spatial", "temporal", "spatiotemporal")) {
    d <- am[am$data_family == f, ]
    w <- reshape(d[, c("combo", "replicate", "approach", "rms_signal")], idvar = c("combo", "replicate"), timevar = "approach",
                 direction = "wide")
    w$diff <- w$rms_signal.independent - w$rms_signal.default
    e <- strat(w)
    rows[[length(rows) + 1]] <- data.frame(hypothesis = "H3-R (v1.1)", data_family = f, contrast = "independent - default",
      prediction = if (f == "spatiotemporal") "+" else "exploratory", delta = 0.10 * mean(w$rms_signal.independent),
      estimate = e[["estimate"]], se = e[["se"]])
  }
  rs <- list(spatial = list(rough = c(1, 2, 3, 4, 17), smooth = c(11, 12, 13, 14, 18)),
             spatiotemporal = list(rough = c(9, 10), smooth = c(15, 16)))
  for (f in names(rs)) {
    d <- am[am$approach == "default", ]
    r <- d[d$base_id %in% rs[[f]]$rough, ]; s <- d[d$base_id %in% rs[[f]]$smooth, ]
    mr <- vapply(split(r$rms_signal, r$combo), mean, 0); ms <- vapply(split(s$rms_signal, s$combo), mean, 0)
    vr <- vapply(split(r$rms_signal, r$combo), function(x) stats::var(x) / length(x), 0)
    vs <- vapply(split(s$rms_signal, s$combo), function(x) stats::var(x) / length(x), 0)
    rows[[length(rows) + 1]] <- data.frame(hypothesis = "H2b (v1.1)", data_family = f, contrast = "rough - smooth",
      prediction = "+", delta = 0.05, estimate = mean(mr) - mean(ms), se = sqrt(sum(vr) / length(mr)^2 + sum(vs) / length(ms)^2))
  }
  tests <- do.call(rbind, rows)
  tests$p <- 2 * stats::pnorm(abs(tests$estimate) / tests$se, lower.tail = FALSE)
  tests$p_holm <- stats::ave(tests$p, tests$hypothesis, FUN = function(p) stats::p.adjust(p, "holm"))
  tests$confirmed <- ifelse(tests$prediction == "exploratory", NA, tests$estimate >= tests$delta & tests$p_holm < 0.05)
  means <- stats::aggregate(cbind(rms_signal, rms) ~ data_family + approach, data = am, FUN = mean)
  write.csv(tests, file.path("results", "stage7_amendment_tests.csv"), row.names = FALSE)
  write.csv(means, file.path("results", "stage7_amendment_means.csv"), row.names = FALSE)
  cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
}
