# Stage 7 — amendment protocol-v1.2 (registered on 2026-10-27, after the stage 7 results):
# a threshold of the Q3 warning for each data family.
#
# Why. The Q3 warning of verify() missed its registered criterion in stage 7 (at most 5% of
# warnings in every truth without contiguous classes, at least 90% in the contiguous one):
# with the single threshold of 0.14 it fired in 9.5-26% of the temporal replicates without
# contiguous classes and in 6.5% of the spatial null ones. Within each family the
# statistic separates the truths, but at different levels: the temporal bases are larger
# and the statistic grows with them. The pre-registered result stays "criterion not met".
#
# Two steps, run in this order (Rscript analysis/07d_q3_recalibration.R calibrate|validate):
#
# calibrate: from the stage 7 statistics (results/stage7/combo_*.rds; no perturbation,
#   gradient id 19 left out), one threshold per family: the midpoint between the largest
#   97.5th percentile of the truths without contiguous classes (null, structured, crossed,
#   between 0.3) and the 5th percentile of the contiguous truth, rounded to 2 decimals; when
#   the two cross, the value that maximises the true positive rate minus the false positive
#   rate (Youden), as in analysis/06c_q3_warning.R.
# validate: new seeds, 9900 + 1000 * combination + replicate (combinations of
#   results/stage7/grid.csv without perturbation; the replicates of stage 7). The same path
#   as analysis/07_full_grid.R up to the statistic: simulate, diagnose, group, the default
#   correction of D7 (space-time: the group bases, as calibrated), verify(nperm = 0). The
#   registered criterion is applied with the threshold of each family, by family and truth
#   with the scenarios pooled; the rates per scenario, those of the gradient and those of
#   the old threshold (0.14) are reported apart.
#
# Outputs (not tracked): results/q3_thresholds_v12.csv, results/stage7_q3v/combo_NNN.rds,
#   results/stage7_q3v/progress.log, results/stage7_q3v_tests.csv,
#   results/stage7_q3v_scenarios.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

grid <- read.csv(file.path("results", "stage7", "grid.csv"), stringsAsFactors = FALSE)
grid <- grid[grid$perturbation == "none", ]
grid$data_family <- ifelse(grid$family == "fields", "spatial", grid$family)
grid$gradient <- grid$base_id == 19
families <- c("spatial", "temporal", "spatiotemporal")
old_threshold <- 0.14
threshold_file <- file.path("results", "q3_thresholds_v12.csv")

# --- Calibration on the stage 7 statistics ----------------------------------------------------------
calibrate_family <- function(d) {
  neg <- d$value[d$truth_label != "contiguous"]
  neg_q <- max(vapply(split(neg, d$truth_label[d$truth_label != "contiguous"]),
                      function(v) unname(stats::quantile(v, 0.975, type = 1)), 0))
  pos <- d$value[d$truth_label == "contiguous"]
  pos_q <- unname(stats::quantile(pos, 0.05, type = 1))
  if (neg_q < pos_q) return(data.frame(negative_q975 = neg_q, contiguous_q05 = pos_q, threshold = round((neg_q + pos_q) / 2, 2),
                                       rule = "midpoint"))
  cand <- sort(unique(c(neg, pos)))
  j <- vapply(cand, function(t) mean(pos > t) - mean(neg > t), 0)
  data.frame(negative_q975 = neg_q, contiguous_q05 = pos_q, threshold = round(cand[which.max(j)], 2), rule = "Youden")
}

calibrate <- function() {
  st7 <- do.call(rbind, lapply(grid$combo, function(k) {
    d <- readRDS(file.path("results", "stage7", sprintf("combo_%03d.rds", k)))
    d[d$method == "verify" & d$metric == "q3_statistic", c("combo", "replicate", "value")]
  }))
  st7 <- merge(st7, grid[, c("combo", "data_family", "truth_label", "gradient")], by = "combo")
  st7 <- st7[!st7$gradient, ]
  out <- do.call(rbind, lapply(families, function(f) cbind(data_family = f, calibrate_family(st7[st7$data_family == f, ]),
                                                         replicates = sum(st7$data_family == f))))
  write.csv(out, threshold_file, row.names = FALSE)
  print(out, digits = 3, row.names = FALSE)
}

# --- Validation on new seeds -----------------------------------------------------------------------
base_seed <- 9900L
out_dir <- file.path("results", "stage7_q3v")

one_q3v <- function(g, r) {
  seed <- base_seed + 1000L * g$combo + r
  s <- simulate(g, seed, g$truth, between = g$between)
  gr <- group_regimes(diagnose_variables(s))
  fit <- if (s$family == "spatiotemporal") correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union") else
    correct(s, gr, seed = seed, data_basis = "union")
  v <- verify(fit, s, nperm = 0, seed = seed)
  kn <- if (!length(fit$bases)) NA_real_ else ncol(.union_of_bases(fit, nrow(s$X))) / nrow(s$X)
  data.frame(combo = g$combo, replicate = r, seed = seed, statistic = v$q3$statistic, kn = kn)
}

validate <- function() {
  thr <- utils::read.csv(threshold_file)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "progress.log")
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "one_q3v"), envir = environment(one_q3v))
  reps_big <- utils::read.csv(file.path("results", "stage7", "timing_st100.csv"))$replicates
  todo <- grid[order(grid$family != "spatiotemporal" | grid$n != 100), ]          # the slow ones first
  cat(format(Sys.time()), "start;", nrow(todo), "combinations\n", file = log_file, append = TRUE)
  for (k in todo$combo) {
    g <- todo[todo$combo == k, ]
    f <- file.path(out_dir, sprintf("combo_%03d.rds", k))
    if (file.exists(f)) next
    t0 <- Sys.time()
    reps <- seq_len(if (g$base_id %in% c(10, 16)) reps_big else g$n_rep)
    saveRDS(do.call(rbind, parLapply(cl, reps, function(r, g) one_q3v(g, r), g = g)), f)
    cat(format(Sys.time()), sprintf("combo %d (%s %s) done in %.1f min\n", k, g$label, g$truth_label,
                                    as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_file, append = TRUE)
  }
  stopCluster(cl)

  q <- do.call(rbind, lapply(grid$combo, function(k) readRDS(file.path(out_dir, sprintf("combo_%03d.rds", k)))))
  q <- merge(q, grid[, c("combo", "label", "base_id", "data_family", "truth_label", "gradient")], by = "combo")
  q$threshold <- thr$threshold[match(q$data_family, thr$data_family)]
  q$warning <- q$statistic > q$threshold
  q$warning_old <- q$statistic > old_threshold
  rate <- function(d) data.frame(replicates = nrow(d), warning_rate = mean(d$warning), warning_rate_old = mean(d$warning_old),
                                 statistic_mean = mean(d$statistic), kn_mean = mean(d$kn, na.rm = TRUE))
  tests <- do.call(rbind, lapply(split(q, list(q$data_family, q$truth_label, q$gradient), drop = TRUE), function(d)
    cbind(d[1, c("data_family", "truth_label", "gradient", "threshold")], rate(d))))
  tests$criterion <- ifelse(tests$gradient, "reported apart", ifelse(tests$truth_label == "contiguous", ">= 0.90", "<= 0.05"))
  tests$met <- ifelse(tests$gradient, NA, ifelse(tests$truth_label == "contiguous", tests$warning_rate >= 0.90, tests$warning_rate <= 0.05))
  tests$met_old <- ifelse(tests$gradient, NA, ifelse(tests$truth_label == "contiguous", tests$warning_rate_old >= 0.90,
                                                     tests$warning_rate_old <= 0.05))
  tests <- tests[order(match(tests$data_family, families), tests$gradient, tests$truth_label), ]
  scen <- do.call(rbind, lapply(split(q, list(q$base_id, q$truth_label), drop = TRUE), function(d)
    cbind(d[1, c("base_id", "label", "data_family", "truth_label", "threshold")], rate(d))))
  scen <- scen[order(scen$base_id, scen$truth_label), ]
  write.csv(tests, file.path("results", "stage7_q3v_tests.csv"), row.names = FALSE)
  write.csv(scen, file.path("results", "stage7_q3v_scenarios.csv"), row.names = FALSE)
  cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
  print(tests, digits = 3, row.names = FALSE)
}

step <- commandArgs(trailingOnly = TRUE)
if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  if (!length(step) || !step[1] %in% c("calibrate", "validate")) stop("Give the step: calibrate or validate.")
  if (step[1] == "calibrate") calibrate() else validate()
}
