# Stage 8 — the Q3 thresholds of amendment protocol-v1.2 at the size of the real case
# (decision D23), before stage 9.
#
# Why. The thresholds (spatial 0.16, temporal 0.19, space-time 0.13) were validated at the
# sizes of the stage 7 grid. The Foz do Amazonas case is outside them: monthly series of
# 93 months (2015-04 to 2022-12) on 628 cells of 0.25 degrees (the marine basin;
# results/fza/coverage_summary.csv). In the pilot of analysis/08a_q3_reference_pilot.R
# one of 20 replicates without contiguous classes passed the temporal threshold at T = 93.
#
# Design. The groups and truths of the stage 7 grid (analysis/d4_scenarios.R), at the
# sizes of the real case:
# * temporal, T = 93: the groups of ids 5 (high heterogeneity) and 7 (low), 100 replicates;
# * spatial, n = 628 (a map of the basin, for example a time mean): the groups of ids 2
#   (high) and 4 (low), 100 replicates;
# * space-time, 628 sites x 93 times: the groups of id 9 (rough fields), 16 replicates.
#   One replicate at this size took 2.4 h on one core (correct 1.9 h, verify 0.5 h), so the
#   space-time part runs on 16 workers, one replicate per worker, about 12 h in all; it is a
#   coarse check (1 warning in 16 is 6%), run after the temporal and spatial parts.
# Five truths each (null, structured, contiguous, crossed, 0.3 between groups). The same
# path as analysis/07d_q3_recalibration.R up to the statistic; the thresholds are those of
# results/q3_thresholds_v12.csv, fixed.
# Criterion (protocol section 6, unchanged): at most 5% of warnings in each truth without
# contiguous classes, at least 90% in the contiguous truth, by family and truth.
# Seeds: 14300 + 1000 * combination + replicate (combinations numbered below).
#
# Outputs (not tracked): results/q3_size_check/combo_NNN.rds, results/q3_size_check.csv,
#   results/q3_size_check_tests.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 14300L
N_REP <- 100L
N_REP_ST <- 16L          # space-time at 628 x 93: one replicate per worker
WORKERS_ST <- 16L        # fewer workers for the large space-time fits (memory)
out_dir <- file.path("results", "q3_size_check")
bases_all <- rbind(scenarios, extra_scenarios, d14_scenarios)
truths <- data.frame(truth = c("null", "structured", "contiguous", "crossed", "structured"),
                     between = c(0, 0, 0, 0, 0.3), truth_label = c("null", "structured", "contiguous", "crossed", "between03"),
                     stringsAsFactors = FALSE)
resize <- function(id, n = NA, n_times = NA, label) {
  b <- bases_all[bases_all$id == id, setdiff(names(bases_all), "between")]
  if (!is.na(n)) b$n <- n
  if (!is.na(n_times)) b$n_times <- n_times
  b$label <- label; b$base_id <- id
  cbind(b[rep(1, nrow(truths)), ], truths, row.names = NULL)
}
grid <- rbind(resize(5, n_times = 93, label = "T T=93 high"), resize(7, n_times = 93, label = "T T=93 low"),
              resize(2, n = 628, label = "S n=628 high"), resize(4, n = 628, label = "S n=628 low"),
              resize(9, n = 628, n_times = 93, label = "ST 628x93"))
grid$data_family <- ifelse(grid$family == "fields", "spatial", grid$family)
grid$combo <- seq_len(nrow(grid))
grid$reps <- ifelse(grid$family == "spatiotemporal", N_REP_ST, N_REP)

one_size <- function(g, r) {
  seed <- base_seed + 1000L * g$combo + r
  invisible(gc(reset = TRUE))
  s <- simulate(g, seed, g$truth, between = g$between)
  gr <- group_regimes(diagnose_variables(s))
  fit <- if (s$family == "spatiotemporal") correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union") else
    correct(s, gr, seed = seed, data_basis = "union")
  v <- verify(fit, s, nperm = 0, seed = seed)
  data.frame(combo = g$combo, replicate = r, seed = seed, statistic = v$q3$statistic,
             kn = if (!length(fit$bases)) NA_real_ else ncol(.union_of_bases(fit, nrow(s$X))) / nrow(s$X),
             mem_mb = sum(gc()[, 6]))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  thr <- utils::read.csv(file.path("results", "q3_thresholds_v12.csv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "progress.log")
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  start_cluster <- function(workers) {
    cl <- makeCluster(workers, setup_timeout = 900)
    clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
      source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
    clusterExport(cl, c("base_seed", "one_size"), envir = globalenv())
    cl
  }
  # summaries of the combinations done so far (written again after each part)
  summarise_size <- function() {
    done <- grid$combo[file.exists(file.path(out_dir, sprintf("combo_%03d.rds", grid$combo)))]
    q <- do.call(rbind, lapply(done, function(k) readRDS(file.path(out_dir, sprintf("combo_%03d.rds", k)))))
    q <- merge(q, grid[, c("combo", "label", "base_id", "data_family", "truth_label")], by = "combo")
    q$threshold <- thr$threshold[match(q$data_family, thr$data_family)]
    q$warning <- q$statistic > q$threshold
    utils::write.csv(q, file.path("results", "q3_size_check.csv"), row.names = FALSE)
    tests <- do.call(rbind, lapply(split(q, list(q$label, q$truth_label), drop = TRUE), function(d)
      data.frame(label = d$label[1], data_family = d$data_family[1], truth_label = d$truth_label[1], threshold = d$threshold[1],
                 replicates = nrow(d), statistic_mean = mean(d$statistic), statistic_q95 = unname(stats::quantile(d$statistic, 0.95)),
                 statistic_q05 = unname(stats::quantile(d$statistic, 0.05)), warning_rate = mean(d$warning), kn_mean = mean(d$kn, na.rm = TRUE),
                 mem_mb_max = max(d$mem_mb))))
    tests$criterion <- ifelse(tests$truth_label == "contiguous", ">= 0.90", "<= 0.05")
    tests$met <- ifelse(tests$truth_label == "contiguous", tests$warning_rate >= 0.90, tests$warning_rate <= 0.05)
    tests <- tests[order(tests$label, tests$truth_label), ]
    utils::write.csv(tests, file.path("results", "q3_size_check_tests.csv"), row.names = FALSE)
    cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
    print(tests, digits = 3, row.names = FALSE)
  }
  cl <- start_cluster(max(1, min(24, detectCores() - 4)))
  cat(format(Sys.time()), "start;", nrow(grid), "combinations\n", file = log_file, append = TRUE)
  # the cheap temporal and spatial parts first, then the space-time part on fewer workers
  for (k in grid$combo[order(grid$family == "spatiotemporal", grid$combo)]) {
    g <- grid[grid$combo == k, ]
    if (g$family == "spatiotemporal" && length(cl) != WORKERS_ST) {
      summarise_size(); stopCluster(cl); cl <- start_cluster(WORKERS_ST)
    }
    f <- file.path(out_dir, sprintf("combo_%03d.rds", k))
    if (file.exists(f)) next
    t0 <- Sys.time()
    saveRDS(do.call(rbind, parLapply(cl, seq_len(g$reps), function(r, g) one_size(g, r), g = g)), f)
    cat(format(Sys.time()), sprintf("combo %d (%s %s) done in %.1f min\n", k, g$label, g$truth_label,
                                    as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_file, append = TRUE)
  }
  stopCluster(cl)
  summarise_size()
}
