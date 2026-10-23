# Stage 6 — Calibration of the Q3 warning of verify(): "signal at the scale of the basis".
#
# In the pilot the correction erased classes contiguous in space or time (regions,
# periods) and nothing warned the user (protocol, question Q3). verify() now measures the
# association between variables carried by the projected part of a conditional
# correction (q3$statistic: the mean over the pairs of standardised variables of
# |cov_large|, on the union of the bases).
# This script computes it in the pilot scenarios (analysis/d4_scenarios.R,
# pilot_scenarios): contiguous classes against the null truth, plus the structured truth
# (classes at random positions, which the correction does not remove) as a check that
# the warning does not fire for classes that are not at the scale of the basis.
# The share of the association in the projected part (q3$share) is recorded too; it did
# not separate the truths in a first look (it is high in the null truth, where no
# association is left in the residuals), so it is not the warning's statistic.
# Approaches: the default of decision D7 (global basis in space and time; in space-time,
# where the default is whitening and has no projected part, one basis per group) and one
# basis per group.
#
# Threshold: the midpoint between the largest null statistic and the smallest contiguous
# statistic when they do not overlap; otherwise the value that maximises the true
# positive rate minus the false positive rate (Youden), contiguous against null.
#
# Seeds: 7700 + 1000 * scenario + replicate (scenario = id in pilot_scenarios),
# replicates 1-5; the same number seeds the bases. Recorded in the outputs.
# Outputs (not tracked by git): results/q3_warning_replicates.csv,
#   results/q3_warning_summary.csv, results/q3_warning_threshold.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 7700L
n_rep <- 5L

one_q3 <- function(sc, r) {
  seed <- base_seed + 1000L * sc$id + r
  s <- simulate(sc, seed, sc$truth)
  gr <- group_regimes(diagnose_variables(s))
  fits <- list(default = if (s$family == "spatiotemporal") NULL else correct(s, gr, seed = seed, data_basis = "union"),
               groups = correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union"))
  if (is.null(fits$default)) fits$default <- fits$groups
  do.call(rbind, lapply(names(fits), function(a) {
    v <- verify(fits[[a]], s, nperm = 0)
    U <- .union_of_bases(fits[[a]], nrow(s$X))
    data.frame(scenario = sc$id, label = sc$label, family = s$family, truth = sc$truth, replicate = r, seed = seed,
               approach = a, scope = fits[[a]]$scope, kn = ncol(U) / nrow(s$X), statistic = v$q3$statistic,
               share = v$q3$share,
               median_variance_share = stats::median(v$variables$large_share[fits[[a]]$group != "none"]),
               residual_structure = mean(v$variables$residual_structure))
  }))
}

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(16, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
clusterExport(cl, c("one_q3", "base_seed"))
jobs <- expand.grid(i = seq_len(nrow(pilot_scenarios)), r = seq_len(n_rep))
jobs <- jobs[order(pilot_scenarios$family[jobs$i] != "spatiotemporal"), ]       # the slow ones first
res <- parLapply(cl, seq_len(nrow(jobs)), function(j, jobs, todo) one_q3(todo[jobs$i[j], ], jobs$r[j]),
                 jobs = jobs, todo = pilot_scenarios)
stopCluster(cl)
res <- do.call(rbind, res)
res <- res[order(res$scenario, res$approach, res$replicate), ]
dir.create("results", showWarnings = FALSE)
write.csv(res, file.path("results", "q3_warning_replicates.csv"), row.names = FALSE)

# threshold: contiguous against null, all families and both approaches together
calibrate <- function(null, pos) {
  if (max(null) < min(pos)) return(c(threshold = (max(null) + min(pos)) / 2, rule = 1))
  cand <- sort(unique(c(null, pos)))
  j <- vapply(cand, function(t) mean(pos > t) - mean(null > t), 0)
  c(threshold = cand[which.max(j)], rule = 2)
}
cal <- calibrate(res$statistic[res$truth == "null"], res$statistic[res$truth == "contiguous"])
threshold <- round(unname(cal["threshold"]), 2)
res$warning <- res$statistic > threshold
summary_tab <- do.call(rbind, lapply(split(res, list(res$label, res$approach), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], truth = d$truth[1],
             approach = d$approach[1], replicates = nrow(d), kn = mean(d$kn), statistic_mean = mean(d$statistic),
             statistic_min = min(d$statistic), statistic_max = max(d$statistic), warning_rate = mean(d$warning))))
summary_tab <- summary_tab[order(summary_tab$scenario, summary_tab$approach), ]
write.csv(summary_tab, file.path("results", "q3_warning_summary.csv"), row.names = FALSE)
write.csv(data.frame(threshold = threshold, rule = if (cal["rule"] == 1) "midpoint, no overlap" else "Youden",
                     null_max = max(res$statistic[res$truth == "null"]),
                     contiguous_min = min(res$statistic[res$truth == "contiguous"]),
                     false_positive_rate = mean(res$warning[res$truth == "null"]),
                     true_positive_rate = mean(res$warning[res$truth == "contiguous"]),
                     structured_rate = mean(res$warning[res$truth == "structured"]),
                     base_seed = base_seed, replicates = n_rep),
          file.path("results", "q3_warning_threshold.csv"), row.names = FALSE)
print(summary_tab, digits = 3, row.names = FALSE)
cat(sprintf("Q3 threshold: %.2f\n", threshold))
