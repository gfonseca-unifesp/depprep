# Stage 8 — pilot of a Q3 statistic normalised by an independence reference (decision D22,
# second part; the choice among the normalisations is decision D23).
#
# Why. Without contiguous classes, the Q3 statistic of verify() grows with the size of the
# basis and is larger in short series (analysis/07e_q3_cause.R): the projected parts of
# independent autocorrelated variables are correlated by chance. Amendment protocol-v1.2
# answered with one threshold per data family, validated only at the sizes of the stage 7
# grid. The real case (93 months) is outside them.
#
# Reference. For each replicate, B = 99 sets of spectral surrogates (each variable keeps its
# spectrum on a complete eigenbasis, the variables become independent; Wagner and Dray
# 2015), and the statistic recomputed on the same projector (verify.R: .q3_reference()).
# Three normalisations of the observed statistic: excess (observed - reference mean), z
# (excess / reference sd) and ratio (observed / reference mean); the raw statistic as the
# comparison. In the truths without classes the variables of a group share large-scale
# fields, so the observed value is above independence there too; what is sought is a
# normalisation whose threshold is the same in every family and size.
#
# Design. The combinations of results/stage7/grid.csv without perturbation, plus two out of
# the grid at the length of the real case: temporal T = 93 (as T=120 high) and space-time
# 50 sites x 93 times (as ST 50x60); five truths each. 5 replicates per combination
# (2 for space-time 100 x 120). Seeds: 10500 + 1000 * combination + replicate; the extra
# combinations are numbered after those of the grid. The same path as
# analysis/07d_q3_recalibration.R up to the statistic.
#
# Outputs (not tracked): results/q3_reference_pilot.csv, results/q3_reference_pilot_summary.csv,
#   results/q3_reference_pilot_thresholds.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 10500L
B <- 99L
grid <- read.csv(file.path("results", "stage7", "grid.csv"), stringsAsFactors = FALSE)
grid <- grid[grid$perturbation == "none", ]
truths <- unique(grid[, c("truth", "between", "truth_label")])
extra <- function(id, n_times, label) {
  b <- grid[grid$base_id == id, ][1, setdiff(names(grid), c("truth", "between", "truth_label", "combo", "n_rep"))]
  b$n_times <- n_times; b$label <- label; b$base_id <- 900 + id
  cbind(b[rep(1, nrow(truths)), ], truths)
}
ext <- rbind(extra(5, 93, "T T=93 high"), extra(9, 93, "ST 50x93"))
ext$combo <- max(grid$combo) + seq_len(nrow(ext))
ext$n_rep <- NA
grid <- rbind(grid[, names(ext)], ext)
grid$data_family <- ifelse(grid$family == "fields", "spatial", grid$family)
grid$reps <- ifelse(grid$family == "spatiotemporal" & !is.na(grid$n) & grid$n == 100, 2L, 5L)

one_ref <- function(g, r) {
  seed <- base_seed + 1000L * g$combo + r
  s <- simulate(g, seed, g$truth, between = g$between)
  gr <- group_regimes(diagnose_variables(s))
  fit <- if (s$family == "spatiotemporal") correct(s, gr, mode = "conditional", scope = "groups", seed = seed, data_basis = "union") else
    correct(s, gr, seed = seed, data_basis = "union")
  v <- verify(fit, s, nperm = 0, seed = seed)
  U <- .union_of_bases(fit, nrow(s$X))
  sb <- .surrogate_bases(s$family, coords = s$coords, time = if (s$family == "temporal") s$time, index = s$index)
  ref <- .q3_reference(s$X, U, sb, B = B, seed = seed)
  data.frame(combo = g$combo, replicate = r, seed = seed, statistic = v$q3$statistic, ref_mean = ref$mean, ref_sd = ref$sd,
             kn = ncol(U) / nrow(s$X))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "B", "one_ref"))
  jobs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) data.frame(i = i, r = seq_len(grid$reps[i]))))
  jobs <- jobs[order(grid$family[jobs$i] != "spatiotemporal"), ]                  # the slow ones first
  res <- parLapply(cl, seq_len(nrow(jobs)), function(j, jobs, grid) tryCatch(one_ref(grid[jobs$i[j], ], jobs$r[j]),
    error = function(e) data.frame(combo = grid$combo[jobs$i[j]], replicate = jobs$r[j], seed = NA, statistic = NA, ref_mean = NA,
                                   ref_sd = NA, kn = NA, error = conditionMessage(e))), jobs = jobs, grid = grid)
  stopCluster(cl)
  res <- do.call(rbind, lapply(res, function(d) { if (is.null(d$error)) d$error <- NA_character_; d }))
  res <- merge(res, grid[, c("combo", "label", "base_id", "data_family", "truth_label")], by = "combo")
  res$excess <- res$statistic - res$ref_mean
  res$z <- res$excess / res$ref_sd
  res$ratio <- res$statistic / res$ref_mean
  dir.create("results", showWarnings = FALSE)
  write.csv(res, file.path("results", "q3_reference_pilot.csv"), row.names = FALSE)

  # per normalisation: the midpoint threshold of amendment v1.2 by family and pooled, and the
  # warning rates of the pooled threshold by family, truth and size (gradient id 19 apart)
  ok <- res[is.finite(res$statistic) & res$base_id != 19, ]
  in_grid <- ok[ok$base_id < 900, ]
  mid <- function(d, m) {
    neg <- d[d$truth_label != "contiguous", ]
    a <- max(vapply(split(neg[[m]], neg$truth_label), function(v) unname(stats::quantile(v, 0.975, type = 1)), 0))
    b <- unname(stats::quantile(d[[m]][d$truth_label == "contiguous"], 0.05, type = 1))
    c(negative_q975 = a, contiguous_q05 = b, threshold = (a + b) / 2, gap = b - a)
  }
  measures <- c("statistic", "excess", "z", "ratio")
  thr <- do.call(rbind, lapply(measures, function(m) rbind(
    data.frame(measure = m, data_family = "pooled", t(mid(in_grid, m))),
    do.call(rbind, lapply(split(in_grid, in_grid$data_family), function(d) data.frame(measure = m, data_family = d$data_family[1], t(mid(d, m))))))))
  write.csv(thr, file.path("results", "q3_reference_pilot_thresholds.csv"), row.names = FALSE)
  summ <- do.call(rbind, lapply(measures, function(m) {
    t0 <- thr$threshold[thr$measure == m & thr$data_family == "pooled"]
    do.call(rbind, lapply(split(ok, list(ok$label, ok$truth_label), drop = TRUE), function(d)
      data.frame(measure = m, label = d$label[1], base_id = d$base_id[1], data_family = d$data_family[1], truth_label = d$truth_label[1],
                 replicates = nrow(d), mean = mean(d[[m]]), warning_rate_pooled = mean(d[[m]] > t0))))
  }))
  write.csv(summ, file.path("results", "q3_reference_pilot_summary.csv"), row.names = FALSE)
  print(thr, digits = 3, row.names = FALSE)
  cat("errors:", sum(!is.na(res$error)), "\n")
}
