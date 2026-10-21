# Stage 3 — Diagnosis of every variable in the D4 grid, 200 replicates per scenario.
#
# For each scenario (family x sample size x heterogeneity, with the values fixed by
# decision D4), every replicate is simulated, every variable is diagnosed with
# diagnose_variables() and classified with classify_regime(). A subsample of
# replicates also gets the parametric bootstrap (one variable per group), to measure
# the coverage of its 90% intervals.
#
# Seeds: base_seed + 1000 * scenario + replicate (protocol section 9).
# Runs in parallel and can be resumed: one file per scenario in results/diagnostics/,
# scenarios already saved are skipped. Progress in results/diagnostics/progress.log.
#
# Outputs (not tracked by git):
#   results/diagnostics/scenario_NN.rds   signatures, truth and seed of every variable
#   results/diagnostics/bootstrap_NN.rds  bootstrap intervals (subsample)
#   results/diagnostics_variables.rds     all scenarios together
#   results/diagnostics_summary.csv       accuracy and scale error by scenario and regime

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)

base_seed <- 3300L
n_rep <- 200L
n_boot_rep <- 20L
B <- 99L
out_dir <- file.path("results", "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

# --- Scenarios (decision D4) -------------------------------------------------------------
spatial_groups <- function(short_range) sim_groups(
  sim_group("long", 8, "exponential", range = 0.5),
  sim_group("short", 8, "exponential", range = short_range),
  sim_group("aniso", 7, "exponential", range = 0.25, anisotropy = 4, angle = 30),
  sim_group("none", 7, "none"))
temporal_groups <- function(weak_phi) sim_groups(
  sim_group("strong", 6, "ar1", phi = 0.9),
  sim_group("weak", 6, "ar1", phi = weak_phi),
  sim_group("seasonal", 6, "seasonal", period = 12),
  sim_group("trend", 6, "trend_shift"),
  sim_group("none", 6, "none"))
st_groups <- function() sim_groups(
  sim_group("sep_long", 6, "separable_ar1", range = 0.5, phi = 0.9),
  sim_group("sep_short", 6, "separable_seasonal", range = 0.25, period = 12),   # 0.1 was not identifiable (D10)
  sim_group("moving", 6, "moving", range = 0.25, velocity = c(0.02, 0)),
  sim_group("common", 6, "temporal_common", phi = 0.9),
  sim_group("none", 6, "none"))

scenarios <- rbind(
  data.frame(family = "spatial", n = c(100, 400, 100, 400), n_times = NA,
             heterogeneity = c("high", "high", "low", "low")),
  data.frame(family = "temporal", n = NA, n_times = c(120, 480, 120, 480),
             heterogeneity = c("high", "high", "low", "low")),
  data.frame(family = "spatiotemporal", n = c(50, 100), n_times = c(60, 120), heterogeneity = "D4"))
scenarios$id <- seq_len(nrow(scenarios))
scenarios$label <- with(scenarios, ifelse(family == "spatial", sprintf("S n=%d %s", n, heterogeneity),
                    ifelse(family == "temporal", sprintf("T T=%d %s", n_times, heterogeneity),
                           sprintf("ST %dx%d", n, n_times))))
write.csv(scenarios, file.path(out_dir, "scenarios.csv"), row.names = FALSE)

# truth of each variable: regime type (as classify_regime() names it) and scale
truth_table <- function(groups) {
  do.call(rbind, lapply(groups, function(g) {
    src <- if (is.null(g$shared_with)) g else groups[[g$shared_with]]
    p <- src$params
    type <- switch(src$regime,
      exponential = if (!is.null(p$anisotropy) && p$anisotropy > 1) "anisotropic" else "isotropic",
      ar1 = "ar", trend_shift = "trend", seasonal = "seasonal", src$regime)
    data.frame(group = g$name, true_type = type,
               true_range = if (is.null(p$range)) NA_real_ else p$range,
               true_phi = if (is.null(p$phi)) NA_real_ else p$phi,
               true_period = if (is.null(p$period)) NA_real_ else p$period,
               share = g$large_scale_share, stringsAsFactors = FALSE)[rep(1, g$n_vars), ]
  }))
}

simulate <- function(sc, seed) {
  switch(sc$family,
    spatial = sim_spatial(sc$n, spatial_groups(if (sc$heterogeneity == "high") 0.1 else 0.25), seed = seed),
    temporal = sim_temporal(sc$n_times, temporal_groups(if (sc$heterogeneity == "high") 0.3 else 0.6), seed = seed),
    spatiotemporal = sim_spatiotemporal(sc$n, sc$n_times, st_groups(), seed = seed))
}

one_replicate <- function(sc, r) {
  seed <- base_seed + 1000L * sc$id + r
  s <- simulate(sc, seed)
  d <- diagnose_variables(s)
  d$type <- classify_regime(d)
  cbind(scenario = sc$id, replicate = r, seed = seed, truth_table(s$groups), d, row.names = NULL)
}

one_bootstrap <- function(sc, r) {
  seed <- base_seed + 1000L * sc$id + r
  s <- simulate(sc, seed)
  first <- !duplicated(s$group)                                   # one variable per group
  d <- if (sc$family == "spatial") diagnose_variables(s$X[, first, drop = FALSE], coords = s$coords, B = B, seed = seed)
       else diagnose_variables(s$X[, first, drop = FALSE], time = s$time, B = B, seed = seed)
  cbind(scenario = sc$id, replicate = r, seed = seed, truth_table(s$groups)[first, ], d, row.names = NULL)
}

# --- Run -------------------------------------------------------------------------------
# Workers start their own R session and activate renv; its sandbox step can take more
# than a minute when many workers start together, so it is switched off for them and
# the connection timeout is raised.
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE)); NULL }, getwd())
clusterExport(cl, c("base_seed", "B", "spatial_groups", "temporal_groups", "st_groups", "truth_table",
                    "simulate", "one_replicate", "one_bootstrap"))
cat(format(Sys.time()), "start\n", file = log_file, append = TRUE)
for (k in scenarios$id) {
  sc <- scenarios[k, ]
  f <- file.path(out_dir, sprintf("scenario_%02d.rds", k))
  if (!file.exists(f)) {
    t0 <- Sys.time()
    res <- parLapply(cl, seq_len(n_rep), function(r, sc) one_replicate(sc, r), sc = sc)
    saveRDS(do.call(rbind, res), f)
    cat(format(Sys.time()), sc$label, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n",
        file = log_file, append = TRUE)
  }
  fb <- file.path(out_dir, sprintf("bootstrap_%02d.rds", k))
  if (sc$family != "spatiotemporal" && !file.exists(fb)) {
    res <- parLapply(cl, seq_len(n_boot_rep), function(r, sc) one_bootstrap(sc, r), sc = sc)
    saveRDS(do.call(rbind, res), fb)
    cat(format(Sys.time()), sc$label, "bootstrap done\n", file = log_file, append = TRUE)
  }
}
stopCluster(cl)

# --- Summary -----------------------------------------------------------------------------
all <- do.call(rbind, lapply(scenarios$id, function(k) readRDS(file.path(out_dir, sprintf("scenario_%02d.rds", k)))))
all$label <- scenarios$label[all$scenario]
all$type <- unname(classify_regime(all))   # classified again, so new thresholds need no new run
all$correct <- all$type == all$true_type
saveRDS(all, file.path("results", "diagnostics_variables.rds"))

summary_tab <- do.call(rbind, lapply(split(all, list(all$scenario, all$true_type), drop = TRUE), function(d) {
  rr <- d$range / d$true_range
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], true_type = d$true_type[1],
             n_variables = nrow(d), accuracy = mean(d$correct),
             range_ratio_median = stats::median(rr, na.rm = TRUE),
             range_within_25 = mean(abs(rr - 1) <= 0.25, na.rm = TRUE),
             range_failed = mean(is.na(d$range) & !is.na(d$true_range)),
             phi_median = stats::median(d$phi, na.rm = TRUE),
             phi_expected = stats::median(d$share * d$true_phi, na.rm = TRUE))
}))
summary_tab <- summary_tab[order(summary_tab$scenario, summary_tab$true_type), ]
write.csv(summary_tab, file.path("results", "diagnostics_summary.csv"), row.names = FALSE)

boot <- do.call(rbind, lapply(scenarios$id[scenarios$family != "spatiotemporal"], function(k)
  readRDS(file.path(out_dir, sprintf("bootstrap_%02d.rds", k)))))
boot$label <- scenarios$label[boot$scenario]
boot$range_cover <- with(boot, ifelse(is.na(true_range) | is.na(range_lo), NA, range_lo <= true_range & true_range <= range_hi))
boot$phi_target <- boot$share * boot$true_phi                       # lag-1 autocorrelation of the observed variable
boot$phi_cover <- with(boot, ifelse(is.na(phi_target) | is.na(phi_lo), NA, phi_lo <= phi_target & phi_target <= phi_hi))
saveRDS(boot, file.path("results", "diagnostics_bootstrap.rds"))
cat(format(Sys.time()), "summary written\n", file = log_file, append = TRUE)
print(summary_tab, digits = 2)
