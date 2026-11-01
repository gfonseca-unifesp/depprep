# Stage 8 — time of correct() on larger space-time grids, on one process and on several
# (decision D25). Not a pre-registered test.
#
# Question. In space-time, correct() spends most of its time selecting the basis of each
# variable. How long does it take as the number of sites grows, and how much do the
# processes of `cores` save? The result must be identical: the bases of the variables are
# built in parallel with the same seeds.
#
# Design. One data set per size: 100, 200 and 300 sites x 93 times (the length of the
# real case), the groups of id 9 (analysis/d4_scenarios.R), structured truth, seed 777.
# Each is corrected with the defaults of correct(), on one process and on
# min(16, cores - 4) processes, and the outputs (R, X and the bases of the variables)
# are compared. Times are wall-clock seconds on this machine; nothing else should run
# at the same time.
#
# Outputs (not tracked): results/correct_time.csv

suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("analysis", "d4_scenarios.R"))
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")

seed <- 777L
n_cores <- max(1L, min(16L, parallel::detectCores() - 4L))
sc <- rbind(scenarios, extra_scenarios, d14_scenarios)
sc <- sc[sc$id == 9, setdiff(names(sc), "between")]
sc$n_times <- 93

res <- list()
for (n in c(100, 200, 300)) {
  g <- sc; g$n <- n
  s <- simulate(g, seed, "structured")
  gr <- group_regimes(diagnose_variables(s))
  t1 <- system.time(c1 <- correct(s, gr, seed = seed))[["elapsed"]]
  tp <- system.time(cp <- correct(s, gr, seed = seed, cores = n_cores))[["elapsed"]]
  k <- vapply(c1$variable_bases, function(x) if (is.null(x)) 0L else ncol(x), 0L)
  res[[length(res) + 1]] <- data.frame(
    n_sites = n, n_times = 93, n_obs = nrow(s$X), variables = ncol(s$X), seed = seed, cores = n_cores,
    seconds_1_core = t1, seconds_cores = tp, speedup = t1 / tp,
    identical = identical(c1$R, cp$R) && identical(c1$X, cp$X) && identical(c1$variable_bases, cp$variable_bases),
    own_k_mean = mean(k), own_k_max = max(k))
  print(res[[length(res)]], row.names = FALSE)
  dir.create("results", showWarnings = FALSE)
  utils::write.csv(do.call(rbind, res), file.path("results", "correct_time.csv"), row.names = FALSE)
}
