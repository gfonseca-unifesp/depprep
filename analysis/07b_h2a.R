# Stage 7 — H2a on new data (protocol-v1, section 2.3): the error of correcting two
# variables with different bases against its theory, rho (k - sum cos^2) / (n - k), and the
# error of a shared basis. The design and the functions are those of
# analysis/05a_h2_mechanism.R (read from that file up to its run section); only the seeds
# change: base 8800, the stage 7 base, so the stage 5 results are not reused.
#
# Outputs (not tracked): results/stage7_h2/h2a.rds, results/stage7_h2a_summary.csv

src <- readLines(file.path("analysis", "05a_h2_mechanism.R"))
eval(parse(text = src[seq_len(grep("^# --- Run", src) - 1)]))   # loads the package and the functions
library(parallel)
base_seed <- 8800L                                       # seed_of() reads it
out_dir <- file.path("results", "stage7_h2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
f <- file.path(out_dir, "h2a.rds")
if (!file.exists(f)) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE)); NULL }, getwd())
  clusterExport(cl, c("base_seed", "seed_of", "field", "sim_residuals", "resid_cor", "capture", "one_h2a"))
  tasks <- merge(expand.grid(j = c(0, 5, 10, 15, 20, 25, 30), rho = c(0, 0.5)), data.frame(r = seq_len(500L)))
  res <- do.call(rbind, parLapply(cl, split(tasks, seq_len(nrow(tasks))), function(t, task_fun) do.call(task_fun, as.list(t)),
                                  task_fun = one_h2a))
  stopCluster(cl)
  saveRDS(res, f)
}
h2a <- readRDS(f)
s <- do.call(rbind, lapply(split(h2a, list(h2a$j, h2a$rho), drop = TRUE), function(x)
  data.frame(j = x$j[1], rho = x$rho[1], n_rep = nrow(x), distance = mean(x$distance),
             error_different_mean = mean(x$error_different), error_shared_mean = mean(x$error_shared),
             theory = -x$rho[1] * x$j[1] / (400 - 40))))
s$within_0.005 <- abs(s$error_different_mean - s$theory) < 0.005
s$shared_below_0.002 <- abs(s$error_shared_mean) < 0.002
write.csv(s[order(s$rho, s$j), ], file.path("results", "stage7_h2a_summary.csv"), row.names = FALSE)
