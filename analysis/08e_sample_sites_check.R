# Stage 8 — descriptive check of the selection of basis terms on a sample of sites
# (decision D25). Not a pre-registered test.
#
# Question. correct() spends most of its time selecting the basis of each variable
# (78% at 100 sites x 93 times, 30 variables). With sample_sites < 1 the terms are
# selected on a seeded share of the sites (all times) and then taken on all rows. How
# much time does that save, and how much does it change the bases and the outputs?
#
# Design. In each replicate the same data are corrected with sample_sites = 1 (all sites,
# the reference), 0.5 and 0.25, with the defaults of correct():
# * space-time, 200 sites x 93 times, the groups of id 9 (analysis/d4_scenarios.R), null
#   and structured truths, 10 replicates each;
# * spatial, 628 sites (the map of the Foz basin), the groups of id 2, structured truth,
#   10 replicates, sample_sites = 1 and 0.5: few sites per candidate, where the tests are
#   expected to lose power.
# Measures, per replicate and share: the time of correct(); per dependent variable, on
# its own basis (the data matrix X), k/n, capture and absorption of the true components
# (.capture_absorption()) and the distance to the basis selected on all sites
# (basis_distance()); on the outputs, the error of R against the correlation of the true
# signal and the k-means ARI on X.
# Seeds: 15400 + 1000 * combination + replicate.
#
# Outputs (not tracked): results/sample_sites_check.csv, results/sample_sites_check_summary.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 15400L
n_rep <- 10L
bases_all <- rbind(scenarios, extra_scenarios, d14_scenarios)
grid <- rbind(
  data.frame(base_id = 9, n = 200, n_times = 93, truth = c("null", "structured"), label = "ST 200x93", shares = "1 0.5 0.25"),
  data.frame(base_id = 2, n = 628, n_times = NA, truth = "structured", label = "S n=628", shares = "1 0.5"))
grid$combo <- seq_len(nrow(grid))

one_check <- function(k, r) {
  g <- grid[grid$combo == k, ]
  sc <- bases_all[bases_all$id == g$base_id, setdiff(names(bases_all), "between")]
  sc$n <- g$n
  if (!is.na(g$n_times)) sc$n_times <- g$n_times
  seed <- base_seed + 1000L * k + r
  s <- simulate(sc, seed, g$truth)
  gr <- group_regimes(diagnose_variables(s))
  latent <- if (!is.null(s$X_latent)) s$X_latent else s$X
  share <- truth_table(s$groups)$share
  signal <- latent - sweep(s$components$large, 2, sqrt(share), "*")
  R_sig <- suppressWarnings(stats::cor(signal)); ut <- upper.tri(R_sig)
  truth <- list(classes = s$classes, signal = signal)
  fits <- lapply(as.numeric(strsplit(g$shares, " ")[[1]]), function(f) {
    t0 <- Sys.time()
    cr <- correct(s, gr, seed = seed, sample_sites = f)
    list(f = f, cr = cr, time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  })
  ref <- fits[[1]]$cr$variable_bases
  do.call(rbind, lapply(fits, function(ft) {
    vb <- ft$cr$variable_bases
    per_var <- do.call(rbind, lapply(names(vb), function(v) {
      j <- match(v, colnames(s$X))
      ca <- .capture_absorption(vb[[v]], s$components$large[, j], s$components$residual[, j])
      data.frame(variable = v, kn = ncol(vb[[v]]) / nrow(s$X), capture = ca[["capture"]], absorption = ca[["absorption"]],
                 distance_to_full = basis_distance(vb[[v]], ref[[v]]))
    }))
    km <- if (is.null(s$classes)) NA_real_ else {
      d <- run_downstream(ft$cr$X, "kmeans", truth = truth, seed = seed)
      d$value[d$metric == "ari_true_k"]
    }
    data.frame(combo = k, label = g$label, truth = g$truth, replicate = r, seed = seed, sample_sites = ft$f, time_correct = ft$time,
               dependent = nrow(per_var), kn = mean(per_var$kn), capture = mean(per_var$capture), absorption = mean(per_var$absorption),
               distance_to_full = mean(per_var$distance_to_full), rms_signal = sqrt(mean((ft$cr$R - R_sig)[ut]^2)), kmeans_ari = km)
  }))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(16, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "bases_all", "grid", "one_check"))
  jobs <- expand.grid(combo = grid$combo, r = seq_len(n_rep))
  jobs <- jobs[order(jobs$combo != 1 & jobs$combo != 2), ]                  # the slow space-time jobs first
  res <- clusterApplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) tryCatch(one_check(jobs$combo[i], jobs$r[i]),
    error = function(e) data.frame(combo = jobs$combo[i], replicate = jobs$r[i], error = conditionMessage(e))), jobs = jobs)
  stopCluster(cl)
  ok <- Filter(function(d) !"error" %in% names(d), res)
  bad <- Filter(function(d) "error" %in% names(d), res)
  res <- do.call(rbind, ok)
  dir.create("results", showWarnings = FALSE)
  utils::write.csv(res, file.path("results", "sample_sites_check.csv"), row.names = FALSE)
  summ <- stats::aggregate(cbind(time_correct, kn, capture, absorption, distance_to_full, rms_signal, kmeans_ari) ~ label + truth + sample_sites,
                           data = transform(res, kmeans_ari = ifelse(is.na(kmeans_ari), -1, kmeans_ari)), FUN = mean)
  summ$kmeans_ari[summ$kmeans_ari < 0] <- NA
  summ$replicates <- stats::aggregate(replicate ~ label + truth + sample_sites, data = res, FUN = length)$replicate
  utils::write.csv(summ, file.path("results", "sample_sites_check_summary.csv"), row.names = FALSE)
  print(summ, digits = 3, row.names = FALSE)
  cat("errors:", length(bad), "\n")
  if (length(bad)) print(do.call(rbind, bad))
}
