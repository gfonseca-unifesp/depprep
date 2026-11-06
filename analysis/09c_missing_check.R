# Stage 9 — descriptive check of the correction with missing values (decision D26). Not a
# pre-registered test.
#
# Question. With missing values each variable is regressed on its basis over its observed
# rows, the bases are selected with the gaps at the variable's mean, and R is computed
# pairwise (correct(), D26). Two variables of a group keep the same basis, but when their
# gaps differ their projections differ: the core rule then holds only approximately. How
# much do the bases, R and the classes change with gaps of the sizes and shapes of the Foz
# cube (analysis/09b_fza_gap_census.R)?
#
# Design. In each replicate the same data are corrected complete (the reference) and with
# gaps, with the defaults of correct() and the grouping of the complete data:
# * space-time, 100 sites x 93 times (the length of the Foz series), the groups of id 9
#   (analysis/d4_scenarios.R), null and structured truths;
# * spatial, 400 sites, the groups of id 2, structured truth;
# 10 replicates each. Gaps: a share f = 0.1, 0.2 or 0.3 of the values of a variable, of
# two kinds:
# * "random": every value of every variable missing with probability f, independently;
# * "patches" (like clouds): half of the variables (drawn) lose the same block, a region
#   of the sites nearest to a drawn centre, in space-time only during a season of 4
#   consecutive calendar months of a 12-month year (the region then holds 3f of the
#   sites, so that the block is a share f of the values); the other half stay complete.
# Measures, per replicate, kind and share: the missing share; on the own basis of each
# dependent variable (the data matrix X), k/n, capture and absorption of the true
# components (.capture_absorption(), over all rows; the capture over the variables with a
# large-scale component in the truth) and the distance to the basis of the
# complete data (basis_distance()); the error of R against the correlation of the true
# signal and its distance to the R of the complete data (RMS over the pairs); the share
# of the regime types (classify_regime() on diagnose_variables(nsim = 199)) equal to those
# of the complete data; and the k-means ARI on the rows of X where every variable has a
# value (NA with fewer than 30 such rows), with the number of those rows.
# Seeds: 16400 + 1000 * combination + replicate for the data; for the gaps, the same plus
# 500 + 100 f (plus 50 for the patches), recorded as gaps_seed.
#
# Outputs (not tracked): results/missing_check.csv, results/missing_check_summary.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 16400L
n_rep <- 10L
shares <- c(0.1, 0.2, 0.3)
kinds <- c("random", "patches")
bases_all <- rbind(scenarios, extra_scenarios, d14_scenarios)
grid <- rbind(
  data.frame(base_id = 9, n = 100, n_times = 93, truth = c("null", "structured"), label = "ST 100x93"),
  data.frame(base_id = 2, n = 400, n_times = NA, truth = "structured", label = "S n=400"))
grid$combo <- seq_len(nrow(grid))

# The gaps of one kind and share, as a logical matrix (TRUE = missing). Uses the caller's
# random stream.
make_gaps <- function(s, kind, f) {
  n <- nrow(s$X); p <- ncol(s$X)
  if (kind == "random") return(matrix(stats::runif(n * p) < f, n, p))
  coords <- as.matrix(s$coords)
  st <- !is.null(s$index)
  site_share <- if (st) min(1, 3 * f) else f
  centre <- coords[sample.int(nrow(coords), 1), ]
  d <- sqrt(colSums((t(coords) - centre)^2))
  region <- rank(d, ties.method = "first") <= round(site_share * nrow(coords))
  rows <- if (st) {
    start <- sample.int(12, 1)
    season <- ((start - 1 + 0:3) %% 12) + 1
    region[s$index$site] & ((s$index$time - 1) %% 12 + 1) %in% season
  } else region
  M <- matrix(FALSE, n, p)
  M[rows, sample.int(p, p %/% 2)] <- TRUE
  M
}

one_check <- function(k, r) {
  g <- grid[grid$combo == k, ]
  sc <- bases_all[bases_all$id == g$base_id, setdiff(names(bases_all), "between")]
  sc$n <- g$n
  if (!is.na(g$n_times)) sc$n_times <- g$n_times
  seed <- base_seed + 1000L * k + r
  s <- simulate(sc, seed, g$truth)
  pos <- if (is.null(s$index)) list(coords = s$coords) else list(coords = s$coords, index = s$index)
  gr <- group_regimes(diagnose_variables(s))
  latent <- if (!is.null(s$X_latent)) s$X_latent else s$X
  share <- truth_table(s$groups)$share
  signal <- latent - sweep(s$components$large, 2, sqrt(share), "*")
  R_sig <- suppressWarnings(stats::cor(signal)); ut <- upper.tri(R_sig)
  types0 <- classify_regime(do.call(diagnose_variables, c(list(s$X, nsim = 199), pos)))
  ref <- correct(s, gr, seed = seed)
  cases <- c(list(list(kind = "complete", f = 0, X = s$X, gaps_seed = NA)),
             unlist(lapply(kinds, function(kd) lapply(shares, function(f) {
               gs <- seed + 500L + round(100 * f) + if (kd == "patches") 50L else 0L
               set.seed(gs)
               X <- s$X; X[make_gaps(s, kd, f)] <- NA
               list(kind = kd, f = f, X = X, gaps_seed = gs)
             })), recursive = FALSE))
  do.call(rbind, lapply(cases, function(cs) {
    t0 <- Sys.time()
    cr <- if (cs$kind == "complete") ref else do.call(correct, c(list(cs$X, gr, seed = seed), pos))
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    vb <- cr$variable_bases
    per_var <- do.call(rbind, lapply(names(vb), function(v) {
      j <- match(v, colnames(s$X))
      ca <- .capture_absorption(vb[[v]], s$components$large[, j], s$components$residual[, j])
      has_large <- sum(s$components$large[, j]^2) > 0
      data.frame(kn = ncol(vb[[v]]) / nrow(s$X), capture = if (has_large) ca[["capture"]] else NA_real_, absorption = ca[["absorption"]],
                 distance_to_complete = basis_distance(vb[[v]], ref$variable_bases[[v]]))
    }))
    types <- if (cs$kind == "complete") types0 else classify_regime(do.call(diagnose_variables, c(list(cs$X, nsim = 199), pos)))
    whole <- stats::complete.cases(cr$X)
    km <- if (is.null(s$classes) || sum(whole) < 30) NA_real_ else {
      d <- run_downstream(cr$X[whole, , drop = FALSE], "kmeans",
                          truth = list(classes = s$classes[whole], signal = signal[whole, , drop = FALSE]), seed = seed)
      d$value[d$metric == "ari_true_k"]
    }
    data.frame(combo = k, label = g$label, truth = g$truth, replicate = r, seed = seed, gaps_seed = cs$gaps_seed,
               kind = cs$kind, share = cs$f, missing = mean(is.na(cs$X)), time_correct = time,
               kn = mean(per_var$kn), capture = mean(per_var$capture, na.rm = TRUE), absorption = mean(per_var$absorption),
               distance_to_complete = mean(per_var$distance_to_complete),
               rms_signal = sqrt(mean((cr$R - R_sig)[ut]^2)), rms_to_complete = sqrt(mean((cr$R - ref$R)[ut]^2)),
               regimes_kept = mean(types == types0), complete_rows = sum(whole), kmeans_ari = km)
  }))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(16, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "bases_all", "grid", "shares", "kinds", "make_gaps", "one_check"))
  jobs <- expand.grid(combo = grid$combo, r = seq_len(n_rep))
  jobs <- jobs[order(jobs$combo == 3), ]                                   # the slow space-time jobs first
  res <- clusterApplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) tryCatch(one_check(jobs$combo[i], jobs$r[i]),
    error = function(e) data.frame(combo = jobs$combo[i], replicate = jobs$r[i], error = conditionMessage(e))), jobs = jobs)
  stopCluster(cl)
  ok <- Filter(function(d) !"error" %in% names(d), res)
  bad <- Filter(function(d) "error" %in% names(d), res)
  res <- do.call(rbind, ok)
  dir.create("results", showWarnings = FALSE)
  utils::write.csv(res, file.path("results", "missing_check.csv"), row.names = FALSE)
  by <- res[c("label", "truth", "kind", "share")]
  measures <- c("missing", "time_correct", "kn", "capture", "absorption", "distance_to_complete", "rms_signal",
                "rms_to_complete", "regimes_kept", "complete_rows", "kmeans_ari")
  summ <- stats::aggregate(res[measures], by, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  summ$ari_replicates <- stats::aggregate(list(n = !is.na(res$kmeans_ari)), by, sum)$n   # replicates with an ARI
  summ$replicates <- stats::aggregate(list(n = res$replicate), by, length)$n
  summ <- summ[order(summ$label, summ$truth, match(summ$kind, c("complete", kinds)), summ$share), ]
  utils::write.csv(summ, file.path("results", "missing_check_summary.csv"), row.names = FALSE)
  print(summ, digits = 3, row.names = FALSE)
  cat("errors:", length(bad), "\n")
  if (length(bad)) print(do.call(rbind, bad))
}
