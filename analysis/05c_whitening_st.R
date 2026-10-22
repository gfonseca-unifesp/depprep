# Stage 5 — Why whitening helps in space-time (support for decision D7).
#
# For the spatio-temporal scenarios (ids 9, 10, 15, 16; the data of analysis/05b), the
# within-group error of the residual correlations under the true grouping, with the
# whitening split into its parts: spatial whitening only, AR filtering of each site
# only, both (the whitening mode), the projection on the group basis (conditional mode),
# and the projection followed by the AR filter. Also the share of the large-scale
# variance at slow frequencies. Error = mean absolute difference from the true (sample)
# residual correlations. 10 replicates per scenario.
#
# Seeds: the data seeds of analysis/05b (3300 + 1000 * scenario + replicate); the basis
# selection uses the replicate number as seed.
#
# Outputs (not tracked by git): results/whitening_st.rds, results/whitening_st_summary.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

one <- function(r, id) {
  all_sc <- rbind(scenarios, extra_scenarios)
  sc <- all_sc[all_sc$id == id, ]
  s <- simulate(sc, 3300L + 1000L * id + r)
  truth <- truth_table(s$groups)
  n_s <- nrow(s$coords); tt <- s$index$time; n_t <- max(tt)
  to_cube <- function(v) { m <- matrix(0, n_s, n_t); m[cbind(s$index$site, tt)] <- v; m }
  ar_filter <- function(M, phi) { p <- length(phi); out <- M[-(1:p), , drop = FALSE]
    for (i in seq_len(p)) out <- out - phi[i] * M[(p + 1 - i):(nrow(M) - i), , drop = FALSE]; out }
  res <- list()
  for (g in setdiff(unique(truth$group), "none")) {
    v <- which(truth$group == g); type <- truth$true_type[v[1]]
    Z <- scale(s$X[, v])
    R_res <- stats::cor(s$components$residual[, v])
    # group AR coefficient from the site-mean series (as in the whitening mode)
    tm <- sapply(seq_along(v), function(j) colMeans(to_cube(Z[, j])))
    phi <- mean(sapply(seq_along(v), function(j) stats::ar(tm[, j] - mean(tm[, j]), aic = FALSE, order.max = 1)$ar))
    # spatial whitening matrix from the site means
    sm <- sapply(seq_along(v), function(j) rowMeans(to_cube(Z[, j])))
    fits <- lapply(seq_along(v), function(j) .fit_variogram(s$coords, sm[, j]))
    ok <- which(vapply(fits, function(f) !is.null(f$fit), TRUE))
    W <- if (length(ok)) { S <- gstat::variogramLine(fits[[ok[1]]]$fit, dist_vector = as.matrix(dist(s$coords)), covariance = TRUE)
      solve(t(chol(S + diag(1e-8 * mean(diag(S)), n_s)))) } else diag(n_s)
    variants <- list(
      none = function(z) z,
      spatial_only = function(z) as.vector(W %*% to_cube(z)),
      temporal_only = function(z) as.vector(t(ar_filter(t(to_cube(z)), phi))),
      both = function(z) as.vector(t(ar_filter(t(W %*% to_cube(z)), phi))))
    b <- build_group_basis(s$X[, v], coords = s$coords, index = s$index, type = type,
                           period = truth$true_period[v[1]], seed = r)
    E <- .resid_on(s$X[, v], b)
    err <- function(M) { R <- stats::cor(M); mean(abs((R - R_res)[upper.tri(R)])) }
    out <- vapply(variants, function(f) err(sapply(seq_along(v), function(j) f(Z[, j] - mean(Z[, j])))), 0)
    out["conditional"] <- err(E)
    out["conditional_then_ar"] <- err(sapply(seq_along(v), function(j) as.vector(t(ar_filter(t(to_cube(E[, j])), phi)))))
    # where the shared structure lives: share of the large-scale variance at slow frequencies
    Lc <- to_cube(s$components$large[, v[1]])
    spec <- Mod(stats::mvfft(t(Lc - rowMeans(Lc))))^2
    slow <- sum(spec[2:max(2, floor(n_t / 10)), ]) / sum(spec[-1, ])
    res[[g]] <- data.frame(scenario = sc$label, replicate = r, group = g, type = type, phi = phi,
                           slow_share = slow, t(out))
  }
  do.call(rbind, res)
}

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
invisible(clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd()))
clusterExport(cl, "one")
tasks <- expand.grid(r = 1:10, id = c(9, 10, 15, 16))
out <- do.call(rbind, parLapply(cl, split(tasks, seq_len(nrow(tasks))), function(t) one(t$r, t$id)))
stopCluster(cl)
saveRDS(out, file.path("results", "whitening_st.rds"))
summ <- aggregate(cbind(none, conditional, temporal_only, spatial_only, both, conditional_then_ar, phi, slow_share) ~ scenario + type,
                  data = out, FUN = mean)
write.csv(summ, file.path("results", "whitening_st_summary.csv"), row.names = FALSE)
print(summ, digits = 2, row.names = FALSE)
