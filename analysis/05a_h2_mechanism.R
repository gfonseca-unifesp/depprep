# Stage 5 — The mechanism (H2a, H2b) and the calibration of the second grouping level (D11).
#
# Variables x_j = sqrt(0.6) L_j + sqrt(0.4) e_j on n random sites; the residuals e have a
# known correlation rho. Error = estimated minus true (sample) residual correlation.
#
#   H2a  the between-basis error. L lies exactly in the span of 10 dbMEM and both bases
#        contain it; each basis adds 30 other dbMEM, of which j differ between the two
#        bases (distance between the bases d = j / 40). Theory: bias ~ -rho j / (n - k).
#        500 replicates per j and rho.
#   H2b  the uncaptured shared structure. One shared field L, rough (exponential) or
#        smooth (Matern 2.5), range 0.5; groups of 2, 4, 8 or 16 variables; white or
#        autocorrelated residuals (70% of a range-0.1 field). For the pair (x1, x2):
#        independent bases (one per variable), the shared basis of the group, and the
#        whitening mode. 200 replicates per cell.
#   D11  pairs with the same field or with independent fields of the same range; the
#        distance between their individually selected bases, by range, roughness, n and
#        rho; ROC, AUC and the threshold that best separates them. 100 replicates.
#   didactic  one data set with two groups sharing one field, for heatmaps of R_true,
#        raw R, R corrected with independent bases and with shared bases.
#
# Seeds: 5500 + 1000 * part + replicate (part 1 = H2a, 2 = H2b, 3 = D11, 4 = didactic).
# Runs in parallel and can be resumed (one file per part in results/h2/).
#
# Outputs (not tracked by git): results/h2/{h2a,h2b,d11,didactic}.rds and
#   results/h2a_summary.csv, results/h2b_summary.csv, results/d11_summary.csv,
#   results/d11_threshold.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)

base_seed <- 5500L
out_dir <- file.path("results", "h2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")
seed_of <- function(part, r) base_seed + 1000L * part + r

# --- Shared pieces --------------------------------------------------------------------------
field <- function(coords, range, smoothness, n_sim = 1) depprep:::.grf_chol(coords, range, n_sim = n_sim, smoothness = smoothness)
sim_residuals <- function(coords, m, rho, kind) {
  R <- matrix(rho, m, m); diag(R) <- 1
  W <- matrix(stats::rnorm(nrow(coords) * m), ncol = m) %*% chol(R)
  if (kind == "white") return(W)
  Z <- field(coords, 0.1, 0.5, n_sim = m) %*% chol(R)
  sqrt(0.7) * Z + sqrt(0.3) * W
}
resid_cor <- function(B, x, y) { E <- depprep:::.resid_on(cbind(x, y), B); stats::cor(E[, 1], E[, 2]) }
capture <- function(B, L) depprep:::.capture_absorption(B, L, cbind(L))[["capture"]]

# --- H2a: controlled difference between bases ----------------------------------------------
one_h2a <- function(r, j, rho) {
  set.seed(seed_of(1, r))
  n <- 400
  coords <- matrix(stats::runif(2 * n), n)
  cand <- depprep:::.spatial_candidates(coords)
  L <- as.numeric(cand[, 1:10] %*% stats::rnorm(10)); L <- L / stats::sd(L)
  e <- sim_residuals(coords, 2, rho, "white")
  x <- sqrt(0.6) * L + sqrt(0.4) * e[, 1]; y <- sqrt(0.6) * L + sqrt(0.4) * e[, 2]
  extra <- sample(11:ncol(cand), 30 + j)
  b1 <- cand[, c(1:10, extra[1:30])]
  b2 <- cand[, c(1:10, extra[seq_len(30 - j)], utils::tail(extra, j))]
  tr <- stats::cor(e[, 1], e[, 2])
  B1 <- structure(b1, class = "depprep_basis"); B2 <- structure(b2, class = "depprep_basis")
  E1 <- depprep:::.resid_on(cbind(x), B1); E2 <- depprep:::.resid_on(cbind(y), B2)
  data.frame(replicate = r, j = j, rho = rho, distance = basis_distance(B1, B2),
             error_different = stats::cor(E1[, 1], E2[, 1]) - tr, error_shared = resid_cor(B1, x, y) - tr,
             theory = -rho * j / (n - 40))
}

# --- H2b: capture of the shared structure ----------------------------------------------------
one_h2b <- function(r, smoothness, m, rho, kind) {
  seed <- seed_of(2, r)
  set.seed(seed)
  n <- 400
  coords <- matrix(stats::runif(2 * n), n)
  L <- as.numeric(field(coords, 0.5, smoothness))
  E <- sim_residuals(coords, m, rho, kind)
  X <- sqrt(0.6) * L + sqrt(0.4) * E
  colnames(X) <- paste0("v", seq_len(m))
  tr <- stats::cor(E[, 1], E[, 2])
  b1 <- build_group_basis(X[, 1, drop = FALSE], coords = coords, type = "isotropic", seed = seed)
  b2 <- build_group_basis(X[, 2, drop = FALSE], coords = coords, type = "isotropic", seed = seed + 1e5)
  bs <- build_group_basis(X, coords = coords, type = "isotropic", seed = seed)
  E1 <- depprep:::.resid_on(X[, 1, drop = FALSE], b1); E2 <- depprep:::.resid_on(X[, 2, drop = FALSE], b2)
  grp <- stats::setNames(rep("g1", m), colnames(X))
  cw <- correct(X, grp, coords = coords, bases = list(g1 = bs), mode = "whiten")
  data.frame(replicate = r, smoothness = smoothness, m = m, rho = rho, residual = kind, seed = seed,
             error_none = stats::cor(X[, 1], X[, 2]) - tr,
             error_independent = stats::cor(E1[, 1], E2[, 1]) - tr,
             error_shared = resid_cor(bs, X[, 1], X[, 2]) - tr,
             error_whiten = cw$R[1, 2] - tr,
             capture_independent = (capture(b1, L) + capture(b2, L)) / 2, capture_shared = capture(bs, L),
             k_independent = (ncol(b1) + ncol(b2)) / 2, k_shared = ncol(bs),
             distance_independent = basis_distance(b1, b2))
}

# --- D11: same field or independent fields? ------------------------------------------------
one_d11 <- function(r, same, range, smoothness, n, rho) {
  seed <- seed_of(3, r)
  set.seed(seed)
  coords <- matrix(stats::runif(2 * n), n)
  L <- field(coords, range, smoothness, n_sim = 2)
  if (same) L[, 2] <- L[, 1]
  E <- sim_residuals(coords, 2, rho, "white")
  X <- sqrt(0.6) * L + sqrt(0.4) * E
  b1 <- build_group_basis(X[, 1, drop = FALSE], coords = coords, type = "isotropic", seed = seed)
  b2 <- build_group_basis(X[, 2, drop = FALSE], coords = coords, type = "isotropic", seed = seed + 1e5)
  data.frame(replicate = r, same_field = same, range = range, smoothness = smoothness, n = n, rho = rho, seed = seed,
             distance = basis_distance(b1, b2), k1 = ncol(b1), k2 = ncol(b2))
}

# --- Didactic example ------------------------------------------------------------------------
didactic <- function() {
  seed <- seed_of(4, 1)
  g <- sim_groups(sim_group("A", 4, "matern", range = 0.5), sim_group("B", 4, "matern", shared_with = "A"),
                  sim_group("C", 4, "matern", range = 0.2))
  R_true <- make_R_true(g, "blocks", within = 0.5, between = 0)
  s <- sim_spatial(400, g, R_true = R_true, seed = seed)
  indep <- vapply(seq_len(ncol(s$X)), function(j) {
    b <- build_group_basis(s$X[, j, drop = FALSE], coords = s$coords, type = "isotropic", seed = seed + j)
    depprep:::.resid_on(s$X[, j, drop = FALSE], b)
  }, numeric(nrow(s$X)))
  grp <- stats::setNames(ifelse(s$group == "C", "g2", "g1"), colnames(s$X))   # A and B share their field
  shared <- correct(s, grp, seed = seed)
  list(seed = seed, group = s$group, R_true = R_true, R_raw = stats::cor(s$X), R_independent = stats::cor(indep),
       R_shared = shared$R, R_residual = stats::cor(s$components$residual))
}

# --- Run -------------------------------------------------------------------------------------
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE)); NULL }, getwd())
clusterExport(cl, c("base_seed", "seed_of", "field", "sim_residuals", "resid_cor", "capture", "one_h2a", "one_h2b", "one_d11"))
run_part <- function(name, grid, fun, reps) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  tasks <- merge(grid, data.frame(r = seq_len(reps)))
  res <- do.call(rbind, parLapply(cl, split(tasks, seq_len(nrow(tasks))), function(t, task_fun) do.call(task_fun, as.list(t)), task_fun = fun))
  saveRDS(res, f)
  cat(format(Sys.time()), name, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n",
      file = log_file, append = TRUE)
  res
}
cat(format(Sys.time()), "start\n", file = log_file, append = TRUE)
h2a <- run_part("h2a", expand.grid(j = c(0, 5, 10, 15, 20, 25, 30), rho = c(0, 0.5)), one_h2a, 500L)
h2b <- run_part("h2b", expand.grid(smoothness = c(0.5, 2.5), m = c(2, 4, 8, 16), rho = c(0, 0.5),
                                   kind = c("white", "autocorrelated"), stringsAsFactors = FALSE), one_h2b, 200L)
d11 <- run_part("d11", expand.grid(same = c(TRUE, FALSE), range = c(0.25, 0.5), smoothness = c(0.5, 2.5),
                                   n = c(100, 400), rho = c(0, 0.5)), one_d11, 100L)
stopCluster(cl)
f_did <- file.path(out_dir, "didactic.rds")
if (!file.exists(f_did)) saveRDS(didactic(), f_did)

# --- Summaries -------------------------------------------------------------------------------
agg <- function(d, by, cols) do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(x)
  cbind(x[1, by, drop = FALSE], n_rep = nrow(x),
        as.data.frame(as.list(unlist(lapply(cols, function(cc) c(stats::setNames(mean(x[[cc]]), paste0(cc, "_mean")),
                                                                    stats::setNames(sqrt(mean(x[[cc]]^2)), paste0(cc, "_rmse"))))))))))
h2a_sum <- agg(h2a, c("j", "rho"), c("error_different", "error_shared"))
h2a_sum$distance <- tapply(h2a$distance, list(h2a$j, h2a$rho), mean)[cbind(as.character(h2a_sum$j), as.character(h2a_sum$rho))]
h2a_sum$theory <- -h2a_sum$rho * h2a_sum$j / (400 - 40)
write.csv(h2a_sum[order(h2a_sum$rho, h2a_sum$j), ], file.path("results", "h2a_summary.csv"), row.names = FALSE)

h2b_sum <- agg(h2b, c("smoothness", "m", "rho", "residual"),
               c("error_none", "error_independent", "error_shared", "error_whiten"))
cap <- aggregate(cbind(capture_independent, capture_shared, k_independent, k_shared, distance_independent) ~
                   smoothness + m + rho + residual, data = h2b, FUN = mean)
h2b_sum <- merge(h2b_sum, cap)
write.csv(h2b_sum[order(h2b_sum$smoothness, h2b_sum$residual, h2b_sum$rho, h2b_sum$m), ],
          file.path("results", "h2b_summary.csv"), row.names = FALSE)

auc <- function(a, b) mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="))
d11_sum <- do.call(rbind, lapply(split(d11, d11[c("range", "smoothness", "n", "rho")], drop = TRUE), function(x)
  data.frame(x[1, c("range", "smoothness", "n", "rho")], dist_same = mean(x$distance[x$same_field]),
             dist_independent = mean(x$distance[!x$same_field]),
             auc = auc(x$distance[!x$same_field], x$distance[x$same_field]))))
write.csv(d11_sum, file.path("results", "d11_summary.csv"), row.names = FALSE)
# threshold: largest Youden index (sensitivity + specificity - 1) over all conditions
ths <- seq(0.2, 0.9, by = 0.01)
youden <- vapply(ths, function(t) mean(d11$distance[d11$same_field] <= t) + mean(d11$distance[!d11$same_field] > t) - 1, 0)
write.csv(data.frame(threshold = ths, youden = youden,
                     same_kept_together = vapply(ths, function(t) mean(d11$distance[d11$same_field] <= t), 0),
                     independent_split = vapply(ths, function(t) mean(d11$distance[!d11$same_field] > t), 0)),
          file.path("results", "d11_threshold.csv"), row.names = FALSE)
cat(format(Sys.time()), "summaries written; best threshold", ths[which.max(youden)], "\n", file = log_file, append = TRUE)
print(h2a_sum, digits = 3); print(h2b_sum, digits = 3); print(d11_sum, digits = 3)
