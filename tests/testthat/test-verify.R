test_that("verify reports removed variance, autocorrelation and Delta R", {
  g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(300, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 1)
  grp <- stats::setNames(c(rep("g1", 4), rep("none", 2)), colnames(s$X))
  b <- list(g1 = build_group_basis(s$X[, 1:4], coords = s$coords, type = "isotropic", seed = 1))
  v <- verify(correct(s, grp, bases = b), s, nperm = 0)
  expect_s3_class(v, "depprep_verification")
  long <- v$variables$group == "g1"
  expect_true(all(v$variables$removed[long] > 0.3))       # the shared field is mostly removed
  expect_true(all(v$variables$autocor_after[long] < v$variables$autocor_before[long]))
  expect_setequal(v$delta_summary$kind, c("within groups", "involving no dependence"))
  # the apparent association inside the group falls
  expect_lt(v$delta_summary$mean_delta[v$delta_summary$kind == "within groups"], 0)
  expect_equal(v$delta_R, t(v$delta_R))
})

# Large scale exactly in the span of `B`, shared by the variables, plus white noise.
exact_field <- function(B, p, seed) {
  set.seed(seed)
  n <- nrow(B)
  L <- as.numeric(B %*% stats::rnorm(ncol(B))); L <- L / stats::sd(L)
  X <- sqrt(0.6) * L + sqrt(0.4) * matrix(stats::rnorm(n * p), n)
  colnames(X) <- paste0("v", seq_len(p))
  X
}
as_basis <- function(B) structure(B, class = "depprep_basis", type = "isotropic", k = ncol(B), kn = ncol(B) / nrow(B),
                                  terms = colnames(B))

test_that("data corrected with the true basis show no residual structure, and too small a basis does", {
  set.seed(11)
  xy <- matrix(stats::runif(400), 200)
  cand <- depprep:::.spatial_candidates(xy)
  X <- exact_field(cand[, 1:6], 6, seed = 12)
  grp <- stats::setNames(rep("g1", 6), colnames(X))
  v <- verify(correct(X, grp, coords = xy, bases = list(g1 = as_basis(cand[, 1:6]))), X, coords = xy, nperm = 99)
  expect_lte(sum(v$variables$residual_structure), 1)
  expect_gt(v$mantel_p, 0.05)
  # the correlogram uses the classes of vegan::mantel.correlog and the moments of residuals
  expect_equal(unique(v$correlograms$distance), v$mantel$distance)
  expect_true(all(v$correlograms$expected[v$correlograms$class == 1] < -1 / 199))
  # one eigenvector leaves most of the field: the correlograms and the Mantel test see it
  w <- verify(correct(X, grp, coords = xy, bases = list(g1 = as_basis(cand[, 1, drop = FALSE]))), X, coords = xy,
              nperm = 499)
  expect_true(all(w$variables$residual_structure))
  expect_true(all(w$variables$extent >= 1))
  expect_lt(w$mantel_p, 0.05)
})

test_that("in time, the Ljung-Box null accounts for the basis that was removed", {
  n <- 240
  time <- seq_len(n)
  cand <- depprep:::.temporal_candidates(time, "seasonal", 12)
  X <- exact_field(cand, 5, seed = 21)
  grp <- stats::setNames(rep("g1", 5), colnames(X))
  v <- verify(correct(X, grp, time = time, bases = list(g1 = as_basis(cand))), X, time = time, nperm = 99)
  expect_equal(unique(v$variables$lb_lag), 24)               # max(2 * 12, 10 * log10(240)), period from the terms
  expect_lte(sum(v$variables$residual_structure), 1)
  # many broad eigenfunctions removed from white noise: the chi-square reference would
  # flag the induced negative autocorrelation, the projected null does not
  set.seed(22)
  W <- matrix(stats::rnorm(n * 8), n, dimnames = list(NULL, paste0("w", 1:8)))
  mem <- depprep:::.temporal_candidates(time, "ar")[, 1:60]
  cw <- correct(W, stats::setNames(rep("g1", 8), colnames(W)), time = time, bases = list(g1 = as_basis(mem)))
  vw <- verify(cw, W, time = time, nperm = 199)
  chi <- verify(cw, W, time = time, nperm = 0)
  expect_lte(sum(vw$variables$residual_structure), 1)
  expect_gt(sum(chi$variables$residual_structure), sum(vw$variables$residual_structure))
  # without the seasonal basis the cycle is left and is detected
  u <- verify(correct(X, grp, time = time, bases = list(g1 = as_basis(mem[, 1:2]))), X, time = time, nperm = 0, period = 12)
  expect_true(all(u$variables$residual_structure))
})

test_that("the projected and residual parts rebuild the covariance", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("b", 3, "matern", range = 0.15),
                  sim_group("none", 2, "none"))
  s <- sim_spatial(250, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0.2), seed = 3)
  grp <- stats::setNames(c(rep("g1", 3), rep("g2", 3), rep("none", 2)), colnames(s$X))
  bases <- list(g1 = build_group_basis(s$X[, 1:3], coords = s$coords, type = "isotropic", seed = 1),
                g2 = build_group_basis(s$X[, 4:6], coords = s$coords, type = "isotropic", seed = 1))
  cr <- correct(s, grp, bases = bases, data_basis = "union")   # the identities hold for the union projector
  v <- verify(cr, s, nperm = 0)
  C <- stats::cov(s$X)
  expect_equal(v$cov_large + v$cov_residual, C, ignore_attr = TRUE)
  # R_large and the residual correlation, rescaled by the standard deviations of each part
  sd_L <- sqrt(diag(v$cov_large)); sd_E <- sqrt(diag(v$cov_residual))
  R_E <- stats::cov2cor(v$cov_residual)
  expect_equal(diag(sd_L) %*% v$R_large %*% diag(sd_L) + diag(sd_E) %*% R_E %*% diag(sd_E), C, ignore_attr = TRUE)
  # the residual part is that of the corrected data matrix
  expect_equal(R_E, stats::cor(cr$X), ignore_attr = TRUE)
  expect_equal(v$R_large, stats::cor(s$X - cr$X - rep(colMeans(s$X), each = nrow(s$X))), ignore_attr = TRUE)
  expect_equal(unname(diag(v$large_share)), v$variables$removed)
  # pairwise shares use the projector of each pair in R: group basis within, union across
  Xc <- scale(s$X, scale = FALSE)
  B1 <- unclass(bases$g1); U <- unclass(union_basis(bases$g1, bases$g2))
  expect_equal(v$large_share_pairwise[1, 2], sum(Xc[, 1] * B1 %*% crossprod(B1, Xc[, 2])) / sum(Xc[, 1] * Xc[, 2]))
  expect_equal(v$large_share_pairwise[1, 4], sum(Xc[, 1] * U %*% crossprod(U, Xc[, 4])) / sum(Xc[, 1] * Xc[, 4]))
  expect_equal(v$large_share_pairwise[7, 8], 0)
  # whitening has no projected part
  vw <- verify(correct(s, grp, bases = bases, mode = "whiten"), s, nperm = 0)
  expect_true(all(is.na(vw$R_large)))
  expect_true(is.na(vw$q3$statistic))
  expect_false(vw$q3$warning)
  expect_match(vw$notes, "no projected part")
})

test_that("the Q3 warning fires for contiguous classes and not for the null truth", {
  skip_on_cran()
  g <- sim_groups(sim_group("long", 5, "exponential", range = 0.5), sim_group("short", 5, "exponential", range = 0.15),
                  sim_group("aniso", 5, "exponential", range = 0.25, anisotropy = 4, angle = 30),
                  sim_group("none", 5, "none"))
  q3 <- function(truth) {
    R_true <- if (truth == "null") NULL else make_R_true(g, "blocks", within = 0.5, between = 0)
    cls <- if (truth == "null") NULL else list(k = 3, type = "contiguous", effect = 1)
    s <- sim_spatial(250, g, R_true = R_true, classes = cls, seed = 1)
    grp <- stats::setNames(ifelse(s$group == "none", "none", s$group), colnames(s$X))
    verify(correct(s, grp, seed = 1), s, nperm = 0)
  }
  v0 <- q3("null"); v1 <- q3("contiguous")
  expect_false(v0$q3$warning)
  expect_true(v1$q3$warning)
  expect_match(v1$warnings, "signal at the scale of the basis", all = FALSE)
  expect_false(any(grepl("signal at the scale of the basis", v0$warnings)))
  expect_gt(v1$q3$statistic, v0$q3$statistic)
})

test_that("the LISA finds a region left by the correction and nothing in white noise", {
  set.seed(7)
  xy <- cbind(runif(200), runif(200))
  patch <- as.numeric((xy[, 1] - 0.8)^2 + (xy[, 2] - 0.8)^2 < 0.02)
  Z <- cbind(noise = rnorm(200), patch = 3 * patch + rnorm(200, sd = 0.5))
  l <- depprep:::.lisa(Z, xy)
  expect_equal(nrow(l), 400)
  expect_lt(mean(l$cluster[l$variable == "noise"] != "none"), 0.03)
  hh <- l$variable == "patch" & l$cluster == "high_high"
  expect_gt(mean(patch[l$site[hh]]), 0.8)                       # the high-high sites are in the patch
  expect_gt(sum(hh), 5)
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(150, g, seed = 3)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  v <- verify(correct(s, grp, seed = 1), s, nperm = 0)
  expect_true(all(c("local_I", "p_adjusted", "cluster") %in% names(v$lisa)))
  expect_true(all(v$variables$lisa_share >= 0 & v$variables$lisa_share <= 1))
  expect_output(print(v), "LISA")
  gt <- sim_groups(sim_group("a", 2, "ar1", phi = 0.8))
  st <- sim_temporal(100, gt, seed = 1)
  grt <- data.frame(variable = colnames(st$X), group = "g1", group_type = "ar")
  class(grt) <- c("depprep_grouping", "data.frame"); attr(grt, "period") <- c(g1 = NA)
  expect_null(verify(correct(st, grt), st, nperm = 0)$lisa)
})

test_that("spectral surrogates keep each variable's spectrum and break the association", {
  set.seed(2)
  coords <- cbind(runif(120), runif(120))
  f <- sin(3 * coords[, 1]) + cos(4 * coords[, 2])
  X <- cbind(a = f + rnorm(120, sd = 0.3), b = f + rnorm(120, sd = 0.3))
  sb <- depprep:::.surrogate_bases("spatial", coords = coords)
  s <- depprep:::.with_stream(1, 11L, depprep:::.spectral_surrogate(X, sb))
  Xc <- scale(X, scale = FALSE)
  expect_equal(crossprod(sb$U, s)^2, crossprod(sb$U, Xc)^2, tolerance = 1e-8)   # same power at every scale
  expect_equal(colSums(s^2), colSums(Xc^2), tolerance = 1e-8)
  r <- depprep:::.with_stream(1, 11L, vapply(1:200, function(b) stats::cor(depprep:::.spectral_surrogate(X, sb))[1, 2], 0))
  expect_gt(stats::cor(X)[1, 2], 0.8)
  expect_lt(abs(mean(r)), 0.1)                                                 # association broken on average
  # temporal: the lag-1 autocorrelation is kept on average
  x <- as.numeric(stats::arima.sim(list(ar = 0.8), 200))
  sbt <- depprep:::.surrogate_bases("temporal", time = 1:200)
  a1 <- function(z) stats::acf(z, lag.max = 1, plot = FALSE)$acf[2]
  st <- depprep:::.with_stream(1, 11L, vapply(1:50, function(b) a1(depprep:::.spectral_surrogate(cbind(x), sbt)), 0))
  expect_lt(abs(mean(st) - a1(x)), 0.1)
  # space-time: the product basis is complete; an incomplete grid stops
  idx <- data.frame(site = rep(1:6, times = 5), time = rep(1:5, each = 6))
  sbs <- depprep:::.surrogate_bases("spatiotemporal", coords = coords[1:6, ], index = idx)
  z <- cbind(rnorm(30))
  zs <- depprep:::.with_stream(1, 11L, depprep:::.spectral_surrogate(z, sbs))
  expect_equal(sum(zs^2), sum(scale(z, scale = FALSE)^2), tolerance = 1e-8)
  expect_error(depprep:::.surrogate_bases("spatiotemporal", coords = coords[1:6, ], index = idx[-1, ]), "complete grid")
})

test_that("the Q3 reference matches the statistic for independent white noise and is reproducible", {
  set.seed(3)
  n <- 150
  X <- matrix(rnorm(n * 8), n)
  U <- qr.Q(qr(scale(matrix(rnorm(n * 15), n), scale = FALSE)))
  sb <- depprep:::.surrogate_bases("temporal", time = seq_len(n))
  obs <- depprep:::.large_association(crossprod(U %*% crossprod(U, scale(X, scale = FALSE))) / (n - 1),
                                      crossprod(scale(X, scale = FALSE) - U %*% crossprod(U, scale(X, scale = FALSE))) / (n - 1))[["mean"]]
  ref <- depprep:::.q3_reference(X, U, sb, B = 49, seed = 5)
  expect_length(ref$values, 49)
  expect_lt(abs(obs / ref$mean - 1), 0.3)
  expect_identical(ref$values, depprep:::.q3_reference(X, U, sb, B = 49, seed = 5)$values)
})

test_that("the default Q3 threshold depends on the data family and an explicit one overrides it", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(150, g, seed = 3)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  cr <- correct(s, grp, seed = 1)
  expect_equal(verify(cr, s, nperm = 0)$q3$threshold, 0.16)
  v <- verify(cr, s, nperm = 0, q3_threshold = 0.01)
  expect_equal(v$q3$threshold, 0.01)
  expect_true(v$q3$warning)
  gt <- sim_groups(sim_group("a", 3, "ar1", phi = 0.9), sim_group("none", 2, "none"))
  st <- sim_temporal(120, gt, seed = 3)
  gr <- data.frame(variable = colnames(st$X), group = c("g1", "g1", "g1", "none", "none"),
                   group_type = c("ar", "ar", "ar", "none", "none"))
  class(gr) <- c("depprep_grouping", "data.frame"); attr(gr, "period") <- c(g1 = NA)
  ct <- correct(st, gr, seed = 1)
  expect_equal(verify(ct, st, nperm = 0)$q3$threshold, 0.19)
  expect_equal(unname(depprep:::.q3_thresholds[c("spatial", "temporal", "spatiotemporal")]), c(0.16, 0.19, 0.13))
})

test_that("verify() tests each column of an own-basis data matrix against its own basis", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(200, g, seed = 6)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  cr <- correct(s, grp, seed = 1)                          # default: data_basis = "own"
  v <- verify(cr, s, nperm = 49)
  expect_equal(nrow(v$variables), 5)
  expect_true(all(is.finite(v$variables$correlogram_p)))
  expect_lte(sum(v$variables$residual_structure), 1)       # each null removes the variable's own basis
  ob <- depprep:::.own_bases(cr, nrow(s$X), colnames(s$X))
  expect_equal(ncol(ob[[1]]), ncol(cr$variable_bases[[1]]))
  expect_equal(ncol(ob[[4]]), 0)
})
