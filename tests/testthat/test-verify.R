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
  cr <- correct(s, grp, bases = bases)
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
