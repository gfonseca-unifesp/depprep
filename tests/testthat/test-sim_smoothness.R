test_that("the Matern correlation has the practical range and reduces to the exponential", {
  for (nu in c(0.5, 1.5, 2.5)) expect_equal(depprep:::.matern_cor(0.3, 0.3, nu), 0.05, tolerance = 0.01)
  h <- seq(0, 1, by = 0.1)
  expect_identical(depprep:::.matern_cor(h, 0.4, 0.5), depprep:::.exp_cor(h, 0.4))
  # smoother fields are more correlated at short distances, for the same practical range
  expect_gt(depprep:::.matern_cor(0.05, 0.4, 2.5), depprep:::.matern_cor(0.05, 0.4, 0.5))
  expect_error(depprep:::.matern_cor(0.1, 0.4, 1), "smoothness")
})

test_that("Cholesky and Fourier-feature fields have the Matern correlation", {
  coords <- cbind(c(0, 0.1, 0.3), c(0, 0, 0))
  set.seed(1)
  z <- depprep:::.grf_chol(coords, range = 0.4, n_sim = 4000, smoothness = 2.5)
  # absolute differences: the standard error of a correlation from 4000 draws is <= 0.016
  expect_lt(abs(cor(z[1, ], z[2, ]) - depprep:::.matern_cor(0.1, 0.4, 2.5)), 0.05)
  expect_lt(abs(cor(z[1, ], z[3, ]) - depprep:::.matern_cor(0.3, 0.4, 2.5)), 0.05)
  set.seed(2)
  f <- replicate(2000, depprep:::.rff_sampler(0.4, smoothness = 2.5)(coords))
  expect_lt(abs(cor(f[1, ], f[2, ]) - depprep:::.matern_cor(0.1, 0.4, 2.5)), 0.06)
})

test_that("the matern regime is smoother and the exponential default is unchanged", {
  g_m <- sim_groups(sim_group("m", 1, "matern", range = 0.3, large_scale_share = 1))
  g_e <- sim_groups(sim_group("e", 1, "exponential", range = 0.3, large_scale_share = 1))
  nn_cor <- function(s) {
    d <- as.matrix(dist(s$coords)); diag(d) <- Inf
    cor(s$X[, 1], s$X[apply(d, 1, which.min), 1])
  }
  sm <- mean(vapply(1:10, function(r) nn_cor(sim_spatial(300, g_m, seed = r)), 0))
  se <- mean(vapply(1:10, function(r) nn_cor(sim_spatial(300, g_e, seed = r)), 0))
  expect_gt(sm, se + 0.05)                                  # neighbours more alike in the smooth field
  g_e2 <- sim_groups(sim_group("e", 1, "exponential", range = 0.3, smoothness = 0.5, large_scale_share = 1))
  expect_identical(sim_spatial(100, g_e, seed = 3)$X, sim_spatial(100, g_e2, seed = 3)$X)
})

test_that("spatio-temporal regimes accept smoothness", {
  g <- sim_groups(sim_group("s", 2, "separable_ar1", range = 0.5, phi = 0.9, smoothness = 2.5),
                  sim_group("m", 2, "moving", range = 0.25, velocity = c(0.02, 0), smoothness = 2.5),
                  sim_group("q", 2, "separable_seasonal", range = 0.25, period = 12, smoothness = 2.5))
  s <- sim_spatiotemporal(30, 20, g, seed = 6)
  expect_equal(dim(s$X), c(600, 6))
  expect_false(anyNA(s$X))
  expect_identical(sim_spatiotemporal(30, 20, g, seed = 6)$X, s$X)
})
