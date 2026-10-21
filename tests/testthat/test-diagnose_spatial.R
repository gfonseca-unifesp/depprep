sp_field <- function(range, seed, n = 400, anisotropy = 1, angle = 0, share = 1) {
  g <- sim_groups(sim_group("f", 1, "exponential", range = range, anisotropy = anisotropy,
                            angle = angle, large_scale_share = share))
  sim_spatial(n, g, seed = seed)
}

test_that(".moran detects spatial structure and not white noise", {
  s <- sp_field(0.25, seed = 1)
  lw <- depprep:::.knn_weights(s$coords, k = 8)
  m_field <- depprep:::.moran(s$X[, 1], lw)
  noise <- depprep:::.with_stream(2, 1, stats::rnorm(400))
  m_noise <- depprep:::.moran(noise, lw)
  expect_gt(m_field[["moran_I"]], 0.3)
  expect_lt(m_field[["moran_p"]], 0.001)
  expect_lt(abs(m_noise[["moran_I"]]), 0.1)
})

test_that("Moran p-values are roughly uniform without dependence", {
  xy <- depprep:::.with_stream(3, 1, matrix(stats::runif(400), ncol = 2))
  lw <- depprep:::.knn_weights(xy, k = 8)
  p <- vapply(1:200, function(i) depprep:::.moran(depprep:::.with_stream(i, 2, stats::rnorm(200)), lw)[["moran_p"]], 0)
  expect_equal(mean(p < 0.05), 0.05, tolerance = 0.6)   # about 5% false positives
})

test_that(".fit_variogram recovers the range of the simulated fields", {
  skip_on_cran()
  for (rng in c(0.1, 0.25, 0.5)) {
    est <- vapply(1:30, function(i) {
      s <- sp_field(rng, seed = 300 + i)
      depprep:::.fit_variogram(s$coords, s$X[, 1])$range
    }, numeric(1))
    expect_lt(abs(mean(est, na.rm = TRUE) / rng - 1), 0.25)
  }
})

test_that(".fit_variogram returns a nugget for noisy variables", {
  s <- sp_field(0.25, seed = 4, share = 0.6)
  f <- depprep:::.fit_variogram(s$coords, s$X[, 1])
  expect_true(f$vario_model %in% c("Exp", "Sph"))
  expect_gt(f$nugget_ratio, 0.15)
  expect_lt(f$nugget_ratio, 0.7)
})

test_that(".anisotropy separates anisotropic from isotropic fields", {
  skip_on_cran()
  ratio <- function(a, seed) {
    s <- sp_field(0.25, seed = seed, anisotropy = a, angle = 30)
    depprep:::.anisotropy(s$coords, s$X[, 1])[["aniso_ratio"]]
  }
  iso <- vapply(1:20, function(i) ratio(1, 500 + i), 0)
  ani <- vapply(1:20, function(i) ratio(4, 500 + i), 0)
  expect_gt(stats::median(ani, na.rm = TRUE), stats::median(iso, na.rm = TRUE) * 1.5)
})

test_that(".fit_variogram returns NA values, not an error, when nothing converges", {
  xy <- cbind(c(0, 1, 0, 1, 0.5), c(0, 0, 1, 1, 0.5))
  f <- depprep:::.fit_variogram(xy, c(1, 1, 1, 1, 1.0001))
  expect_type(f, "list")
  expect_true(is.na(f$range))
})
