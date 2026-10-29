tm_series <- function(regime, seed, n_times = 480, share = 1, ...) {
  g <- sim_groups(sim_group("f", 1, regime, ..., large_scale_share = share))
  sim_temporal(n_times, g, seed = seed)$X[, 1]
}

test_that(".temporal_stats recovers phi (mean of 30 replicates within 0.05)", {
  skip_on_cran()
  for (phi in c(0.3, 0.9)) {
    est <- vapply(1:30, function(i) depprep:::.temporal_stats(tm_series("ar1", 600 + i, phi = phi))[["phi"]], 0)
    expect_lt(abs(mean(est) - phi), 0.05)
  }
})

test_that(".temporal_stats finds the seasonal period and the trend", {
  s <- depprep:::.temporal_stats(tm_series("seasonal", 7, period = 12, share = 0.6))
  expect_equal(s[["period"]], 12)
  expect_gt(s[["period_strength"]], 0.2)
  tr <- depprep:::.temporal_stats(tm_series("trend_shift", 8, share = 0.6))
  ar <- depprep:::.temporal_stats(tm_series("ar1", 8, phi = 0.3, share = 0.6))
  expect_gt(tr[["trend_r2"]], 0.4)
  expect_lt(ar[["trend_r2"]], 0.1)
})

test_that(".temporal_stats finds a period off the Fourier grid and accepts a known period", {
  x <- 2 * sin(2 * pi * (1:100) / 12.5) + depprep:::.with_stream(3, 1, stats::rnorm(100))
  expect_lt(abs(depprep:::.temporal_stats(x)[["period"]] - 12.5), 0.2)
  expect_equal(depprep:::.temporal_stats(x, period = 12)[["period"]], 12)
})

test_that(".temporal_stats measures persistence after removing trend and cycle", {
  e <- depprep:::.with_stream(5, 1, as.numeric(stats::arima.sim(list(ar = 0.5), 480)))
  y <- 1.5 * sin(2 * pi * (1:480) / 12) + 0.01 * (1:480) + e
  s <- depprep:::.temporal_stats(y)
  expect_lt(abs(s[["period"]] - 12), 0.1)
  expect_lt(abs(s[["phi"]] - 0.5), 0.15)      # the trend does not push phi towards 1
  expect_gt(s[["trend_r2"]], 0.4)             # the trend is seen through the cycle
})

test_that("the period test rarely flags white noise", {
  p <- vapply(1:100, function(i) depprep:::.temporal_stats(depprep:::.with_stream(i, 1, stats::rnorm(120)))[["period_p"]], 0)
  expect_lt(mean(p < 0.05), 0.12)
})

test_that("Ljung-Box flags persistence and rarely flags white noise", {
  expect_lt(depprep:::.temporal_stats(tm_series("ar1", 9, phi = 0.9, share = 0.6))[["ljung_p"]], 0.001)
  p <- vapply(1:200, function(i) depprep:::.temporal_stats(depprep:::.with_stream(i, 1, stats::rnorm(120)))[["ljung_p"]], 0)
  expect_lt(mean(p < 0.05), 0.1)
})

st_sim <- function(regime, seed, ...) {
  g <- sim_groups(sim_group("f", 1, regime, ..., large_scale_share = 0.6))
  sim_spatiotemporal(100, 60, g, seed = seed)
}

test_that(".st_stats separates moving, separable, temporal-only and no dependence", {
  skip_on_cran()
  st <- function(s) depprep:::.st_stats(s$X[, 1], s$coords, s$index)
  mov <- st(st_sim("moving", 1, range = 0.25, velocity = c(0.02, 0)))
  sep <- st(st_sim("separable_ar1", 1, range = 0.25, phi = 0.9))
  com <- st(st_sim("temporal_common", 1, phi = 0.9))
  non <- st(st_sim("none", 1))
  expect_gt(mov[["shift_indicator"]], sep[["shift_indicator"]])
  expect_gt(mov[["shift_indicator"]], 0)
  expect_equal(mov[["shift_direction"]], 0)               # moving towards +x
  expect_lt(sep[["moran_p"]], 0.05)
  expect_gt(com[["moran_p"]], 0.05)                        # no spatial structure in a slice
  expect_lt(com[["ljung_p"]], 0.05)
  expect_gt(non[["moran_p"]], 0.05)
  expect_gt(non[["ljung_p"]], 0.05)
})
