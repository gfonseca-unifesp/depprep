groups_t <- function() {
  sim_groups(sim_group("strong", 3, "ar1", phi = 0.9),
             sim_group("weak", 3, "ar1", phi = 0.3),
             sim_group("season", 2, "seasonal", period = 12),
             sim_group("trend", 2, "trend_shift", shift_at = 0.5),
             sim_group("none", 2, "none"))
}

test_that("sim_temporal returns the documented structure", {
  s <- sim_temporal(120, groups_t(), seed = 1)
  expect_s3_class(s, "depprep_sim")
  expect_equal(dim(s$X), c(120L, 12L))
  expect_equal(s$time, 1:120)
  expect_equal(s$family, "temporal")
  rebuilt <- sweep(s$components$large, 2, sqrt(s$large_scale_share), "*") +
    sweep(s$components$residual, 2, sqrt(1 - s$large_scale_share), "*")
  expect_equal(s$X, rebuilt)
})

test_that("same seed gives identical data and invalid inputs fail clearly", {
  expect_identical(sim_temporal(60, groups_t(), seed = 3)$X, sim_temporal(60, groups_t(), seed = 3)$X)
  expect_error(sim_temporal(60, sim_groups(sim_group("a", 1, "ar1")), seed = 1), "needs `phi`")
  expect_error(sim_temporal(60, sim_groups(sim_group("s", 1, "seasonal")), seed = 1), "needs `period`")
  expect_error(sim_temporal(60, sim_groups(sim_group("e", 1, "exponential", range = 1)), seed = 1), "not available")
})

test_that("phi of the AR(1) groups is recovered (mean of 50 replicates within 0.05)", {
  est <- t(vapply(1:50, function(i) {
    L <- sim_temporal(480, groups_t(), seed = 200 + i)$components$large
    c(stats::acf(L[, "strong_1"], 1, plot = FALSE)$acf[2], stats::acf(L[, "weak_1"], 1, plot = FALSE)$acf[2])
  }, numeric(2)))
  expect_lt(abs(mean(est[, 1]) - 0.9), 0.05)
  expect_lt(abs(mean(est[, 2]) - 0.3), 0.05)
})

test_that("the seasonal group peaks at its period and the trend jumps where specified", {
  s <- sim_temporal(480, groups_t(), seed = 4)
  sp <- stats::spec.pgram(s$components$large[, "season_1"], taper = 0, plot = FALSE)
  expect_equal(1 / sp$freq[which.max(sp$spec)], 12)
  expect_equal(which.max(diff(s$components$large[, "trend_1"])), 240L)
})

test_that("contiguous classes are blocks of time", {
  s <- sim_temporal(90, groups_t(), classes = list(k = 3, type = "contiguous"), seed = 5)
  expect_equal(s$classes, rep(1:3, each = 30))
})

test_that("residual components recover R_true", {
  g <- groups_t()
  r <- make_R_true(g, "blocks", within = 0.4, between = 0.1)
  reps <- lapply(1:50, function(i) stats::cor(sim_temporal(480, g, R_true = r, seed = i)$components$residual))
  expect_lt(max(abs(Reduce(`+`, reps) / 50 - r)), 0.05)
})
