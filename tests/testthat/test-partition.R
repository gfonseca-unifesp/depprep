known_partition <- function(n = 2000, seed = 1) {
  # orthonormal patterns: q1 in both tables (shared), q2 in the position only, q6 in the
  # drivers only; unit-variance parts 0.2, 0.3 and 0.1, and noise 0.4
  set.seed(seed)
  Q <- qr.Q(qr(scale(matrix(stats::rnorm(n * 6), n), scale = FALSE))) * sqrt(n - 1)
  P <- Q[, 1:5]
  D <- cbind(d1 = Q[, 1] + 0.5 * Q[, 6], d2 = Q[, 6])
  x <- sqrt(0.2) * Q[, 1] + sqrt(0.3) * Q[, 2] + sqrt(0.1) * Q[, 6] + sqrt(0.4) * stats::rnorm(n)
  list(x = cbind(x = x, y = stats::rnorm(n)), P = P, D = D)
}

test_that("the partition recovers a known split of the variance", {
  d <- known_partition()
  vp <- variation_partition(d$x, d$P, d$D)
  expect_s3_class(vp, "data.frame")
  expect_identical(vp$variable, c("x", "y"))
  expect_identical(vp$k_position, c(5L, 5L))
  expect_identical(vp$k_drivers, c(2L, 2L))
  truth <- c(position_only = 0.3, shared = 0.2, drivers_only = 0.1, residual = 0.4)
  expect_lt(max(abs(unlist(vp[1, names(truth)]) - truth)), 0.03)
  expect_lt(max(abs(unlist(vp[2, names(truth)[1:3]]))), 0.01)       # pure noise: nothing explained
  expect_equal(rowSums(vp[, names(truth)]), c(1, 1))
})

test_that("the fractions come from the adjusted R-squared of the three regressions", {
  d <- known_partition(n = 300, seed = 2)
  vp <- variation_partition(d$x[, 1, drop = FALSE], d$P, d$D)
  adj <- function(M) summary(stats::lm(d$x[, 1] ~ M))$adj.r.squared
  a <- adj(d$P); b <- adj(d$D); ab <- adj(cbind(d$P, d$D[, 2]))   # d1 lies in the span of the others
  expect_equal(unlist(vp[1, c("position_only", "shared", "drivers_only", "residual")]),
               c(ab - b, a + b - ab, ab - a, 1 - ab), ignore_attr = TRUE)
})

test_that("one basis per variable, or none, is accepted", {
  d <- known_partition(n = 300, seed = 3)
  both <- variation_partition(d$x, list(y = d$P[, 1:2], x = d$P), d$D)
  expect_equal(both[1, ], variation_partition(d$x[, "x", drop = FALSE], d$P, d$D), ignore_attr = TRUE)
  expect_identical(both$k_position, c(5L, 2L))
  none <- variation_partition(d$x, NULL, as.data.frame(d$D))
  expect_equal(none$position_only, c(0, 0))
  expect_equal(none$shared, c(0, 0))
  expect_error(variation_partition(d$x, list(d$P), d$D), "one basis per column")
  expect_error(variation_partition(d$x, d$P, d$D[-1, ]), "one row per observation")
  expect_error(variation_partition(d$x, d$P, replace(d$D, 1, NA)), "missing")
})

test_that("a variable with gaps is partitioned on its observed rows", {
  set.seed(11)
  t <- 1:120
  position <- cbind(sin(2 * pi * t / 12), cos(2 * pi * t / 12), t / 120)
  driver <- cbind(index = stats::rnorm(120))
  x <- 2 * position[, 1] + driver[, 1] + stats::rnorm(120)
  gone <- c(3, 20:31, 77)
  y <- x
  y[gone] <- NA
  obs <- -gone
  with_gaps <- variation_partition(cbind(x = y, full = x), position, driver)
  alone <- variation_partition(cbind(x = x[obs]), position[obs, ], driver[obs, , drop = FALSE])
  cols <- c("position_only", "shared", "drivers_only", "residual")
  expect_equal(unlist(with_gaps[1, cols]), unlist(alone[1, cols]))
  expect_equal(unlist(with_gaps[2, cols]), unlist(variation_partition(cbind(full = x), position, driver)[1, cols]))
  expect_equal(sum(with_gaps[1, cols]), 1)
})
