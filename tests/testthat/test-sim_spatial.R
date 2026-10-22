groups_s <- function() {
  sim_groups(sim_group("long", 3, "exponential", range = 0.25),
             sim_group("short", 3, "exponential", range = 0.1, anisotropy = 2, angle = 30),
             sim_group("copy", 2, "exponential", shared_with = "long", large_scale_share = 0.4),
             sim_group("none", 2, "none"))
}

test_that("sim_spatial returns the documented structure", {
  s <- sim_spatial(60, groups_s(), seed = 1)
  expect_s3_class(s, "depprep_sim")
  expect_equal(dim(s$X), c(60L, 10L))
  expect_equal(dim(s$coords), c(60L, 2L))
  expect_true(all(s$coords >= 0 & s$coords <= 1))
  expect_named(s$components, c("large", "residual", "class"))
  expect_equal(unname(s$group), attr(groups_s(), "var_group"))
  expect_equal(s$seed, 1)
  expect_output(print(s), "spatial")
})

test_that("X is the sum of its components with the right shares", {
  s <- sim_spatial(50, groups_s(), seed = 2)
  share <- s$large_scale_share
  rebuilt <- sweep(s$components$large, 2, sqrt(share), "*") +
    sweep(s$components$residual, 2, sqrt(1 - share), "*") + s$components$class
  expect_equal(s$X, rebuilt)
  expect_equal(unname(share[c("long_1", "copy_1", "none_1")]), c(0.6, 0.4, 0))
  expect_true(all(s$components$large[, "none_1"] == 0))
})

test_that("variables of a group, and groups that share a field, have the same large-scale component", {
  s <- sim_spatial(50, groups_s(), seed = 3)
  L <- s$components$large
  expect_equal(L[, "long_1"], L[, "long_3"])
  expect_equal(L[, "long_1"], L[, "copy_2"])
  expect_false(isTRUE(all.equal(L[, "long_1"], L[, "short_1"])))
})

test_that("same seed gives identical data and different seeds do not", {
  a <- sim_spatial(40, groups_s(), seed = 10)
  b <- sim_spatial(40, groups_s(), seed = 10)
  c <- sim_spatial(40, groups_s(), seed = 11)
  expect_identical(a$X, b$X)
  expect_false(isTRUE(all.equal(a$X, c$X)))
})

test_that("given coordinates are kept and invalid inputs fail clearly", {
  xy <- cbind(seq(0, 1, length.out = 30), 0.5)
  s <- sim_spatial(groups = groups_s(), coords = xy, seed = 4)
  expect_equal(unname(s$coords), xy)
  expect_error(sim_spatial(30, groups_s(), R_true = diag(3), seed = 1), "10 x 10")
  bad <- sim_groups(sim_group("t", 2, "ar1", phi = 0.5))
  expect_error(sim_spatial(30, bad, seed = 1), "not available")
  nr <- sim_groups(sim_group("x", 2, "exponential"))
  expect_error(sim_spatial(30, nr, seed = 1), "needs `range`")
})

test_that("residual components recover R_true", {
  g <- groups_s()
  r <- make_R_true(g, "blocks", within = 0.5, between = 0.2)
  reps <- lapply(1:50, function(i) stats::cor(sim_spatial(400, g, R_true = r, seed = i)$components$residual))
  mean_cor <- Reduce(`+`, reps) / length(reps)
  expect_lt(max(abs(mean_cor - r)), 0.05)
  # one replicate: Frobenius error within 1.5 times its sampling expectation,
  # E||R_hat - R||_F^2 ~ sum_{i != j} (1 - rho_ij^2)^2 / n
  expected <- sqrt(sum(((1 - r^2)^2)[upper.tri(r)]) * 2 / 400)
  expect_lt(norm(reps[[1]] - r, "F"), 1.5 * expected)
})

test_that("classes are attached as specified", {
  s <- sim_spatial(200, groups_s(), classes = list(k = 3, type = "contiguous", effect = 2, vars = 1:3), seed = 5)
  expect_setequal(unique(s$classes), 1:3)
  expect_true(all(s$components$class[, 4:10] == 0))
  expect_gt(stats::sd(s$components$class[, 1]), 0)
})

test_that("variogram range of the simulated fields matches the simulated range", {
  skip_on_cran()
  for (rng in c(0.1, 0.25)) {
    g <- sim_groups(sim_group("f", 1, "exponential", range = rng, large_scale_share = 1))
    est <- vapply(1:50, function(i) {
      s <- sim_spatial(400, g, seed = 100 + i)
      d <- data.frame(s$coords, z = s$components$large[, 1])
      v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = d, cutoff = 0.5)
      fit <- suppressWarnings(gstat::fit.variogram(v, gstat::vgm(1, "Exp", rng / 3)))
      3 * fit$range[nrow(fit)]
    }, numeric(1))
    expect_lt(abs(mean(est) / rng - 1), 0.25)
  }
})

test_that("the gradient regime is the specified surface, exactly without a field", {
  g <- sim_groups(sim_group("plane", 2, "gradient", direction = 30),
                  sim_group("bowl", 1, "gradient", direction = 0, shape = "quadratic"),
                  sim_group("none", 1, "none"))
  s <- sim_spatial(100, g, seed = 6)
  L <- s$components$large
  expect_equal(L[, "plane_1"], L[, "plane_2"])
  fit <- stats::lm(L[, "plane_1"] ~ s$coords)
  expect_lt(max(abs(stats::residuals(fit))), 1e-10)
  b <- stats::coef(fit)[-1]
  expect_equal(atan2(b[2], b[1]) * 180 / pi, 30, ignore_attr = TRUE)
  expect_equal(c(mean(L[, "bowl_1"]), stats::sd(L[, "bowl_1"])), c(0, 1), tolerance = 1e-12)
  expect_equal(stats::cor(L[, "bowl_1"], (s$coords[, "x"] - mean(range(s$coords[, "x"])))^2), 1)
})

test_that("a gradient plus a field recovers direction, strength and the field's range", {
  skip_on_cran()
  # mean of 50 replicates: direction within 5 degrees, strength within 15% and the
  # variogram range of the detrended component within 25% (protocol, section 7)
  g <- sim_groups(sim_group("grad", 1, "gradient", direction = 60, strength = 1, range = 0.25,
                            large_scale_share = 1))
  est <- t(vapply(1:50, function(i) {
    s <- sim_spatial(400, g, seed = 300 + i)
    d <- data.frame(s$coords, z = s$components$large[, 1])
    fit <- stats::lm(z ~ x + y, d)
    b <- stats::coef(fit)[c("x", "y")]
    d$r <- stats::residuals(fit)
    v <- gstat::variogram(r ~ 1, locations = ~ x + y, data = d, cutoff = 0.5)
    vf <- suppressWarnings(gstat::fit.variogram(v, gstat::vgm(1, "Exp", 0.25 / 3)))
    c(atan2(b[2], b[1]) * 180 / pi, stats::sd(stats::fitted(fit)) / stats::sd(d$r), 3 * vf$range[nrow(vf)])
  }, numeric(3)))
  m <- colMeans(est)
  expect_lt(abs(m[1] - 60), 5)
  expect_lt(abs(m[2] - 1), 0.15)
  expect_lt(abs(m[3] / 0.25 - 1), 0.25)
})
