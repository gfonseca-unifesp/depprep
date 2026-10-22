is_orthonormal <- function(b) isTRUE(all.equal(crossprod(unclass(b)), diag(ncol(b)), tolerance = 1e-8,
                                                 check.attributes = FALSE))

test_that("spatial basis is orthonormal, bounded and captures a long-range field", {
  g <- sim_groups(sim_group("long", 6, "exponential", range = 0.5), sim_group("none", 4, "none"))
  s <- sim_spatial(400, g, seed = 11)
  long <- s$group == "long"
  b <- build_group_basis(s$X[, long], coords = s$coords, type = "isotropic", seed = 1)
  expect_s3_class(b, "depprep_basis")
  expect_true(is_orthonormal(b))
  expect_gt(attr(b, "k"), 0)
  expect_lte(attr(b, "kn"), 0.25)
  ca <- depprep:::.capture_absorption(b, s$components$large[, long], s$components$residual[, long])
  expect_gt(ca[["capture"]], 0.7)
  expect_lt(ca[["absorption"]], 0.1)
})

test_that("one basis per group: reproducible and independent of the order of the variables", {
  g <- sim_groups(sim_group("long", 5, "exponential", range = 0.5))
  s <- sim_spatial(200, g, seed = 3)
  b1 <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 7)
  b2 <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 7)
  b3 <- build_group_basis(s$X[, 5:1], coords = s$coords, type = "isotropic", seed = 7)
  expect_identical(b1, b2)
  # same subspace (near ties between candidates may swap the order of the columns)
  proj <- function(b) unclass(b) %*% t(unclass(b))
  expect_equal(proj(b1), proj(b3), tolerance = 1e-8)
  # every variable of the group is projected on the very same matrix
  P <- unclass(b1) %*% t(unclass(b1))
  one_by_one <- vapply(seq_len(ncol(s$X)), function(j) P %*% s$X[, j], numeric(nrow(s$X)))
  expect_equal(one_by_one, P %*% s$X, ignore_attr = TRUE)
})

test_that("the none group has an empty basis and kn_max bounds k", {
  g <- sim_groups(sim_group("long", 4, "exponential", range = 0.5))
  s <- sim_spatial(200, g, seed = 5)
  expect_equal(ncol(build_group_basis(s$X, coords = s$coords, type = "none")), 0)
  b <- build_group_basis(s$X, coords = s$coords, type = "isotropic", kn_max = 0.02, seed = 1)
  expect_lte(attr(b, "k"), 4)
  expect_error(build_group_basis(s$X, coords = s$coords, type = "wave"), "Unknown")
})

test_that("temporal bases: harmonics for seasonal groups, slow terms for AR groups", {
  g <- sim_groups(sim_group("seas", 5, "seasonal", period = 12), sim_group("ar", 5, "ar1", phi = 0.9))
  s <- sim_temporal(240, g, seed = 2)
  seas <- s$group == "seas"
  b <- build_group_basis(s$X[, seas], time = s$time, type = "seasonal", period = 12, seed = 1)
  expect_true(all(grepl("_p12$|^trend", attr(b, "terms"))))
  expect_true(any(attr(b, "terms") %in% c("sin1_p12", "cos1_p12")))
  expect_true(is_orthonormal(b))
  expect_error(build_group_basis(s$X[, seas], time = s$time, type = "seasonal"), "period")
  a <- build_group_basis(s$X[, !seas], time = s$time, type = "ar", seed = 1)
  expect_gt(attr(a, "k"), 0)
  expect_true(all(grepl("^tMEM", attr(a, "terms"))))
  # rows in a different order give the same basis for the same rows
  o <- sample.int(240)
  a2 <- build_group_basis(s$X[o, !seas], time = s$time[o], type = "ar", seed = 1)
  expect_equal(unclass(a2) %*% t(unclass(a2)), (unclass(a) %*% t(unclass(a)))[o, o], tolerance = 1e-8)
})

test_that("spatio-temporal bases are products of spatial and temporal terms", {
  g <- sim_groups(sim_group("sep", 4, "separable_ar1", range = 0.5, phi = 0.9),
                  sim_group("common", 4, "temporal_common", phi = 0.9))
  s <- sim_spatiotemporal(50, 60, g, seed = 4)
  sep <- s$group == "sep"
  b <- build_group_basis(s$X[, sep], coords = s$coords, index = s$index, type = "separable_ar1", seed = 1)
  expect_equal(nrow(b), 50 * 60)
  expect_true(is_orthonormal(b))
  expect_true(all(grepl(":", attr(b, "terms"))))
  ca <- depprep:::.capture_absorption(b, s$components$large[, sep], s$components$residual[, sep])
  # an AR(1) field with phi = 0.9 keeps much variance at short lags, out of reach of a
  # slow basis: all 449 products of the candidates capture about 0.57 here
  expect_gt(ca[["capture"]], 0.35)
  expect_lt(ca[["absorption"]], 0.05)
  cm <- build_group_basis(s$X[, !sep], coords = s$coords, index = s$index, type = "temporal_common", seed = 1)
  expect_true(all(startsWith(attr(cm, "terms"), "const:")))
  expect_error(build_group_basis(s$X[-1, sep], coords = s$coords, index = s$index[-1, ], type = "moving"),
               "complete grid")
})

test_that("the union basis contains each basis and drops redundant directions", {
  g <- sim_groups(sim_group("long", 3, "exponential", range = 0.5), sim_group("short", 3, "exponential", range = 0.1))
  s <- sim_spatial(300, g, seed = 8)
  b1 <- build_group_basis(s$X[, 1:3], coords = s$coords, type = "isotropic", seed = 1)
  b2 <- build_group_basis(s$X[, 4:6], coords = s$coords, type = "isotropic", seed = 1)
  u <- union_basis(b1, b2)
  expect_true(is_orthonormal(u))
  expect_equal(attr(u, "k"), qr(cbind(unclass(b1), unclass(b2)))$rank)
  inside <- function(b) max(abs(unclass(b) - unclass(u) %*% crossprod(unclass(u), unclass(b))))
  expect_lt(inside(b1), 1e-8)
  expect_lt(inside(b2), 1e-8)
  expect_equal(attr(union_basis(b1, b1), "k"), attr(b1, "k"))
  expect_warning(union_basis(b1, b2, kn_max = 0.001), "kn_max")
})

test_that("trend terms come first among the candidates and keep their shape", {
  set.seed(1)
  xy <- matrix(stats::runif(400), 200)
  cand <- depprep:::.spatial_candidates(xy)
  expect_identical(colnames(cand)[1:5], c("trend_x", "trend_y", "trend_x2", "trend_xy", "trend_y2"))
  expect_equal(abs(stats::cor(cand[, 1], xy[, 1])), 1)
  expect_true(startsWith(colnames(cand)[6], "MEM"))
  expect_false(any(grepl("^trend", colnames(depprep:::.spatial_candidates(xy, trend = FALSE)))))
  tc <- depprep:::.temporal_candidates(1:120, "seasonal", 12)
  expect_identical(colnames(tc)[1:4], c("trend1", "trend2", "trend3", "sin1_p12"))
  expect_false(any(grepl("^trend", colnames(depprep:::.temporal_candidates(1:120, "ar")))))
})

test_that("a spatial gradient plus a field is captured with fewer terms than without trend terms", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.2))
  s <- sim_spatial(300, g, seed = 1)
  grad <- 2 * as.numeric(scale(s$coords[, 1] + 0.5 * s$coords[, 2]))
  Y <- scale(s$X + grad)
  set.seed(1); with_trend <- depprep:::.forward_select(Y, depprep:::.spatial_candidates(s$coords), 0.25, 0.05, 199)
  set.seed(1); without <- depprep:::.forward_select(Y, depprep:::.spatial_candidates(s$coords, trend = FALSE), 0.25, 0.05, 199)
  expect_true(all(c("trend_x", "trend_y") %in% with_trend$terms))
  expect_lt(length(with_trend$terms), length(without$terms))
  capture <- function(B) sum(crossprod(B, grad - mean(grad))^2) / sum((grad - mean(grad))^2)
  expect_gt(capture(with_trend$basis), 0.999)
  expect_gt(capture(with_trend$basis), capture(without$basis))
  # through build_group_basis() as well
  b <- build_group_basis(s$X + grad, coords = s$coords, type = "isotropic", seed = 1)
  expect_true("trend_x" %in% attr(b, "terms"))
})

test_that("a seasonal series with a trend gets a trend term", {
  set.seed(2)
  t <- 1:240
  x <- cbind(sin(2 * pi * t / 12) + 0.006 * t + stats::rnorm(240, sd = 0.7),
             cos(2 * pi * t / 12) + 0.006 * t + stats::rnorm(240, sd = 0.7))
  b <- build_group_basis(x, time = t, type = "seasonal", period = 12, seed = 1)
  expect_true("trend1" %in% attr(b, "terms"))
  expect_true(all(c("sin1_p12", "cos1_p12") %in% attr(b, "terms")))
  expect_true(is_orthonormal(b))
})

test_that("the basis keeps the scalogram of every candidate", {
  g <- sim_groups(sim_group("long", 4, "exponential", range = 0.5))
  s <- sim_spatial(200, g, seed = 3)
  b <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 1)
  sg <- attr(b, "scalogram")
  expect_identical(names(sg), colnames(depprep:::.spatial_candidates(s$coords)))
  expect_true(all(attr(b, "terms") %in% names(sg)))
  expect_true(all(sg >= 0) && sum(sg) <= 1 + 1e-8)
  expect_equal(sum(sg[attr(b, "terms")]), attr(b, "r2"))
  expect_length(attr(build_group_basis(s$X, coords = s$coords, type = "none"), "scalogram"), 0)
})

test_that("capture and absorption measure what a basis removes", {
  n <- 100
  B <- qr.Q(qr(cbind(scale(1:n), scale((1:n)^2))))
  large <- B %*% matrix(c(1, 2, 3, 4), 2)
  res <- matrix(stats::rnorm(n * 2), n)
  ca <- depprep:::.capture_absorption(structure(B, class = "depprep_basis"), large, res)
  expect_equal(ca[["capture"]], 1)
  expect_lt(abs(ca[["absorption"]]), 0.1)
})
