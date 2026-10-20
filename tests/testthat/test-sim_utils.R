test_that(".with_stream is reproducible and restores the caller's RNG", {
  a <- depprep:::.with_stream(1, 2, stats::rnorm(3))
  b <- depprep:::.with_stream(1, 2, stats::rnorm(3))
  c <- depprep:::.with_stream(1, 3, stats::rnorm(3))
  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, c)))
  set.seed(42); u1 <- stats::runif(1)
  set.seed(42); depprep:::.with_stream(7, 1, stats::rnorm(10)); u2 <- stats::runif(1)
  expect_identical(u1, u2)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
})

test_that(".aniso_coords stretches the minor axis only", {
  xy <- rbind(c(0, 0), c(1, 0), c(0, 1))
  out <- depprep:::.aniso_coords(xy, ratio = 2, angle = 0)
  expect_equal(out[2, ], c(1, 0))
  expect_equal(out[3, ], c(0, 2))
  out90 <- depprep:::.aniso_coords(xy, ratio = 2, angle = 90)
  expect_equal(abs(out90[3, ]), c(1, 0), tolerance = 1e-12)
})

test_that(".grf_chol has the exponential correlation", {
  xy <- rbind(c(0, 0), c(0.1, 0))
  z <- depprep:::.with_stream(3, 1, depprep:::.grf_chol(xy, range = 0.3, n_sim = 4000))
  expect_equal(stats::cor(z[1, ], z[2, ]), exp(-1), tolerance = 0.05)
  expect_equal(mean(apply(z, 1, stats::var)), 1, tolerance = 0.05)
})

test_that(".rff_sampler approximates the exponential correlation", {
  xy <- rbind(c(0.2, 0.2), c(0.3, 0.2))
  z <- depprep:::.with_stream(4, 1, replicate(3000, depprep:::.rff_sampler(range = 0.3, n_features = 500)(xy)))
  expect_equal(stats::cor(z[1, ], z[2, ]), exp(-1), tolerance = 0.06)
  expect_equal(mean(apply(z, 1, stats::var)), 1, tolerance = 0.06)
})

test_that(".ar1, .seasonal and .trend_shift behave as specified", {
  x <- depprep:::.with_stream(5, 1, depprep:::.ar1(2000, phi = 0.7, n_series = 3))
  expect_equal(dim(x), c(2000L, 3L))
  expect_equal(mean(apply(x, 2, function(s) stats::acf(s, 1, plot = FALSE)$acf[2])), 0.7, tolerance = 0.05)
  s <- depprep:::.seasonal(120, period = 12)
  expect_equal(stats::var(s), 1, tolerance = 0.02)
  expect_equal(s[1:12], s[13:24])
  tr <- depprep:::.trend_shift(100, slope = 0, jump = 1, shift_at = 0.5)
  expect_equal(which.max(diff(tr)), 50L)
  expect_equal(mean(tr), 0, tolerance = 1e-12)
})

test_that(".residual_mvn has correlation R and rejects invalid matrices", {
  r <- matrix(c(1, 0.6, 0.6, 1), 2)
  e <- depprep:::.with_stream(6, 1, depprep:::.residual_mvn(20000, r))
  expect_equal(stats::cor(e)[1, 2], 0.6, tolerance = 0.02)
  expect_error(depprep:::.residual_mvn(10, matrix(c(1, 2, 2, 1), 2)), "positive definite")
})

test_that("contiguous classes are spatially clustered and random ones are not", {
  xy <- depprep:::.with_stream(8, 1, matrix(stats::runif(800), ncol = 2))
  nn <- apply(as.matrix(stats::dist(xy)) + diag(Inf, 400), 1, which.min)
  cont <- depprep:::.with_stream(8, 2, depprep:::.contiguous_classes(3, coords = xy))
  rand <- depprep:::.with_stream(8, 2, depprep:::.random_classes(3, 400))
  expect_gt(mean(cont == cont[nn]), 0.8)
  expect_lt(mean(rand == rand[nn]), 0.5)
  expect_equal(depprep:::.contiguous_classes(3, n_times = 9), rep(1:3, each = 3))
})

test_that("sim_groups validates groups and names variables", {
  g <- sim_groups(sim_group("long", 2, "exponential", range = 0.5),
                  sim_group("copy", 1, "exponential", shared_with = "long"),
                  sim_group("none", 2, "none", large_scale_share = 0.9))
  expect_equal(attr(g, "var_group"), c("long", "long", "copy", "none", "none"))
  expect_equal(attr(g, "var_names")[1:3], c("long_1", "long_2", "copy_1"))
  expect_equal(g$none$large_scale_share, 0)
  expect_error(sim_groups(sim_group("a", 1, "none"), sim_group("a", 1, "none")), "unique")
  expect_error(sim_groups(sim_group("a", 1, "ar1", shared_with = "b")), "unknown group")
  expect_error(sim_groups(sim_group("n", 1, "none"), sim_group("a", 1, "ar1", shared_with = "n")), "without dependence")
})

test_that("make_R_true builds valid identity and block matrices", {
  g <- sim_groups(sim_group("a", 2, "none"), sim_group("b", 3, "none"))
  expect_equal(unname(make_R_true(g)), diag(5))
  r <- make_R_true(g, "blocks", within = 0.5, between = 0.2)
  expect_equal(r["a_1", "a_2"], 0.5)
  expect_equal(r["a_1", "b_1"], 0.2)
  expect_true(isSymmetric(r))
  expect_error(make_R_true(g, "blocks", within = -0.9, between = 0.9), "positive definite")
})
