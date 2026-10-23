test_that("counts: a share of the variables becomes non-negative integers that follow the latent", {
  g <- sim_groups(sim_group("a", 6, "exponential", range = 0.3), sim_group("none", 3, "none"))
  s <- sim_spatial(300, g, seed = 3, counts = 1 / 3)
  expect_identical(s$count_vars, colnames(s$X)[c(3, 6, 9)])
  C <- s$X[, s$count_vars]
  expect_true(all(C >= 0 & C == round(C)))
  expect_lt(abs(mean(C) - 5), 0.6)                         # mean count about 5
  rho <- vapply(s$count_vars, function(v) stats::cor(s$X[, v], s$X_latent[, v], method = "spearman"), 0)
  expect_true(all(rho > 0.6))
  expect_equal(s$X[, -c(3, 6, 9)], s$X_latent[, -c(3, 6, 9)])   # the other variables are untouched
  t1 <- sim_temporal(200, sim_groups(sim_group("ar", 3, "ar1", phi = 0.5)), seed = 2, counts = 1 / 3)
  expect_identical(t1$count_vars, "ar_3")
})

test_that("without perturbations the data are unchanged", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.3), sim_group("none", 2, "none"))
  s0 <- sim_spatial(150, g, seed = 5)
  s1 <- sim_spatial(150, g, seed = 5, counts = 0)
  expect_identical(s0$X, s1$X)
  expect_null(s0$X_latent)
  expect_length(s0$count_vars, 0)
})

test_that("non-stationary range: the field is shorter on the left and longer on the right", {
  g <- sim_groups(sim_group("ns", 1, "exponential", range = 0.3, nonstationary = c(0.3, 2),
                            large_scale_share = 1))
  ratio <- vapply(1:6, function(i) {
    s <- sim_spatial(800, g, seed = 40 + i)
    z <- s$components$large[, 1]
    side <- function(keep) {                               # mean squared difference of near neighbours
      xy <- s$coords[keep, ]; zz <- z[keep]
      nn <- spdep::knearneigh(xy, 3)$nn
      mean((zz - zz[nn[, 1]])^2)
    }
    side(s$coords[, 1] < 0.3) / side(s$coords[, 1] > 0.7)
  }, 0)
  expect_gt(stats::median(ratio), 1.3)                    # neighbours differ more where the range is short
  s <- sim_spatial(400, g, seed = 1)
  expect_equal(stats::var(s$components$large[, 1]), 1, tolerance = 0.35)
})
