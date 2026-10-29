g_sp <- function() sim_groups(sim_group("long", 2, "exponential", range = 0.5),
                              sim_group("short", 2, "exponential", range = 0.1),
                              sim_group("none", 2, "none"))

test_that("diagnose_variables returns one signature row per variable", {
  skip_on_cran()
  s <- sim_spatial(400, g_sp(), seed = 1)
  d <- diagnose_variables(s)
  expect_s3_class(d, "depprep_signature")
  expect_equal(nrow(d), 6L)
  expect_equal(d$variable, colnames(s$X))
  expect_true(all(d$family == "spatial"))
  expect_true(all(is.na(d$phi)))
  expect_true(all(d$moran_p[1:4] < 0.05))
  expect_gt(stats::median(d$range[1:2], na.rm = TRUE), stats::median(d$range[3:4], na.rm = TRUE))
})

test_that("a variable's signature does not depend on the other variables", {
  skip_on_cran()
  s <- sim_spatial(200, g_sp(), seed = 2)
  d1 <- diagnose_variables(s$X, coords = s$coords)
  x2 <- s$X
  x2[, -1] <- depprep:::.with_stream(9, 1, matrix(stats::rnorm(length(x2[, -1])), nrow(x2)))  # replace the others
  d2 <- diagnose_variables(x2, coords = s$coords)
  expect_equal(d1[1, ], d2[1, ], ignore_attr = TRUE)
  x3 <- s$X[, c(1, 6, 5, 4, 3, 2)]                                                               # reorder the others
  d3 <- diagnose_variables(x3, coords = s$coords)
  expect_equal(d1[1, -1], d3[1, -1], ignore_attr = TRUE)
})

test_that("temporal and spatio-temporal families are diagnosed", {
  gt <- sim_groups(sim_group("ar", 2, "ar1", phi = 0.9), sim_group("seas", 2, "seasonal", period = 12),
                   sim_group("none", 2, "none"))
  dt <- diagnose_variables(sim_temporal(240, gt, seed = 3))
  expect_true(all(dt$family == "temporal"))
  expect_true(all(dt$period[3:4] == 12))
  expect_true(all(is.na(dt$moran_I)))
  gst <- sim_groups(sim_group("move", 1, "moving", range = 0.25, velocity = c(0.02, 0)),
                    sim_group("none", 1, "none"))
  dst <- diagnose_variables(sim_spatiotemporal(60, 40, gst, seed = 4))
  expect_true(all(dst$family == "spatiotemporal"))
  expect_gt(dst$shift_indicator[1], dst$shift_indicator[2])
})

test_that("time order is respected and inputs are checked", {
  gt <- sim_groups(sim_group("ar", 1, "ar1", phi = 0.9))
  s <- sim_temporal(200, gt, seed = 5)
  shuffled <- sample(200)
  expect_equal(diagnose_variables(s$X[shuffled, , drop = FALSE], time = s$time[shuffled])$phi,
               diagnose_variables(s)$phi)
  expect_error(diagnose_variables(s$X), "Give `coords`")
  expect_error(diagnose_variables(s$X, time = s$time, B = 10), "seed")
})

test_that("the bootstrap is reproducible and gives a positive standard error", {
  s <- sim_spatial(120, sim_groups(sim_group("f", 1, "exponential", range = 0.25)), seed = 6)
  a <- diagnose_variables(s, B = 15, seed = 1)
  b <- diagnose_variables(s, B = 15, seed = 1)
  expect_equal(a$range_se, b$range_se)
  expect_gt(a$range_se, 0)
  expect_lt(a$range_lo, a$range_hi)
  st <- sim_temporal(200, sim_groups(sim_group("f", 1, "ar1", phi = 0.6)), seed = 7)
  ph <- diagnose_variables(st, B = 30, seed = 2)
  expect_gt(ph$phi_se, 0)
})

# The parametric bootstrap slightly underestimates the spread of the AR estimate in
# series of this length (about 80-85% actual coverage for a nominal 90%); the test
# guards against gross errors, and the coverage is reported by the stage 3 analysis.
test_that("the 90% bootstrap interval of phi has roughly nominal coverage", {
  skip_on_cran()
  g <- sim_groups(sim_group("f", 1, "ar1", phi = 0.6, large_scale_share = 1))
  cover <- vapply(1:60, function(i) {
    d <- diagnose_variables(sim_temporal(240, g, seed = 700 + i), B = 99, seed = i)
    d$phi_lo <= 0.6 && 0.6 <= d$phi_hi
  }, logical(1))
  expect_gt(mean(cover), 0.75)
  expect_lt(mean(cover), 0.99)
})
