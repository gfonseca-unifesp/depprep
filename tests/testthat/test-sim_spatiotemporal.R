groups_st <- function() {
  sim_groups(sim_group("sep_long", 3, "separable_ar1", range = 0.5, phi = 0.9),
             sim_group("sep_short", 3, "separable_seasonal", range = 0.1, period = 12),
             sim_group("move", 2, "moving", range = 0.25, velocity = c(0.02, 0)),
             sim_group("common", 2, "temporal_common", phi = 0.6),
             sim_group("none", 2, "none"))
}

field <- function(s, var) matrix(s$components$large[, var], nrow(s$coords))  # sites x times

test_that("sim_spatiotemporal returns the documented structure and order", {
  s <- sim_spatiotemporal(20, 15, groups_st(), seed = 1)
  expect_equal(dim(s$X), c(300L, 12L))
  expect_equal(nrow(s$index), 300L)
  expect_equal(s$index[1:3, "site"], 1:3)
  expect_equal(s$index[21, ], data.frame(site = 1L, time = 2L), ignore_attr = TRUE)
  rebuilt <- sweep(s$components$large, 2, sqrt(s$large_scale_share), "*") +
    sweep(s$components$residual, 2, sqrt(1 - s$large_scale_share), "*")
  expect_equal(s$X, rebuilt)
  expect_identical(s$X, sim_spatiotemporal(20, 15, groups_st(), seed = 1)$X)
})

test_that("separable AR(1) fields have temporal phi and the spatial exponential correlation", {
  xy <- rbind(c(0.2, 0.2), c(0.3, 0.2))
  g <- sim_groups(sim_group("f", 1, "separable_ar1", range = 0.3, phi = 0.5, large_scale_share = 1))
  s <- sim_spatiotemporal(groups = g, n_times = 6000, coords = xy, seed = 2)
  z <- field(s, "f_1")
  expect_equal(stats::acf(z[1, ], 1, plot = FALSE)$acf[2], 0.5, tolerance = 0.05)
  expect_lt(abs(stats::cor(z[1, ], z[2, ]) - exp(-1)), 0.05)  # absolute tolerance
})

test_that("separable seasonal fields are an amplitude map times one cycle", {
  s <- sim_spatiotemporal(10, 48, groups_st(), seed = 3)
  z <- field(s, "sep_short_1")
  expect_equal(abs(stats::cor(z[1, ], z[2, ])), 1, tolerance = 1e-10)
  sp <- stats::spec.pgram(z[1, ], taper = 0, plot = FALSE)
  expect_equal(1 / sp$freq[which.max(sp$spec)], 12)
})

test_that("the moving pattern is the same field displaced by the velocity at each step", {
  # grid with spacing equal to the velocity, so s - v is also a grid point
  xy <- as.matrix(expand.grid(x = seq(0, 1, by = 0.02), y = c(0.3, 0.6)))
  g <- sim_groups(sim_group("m", 1, "moving", range = 0.25, velocity = c(0.02, 0), large_scale_share = 1))
  s <- sim_spatiotemporal(groups = g, n_times = 5, coords = xy, seed = 4)
  z <- field(s, "m_1")
  left <- match(paste(round(xy[, 1] - 0.02, 6), xy[, 2]), paste(round(xy[, 1], 6), xy[, 2]))
  ok <- !is.na(left)
  expect_equal(z[ok, 3], z[left[ok], 2], tolerance = 1e-8)
  # a frozen pattern would not decorrelate; the moving one does
  expect_lt(stats::cor(z[, 1], z[, 5]), stats::cor(z[, 1], z[, 2]))
})

test_that("the temporal-only regime is identical at every site", {
  s <- sim_spatiotemporal(8, 30, groups_st(), seed = 5)
  z <- field(s, "common_1")
  expect_true(all(apply(z, 2, function(col) diff(range(col))) < 1e-12))
})

test_that("contiguous classes keep their class at every time", {
  s <- sim_spatiotemporal(60, 10, groups_st(), classes = list(k = 3, type = "contiguous"), seed = 6)
  by_site <- tapply(s$classes, s$index$site, function(x) length(unique(x)))
  expect_true(all(by_site == 1))
})

test_that("invalid regimes and missing parameters fail clearly", {
  expect_error(sim_spatiotemporal(10, 5, sim_groups(sim_group("a", 1, "ar1", phi = 0.5)), seed = 1), "not available")
  expect_error(sim_spatiotemporal(10, 5, sim_groups(sim_group("m", 1, "moving", range = 0.2)), seed = 1), "velocity")
})
