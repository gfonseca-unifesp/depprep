test_that("the checks stop on infinite values, on incomplete positions and on columns with too few values", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(80, g, seed = 1)
  X <- s$X
  X[3, 2] <- Inf
  expect_error(diagnose_variables(X, coords = s$coords, nsim = 9), "only `X` may have missing values")
  co <- s$coords
  co[1, 1] <- NaN
  expect_error(diagnose_variables(s$X, coords = co, nsim = 9), "`coords`")
  expect_error(diagnose_variables(s$X[, 1:2], time = c(NA, seq_len(79))), "`time`")
  X <- s$X
  X[-(1:2), 4] <- NA
  expect_error(diagnose_variables(X, coords = s$coords, nsim = 9), "Fewer than 3 observed values in none_1")
  expect_true(depprep:::.check_values(s$X, s$coords))
  X[1, 1] <- NA
  expect_error(depprep:::.check_values(X, s$coords, missing = FALSE), "not supported here")
})

test_that("a spatial variable with missing values is diagnosed at its observed sites", {
  g <- sim_groups(sim_group("a", 2, "matern", range = 0.4), sim_group("none", 1, "none"))
  s <- sim_spatial(150, g, seed = 2)
  X <- s$X
  gone <- c(3, 17, 40:52, 99)
  X[gone, 2] <- NA
  full <- diagnose_variables(X, coords = s$coords, nsim = 19, perm_seed = 1)
  obs <- setdiff(seq_len(nrow(X)), gone)
  # the same statistics as on the observed sites alone (same permutation stream: seed + column)
  alone <- diagnose_variables(X[obs, 2, drop = FALSE], coords = s$coords[obs, ], nsim = 19, perm_seed = 2)
  cols <- c("moran_I", "moran_p", "correlogram_p", "trend_r2", "range", "nugget_ratio", "aniso_ratio")
  expect_equal(unlist(full[2, cols]), unlist(alone[1, cols]))
  # the complete variables are not affected by the gaps of another one
  ref <- diagnose_variables(s$X, coords = s$coords, nsim = 19, perm_seed = 1)
  expect_equal(full[c(1, 3), cols], ref[c(1, 3), cols], ignore_attr = TRUE)
})

test_that("a series with gaps keeps its persistence and its cycle", {
  g <- sim_groups(sim_group("s", 2, "seasonal", period = 12))
  s <- sim_temporal(240, g, seed = 3)
  X <- s$X
  set.seed(4)
  X[sample.int(240, 24), 1] <- NA                         # 10% at random
  X[100:111, 2] <- NA                                     # a whole year
  sig <- diagnose_variables(X, time = s$time, nsim = 0)
  ref <- diagnose_variables(s$X, time = s$time, nsim = 0)
  expect_true(all(is.finite(as.matrix(sig[, c("phi", "acf1", "period", "period_strength", "trend_r2", "ljung_p")]))))
  expect_equal(sig$period, ref$period, tolerance = 0.05)
  expect_equal(sig$phi, ref$phi, tolerance = 0.25)
  expect_true(all(sig$period_p < 0.05))
  # without gaps nothing changes
  expect_equal(depprep:::.temporal_stats(s$X[, 1]), depprep:::.temporal_stats(replace(s$X[, 1], integer(0), NA)))
})

test_that("a spatio-temporal variable with gaps gets a full signature", {
  g <- sim_groups(sim_group("tc", 2, "temporal_common", phi = 0.7), sim_group("none", 1, "none"))
  s <- sim_spatiotemporal(60, 40, g, seed = 5)
  X <- s$X
  set.seed(6)
  X[sample.int(nrow(X), round(0.15 * nrow(X))), 1] <- NA  # 15% at random
  X[s$index$time == 7, 2] <- NA                           # one whole time missing
  sig <- diagnose_variables(X, coords = s$coords, index = s$index, nsim = 0)
  ref <- diagnose_variables(s$X, coords = s$coords, index = s$index, nsim = 0)
  cols <- c("moran_I", "moran_p", "phi", "acf1", "ljung_p", "shift_indicator")
  expect_true(all(is.finite(as.matrix(sig[, cols]))))
  expect_equal(sig$phi[1:2], ref$phi[1:2], tolerance = 0.2)
  expect_identical(classify_regime(sig), classify_regime(ref))
})

test_that("the basis of a variable with gaps is the basis with its gaps at its mean", {
  g <- sim_groups(sim_group("a", 1, "matern", range = 0.4))
  s <- sim_spatial(200, g, seed = 7)
  x <- s$X
  x[c(5, 30:45, 120), 1] <- NA
  b_na <- build_group_basis(x, coords = s$coords, type = "isotropic", seed = 1)
  filled <- x
  filled[is.na(filled)] <- mean(x, na.rm = TRUE)
  # filling with the mean only rescales the standardised column: the same tests, the same terms
  b_mean <- build_group_basis(filled, coords = s$coords, type = "isotropic", seed = 1)
  expect_gt(attr(b_na, "k"), 0)
  expect_identical(attr(b_na, "terms"), attr(b_mean, "terms"))
  expect_equal(nrow(b_na), 200)
})

test_that("a spatio-temporal group with gaps selects nearly the terms of the complete group", {
  g <- sim_groups(sim_group("tc", 3, "temporal_common", phi = 0.9))
  s <- sim_spatiotemporal(40, 30, g, seed = 8)
  X <- s$X
  set.seed(9)
  X[sample.int(length(X), round(0.1 * length(X)))] <- NA
  b <- build_group_basis(X, coords = s$coords, index = s$index, type = "temporal_common", seed = 1)
  ref <- build_group_basis(s$X, coords = s$coords, index = s$index, type = "temporal_common", seed = 1)
  expect_gt(attr(b, "k"), 0)
  expect_gt(length(intersect(attr(b, "terms"), attr(ref, "terms"))) / attr(ref, "k"), 0.7)
})

test_that("correct() projects each variable on its basis over its observed rows", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.4), sim_group("none", 2, "none"))
  s <- sim_spatial(200, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 10)
  X <- s$X
  set.seed(11)
  gone <- sample.int(200, 20)
  X[gone, 1] <- NA
  X[1:15, 2] <- NA
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(X))
  cr <- correct(X, grp, coords = s$coords, scope = "groups", seed = 1)
  b <- unclass(cr$bases$g1)
  o <- !is.na(X[, 1])
  expect_equal(unname(cr$X[o, 1]), unname(stats::lm.fit(cbind(1, unclass(cr$variable_bases[[1]])[o, ]), X[o, 1])$residuals))
  expect_true(all(is.na(cr$X[!o, 1])))
  # R pairwise: both residuals on the group basis, on the rows where both have a value
  e1 <- stats::lm.fit(cbind(1, b[o, ]), X[o, 1])$residuals
  o2 <- !is.na(X[, 2])
  e2 <- stats::lm.fit(cbind(1, b[o2, ]), X[o2, 2])$residuals
  both <- o & o2
  expect_equal(cr$R_pairwise[1, 2], stats::cor(e1[both[o]], e2[both[o2]]))
  # small gaps change R little
  ref <- correct(s$X, grp, coords = s$coords, scope = "groups", seed = 1)
  expect_lt(max(abs(cr$R - ref$R)), 0.1)
})

test_that("the whitening mode accepts missing values in every family", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.4), sim_group("none", 1, "none"))
  s <- sim_spatial(150, g, seed = 12)
  X <- s$X; X[c(2, 30:40), 1] <- NA
  w <- correct(X, stats::setNames(c("g1", "g1", "g1", "none"), colnames(X)), coords = s$coords, mode = "whiten", seed = 1)
  expect_true(all(is.finite(w$R)))
  expect_true(all(is.na(w$X[c(2, 30:40), 1])))

  gt <- sim_groups(sim_group("ar", 3, "ar1", phi = 0.8))
  st <- sim_temporal(150, gt, seed = 13)
  X <- st$X; X[50:55, 2] <- NA
  wt <- correct(X, group_regimes(diagnose_variables(X, time = st$time, nsim = 0)), time = st$time, mode = "whiten", seed = 1)
  expect_true(all(is.finite(wt$R)))
  expect_gt(sum(is.na(wt$X[, 2])), 5)                     # the gap and the lags that need it

  gs <- sim_groups(sim_group("tc", 3, "temporal_common", phi = 0.8), sim_group("none", 1, "none"))
  ss <- sim_spatiotemporal(40, 25, gs, seed = 14)
  X <- ss$X
  set.seed(15)
  X[sample.int(nrow(X), 100), 1] <- NA
  X[ss$index$time == 4, 2] <- NA
  gr <- group_regimes(diagnose_variables(ss$X, coords = ss$coords, index = ss$index, nsim = 0), method = "type")
  ws <- correct(X, gr, coords = ss$coords, index = ss$index, seed = 1)
  expect_identical(ws$mode, "whiten")
  expect_true(all(is.finite(ws$R)))
  ref <- correct(ss$X, gr, coords = ss$coords, index = ss$index, seed = 1)
  expect_lt(max(abs(ws$R - ref$R)), 0.15)
})

test_that("correct() stops on missing values in the scale mode and on pairs never observed together", {
  g <- sim_groups(sim_group("a", 2, "matern", range = 0.4), sim_group("none", 2, "none"))
  s <- sim_spatial(100, g, seed = 16)
  X <- s$X
  X[1:50, 3] <- NA
  X[51:100, 4] <- NA
  grp <- stats::setNames(c("g1", "g1", "none", "none"), colnames(X))
  expect_error(correct(X, grp, coords = s$coords, mode = "scale"), "scale mode needs complete data")
  expect_error(correct(X, grp, coords = s$coords, scope = "groups"), "none_1 and none_2")
})

test_that("verify() describes every family with missing values on the observed rows", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.4), sim_group("none", 1, "none"))
  s <- sim_spatial(150, g, seed = 17)
  X <- s$X
  X[c(4, 60:75), 2] <- NA
  grp <- stats::setNames(c("g1", "g1", "g1", "none"), colnames(X))
  cr <- correct(X, grp, coords = s$coords, scope = "groups", seed = 1)
  v <- verify(cr, X, coords = s$coords, nperm = 19)
  expect_true(all(is.finite(as.matrix(v$variables[, c("removed", "autocor_before", "autocor_after", "correlogram_p")]))))
  expect_false(any(v$lisa$site[v$lisa$variable == colnames(X)[2]] %in% c(4, 60:75)))
  expect_true(is.finite(v$q3$statistic))
  expect_false(is.null(v$mantel))
  ref <- verify(correct(s$X, grp, coords = s$coords, scope = "groups", seed = 1), s$X, coords = s$coords, nperm = 19)
  expect_equal(v$q3$statistic, ref$q3$statistic, tolerance = 0.2)

  gt <- sim_groups(sim_group("ar", 3, "ar1", phi = 0.8))
  st <- sim_temporal(150, gt, seed = 18)
  X <- st$X; X[c(10, 80:90), 1] <- NA
  grt <- group_regimes(diagnose_variables(X, time = st$time, nsim = 0))
  vt <- verify(correct(X, grt, time = st$time, seed = 1), X, time = st$time, nperm = 19)
  expect_true(all(is.finite(as.matrix(vt$variables[, c("removed", "autocor_after", "lb_p", "lb_p_before")]))))

  gs <- sim_groups(sim_group("tc", 3, "temporal_common", phi = 0.8), sim_group("none", 1, "none"))
  ss <- sim_spatiotemporal(40, 25, gs, seed = 19)
  X <- ss$X
  set.seed(20)
  X[sample.int(nrow(X), 100), 1] <- NA
  X[ss$index$site == 3, 2] <- NA                          # a site never observed in one variable
  grs <- group_regimes(diagnose_variables(ss$X, coords = ss$coords, index = ss$index, nsim = 0), method = "type")
  for (m in c("conditional", "whiten")) {
    vs <- verify(correct(X, grs, coords = ss$coords, index = ss$index, mode = m, seed = 1), X,
                 coords = ss$coords, index = ss$index, nperm = 19)
    expect_true(all(is.finite(as.matrix(vs$variables[, c("removed", "autocor_before", "autocor_after", "correlogram_p")]))))
    expect_false(3 %in% vs$lisa$site[vs$lisa$variable == colnames(X)[2]])
  }
})
