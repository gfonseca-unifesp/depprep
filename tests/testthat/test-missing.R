test_that("missing values stop the diagnosis, the bases and the correction with a clear message", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(80, g, seed = 1)
  X <- s$X
  X[3, 2] <- NA
  expect_error(diagnose_variables(X, coords = s$coords, nsim = 9), "missing values are not supported yet")
  expect_error(build_group_basis(X[, 1:3], coords = s$coords, type = "isotropic"), "`X`")
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  expect_error(correct(X, grp, coords = s$coords), "missing values are not supported yet")
  co <- s$coords
  co[1, 1] <- NaN
  expect_error(correct(s$X, grp, coords = co), "`coords`")
  expect_error(diagnose_variables(s$X[, 1:2], time = c(NA, seq_len(79))), "`time`")
  expect_true(depprep:::.stop_if_missing(s$X, s$coords))
})
