test_that("verify reports removed variance, autocorrelation and Delta R", {
  g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(300, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 1)
  grp <- stats::setNames(c(rep("g1", 4), rep("none", 2)), colnames(s$X))
  b <- list(g1 = build_group_basis(s$X[, 1:4], coords = s$coords, type = "isotropic", seed = 1))
  v <- verify(correct(s, grp, bases = b), s)
  expect_s3_class(v, "depprep_verification")
  long <- v$variables$group == "g1"
  expect_true(all(v$variables$removed[long] > 0.3))       # the shared field is mostly removed
  expect_true(all(v$variables$autocor_after[long] < v$variables$autocor_before[long]))
  expect_setequal(v$delta_summary$kind, c("within groups", "involving no dependence"))
  # the apparent association inside the group falls
  expect_lt(v$delta_summary$mean_delta[v$delta_summary$kind == "within groups"], 0)
  expect_equal(v$delta_R, t(v$delta_R))
})
