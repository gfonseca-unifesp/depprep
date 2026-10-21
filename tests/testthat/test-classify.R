sig <- function(family, ...) {
  row <- depprep:::.empty_signature()
  row[names(list(...))] <- list(...)
  row$variable <- "v"
  d <- as.data.frame(row, stringsAsFactors = FALSE)
  d$family <- family
  d
}

test_that("classify_regime applies the rules of each family", {
  expect_equal(unname(classify_regime(sig("spatial", moran_p = 0.4))), "none")
  expect_equal(unname(classify_regime(sig("spatial", moran_p = 0.001, aniso_ratio = 1.3))), "isotropic")
  expect_equal(unname(classify_regime(sig("spatial", moran_p = 0.001, aniso_ratio = 3.5))), "anisotropic")
  expect_equal(unname(classify_regime(sig("temporal", ljung_p = 0.3))), "none")
  expect_equal(unname(classify_regime(sig("temporal", ljung_p = 0.001, period_strength = 0.4, trend_r2 = 0))), "seasonal")
  expect_equal(unname(classify_regime(sig("temporal", ljung_p = 0.001, period_strength = 0.02, trend_r2 = 0.7))), "trend")
  expect_equal(unname(classify_regime(sig("temporal", ljung_p = 0.001, period_strength = 0.02, trend_r2 = 0.05))), "ar")
  expect_equal(unname(classify_regime(sig("spatiotemporal", moran_p = 0.6, ljung_p = 0.6))), "none")
  expect_equal(unname(classify_regime(sig("spatiotemporal", moran_p = 0.6, ljung_p = 0.001))), "temporal_common")
  expect_equal(unname(classify_regime(sig("spatiotemporal", moran_p = 0.001, ljung_p = 0.001, shift_indicator = 0.2))), "moving")
  expect_equal(unname(classify_regime(sig("spatiotemporal", moran_p = 0.001, ljung_p = 0.001, shift_indicator = 0.01,
                                          period_strength = 0.5))), "separable_seasonal")
  expect_equal(unname(classify_regime(sig("spatiotemporal", moran_p = 0.001, ljung_p = 0.001, shift_indicator = 0.01,
                                          period_strength = 0.01))), "separable_ar1")
})

test_that("thresholds change the labels and default_thresholds is complete", {
  s <- sig("spatial", moran_p = 0.03, aniso_ratio = 1.5)
  expect_equal(unname(classify_regime(s)), "isotropic")
  expect_equal(unname(classify_regime(s, list(p_dependence = 0.01))), "none")
  expect_equal(unname(classify_regime(s, list(aniso_ratio = 1.2))), "anisotropic")
  expect_named(default_thresholds(), c("p_dependence", "aniso_ratio", "seasonal_strength", "trend_r2", "shift"))
})

test_that("clean simulated temporal regimes are classified correctly", {
  g <- sim_groups(sim_group("ar", 3, "ar1", phi = 0.6), sim_group("seasonal", 3, "seasonal", period = 12),
                  sim_group("trend", 3, "trend_shift"), sim_group("none", 3, "none"))
  s <- sim_temporal(480, g, seed = 11)
  lab <- classify_regime(diagnose_variables(s))
  expect_gte(mean(lab == rep(c("ar", "seasonal", "trend", "none"), each = 3)), 0.9)
})
