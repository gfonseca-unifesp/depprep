test_that("the codependence spectrum sums to the correlation (Parseval)", {
  set.seed(1)
  coords <- cbind(stats::runif(80), stats::runif(80))
  x <- coords[, 1] + stats::rnorm(80)
  y <- coords[, 1] + stats::rnorm(80)
  cd <- codependence(x, y, coords = coords)
  expect_equal(nrow(cd), 79)
  expect_equal(utils::tail(cd$cumulative, 1), stats::cor(x, y), tolerance = 1e-8)
  expect_equal(sum(cd$r_x^2), 1, tolerance = 1e-8)
  U <- attr(cd, "basis")
  expect_equal(crossprod(U), diag(79), tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(all(diff(attr(U, "eigenvalue")) <= 1e-10))           # broad to fine
})

test_that("a shared broad pattern puts the codependence at broad scales", {
  set.seed(2)
  tt <- 1:120
  common <- sin(2 * pi * tt / 120)
  x <- common + stats::rnorm(120, sd = 0.3)
  y <- common + stats::rnorm(120, sd = 0.3)
  cd <- codependence(x, y, time = tt)
  expect_gt(sum(cd$codependence[1:5]), 0.8 * stats::cor(x, y))
  cd2 <- codependence(x, y, basis = attr(cd, "basis"))
  expect_equal(cd2$codependence, cd$codependence)
  # residual correlation outside the first 5 functions = renormalised partial sum
  B <- attr(cd, "basis")[, 1:5]
  rx <- x - mean(x) - B %*% crossprod(B, x - mean(x))
  ry <- y - mean(y) - B %*% crossprod(B, y - mean(y))
  out <- -(1:5)
  expect_equal(stats::cor(as.vector(rx), as.vector(ry)),
               sum(cd$codependence[out]) / sqrt(sum(cd$r_x[out]^2) * sum(cd$r_y[out]^2)), tolerance = 1e-8)
})
