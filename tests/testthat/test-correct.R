exact_case <- function(rho, n = 400, seed = 1) {
  # large scale exactly in the span of 10 dbMEM, shared by two variables
  set.seed(seed)
  xy <- matrix(stats::runif(2 * n), n)
  B <- depprep:::.spatial_candidates(xy)[, 1:10]
  L <- as.numeric(B %*% stats::rnorm(10)); L <- L / stats::sd(L)
  e <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  X <- cbind(a = sqrt(0.6) * L + sqrt(0.4) * e[, 1], b = sqrt(0.6) * L + sqrt(0.4) * e[, 2])
  list(X = X, xy = xy, B = structure(B, class = "depprep_basis", type = "isotropic", k = 10, kn = 10 / n), e = e)
}

test_that("conditional correction with the true shared basis recovers the residual correlation", {
  d <- exact_case(0.5)
  cr <- correct(d$X, c(a = "g1", b = "g1"), coords = d$xy, bases = list(g1 = d$B))
  expect_s3_class(cr, "depprep_corrected")
  # (I - P) L = 0, so the result is the correlation of the projected true residuals
  E <- depprep:::.resid_on(d$e, d$B)
  expect_equal(cr$R[1, 2], stats::cor(E[, 1], E[, 2]))
  expect_lt(abs(cr$R[1, 2] - stats::cor(d$e[, 1], d$e[, 2])), 0.03)
  expect_lt(abs(stats::cor(d$X)[1, 2] - stats::cor(d$e[, 1], d$e[, 2])) - 0.2, 0.1)  # raw: inflated
})

test_that("both outputs are valid: a correlation matrix and residuals orthogonal to the union", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("b", 3, "matern", range = 0.15),
                  sim_group("none", 2, "none"))
  s <- sim_spatial(300, g, seed = 2)
  grp <- c(rep("g1", 3), rep("g2", 3), rep("none", 2)); names(grp) <- colnames(s$X)
  bases <- list(g1 = build_group_basis(s$X[, 1:3], coords = s$coords, type = "isotropic", seed = 1),
                g2 = build_group_basis(s$X[, 4:6], coords = s$coords, type = "isotropic", seed = 1))
  cr <- correct(s, grp, bases = bases)
  expect_equal(cr$R, t(cr$R))
  expect_equal(unname(diag(cr$R)), rep(1, 8))
  expect_gt(min(eigen(cr$R, only.values = TRUE)$values), -1e-8)
  expect_equal(dim(cr$X), dim(s$X))
  U <- unclass(union_basis(bases$g1, bases$g2))
  expect_lt(max(abs(crossprod(U, cr$X))), 1e-8)
  expect_equal(nrow(cr$union_kn), 3)                      # g1-g2, g1-none, g2-none
  expect_identical(correct(s, grp, bases = bases), cr)
  expect_error(correct(s, c(x = "g1")), "every column")
})

test_that("without any basis the correction changes nothing", {
  g <- sim_groups(sim_group("none", 4, "none"))
  s <- sim_spatial(200, g, seed = 3)
  grp <- stats::setNames(rep("none", 4), colnames(s$X))
  cr <- correct(s, grp)
  expect_equal(cr$R, stats::cor(s$X), tolerance = 1e-10, ignore_attr = TRUE)
  v <- verify(cr, s)
  expect_lt(max(abs(v$delta_R)), 1e-10)
})

test_that("whitening removes the spatial autocorrelation of the group", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.3))
  s <- sim_spatial(300, g, seed = 4)
  grp <- stats::setNames(rep("g1", 3), colnames(s$X))
  b <- list(g1 = build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 1))
  cw <- correct(s, grp, bases = b, mode = "whiten")
  lw <- depprep:::.knn_weights(s$coords)
  before <- depprep:::.moran(s$X[, 1], lw)[["moran_I"]]
  after <- depprep:::.moran(cw$X[, 1], lw)[["moran_I"]]
  expect_gt(before, 0.2)
  expect_lt(abs(after), 0.08)
})

test_that("temporal whitening filters the shared persistence", {
  g <- sim_groups(sim_group("ar", 3, "ar1", phi = 0.9), sim_group("seas", 3, "seasonal", period = 12))
  s <- sim_temporal(300, g, seed = 5)
  gr <- data.frame(variable = colnames(s$X), group = rep(c("g1", "g2"), each = 3),
                   group_type = rep(c("ar", "seasonal"), each = 3))
  class(gr) <- c("depprep_grouping", "data.frame"); attr(gr, "period") <- c(g1 = NA, g2 = 12)
  cw <- correct(s, gr, mode = "whiten")
  acf1 <- function(z) stats::acf(z, lag.max = 1, plot = FALSE)$acf[2]
  expect_lt(abs(acf1(cw$X[, 1])), 0.15)
  expect_lt(abs(acf1(cw$X[, 4])), 0.15)
  expect_true(all(is.finite(cw$R)))
})

test_that("the scale mode gives one valid correlation matrix per band", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5))
  s <- sim_spatial(200, g, seed = 6)
  grp <- stats::setNames(rep("g1", 3), colnames(s$X))
  cs <- correct(s, grp, bases = list(g1 = build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 1)),
                mode = "scale")
  expect_named(cs$R_scale, c("broad", "medium", "fine"))
  for (r in cs$R_scale) expect_equal(unname(diag(r)), rep(1, 3))
  expect_error(verify(cs, s), "conditional")
})
