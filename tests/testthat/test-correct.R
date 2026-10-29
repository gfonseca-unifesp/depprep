exact_case <- function(rho, n = 400, seed = 1) {
  # large scale exactly in the span of 10 dbMEM, shared by two variables
  set.seed(seed)
  xy <- matrix(stats::runif(2 * n), n)
  B <- depprep:::.spatial_candidates(xy, trend = FALSE)[, 1:10]
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
  cr <- correct(s, grp, bases = bases, data_basis = "union")
  expect_equal(cr$R, t(cr$R))
  expect_equal(unname(diag(cr$R)), rep(1, 8))
  expect_gt(min(eigen(cr$R, only.values = TRUE)$values), -1e-8)
  expect_equal(dim(cr$X), dim(s$X))
  U <- unclass(union_basis(bases$g1, bases$g2))
  expect_lt(max(abs(crossprod(U, cr$X))), 1e-8)
  expect_equal(nrow(cr$union_kn), 3)                      # g1-g2, g1-none, g2-none
  expect_identical(correct(s, grp, bases = bases, data_basis = "union"), cr)
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

test_that("the scale mode gives one valid correlation matrix per band of the selected functions", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5))
  s <- sim_spatial(200, g, seed = 6)
  grp <- stats::setNames(rep("g1", 3), colnames(s$X))
  b <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 1)
  cs <- correct(s, grp, bases = list(g1 = b), mode = "scale")
  expect_identical(names(cs$R_scale), cs$bands$band)
  expect_true(all(cs$bands$band %in% c("trend", "all", "broad", "medium", "fine")))
  for (r in cs$R_scale) expect_equal(unname(diag(r)), rep(1, 3))
  # every selected function is in exactly one band
  expect_setequal(unlist(strsplit(cs$bands$terms, ", ")), attr(b, "terms"))
  expect_equal(dim(cs$r2_scale), c(3, nrow(cs$bands)))
  expect_true(all(cs$r2_scale >= 0 & cs$r2_scale <= 1))
  expect_identical(cs$scalogram$g1, attr(b, "scalogram"))
  expect_error(verify(cs, s), "conditional")
})

test_that("a spatial pattern at a known scale lands in the right band, and an empty band is flagged", {
  set.seed(31)
  xy <- matrix(stats::runif(800), 400)
  C <- depprep:::.spatial_candidates(xy)
  Y <- cbind(c1 = 3 * C[, "MEM2"], c2 = 3 * C[, "MEM2"], f1 = 3 * C[, "MEM60"], f2 = 3 * C[, "MEM60"]) +
    matrix(stats::rnorm(1600, sd = 0.1), 400)
  grp <- c(c1 = "g", c2 = "g", f1 = "h", f2 = "h")
  cs <- correct(Y, grp, coords = xy, mode = "scale", scope = "groups")
  # selection may add a few noise terms; the two patterns must still fall in distinct bands
  band_of <- function(term) cs$bands$band[grepl(paste0("\\b", term, "\\b"), cs$bands$terms)]
  best <- colnames(cs$r2_scale)[apply(cs$r2_scale, 1, which.max)]
  expect_identical(best, rep(c(band_of("MEM2"), band_of("MEM60")), each = 2))
  expect_lt(match(band_of("MEM2"), cs$bands$band), match(band_of("MEM60"), cs$bands$band))
  expect_gt(min(cs$r2_scale[c("c1", "c2"), band_of("MEM2")]), 0.5)
  expect_gt(min(cs$r2_scale[c("f1", "f2"), band_of("MEM60")]), 0.5)
  # a band with a function that no variable follows is flagged
  B <- structure(C[, c("MEM2", "MEM60", "MEM120")], class = "depprep_basis", k = 3, kn = 3 / 400)
  Y0 <- Y - C[, "MEM120"] %*% crossprod(C[, "MEM120"], Y)            # nothing at the scale of MEM120
  cf <- correct(Y0, grp, coords = xy, mode = "scale", bases = list(g = B, h = B))
  expect_identical(cf$bands$band, c("broad", "medium", "fine"))
  expect_match(cf$warnings, "band fine: no variable has signal", all = FALSE)
  expect_false(any(grepl("band (broad|medium)", cf$warnings)))
})

test_that("a temporal pattern at a known wavelength lands in the right band", {
  set.seed(32)
  t <- 1:240
  slow <- 8 * depprep:::.temporal_candidates(t, "ar")[, "tMEM4"]
  seas <- sin(2 * pi * t / 12)
  X <- cbind(a1 = slow, a2 = slow, b1 = seas, b2 = seas) + matrix(stats::rnorm(960, sd = 0.5), 240)
  bs <- list(ga = build_group_basis(X[, 1:2], time = t, type = "ar", seed = 1),
             gb = build_group_basis(X[, 3:4], time = t, type = "seasonal", period = 12, seed = 1))
  cs <- correct(X, c(a1 = "ga", a2 = "ga", b1 = "gb", b2 = "gb"), time = t, bases = bs, mode = "scale")
  best <- colnames(cs$r2_scale)[apply(cs$r2_scale, 1, which.max)]
  expect_identical(best[1], best[2])
  expect_identical(best[3], best[4])
  ba <- cs$bands[cs$bands$band == best[1], ]
  bb <- cs$bands[cs$bands$band == best[3], ]
  expect_true(ba$wavelength_min > 60 && ba$wavelength_max < 150)   # tMEM4: about 98 time steps
  expect_true(bb$wavelength_min <= 12 && bb$wavelength_max >= 12)
  expect_lt(match(best[1], cs$bands$band), match(best[3], cs$bands$band))
})

test_that("scope variables: the same basis for every variable gives the same R as one group", {
  d <- exact_case(0.5)
  grp <- c(a = "g1", b = "g1")
  cv <- correct(d$X, grp, coords = d$xy, scope = "variables", bases = list(a = d$B, b = d$B), data_basis = "union")
  cg <- correct(d$X, grp, coords = d$xy, bases = list(g1 = d$B), data_basis = "union")
  expect_identical(cv$scope, "variables")
  expect_named(cv$bases, "union")
  expect_equal(attr(cv$bases$union, "k"), 10)
  expect_equal(cv$R, cg$R)
  expect_equal(cv$X, cg$X)
  expect_equal(correct(d$X, grp, coords = d$xy, scope = "variables", bases = list(a = d$B, b = d$B),
                       data_basis = "own")$X, cg$X, ignore_attr = TRUE)
  expect_error(correct(d$X, grp, coords = d$xy, scope = "variables", bases = list(a = d$B)), "one basis per dependent")
  # D17: in whitening mode the data matrix can still hold each variable on its own basis
  cw <- correct(d$X, grp, coords = d$xy, bases = list(g1 = d$B), mode = "whiten", data_basis = "own")
  expect_null(cw$variable_bases)                          # bases given: each variable uses its group's
  expect_equal(cw$X, depprep:::.resid_on(d$X, d$B), ignore_attr = TRUE)
})

test_that("scope variables: the union contains each variable's basis and every pair uses it", {
  g <- sim_groups(sim_group("a", 2, "matern", range = 0.5), sim_group("b", 2, "matern", range = 0.15),
                  sim_group("none", 2, "none"))
  s <- sim_spatial(300, g, seed = 2)
  grp <- stats::setNames(c("g1", "g1", "g2", "g2", "none", "none"), colnames(s$X))
  cv <- correct(s, grp, scope = "variables", seed = 3, data_basis = "union")
  expect_named(cv$variable_bases, colnames(s$X)[1:4])
  expect_true(all(cv$group[1:4] == "union") && all(cv$group[5:6] == "none"))
  U <- unclass(cv$bases$union)
  expect_true(isTRUE(all.equal(crossprod(U), diag(ncol(U)), check.attributes = FALSE)))
  inside <- function(b) max(abs(unclass(b) - U %*% crossprod(U, unclass(b))))
  for (b in cv$variable_bases) expect_lt(inside(b), 1e-8)
  # the between-basis error is zero: each pair is the correlation of residuals on the union
  E <- depprep:::.resid_on(s$X[, 1:4], U)
  expect_equal(cv$R_pairwise[1:4, 1:4], stats::cor(E), ignore_attr = TRUE)
  expect_lt(max(abs(crossprod(U, cv$X))), 1e-8)
  # data_basis = "own": each column on its own basis only
  co <- correct(s, grp, scope = "variables", bases = cv$variable_bases, data_basis = "own")
  expect_equal(co$R, cv$R)
  for (j in 1:4) expect_equal(co$X[, j], depprep:::.resid_on(s$X[, j, drop = FALSE], cv$variable_bases[[j]])[, 1],
                              ignore_attr = TRUE)
  expect_equal(co$X[, 5], s$X[, 5] - mean(s$X[, 5]), ignore_attr = TRUE)
  # kn_max caps the union
  cc <- correct(s, grp, scope = "variables", bases = cv$variable_bases, kn_max = 0.02)
  expect_lte(attr(cc$bases$union, "k"), 6)
  expect_match(cc$warnings, "capped", all = FALSE)
})

test_that("the period argument sets the seasonal period; a missing period warns", {
  g <- sim_groups(sim_group("seas", 3, "seasonal", period = 12))
  s <- sim_temporal(240, g, seed = 9)
  gr <- data.frame(variable = colnames(s$X), group = "g1", group_type = "seasonal")
  expect_warning(cr <- correct(s, gr, scope = "groups"), "No seasonal period known for group g1")
  expect_true(all(grepl("_p12$|^trend", attr(cr$bases$g1, "terms"))))
  expect_no_warning(c10 <- correct(s, gr, scope = "groups", period = 10))
  expect_true(all(grepl("_p10$|^trend", attr(c10$bases$g1, "terms"))))
  expect_no_warning(correct(s, gr, scope = "groups", period = c(g1 = 12)))
  expect_warning(correct(s, gr), "the global basis")
  expect_null(depprep:::.period_of("ar", "g1", NULL, gr))
  expect_equal(depprep:::.period_of("seasonal", NULL, c(a = 6, b = 12, c = 12), gr), 12)
})

test_that("defaults follow decision D7: global basis in space and time, whitening in space-time", {
  skip_on_cran()
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("b", 3, "matern", range = 0.15),
                  sim_group("none", 2, "none"))
  s <- sim_spatial(300, g, seed = 7)
  gr <- group_regimes(diagnose_variables(s, nsim = 99))
  cr <- correct(s, gr)
  expect_identical(cr$mode, "conditional")
  expect_identical(cr$scope, "global")
  expect_named(cr$bases, "global")
  dep <- gr$group != "none"
  expect_true(all(cr$group[dep] == "global"))
  # variables without dependence are not corrected among themselves
  nn <- which(!dep)
  if (length(nn) > 1) expect_equal(cr$R[nn[1], nn[2]], stats::cor(s$X)[nn[1], nn[2]], tolerance = 1e-8)
  expect_identical(correct(s, gr, scope = "groups")$scope, "groups")
  g2 <- sim_groups(sim_group("sep", 3, "separable_ar1", range = 0.5, phi = 0.9), sim_group("none", 2, "none"))
  s2 <- sim_spatiotemporal(30, 40, g2, seed = 8)
  gr2 <- group_regimes(diagnose_variables(s2, nsim = 99))
  cw <- correct(s2, gr2)
  expect_identical(cw$mode, "whiten")
  expect_identical(cw$scope, "groups")
})

test_that("the union of the variable bases works when it has more columns than sites", {
  g <- sim_groups(sim_group("a", 12, "exponential", range = 0.5), sim_group("b", 12, "exponential", range = 0.1))
  s <- sim_spatial(60, g, seed = 11)
  grp <- stats::setNames(rep("g1", ncol(s$X)), colnames(s$X))
  cand <- depprep:::.spatial_candidates(s$coords)
  vb <- lapply(seq_len(ncol(s$X)), function(j) {             # overlapping selections, as variables make
    cols <- depprep:::.with_stream(j, 1, sort(sample(ncol(cand), 15)))
    structure(cand[, cols], class = "depprep_basis", type = "isotropic", terms = colnames(cand)[cols])
  })
  names(vb) <- colnames(s$X)
  expect_gt(sum(vapply(vb, ncol, 0)), nrow(s$X))          # wider than tall before the union
  cr <- correct(s$X, grp, coords = s$coords, bases = vb, scope = "variables", kn_max = 0.3)
  u <- cr$bases$union
  expect_true(all(is.finite(unclass(u))))
  expect_lte(ncol(u), floor(0.3 * nrow(s$X)))
  expect_equal(crossprod(unclass(u)), diag(ncol(u)), tolerance = 1e-8, ignore_attr = TRUE)
  M <- cbind(diag(3), diag(3)); colnames(M) <- c("a", "b", "c", "a", "b", "x")
  expect_identical(depprep:::.distinct_columns(M, colnames(M)), c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))
})


test_that("by default the data matrix holds each variable on a basis selected on it alone (D17)", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(200, g, seed = 4)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  cr <- correct(s, grp, seed = 2)
  expect_identical(cr$data_basis, "own")
  expect_named(cr$variable_bases, colnames(s$X)[1:3])
  for (j in 1:3) {
    b <- build_group_basis(s$X[, j, drop = FALSE], coords = s$coords, type = "isotropic", seed = 2 + j)
    expect_equal(cr$X[, j], depprep:::.resid_on(s$X[, j, drop = FALSE], b)[, 1], ignore_attr = TRUE)
  }
  expect_equal(cr$X[, 4], s$X[, 4] - mean(s$X[, 4]), ignore_attr = TRUE)   # no dependence: centred only
})
