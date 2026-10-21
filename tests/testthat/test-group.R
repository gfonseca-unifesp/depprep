spatial_case <- function(seed = 21) {
  g <- sim_groups(sim_group("long", 8, "exponential", range = 0.5),
                  sim_group("short", 8, "exponential", range = 0.1),
                  sim_group("none", 7, "none"))
  sim_spatial(400, g, seed = seed)
}

test_that("grouping recovers the true groups and keeps the none group apart", {
  s <- spatial_case()
  sig <- diagnose_variables(s)
  gr <- group_regimes(sig, method = "hclust")
  expect_s3_class(gr, "depprep_grouping")
  expect_gte(mclust::adjustedRandIndex(gr$group, s$group), 0.8)
  expect_true(all(gr$group[gr$type == "none"] == "none"))
  expect_true(all(gr$group_type[gr$group == "none"] == "none"))
  expect_true(all(is.na(gr$stability)))
  gm <- group_regimes(sig)                                   # default: mclust (decision D6)
  expect_identical(attr(gm, "method"), "mclust")
  expect_gte(mclust::adjustedRandIndex(gm$group, s$group), 0.4)
  expect_true(all(gm$group[gm$type == "none"] == "none"))
  expect_equal(attr(group_regimes(sig, method = "hclust", n_groups = 2), "n_groups"), 2)
})

test_that("grouping uses only the signatures, not the order or correlations of the variables", {
  s <- spatial_case()
  sig <- diagnose_variables(s)
  gr <- group_regimes(sig)
  o <- sample.int(nrow(sig))
  gr_o <- group_regimes(sig[o, ])
  expect_equal(mclust::adjustedRandIndex(gr_o$group, gr$group[o]), 1)
  # the partition is a function of the signature table alone
  expect_identical(group_regimes(sig)$group, gr$group)
})

test_that("stability: co-association is a symmetric share and defines the stability", {
  s <- spatial_case(22)
  sig <- diagnose_variables(s, nsim = 99)
  # hclust: in spatial data mclust splits more and its stability is lower (stage 4)
  gr <- group_regimes(sig, data = s, method = "hclust", B = 8, seed = 3)
  co <- attr(gr, "coassociation")
  expect_equal(co, t(co))
  expect_true(all(diag(co) == 1))
  expect_true(all(co >= 0 & co <= 1))
  dep <- gr[gr$group != "none", ]
  i <- 1
  same <- setdiff(which(dep$group == dep$group[i]), i)
  expect_equal(dep$stability[i], mean(co[i, same]))
  expect_identical(gr$ambiguous, !is.na(gr$stability) & gr$stability < 0.7)
  expect_gte(stats::median(dep$stability), 0.7)                  # clean case: stable groups
  expect_identical(group_regimes(sig, data = s, method = "hclust", B = 8, seed = 3), gr)
  expect_error(group_regimes(sig, B = 5), "required")
})

test_that("temporal and spatio-temporal grouping run and report periods", {
  g <- sim_groups(sim_group("seas", 5, "seasonal", period = 12), sim_group("ar", 5, "ar1", phi = 0.9),
                  sim_group("none", 4, "none"))
  s <- sim_temporal(240, g, seed = 6)
  gr <- group_regimes(diagnose_variables(s), data = s, B = 4, seed = 1)
  expect_gte(mclust::adjustedRandIndex(gr$group, s$group), 0.8)
  seas_group <- gr$group[gr$group_type == "seasonal"][1]
  expect_equal(unname(attr(gr, "period")[seas_group]), 12)
  g2 <- sim_groups(sim_group("sep", 4, "separable_ar1", range = 0.5, phi = 0.9),
                   sim_group("common", 4, "temporal_common", phi = 0.9))
  s2 <- sim_spatiotemporal(50, 60, g2, seed = 4)
  gr2 <- group_regimes(diagnose_variables(s2, nsim = 99), data = s2, B = 2, seed = 1)
  expect_equal(nrow(gr2), 8)
  shared <- gr2$group != "none" & (duplicated(gr2$group) | duplicated(gr2$group, fromLast = TRUE))
  expect_false(anyNA(gr2$stability[shared]))
})

test_that("basis_distance is 0 for the same subspace and 1 for orthogonal ones", {
  Q <- qr.Q(qr(matrix(stats::rnorm(400), 100)))
  b1 <- structure(Q[, 1:2], class = "depprep_basis")
  b2 <- structure(Q[, 3:4], class = "depprep_basis")
  b3 <- structure(Q[, 2:1], class = "depprep_basis")
  expect_equal(basis_distance(b1, b3), 0)
  expect_equal(basis_distance(b1, b2), 1)
  expect_equal(basis_distance(b1, b2), basis_distance(b2, b1))
  empty <- structure(matrix(0, 100, 0), class = "depprep_basis")
  expect_equal(basis_distance(empty, empty), 0)
  expect_equal(basis_distance(b1, empty), 1)
})

test_that("the subspace level separates variables of the same regime with independent fields", {
  g <- sim_groups(sim_group("a1", 4, "matern", range = 0.5), sim_group("a2", 4, "matern", range = 0.5),
                  sim_group("none", 3, "none"))
  s <- sim_spatial(400, g, seed = 1)
  sig <- diagnose_variables(s, nsim = 99)
  gr <- group_regimes(sig, data = s, method = "type", split = "subspace", split_threshold = 0.55, seed = 1)
  expect_equal(mclust::adjustedRandIndex(gr$group, s$group), 1)
  D <- attr(gr, "basis_distance")[[1]]
  expect_equal(D, t(D))
  expect_true(all(diag(D) == 0))
  # same field: closer than independent fields (on average; single pairs can overlap)
  f <- sub("_[0-9]+$", "", rownames(D))
  expect_lt(mean(D[outer(f, f, "==") & upper.tri(D)]), mean(D[outer(f, f, "!=")]))
  o <- sample.int(nrow(sig))
  s_o <- s; s_o$X <- s$X[, o]
  gr_o <- group_regimes(sig[o, ], data = s_o, method = "type", split = "subspace", split_threshold = 0.55, seed = 1)
  expect_equal(mclust::adjustedRandIndex(gr_o$group, gr$group[o]), 1)
  expect_error(group_regimes(sig, split = "subspace"), "required")
})
