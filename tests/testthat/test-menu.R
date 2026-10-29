test_that("menus are built from terms and checked", {
  m <- basis_menu(term_trend(), term_cycle(period = 12), term_mem(margin = "time"))
  expect_s3_class(m, "depprep_menu")
  expect_length(m, 3)
  expect_output(print(m), "cycle\\(period = 12\\)")
  expect_error(basis_menu(), "one or more terms")
  expect_error(basis_menu(list(kind = "x")), "one or more terms")
  expect_error(term_mem(W = "voronoi"), "Unknown `W`")
  expect_error(term_cycle(period = 1), "at least 2")
  # a period of 2 (the alternation of a unit-step series) is accepted, as before the menu
  expect_silent(term_cycle(period = 2))
  expect_error(term_aem("north"), "angle")
  expect_identical(depprep:::.term_margins(m, "spatiotemporal"), c("both", "time", "time"))
  expect_identical(depprep:::.term_margins(basis_menu(term_mem()), "temporal"), "time")
})

test_that("the menu of a type gives the same basis as the type", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.5))
  s <- sim_spatial(120, g, seed = 3)
  b_type <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 2)
  b_menu <- build_group_basis(s$X, coords = s$coords, menu = basis_menu(term_trend(), term_mem()), seed = 2)
  expect_equal(unclass(b_menu)[, ], unclass(b_type)[, ])
  expect_identical(attr(b_menu, "terms"), attr(b_type, "terms"))
  expect_identical(attr(b_menu, "type"), "menu")
  expect_s3_class(attr(b_menu, "menu"), "depprep_menu")
  gt <- sim_groups(sim_group("s", 3, "seasonal", period = 12))
  st <- sim_temporal(144, gt, seed = 3)
  expect_equal(unclass(build_group_basis(st$X, time = st$time, menu = basis_menu(term_trend(), term_cycle(12)), seed = 1))[, ],
               unclass(build_group_basis(st$X, time = st$time, type = "seasonal", period = 12, seed = 1))[, ])
})

test_that("term_cycle() finds the period, and the terms are restricted to the menu", {
  gt <- sim_groups(sim_group("s", 3, "seasonal", period = 12))
  st <- sim_temporal(144, gt, seed = 4)
  b <- build_group_basis(st$X, time = st$time, menu = basis_menu(term_cycle()), seed = 1)
  expect_true(all(grepl("^(sin|cos)[0-9]+_p1[12]", attr(b, "terms"))))
  expect_gt(attr(b, "k"), 0)
  expect_error(build_group_basis(st$X, time = st$time, menu = basis_menu(term_aem(0))), "Spatial terms")
  expect_error(depprep:::.menu_period(matrix(0, 50, 2), 1:50), "give `period`")
})

test_that("term_mem() chooses the graph and term_aem() gives directional maps", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.4))
  s <- sim_spatial(100, g, seed = 5)
  # six graphs: the global test runs at the Sidak level 0.0085, which needs more than 117 permutations
  b <- build_group_basis(s$X, coords = s$coords, menu = basis_menu(term_mem(W = "auto")), seed = 1, nperm = 199)
  expect_true(attr(b, "W") %in% c("del_binary", "del_flin", "gab_binary", "gab_flin", "rel_binary", "rel_flin"))
  expect_true(all(grepl(paste0("_", attr(b, "W"), "$"), attr(b, "terms"))))
  bg <- build_group_basis(s$X, coords = s$coords, menu = basis_menu(term_mem(W = "gab")), seed = 1, nperm = 99)
  expect_true(all(grepl("_gab_binary$", attr(bg, "terms"))))
  a <- depprep:::.aem_vectors(s$coords, 30)
  expect_equal(nrow(a), 100)
  expect_true(all(grepl("^AEM[0-9]+_d30$", colnames(a))))
  # only the positively autocorrelated maps are kept, so the global test keeps its power
  expect_lt(ncol(a), 60)
  W <- spdep::nb2mat(spdep::tri2nb(s$coords), style = "B")
  moran <- 100 / sum(W) * colSums(a * (W %*% a)) / colSums(a^2)
  expect_true(all(moran > -1 / 99))
  ba <- build_group_basis(s$X, coords = s$coords, menu = basis_menu(term_aem(30)), seed = 1, nperm = 99)
  expect_gt(length(attr(ba, "terms")), 2)
  expect_true(all(grepl("^AEM", attr(ba, "terms"))))
  expect_lte(attr(ba, "kn"), 0.25)
})

test_that("correct() takes one menu or a menu per group", {
  g <- sim_groups(sim_group("a", 3, "exponential", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(120, g, seed = 6)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  m <- basis_menu(term_trend(degree = 1))
  cr <- correct(s, grp, menu = m, scope = "groups")
  expect_true(all(attr(cr$bases$g1, "terms") %in% c("trend_x", "trend_y")))
  cr2 <- correct(s, grp, menu = list(g1 = m))
  expect_true(all(attr(cr2$bases$global, "terms") %in% c("trend_x", "trend_y")))
  expect_error(correct(s, grp, menu = list(m)), "named by group")
  expect_error(correct(s, grp, menu = m, mode = "scale"), "scale mode")
})

test_that("a seasonal basis with a diagnosed period of 2 is built as before the menu", {
  set.seed(3)
  tt <- 1:93
  x <- matrix(cos(pi * tt) + stats::rnorm(93, sd = 0.5), ncol = 1)
  b <- build_group_basis(x, time = tt, type = "seasonal", period = 2, seed = 1)
  expect_s3_class(b, "depprep_basis")
  expect_true(all(grepl("_p2$|^trend", colnames(unclass(b)))))
})

test_that("the AEM maps of the same positions and direction are built once", {
  s <- sim_spatial(60, sim_groups(sim_group("a", 2, "exponential", range = 0.4)), seed = 4)
  cache <- depprep:::.aem_cache
  cache$entries <- NULL
  a1 <- depprep:::.aem_vectors(s$coords, 30)
  expect_length(cache$entries, 1)
  expect_identical(depprep:::.aem_vectors(s$coords, 30), a1)
  expect_length(cache$entries, 1)
  a2 <- depprep:::.aem_vectors(s$coords, 120)
  expect_length(cache$entries, 2)
  expect_false(isTRUE(all.equal(a1, a2)))
  expect_identical(depprep:::.aem_build(s$coords, 30), a1)
})
