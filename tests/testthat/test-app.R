test_that("the preview of each technique gives scores or classes for the raw and corrected data", {
  df <- depprep:::.app_example()
  d <- depprep:::.app_data(df, "spatial", x = "x", y = "y")
  cr <- correct(d$X, group_regimes(diagnose_variables(d$X, coords = d$coords)), coords = d$coords, seed = 1)
  for (tq in c("pca", "mds", "kmeans", "som", "hdbscan")) {
    p <- depprep:::.technique_preview(tq, d$X, cr$X, cr$R, k = 3)
    expect_identical(p$technique, tq)
    expect_identical(p$rows, seq_len(nrow(d$X)))
    for (part in c("raw", "corrected")) {
      if (tq %in% c("pca", "mds")) expect_equal(dim(p[[part]]$xy), c(nrow(d$X), 2)) else expect_length(p[[part]]$cl, nrow(d$X))
    }
    if (tq %in% c("kmeans", "som")) expect_setequal(unique(p$raw$cl), 1:3)
    f <- tempfile(fileext = ".pdf")
    expect_silent(depprep:::.save_figure(f, function() depprep:::.draw_preview(p, "spatial", d$coords, cex = 1.3, pt = 0.8), "pdf"))
    expect_gt(file.size(f), 1000)
  }
  # seeded: the same preview twice; a subsample above max_n
  expect_identical(depprep:::.technique_preview("kmeans", d$X, cr$X, k = 3), depprep:::.technique_preview("kmeans", d$X, cr$X, k = 3))
  ps <- depprep:::.technique_preview("mds", d$X, cr$X, max_n = 50)
  expect_length(ps$rows, 50)
  expect_equal(nrow(ps$raw$xy), 50)
  # the PCA of the corrected data uses the axes of the corrected R
  pp <- depprep:::.technique_preview("pca", d$X, cr$X, cr$R)
  expect_equal(pp$corrected$share, eigen(cr$R, symmetric = TRUE)$values[1:2] / ncol(cr$R))
})

test_that("figures are saved in every format at the size asked", {
  draw <- function() depprep:::.draw_autocor(data.frame(variable = c("a", "b"), autocor_before = c(0.8, 0.5), autocor_after = c(0.1, -0.05)),
                                             cex = 1.2, pt = 1.5)
  for (fmt in c("png", "pdf", "tiff", if (isTRUE(capabilities("cairo"))) "svg")) {
    f <- tempfile(fileext = paste0(".", fmt))
    depprep:::.save_figure(f, draw, fmt, width_cm = 12, height_cm = 8, dpi = 150)
    expect_gt(file.size(f), 500)
  }
  f <- tempfile(fileext = ".png")
  depprep:::.save_figure(f, draw, "png", width_cm = 12, height_cm = 8, dpi = 150)
  sig <- readBin(f, "raw", 24)
  expect_equal(c(readBin(sig[17:20], "integer", size = 4, endian = "big"), readBin(sig[21:24], "integer", size = 4, endian = "big")),
               c(round(12 / 2.54 * 150), round(8 / 2.54 * 150)), tolerance = 1)
  expect_identical(depprep:::.class_colours(c(0, 1, 9)), c("grey80", "#E69F00", "#E69F00"))
})

test_that("every plot of the diagnosis is drawn for its family", {
  set.seed(1)
  sp <- depprep:::.app_data(depprep:::.app_example(), "spatial", x = "x", y = "y")
  tm <- depprep:::.app_data(data.frame(t = 1:60, a = cumsum(stats::rnorm(60)), b = sin(1:60 / 2) + stats::rnorm(60)), "temporal", time = "t")
  st_df <- data.frame(site = rep(1:20, times = 12), t = rep(1:12, each = 20), x = rep(stats::runif(20), 12), y = rep(stats::runif(20), 12))
  st_df$a <- st_df$x + st_df$t / 6 + stats::rnorm(240); st_df$b <- stats::rnorm(240)
  st <- depprep:::.app_data(st_df, "spatiotemporal", x = "x", y = "y", time = "t", site = "site")
  expect_identical(unname(depprep:::.diagnosis_plots("spatial")), c("variogram", "correlogram", "map"))
  for (d in list(sp, tm, st)) for (w in depprep:::.diagnosis_plots(d$family)) {
    f <- tempfile(fileext = ".pdf")
    expect_silent(depprep:::.save_figure(f, function() depprep:::.draw_diagnosis(d, colnames(d$X)[1], 1.1, 0.8, w), "pdf"))
    expect_gt(file.size(f), 1000)
    expect_true(nzchar(depprep:::.diagnosis_plot_help[[w]]))
  }
  expect_error(depprep:::.draw_diagnosis(tm, "a", what = "map"), "should be one of")
})

test_that("the graphs of the menu are drawn and the menus put in words", {
  set.seed(3)
  xy <- cbind(lon = stats::runif(40), lat = stats::runif(40))
  links <- lapply(c(dbmem = "dbmem", del = "del", gab = "gab", rel = "rel"), function(w) depprep:::.menu_graph_links(xy, w))
  key <- function(e) paste(e[, 1], e[, 2])
  expect_true(all(vapply(links, function(e) all(e[, 1] < e[, 2]), TRUE)))
  # the relative-neighbourhood graph is inside the Gabriel graph, which is inside the Delaunay graph
  expect_true(all(key(links$rel) %in% key(links$gab)))
  expect_true(all(key(links$gab) %in% key(links$del)))
  # the distance-based graph links every site
  expect_setequal(unique(c(links$dbmem)), 1:40)
  expect_identical(depprep:::.menu_graph_links(xy, "auto"), links$del)
  expect_identical(depprep:::.menu_graph_links(xy, "gab_flin"), links$gab)
  expect_true(all(names(depprep:::.graph_help) %in% depprep:::.graph_choices))
  f <- tempfile(fileext = ".pdf")
  expect_silent(depprep:::.save_figure(f, function() depprep:::.draw_menu_graph(xy, "gab", direction = 45), "pdf"))
  expect_gt(file.size(f), 1000)
  expect_identical(depprep:::.menu_words(NULL, "spatial"), "automatic choice")
  expect_identical(depprep:::.menu_words(basis_menu(term_trend(), term_mem(W = "gab"), term_aem(30)), "spatial"),
                   "trend + spatial MEM (gabriel graph) + directional maps toward 30 degrees")
  expect_identical(depprep:::.menu_words(basis_menu(term_cycle(), term_mem()), "temporal"), "cycle (period from the periodogram) + temporal MEM")
  expect_identical(depprep:::.menu_words(depprep:::.menu_of_type("isotropic"), "spatial"), "trend + spatial MEM (distance-based)")
})

test_that(".format_table() shows 3 significant digits and leaves missing values blank", {
  d <- data.frame(a = c("x", NA), p = c(0.005, NA), big = c(35100, 0.0273), ok = c(TRUE, NA), stringsAsFactors = FALSE)
  f <- depprep:::.format_table(d)
  expect_identical(f$a, c("x", ""))
  expect_identical(f$p, c("0.005", ""))
  expect_identical(f$ok, c("yes", ""))
  expect_false(any(grepl("NA", unlist(f))))
  expect_identical(dim(f), dim(d))
})

test_that(".coord_like() finds the coordinate columns, with or without a unit", {
  nm <- c("x_km", "y_km", "lon", "lat", "Longitude", "x", "easting", "chl", "sst", "lat2", "max")
  expect_identical(nm[depprep:::.coord_like(nm)], c("x_km", "y_km", "lon", "lat", "Longitude", "x", "easting"))
  expect_identical(nm[depprep:::.coord_like(nm, "x")], c("x_km", "lon", "Longitude", "x", "easting"))
  expect_identical(nm[depprep:::.coord_like(nm, "y")], c("y_km", "lat"))
})

test_that(".app_data() reads the positions of each family and checks the data", {
  df <- depprep:::.app_example()
  d <- depprep:::.app_data(df, "spatial", x = "x", y = "y")
  expect_equal(ncol(d$X), 10)
  expect_equal(dim(d$coords), c(200, 2))
  expect_false(any(c("x", "y") %in% colnames(d$X)))
  expect_length(d$notes, 0)
  deg <- df; deg$x <- deg$x * 50 - 50; deg$y <- deg$y * 10
  expect_match(depprep:::.app_data(deg, "spatial", x = "x", y = "y")$notes, "degrees")
  # near the equator degrees are accepted; far from it the app asks to project them
  expect_match(depprep:::.app_data(deg, "spatial", x = "x", y = "y")$notes, "Near the equator this is fine")
  deg$y <- deg$y + 50
  expect_match(depprep:::.app_data(deg, "spatial", x = "x", y = "y")$notes, "Project them")
  gappy <- df; gappy[3, 4] <- NA
  expect_match(depprep:::.app_data(gappy, "spatial", x = "x", y = "y")$notes, "of the values are missing")
  bad <- df; bad$x[2] <- NA
  expect_error(depprep:::.app_data(bad, "spatial", x = "x", y = "y"), "`coords`")
  expect_error(depprep:::.app_data(df, "temporal"), "time column")
  expect_error(depprep:::.app_data(df, "spatial", x = "x", y = "y", variables = "long_1"), "at least two")
  st <- data.frame(site = rep(c("a", "b", "c"), times = 4), t = rep(1:4, each = 3), x = rep(c(0, 1, 2), 4), y = rep(c(0, 0, 1), 4),
                   v1 = rnorm(12), v2 = rnorm(12))
  ds <- depprep:::.app_data(st, "spatiotemporal", x = "x", y = "y", time = "t", site = "site")
  expect_equal(nrow(ds$coords), 3)
  expect_equal(ds$index$site, rep(1:3, 4))
  expect_equal(unname(ds$coords[2, ]), c(1, 0))
})

test_that("regimes changed by hand and menus built in the app", {
  g <- data.frame(variable = c("a", "b", "c"), type = "isotropic", group = c("g1", "g1", "g2"), group_type = "isotropic")
  class(g) <- c("depprep_grouping", "data.frame"); attr(g, "period") <- c(g1 = NA, g2 = NA)
  o <- depprep:::.override_regime(g, c("a", "c"), "anisotropic")
  expect_equal(o$group, c("manual_anisotropic", "g1", "manual_anisotropic"))
  expect_equal(o$group_type[1], "anisotropic")
  expect_equal(depprep:::.override_regime(g, "b", "none")$group[2], "none")
  expect_error(depprep:::.override_regime(g, "z", "none"), "Unknown variable")
  expect_error(depprep:::.override_regime(g, "a", "seasonal"), "period")
  expect_equal(unname(attr(depprep:::.override_regime(g, "a", "seasonal", 12), "period")["manual_seasonal"]), 12)
  expect_null(depprep:::.app_menu())
  m <- depprep:::.app_menu(trend = TRUE, mem = TRUE, W = "gab", aem = TRUE, direction = 30)
  expect_s3_class(m, "depprep_menu")
  expect_equal(vapply(m, `[[`, "", "kind"), c("trend", "mem", "aem"))
  expect_equal(depprep:::.menu_code(m), "basis_menu(term_trend(), term_mem(W = \"gab\"), term_aem(direction = 30))")
})

test_that("the R code of the app parses and reproduces the correction", {
  skip_on_cran()
  code <- depprep:::.app_code("data.csv", "spatial", "x", "y", NULL, NULL, c("long_1", "long_2", "none_1"),
                              list(list(variables = "none_1", type = "none", period = NULL)),
                              list(g1 = basis_menu(term_trend())), "conditional", "groups")
  expect_silent(parse(text = code))
  expect_true(any(grepl("menu = menus", code)))
  both <- depprep:::.app_code("data.csv", "spatiotemporal", "x", "y", "t", "s", c("a", "b"), list(), list(), "conditional_ar", "groups",
                              list(mode = "whiten", scope = "groups"))
  expect_silent(parse(text = both))
  expect_true(any(grepl("mode = \"whiten\", scope = \"groups\", data_basis = \"union\", bases = cr$bases", both, fixed = TRUE)))
  expect_true(any(grepl("grouping$group[i] <- \"none\"", code, fixed = TRUE)))
  skip_on_cran()
  dir <- tempfile(); dir.create(dir); dir <- normalizePath(dir, winslash = "/")
  df <- depprep:::.app_example()
  utils::write.csv(df, file.path(dir, "data.csv"), row.names = FALSE)
  code <- depprep:::.app_code(file.path(dir, "data.csv"), "spatial", "x", "y", NULL, NULL, setdiff(names(df), c("x", "y")),
                              list(), list(), "conditional", "global")
  code <- sub("depprep_report.html", file.path(dir, "r.html"), code, fixed = TRUE)
  code <- gsub("\"corrected_", paste0("\"", dir, "/corrected_"), code, fixed = TRUE)
  code <- gsub("\"diagnosis.csv\"", paste0("\"", dir, "/diagnosis.csv\""), code, fixed = TRUE)
  env <- new.env()
  eval(parse(text = code), envir = env)
  expect_true(file.exists(file.path(dir, "r.html")))
  expect_true(file.exists(file.path(dir, "diagnosis.csv")))
  expect_equal(dim(env$cr$X), c(200, 10))
})

test_that("the question of step 2 gives the mode of each family", {
  expect_setequal(unname(depprep:::.question_choices), c("association", "innovations", "unsure"))
  for (q in c("association", "unsure")) {
    expect_identical(depprep:::.question_mode(q, "spatial"), "conditional")
    expect_identical(depprep:::.question_mode(q, "temporal"), "conditional")
    expect_identical(depprep:::.question_mode(q, "spatiotemporal"), "conditional_ar")
  }
  for (fam in c("spatial", "temporal", "spatiotemporal")) expect_identical(depprep:::.question_mode("innovations", fam), "whiten")
})

test_that("the two modes of the app are compared pair by pair", {
  R1 <- diag(3); dimnames(R1) <- list(c("a", "b", "c"), c("a", "b", "c"))
  R2 <- R1
  R1["a", "b"] <- R1["b", "a"] <- 0.1; R2["a", "b"] <- R2["b", "a"] <- 0.25
  R1["a", "c"] <- R1["c", "a"] <- -0.5; R2["a", "c"] <- R2["c", "a"] <- 0.1
  R1["b", "c"] <- R1["c", "b"] <- 0.3; R2["b", "c"] <- R2["c", "b"] <- 0.5
  p <- depprep:::.mode_disagreement(R1, R2)
  expect_identical(p$a, c("a", "b"))
  expect_identical(p$b, c("c", "c"))
  expect_equal(p$difference, c(0.6, 0.2))
  expect_identical(attr(p, "n_pairs"), 3L)
  expect_identical(nrow(depprep:::.mode_disagreement(R1, R1)), 0L)
  expect_identical(depprep:::.other_mode("whiten", "spatiotemporal"), "conditional_ar")
  expect_identical(depprep:::.other_mode("conditional_ar", "spatiotemporal"), "whiten")
  expect_identical(depprep:::.other_mode("whiten", "temporal"), "conditional")
  expect_identical(depprep:::.default_scope("conditional", "spatial"), "global")
  expect_identical(depprep:::.default_scope("conditional", "spatial", list(g1 = 1)), "groups")
  expect_identical(depprep:::.default_scope("whiten", "spatial"), "groups")
  expect_identical(depprep:::.default_scope("conditional_ar", "spatiotemporal"), "groups")
})

test_that("the other mode gives the R of a correction run alone, with the bases reused in space-time", {
  skip_on_cran()
  g <- sim_groups(sim_group("sep", 3, "separable_ar1", range = 0.5, phi = 0.9), sim_group("none", 2, "none"))
  s <- sim_spatiotemporal(30, 40, g, seed = 8)
  gr <- data.frame(variable = colnames(s$X), group = rep(c("g1", "none"), c(3, 2)), group_type = rep(c("separable_ar1", "none"), c(3, 2)))
  pos <- list(coords = s$coords, index = s$index)
  first <- correct(s$X, gr, coords = s$coords, index = s$index, scope = "groups", seed = 1)
  other <- depprep:::.run_other_mode(first, s$X, gr, pos, "spatiotemporal")
  expect_identical(other$mode, "whiten")
  expect_equal(other$R, correct(s$X, gr, coords = s$coords, index = s$index, mode = "whiten", scope = "groups", seed = 1)$R)
  # a spatial pair of modes has different scopes, so the whitening builds its own bases
  d <- depprep:::.app_data(depprep:::.app_example(), "spatial", x = "x", y = "y")
  grs <- group_regimes(diagnose_variables(d$X, coords = d$coords, nsim = 0))
  f <- correct(d$X, grs, coords = d$coords, seed = 1)
  o <- depprep:::.run_other_mode(f, d$X, grs, list(coords = d$coords), "spatial")
  expect_identical(o$scope, "groups")
  expect_identical(o$mode, "whiten")
})

test_that("the predicted time holds the times measured on real and simulated data", {
  inside <- function(rows, cores, seconds) { m <- depprep:::.predicted_minutes(rows, cores); seconds / 60 >= m[["low"]] && seconds / 60 <= m[["high"]] }
  expect_true(inside(20925, 16, 525))     # the shelf of the Foz, default mode (analysis/09h)
  expect_true(inside(37479, 16, 223))     # the slope, conditional + AR (analysis/09g)
  expect_true(inside(58404, 16, 3349))    # the whole basin, conditional + AR (analysis/09d)
  expect_true(inside(9300, 1, 196))       # analysis/08f
  expect_true(inside(27900, 16, 582))
  w <- depprep:::.time_words
  expect_identical(w(c(low = 0.1, high = 0.5)), "under a minute")
  expect_identical(w(c(low = 2, high = 13)), "about 2 to 13 minutes")
  expect_identical(w(c(low = 60, high = 300)), "about 1 to 5 hours")
  expect_match(depprep:::.time_message(20925, "spatiotemporal", "conditional_ar", 16), "reuses the bases")
  expect_match(depprep:::.time_message(20925, "spatiotemporal", "conditional_ar", 16, "verify"), "pair by pair")
  expect_match(depprep:::.time_message(20925, "spatiotemporal", "whiten", 16, "verify"), "no projected part")
  expect_match(depprep:::.time_message(500, "spatial", "conditional"), "seconds to a few minutes")
})

test_that("the app server loads the example and builds the data", {
  skip_if_not_installed("shiny")
  shiny::testServer(depprep:::.app_server, {
    session$setInputs(example = 1, family = "spatial")
    expect_equal(nrow(rv$df), 200)
    session$setInputs(x = "x", y = "y", variables = setdiff(names(rv$df), c("x", "y")), technique = "x", regions = "yes",
                      question = "association", mode = "default", scope = "default")
    expect_equal(ncol(dat()$X), 10)
    expect_equal(settings()$scope, "global")
    expect_equal(settings()$mode, "conditional")
    session$setInputs(question = "innovations")
    expect_equal(settings()$mode, "whiten")
    session$setInputs(mode = "conditional")           # step 4 can still overrule the answer
    expect_equal(settings()$mode, "conditional")
  })
  expect_s3_class(depprep:::.app_ui(), "shiny.tag.list")
})

test_that("the previews and the diagnosis plots of the app accept missing values", {
  df <- depprep:::.app_example()
  df[c(3, 40:49), 4] <- NA
  d <- depprep:::.app_data(df, "spatial", x = "x", y = "y")
  gr <- group_regimes(diagnose_variables(d$X, coords = d$coords, nsim = 0))
  cr <- correct(d$X, gr, coords = d$coords, seed = 1)
  for (tq in c("pca", "kmeans")) {
    p <- depprep:::.technique_preview(tq, d$X, cr$X, cr$R, k = 3)
    expect_false(any(p$rows %in% c(3, 40:49)))
  }
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f))
  v <- colnames(d$X)[3]
  for (w in depprep:::.diagnosis_plots("spatial")) {
    expect_silent(depprep:::.save_figure(f, function() depprep:::.draw_diagnosis(d, v, 1, 1, w), "pdf"))
  }
  tm <- data.frame(t = 1:60, a = sin(1:60 / 3) + stats::rnorm(60, sd = 0.1), b = stats::rnorm(60))
  tm$a[10:14] <- NA
  dt <- depprep:::.app_data(tm, "temporal", time = "t")
  for (w in depprep:::.diagnosis_plots("temporal")) {
    expect_silent(depprep:::.save_figure(f, function() depprep:::.draw_diagnosis(dt, "a", 1, 1, w), "pdf"))
  }
})
