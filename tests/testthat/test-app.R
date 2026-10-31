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
  bad <- df; bad[3, 4] <- NA
  expect_error(depprep:::.app_data(bad, "spatial", x = "x", y = "y"), "missing values")
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
  expect_true(any(grepl("grouping$group[i] <- \"none\"", code, fixed = TRUE)))
  skip_on_cran()
  dir <- tempfile(); dir.create(dir); dir <- normalizePath(dir, winslash = "/")
  df <- depprep:::.app_example()
  utils::write.csv(df, file.path(dir, "data.csv"), row.names = FALSE)
  code <- depprep:::.app_code(file.path(dir, "data.csv"), "spatial", "x", "y", NULL, NULL, setdiff(names(df), c("x", "y")),
                              list(), list(), "conditional", "global")
  code <- sub("depprep_report.html", file.path(dir, "r.html"), code, fixed = TRUE)
  code <- gsub("\"corrected_", paste0("\"", dir, "/corrected_"), code, fixed = TRUE)
  env <- new.env()
  eval(parse(text = code), envir = env)
  expect_true(file.exists(file.path(dir, "r.html")))
  expect_equal(dim(env$cr$X), c(200, 10))
})

test_that("the app server loads the example and builds the data", {
  skip_if_not_installed("shiny")
  shiny::testServer(depprep:::.app_server, {
    session$setInputs(example = 1, family = "spatial")
    expect_equal(nrow(rv$df), 200)
    session$setInputs(x = "x", y = "y", variables = setdiff(names(rv$df), c("x", "y")), technique = "x", regions = "yes",
                      mode = "default", scope = "default")
    expect_equal(ncol(dat()$X), 10)
    expect_equal(settings()$scope, "global")
  })
  expect_s3_class(depprep:::.app_ui(), "shiny.tag.list")
})
