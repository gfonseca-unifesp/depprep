test_that("report_correction() writes one self-contained HTML file", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(120, g, seed = 2)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  cr <- correct(s, grp, seed = 1)
  v <- verify(cr, s, nperm = 0)
  f <- report_correction(cr, s, v, file = tempfile(fileext = ".html"))
  expect_true(file.exists(f))
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (s_ in c("<h2>Bases</h2>", "What was removed", "Delta R", "(Q3)", "LISA", "Reproduce", "correct(X, grouping")) {
    expect_true(grepl(s_, h, fixed = TRUE), info = s_)
  }
  expect_false(grepl("<script src=|<link ", h))                     # nothing external
  if (isTRUE(capabilities("cairo"))) expect_true(grepl("<svg", h, fixed = TRUE))
  # temporal data: no LISA section; given code is shown escaped
  gt <- sim_groups(sim_group("a", 2, "ar1", phi = 0.8), sim_group("none", 1, "none"))
  st <- sim_temporal(100, gt, seed = 1)
  grt <- data.frame(variable = colnames(st$X), group = c("g1", "g1", "none"), group_type = c("ar", "ar", "none"))
  class(grt) <- c("depprep_grouping", "data.frame"); attr(grt, "period") <- c(g1 = NA)
  ft <- report_correction(correct(st, grt), st, file = tempfile(fileext = ".html"), code = "x <- 1 < 2")
  ht <- paste(readLines(ft, warn = FALSE), collapse = "\n")
  expect_false(grepl("LISA", ht, fixed = TRUE))
  expect_true(grepl("x &lt;- 1 &lt; 2", ht, fixed = TRUE))
})
