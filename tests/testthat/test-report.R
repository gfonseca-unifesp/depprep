test_that("report_correction() writes one self-contained HTML file", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(120, g, seed = 2)
  grp <- stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X))
  cr <- correct(s, grp, seed = 1)
  v <- verify(cr, s, nperm = 0)
  f <- report_correction(cr, s, v, file = tempfile(fileext = ".html"))
  expect_true(file.exists(f))
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (s_ in c("Diagnosis before the correction", "Moran&#39;s I|Moran's I", "<h2>Bases</h2>", "What was removed", "Correlogram p before",
               "LISA before", "Moran correlograms before", "Delta R", "(Q3)", "Local structure before and after", "Reproduce",
               "correct(X, grouping", "signature = sig")) {
    if (grepl("|", s_, fixed = TRUE)) { expect_true(grepl(s_, h), info = s_); next }
    expect_true(grepl(s_, h, fixed = TRUE), info = s_)
  }
  expect_false(grepl("<script src=|<link ", h))                     # nothing external
  if (isTRUE(capabilities("cairo"))) {
    expect_true(grepl("<svg", h, fixed = TRUE))
    # every figure names its glyphs apart, so no id is defined twice in the page
    ids <- regmatches(h, gregexpr('id="[^"]+"', h))[[1]]
    expect_gt(length(ids), 0)
    expect_false(anyDuplicated(ids) > 0)
  }
  # temporal data: no LISA section; given code is shown escaped
  gt <- sim_groups(sim_group("a", 2, "ar1", phi = 0.8), sim_group("none", 1, "none"))
  st <- sim_temporal(100, gt, seed = 1)
  grt <- data.frame(variable = colnames(st$X), group = c("g1", "g1", "none"), group_type = c("ar", "ar", "none"))
  class(grt) <- c("depprep_grouping", "data.frame"); attr(grt, "period") <- c(g1 = NA)
  ft <- report_correction(correct(st, grt), st, file = tempfile(fileext = ".html"), code = "x <- 1 < 2")
  ht <- paste(readLines(ft, warn = FALSE), collapse = "\n")
  expect_false(grepl("LISA", ht, fixed = TRUE))
  expect_true(grepl("x &lt;- 1 &lt; 2", ht, fixed = TRUE))
  expect_true(grepl("Ljung-Box p before", ht, fixed = TRUE))
  expect_true(grepl("Diagnosis before the correction", ht, fixed = TRUE))
})

test_that(".before_after() lays out each variable's tests before and after, by family", {
  g <- sim_groups(sim_group("a", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(120, g, seed = 2)
  v <- verify(correct(s, stats::setNames(c("g1", "g1", "g1", "none", "none"), colnames(s$X)), seed = 1), s, nperm = 0)
  b <- depprep:::.before_after(v, "spatial")
  expect_identical(b$variable, colnames(s$X))
  expect_identical(b$lisa_before, signif(v$variables$lisa_share_before, 3))
  expect_identical(b$correlogram_p_after, signif(v$variables$correlogram_p, 3))
  expect_length(attr(b, "labels"), ncol(b))
  gt <- sim_groups(sim_group("a", 2, "ar1", phi = 0.8))
  st <- sim_temporal(100, gt, seed = 1)
  grt <- data.frame(variable = colnames(st$X), group = "g1", group_type = "ar")
  class(grt) <- c("depprep_grouping", "data.frame"); attr(grt, "period") <- c(g1 = NA)
  bt <- depprep:::.before_after(verify(correct(st, grt), st, nperm = 0), "temporal")
  expect_true(all(c("lb_p_before", "lb_p_after") %in% names(bt)))
  expect_false("lisa_before" %in% names(bt))
})

test_that("diagnosis_table() joins the signature of the family with the grouping", {
  g <- sim_groups(sim_group("long", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
  s <- sim_spatial(100, g, seed = 1)
  sig <- diagnose_variables(s)
  gr <- group_regimes(sig)
  d <- diagnosis_table(sig, gr)
  expect_identical(d$variable, sig$variable)
  expect_identical(d$group, gr$group[match(sig$variable, gr$variable)])
  expect_true(all(c("type", "stability", "moran_I", "moran_p", "range", "nugget_ratio") %in% names(d)))
  expect_false(any(c("phi", "period", "ljung_p") %in% names(d)))           # temporal columns hold no value in space
  expect_equal(d$moran_I, signif(sig$moran_I, 3))
  expect_identical(diagnosis_table(sig, stats::setNames(gr$group, gr$variable))$group, d$group)
  expect_identical(names(diagnosis_table(sig))[1], "variable")
  st <- sim_temporal(120, sim_groups(sim_group("s", 2, "seasonal", period = 12), sim_group("none", 1, "none")), seed = 1)
  dt <- diagnosis_table(diagnose_variables(st))
  expect_true(all(c("acf1", "phi", "ljung_p", "period", "period_strength") %in% names(dt)))
  expect_false(any(c("moran_I", "range") %in% names(dt)))
})
