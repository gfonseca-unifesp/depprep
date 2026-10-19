test_that("the package loads and declares its core dependencies", {
  expect_true(isNamespaceLoaded("depprep"))
  imports <- read.dcf(system.file("DESCRIPTION", package = "depprep"), fields = "Imports")
  for (p in c("adespatial", "gstat", "spdep", "mclust", "kohonen")) {
    expect_match(imports, p)
  }
})
