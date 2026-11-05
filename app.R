# Opens the app from a copy of the source, without installing the package, as in
#   shiny::runGitHub("depprep", "gfonseca-unifesp", ref = "main")
# It needs shiny, pkgload and the packages listed under Imports in DESCRIPTION. This file
# is not part of the package; with the package installed, use depprep::run_app().
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The app needs the pkgload package: install.packages(\"pkgload\").", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
shiny::shinyApp(depprep:::.app_ui(), depprep:::.app_server)
