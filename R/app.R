#' Run the depprep app
#'
#' Opens a point-and-click interface to the whole workflow, in six steps:
#' 1. load a CSV file and mark the position columns (coordinates, time, site);
#' 2. say which technique the corrected data are for and whether the classes sought
#'    could be regions or periods;
#' 3. diagnose the dependence of each variable, with its plots (variogram, ACF,
#'    periodogram), and change a regime by hand or build the menu of structures of a
#'    group ([basis_menu()]);
#' 4. correct, with the default explained (the data matrix with one basis per variable,
#'    the correlation matrix with the shared basis) and the alternatives;
#' 5. verify what the correction removed and kept, with a preview of the technique
#'    before and after;
#' 6. download the corrected data, the correlation matrix, the report
#'    ([report_correction()]) and the R code that reproduces everything.
#'
#' Limits shown in the app: missing values are not supported yet; coordinates in degrees
#' should be projected first; large data (over a few thousand observations) are slow.
#' Needs the shiny package.
#'
#' @param launch.browser Open the app in the browser.
#' @return Runs the app; returns nothing useful.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(launch.browser = interactive()) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The app needs the shiny package: install.packages(\"shiny\").", call. = FALSE)
  }
  shiny::runApp(shiny::shinyApp(.app_ui(), .app_server), launch.browser = launch.browser)
}

# --- Pure helpers (tested without the app) --------------------------------------------

# The data of the app from a data frame and the chosen columns: the numeric variables,
# the positions of the family, and the checks the app shows (missing values stop;
# coordinates that look like degrees and large data give notes).
.app_data <- function(df, family, x = NULL, y = NULL, time = NULL, site = NULL, variables = NULL) {
  family <- match.arg(family, c("spatial", "temporal", "spatiotemporal"))
  pos_cols <- c(x, y, time, site)
  if (is.null(variables)) variables <- setdiff(names(df)[vapply(df, is.numeric, TRUE)], pos_cols)
  variables <- setdiff(variables, pos_cols)
  if (length(variables) < 2) stop("Choose at least two numeric variables.", call. = FALSE)
  X <- as.matrix(df[, variables, drop = FALSE])
  if (!is.numeric(X)) stop("The variables must be numeric.", call. = FALSE)
  notes <- character(0)
  out <- list(X = X, coords = NULL, time = NULL, index = NULL, family = family)
  if (family %in% c("spatial", "spatiotemporal")) {
    if (is.null(x) || is.null(y)) stop("Choose the x and y columns.", call. = FALSE)
    xy <- as.matrix(df[, c(x, y)])
    lonlat_names <- any(grepl("^(lon|long|longitude|lat|latitude)$", c(x, y), ignore.case = TRUE))
    in_bounds <- all(abs(xy[, 1]) <= 180, na.rm = TRUE) && all(abs(xy[, 2]) <= 90, na.rm = TRUE)
    if (in_bounds && (lonlat_names || diff(range(xy[, 1], na.rm = TRUE)) > 5)) {
      notes <- c(notes, "The coordinates look like degrees: distances in degrees are distorted away from the equator; project them (for example to UTM, in km) first.")
    }
  }
  if (family == "spatial") out$coords <- xy
  if (family == "temporal") {
    if (is.null(time)) stop("Choose the time column.", call. = FALSE)
    out$time <- as.numeric(df[[time]])
  }
  if (family == "spatiotemporal") {
    if (is.null(time) || is.null(site)) stop("Choose the time and site columns.", call. = FALSE)
    s <- match(df[[site]], unique(df[[site]]))
    first <- !duplicated(s)
    out$coords <- xy[first, , drop = FALSE][order(s[first]), , drop = FALSE]
    out$index <- data.frame(site = s, time = as.numeric(df[[time]]))
  }
  .stop_if_missing(X, if (family == "spatial") out$coords else NULL, out$time)
  if (nrow(X) > 3000) notes <- c(notes, sprintf("%d observations: the diagnosis and the correction can take minutes.", nrow(X)))
  out$notes <- notes
  out
}

# Change the regime of some variables by hand: each gets the group "manual_<type>" (or
# "none"), with that type; a seasonal type takes `period`.
.override_regime <- function(grouping, variables, type, period = NULL) {
  i <- match(variables, grouping$variable)
  if (anyNA(i)) stop("Unknown variable: ", paste(variables[is.na(i)], collapse = ", "), call. = FALSE)
  g <- if (type == "none") "none" else paste0("manual_", type)
  grouping$type[i] <- type
  grouping$group[i] <- g
  grouping$group_type[i] <- type
  if (type %in% c("seasonal", "separable_seasonal")) {
    if (is.null(period)) stop("A seasonal regime needs a period.", call. = FALSE)
    p <- attr(grouping, "period")
    attr(grouping, "period") <- c(p[setdiff(names(p), g)], stats::setNames(period, g))
  }
  grouping
}

# The menu built in the app (NULL when nothing is ticked).
.app_menu <- function(trend = FALSE, cycle = FALSE, period = NULL, mem = FALSE, W = "dbmem", time_mem = FALSE,
                      aem = FALSE, direction = 0) {
  terms <- c(if (trend) list(term_trend()), if (cycle) list(term_cycle(period = if (is.null(period) || is.na(period)) NULL else period)),
             if (mem) list(term_mem(W = W)), if (time_mem) list(term_mem(margin = "time")), if (aem) list(term_aem(direction)))
  if (!length(terms)) NULL else do.call(basis_menu, terms)
}

# The R code that reproduces what was done in the app.
.app_code <- function(file, family, x, y, time, site, variables, overrides, menus, mode, scope) {
  q <- function(v) encodeString(v, quote = "\"")
  vec <- function(v) paste0("c(", paste(q(v), collapse = ", "), ")")
  pos <- switch(family, spatial = "coords = coords", temporal = "time = time", spatiotemporal = "coords = coords, index = index")
  positions <- switch(family,
    spatial = sprintf("coords <- as.matrix(df[, %s])", vec(c(x, y))),
    temporal = sprintf("time <- df[[%s]]", q(time)),
    spatiotemporal = c(sprintf("site <- match(df[[%s]], unique(df[[%s]]))", q(site), q(site)),
                       sprintf("index <- data.frame(site = site, time = df[[%s]])", q(time)),
                       sprintf("coords <- as.matrix(df[!duplicated(site), %s])[order(site[!duplicated(site)]), ]", vec(c(x, y)))))
  override <- function(o) {
    g <- if (o$type == "none") "none" else paste0("manual_", o$type)
    c(sprintf("i <- match(%s, grouping$variable)", vec(o$variables)),
      sprintf("grouping$type[i] <- %s; grouping$group[i] <- %s; grouping$group_type[i] <- %s", q(o$type), q(g), q(o$type)),
      if (!is.null(o$period)) sprintf("attr(grouping, \"period\")[%s] <- %s", q(g), o$period))
  }
  c("library(depprep)",
    sprintf("df <- read.csv(%s)", q(file)),
    sprintf("X <- as.matrix(df[, %s])", vec(variables)),
    positions,
    sprintf("sig <- diagnose_variables(X, %s)", pos),
    "grouping <- group_regimes(sig)",
    if (length(overrides)) c("# regimes changed by hand", unlist(lapply(overrides, override))),
    if (length(menus)) c("menus <- list(", paste0("  ", names(menus), " = ", vapply(menus, function(m) .menu_code(m), ""), collapse = ",\n"), ")"),
    sprintf("cr <- correct(X, grouping, %s, mode = %s, scope = %s%s)", pos, q(mode), q(scope), if (length(menus)) ", menu = menus" else ""),
    sprintf("v <- verify(cr, X, %s)", pos),
    "write.csv(cr$X, \"corrected_data.csv\", row.names = FALSE)",
    "write.csv(cr$R, \"corrected_R.csv\")",
    sprintf("report_correction(cr, X, v, file = \"depprep_report.html\", %s)", pos))
}

.menu_code <- function(m) {
  one <- function(t) {
    par <- t[setdiff(names(t), c("kind", "margin"))]
    par <- par[!vapply(par, is.null, TRUE)]
    if (!is.null(t$margin) && t$kind %in% c("mem", "trend")) par$margin <- t$margin
    if (t$kind == "aem") par$margin <- NULL
    args <- vapply(names(par), function(k) sprintf("%s = %s", k, if (is.character(par[[k]])) paste0("\"", par[[k]], "\"") else par[[k]]), "")
    sprintf("term_%s(%s)", t$kind, paste(args, collapse = ", "))
  }
  sprintf("basis_menu(%s)", paste(vapply(m, one, ""), collapse = ", "))
}

# Example data for the app: a small spatial simulation with three regimes.
.app_example <- function() {
  g <- sim_groups(sim_group("long", 4, "exponential", range = 0.5), sim_group("short", 3, "exponential", range = 0.15),
                  sim_group("none", 3, "none"))
  s <- sim_spatial(200, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0),
                   classes = list(k = 3, type = "random", effect = 1), seed = 1)
  data.frame(x = s$coords[, 1], y = s$coords[, 2], s$X)
}

# --- User interface -------------------------------------------------------------------

.app_ui <- function() {
  shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(paste(
      "body{max-width:1200px;margin:auto}.note{background:#fff4e5;border-left:4px solid #f39c12;padding:.5rem;margin:.5rem 0}",
      ".ok{background:#eef7ee;border-left:4px solid #2e8b57;padding:.5rem;margin:.5rem 0}.help{color:#555;font-size:.95em}")))),
    shiny::titlePanel("depprep: dependence-aware preprocessing"),
    shiny::navlistPanel(id = "step", widths = c(2, 10),
      shiny::tabPanel("1. Data",
        shiny::p(class = "help", "Load a CSV with one row per observation. Mark the columns of the positions; the other numeric columns are the variables."),
        shiny::fileInput("file", "CSV file", accept = ".csv"),
        shiny::actionButton("example", "Use example data (spatial simulation)"),
        shiny::radioButtons("family", "Data", c("Spatial" = "spatial", "Temporal" = "temporal", "Space-time" = "spatiotemporal"), inline = TRUE),
        shiny::uiOutput("cols"),
        shiny::uiOutput("data_msg"),
        shiny::tableOutput("head")),
      shiny::tabPanel("2. Context",
        shiny::radioButtons("technique", "The corrected data are for",
          c("PCA (uses the corrected correlation matrix R)" = "pca", "MDS, k-means, SOM or HDBSCAN (use the corrected data matrix X)" = "x")),
        shiny::radioButtons("regions", "Could the classes you look for be regions or periods (contiguous in space or time)?",
          c("No" = "no", "Yes" = "yes", "Not sure" = "unsure")),
        shiny::uiOutput("context_msg")),
      shiny::tabPanel("3. Diagnosis",
        shiny::p(class = "help", "Each variable is diagnosed alone (variogram or ACF), never through its correlations with the others. Variables with similar signatures form a group that shares one basis."),
        shiny::actionButton("diagnose", "Diagnose"),
        shiny::tableOutput("signatures"),
        shiny::fluidRow(shiny::column(4, shiny::uiOutput("plot_var")), shiny::column(8, shiny::plotOutput("diag_plot", height = 300))),
        shiny::h4("Change a regime by hand"),
        shiny::uiOutput("override_ui"),
        shiny::h4("Menu of structures for a group"),
        shiny::p(class = "help", "Restrict the structures a group's basis is built from; the parameters are still searched with the usual limits (tests, k/n cap). Empty: the automatic choice."),
        shiny::uiOutput("menu_ui"),
        shiny::verbatimTextOutput("menus")),
      shiny::tabPanel("4. Correction",
        shiny::uiOutput("correct_help"),
        shiny::selectInput("mode", "Mode", c("Default for the data" = "default", "Conditional (residuals on the bases)" = "conditional",
                                         "Whitening (innovations)" = "whiten")),
        shiny::selectInput("scope", "Bases for R", c("Default for the data" = "default", "One per group" = "groups", "One global basis" = "global",
                                                  "Union of the variables' bases" = "variables")),
        shiny::actionButton("correct", "Correct"),
        shiny::tableOutput("bases"),
        shiny::uiOutput("correct_msg")),
      shiny::tabPanel("5. Verification",
        shiny::actionButton("verify", "Verify"),
        shiny::uiOutput("q3_msg"),
        shiny::tableOutput("delta"),
        shiny::plotOutput("autocor_plot", height = 320),
        shiny::h4("Preview of the technique, before and after"),
        shiny::numericInput("k", "Number of groups (k-means)", 3, min = 2, max = 20),
        shiny::plotOutput("preview", height = 360)),
      shiny::tabPanel("6. Export",
        shiny::p(class = "help", "Download the results and the R code that reproduces them."),
        shiny::downloadButton("dl_x", "Corrected data (CSV)"),
        shiny::downloadButton("dl_r", "Corrected correlation matrix (CSV)"),
        shiny::downloadButton("dl_report", "Report (HTML)"),
        shiny::downloadButton("dl_code", "R code"),
        shiny::hr(),
        shiny::p(class = "help", "Limits: missing values are not supported yet; project coordinates in degrees first; filters such as SAR/CAR or STAR are not in the menu."))))
}

# --- Server ---------------------------------------------------------------------------

.app_server <- function(input, output, session) {
  rv <- shiny::reactiveValues(df = NULL, file = "data.csv", sig = NULL, grouping = NULL, overrides = list(), menus = list(),
                          cr = NULL, v = NULL)
  shiny::observeEvent(input$file, { rv$df <- utils::read.csv(input$file$datapath); rv$file <- input$file$name })
  shiny::observeEvent(input$example, {
    rv$df <- .app_example(); rv$file <- "example.csv"
    shiny::updateRadioButtons(session, "family", selected = "spatial")
  })
  num_cols <- shiny::reactive({ shiny::req(rv$df); names(rv$df)[vapply(rv$df, is.numeric, TRUE)] })
  output$cols <- shiny::renderUI({
    shiny::req(rv$df)
    nm <- names(rv$df)
    guess <- function(p) { m <- grep(p, nm, ignore.case = TRUE, value = TRUE); if (length(m)) m[1] else NULL }
    shiny::tagList(
      if (input$family != "temporal") shiny::fluidRow(shiny::column(3, shiny::selectInput("x", "x column", nm, guess("^(x|lon|long|longitude|easting)$"))),
                                                   shiny::column(3, shiny::selectInput("y", "y column", nm, guess("^(y|lat|latitude|northing)$")))),
      if (input$family != "spatial") shiny::fluidRow(shiny::column(3, shiny::selectInput("time", "time column", nm, guess("^(t|time|date|month|year)$"))),
        if (input$family == "spatiotemporal") shiny::column(3, shiny::selectInput("site", "site column", nm, guess("^(site|station|pixel|id)$")))),
      shiny::selectizeInput("variables", "Variables", num_cols(), setdiff(num_cols(), c(input$x, input$y, input$time, input$site)), multiple = TRUE))
  })
  dat <- shiny::reactive({
    shiny::req(rv$df, input$variables)
    .app_data(rv$df, input$family, x = if (input$family != "temporal") input$x, y = if (input$family != "temporal") input$y,
              time = if (input$family != "spatial") input$time, site = if (input$family == "spatiotemporal") input$site,
              variables = input$variables)
  })
  output$data_msg <- shiny::renderUI({
    d <- tryCatch(dat(), error = function(e) e)
    if (inherits(d, "error")) return(shiny::div(class = "note", conditionMessage(d)))
    shiny::tagList(shiny::div(class = "ok", sprintf("%d observations, %d variables.", nrow(d$X), ncol(d$X))),
               lapply(d$notes, function(n) shiny::div(class = "note", n)))
  })
  output$head <- shiny::renderTable(utils::head(rv$df, 6))
  output$context_msg <- shiny::renderUI({
    out <- list(shiny::p(if (input$technique == "pca") "PCA reads the corrected correlation matrix R, built with the shared basis of each group." else
      "These techniques read the corrected data matrix X, in which each variable is corrected with its own basis: it keeps the classes best."))
    if (input$regions != "no") out <- c(out, list(shiny::div(class = "note",
      "Classes that are regions or periods live at the scale of the basis, and the correction removes them together with the dependence. Step 5 warns when much of the association is at that scale (Q3); in that case use the raw data, or the scale mode of correct(), for those classes.")))
    shiny::tagList(out)
  })
  pos_args <- function(d) switch(d$family, spatial = list(coords = d$coords), temporal = list(time = d$time),
                                 spatiotemporal = list(coords = d$coords, index = d$index))
  shiny::observeEvent(input$diagnose, {
    d <- dat()
    shiny::withProgress(message = "Diagnosing each variable", {
      rv$sig <- do.call(diagnose_variables, c(list(d$X, nsim = 199), pos_args(d)))
      rv$grouping <- group_regimes(rv$sig)
      rv$overrides <- list(); rv$menus <- list(); rv$cr <- NULL; rv$v <- NULL
    })
  })
  output$signatures <- shiny::renderTable({
    shiny::req(rv$grouping, rv$sig)
    s <- rv$sig; g <- rv$grouping
    cols <- intersect(c("dependence_p", "range", "nugget_ratio", "aniso_ratio", "phi", "period", "trend_r2", "shift_indicator"), names(s))
    keep <- cols[vapply(cols, function(k) any(is.finite(s[[k]])), TRUE)]
    cbind(g[, c("variable", "type", "group")], signif(s[match(g$variable, s$variable), keep, drop = FALSE], 3))
  })
  output$plot_var <- shiny::renderUI({ shiny::req(input$variables); shiny::selectInput("pvar", "Plot of", input$variables) })
  output$diag_plot <- shiny::renderPlot({
    d <- dat(); shiny::req(input$pvar %in% colnames(d$X))
    z <- d$X[, input$pvar]
    if (d$family == "spatial") {
      v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(z = z, x = d$coords[, 1], y = d$coords[, 2]))
      graphics::plot(v$dist, v$gamma, pch = 19, xlab = "distance", ylab = "semivariance", main = paste("Variogram of", input$pvar))
    } else if (d$family == "temporal") {
      o <- order(d$time)
      graphics::par(mfrow = c(1, 2))
      stats::acf(z[o], main = paste("ACF of", input$pvar))
      sp <- stats::spec.pgram(z[o], taper = 0, detrend = TRUE, plot = FALSE)
      graphics::plot(1 / sp$freq, sp$spec, type = "h", log = "x", xlab = "period", ylab = "power", main = "Periodogram")
    } else {
      graphics::par(mfrow = c(1, 2))
      m <- tapply(z, d$index$site, mean)
      graphics::plot(d$coords, pch = 19, col = grDevices::hcl.colors(20)[cut(m, 20)], asp = 1, main = "Site means", xlab = "x", ylab = "y")
      stats::acf(tapply(z, d$index$time, mean), main = "ACF of the mean series")
    }
  })
  types_of <- function(f) switch(f, spatial = c("isotropic", "anisotropic", "none"), temporal = c("ar", "seasonal", "trend", "none"),
    spatiotemporal = c("temporal_common", "spatial_only", "moving", "separable_seasonal", "separable_ar1", "none"))
  output$override_ui <- shiny::renderUI({
    shiny::req(rv$grouping)
    shiny::fluidRow(shiny::column(4, shiny::selectizeInput("ov_vars", "Variables", rv$grouping$variable, multiple = TRUE)),
                shiny::column(3, shiny::selectInput("ov_type", "New regime", types_of(input$family))),
                shiny::column(2, shiny::numericInput("ov_period", "Period (seasonal)", NA)),
                shiny::column(2, shiny::br(), shiny::actionButton("ov_apply", "Apply")))
  })
  shiny::observeEvent(input$ov_apply, {
    shiny::req(rv$grouping, input$ov_vars)
    per <- if (is.na(input$ov_period)) NULL else input$ov_period
    g <- tryCatch(.override_regime(rv$grouping, input$ov_vars, input$ov_type, per), error = function(e) e)
    if (inherits(g, "error")) { shiny::showNotification(conditionMessage(g), type = "error"); return() }
    rv$grouping <- g
    rv$overrides <- c(rv$overrides, list(list(variables = input$ov_vars, type = input$ov_type, period = per)))
  })
  output$menu_ui <- shiny::renderUI({
    shiny::req(rv$grouping)
    gs <- setdiff(unique(rv$grouping$group), "none")
    shiny::tagList(
      shiny::fluidRow(shiny::column(3, shiny::selectInput("mn_group", "Group", gs)),
                  shiny::column(2, shiny::checkboxInput("mn_trend", "Trend")),
                  if (input$family != "spatial") shiny::column(2, shiny::checkboxInput("mn_cycle", "Cycle"), shiny::numericInput("mn_period", "Period (empty: searched)", NA)),
                  if (input$family != "temporal") shiny::column(2, shiny::checkboxInput("mn_mem", "MEM"), shiny::selectInput("mn_W", "Graph", c("dbmem", "auto", "del", "gab", "rel"))),
                  if (input$family != "spatial") shiny::column(1, shiny::checkboxInput("mn_tmem", "Temporal MEM")),
                  if (input$family != "temporal") shiny::column(2, shiny::checkboxInput("mn_aem", "Directional (AEM)"), shiny::numericInput("mn_dir", "Direction (degrees)", 0))),
      shiny::actionButton("mn_apply", "Set the menu of this group"), shiny::actionButton("mn_clear", "Back to automatic"))
  })
  shiny::observeEvent(input$mn_apply, {
    m <- tryCatch(.app_menu(trend = isTRUE(input$mn_trend), cycle = isTRUE(input$mn_cycle), period = input$mn_period, mem = isTRUE(input$mn_mem),
                            W = if (is.null(input$mn_W)) "dbmem" else input$mn_W, time_mem = isTRUE(input$mn_tmem), aem = isTRUE(input$mn_aem),
                            direction = if (is.null(input$mn_dir)) 0 else input$mn_dir), error = function(e) e)
    if (inherits(m, "error")) { shiny::showNotification(conditionMessage(m), type = "error"); return() }
    menus <- rv$menus; menus[[input$mn_group]] <- m; rv$menus <- Filter(Negate(is.null), menus)
  })
  shiny::observeEvent(input$mn_clear, { menus <- rv$menus; menus[[input$mn_group]] <- NULL; rv$menus <- menus })
  output$menus <- shiny::renderPrint({ if (!length(rv$menus)) cat("All groups: automatic choice.\n") else for (g in names(rv$menus)) { cat(g, ": "); print(rv$menus[[g]]) } })
  settings <- shiny::reactive({
    d <- dat()
    mode <- if (input$mode == "default") (if (d$family == "spatiotemporal") "whiten" else "conditional") else input$mode
    scope <- if (input$scope == "default") (if (mode != "whiten" && d$family != "spatiotemporal" && !length(rv$menus)) "global" else "groups") else input$scope
    list(mode = mode, scope = scope)
  })
  output$correct_help <- shiny::renderUI(shiny::p(class = "help",
    "Default: in space and in time, R from one global basis of the dependent variables, conditional mode; in space-time, whitening by group. The data matrix X always has each variable corrected with its own basis, which keeps the classes. With a menu, the bases are built per group."))
  shiny::observeEvent(input$correct, {
    shiny::req(rv$grouping)
    d <- dat(); st <- settings()
    shiny::withProgress(message = "Correcting", {
      rv$cr <- tryCatch(do.call(correct, c(list(d$X, rv$grouping, mode = st$mode, scope = st$scope, menu = if (length(rv$menus)) rv$menus),
                                           pos_args(d))), error = function(e) e)
      rv$v <- NULL
    })
  })
  output$bases <- shiny::renderTable({
    shiny::req(rv$cr, !inherits(rv$cr, "error"))
    b <- rv$cr$bases
    data.frame(basis = names(b), k = vapply(b, ncol, 0L), kn = round(vapply(b, ncol, 0L) / nrow(dat()$X), 3),
               terms = vapply(b, function(x) paste(utils::head(attr(x, "terms"), 6), collapse = ", "), ""))
  })
  output$correct_msg <- shiny::renderUI({
    shiny::req(rv$cr)
    if (inherits(rv$cr, "error")) return(shiny::div(class = "note", conditionMessage(rv$cr)))
    shiny::tagList(shiny::div(class = "ok", sprintf("Corrected: mode %s, scope %s.", rv$cr$mode, rv$cr$scope)),
               lapply(rv$cr$warnings, function(w) shiny::div(class = "note", w)))
  })
  shiny::observeEvent(input$verify, {
    shiny::req(rv$cr, !inherits(rv$cr, "error"))
    d <- dat()
    shiny::withProgress(message = "Verifying", rv$v <- do.call(verify, c(list(rv$cr, d$X, nperm = 99), pos_args(d))))
  })
  output$q3_msg <- shiny::renderUI({
    shiny::req(rv$v); q <- rv$v$q3
    if (!is.finite(q$statistic)) return(shiny::div(class = "note", "Q3 not computed: whitening has no projected part."))
    if (isTRUE(q$warning)) shiny::div(class = "note", sprintf("Q3 warning: much of the association between variables (%.2f, threshold %.2f) was at the scale of the basis. Classes that are regions or periods would be removed.", q$statistic, q$threshold)) else
      shiny::div(class = "ok", sprintf("No Q3 warning (%.2f, threshold %.2f).", q$statistic, q$threshold))
  })
  output$delta <- shiny::renderTable({ shiny::req(rv$v); rv$v$delta_summary }, digits = 3)
  output$autocor_plot <- shiny::renderPlot({
    shiny::req(rv$v); vv <- rv$v$variables; o <- seq_len(nrow(vv))
    graphics::par(mar = c(4, 8, 1, 1), las = 1)
    graphics::plot(vv$autocor_before, o, xlim = range(c(vv$autocor_before, vv$autocor_after, 0)), yaxt = "n", pch = 19, col = "grey55",
                   xlab = "autocorrelation (grey: before; blue: after)", ylab = "")
    graphics::points(vv$autocor_after, o, pch = 19, col = "#2a6fb0")
    graphics::axis(2, at = o, labels = vv$variable, cex.axis = 0.8)
    graphics::abline(v = 0, lty = 3)
  })
  output$preview <- shiny::renderPlot({
    shiny::req(rv$cr, !inherits(rv$cr, "error"))
    d <- dat()
    graphics::par(mfrow = c(1, 2))
    if (input$technique == "pca") {
      for (k in 1:2) {
        R <- if (k == 1) stats::cor(d$X) else rv$cr$R
        e <- eigen(R, symmetric = TRUE)
        graphics::plot(e$vectors[, 1:2], type = "n", xlab = sprintf("axis 1 (%.0f%%)", 100 * e$values[1] / sum(e$values)),
                       ylab = sprintf("axis 2 (%.0f%%)", 100 * e$values[2] / sum(e$values)), main = c("PCA, raw", "PCA, corrected")[k])
        graphics::text(e$vectors[, 1:2], colnames(d$X), cex = 0.7)
      }
    } else {
      Xc <- rv$cr$X
      for (k in 1:2) {
        Z <- scale(if (k == 1) d$X[utils::tail(seq_len(nrow(d$X)), nrow(Xc)), ] else Xc)
        cl <- stats::kmeans(Z, input$k, nstart = 10)$cluster
        pc <- stats::prcomp(Z)$x[, 1:2]
        graphics::plot(pc, col = cl, pch = 19, cex = 0.6, main = c("k-means, raw", "k-means, corrected")[k])
      }
    }
  })
  output$dl_x <- shiny::downloadHandler("corrected_data.csv", function(f) utils::write.csv(rv$cr$X, f, row.names = FALSE))
  output$dl_r <- shiny::downloadHandler("corrected_R.csv", function(f) utils::write.csv(rv$cr$R, f))
  output$dl_report <- shiny::downloadHandler("depprep_report.html", function(f) {
    d <- dat()
    do.call(report_correction, c(list(rv$cr, d$X, rv$v, file = f, code = code_text()), pos_args(d)))
  })
  code_text <- shiny::reactive({
    d <- dat(); st <- settings()
    .app_code(rv$file, d$family, if (d$family != "temporal") input$x, if (d$family != "temporal") input$y,
              if (d$family != "spatial") input$time, if (d$family == "spatiotemporal") input$site, colnames(d$X),
              rv$overrides, rv$menus, st$mode, st$scope)
  })
  output$dl_code <- shiny::downloadHandler("depprep_analysis.R", function(f) writeLines(code_text(), f))
  invisible(rv)
}
