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
#' Coordinates can be projected (for example in km) or longitude and latitude in degrees.
#' Degrees treat a degree of longitude as a degree of latitude, which holds near the
#' equator; the app says how much shorter a degree of longitude is at the largest latitude
#' of the data, and asks to project them when that exceeds 2% (beyond about 11 degrees of
#' latitude).
#'
#' Limits shown in the app: missing values are not supported yet; large data (over a few
#' thousand observations) are slow.
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
  if (is.null(variables)) variables <- setdiff(names(df)[vapply(df, is.numeric, TRUE) & !.coord_like(names(df))], pos_cols)
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
      # one degree of longitude is cos(latitude) times one degree of latitude
      lat_max <- max(abs(xy[, 2]), na.rm = TRUE)
      shrink <- 1 - cos(lat_max * pi / 180)
      notes <- c(notes, if (shrink <= 0.02) {
        sprintf("The coordinates look like longitude and latitude in degrees. Near the equator this is fine: up to latitude %.0f, a degree of longitude is at most %.1f%% shorter than a degree of latitude.",
                lat_max, 100 * shrink)
      } else {
        sprintf("The coordinates look like longitude and latitude in degrees: at latitude %.0f, a degree of longitude is %.0f%% shorter than a degree of latitude, so east-west distances would be overstated. Project them (for example to UTM, in km) first.",
                lat_max, 100 * shrink)
      })
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

# --- Figures of the app: drawn on screen and saved for publication -------------------
# Every figure is a function of `cex` (text size) and `pt` (point size), so the screen
# and the downloaded file are drawn by the same code.

# Save a figure in a file: PNG or TIFF at `dpi`, or PDF or SVG (vector), with the size in cm.
.save_figure <- function(file, draw, format = c("png", "pdf", "svg", "tiff"), width_cm = 16, height_cm = 10, dpi = 300) {
  format <- match.arg(format)
  w <- width_cm / 2.54; h <- height_cm / 2.54
  switch(format,
    png = grDevices::png(file, width = w, height = h, units = "in", res = dpi),
    tiff = grDevices::tiff(file, width = w, height = h, units = "in", res = dpi, compression = "lzw"),
    pdf = grDevices::pdf(file, width = w, height = h),
    svg = grDevices::svg(file, width = w, height = h))
  on.exit(grDevices::dev.off())
  draw()
  invisible(file)
}

# A plot with its download button underneath.
.figure_ui <- function(id, height) {
  shiny::tagList(shiny::plotOutput(id, height = height), shiny::downloadButton(paste0("dl_", id), "Download figure"))
}

# Names of the two coordinate columns, for the axes.
.axis_names <- function(m) { n <- colnames(m); if (is.null(n)) c("x", "y") else n }

# Colours of the classes (colour-blind safe, Okabe-Ito); noise (class 0) in grey.
.class_colours <- function(cl) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  ifelse(cl == 0, "grey80", pal[(cl - 1) %% length(pal) + 1])
}

# The technique of step 2 on the raw and on the corrected data, for the preview of step
# 5: the scores of the observations on the first two PCA axes (on the axes of the
# corrected R for the corrected data) or the MDS configuration, or the classes of
# k-means, SOM (units grouped into k) or HDBSCAN (0 = noise). With more than `max_n`
# observations a seeded subsample is shown. `rows` are the rows of X shown.
.technique_preview <- function(technique = c("pca", "mds", "kmeans", "som", "hdbscan"), X, Xc, Rc = NULL, k = 3,
                               seed = 1L, max_n = 2000) {
  technique <- match.arg(technique)
  rows <- utils::tail(seq_len(nrow(X)), nrow(Xc))          # the whitening drops the first times
  X <- as.matrix(X)[rows, , drop = FALSE]
  n <- nrow(X)
  sub <- if (n > max_n) sort(.with_stream(seed, .streams[["subsample"]], sample(n, max_n))) else seq_len(n)
  one <- function(M, R) {
    Z <- scale(M); Z[!is.finite(Z)] <- 0
    Z <- Z[sub, , drop = FALSE]
    switch(technique,
      pca = {
        e <- eigen(if (is.null(R)) stats::cor(M) else R, symmetric = TRUE)
        list(xy = Z %*% e$vectors[, 1:2], share = e$values[1:2] / sum(e$values), loadings = e$vectors[, 1:2, drop = FALSE])
      },
      mds = list(xy = stats::cmdscale(stats::dist(Z), k = 2)),
      kmeans = list(cl = stats::kmeans(Z, k, nstart = 10, iter.max = 50)$cluster),
      som = {
        sm <- kohonen::som(Z, grid = kohonen::somgrid(5, 5, "hexagonal"), rlen = 100)
        list(cl = stats::cutree(stats::hclust(stats::dist(sm$codes[[1]]), method = "ward.D2"), k)[sm$unit.classif])
      },
      hdbscan = list(cl = dbscan::hdbscan(Z, minPts = max(10L, as.integer(ceiling(0.02 * nrow(Z)))))$cluster))
  }
  out <- .with_stream(seed, .streams[["downstream"]], list(raw = one(X, NULL), corrected = one(as.matrix(Xc), Rc)))
  c(out, list(technique = technique, rows = rows[sub], variables = colnames(X)))
}

# Draw a preview: two panels, raw and corrected. Scores are coloured by position (the
# first coordinate, or the time), so a technique that only redraws the map shows a
# gradient; classes are drawn where the observations are (a map, the time line, or the
# most frequent class of each site in space-time).
.draw_preview <- function(p, family, coords = NULL, time = NULL, index = NULL, cex = 1, pt = 1) {
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1.5), cex = cex)
  nm <- c(pca = "PCA", mds = "MDS", kmeans = "k-means", som = "SOM", hdbscan = "HDBSCAN")[[p$technique]]
  pos1 <- switch(family, spatial = as.matrix(coords)[p$rows, 1], temporal = time[p$rows], spatiotemporal = index$time[p$rows])
  pos_name <- switch(family, spatial = "first coordinate", temporal = "time", spatiotemporal = "time")
  for (k in 1:2) {
    r <- p[[c("raw", "corrected")[k]]]
    title <- sprintf("%s, %s", nm, c("raw", "corrected")[k])
    if (!is.null(r$xy)) {
      col <- grDevices::hcl.colors(50, "viridis")[cut(rank(pos1, ties.method = "first"), 50, labels = FALSE)]
      lab <- if (p$technique == "pca") sprintf("axis %d (%.0f%%)", 1:2, 100 * r$share) else c("MDS 1", "MDS 2")
      graphics::plot(r$xy, pch = 19, cex = 0.6 * pt, col = col, xlab = lab[1], ylab = lab[2], main = title)
      if (p$technique == "pca") {
        s <- 0.9 * max(abs(r$xy)) / max(abs(r$loadings))
        graphics::arrows(0, 0, s * r$loadings[, 1], s * r$loadings[, 2], length = 0.06, col = "grey30")
        graphics::text(1.08 * s * r$loadings[, 1], 1.08 * s * r$loadings[, 2], p$variables, col = "grey20", cex = 0.85, xpd = NA)
      }
      graphics::mtext(sprintf("colour: %s (dark = low, yellow = high)", pos_name), side = 3, line = 0.25, cex = 0.7 * cex)
    } else {
      col <- .class_colours(r$cl)
      if (family == "spatial") {
        xy <- as.matrix(coords)[p$rows, , drop = FALSE]
        graphics::plot(xy, pch = 19, cex = 0.8 * pt, col = col, asp = 1, xlab = .axis_names(coords)[1], ylab = .axis_names(coords)[2], main = title)
      } else if (family == "temporal") {
        graphics::plot(time[p$rows], r$cl, pch = 15, cex = 0.9 * pt, col = col, xlab = "time", ylab = "class", main = title, yaxt = "n")
        graphics::axis(2, at = sort(unique(r$cl)))
      } else {
        site <- index$site[p$rows]
        modal <- tapply(r$cl, site, function(z) as.integer(names(which.max(table(z)))))
        xy <- as.matrix(coords)[as.integer(names(modal)), , drop = FALSE]
        graphics::plot(xy, pch = 19, cex = 0.8 * pt, col = .class_colours(modal), asp = 1, xlab = "x", ylab = "y",
                       main = paste(title, "(most frequent class)"))
      }
      cls <- sort(unique(r$cl))
      graphics::legend("topright", legend = ifelse(cls == 0, "noise", paste("class", cls)), pch = 19, col = .class_colours(cls),
                       bty = "n", cex = 0.75)
    }
  }
}

# Autocorrelation of each variable before (grey) and after (blue), joined by a dotted line.
.draw_autocor <- function(vv, cex = 1, pt = 1) {
  o <- seq_len(nrow(vv))
  graphics::par(mar = c(4.2, 8, 1, 1), las = 1, cex = cex)
  graphics::plot(vv$autocor_before, o, xlim = range(c(vv$autocor_before, vv$autocor_after, 0), na.rm = TRUE), yaxt = "n", type = "n",
                 xlab = "autocorrelation (grey: before; blue: after)", ylab = "")
  graphics::segments(vv$autocor_before, o, vv$autocor_after, o, lty = 3, col = "grey40")
  graphics::points(vv$autocor_before, o, pch = 19, cex = pt, col = "grey55")
  graphics::points(vv$autocor_after, o, pch = 19, cex = pt, col = "#2a6fb0")
  graphics::axis(2, at = o, labels = vv$variable, cex.axis = 0.85)
  graphics::abline(v = 0, lty = 3)
}

# Local structure (LISA quadrants) of one variable before and after the correction.
.draw_lisa <- function(v, variable, coords, cex = 1, pt = 1) {
  pos <- as.matrix(coords)
  cols <- c(none = "grey85", high_high = "#c0392b", low_low = "#2a6fb0", high_low = "#f39c12", low_high = "#8e44ad")
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), cex = cex)
  for (k in 1:2) {
    l <- if (k == 1) v$lisa_before else v$lisa
    l <- l[l$variable == variable, ]
    graphics::plot(pos[l$site, 1], pos[l$site, 2], pch = 19, cex = 0.7 * pt, col = cols[l$cluster], asp = 1,
                   xlab = .axis_names(pos)[1], ylab = .axis_names(pos)[2],
                   main = sprintf("%s, %s (%.0f%% of the sites)", variable, c("before", "after")[k], 100 * mean(l$cluster != "none")))
  }
  graphics::legend("topright", names(cols), pch = 19, col = cols, bty = "n", cex = 0.75)
}

# The plots of the diagnosis of one variable: variogram (space), ACF and periodogram
# (time), site means and the ACF of the mean series (space-time).
.draw_diagnosis <- function(d, variable, cex = 1, pt = 1) {
  z <- d$X[, variable]
  if (d$family == "spatial") {
    graphics::par(cex = cex, mar = c(4.2, 4.2, 2.5, 1))
    v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(z = z, x = d$coords[, 1], y = d$coords[, 2]))
    graphics::plot(v$dist, v$gamma, pch = 19, cex = pt, xlab = "distance", ylab = "semivariance", main = paste("Variogram of", variable))
  } else if (d$family == "temporal") {
    o <- order(d$time)
    graphics::par(mfrow = c(1, 2), cex = cex, mar = c(4.2, 4.2, 2.5, 1))
    stats::acf(z[o], main = paste("ACF of", variable))
    sp <- stats::spec.pgram(z[o], taper = 0, detrend = TRUE, plot = FALSE)
    graphics::plot(1 / sp$freq, sp$spec, type = "h", log = "x", xlab = "period", ylab = "power", main = "Periodogram")
  } else {
    graphics::par(mfrow = c(1, 2), cex = cex, mar = c(4.2, 4.2, 2.5, 1))
    m <- tapply(z, d$index$site, mean)
    graphics::plot(d$coords, pch = 19, cex = pt, col = grDevices::hcl.colors(20)[cut(m, 20)], asp = 1, main = "Site means", xlab = "x", ylab = "y")
    stats::acf(tapply(z, d$index$time, mean), main = "ACF of the mean series")
  }
}

# A table for the screen: numbers with 3 significant digits, TRUE/FALSE as yes/no, and
# missing values (a stability not computed, a test that does not apply) left blank.
.format_table <- function(d) {
  d[] <- lapply(d, function(x) {
    out <- if (is.logical(x)) ifelse(x, "yes", "no") else if (is.numeric(x)) format(x, digits = 3, trim = TRUE, drop0trailing = TRUE) else as.character(x)
    out[is.na(x)] <- ""
    out
  })
  d
}

# Columns whose names look like coordinates (x, y, lon, lat, easting, northing, with an
# optional unit such as x_km or lon_deg): the app guesses x and y from them and leaves
# them out of the variables it proposes.
.coord_like <- function(nm, axis = c("any", "x", "y")) {
  axis <- match.arg(axis)
  stem <- switch(axis, x = "x|lon|long|longitude|easting", y = "y|lat|latitude|northing",
                 any = "x|y|lon|long|longitude|lat|latitude|easting|northing")
  grepl(sprintf("^(%s)([_.]?(km|m|deg))?$", stem), nm, ignore.case = TRUE)
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
    "write.csv(diagnosis_table(sig, grouping), \"diagnosis.csv\", row.names = FALSE)",
    sprintf("report_correction(cr, X, v, file = \"depprep_report.html\", signature = sig, grouping = grouping, %s)", pos))
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
        shiny::p(class = "help", "Load a CSV with one row per observation. Mark the columns of the positions (coordinates in km or m, or longitude and latitude in degrees); the other numeric columns are the variables."),
        shiny::fileInput("file", "CSV file", accept = ".csv"),
        shiny::actionButton("example", "Use example data (spatial simulation)"),
        shiny::radioButtons("family", "Data", c("Spatial" = "spatial", "Temporal" = "temporal", "Space-time" = "spatiotemporal"), inline = TRUE),
        shiny::uiOutput("cols"),
        shiny::uiOutput("data_msg"),
        shiny::tableOutput("head")),
      shiny::tabPanel("2. Context",
        shiny::radioButtons("technique", "The corrected data are for",
          c("PCA (reads the corrected correlation matrix R)" = "pca", "MDS" = "mds", "k-means" = "kmeans", "SOM" = "som",
            "HDBSCAN" = "hdbscan")),
        shiny::p(class = "help", "MDS, k-means, SOM and HDBSCAN read the corrected data matrix X. Step 5 previews the technique chosen here, before and after the correction."),
        shiny::radioButtons("regions", "Could the classes you look for be regions or periods (contiguous in space or time)?",
          c("No" = "no", "Yes" = "yes", "Not sure" = "unsure")),
        shiny::uiOutput("context_msg")),
      shiny::tabPanel("3. Diagnosis",
        shiny::p(class = "help", "Each variable is diagnosed alone (variogram or ACF), never through its correlations with the others. Variables with similar signatures form a group that shares one basis."),
        shiny::actionButton("diagnose", "Diagnose"),
        shiny::div(style = "overflow-x:auto", shiny::tableOutput("signatures")),
        shiny::downloadButton("dl_diag", "Diagnosis (CSV)"),
        shiny::fluidRow(shiny::column(4, shiny::uiOutput("plot_var"), shiny::p(class = "help", "Text size, size and format of the figures: Figures, in step 5.")),
                        shiny::column(8, .figure_ui("diag_plot", 300))),
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
        shiny::wellPanel(shiny::strong("Figures (on screen and in the downloads)"),
          shiny::fluidRow(
            shiny::column(2, shiny::numericInput("fig_cex", "Text size", 1, min = 0.5, max = 3, step = 0.1)),
            shiny::column(2, shiny::numericInput("fig_pt", "Point size", 1, min = 0.2, max = 4, step = 0.1)),
            shiny::column(2, shiny::numericInput("fig_w", "Width (cm)", 17, min = 5, max = 60)),
            shiny::column(2, shiny::numericInput("fig_h", "Height (cm)", 9, min = 5, max = 60)),
            shiny::column(2, shiny::selectInput("fig_fmt", "Format", c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"))),
            shiny::column(2, shiny::numericInput("fig_dpi", "Resolution (dpi, PNG and TIFF)", 300, min = 72, max = 1200, step = 50)))),
        shiny::actionButton("verify", "Verify"),
        shiny::uiOutput("q3_msg"),
        shiny::h4("Each variable before and after the correction"),
        shiny::p(class = "help", "Variance removed, autocorrelation, and the tests of structure on the raw data and on the corrected data: the Moran correlogram (Bonferroni over distance classes) and the share of sites with local structure (LISA, FDR) in space, the Ljung-Box test in time."),
        shiny::div(style = "overflow-x:auto", shiny::tableOutput("before_after")),
        shiny::h4("Autocorrelation of each variable before and after"),
        .figure_ui("autocor_plot", 320),
        shiny::h4("How the correlations changed (Delta R)"),
        shiny::tableOutput("delta"),
        shiny::uiOutput("lisa_ui"),
        shiny::h4("Preview of the technique, before and after"),
        shiny::uiOutput("preview_help"),
        shiny::numericInput("k", "Number of groups (k-means, SOM)", 3, min = 2, max = 20),
        .figure_ui("preview", 420)),
      shiny::tabPanel("6. Export",
        shiny::p(class = "help", "Download the results and the R code that reproduces them."),
        shiny::downloadButton("dl_x", "Corrected data (CSV)"),
        shiny::downloadButton("dl_r", "Corrected correlation matrix (CSV)"),
        shiny::downloadButton("dl_report", "Report (HTML)"),
        shiny::downloadButton("dl_code", "R code"),
        shiny::downloadButton("dl_diag6", "Diagnosis before the correction (CSV)"),
        shiny::downloadButton("dl_lisa", "Local Moran's I before and after (CSV, space)"),
        shiny::hr(),
        shiny::p(class = "help", "Limits: missing values are not supported yet; project coordinates in degrees beyond about 11 degrees of latitude; filters such as SAR/CAR or STAR are not in the menu."))))
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
    guess_xy <- function(axis) { m <- nm[.coord_like(nm, axis)]; if (length(m)) m[1] else NULL }
    # keep what the user chose; this block is drawn again only for new data or a new family
    keep <- function(id, fallback) { v <- shiny::isolate(input[[id]]); if (!is.null(v) && v %in% nm) v else fallback }
    sel <- list(x = keep("x", guess_xy("x")), y = keep("y", guess_xy("y")),
                time = keep("time", guess("^(t|time|date|month|year)$")), site = keep("site", guess("^(site|station|pixel|id)$")))
    used <- unlist(c(if (input$family != "temporal") sel[c("x", "y")], if (input$family != "spatial") sel["time"],
                     if (input$family == "spatiotemporal") sel["site"]))
    shiny::tagList(
      if (input$family != "temporal") shiny::fluidRow(shiny::column(3, shiny::selectInput("x", "x column (km, m or longitude)", nm, sel$x)),
                                                   shiny::column(3, shiny::selectInput("y", "y column (km, m or latitude)", nm, sel$y))),
      if (input$family != "spatial") shiny::fluidRow(shiny::column(3, shiny::selectInput("time", "time column", nm, sel$time)),
        if (input$family == "spatiotemporal") shiny::column(3, shiny::selectInput("site", "site column", nm, sel$site))),
      shiny::selectizeInput("variables", "Variables", num_cols(), setdiff(num_cols()[!.coord_like(num_cols())], used), multiple = TRUE))
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
  # a figure on screen and its download, drawn by the same function of text and point size
  figure <- function(id, name, draw) {
    output[[id]] <- shiny::renderPlot(draw(input$fig_cex, input$fig_pt))
    output[[paste0("dl_", id)]] <- shiny::downloadHandler(
      function() paste0(name(), ".", input$fig_fmt),
      function(f) .save_figure(f, function() draw(input$fig_cex, input$fig_pt), input$fig_fmt, input$fig_w, input$fig_h, input$fig_dpi))
  }
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
    .format_table(diagnosis_table(rv$sig, rv$grouping))
  })
  output$plot_var <- shiny::renderUI({ shiny::req(input$variables); shiny::selectInput("pvar", "Plot of", input$variables) })
  figure("diag_plot", function() paste0("diagnosis_", input$pvar), function(cex, pt) {
    d <- dat(); shiny::req(input$pvar %in% colnames(d$X))
    .draw_diagnosis(d, input$pvar, cex, pt)
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
  output$before_after <- shiny::renderTable({
    shiny::req(rv$v)
    tab <- .before_after(rv$v, dat()$family)
    names(tab) <- attr(tab, "labels")
    .format_table(tab)
  })
  output$lisa_ui <- shiny::renderUI({
    shiny::req(rv$v, !is.null(rv$v$lisa))
    shiny::tagList(shiny::h4("Local structure before and after (LISA)"),
                   shiny::selectInput("lvar", "Variable", rv$v$variables$variable,
                                      rv$v$variables$variable[which.max(rv$v$variables$lisa_share_before)]),
                   .figure_ui("lisa_plot", 360))
  })
  figure("lisa_plot", function() paste0("lisa_", input$lvar), function(cex, pt) {
    shiny::req(rv$v, !is.null(rv$v$lisa), input$lvar)
    .draw_lisa(rv$v, input$lvar, dat()$coords, cex, pt)
  })
  figure("autocor_plot", function() "autocorrelation_before_after", function(cex, pt) {
    shiny::req(rv$v)
    .draw_autocor(rv$v$variables, cex, pt)
  })
  output$preview_help <- shiny::renderUI(shiny::p(class = "help", switch(input$technique,
    pca = "PCA scores of the observations on the first two axes, raw (axes of the raw correlations) and corrected (axes of the corrected R), coloured by position; the arrows are the variables. A gradient of colour along an axis means the axis follows the position.",
    mds = "MDS of the observations (Euclidean distances of the standardised data), raw and corrected, coloured by position. A gradient of colour means the configuration redraws the map or the time line.",
    kmeans = "k-means classes of the observations, raw and corrected, drawn where the observations are.",
    som = "SOM (5 x 5 units, grouped into k classes), raw and corrected, drawn where the observations are.",
    hdbscan = "HDBSCAN classes of the observations (grey: noise), raw and corrected, drawn where the observations are.")))
  preview <- shiny::reactive({
    shiny::req(rv$cr, !inherits(rv$cr, "error"))
    .technique_preview(input$technique, dat()$X, rv$cr$X, rv$cr$R, k = input$k)
  })
  figure("preview", function() paste0("preview_", input$technique), function(cex, pt) {
    d <- dat()
    .draw_preview(preview(), d$family, d$coords, d$time, d$index, cex, pt)
  })
  output$dl_x <- shiny::downloadHandler("corrected_data.csv", function(f) utils::write.csv(rv$cr$X, f, row.names = FALSE))
  output$dl_r <- shiny::downloadHandler("corrected_R.csv", function(f) utils::write.csv(rv$cr$R, f))
  output$dl_report <- shiny::downloadHandler("depprep_report.html", function(f) {
    d <- dat()
    do.call(report_correction, c(list(rv$cr, d$X, rv$v, file = f, code = code_text(), signature = rv$sig, grouping = rv$grouping),
                                 pos_args(d)))
  })
  dl_diagnosis <- shiny::downloadHandler("diagnosis.csv", function(f) utils::write.csv(diagnosis_table(rv$sig, rv$grouping), f, row.names = FALSE))
  output$dl_diag <- dl_diagnosis
  output$dl_diag6 <- dl_diagnosis
  output$dl_lisa <- shiny::downloadHandler("lisa_before_after.csv", function(f) {
    shiny::req(rv$v, !is.null(rv$v$lisa))
    utils::write.csv(rbind(cbind(when = "before", rv$v$lisa_before), cbind(when = "after", rv$v$lisa)), f, row.names = FALSE)
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
