#' Run the depprep app
#'
#' Opens a point-and-click interface to the whole workflow, in six steps:
#' 1. load a CSV file and mark the position columns (coordinates, time, site);
#' 2. say which technique the corrected data are for and whether the classes sought
#'    could be regions or periods;
#' 3. diagnose the dependence of each variable, with its plots (variogram, Moran
#'    correlogram and map in space; series with trends, ACF and periodogram in time), and
#'    change a regime by hand or build the menu of structures of a group
#'    ([basis_menu()]), each structure explained;
#' 4. correct, with the default explained (the data matrix with one basis per variable,
#'    the correlation matrix with the shared basis) and the alternatives;
#' 5. verify what the correction removed and kept, with a preview of the technique
#'    before and after;
#' 6. download the corrected data, the correlation matrix, the report
#'    ([report_correction()]) and the R code that reproduces everything.
#'
#' Coordinates are longitude (x) and latitude (y) in decimal degrees. Degrees treat a degree of longitude as a degree of latitude, which holds near the
#' equator; the app says how much shorter a degree of longitude is at the largest latitude
#' of the data, and asks to project them when that exceeds 2% (beyond about 11 degrees of
#' latitude).
#'
#' Missing values (empty cells) are accepted (decision D26): each variable is diagnosed,
#' corrected and correlated on its observed rows, and the previews of the techniques on
#' observations use the rows where every variable has a value. Large data (over a few
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
# the positions of the family, and the checks the app shows (missing positions stop;
# missing values, coordinates that look like degrees and large data give notes).
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
  .check_values(X, if (family == "spatial") out$coords else NULL, out$time)
  if (anyNA(X)) {
    notes <- c(notes, sprintf(paste("%.1f%% of the values are missing (%s): each variable uses its observed rows, and the",
                                    "previews of the techniques on observations use the %d rows where every variable has a value."),
                              100 * mean(is.na(X)), paste(colnames(X)[colSums(is.na(X)) > 0], collapse = ", "),
                              sum(stats::complete.cases(X))))
  }
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
# observations a seeded subsample is shown. `rows` are the rows of X shown. With missing
# values only the rows where every variable has a value, raw and corrected, are shown
# (the raw correlations of the PCA are then pairwise, as in correct()).
.technique_preview <- function(technique = c("pca", "mds", "kmeans", "som", "hdbscan"), X, Xc, Rc = NULL, k = 3,
                               seed = 1L, max_n = 2000) {
  technique <- match.arg(technique)
  rows <- utils::tail(seq_len(nrow(X)), nrow(Xc))          # the whitening drops the first times
  X <- as.matrix(X)[rows, , drop = FALSE]
  Xc <- as.matrix(Xc)
  whole <- stats::complete.cases(X) & stats::complete.cases(Xc)
  if (!all(whole)) {
    if (sum(whole) < max(k + 1, 10)) stop("Too few rows where every variable has a value for this preview.", call. = FALSE)
    rows <- rows[whole]
  }
  n <- sum(whole)
  sub <- if (n > max_n) sort(.with_stream(seed, .streams[["subsample"]], sample(n, max_n))) else seq_len(n)
  one <- function(M, R) {
    C <- if (anyNA(M)) stats::cor(M, use = "pairwise.complete.obs") else stats::cor(M)
    M <- M[whole, , drop = FALSE]
    Z <- scale(M); Z[!is.finite(Z)] <- 0
    Z <- Z[sub, , drop = FALSE]
    switch(technique,
      pca = {
        e <- eigen(if (is.null(R)) C else R, symmetric = TRUE)
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
  out <- .with_stream(seed, .streams[["downstream"]], list(raw = one(X, NULL), corrected = one(Xc, Rc)))
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

# The plots of the diagnosis of one variable, by family: the variogram, the Moran
# correlogram and the map of the values in space; the series with its trends and the ACF
# and periodogram in time; in space-time, the same on the site means and on the series
# of the sites.
.diagnosis_plots <- function(family) switch(family,
  spatial = c("Variogram" = "variogram", "Moran correlogram" = "correlogram", "Map of the values" = "map"),
  temporal = c("Series and trends" = "series", "ACF and periodogram" = "acf"),
  spatiotemporal = c("Map of the site means" = "map", "Moran correlogram of the site means" = "correlogram",
                     "Series of the sites and trends" = "series", "ACF of the mean series" = "acf"))

.draw_diagnosis <- function(d, variable, cex = 1, pt = 1, what = .diagnosis_plots(d$family)[[1]]) {
  what <- match.arg(what, unname(.diagnosis_plots(d$family)))
  z <- d$X[, variable]
  st <- d$family == "spatiotemporal"
  site_z <- if (st) as.vector(tapply(z, d$index$site, mean, na.rm = TRUE)) else z
  site_z[is.nan(site_z)] <- NA                               # a site never observed
  o <- !is.na(z)                                           # missing values are left out of the plots
  lab <- .axis_names(d$coords)
  graphics::par(cex = cex, mar = c(4.2, 4.2, 2.5, 1))
  switch(what,
    variogram = {
      v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(z = z, x = d$coords[, 1], y = d$coords[, 2])[o, ])
      graphics::plot(v$dist, v$gamma, pch = 19, cex = pt, xlab = "distance", ylab = "semivariance", main = paste("Variogram of", variable))
    },
    correlogram = {
      cg <- .moran_correlogram(matrix(site_z, ncol = 1, dimnames = list(NULL, variable)), d$coords,
                               .fixed_space(matrix(0, length(site_z), 0)))
      sig <- !is.na(cg$p_progressive) & cg$p_progressive <= 0.05
      graphics::plot(cg$distance, cg$moran_I, type = "b", pch = ifelse(sig, 19, 1), cex = pt, col = "#2a6fb0",
                     ylim = range(c(cg$moran_I, cg$expected, 0), na.rm = TRUE), xlab = "distance (centre of the class)",
                     ylab = "Moran's I", main = paste0("Moran correlogram of ", variable, if (st) " (site means)"))
      graphics::abline(h = cg$expected[1], lty = 3)
      graphics::legend("topright", c("significant (5%, progressive Bonferroni)", "not significant", "expected without dependence"),
                       pch = c(19, 1, NA), lty = c(NA, NA, 3), col = c("#2a6fb0", "#2a6fb0", "black"), bty = "n", cex = 0.8)
    },
    map = {
      graphics::par(mar = c(4.2, 4.2, 2.5, 7))
      pal <- grDevices::hcl.colors(20)
      graphics::plot(d$coords, pch = 19, cex = pt, col = pal[cut(site_z, 20)], asp = 1, xlab = lab[1], ylab = lab[2],
                     main = paste0(variable, if (st) " (site means)"))
      at <- seq(min(site_z, na.rm = TRUE), max(site_z, na.rm = TRUE), length.out = 5)
      graphics::legend("topleft", inset = c(1.02, 0), legend = format(signif(at, 3)), pch = 15, pt.cex = 1.6,
                       col = pal[c(1, 5, 10, 15, 20)], bty = "n", xpd = NA, cex = 0.8)
    },
    series = {
      if (st) {
        tm <- sort(unique(d$index$time))
        M <- tapply(z, list(d$index$site, d$index$time), mean)
        shown <- unique(round(seq(1, nrow(M), length.out = min(60, nrow(M)))))
        y <- colMeans(M, na.rm = TRUE)
        graphics::matplot(tm, t(M[shown, , drop = FALSE]), type = "l", lty = 1, col = grDevices::adjustcolor("grey50", 0.35),
                          xlab = "time", ylab = variable, main = sprintf("%s: %d of %d sites (grey) and their mean (black)", variable, length(shown), nrow(M)))
        graphics::lines(tm, y, lwd = 2)
      } else {
        o <- order(d$time); tm <- d$time[o]; y <- z[o]
        graphics::plot(tm, y, type = "l", col = "grey40", xlab = "time", ylab = variable, main = paste("Series of", variable))
      }
      lin <- stats::lm(y ~ tm, na.action = stats::na.exclude)
      cub <- stats::lm(y ~ stats::poly(tm, min(3, length(unique(tm)) - 1)), na.action = stats::na.exclude)
      graphics::lines(tm, stats::fitted(lin), lty = 2, lwd = 2, col = "#c0392b")
      graphics::lines(tm, stats::fitted(cub), lwd = 2, col = "#2a6fb0")
      r2 <- c(summary(lin)$r.squared, summary(cub)$r.squared)
      graphics::legend("topleft", c(sprintf("linear trend, R2 = %.2f (a trend regime above %.2f)", r2[1], default_thresholds()$trend_r2),
                                    sprintf("cubic trend (the menu's trend), R2 = %.2f", r2[2])),
                       lty = c(2, 1), lwd = 2, col = c("#c0392b", "#2a6fb0"), bty = "n", cex = 0.8)
    },
    acf = if (st) {
      m <- tapply(z, d$index$time, mean, na.rm = TRUE); m[is.nan(m)] <- NA
      stats::acf(m, main = paste("ACF of the mean series of", variable), na.action = stats::na.pass)
    } else {
      o <- order(d$time)
      graphics::par(mfrow = c(1, 2))
      stats::acf(z[o], main = paste("ACF of", variable), na.action = stats::na.pass)
      if (anyNA(z)) {                                      # the FFT periodogram needs a complete series
        graphics::plot.new()
        graphics::title("Periodogram")
        graphics::text(0.5, 0.5, "The series has gaps:\nsee the period and its strength\nin the table of the diagnosis.", cex = 0.85)
      } else {
        sp <- stats::spec.pgram(z[o], taper = 0, detrend = TRUE, plot = FALSE)
        graphics::plot(1 / sp$freq, sp$spec, type = "h", log = "x", xlab = "period", ylab = "power", main = "Periodogram")
      }
    })
  invisible(NULL)
}

# How to read each plot of the diagnosis.
.diagnosis_plot_help <- c(
  variogram = "Semivariance (half the mean squared difference) of pairs of sites by distance: it rises while near sites are more alike than distant ones, and levels off at the range of the dependence.",
  correlogram = "Moran's I by distance class: positive when sites at that distance have similar values, near the dotted line when they are unrelated. Filled points are significant (5%, Bonferroni over the classes in order). Step 5 shows the same test after the correction.",
  map = "Where the high (yellow) and low (dark) values are. A smooth gradient or large patches point to a large-scale structure; a salt-and-pepper map to little spatial dependence.",
  series = "The series over time (in space-time, a sample of sites in grey and their mean in black) with two trends: the linear one (red, dashed), whose R2 classifies a variable as a trend regime, and the cubic one (blue), the trend the menu below would remove.",
  acf = "The ACF is the correlation of the series with itself k steps later: bars beyond the dashed band are significant, and a slow decay means persistence. The periodogram shows the power of each period: a clear peak marks a cycle.")

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

# The neighbourhood graphs of the menu, with how each links the sites.
.graph_choices <- c("Distance-based (dbMEM, the default)" = "dbmem", "Best of the three graphs below, chosen from the data" = "auto",
                    "Delaunay triangulation" = "del", "Gabriel graph" = "gab", "Relative neighbourhood graph" = "rel")
.graph_help <- c(
  dbmem = "Sites closer than a threshold are neighbours. The threshold is the longest link needed to keep every site connected (the longest edge of the minimum spanning tree). This is the graph of the automatic choice.",
  auto = "Tries the Delaunay, Gabriel and relative-neighbourhood graphs, each with equal or distance-decreasing weights, and keeps the one whose maps explain the group best (adjusted R2); the test is made stricter for the six tries (Bauman et al. 2018). The figure shows the Delaunay graph.",
  del = "Each site is linked to its natural neighbours (the triangles with no other site inside their circle): many links, some long ones at the edges of the area.",
  gab = "The Delaunay links whose circle with the link as diameter holds no other site: fewer and shorter links.",
  rel = "Two sites are linked only when no third site is closer to both: the sparsest of the three, mostly the nearest neighbours.")

# The links (pairs of sites, i < j) of a neighbourhood graph of the menu.
.menu_graph_links <- function(coords, W = "dbmem") {
  coords <- as.matrix(coords)
  g <- sub("_.*$", "", W)
  nb <- switch(if (g == "auto") "del" else g,
    dbmem = spdep::dnearneigh(coords, 0, adespatial::give.thresh(stats::dist(coords)) * (1 + 1e-9)),
    del = spdep::tri2nb(coords),
    gab = spdep::graph2nb(spdep::gabrielneigh(coords), sym = TRUE),
    rel = spdep::graph2nb(spdep::relativeneigh(coords), sym = TRUE),
    stop("Unknown graph: ", W, call. = FALSE))
  e <- do.call(rbind, lapply(seq_along(nb), function(i) { j <- nb[[i]][nb[[i]] > i]; if (length(j)) cbind(i, j) }))
  if (is.null(e)) matrix(integer(0), 0, 2) else unname(e)
}

# The sites and the links of the graph chosen for the MEM (the Delaunay graph when only
# the directional maps are ticked), with an arrow the way a directional process flows.
.draw_menu_graph <- function(coords, W = "dbmem", direction = NULL, cex = 1, pt = 1, max_sites = 3000) {
  coords <- as.matrix(coords); lab <- .axis_names(coords)
  graphics::par(cex = cex, mar = c(4.2, 4.2, 2.5, 1))
  draw_links <- nrow(coords) <= max_sites
  title <- if (!draw_links) sprintf("%d sites (too many to draw the links)", nrow(coords)) else
    paste("Links:", names(.graph_choices)[match(if (W == "auto") "del" else W, .graph_choices)])
  graphics::plot(coords, type = "n", asp = 1, xlab = lab[1], ylab = lab[2], main = title, cex.main = 0.9)
  if (draw_links) {
    e <- .menu_graph_links(coords, W)
    graphics::segments(coords[e[, 1], 1], coords[e[, 1], 2], coords[e[, 2], 1], coords[e[, 2], 2], col = "grey65", lwd = 0.7)
  }
  graphics::points(coords, pch = 19, cex = 0.5 * pt, col = "grey20")
  if (!is.null(direction) && is.finite(direction)) {
    th <- direction * pi / 180
    ctr <- colMeans(apply(coords, 2, range))
    len <- 0.2 * sqrt(sum(apply(coords, 2, function(v) diff(range(v)))^2))
    graphics::arrows(ctr[1] - len * cos(th), ctr[2] - len * sin(th), ctr[1] + len * cos(th), ctr[2] + len * sin(th),
                     lwd = 4, length = 0.15, col = "#c0392b")
    graphics::mtext(sprintf("red arrow: the process flows toward %g degrees", direction), side = 3, line = 0.2, cex = 0.75 * cex)
  }
  invisible(NULL)
}

# A menu in words, for the screen.
.menu_words <- function(m, family) {
  if (is.null(m)) return("automatic choice")
  one <- function(t) switch(t$kind,
    trend = "trend",
    cycle = if (is.null(t$period)) "cycle (period from the periodogram)" else sprintf("cycle of period %g", t$period),
    mem = if (identical(t$margin, "time") || (is.null(t$margin) && family == "temporal")) "temporal MEM" else
      sprintf("spatial MEM (%s)", tolower(sub(" \\(.*$", "", c(names(.graph_choices)[match(sub("_.*$", "", t$W), .graph_choices)], t$W)[1]))),
    aem = sprintf("directional maps toward %g degrees", t$direction))
  paste(vapply(m, one, ""), collapse = " + ")
}

# The R code that reproduces what was done in the app.
.app_code <- function(file, family, x, y, time, site, variables, overrides, menus, mode, scope, other = NULL) {
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
    if (!is.null(other)) sprintf("other <- correct(X, grouping, %s, mode = %s, scope = %s, data_basis = \"union\"%s%s)", pos, q(other$mode), q(other$scope),
                                 if (other$scope == scope) ", bases = cr$bases" else "", if (length(menus)) ", menu = menus" else ""),
    sprintf("v <- verify(cr, X, %s)", pos),
    "write.csv(cr$X, \"corrected_data.csv\", row.names = FALSE)",
    "write.csv(cr$R, \"corrected_R.csv\")",
    if (!is.null(other)) "write.csv(other$R, \"corrected_R_other_mode.csv\")",
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

.default_mode <- function(family) if (family == "spatiotemporal") "conditional_ar" else "conditional"

.default_scope <- function(mode, family, menus = NULL) {
  if (mode != "whiten" && family != "spatiotemporal" && !length(menus)) "global" else "groups"
}

# The mode that R is compared with (decision D30): what is being estimated cannot be answered
# before the analysis, so both run and the user reads where they differ.
.other_mode <- function(mode, family) if (mode == "whiten") .default_mode(family) else "whiten"

# R by the other mode. Only R is used, since X does not depend on the mode (D17), and the
# bases already built are reused when both modes share the scope.
.run_other_mode <- function(first, X, grouping, pos, family, menus = NULL, cores = 1L) {
  mode <- .other_mode(first$mode, family)
  scope <- .default_scope(mode, family, menus)
  do.call(correct, c(list(X, grouping, mode = mode, scope = scope, data_basis = "union", cores = cores,
                          bases = if (scope == first$scope) first$bases, menu = if (length(menus)) menus), pos))
}

# The pairs whose correlation differs by `threshold` or more between two matrices, largest first.
.mode_disagreement <- function(R1, R2, threshold = 0.2) {
  low <- which(lower.tri(R1), arr.ind = TRUE)
  d <- data.frame(a = rownames(R1)[low[, 2]], b = rownames(R1)[low[, 1]], first = R1[low], second = R2[low], stringsAsFactors = FALSE)
  d$difference <- d$second - d$first
  d <- d[abs(d$difference) >= threshold, ]
  d <- d[order(-abs(d$difference)), ]
  rownames(d) <- NULL
  attr(d, "n_pairs") <- nrow(low)
  d
}

# Predicted minutes of correct() in space-time (decision D30). The time grows with the square
# of the rows and falls slowly with the processes, as measured on the Foz and in the timing
# runs (analysis/09d, 09g, 09h and 08f). The same size can take eight times less or 1.5 times
# more, since the regimes count too: the shelf took longer than the larger slope.
.predicted_minutes <- function(rows, cores = 1L) {
  typical <- 6 * (rows / 1e4)^2 / cores^0.4
  c(low = typical / 8, high = typical * 1.5)
}

.time_words <- function(minutes) {
  if (minutes[["high"]] < 1) return("under a minute")
  round2 <- function(x) as.numeric(formatC(x, digits = 2, format = "g"))
  if (minutes[["low"]] < 1) return(sprintf("up to %g minutes", ceiling(round2(minutes[["high"]]))))
  if (minutes[["high"]] < 120) return(sprintf("about %g to %g minutes", round2(minutes[["low"]]), round2(minutes[["high"]])))
  sprintf("about %g to %g hours", round2(minutes[["low"]] / 60), round2(minutes[["high"]] / 60))
}

# The messages of steps 4 and 5: `step` is "correct" or "verify".
.time_message <- function(rows, family, mode, cores = 1L, step = c("correct", "verify")) {
  step <- match.arg(step)
  if (family != "spatiotemporal") {
    return("Expected time: seconds to a few minutes with a few thousand observations; larger data were not timed.")
  }
  m <- .predicted_minutes(rows, cores)
  if (step == "correct") {
    return(sprintf("Expected time with %d observations and %d process(es): %s. It grows with the square of the observations and depends on the regimes. Running the other mode too adds seconds to minutes, since it reuses the bases.",
                   rows, cores, .time_words(m)))
  }
  if (mode == "whiten") return("Expected time: seconds to a minute, since whitening has no projected part.")
  sprintf("Expected time: %s, because the projected part is computed pair by pair.", .time_words(m * 0.4))
}

# What steps 4 and 5 say about the comparison of the two modes.
.modes_text <- function(p, mode, other) {
  if (!nrow(p)) return(sprintf("Modes %s and %s agree: no pair differs by 0.2 or more, so the correlations do not depend on the mode.", mode, other))
  sprintf("Modes %s and %s disagree on %d of %d pairs (difference of 0.2 or more, in the table of step 4). The conditional modes keep associations that change slowly over sites and time; whitening keeps only how changes go together from one site or time step to the next. For these pairs the answer depends on that, so report both values.",
          mode, other, nrow(p), attr(p, "n_pairs"))
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
        shiny::p(class = "help", "Load a CSV with one row per observation. Mark the columns of the positions (longitude and latitude in decimal degrees, the time, the site); the other numeric columns are the variables."),
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
        shiny::h4("Plots of each variable"),
        shiny::fluidRow(shiny::column(4, shiny::uiOutput("plot_var"), shiny::uiOutput("diag_help"),
                                      shiny::p(class = "help", "Text size, size and format of the figures: Figures, in step 5.")),
                        shiny::column(8, .figure_ui("diag_plot", 360))),
        shiny::h4("Change a regime by hand"),
        shiny::uiOutput("override_ui"),
        shiny::h4("Menu of structures for a group"),
        shiny::p(class = "help", "The basis of a group holds the large-scale structures its variables share, and the correction removes them. By default they are chosen automatically from the group's regime. When you know the system (the period of a cycle, the direction of a current), tick the structures it can have: the basis is then built only from them. The data still decide how many terms enter: each term is tested by permutation and the size of the basis is capped (k/n), so a ticked structure the variables do not follow is left out. In space-time, the spatial and the temporal structures ticked are combined (their products)."),
        shiny::uiOutput("menu_ui"),
        shiny::strong("Menus set"),
        shiny::verbatimTextOutput("menus")),
      shiny::tabPanel("4. Correction",
        shiny::uiOutput("correct_help"),
        shiny::selectInput("mode", "Mode", c("Default for the data" = "default", "Conditional (residuals on the bases)" = "conditional",
                                         "Conditional + AR, space-time (residuals filtered in time)" = "conditional_ar",
                                         "Whitening (innovations)" = "whiten")),
        shiny::selectInput("scope", "Bases for R", c("Default for the data" = "default", "One per group" = "groups", "One global basis" = "global",
                                                  "Union of the variables' bases" = "variables")),
        shiny::uiOutput("cores_ui"),
        shiny::checkboxInput("both_modes", "Compare with the other mode and show where the two disagree", value = TRUE),
        shiny::uiOutput("time_correct"),
        shiny::actionButton("correct", "Correct"),
        shiny::tableOutput("bases"),
        shiny::uiOutput("correct_msg"),
        shiny::uiOutput("modes_msg"),
        shiny::tableOutput("modes_table")),
      shiny::tabPanel("5. Verification",
        shiny::wellPanel(shiny::strong("Figures (on screen and in the downloads)"),
          shiny::fluidRow(
            shiny::column(2, shiny::numericInput("fig_cex", "Text size", 1, min = 0.5, max = 3, step = 0.1)),
            shiny::column(2, shiny::numericInput("fig_pt", "Point size", 1, min = 0.2, max = 4, step = 0.1)),
            shiny::column(2, shiny::numericInput("fig_w", "Width (cm)", 17, min = 5, max = 60)),
            shiny::column(2, shiny::numericInput("fig_h", "Height (cm)", 9, min = 5, max = 60)),
            shiny::column(2, shiny::selectInput("fig_fmt", "Format", c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"))),
            shiny::column(2, shiny::numericInput("fig_dpi", "Resolution (dpi, PNG and TIFF)", 300, min = 72, max = 1200, step = 50)))),
        shiny::uiOutput("time_verify"),
        shiny::actionButton("verify", "Verify"),
        shiny::uiOutput("q3_msg"),
        shiny::uiOutput("modes_msg5"),
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
        shiny::uiOutput("dl_modes_ui"),
        shiny::hr(),
        shiny::p(class = "help", "Limits: with missing values, the previews of the techniques on observations use only the rows where every variable has a value; project coordinates in degrees beyond about 11 degrees of latitude; filters such as SAR/CAR or STAR are not in the menu."))))
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
      if (input$family != "temporal") shiny::fluidRow(shiny::column(3, shiny::selectInput("x", "x column (longitude)", nm, sel$x)),
                                                   shiny::column(3, shiny::selectInput("y", "y column (latitude)", nm, sel$y))),
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
    out <- list(shiny::p(sprintf("Step 4 estimates R with the mode %s and compares it with %s.", .default_mode(input$family), .other_mode(.default_mode(input$family), input$family))),
      shiny::p(class = "help", "The design is yours: the correction removes what is large scale for the grid and the period you loaded. Another extent can change what counts as large scale, and with it the result."),
      shiny::p(if (input$technique == "pca") "PCA reads the corrected correlation matrix R, built with the shared basis of each group." else
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
  output$plot_var <- shiny::renderUI({
    shiny::req(input$variables)
    shiny::tagList(shiny::selectInput("pvar", "Variable", input$variables, shiny::isolate(input$pvar)),
                   shiny::radioButtons("pwhat", "Plot", .diagnosis_plots(input$family)))
  })
  output$diag_help <- shiny::renderUI({ shiny::req(input$pwhat); shiny::p(class = "help", .diagnosis_plot_help[[input$pwhat]]) })
  figure("diag_plot", function() paste0("diagnosis_", input$pwhat, "_", input$pvar), function(cex, pt) {
    d <- dat(); shiny::req(input$pvar %in% colnames(d$X), input$pwhat %in% .diagnosis_plots(d$family))
    .draw_diagnosis(d, input$pvar, cex, pt, input$pwhat)
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
    sp <- input$family != "temporal"; tm <- input$family != "spatial"
    # a structure: its box, what it is, and its settings
    term <- function(id, label, help, ...) shiny::div(style = "margin-bottom:.8rem",
      shiny::checkboxInput(id, shiny::strong(label)), shiny::div(class = "help", style = "margin:-.5rem 0 .3rem 1.4rem", help),
      shiny::div(style = "margin-left:1.4rem", ...))
    trend_help <- if (sp && tm) "A smooth surface of the coordinates and a smooth curve of time: broad gradients across the area and over the whole period." else
      if (sp) "A smooth surface of the coordinates (x, y, x2, xy, y2): one broad gradient across the area, such as from the coast to the open ocean." else
      "A smooth curve of time (a cubic polynomial): a long-term rise or fall over the whole period."
    shiny::tagList(
      shiny::fluidRow(shiny::column(4, shiny::selectInput("mn_group", "Group", gs)), shiny::column(8, shiny::uiOutput("mn_group_info"))),
      shiny::fluidRow(
        shiny::column(if (sp) 6 else 12,
          term("mn_trend", "Trend", trend_help),
          if (tm) term("mn_cycle", "Cycle",
            "Sine and cosine waves of one period, such as 12 for an annual cycle in monthly data. For a regular, repeating cycle.",
            shiny::numericInput("mn_period", "Period, in time steps (empty: taken from the periodogram of the group's variables)", NA, width = "100%")),
          if (tm) term("mn_tmem", "Temporal MEM (Moran's eigenvector maps of time)",
            "Waves along time, from slow to fast, without assuming a period: for persistence and irregular fluctuations at several time scales."),
          if (sp) term("mn_mem", "Spatial MEM (Moran's eigenvector maps)",
            "Map patterns built from which sites are neighbours, from broad (a few large patches) to fine (many small patches); only the patterns the group's variables follow are kept. For patchy structure at several scales rather than one gradient.",
            shiny::selectInput("mn_W", "Neighbourhood graph: which sites count as neighbours", .graph_choices, width = "100%"),
            shiny::uiOutput("mn_graph_help")),
          if (sp) term("mn_aem", "Directional maps (AEM, asymmetric eigenvector maps)",
            "Patterns that spread one way only: upstream sites influence downstream ones, not the reverse. For a current, a river plume or wind transport. Built on the Delaunay graph, oriented by the direction.",
            shiny::numericInput("mn_dir", "Direction the process flows toward, in degrees counter-clockwise from east (0 east, 90 north, 180 west, 270 south)", 0, width = "100%"),
            shiny::uiOutput("mn_dir_help"))),
        if (sp) shiny::column(6, .figure_ui("menu_graph", 380),
          shiny::p(class = "help", "The sites and the links of the graph chosen for the spatial MEM (the Delaunay graph when only the directional maps are ticked). With the directional maps ticked, the red arrow is the way the process flows."))),
      shiny::actionButton("mn_apply", "Set the menu of this group"), shiny::actionButton("mn_clear", "Back to automatic"))
  })
  output$mn_group_info <- shiny::renderUI({
    shiny::req(rv$grouping, input$mn_group)
    g <- rv$grouping[rv$grouping$group == input$mn_group, ]
    per <- attr(rv$grouping, "period")[input$mn_group]
    auto <- tryCatch(.menu_words(.menu_of_type(g$group_type[1], if (is.null(per) || is.na(per)) NULL else unname(per)), input$family),
                     error = function(e) "the structures of its regime")
    shiny::div(class = "ok", sprintf("Variables: %s. Regime: %s. Automatic choice for this regime: %s. Set now: %s.",
                                    paste(g$variable, collapse = ", "), g$group_type[1], auto, .menu_words(rv$menus[[input$mn_group]], input$family)))
  })
  output$mn_graph_help <- shiny::renderUI({ shiny::req(input$mn_W); shiny::div(class = "help", .graph_help[[input$mn_W]]) })
  output$mn_dir_help <- shiny::renderUI({
    shiny::req(rv$sig, rv$grouping, input$mn_group)
    v <- rv$grouping$variable[rv$grouping$group == input$mn_group]
    s <- rv$sig[rv$sig$variable %in% v, , drop = FALSE]
    med <- function(x) if (is.null(x) || all(is.na(x))) NA else stats::median(x, na.rm = TRUE)
    # only the directions of a clear anisotropy or movement (the thresholds of the regimes)
    th <- default_thresholds()
    an <- if (isTRUE(med(s$aniso_ratio) >= th$aniso_ratio)) med(s$aniso_direction) else NA
    sh <- if (isTRUE(med(s$shift_indicator) >= th$shift)) med(s$shift_direction) else NA
    hints <- c(if (!is.na(an)) sprintf("the longest range of the variogram is along the axis %.0f / %.0f degrees (the anisotropy direction of the table, which is clockwise from north, in this convention)",
                                       (90 - an) %% 180, (90 - an) %% 180 + 180),
               if (!is.na(sh)) sprintf("the pattern moves toward about %.0f degrees (shift direction)", sh %% 360))
    if (length(hints)) shiny::div(class = "help", paste0("Hints from the diagnosis of this group: ", paste(hints, collapse = "; "), "."))
  })
  figure("menu_graph", function() paste0("graph_", input$mn_W), function(cex, pt) {
    d <- dat(); shiny::req(d$family != "temporal", input$mn_W)
    .draw_menu_graph(d$coords, if (isTRUE(input$mn_mem) || !isTRUE(input$mn_aem)) input$mn_W else "del",
                     if (isTRUE(input$mn_aem)) input$mn_dir, cex, pt)
  })
  shiny::observeEvent(input$mn_apply, {
    m <- tryCatch(.app_menu(trend = isTRUE(input$mn_trend), cycle = isTRUE(input$mn_cycle), period = input$mn_period, mem = isTRUE(input$mn_mem),
                            W = if (is.null(input$mn_W)) "dbmem" else input$mn_W, time_mem = isTRUE(input$mn_tmem), aem = isTRUE(input$mn_aem),
                            direction = if (is.null(input$mn_dir)) 0 else input$mn_dir), error = function(e) e)
    if (inherits(m, "error")) { shiny::showNotification(conditionMessage(m), type = "error"); return() }
    menus <- rv$menus; menus[[input$mn_group]] <- m; rv$menus <- Filter(Negate(is.null), menus)
  })
  shiny::observeEvent(input$mn_clear, { menus <- rv$menus; menus[[input$mn_group]] <- NULL; rv$menus <- menus })
  output$menus <- shiny::renderPrint({
    if (!length(rv$menus)) cat("All groups: automatic choice.\n") else
      for (g in names(rv$menus)) cat(sprintf("%s: %s\n   %s\n", g, .menu_words(rv$menus[[g]], input$family), .menu_code(rv$menus[[g]])))
  })
  settings <- shiny::reactive({
    d <- dat()
    mode <- if (input$mode == "default") .default_mode(d$family) else input$mode
    scope <- if (input$scope == "default") .default_scope(mode, d$family, rv$menus) else input$scope
    list(mode = mode, scope = scope)
  })
  # several processes pay off only on large data (each takes seconds to start on Windows)
  output$cores_ui <- shiny::renderUI({
    n <- tryCatch(nrow(dat()$X), error = function(e) 0L)
    shiny::numericInput("cores", "Processes that build the bases at the same time (the result is the same with any number; more help above a few thousand observations)",
                        if (n > 5000) min(8L, max(1L, parallel::detectCores() - 1L)) else 1L, min = 1, max = 64, step = 1)
  })
  output$correct_help <- shiny::renderUI(shiny::p(class = "help",
    "Default: in space and in time, R from one global basis of the dependent variables, conditional mode; in space-time, conditional + AR by group. The data matrix X always has each variable corrected with its own basis, which keeps the classes. With a menu, the bases are built per group."))
  time_note <- function(step, mode) {
    d <- tryCatch(dat(), error = function(e) NULL); shiny::req(d)
    cores <- if (is.numeric(input$cores) && input$cores >= 1) input$cores else 1L
    shiny::p(class = "help", .time_message(nrow(d$X), d$family, mode, cores, step))
  }
  output$time_correct <- shiny::renderUI(time_note("correct", settings()$mode))
  output$time_verify <- shiny::renderUI({ shiny::req(rv$cr, !inherits(rv$cr, "error")); time_note("verify", rv$cr$mode) })
  shiny::observeEvent(input$correct, {
    shiny::req(rv$grouping)
    d <- dat(); st <- settings()
    cores <- if (is.numeric(input$cores) && input$cores >= 1) input$cores else 1L
    shiny::withProgress(message = "Correcting", {
      rv$cr <- tryCatch(do.call(correct, c(list(d$X, rv$grouping, mode = st$mode, scope = st$scope, menu = if (length(rv$menus)) rv$menus, cores = cores),
                                           pos_args(d))), error = function(e) e)
      rv$v <- NULL; rv$other <- NULL
      if (isTRUE(input$both_modes) && !inherits(rv$cr, "error")) {
        rv$other <- tryCatch(.run_other_mode(rv$cr, d$X, rv$grouping, pos_args(d), d$family, rv$menus, cores), error = function(e) e)
      }
    })
  })
  disagreement <- shiny::reactive({
    shiny::req(rv$cr, rv$other, !inherits(rv$cr, "error"), !inherits(rv$other, "error"))
    .mode_disagreement(rv$cr$R, rv$other$R)
  })
  modes_note <- function() {
    if (inherits(rv$other, "error")) return(shiny::div(class = "note", paste("The other mode failed:", conditionMessage(rv$other))))
    p <- disagreement()
    shiny::div(class = if (nrow(p)) "note" else "ok", .modes_text(p, rv$cr$mode, rv$other$mode))
  }
  output$modes_msg <- shiny::renderUI({ shiny::req(rv$other); modes_note() })
  output$modes_msg5 <- shiny::renderUI({ shiny::req(rv$v, rv$other); modes_note() })
  output$modes_table <- shiny::renderTable({
    p <- disagreement(); shiny::req(nrow(p) > 0)
    names(p)[3:4] <- c(rv$cr$mode, rv$other$mode)
    p
  }, digits = 3)
  output$dl_modes_ui <- shiny::renderUI({ shiny::req(rv$other, !inherits(rv$other, "error")); shiny::downloadButton("dl_modes", "Pairs on which the two modes disagree (CSV)") })
  output$dl_modes <- shiny::downloadHandler("modes_disagreement.csv", function(f) utils::write.csv(disagreement(), f, row.names = FALSE))
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
              rv$overrides, rv$menus, st$mode, st$scope, if (!is.null(rv$other) && !inherits(rv$other, "error")) list(mode = rv$other$mode, scope = rv$other$scope))
  })
  output$dl_code <- shiny::downloadHandler("depprep_analysis.R", function(f) writeLines(code_text(), f))
  invisible(rv)
}
