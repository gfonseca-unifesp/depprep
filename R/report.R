#' Write a report of a correction for the end user
#'
#' Writes one self-contained HTML file (no external files, no pandoc) that says what the
#' data showed before the correction and what the correction did, in the order a user
#' checks it: the diagnosis of each variable (its dependence signature, regime, group and
#' stability), the bases of each group and the degrees of freedom they used, how much of
#' each variable was removed, the autocorrelation and the tests of structure before and
#' after (Moran correlograms and local Moran's I in space, Ljung-Box in time), how the
#' correlations between variables changed (Delta R), the warning of signal at the scale
#' of the basis (Q3), the warnings, and the R code that reproduces the correction.
#'
#' @param corrected Output of [correct()].
#' @param X The raw data (matrix or `depprep_sim`), as given to [correct()].
#' @param verification Output of [verify()]; computed (with `nperm = 99`) when `NULL`.
#' @param signature Output of [diagnose_variables()] used for the correction; computed
#'   again (without bootstrap) when `NULL`.
#' @param grouping Output of [group_regimes()] (or the grouping given to [correct()]),
#'   for the regime, group and stability of each variable; when `NULL`, the groups of
#'   `corrected`.
#' @param file Path of the HTML file.
#' @param coords,time,index As in [correct()] (not needed with a `depprep_sim`).
#' @param code R code that reproduces the correction, shown in the report; `NULL` for a
#'   template built from the settings of `corrected`.
#' @param title Title of the report.
#' @return The path of the file, invisibly.
#' @examples
#' \donttest{
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(150, g, seed = 1)
#' cr <- correct(s, group_regimes(diagnose_variables(s, nsim = 99)))
#' f <- report_correction(cr, s, file = tempfile(fileext = ".html"))
#' file.exists(f)
#' }
#' @export
report_correction <- function(corrected, X, verification = NULL, file = "depprep_report.html", coords = NULL, time = NULL,
                              index = NULL, code = NULL, title = "depprep: correction report", signature = NULL, grouping = NULL) {
  stopifnot(inherits(corrected, "depprep_corrected"))
  if (inherits(X, "depprep_sim")) {
    d <- .data_parts(X); Xm <- d$X; coords <- d$coords; time <- d$time; index <- d$index
  } else Xm <- as.matrix(X)
  if (is.null(colnames(Xm))) colnames(Xm) <- paste0("V", seq_len(ncol(Xm)))
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else if (!is.null(coords)) "spatial" else "temporal"
  v <- if (!is.null(verification)) verification else if (corrected$mode == "scale") NULL else
    verify(corrected, Xm, coords = coords, time = time, index = index, nperm = 99)
  esc <- .html_escape
  num <- function(x, d = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "")
  sec <- character(0)

  # settings
  sec <- c(sec, sprintf('<p class="meta">%s data &middot; %d observations &middot; %d variables &middot; mode <b>%s</b>, scope <b>%s</b>, data matrix <b>%s</b> &middot; depprep %s &middot; %s</p>',
                        family, nrow(Xm), ncol(Xm), esc(corrected$mode), esc(corrected$scope), esc(if (is.null(corrected$data_basis)) "" else corrected$data_basis),
                        esc(as.character(utils::packageVersion("depprep"))), format(Sys.time(), "%Y-%m-%d %H:%M")))

  # the diagnosis before the correction: what chose the regime and the group of each variable
  if (is.null(signature)) signature <- diagnose_variables(Xm, coords = coords, time = time, index = index)
  if (is.null(grouping)) grouping <- stats::setNames(corrected$group, colnames(Xm))
  dt <- diagnosis_table(signature, grouping)
  labels <- c(type = "Regime", group = "Group", stability = "Stability", ambiguous = "Ambiguous", .diagnosis_columns)
  sec <- c(sec, "<h2>Diagnosis before the correction</h2>",
           "<p>Each variable was diagnosed alone, from its own variogram, autocorrelation and periodogram, never from its correlations with the other variables. These signatures chose the regime of each variable and its group; explore them before trusting the correction. <i>dependence p</i>: the test of any dependence (Moran's I and its correlogram in space, Ljung-Box in time); <i>range</i>: the distance at which the variogram levels off; <i>nugget / sill</i>: the share of variation without spatial structure; <i>anisotropy ratio</i>: the longest over the shortest range by direction; <i>AR phi</i>: persistence from one time to the next; <i>seasonal strength</i>: the share of variance in the cycle of that period; <i>stability</i>: how often the variable stays in its group over subsamples.</p>",
           .html_df(dt, c("Variable", unname(labels[names(dt)[-1]]))))

  # bases
  grp <- corrected$group
  gs <- setdiff(unique(grp), "none")
  if (length(gs)) {
    rows <- vapply(gs, function(g) {
      b <- corrected$bases[[g]]
      tm <- attr(b, "terms")
      sprintf("<tr><td>%s</td><td>%d</td><td>%s</td><td class='n'>%d</td><td class='n'>%s</td><td>%s</td></tr>", esc(g), sum(grp == g),
              esc(if (is.null(attr(b, "type"))) "" else attr(b, "type")), ncol(b), num(ncol(b) / nrow(Xm)),
              esc(paste(c(utils::head(tm, 8), if (length(tm) > 8) sprintf("... (%d more)", length(tm) - 8)), collapse = ", ")))
    }, "")
    sec <- c(sec, "<h2>Bases</h2><p>One basis per group; every variable of a group is corrected with the same basis. k/n is the share of degrees of freedom the basis uses.</p>",
             .html_table(c("Group", "Variables", "Type", "k", "k/n", "Terms"), rows),
             sprintf("<p>Variables without dependence (not corrected): %s.</p>", if (any(grp == "none")) esc(paste(colnames(Xm)[grp == "none"], collapse = ", ")) else "none"))
  }

  if (!is.null(v)) {
    vv <- v$variables
    tab <- data.frame(variable = vv$variable, group = vv$group, removed = signif(vv$removed, 3),
                      autocor_before = signif(vv$autocor_before, 3), autocor_after = signif(vv$autocor_after, 3), stringsAsFactors = FALSE)
    head_t <- c("Variable", "Group", "Removed", "Autocor. before", "Autocor. after")
    if (family == "temporal") {
      tab$lb_before <- signif(vv$lb_p_before, 3); tab$lb_after <- signif(vv$lb_p, 3)
      head_t <- c(head_t, "Ljung-Box p before", "Ljung-Box p after")
    } else {
      tab$cp_before <- signif(vv$correlogram_p_before, 3); tab$cp_after <- signif(vv$correlogram_p, 3)
      tab$lisa_before <- signif(vv$lisa_share_before, 3); tab$lisa_after <- signif(vv$lisa_share, 3)
      head_t <- c(head_t, "Correlogram p before", "Correlogram p after", "LISA before", "LISA after")
    }
    tab$left <- vv$residual_structure
    head_t <- c(head_t, "Structure left")
    sec <- c(sec, "<h2>What was removed from each variable</h2>",
             sprintf("<p>Share of the variance removed, %s before and after, %s and whether structure is left (tested against white noise with the same basis removed).</p>",
                     if (family == "temporal") "lag-1 autocorrelation" else "Moran's I (8 nearest neighbours)",
                     if (family == "temporal") "the Ljung-Box test of the raw series and of the corrected one," else
                       "the Moran correlogram test (Bonferroni over distance classes) and the share of sites with local structure (LISA, FDR) before and after,"),
             .html_df(tab, head_t),
             .svg_figure(function() {
               o <- seq_len(nrow(vv))
               graphics::par(mar = c(4, 7, 1, 1), las = 1)
               r <- range(c(vv$autocor_before, vv$autocor_after, 0), na.rm = TRUE)
               graphics::plot(vv$autocor_before, o, xlim = r, yaxt = "n", pch = 19, col = "grey55", xlab = "autocorrelation", ylab = "")
               graphics::segments(vv$autocor_before, o, vv$autocor_after, o, col = "grey75")
               graphics::points(vv$autocor_after, o, pch = 19, col = "#2a6fb0")
               graphics::axis(2, at = o, labels = vv$variable, cex.axis = 0.7)
               graphics::abline(v = 0, lty = 3)
               graphics::legend("bottomright", c("before", "after"), pch = 19, col = c("grey55", "#2a6fb0"), bty = "n", cex = 0.8)
             }, 6, max(3, 0.18 * nrow(vv) + 1), "Autocorrelation of each variable before (grey) and after (blue) the correction."))
    if (!is.null(v$correlograms_before)) {
      show <- utils::head(vv$variable, 12)
      sec <- c(sec, .svg_figure(function() {
        nr <- ceiling(length(show) / 3)
        graphics::par(mfrow = c(nr, 3), mar = c(3.5, 3.5, 1.5, 0.5), mgp = c(2, 0.6, 0), cex = 0.7)
        for (w in show) {
          b <- v$correlograms_before[v$correlograms_before$variable == w, ]
          a <- v$correlograms[v$correlograms$variable == w, ]
          r <- range(c(b$moran_I, a$moran_I, 0), na.rm = TRUE)
          graphics::plot(b$distance, b$moran_I, type = "b", pch = 19, col = "grey55", ylim = r, xlab = "distance", ylab = "Moran's I", main = w)
          graphics::lines(a$distance, a$moran_I, type = "b", pch = 19, col = "#2a6fb0")
          graphics::abline(h = 0, lty = 3)
        }
      }, 7.5, 2.3 * ceiling(length(show) / 3) + 0.3,
      sprintf("Moran correlograms before (grey) and after (blue) the correction%s: Moran's I by distance class%s.",
              if (length(show) < nrow(vv)) sprintf(", first %d variables", length(show)) else "",
              if (family == "spatiotemporal") ", on the mean of each site" else "")))
    }

    ds <- v$delta_summary
    rows <- vapply(seq_len(nrow(ds)), function(i) sprintf("<tr><td>%s</td><td class='n'>%d</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td></tr>",
      esc(ds$kind[i]), ds$pairs[i], num(ds$mean_delta[i], 3), num(ds$mean_abs_delta[i], 3), num(ds$max_abs_delta[i], 3)), "")
    sec <- c(sec, "<h2>How the correlations changed (Delta R)</h2><p>Corrected minus raw correlation, by kind of pair. Large changes within groups are the shared large-scale structure removed.</p>",
             .html_table(c("Pairs", "Number", "Mean change", "Mean |change|", "Largest |change|"), rows),
             .svg_figure(function() {
               d <- v$delta_R; p <- ncol(d)
               graphics::par(mar = c(5, 5, 1, 1))
               lim <- max(abs(d), 1e-6)
               graphics::image(seq_len(p), seq_len(p), d[, p:1], zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "",
                               col = grDevices::hcl.colors(41, "Blue-Red 3"))
               graphics::axis(1, at = seq_len(p), labels = colnames(d), las = 2, cex.axis = 0.6)
               graphics::axis(2, at = seq_len(p), labels = rev(colnames(d)), las = 1, cex.axis = 0.6)
             }, 6, 6, sprintf("Delta R, corrected minus raw (blue: lower after the correction; red: higher; scale &plusmn;%s).", num(max(abs(v$delta_R)), 2))))

    q3 <- v$q3
    sec <- c(sec, "<h2>Signal at the scale of the basis (Q3)</h2>",
             if (!is.finite(q3$statistic)) "<p>Not computed: the whitening mode has no projected part.</p>" else
               sprintf("<p class='%s'>Mean absolute correlation carried by the projected part: <b>%s</b> (threshold %s). %s</p>",
                       if (isTRUE(q3$warning)) "warn" else "ok", num(q3$statistic, 3), num(q3$threshold, 2),
                       if (isTRUE(q3$warning)) "Much of the association between variables was at the scale of the basis. If the classes you look for could be regions or periods, the correction removes them: check before using the corrected data for them." else
                         "Below the threshold: the association between variables was not mainly at the scale of the basis."))

    if (!is.null(v$lisa)) {
      left <- max(vv$lisa_share) > 0
      worst <- if (left) vv$variable[which.max(vv$lisa_share)] else vv$variable[which.max(vv$lisa_share_before)]
      l <- v$lisa[v$lisa$variable == worst, ]
      l0 <- v$lisa_before[v$lisa_before$variable == worst, ]
      pos <- as.matrix(coords)
      cols <- c(none = "grey85", high_high = "#c0392b", low_low = "#2a6fb0", high_low = "#f39c12", low_high = "#8e44ad")
      sec <- c(sec, "<h2>Local structure before and after (LISA)</h2>",
               sprintf("<p>%s, <b>%s</b>: local structure at %s of the sites before the correction, %s after. A region that stands out after the correction may behave differently from the rest; the LISA describes, it does not correct. The LISA of every variable is in <code>v$lisa_before</code> and <code>v$lisa</code>.</p>",
                       if (left) "Local Moran's I of the variable with the most local structure left after the correction" else
                         "No variable keeps local structure after the correction. Local Moran's I of the variable with the most local structure before it",
                       esc(worst), num(vv$lisa_share_before[vv$variable == worst], 2), num(vv$lisa_share[vv$variable == worst], 2)),
               .svg_figure(function() {
                 graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
                 for (k in 1:2) {
                   d <- if (k == 1) l0 else l
                   graphics::plot(pos[d$site, 1], pos[d$site, 2], pch = 19, cex = 0.6, col = cols[d$cluster], asp = 1, xlab = "x", ylab = "y",
                                  main = c("before", "after")[k])
                 }
                 graphics::legend("topright", names(cols), pch = 19, col = cols, bty = "n", cex = 0.6)
               }, 9, 4.8, "Sites by LISA quadrant (FDR-adjusted) before and after the correction; grey: no local structure."))
    }
    warn <- v$warnings
    notes <- v$notes
  } else {
    warn <- corrected$warnings; notes <- character(0)
  }
  if (length(warn) || length(notes)) {
    sec <- c(sec, "<h2>Warnings and notes</h2><ul>", paste0("<li>", esc(c(warn, notes)), "</li>"), "</ul>")
  }
  if (is.null(code)) code <- .correction_code(corrected, family)
  sec <- c(sec, "<h2>Reproduce</h2><pre>", esc(code), "</pre>")

  html <- c("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
            sprintf("<title>%s</title>", esc(title)),
            "<style>body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.5;background:#fff}",
            "h1{font-size:1.6rem}h2{font-size:1.2rem;margin-top:2rem;border-bottom:1px solid #ddd;padding-bottom:.2rem}",
            "table{border-collapse:collapse;font-size:.85rem;width:100%;display:block;overflow-x:auto}th,td{padding:.25rem .5rem;border-bottom:1px solid #eee;text-align:left}td.n{text-align:right;font-variant-numeric:tabular-nums}",
            ".meta{color:#666;font-size:.85rem}.warn{background:#fff4e5;border-left:4px solid #f39c12;padding:.5rem}.ok{background:#eef7ee;border-left:4px solid #2e8b57;padding:.5rem}",
            "pre{background:#f6f6f6;padding:.75rem;overflow-x:auto;font-size:.8rem}figure{margin:1rem 0}figcaption{color:#555;font-size:.85rem}svg{max-width:100%;height:auto}</style></head><body>",
            sprintf("<h1>%s</h1>", esc(title)), sec, "</body></html>")
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(file)
}

# Readable names of the signature columns, in the order they are shown.
.diagnosis_columns <- c(
  moran_I = "Moran's I", moran_p = "Moran p", correlogram_p = "correlogram p", dependence_p = "dependence p",
  vario_model = "variogram model", nugget = "nugget", sill = "sill", range = "range", range_lo = "range, low",
  range_hi = "range, high", nugget_ratio = "nugget / sill", aniso_ratio = "anisotropy ratio",
  aniso_direction = "anisotropy direction", acf1 = "ACF lag 1", phi = "AR phi", phi_lo = "phi, low", phi_hi = "phi, high",
  ar_order = "AR order", ljung_p = "Ljung-Box p", period = "period", period_strength = "seasonal strength",
  period_p = "period p", trend_r2 = "trend R2", shift_indicator = "shift indicator", shift_direction = "shift direction")

#' The diagnosis of each variable, as one table
#'
#' Joins what the diagnosis found for each variable before the correction, for the user
#' to read, explore or save: the regime, the group and its stability (from the
#' grouping), and the columns of the dependence signature that hold a value for the data
#' family (Moran's I and its p-value, the variogram model, nugget, sill, range and
#' anisotropy in space; the lag-1 autocorrelation, the AR coefficient, the Ljung-Box
#' p-value, the period and its strength and the trend in time; the shift of a moving
#' pattern in space-time; the bootstrap intervals when computed). The same table opens
#' the report of [report_correction()] and step 3 of [run_app()].
#'
#' @param signature Output of [diagnose_variables()].
#' @param grouping Output of [group_regimes()], or a vector of groups named by variable;
#'   `NULL` for the signature only.
#' @return A data frame, one row per variable (numbers rounded to 3 significant digits).
#' @examples
#' g <- sim_groups(sim_group("long", 3, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(100, g, seed = 1)
#' sig <- diagnose_variables(s)
#' diagnosis_table(sig, group_regimes(sig))
#' @export
diagnosis_table <- function(signature, grouping = NULL) {
  s <- as.data.frame(signature)
  cols <- names(.diagnosis_columns)[names(.diagnosis_columns) %in% names(s)]
  cols <- cols[vapply(cols, function(k) if (is.numeric(s[[k]])) any(is.finite(s[[k]])) else any(!is.na(s[[k]])), TRUE)]
  out <- data.frame(variable = s$variable, stringsAsFactors = FALSE)
  if (is.character(grouping)) grouping <- data.frame(variable = names(grouping), group = unname(grouping), stringsAsFactors = FALSE)
  if (!is.null(grouping)) {
    g <- as.data.frame(grouping)[match(s$variable, grouping$variable), , drop = FALSE]
    for (k in intersect(c("type", "group", "stability", "ambiguous"), names(g))) out[[k]] <- g[[k]]
  }
  for (k in cols) out[[k]] <- if (is.numeric(s[[k]])) signif(s[[k]], 3) else s[[k]]
  out
}

.html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

.html_table <- function(head, rows) {
  paste0("<table><thead><tr>", paste0("<th>", .html_escape(head), "</th>", collapse = ""), "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

# A data frame as an HTML table: numbers to the right, TRUE/FALSE as yes/no, NA empty.
.html_df <- function(d, head = names(d)) {
  cell <- function(x) if (is.logical(x)) ifelse(is.na(x), "<td></td>", ifelse(x, "<td>yes</td>", "<td>no</td>")) else
    if (is.numeric(x)) ifelse(is.finite(x), sprintf("<td class='n'>%s</td>", format(x, trim = TRUE, drop0trailing = TRUE)), "<td></td>") else
      ifelse(is.na(x), "<td></td>", sprintf("<td>%s</td>", .html_escape(x)))
  cells <- vapply(d, cell, character(nrow(d)))
  if (!is.matrix(cells)) cells <- matrix(cells, nrow = nrow(d))
  .html_table(head, paste0("<tr>", apply(cells, 1, paste, collapse = ""), "</tr>"))
}

# A base-graphics figure as inline SVG (needs cairo); a note when it cannot be drawn.
.svg_figure <- function(draw, width, height, caption) {
  if (!isTRUE(capabilities("cairo"))) return(sprintf("<p><i>Figure not drawn (no cairo device): %s</i></p>", caption))
  f <- tempfile(fileext = ".svg")
  grDevices::svg(f, width = width, height = height)
  ok <- tryCatch({ draw(); TRUE }, error = function(e) FALSE, finally = grDevices::dev.off())
  if (!ok) return(sprintf("<p><i>Figure not drawn: %s</i></p>", caption))
  s <- readLines(f, warn = FALSE)
  unlink(f)
  s <- s[!grepl("^<\\?xml", s)]
  # cairo names the glyphs and clip paths of every figure alike (glyph0-1, clip1, ...); in
  # one page the later figures would draw with the glyphs of the first, so each figure
  # gets its own prefix
  p <- paste0(sub("^file", "f", basename(tools::file_path_sans_ext(f))), "-")
  s <- gsub('id="', paste0('id="', p), s, fixed = TRUE)
  s <- gsub('href="#', paste0('href="#', p), s, fixed = TRUE)
  s <- gsub("url(#", paste0("url(#", p), s, fixed = TRUE)
  paste0("<figure>", paste(s, collapse = "\n"), "<figcaption>", caption, "</figcaption></figure>")
}

# Template of the R code that reproduces a correction, from its settings.
.correction_code <- function(corrected, family) {
  pos <- switch(family, spatial = "coords = coords", temporal = "time = time", spatiotemporal = "coords = coords, index = index")
  menus <- unique(unlist(lapply(corrected$bases, function(b) if (!is.null(attr(b, "menu"))) paste(utils::capture.output(print(attr(b, "menu")))[-1], collapse = "; "))))
  c("library(depprep)",
    "# X: the data (observations x variables), with its positions",
    sprintf("sig <- diagnose_variables(X, %s)", pos),
    "grouping <- group_regimes(sig)",
    sprintf("cr <- correct(X, grouping, %s, mode = \"%s\", scope = \"%s\"%s)", pos, corrected$mode, corrected$scope,
            if (!is.null(corrected$data_basis)) sprintf(", data_basis = \"%s\"", corrected$data_basis) else ""),
    if (length(menus)) paste("# built with the menu:", menus),
    sprintf("v <- verify(cr, X, %s)", pos),
    "report_correction(cr, X, v, file = \"depprep_report.html\", signature = sig, grouping = grouping,",
    sprintf("                  %s)", pos))
}
