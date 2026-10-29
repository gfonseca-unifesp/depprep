#' Write a report of a correction for the end user
#'
#' Writes one self-contained HTML file (no external files, no pandoc) that says what the
#' correction did, in the order a user checks it: the bases of each group and the
#' degrees of freedom they used, how much of each variable was removed, the
#' autocorrelation before and after, how the correlations between variables changed
#' (Delta R), the warning of signal at the scale of the basis (Q3), where local structure
#' is left (LISA, in space), the warnings, and the R code that reproduces the correction.
#'
#' @param corrected Output of [correct()].
#' @param X The raw data (matrix or `depprep_sim`), as given to [correct()].
#' @param verification Output of [verify()]; computed (with `nperm = 99`) when `NULL`.
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
                              index = NULL, code = NULL, title = "depprep: correction report") {
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
    rows <- vapply(seq_len(nrow(vv)), function(i) sprintf(
      "<tr><td>%s</td><td>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td>%s</td>%s</tr>",
      esc(vv$variable[i]), esc(vv$group[i]), num(vv$removed[i]), num(vv$autocor_before[i]), num(vv$autocor_after[i]),
      if (isTRUE(vv$residual_structure[i])) "yes" else "no", if (!is.null(vv$lisa_share)) sprintf("<td class='n'>%s</td>", num(vv$lisa_share[i])) else ""), "")
    sec <- c(sec, "<h2>What was removed from each variable</h2>",
             sprintf("<p>Share of the variance removed, %s before and after, and whether structure is left (tested against white noise with the same basis removed)%s.</p>",
                     if (family == "temporal") "lag-1 autocorrelation" else "Moran's I (8 nearest neighbours)",
                     if (!is.null(vv$lisa_share)) "; LISA: share of sites with local structure left" else ""),
             .html_table(c("Variable", "Group", "Removed", "Autocor. before", "Autocor. after", "Structure left", if (!is.null(vv$lisa_share)) "LISA"), rows),
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
      worst <- vv$variable[which.max(vv$lisa_share)]
      l <- v$lisa[v$lisa$variable == worst, ]
      pos <- as.matrix(coords)
      cols <- c(none = "grey85", high_high = "#c0392b", low_low = "#2a6fb0", high_low = "#f39c12", low_high = "#8e44ad")
      sec <- c(sec, "<h2>Where local structure is left (LISA)</h2>",
               sprintf("<p>Local Moran's I of the corrected variable with the most local structure, <b>%s</b> (%s of the sites). A region that stands out may behave differently from the rest; the LISA describes, it does not correct.</p>",
                       esc(worst), num(max(vv$lisa_share), 2)),
               .svg_figure(function() {
                 graphics::par(mar = c(4, 4, 1, 1))
                 graphics::plot(pos[l$site, 1], pos[l$site, 2], pch = 19, cex = 0.8, col = cols[l$cluster], asp = 1, xlab = "x", ylab = "y")
                 graphics::legend("topright", names(cols), pch = 19, col = cols, bty = "n", cex = 0.7)
               }, 6, 5.5, "Sites by LISA quadrant (FDR-adjusted); grey: no local structure."))
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

.html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

.html_table <- function(head, rows) {
  paste0("<table><thead><tr>", paste0("<th>", .html_escape(head), "</th>", collapse = ""), "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
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
    "report_correction(cr, X, v, file = \"depprep_report.html\", ", sprintf("                  %s)", pos))
}
