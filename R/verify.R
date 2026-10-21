#' Verify what a correction removed and what it kept
#'
#' Summarises a [correct()] result against the raw data:
#' * for each variable, the share of its variance that was removed and its
#'   autocorrelation before and after (Moran's I with 8 nearest neighbours in space,
#'   on the site means in space-time; lag-1 autocorrelation in time);
#' * Delta R = R corrected - R raw for every pair, summarised within groups, across
#'   groups and for pairs involving variables without dependence;
#' * the share of degrees of freedom used by each group basis and each union, with the
#'   warnings of [correct()].
#'
#' @param corrected Output of [correct()] (conditional or whiten mode).
#' @param X The raw data (matrix or `depprep_sim`), as given to [correct()].
#' @param coords,time,index As in [correct()].
#' @return A list of class `depprep_verification` with `variables` (data frame),
#'   `delta_R` (matrix), `delta_summary` (data frame), `kn`, `union_kn` and `warnings`.
#' @examples
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(300, g, seed = 1)
#' cr <- correct(s, group_regimes(diagnose_variables(s, nsim = 99)))
#' verify(cr, s)$delta_summary
#' @export
verify <- function(corrected, X, coords = NULL, time = NULL, index = NULL) {
  stopifnot(inherits(corrected, "depprep_corrected"))
  if (corrected$mode == "scale") stop("`verify()` needs the conditional or whiten mode.", call. = FALSE)
  if (inherits(X, "depprep_sim")) {
    d <- .data_parts(X); X <- d$X; coords <- d$coords; time <- d$time; index <- d$index
  }
  X <- as.matrix(X)
  E <- corrected$X
  grp <- corrected$group
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else if (!is.null(coords)) "spatial" else "temporal"
  # rows of the corrected data: AR filters drop the first times (time-ordered output)
  ord <- if (family == "temporal") order(time) else seq_len(nrow(X))
  Xo <- X[ord, , drop = FALSE]
  keep <- utils::tail(seq_len(nrow(Xo)), nrow(E))
  autocor <- switch(family,
    spatial = { lw <- .knn_weights(coords); function(z, rows) .moran(z, lw)[["moran_I"]] },
    temporal = function(z, rows) stats::acf(z, lag.max = 1, plot = FALSE)$acf[2],
    spatiotemporal = {
      lw <- .knn_weights(coords)
      function(z, rows) .moran(tapply(z, index$site[rows], mean), lw)[["moran_I"]]
    })
  vars <- data.frame(
    variable = colnames(X), group = grp,
    removed = vapply(seq_len(ncol(X)), function(j) 1 - stats::var(E[, j]) / stats::var(Xo[keep, j]), 0),
    autocor_before = vapply(seq_len(ncol(X)), function(j) autocor(Xo[, j], seq_len(nrow(Xo))), 0),
    autocor_after = vapply(seq_len(ncol(X)), function(j) autocor(E[, j], keep), 0),
    stringsAsFactors = FALSE)
  R_raw <- stats::cor(X)
  dR <- corrected$R - R_raw
  ut <- upper.tri(dR)
  kind <- outer(grp, grp, function(a, b) ifelse(a == "none" | b == "none", "involving no dependence",
                                                 ifelse(a == b, "within groups", "across groups")))
  delta_summary <- do.call(rbind, lapply(split(seq_len(sum(ut)), kind[ut]), function(i) {
    v <- dR[ut][i]
    data.frame(pairs = length(v), mean_delta = mean(v), mean_abs_delta = mean(abs(v)), max_abs_delta = max(abs(v)))
  }))
  delta_summary <- cbind(kind = rownames(delta_summary), delta_summary, row.names = NULL)
  structure(list(variables = vars, delta_R = dR, delta_summary = delta_summary, kn = corrected$kn,
                 union_kn = corrected$union_kn, warnings = corrected$warnings),
            class = "depprep_verification")
}

#' @export
print.depprep_verification <- function(x, ...) {
  cat("<depprep_verification>\n")
  print(x$delta_summary, digits = 3, row.names = FALSE)
  cat(sprintf("Median share of variance removed: %.2f; median autocorrelation before %.2f, after %.2f\n",
              stats::median(x$variables$removed), stats::median(x$variables$autocor_before),
              stats::median(x$variables$autocor_after)))
  if (length(x$warnings)) cat("Warnings:", paste(x$warnings, collapse = "; "), "\n")
  invisible(x)
}
