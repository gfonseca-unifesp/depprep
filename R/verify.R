#' Verify what a correction removed and what it kept
#'
#' Summarises a [correct()] result against the raw data in three parts.
#'
#' **Residual structure** (what the correction left; Legendre and Legendre 2012,
#' sections 12.3 and 13.1):
#' * for each variable, the share of its variance that was removed and its
#'   autocorrelation before and after (Moran's I with 8 nearest neighbours in space,
#'   on the site means in space-time; lag-1 autocorrelation in time);
#' * in space (and on the site means in space-time), a Moran correlogram of each
#'   corrected variable by distance class. The classes are those of
#'   [vegan::mantel.correlog()]: equal widths, Sturges' rule for their number, and the
#'   far classes dropped when some site has no partner in them. Each class is tested
#'   two-sided. The correlogram is significant when its smallest p-value passes a
#'   Bonferroni correction over the classes; `extent` counts the first classes that stay
#'   significant with the progressive Bonferroni correction (class k tested at
#'   `alpha / k`), that is, how far the structure that is left extends;
#' * in time, a Ljung-Box test of each corrected variable up to lag
#'   `max(2 * period, 10 * log10(n))` (at most `n / 3`), where the period is `period`,
#'   or the one stored in `corrected`, or the one of the seasonal terms of its bases;
#'   without a period, up to `10 * log10(n)`;
#' * on the whole corrected matrix (standardised variables, Euclidean distance), a
#'   multivariate Mantel correlogram (classes and statistics of
#'   [vegan::mantel.correlog()], progressive Bonferroni correction), against the
#'   distance between sites in space and space-time (on the site means) and between times
#'   in time.
#'
#' Residuals on a basis are not exchangeable: removing the broad scales leaves a
#' negative autocorrelation at short distances and lags, which the usual tests would
#' report as structure. The null is therefore white noise with the same basis (and the
#' mean) removed: the exact moments of Moran's I for residuals (Cliff and Ord, as in
#' [spdep::lm.morantest()]; the `expected` column), and, for the Ljung-Box and Mantel
#' tests, `nperm` permutations of the rows projected off the basis again (Freedman and
#' Lane). In the whiten mode only the mean is removed, so these are the usual tests.
#'
#' **What was removed** (the logic of variation partitioning). In the conditional
#' mode each variable splits into a projected (large-scale) part `P x` and a residual
#' part `(I - P) x`, and so does the covariance of every pair, exactly:
#' `x'y = x'P y + x'(I - P) y`. On the union of all bases (the projector of the data
#' matrix when `data_basis = "union"`; with the default `"own"` each column of
#' `corrected$X` has its own basis, and these parts still describe the union) this gives
#' `cov_large + cov_residual = cov(X)`, `R_large`
#' (the correlation of the projected parts) and `large_share` (the share of each
#' covariance, and on the diagonal of each variance, that lay in the projected part).
#' `large_share_pairwise` gives the same share with the projector `correct()` uses for
#' each pair in `R`: the group basis within a group, the union of the two group bases
#' across groups, none between variables without dependence. Whitening has no
#' projected part: these elements are then `NA`, and `notes` says so.
#'
#' **Warning for question Q3** ("signal at the scale of the basis"). Classes that are
#' contiguous in space or time (regions, periods) live at the scale of the basis, and the
#' correction removes them together with the dependence. `q3$statistic` is the mean
#' absolute association between variables carried by the projected part: the mean over
#' all pairs of `|cov_large| / (sd_x sd_y)`, the projected part of each correlation on
#' the union of the bases. The warning is raised when it exceeds `q3_threshold`.
#' `q3$share` gives that association as a share of the total, `sum |cov_large| /
#' sum (|cov_large| + |cov_residual|)`, for information only: it is high whenever little
#' association is left in the residuals, also in the null truth, so it does not tell
#' classes from shared fields. The statistic grows with the size of the basis also without
#' contiguous classes (the projected parts of independent autocorrelated variables are
#' correlated by chance), so the default threshold depends on the data family: 0.16 in
#' space, 0.19 in time, 0.13 in space-time. They were calibrated in
#' `analysis/07d_q3_recalibration.R` on the stage 7 grid (the midpoint between the 97.5th
#' percentile of the truths without contiguous classes and the 5th percentile of the
#' contiguous truth) and validated there on new seeds: at most 1.2% of warnings without
#' contiguous classes, at least 98.8% with them. In those designs the large-scale fields are shared within groups of
#' 6-8 of the 30 variables, while contiguous classes reach every variable. The warning
#' cannot tell a large-scale field shared by most variables from regional classes: it
#' says that much of the association between variables is at the scale of the basis.
#' When it fires, check whether the classes sought could be regions or periods before
#' using the corrected data for them.
#'
#' @param corrected Output of [correct()] (conditional or whiten mode).
#' @param X The raw data (matrix or `depprep_sim`), as given to [correct()].
#' @param coords,time,index As in [correct()].
#' @param period Known seasonal period, for the Ljung-Box lag in time; `NULL` to take
#'   it from `corrected` when it records one.
#' @param alpha Level of the residual-structure tests.
#' @param nperm Permutations of the Mantel and Ljung-Box tests; 0 skips the Mantel
#'   correlogram and takes the Ljung-Box p-values from the chi-square distribution.
#' @param q3_threshold Mean absolute association carried by the projected part above
#'   which the Q3 warning is raised; `NULL` takes the threshold of the data family (see
#'   Details).
#' @param seed Seed of the permutations.
#' @return A list of class `depprep_verification` with
#'   * `variables` (data frame: group, `removed`, `autocor_before`, `autocor_after`,
#'     `large_share`, then `correlogram_p` and `extent` in space and space-time or
#'     `lb_lag` and `lb_p` in time, and `residual_structure`);
#'   * `correlograms` (data frame, one row per variable and distance class: `class`,
#'     `distance`, `pairs`, `moran_I`, `expected`, `p`, `p_progressive`; `NULL` in time);
#'   * `mantel` (data frame, one row per distance class: `class`, `distance`, `pairs`,
#'     `mantel_r`, `p`, `p_progressive`; `NULL` when `nperm = 0`) and `mantel_p` (its
#'     Bonferroni global test);
#'   * `R_large`, `cov_large`, `cov_residual`, `large_share`, `large_share_pairwise`
#'     (matrices, see above);
#'   * `q3` (list: `statistic`, `share`, `threshold`, `warning`);
#'   * `delta_R` (corrected R minus raw R), `delta_summary` (by kind of pair), `kn`,
#'     `union_kn`, `warnings` (those of [correct()] plus the Q3 warning) and `notes`.
#' @references Legendre, P. and Legendre, L. (2012). *Numerical Ecology*, 3rd ed.
#'   Elsevier. Sections 12.3, 13.1.1, 13.1.6 and 13.5.
#' @examples
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(300, g, seed = 1)
#' cr <- correct(s, group_regimes(diagnose_variables(s, nsim = 99)))
#' v <- verify(cr, s, nperm = 99)
#' v$delta_summary
#' v$variables[, c("variable", "large_share", "correlogram_p", "extent")]
#' v$q3    # fires: one large-scale field is shared by 4 of the 6 variables (see Details)
#' @importFrom vegan mantel.correlog
#' @export
verify <- function(corrected, X, coords = NULL, time = NULL, index = NULL, period = NULL, alpha = 0.05,
                   nperm = 199, q3_threshold = NULL, seed = 1L) {
  stopifnot(inherits(corrected, "depprep_corrected"))
  if (corrected$mode == "scale") stop("`verify()` needs the conditional or whiten mode.", call. = FALSE)
  if (inherits(X, "depprep_sim")) {
    d <- .data_parts(X); X <- d$X; coords <- d$coords; time <- d$time; index <- d$index
  }
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  E <- corrected$X
  grp <- corrected$group
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else if (!is.null(coords)) "spatial" else "temporal"
  if (is.null(q3_threshold)) q3_threshold <- .q3_thresholds[[family]]
  # rows in time order in time: the conditional mode keeps the rows of X, the AR filters of
  # the whitening return them in time order without the first times
  ord <- if (family == "temporal") order(time) else seq_len(nrow(X))
  Xo <- X[ord, , drop = FALSE]
  # with data_basis = "own" (the default, D17) every column of the data matrix holds that
  # variable's residual on its own basis, in the rows of X, in any mode
  own <- identical(corrected$data_basis, "own")
  if (corrected$mode == "conditional" || own) E <- E[ord, , drop = FALSE]
  keep <- utils::tail(seq_len(nrow(Xo)), nrow(E))
  autocor <- switch(family,
    spatial = { lw <- .knn_weights(coords); function(z, rows) .moran(z, lw)[["moran_I"]] },
    temporal = function(z, rows) stats::acf(z, lag.max = 1, plot = FALSE)$acf[2],
    spatiotemporal = {
      lw <- .knn_weights(coords)
      function(z, rows) .moran(tapply(z, index$site[rows], mean), lw)[["moran_I"]]
    })
  U <- if (corrected$mode == "conditional") .union_of_bases(corrected, nrow(X)) else matrix(0, nrow(X), 0)
  parts <- if (corrected$mode == "conditional") .large_scale_parts(X, corrected, U) else .no_large_parts(X)
  vars <- data.frame(
    variable = colnames(X), group = grp,
    removed = vapply(seq_len(ncol(X)), function(j) 1 - stats::var(E[, j]) / stats::var(Xo[keep, j]), 0),
    autocor_before = vapply(seq_len(ncol(X)), function(j) autocor(Xo[, j], seq_len(nrow(Xo))), 0),
    autocor_after = vapply(seq_len(ncol(X)), function(j) autocor(E[, j], keep), 0),
    large_share = unname(diag(parts$large_share)),
    stringsAsFactors = FALSE)

  # residual structure: per variable, then on the whole matrix. Residuals on a basis are
  # not exchangeable (removing the broad scales leaves a negative autocorrelation at short
  # distances), so the null is white noise with the same basis removed: exact moments of
  # Moran's I for residuals, and permutations of the rows projected again off the basis
  # (Freedman and Lane) for the Ljung-Box and Mantel tests.
  site_means <- function(M) {
    if (!ncol(M)) return(matrix(0, nrow(as.matrix(coords)), 0))
    apply(M, 2, function(z) tapply(z, index$site[keep], mean))
  }
  Z <- switch(family, spatiotemporal = site_means(E), E)
  fixed <- function(B) .fixed_space(switch(family, temporal = B[ord, , drop = FALSE][keep, , drop = FALSE],
                                            spatiotemporal = site_means(B[keep, , drop = FALSE]), B))
  Fq <- fixed(U)
  # the null of each variable removes the basis that variable was corrected with
  own_b <- if (own) .own_bases(corrected, nrow(X), colnames(X)) else NULL
  Fq_of <- function(j) if (is.null(own_b)) Fq else fixed(own_b[[j]])
  Fq_all <- if (is.null(own_b)) Fq else {
    nz <- Filter(ncol, own_b)
    fixed(if (length(nz)) unclass(suppressWarnings(do.call(union_basis, c(nz, kn_max = 1)))) else matrix(0, nrow(X), 0))
  }
  correlograms <- NULL
  if (family == "temporal") {
    T_known <- .known_period(period, corrected)
    n_t <- nrow(E)
    lag <- min(ceiling(max(if (is.finite(T_known)) 2 * T_known else 0, 10 * log10(n_t))), floor(n_t / 3))
    vars$lb_lag <- lag
    vars$lb_p <- if (is.null(own_b)) .ljung_box(Z, lag, Fq, nperm, seed) else
      vapply(seq_len(ncol(Z)), function(j) .ljung_box(Z[, j, drop = FALSE], lag, Fq_of(j), nperm, seed + j), 0)
    vars$residual_structure <- vars$lb_p <= alpha
    position <- matrix(time[ord][keep])
  } else {
    position <- coords
    correlograms <- if (is.null(own_b)) .moran_correlogram(Z, position, Fq) else
      do.call(rbind, lapply(seq_len(ncol(Z)), function(j) .moran_correlogram(Z[, j, drop = FALSE], position, Fq_of(j))))
    by_var <- split(correlograms, factor(correlograms$variable, levels = colnames(X)))
    vars$correlogram_p <- vapply(by_var, function(d) if (all(is.na(d$p))) NA_real_ else
      min(1, min(d$p, na.rm = TRUE) * sum(!is.na(d$p))), 0, USE.NAMES = FALSE)
    vars$extent <- vapply(by_var, function(d) {
      s <- !is.na(d$p_progressive) & d$p_progressive <= alpha
      as.integer(min(c(which(!s), length(s) + 1)) - 1)
    }, 0L, USE.NAMES = FALSE)
    vars$residual_structure <- !is.na(vars$correlogram_p) & vars$correlogram_p <= alpha
  }
  mantel <- .mantel_correlogram(Z, position, Fq_all, nperm, seed)
  mantel_p <- if (is.null(mantel)) NA_real_ else min(1, min(mantel$p, na.rm = TRUE) * sum(!is.na(mantel$p)))

  # what the correction did to the correlations
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

  # question Q3: was the association at the scale of the basis?
  q3_stat <- .large_association(parts$cov_large, parts$cov_residual)
  q3_warn <- is.finite(q3_stat[["mean"]]) && q3_stat[["mean"]] > q3_threshold
  warn <- corrected$warnings
  if (q3_warn) {
    warn <- c(warn, sprintf(paste("signal at the scale of the basis: the projected part carried a mean absolute",
                                  "correlation of %.2f between variables (threshold %.2f); classes that are regions or",
                                  "periods would be removed by the correction"), q3_stat[["mean"]], q3_threshold))
  }
  structure(list(variables = vars, correlograms = correlograms, mantel = mantel, mantel_p = mantel_p,
                 R_large = parts$R_large, cov_large = parts$cov_large, cov_residual = parts$cov_residual,
                 large_share = parts$large_share, large_share_pairwise = parts$large_share_pairwise,
                 q3 = list(statistic = q3_stat[["mean"]], share = q3_stat[["share"]], threshold = q3_threshold,
                           warning = q3_warn),
                 delta_R = dR, delta_summary = delta_summary, kn = corrected$kn, union_kn = corrected$union_kn,
                 warnings = warn, notes = parts$notes),
            class = "depprep_verification")
}

#' @export
print.depprep_verification <- function(x, ...) {
  cat("<depprep_verification>\n")
  print(x$delta_summary, digits = 3, row.names = FALSE)
  cat(sprintf("Median share of variance removed: %.2f; median autocorrelation before %.2f, after %.2f\n",
              stats::median(x$variables$removed), stats::median(x$variables$autocor_before),
              stats::median(x$variables$autocor_after)))
  cat(sprintf("Variables with residual structure: %d of %d; Mantel correlogram p = %s\n",
              sum(x$variables$residual_structure), nrow(x$variables),
              if (is.na(x$mantel_p)) "not computed" else format.pval(x$mantel_p, digits = 2)))
  if (is.finite(x$q3$statistic)) {
    cat(sprintf("Mean absolute correlation carried by the projected part: %.2f (Q3 threshold %.2f; %.0f%% of the association)\n",
                x$q3$statistic, x$q3$threshold, 100 * x$q3$share))
  }
  if (length(x$warnings)) cat("Warnings:", paste(x$warnings, collapse = "; "), "\n")
  if (length(x$notes)) cat("Notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

# Default thresholds of the Q3 warning by data family (amendment protocol-v1.2,
# analysis/07d_q3_recalibration.R).
.q3_thresholds <- c(spatial = 0.16, temporal = 0.19, spatiotemporal = 0.13)

# Basis each column of the data matrix was corrected with when data_basis = "own": the
# variable's own selection, or its group's basis when the bases were given; empty for the
# variables without dependence.
.own_bases <- function(corrected, n, vars) {
  empty <- matrix(0, n, 0)
  lapply(seq_along(vars), function(j) {
    g <- corrected$group[[j]]
    if (g == "none") return(empty)
    b <- if (!is.null(corrected$variable_bases)) corrected$variable_bases[[vars[j]]] else corrected$bases[[g]]
    if (is.null(b)) empty else unclass(b)
  })
}

# Union of the bases of a correction (the projector of its data matrix).
.union_of_bases <- function(corrected, n) {
  groups <- setdiff(unique(corrected$group), "none")
  b <- lapply(groups, function(g) unclass(corrected$bases[[g]]))
  if (!length(b) || !sum(vapply(b, ncol, 0L))) return(matrix(0, n, 0))
  unclass(suppressWarnings(do.call(union_basis, c(b, kn_max = 1))))
}

# Projected (large-scale) and residual parts of the covariance of X under a conditional
# correction: on U, the union of all its bases (the projector of the data matrix), and
# pair by pair with the projector correct() uses for R (group basis within a group, union
# of the two bases across groups, none between variables without dependence).
.large_scale_parts <- function(X, corrected, U) {
  Xc <- scale(as.matrix(X), scale = FALSE)
  n <- nrow(Xc)
  grp <- corrected$group
  groups <- setdiff(unique(grp), "none")
  empty <- matrix(0, n, 0)
  basis_of <- function(g) if (g == "none") empty else unclass(corrected$bases[[g]])
  project <- function(B, M) if (ncol(B)) B %*% crossprod(B, M) else M * 0
  L <- project(U, Xc)
  C <- crossprod(Xc) / (n - 1)
  CL <- crossprod(L) / (n - 1)
  CE <- crossprod(Xc - L) / (n - 1)
  sd_L <- sqrt(diag(CL))
  R_large <- CL / outer(sd_L, sd_L)
  R_large[!is.finite(R_large)] <- NA
  ug <- c(groups, if (any(grp == "none")) "none")
  S <- matrix(NA_real_, ncol(Xc), ncol(Xc))
  for (a in seq_along(ug)) for (b in a:length(ug)) {
    ia <- which(grp == ug[a]); ib <- which(grp == ug[b])
    B <- if (a == b) basis_of(ug[a]) else {
      bb <- list(basis_of(ug[a]), basis_of(ug[b]))
      if (sum(vapply(bb, ncol, 0L))) unclass(suppressWarnings(do.call(union_basis, c(bb, kn_max = 1)))) else empty
    }
    S[ia, ib] <- crossprod(Xc[, ia, drop = FALSE], project(B, Xc[, ib, drop = FALSE])) /
      crossprod(Xc[, ia, drop = FALSE], Xc[, ib, drop = FALSE])
    S[ib, ia] <- t(S[ia, ib])
  }
  share <- CL / C
  dn <- list(colnames(X), colnames(X))
  dimnames(CL) <- dimnames(CE) <- dimnames(R_large) <- dimnames(S) <- dimnames(share) <- dn
  list(R_large = R_large, cov_large = CL, cov_residual = CE, large_share = share, large_share_pairwise = S,
       notes = character(0))
}

.no_large_parts <- function(X) {
  p <- ncol(X)
  na <- matrix(NA_real_, p, p, dimnames = list(colnames(X), colnames(X)))
  list(R_large = na, cov_large = na, cov_residual = na, large_share = na, large_share_pairwise = na,
       notes = paste("whitening has no projected part: R_large, the covariance shares and the Q3 statistic are NA",
                     "(use the conditional mode to see what was removed)"))
}

# Association between variables carried by the projected part, over all pairs of
# standardised variables: the mean of |cov_large| (the statistic of the Q3 warning) and
# its share of the total, sum |cov_large| / sum (|cov_large| + |cov_residual|).
.large_association <- function(CL, CE) {
  if (anyNA(CL) || nrow(CL) < 2) return(c(mean = NA_real_, share = NA_real_))
  s <- sqrt(diag(CL) + diag(CE))
  ut <- upper.tri(CL)
  aL <- abs(CL / outer(s, s))[ut]; aE <- abs(CE / outer(s, s))[ut]
  c(mean = mean(aL), share = sum(aL) / sum(aL + aE))
}

# Seasonal period known to the verification: the argument, the one recorded in the
# corrected object, or the largest period of the seasonal terms of its bases.
.known_period <- function(period, corrected) {
  if (!is.null(period)) return(as.numeric(period))
  if (!is.null(corrected$period) && any(is.finite(corrected$period))) return(max(corrected$period, na.rm = TRUE))
  terms <- unlist(lapply(corrected$bases, function(b) attr(b, "terms")))
  p <- suppressWarnings(as.numeric(sub(".*_p", "", grep("_p[0-9.]+$", terms, value = TRUE))))
  if (length(p) && any(is.finite(p))) max(p, na.rm = TRUE) else NA_real_
}

# Orthonormal basis of the space a correction fixes: the constant (the centring) and the
# columns of B.
.fixed_space <- function(B) {
  q <- qr(cbind(1, B), tol = 1e-7)
  qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

# One draw of the null of residual structure: the rows of Z permuted, then projected off
# the fixed space again (Freedman and Lane).
.null_draw <- function(Z, Fq) {
  Zp <- Z[sample.int(nrow(Z)), , drop = FALSE]
  Zp - Fq %*% crossprod(Fq, Zp)
}

# Ljung-Box statistic of every column of Z (rows in time order) up to `lag`.
.ljung_box_q <- function(Z, lag) {
  Z <- scale(as.matrix(Z), scale = FALSE)
  n <- nrow(Z)
  den <- colSums(Z^2)
  q <- 0
  for (h in seq_len(lag)) q <- q + (colSums(Z[-(1:h), , drop = FALSE] * Z[1:(n - h), , drop = FALSE]) / den)^2 / (n - h)
  n * (n + 2) * q
}

# Ljung-Box p-value of every column of Z: from the permutation null of residuals on the
# fixed space, or from the chi-square distribution when nperm = 0.
.ljung_box <- function(Z, lag, Fq, nperm, seed) {
  q <- .ljung_box_q(Z, lag)
  if (nperm < 1) return(unname(stats::pchisq(q, lag, lower.tail = FALSE)))
  q0 <- .with_stream(seed, .streams[["permutation"]],
                     vapply(seq_len(nperm), function(b) .ljung_box_q(.null_draw(Z, Fq), lag), numeric(ncol(Z))))
  q0 <- matrix(q0, ncol(Z))
  unname((1 + rowSums(q0 >= q)) / (nperm + 1))
}

# Distance classes of vegan::mantel.correlog(): equal widths between the smallest and the
# largest distance, ceiling(1 + log2(number of pairs)) classes (Sturges' rule), and the
# classes past the first half dropped when some site has no partner in them.
.distance_classes <- function(D) {
  d <- D[upper.tri(D)]
  n_class <- ceiling(1 + log2(length(d)))
  lo <- min(d); width <- (max(d) - lo) / n_class
  cls <- pmin(pmax(ceiling((D - lo) / width), 1L), n_class)
  diag(cls) <- 0L
  mid <- lo + width * (seq_len(n_class) - 0.5)
  tested <- vapply(seq_len(n_class), function(k) {
    W <- cls == k
    any(W) && (k <= n_class %/% 2 || all(rowSums(W) > 0))
  }, TRUE)
  list(class = cls, mid = mid, tested = which(tested))
}

# Moran correlogram of every column of Z (rows = sites at `coords`) after removing the
# fixed space Fq: Moran's I in each tested distance class (binary weights), two-sided
# normal test with the exact moments of I for residuals (Cliff and Ord; as in
# spdep::lm.morantest), and the progressive Bonferroni correction (class k at alpha / k).
.moran_correlogram <- function(Z, coords, Fq) {
  Z <- as.matrix(Z)
  if (is.null(colnames(Z))) colnames(Z) <- paste0("V", seq_len(ncol(Z)))
  Z <- Z - Fq %*% crossprod(Fq, Z)
  D <- as.matrix(stats::dist(as.matrix(coords)))
  n <- nrow(D); k <- ncol(Fq)
  cl <- .distance_classes(D)
  m2 <- colSums(Z^2)
  out <- do.call(rbind, lapply(seq_along(cl$tested), function(j) {
    W <- (cl$class == cl$tested[j]) * 1
    S0 <- sum(W)
    MW <- W - Fq %*% crossprod(Fq, W)
    tr1 <- sum(diag(MW)); tr2 <- sum(MW * t(MW))
    EI <- n / S0 * tr1 / (n - k)
    VI <- (n / S0)^2 * (2 * tr2 + tr1^2) / ((n - k) * (n - k + 2)) - EI^2
    I <- n / S0 * colSums(Z * (W %*% Z)) / m2
    p <- 2 * stats::pnorm(-abs(I - EI) / sqrt(VI))
    p[m2 < 1e-12] <- NA
    data.frame(variable = colnames(Z), class = j, distance = cl$mid[cl$tested[j]], pairs = S0 / 2,
               moran_I = unname(I), expected = EI, p = unname(p), p_progressive = pmin(1, j * unname(p)),
               stringsAsFactors = FALSE)
  }))
  out <- out[order(match(out$variable, colnames(Z)), out$class), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Multivariate Mantel correlogram of the standardised columns of Z against the distance
# between positions: classes and statistics from vegan::mantel.correlog(), p-values from
# the permutation null of residuals on the fixed space (one-sided in the direction of the
# statistic, as in vegan), progressive Bonferroni correction.
.mantel_correlogram <- function(Z, position, Fq, nperm, seed) {
  if (nperm < 1) return(NULL)
  Z <- as.matrix(Z)
  Z <- Z - Fq %*% crossprod(Fq, Z)
  dz <- function(M) { M <- scale(M); M[!is.finite(M)] <- 0; stats::dist(M) }
  D_geo <- stats::dist(as.matrix(position))
  mc <- vegan::mantel.correlog(dz(Z), D.geo = D_geo, nperm = 0, mult = "none")
  res <- as.data.frame(mc$mantel.res)[, 1:3]
  names(res) <- c("distance", "pairs", "mantel_r")
  res$class <- seq_len(nrow(res))
  res <- res[!is.na(res$mantel_r), , drop = FALSE]
  cls <- .distance_classes(as.matrix(D_geo))$class
  ind <- vapply(res$class, function(k) as.numeric(cls[lower.tri(cls)] == k), numeric(length(D_geo)))
  r_of <- function(M) -as.vector(stats::cor(as.vector(dz(M)), ind))
  r0 <- .with_stream(seed, .streams[["permutation"]],
                     vapply(seq_len(nperm), function(b) r_of(.null_draw(Z, Fq)), numeric(nrow(res))))
  r0 <- matrix(r0, nrow(res))
  res$pairs <- res$pairs / 2
  res$p <- (1 + rowSums(sign(res$mantel_r) * r0 >= abs(res$mantel_r))) / (nperm + 1)
  res$p_progressive <- pmin(1, seq_len(nrow(res)) * res$p)
  rownames(res) <- NULL
  res[, c("class", "distance", "pairs", "mantel_r", "p", "p_progressive")]
}
