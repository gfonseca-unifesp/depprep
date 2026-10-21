#' Diagnose the dependence of each variable
#'
#' Computes a dependence signature for every column of `X`, using only that column
#' (its variogram, autocorrelation and periodogram). Correlations between variables
#' are never used, so the signatures can later choose a correction without looking
#' at the associations the correction is meant to estimate.
#'
#' The family is set by the arguments: `coords` only gives a spatial diagnosis,
#' `time` only a temporal one, and `coords` with `index` (site and time of each row,
#' as returned by [sim_spatiotemporal()]) a spatio-temporal one. A `depprep_sim`
#' object can be passed directly as `X`.
#'
#' Columns of the result (`NA` when not applicable to the family):
#' * spatial: `moran_I`, `moran_p` (k nearest neighbours), `vario_model` (`"Exp"` or
#'   `"Sph"`), `nugget`, `sill`, `range` (practical range), `nugget_ratio`,
#'   `aniso_ratio` and `aniso_direction` (directional variograms; gstat angle,
#'   clockwise from north);
#' * temporal: `acf1`, `phi` (AR(1) coefficient), `ar_order` (AIC), `ljung_p`,
#'   `period` and `period_strength` (dominant seasonal period and its share of the
#'   variance), `trend_r2`;
#' * spatio-temporal: the spatial statistics on time slices, the temporal ones on
#'   sites, `shift_indicator` and `shift_direction` (a pattern that moves);
#' * with `B > 0`: parametric bootstrap of `range` (spatial) and `phi` (temporal):
#'   `range_se`, `range_lo`, `range_hi` (percentile 90% interval), `phi_se`, `phi_lo`,
#'   `phi_hi` (basic 90% interval, which corrects the bias of the AR estimate).
#'
#' @param X Numeric matrix (observations x variables) or a `depprep_sim` object.
#' @param coords Two-column matrix of site coordinates (one row per site).
#' @param time Time of each row (temporal family); rows are sorted by it.
#' @param index Data frame with `site` and `time` of each row (spatio-temporal family).
#' @param B Number of parametric bootstrap replicates (0 = none). Not available for
#'   the spatio-temporal family.
#' @param seed Seed for the bootstrap (required when `B > 0`).
#' @param k Number of nearest neighbours for Moran's I.
#' @return A data frame of class `depprep_signature`, one row per variable.
#' @examples
#' g <- sim_groups(sim_group("long", 2, "exponential", range = 0.5),
#'                 sim_group("none", 2, "none"))
#' s <- sim_spatial(200, g, seed = 1)
#' diagnose_variables(s)[, c("variable", "moran_p", "range", "nugget_ratio")]
#' @export
diagnose_variables <- function(X, coords = NULL, time = NULL, index = NULL, B = 0, seed = NULL, k = 8) {
  if (inherits(X, "depprep_sim")) {
    sim <- X
    X <- sim$X
    if (sim$family %in% c("spatial", "spatiotemporal")) coords <- sim$coords
    if (sim$family == "temporal") time <- sim$time
    if (sim$family == "spatiotemporal") index <- sim$index
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be numeric.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else
    if (!is.null(coords)) "spatial" else if (!is.null(time)) "temporal" else
      stop("Give `coords` (spatial), `time` (temporal) or `coords` and `index` (spatio-temporal).", call. = FALSE)
  if (B > 0 && is.null(seed)) stop("`seed` is required when `B > 0`.", call. = FALSE)
  if (B > 0 && family == "spatiotemporal") {
    warning("The bootstrap is not available for spatio-temporal data; `B` is ignored.", call. = FALSE)
    B <- 0
  }

  if (family == "temporal") {
    ord <- order(time)
    X <- X[ord, , drop = FALSE]
  }
  lw <- if (family != "temporal") .knn_weights(coords, k) else NULL   # depends on the sites only

  rows <- lapply(seq_len(ncol(X)), function(j) {
    z <- X[, j]
    out <- .empty_signature()
    if (family == "spatial") {
      m <- .moran(z, lw)
      f <- .fit_variogram(coords, z)
      out[c("moran_I", "moran_p")] <- as.list(m)
      out[c("vario_model", "nugget", "sill", "range", "nugget_ratio")] <- f[c("vario_model", "nugget", "sill", "range", "nugget_ratio")]
      out[c("aniso_ratio", "aniso_direction")] <- as.list(.anisotropy(coords, z, f))
      if (B > 0 && !is.null(f$fit)) {
        out[c("range_se", "range_lo", "range_hi")] <- as.list(.with_stream(seed + j, .streams[["bootstrap"]],
                                                                      .boot_range(coords, f$fit, B)))
      }
    } else if (family == "temporal") {
      ts <- .temporal_stats(z)
      out[names(ts)] <- as.list(ts)
      if (B > 0) {
        out[c("phi_se", "phi_lo", "phi_hi")] <- as.list(.with_stream(seed + j, .streams[["bootstrap"]],
                                                                .boot_phi(length(z), ts[["phi"]], B)))
      }
    } else {
      st <- .st_stats(z, coords, index, lw = lw, k = k)
      out[names(st)] <- as.list(st)
    }
    out$variable <- colnames(X)[j]
    out
  })
  res <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  res$family <- family
  res <- res[, c("variable", "family", setdiff(names(res), c("variable", "family")))]
  structure(res, class = c("depprep_signature", "data.frame"), B = B, seed = seed)
}

# One row with every signature column, all missing.
.empty_signature <- function() {
  num <- c("moran_I", "moran_p", "nugget", "sill", "range", "nugget_ratio", "aniso_ratio", "aniso_direction",
           "acf1", "phi", "ar_order", "ljung_p", "period", "period_strength", "trend_r2",
           "shift_indicator", "shift_direction",
           "range_se", "range_lo", "range_hi", "phi_se", "phi_lo", "phi_hi")
  c(list(variable = NA_character_, vario_model = NA_character_),
    stats::setNames(as.list(rep(NA_real_, length(num))), num))
}

# Parametric bootstrap of the practical range: simulate the fitted model at the same
# sites, refit, and summarise the refitted ranges (standard error, 90% interval).
.boot_range <- function(coords, fit, B) {
  d <- as.matrix(stats::dist(as.matrix(coords)))
  cmat <- gstat::variogramLine(fit, dist_vector = d, covariance = TRUE)
  u <- chol(cmat + diag(1e-8, nrow(cmat)))
  sims <- t(u) %*% matrix(stats::rnorm(nrow(d) * B), nrow(d))
  r <- apply(sims, 2, function(z) .fit_variogram(coords, z)$range)
  .boot_summary(r)
}

# Parametric bootstrap of the AR(1) coefficient.
.boot_phi <- function(n_times, phi, B) {
  phi <- max(min(phi, 0.99), -0.99)
  sims <- .ar1(n_times, phi, n_series = B)
  r <- apply(sims, 2, function(x) .temporal_stats(x)[["phi"]])
  .boot_summary(r, estimate = phi, type = "basic")
}

# Standard error and 90% interval. "percentile" uses the bootstrap quantiles; "basic"
# reflects them around the estimate (2 * estimate - quantiles), which removes the
# estimator's bias (AR coefficients are biased towards zero in short series).
.boot_summary <- function(r, estimate = NULL, type = c("percentile", "basic")) {
  type <- match.arg(type)
  r <- r[is.finite(r)]
  if (length(r) < 2) return(c(se = NA_real_, lo = NA_real_, hi = NA_real_))
  q <- stats::quantile(r, c(0.05, 0.95), names = FALSE)
  if (type == "basic") q <- rev(2 * estimate - q)
  c(se = stats::sd(r), lo = q[1], hi = q[2])
}
