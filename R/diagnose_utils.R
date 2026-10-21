# Building blocks of diagnose_variables(). Every function here looks at ONE variable
# at a time: dependence signatures use only marginal information (variogram, ACF),
# never correlations between variables.

# --- Spatial ---------------------------------------------------------------------------

# Row-standardised weights from the k nearest neighbours (made symmetric).
.knn_weights <- function(coords, k = 8) {
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k))
  spdep::nb2listw(spdep::make.sym.nb(nb), style = "W")
}

# Moran's I with its one-sided p-value (positive autocorrelation). With nsim > 0 the
# p-value comes from a permutation test (exact; uses the current random stream). The
# normal approximation (nsim = 0) is optimistic in small samples: in stage 3 it gave
# 5-13% false positives at a nominal 5%.
.moran <- function(z, lw, nsim = 0) {
  if (nsim > 0) {
    t <- spdep::moran.mc(z, lw, nsim = nsim, alternative = "greater")
    return(c(moran_I = unname(t$statistic), moran_p = t$p.value))
  }
  t <- spdep::moran.test(z, lw, randomisation = TRUE, alternative = "greater")
  c(moran_I = unname(t$estimate[["Moran I statistic"]]), moran_p = t$p.value)
}

# gstat::fit.variogram() prints advice when a fit does not converge; those cases are
# handled as missing values, so the printout is silenced.
.quiet_fit <- function(v, model) {
  out <- NULL
  utils::capture.output(out <- tryCatch(suppressWarnings(gstat::fit.variogram(v, model)),
                                        error = function(e) NULL))
  out
}

# Practical range of a fitted model: 3a for the exponential (correlation 0.05), a for
# the spherical (correlation exactly 0).
.practical_range <- function(model, a) ifelse(model == "Exp", 3 * a, a)

# Fit exponential and spherical variograms with a nugget; keep the one with the
# smaller weighted sum of squared errors. Returns NA values when no fit converges.
.fit_variogram <- function(coords, z, cutoff = NULL, models = c("Exp", "Sph")) {
  d <- data.frame(as.matrix(coords), z = z)
  names(d)[1:2] <- c("x", "y")
  if (is.null(cutoff)) cutoff <- 0.5 * sqrt(sum(apply(d[, 1:2], 2, function(v) diff(range(v)))^2))
  v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = d, cutoff = cutoff)
  vz <- stats::var(z)
  best <- NULL
  for (m in models) {
    init <- gstat::vgm(0.7 * vz, m, cutoff / 3, nugget = 0.3 * vz)
    fit <- .quiet_fit(v, init)
    if (is.null(fit) || isTRUE(attr(fit, "singular"))) next
    a <- fit$range[fit$model == m]
    if (!length(a) || !is.finite(a) || a <= 0 || any(fit$psill < 0)) next
    sse <- attr(fit, "SSErr")
    if (is.null(best) || sse < best$sse) best <- list(fit = fit, model = m, sse = sse)
  }
  if (is.null(best)) {
    return(list(vario_model = NA_character_, nugget = NA_real_, sill = NA_real_, range = NA_real_,
                nugget_ratio = NA_real_, fit = NULL))
  }
  f <- best$fit
  nug <- sum(f$psill[f$model == "Nug"])
  sill <- sum(f$psill)
  list(vario_model = best$model, nugget = nug, sill = sill,
       range = .practical_range(best$model, f$range[f$model == best$model]),
       nugget_ratio = if (sill > 0) nug / sill else NA_real_, fit = f)
}

# Anisotropy from directional variograms, without fitting one model per direction
# (those fits are unstable). The nugget and range come from the omnidirectional fit
# `fit` (output of .fit_variogram()). In each direction, the structural semivariance
# (semivariance minus nugget) is averaged over the short lags, up to half the range,
# where the variogram rises roughly in proportion to distance / range: the direction
# with the slowest rise has the longest range. The ratio between the largest and the
# smallest rise estimates the anisotropy ratio. gstat angles are clockwise from north
# (the y axis); `aniso_direction` is the direction of the longest range.
.anisotropy <- function(coords, z, fit, angles = c(0, 45, 90, 135)) {
  if (is.null(fit$fit) || !is.finite(fit$range)) return(c(aniso_ratio = NA_real_, aniso_direction = NA_real_))
  d <- data.frame(as.matrix(coords), z = z)
  names(d)[1:2] <- c("x", "y")
  short <- fit$range / 2
  v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = d, cutoff = short, width = short / 5,
                        alpha = angles, tol.hor = 180 / length(angles) / 2)
  if (is.null(v) || !NROW(v)) return(c(aniso_ratio = NA_real_, aniso_direction = NA_real_))
  rise <- vapply(angles, function(a) {
    va <- v[v$dir.hor == a & v$np > 0, ]
    if (!NROW(va)) return(NA_real_)
    stats::weighted.mean(pmax(va$gamma - fit$nugget, 0) / va$dist, va$np)
  }, numeric(1))
  ok <- is.finite(rise) & rise > 0
  if (sum(ok) < 2) return(c(aniso_ratio = NA_real_, aniso_direction = NA_real_))
  c(aniso_ratio = max(rise[ok]) / min(rise[ok]), aniso_direction = angles[ok][which.min(rise[ok])])
}

# --- Temporal --------------------------------------------------------------------------

# Temporal statistics of one regular series. The seasonal period is searched among
# periods up to a third of the series, so that slow persistence or trends (whose
# peak sits at the lowest frequencies) are not read as seasonality.
.temporal_stats <- function(x, max_order = 10) {
  n <- length(x)
  x <- x - mean(x)
  t <- seq_len(n)
  acf1 <- stats::acf(x, 1, plot = FALSE)$acf[2]
  ar_order <- tryCatch(stats::ar(x, aic = TRUE, order.max = max_order, method = "yule-walker")$order,
                       error = function(e) NA_integer_)
  phi <- tryCatch(unname(stats::arima(x, order = c(1, 0, 0), include.mean = FALSE)$coef[1]),
                  error = function(e) acf1)
  lb_p <- stats::Box.test(x, lag = min(10, floor(n / 5)), type = "Ljung-Box")$p.value
  sp <- stats::spec.pgram(x, taper = 0, detrend = FALSE, demean = TRUE, fast = FALSE, plot = FALSE)
  cand <- which(1 / sp$freq <= n / 3)
  k <- cand[which.max(sp$spec[cand])]
  c(acf1 = acf1, phi = phi, ar_order = ar_order, ljung_p = lb_p,
    period = 1 / sp$freq[k], period_strength = sp$spec[k] / sum(sp$spec),
    trend_r2 = summary(stats::lm(x ~ t))$r.squared)
}

# --- Spatio-temporal -------------------------------------------------------------------

# Statistics of one spatio-temporal variable `z` (rows ordered as in
# sim_spatiotemporal(): sites vary fastest within each time).
# - spatial: Moran and variogram on `n_slices` time slices. Moran's p-value is the
#   smallest one times the number of slices (Bonferroni), because a seasonal amplitude
#   map is invisible in slices where the cycle crosses zero; the range is the median
#   over the slices with structure;
# - temporal: .temporal_stats() on up to `n_sites` sites (medians);
# - movement, between neighbours i, j of the same variable: |cor(z_i(t + L), z_j(t))|
#   for lags L = 1..max_lag against the same-time |cor(z_i, z_j)|. A pattern that
#   moves brings the upstream neighbour's past to the site, but only from one
#   direction; noise and separable fields give the same (biased) gain in every
#   direction. shift_indicator is the gain of the best direction minus the median
#   gain over directions; shift_direction is that direction (degrees, from the
#   neighbour to the site, counter-clockwise from the x axis).
.st_stats <- function(z, coords, index, lw = NULL, k = 8, n_slices = 6, n_sites = 30, max_lag = 10, nsim = 0) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  nt <- max(index$time)
  zm <- matrix(NA_real_, n, nt)
  zm[cbind(index$site, index$time)] <- z
  if (is.null(lw)) lw <- .knn_weights(coords, k)
  med <- function(x) stats::median(x, na.rm = TRUE)
  slices <- unique(round(seq(1, nt, length.out = n_slices)))
  sp <- t(vapply(slices, function(s) {
    m <- .moran(zm[, s], lw, nsim)
    f <- .fit_variogram(coords, zm[, s])
    c(m, range = f$range, nugget_ratio = f$nugget_ratio)
  }, numeric(4)))
  structured <- sp[, "moran_p"] < 0.05
  use <- if (any(structured)) structured else rep(TRUE, nrow(sp))
  sites <- unique(round(seq(1, n, length.out = min(n, n_sites))))
  tm <- t(vapply(sites, function(i) .temporal_stats(zm[i, ]), numeric(7)))
  nb <- spdep::knearneigh(coords, k = k)$nn
  pairs <- data.frame(i = rep(sites, each = k), j = as.vector(t(nb[sites, , drop = FALSE])))
  lags <- seq_len(max(1, min(max_lag, floor(nt / 4))))
  lagcor <- function(x, y, L) stats::cor(x[(L + 1):nt], y[1:(nt - L)])  # x(t + L) vs y(t)
  r_same <- abs(mapply(function(i, j) stats::cor(zm[i, ], zm[j, ]), pairs$i, pairs$j))
  r_lag <- abs(matrix(unlist(mapply(function(i, j) vapply(lags, function(L) lagcor(zm[i, ], zm[j, ], L), 0),
                                    pairs$i, pairs$j, SIMPLIFY = FALSE)), ncol = length(lags), byrow = TRUE))
  gain <- apply(r_lag, 1, max, na.rm = TRUE) - r_same
  ang <- atan2(coords[pairs$i, 2] - coords[pairs$j, 2], coords[pairs$i, 1] - coords[pairs$j, 1]) * 180 / pi
  bin <- (round(ang / 45) %% 8) * 45
  gain_bin <- tapply(gain, bin, mean, na.rm = TRUE)
  c(moran_I = med(sp[use, "moran_I"]), moran_p = min(1, min(sp[, "moran_p"]) * nrow(sp)),
    range = med(sp[use, "range"]), nugget_ratio = med(sp[use, "nugget_ratio"]),
    acf1 = med(tm[, "acf1"]), phi = med(tm[, "phi"]), ar_order = med(tm[, "ar_order"]),
    ljung_p = med(tm[, "ljung_p"]),
    period = as.numeric(names(which.max(table(round(tm[, "period"]))))),
    period_strength = med(tm[, "period_strength"]), trend_r2 = med(tm[, "trend_r2"]),
    shift_indicator = max(gain_bin, na.rm = TRUE) - med(gain_bin),
    shift_direction = as.numeric(names(gain_bin)[which.max(gain_bin)]))
}
