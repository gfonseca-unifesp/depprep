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

# Moran correlogram (Legendre & Legendre 2012, §13.1.1): Moran's I in distance classes
# of equal frequency up to half the largest distance, each tested by permutation
# (two-sided: structure can be positive or negative at a given distance). The global
# p-value is the Bonferroni one, the number of classes times the smallest p-value
# (Oden 1984). It sees structure the nearest neighbours miss: patches smaller than the
# sampling interval, or structure at intermediate distances.
.correlogram_test <- function(coords, z, n_classes = 8, nsim = 499) {
  coords <- as.matrix(coords)
  D <- stats::dist(coords)
  d <- as.vector(D)
  br <- unique(stats::quantile(d[d <= max(d) / 2], seq(0, 1, length.out = n_classes + 1), names = FALSE))
  br[1] <- 0
  res <- t(vapply(seq_len(length(br) - 1), function(i) {
    # distance-class graphs are often split in sub-graphs; spdep warns, harmlessly
    nb <- suppressWarnings(spdep::dnearneigh(coords, br[i], br[i + 1], bounds = c("GT", "LE")))
    if (sum(spdep::card(nb)) == 0) return(c(I = NA_real_, p = NA_real_))
    lw <- suppressWarnings(spdep::nb2listw(nb, style = "B", zero.policy = TRUE))
    m <- spdep::moran.mc(z, lw, nsim = max(nsim, 1), zero.policy = TRUE, alternative = "two.sided")
    c(I = unname(m$statistic), p = m$p.value)
  }, numeric(2)))
  ok <- is.finite(res[, "p"])
  list(classes = data.frame(lower = utils::head(br, -1), upper = br[-1], moran_I = res[, "I"], p = res[, "p"]),
       p = if (any(ok)) min(1, sum(ok) * min(res[ok, "p"])) else NA_real_)
}

# Linear trend surface in the coordinates (Legendre & Legendre 2012, §13.1): its
# R-squared and the residuals. A gradient has no sill, so a variogram computed on it
# reads the trend as a long range or as anisotropy.
.spatial_trend <- function(coords, z) {
  X <- scale(as.matrix(coords), scale = FALSE)
  zc <- z - mean(z)
  e <- stats::lm.fit(X, zc)$residuals
  list(r2 = if (sum(zc^2) > 0) 1 - sum(e^2) / sum(zc^2) else 0, resid = e)
}

# Normal scores (ranks mapped to normal quantiles), used only for the diagnosis of
# skewed variables such as counts: kurtosis lowers the power of correlograms and
# variograms (Legendre & Legendre 2012, §13.1.1).
.normal_scores <- function(z) stats::qnorm((rank(z) - 0.5) / length(z))
.is_count <- function(z) all(is.finite(z)) && all(z >= 0) && all(abs(z - round(z)) < 1e-8) && length(unique(z)) > 2

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

# Temporal statistics of one regular series, following Legendre & Legendre (2012,
# ch. 12): the trend and the cycle are removed before the persistence is measured.
# 1. Period: least-squares periodogram (one sine-cosine pair fitted at each period of a
#    fine grid, Dutilleul's approach) on the linearly detrended series, so periods off
#    the Fourier grid (12.5 in 100 months, 52.18 weeks) are found. Periods up to a
#    third of the series only, so slow persistence is not read as a cycle. With
#    `period` given, the search is skipped. `period_strength` is the share of the
#    detrended variance explained by that pair; `period_p` tests the peak against the
#    local background spectrum (red noise has large low-frequency ordinates).
# 2. Trend: `trend_r2` is the partial R-squared of a linear trend given the cycle
#    (the cycle enters only when significant).
# 3. Persistence: `acf1`, `phi` and `ar_order` on the residuals of trend and cycle.
#    `ljung_p` stays on the raw series: it detects dependence of any kind.
.temporal_stats <- function(x, max_order = 10, period = NULL, alpha = 0.05) {
  n <- length(x)
  x <- x - mean(x)
  t <- seq_len(n)
  lb_p <- stats::Box.test(x, lag = min(10, floor(n / 5)), type = "Ljung-Box")$p.value
  tc <- t - mean(t)
  xd <- x - tc * sum(tc * x) / sum(tc^2)                          # linear trend removed
  pg <- .ls_periodogram(xd, t, period = period)
  harm <- if (is.finite(pg$period) && (!is.null(period) || pg$p < alpha)) {
    h <- seq_len(max(1, min(floor(pg$period / 2), 6)))
    m <- do.call(cbind, lapply(h, function(i) cbind(sin(2 * pi * i * t / pg$period), cos(2 * pi * i * t / pg$period))))
    scale(m, scale = FALSE)
  } else matrix(0, n, 0)
  rss <- function(M) if (ncol(M)) sum(stats::lm.fit(M, x)$residuals^2) else sum(x^2)
  rss_h <- rss(harm)
  full <- cbind(tc, harm)
  e <- stats::lm.fit(full, x)$residuals
  trend_r2 <- if (rss_h > 0) max(0, (rss_h - sum(e^2)) / rss_h) else 0
  acf1 <- stats::acf(e, 1, plot = FALSE)$acf[2]
  ar_order <- tryCatch(stats::ar(e, aic = TRUE, order.max = max_order, method = "yule-walker")$order,
                       error = function(err) NA_integer_)
  phi <- tryCatch(unname(stats::arima(e, order = c(1, 0, 0), include.mean = FALSE)$coef[1]),
                  error = function(err) acf1)
  c(acf1 = acf1, phi = phi, ar_order = ar_order, ljung_p = lb_p, period = pg$period,
    period_strength = pg$strength, period_p = pg$p, trend_r2 = trend_r2)
}

# Least-squares periodogram of a centred series: for every frequency of a grid
# oversampled `oversample` times relative to the Fourier spacing 1/n, the share of the
# variance explained by a sine-cosine pair (closed form of the two-column regression).
# The peak is tested against the local background spectrum, so red noise (large
# low-frequency ordinates) is not read as a cycle: under the null every ordinate is
# about s(f) * chi-squared(2), and the peak divided by the mean of the K ordinates
# around it follows about F(2, 2K). The p-value accounts for the frequencies searched
# (twice the Fourier frequencies of the band, for the oversampled grid).
.ls_periodogram <- function(x, t = seq_along(x), period = NULL, min_period = 2, max_period = length(x) / 3,
                            oversample = 10, window = 8) {
  n <- length(x)
  ss <- sum(x^2)
  na <- list(period = NA_real_, strength = NA_real_, p = NA_real_)
  if (n < 8 || ss <= 0) return(na)
  f <- if (is.null(period)) seq(1 / max_period, 1 / min_period, by = 1 / (oversample * n)) else 1 / period
  ang <- 2 * pi * outer(t, f)
  S <- sin(ang); C <- cos(ang)
  S <- sweep(S, 2, colMeans(S)); C <- sweep(C, 2, colMeans(C))
  sss <- colSums(S^2); scc <- colSums(C^2); ssc <- colSums(S * C)
  bs <- colSums(S * x); bc <- colSums(C * x)
  det <- sss * scc - ssc^2
  ess <- ifelse(det > 1e-8 * pmax(sss * scc, 1e-12), (scc * bs^2 - 2 * ssc * bs * bc + sss * bc^2) / det,
                ifelse(scc > 1e-12, bc^2 / scc, 0))
  k <- which.max(ess)
  # local background: mean ordinate at Fourier-spaced frequencies within `window` steps
  # of the peak (the peak and its immediate neighbours excluded), on the full grid
  fb <- seq(max(1 / max_period, f[k] - window / n), min(1 / min_period, f[k] + window / n), by = 1 / n)
  fb <- fb[abs(fb - f[k]) > 1.5 / n]
  ab <- 2 * pi * outer(t, fb)
  Sb <- sweep(sin(ab), 2, colMeans(sin(ab))); Cb <- sweep(cos(ab), 2, colMeans(cos(ab)))
  s1 <- colSums(Sb^2); c1 <- colSums(Cb^2); sc1 <- colSums(Sb * Cb); b1 <- colSums(Sb * x); b2 <- colSums(Cb * x)
  d1 <- s1 * c1 - sc1^2
  eb <- ifelse(d1 > 1e-8 * pmax(s1 * c1, 1e-12), (c1 * b1^2 - 2 * sc1 * b1 * b2 + s1 * b2^2) / d1,
               ifelse(c1 > 1e-12, b2^2 / c1, 0))
  K <- length(eb)
  if (K < 2 || mean(eb) <= 0) return(list(period = 1 / f[k], strength = ess[k] / ss, p = NA_real_))
  g <- ess[k] / mean(eb)                                   # about F(2, 2K) under the null
  p1 <- (1 + g / K)^(-K)
  m <- if (is.null(period)) max(1, 2 * (floor(n / min_period) - floor(n / max_period))) else 1
  list(period = 1 / f[k], strength = ess[k] / ss, p = 1 - (1 - p1)^m)
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
.st_stats <- function(z, coords, index, lw = NULL, k = 8, n_slices = 6, n_sites = 30, max_lag = 10, nsim = 0,
                      period = NULL) {
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
  tm <- t(vapply(sites, function(i) .temporal_stats(zm[i, ], period = period), numeric(8)))
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
    period_strength = med(tm[, "period_strength"]), period_p = med(tm[, "period_p"]), trend_r2 = med(tm[, "trend_r2"]),
    shift_indicator = max(gain_bin, na.rm = TRUE) - med(gain_bin),
    shift_direction = as.numeric(names(gain_bin)[which.max(gain_bin)]))
}
