# Building blocks of diagnose_variables(). Every function here looks at ONE variable
# at a time: dependence signatures use only marginal information (variogram, ACF),
# never correlations between variables.

# --- Spatial ---------------------------------------------------------------------------

# Row-standardised weights from the k nearest neighbours (made symmetric).
.knn_weights <- function(coords, k = 8) {
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k))
  spdep::nb2listw(spdep::make.sym.nb(nb), style = "W")
}

# Moran's I with its p-value (randomisation, one-sided: positive autocorrelation).
.moran <- function(z, lw) {
  t <- spdep::moran.test(z, lw, randomisation = TRUE, alternative = "greater")
  c(moran_I = unname(t$estimate[["Moran I statistic"]]), moran_p = t$p.value)
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
    fit <- tryCatch(suppressWarnings(gstat::fit.variogram(v, init)), error = function(e) NULL)
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

# Anisotropy from directional variograms. gstat angles are clockwise from north
# (the y axis). Returns the ratio between the largest and smallest fitted range and
# the direction (gstat angle) of the largest one.
.anisotropy <- function(coords, z, angles = c(0, 45, 90, 135), cutoff = NULL) {
  d <- data.frame(as.matrix(coords), z = z)
  names(d)[1:2] <- c("x", "y")
  if (is.null(cutoff)) cutoff <- 0.5 * sqrt(sum(apply(d[, 1:2], 2, function(v) diff(range(v)))^2))
  v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = d, cutoff = cutoff,
                        alpha = angles, tol.hor = 180 / length(angles) / 2)
  vz <- stats::var(z)
  rng <- vapply(angles, function(a) {
    va <- v[v$dir.hor == a, ]
    if (nrow(va) < 4) return(NA_real_)
    fit <- tryCatch(suppressWarnings(gstat::fit.variogram(va, gstat::vgm(0.7 * vz, "Exp", cutoff / 3, nugget = 0.3 * vz))),
                    error = function(e) NULL)
    if (is.null(fit) || isTRUE(attr(fit, "singular"))) return(NA_real_)
    r <- 3 * fit$range[fit$model == "Exp"]
    if (!is.finite(r) || r <= 0) NA_real_ else min(r, 10 * cutoff)  # cap runaway fits
  }, numeric(1))
  ok <- is.finite(rng)
  if (sum(ok) < 2) return(c(aniso_ratio = NA_real_, aniso_direction = NA_real_))
  c(aniso_ratio = max(rng[ok]) / min(rng[ok]), aniso_direction = angles[ok][which.max(rng[ok])])
}
