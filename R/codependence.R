#' Multiscale codependence of two variables
#'
#' Splits the correlation between two variables by spatial or temporal scale
#' (multiscale codependence analysis; Guénard et al. 2010; Legendre & Legendre 2012,
#' §14.5.2). The positions give a complete orthonormal set of eigenfunctions of the
#' centred space (distance-based Moran's eigenvector maps, from the broadest to the
#' finest scale, completed with the null directions). For standardised `x` and `y`,
#' the correlation is exactly the sum over eigenfunctions of
#' `r(x, u_j) * r(y, u_j)`, the codependence coefficient of scale `j`.
#'
#' This gives the correction's estimand a spectral reading: the residual correlation
#' after projection on a basis made of some of these functions is the sum of the
#' coefficients outside the basis, renormalised by the residual variances. Applied to
#' the residuals of a correction, the spectrum shows at which scales the shared
#' structure that the basis did not capture remains.
#'
#' @param x,y Numeric vectors of the same length (one value per site or time).
#' @param coords Site coordinates (spatial), or
#' @param time Times (temporal; a regular series).
#' @param basis Optional complete orthonormal basis (n x (n - 1)), for instance from a
#'   previous call (`attr(result, "basis")`), to avoid recomputing it.
#' @return A data frame, one row per eigenfunction in order of scale: `rank`,
#'   `eigenvalue` (proportional to its Moran's I; 0 for the null directions),
#'   `r_x`, `r_y` (correlation of each variable with the function), `codependence`
#'   (their product) and `cumulative` (running sum; the last value is `cor(x, y)`).
#'   The basis is attached as `attr(, "basis")`.
#' @examples
#' set.seed(1)
#' tt <- 1:120
#' common <- sin(2 * pi * tt / 60)
#' x <- common + rnorm(120, sd = 0.5)
#' y <- common + rnorm(120, sd = 0.5)
#' cd <- codependence(x, y, time = tt)
#' head(cd[order(-abs(cd$codependence)), ], 3)   # the broad scale carries the association
#' @export
codependence <- function(x, y, coords = NULL, time = NULL, basis = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("`x` and `y` must have the same length.", call. = FALSE)
  if (is.null(basis)) {
    pos <- if (!is.null(coords)) as.matrix(coords) else if (!is.null(time)) matrix(as.numeric(time)) else
      stop("Give `coords` or `time`.", call. = FALSE)
    basis <- .complete_basis(pos)
  }
  U <- as.matrix(basis)
  unit <- function(v) { v <- v - mean(v); v / sqrt(sum(v^2)) }
  rx <- as.vector(crossprod(U, unit(x)))
  ry <- as.vector(crossprod(U, unit(y)))
  ev <- attr(basis, "eigenvalue")
  if (is.null(ev)) ev <- rep(NA_real_, ncol(U))
  out <- data.frame(rank = seq_len(ncol(U)), eigenvalue = ev, r_x = rx, r_y = ry, codependence = rx * ry,
                    cumulative = cumsum(rx * ry))
  structure(out, basis = basis)
}

# Complete orthonormal basis of the centred space from the positions: dbMEM with
# positive and negative eigenvalues, completed with the null directions, ordered from
# the broadest (largest eigenvalue) to the finest (most negative) scale.
.complete_basis <- function(pos) {
  m <- adespatial::dbmem(pos, MEM.autocor = "all", silent = TRUE)
  ev <- attr(m, "values")
  U <- as.matrix(m)
  U <- sweep(U, 2, sqrt(colSums(U^2)), "/")
  n <- nrow(U)
  Q <- qr.Q(qr(cbind(1, U)), complete = TRUE)
  extra <- if (ncol(U) + 1 < n) Q[, (ncol(U) + 2):n, drop = FALSE] else matrix(0, n, 0)
  B <- cbind(U, extra)
  e <- c(ev, rep(0, ncol(extra)))
  o <- order(-e)
  structure(B[, o, drop = FALSE], eigenvalue = e[o])
}
