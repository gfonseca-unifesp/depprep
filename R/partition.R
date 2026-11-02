#' Variation partition between position and drivers
#'
#' Splits the variance of each variable into the part explained by its position only
#' (the spatial, temporal or spatio-temporal basis), by the drivers only, by both
#' together (shared), and the residual, with adjusted R-squared (Peres-Neto et al.
#' 2006; the same fractions as `vegan::varpart()` with two tables). With `A` the
#' position basis, `B` the drivers and `AB` the two together, each with its adjusted
#' R-squared:
#' * position only = R2adj(AB) - R2adj(B);
#' * drivers only = R2adj(AB) - R2adj(A);
#' * shared = R2adj(A) + R2adj(B) - R2adj(AB);
#' * residual = 1 - R2adj(AB).
#'
#' The four add up to 1. The shared fraction is not a test and can be negative (when the
#' two tables explain more together than the sum of their parts). When the position
#' basis was selected on the same data (as the bases of [build_group_basis()] are), its
#' adjusted R-squared is somewhat optimistic (Blanchet et al. 2008), and so is the
#' position-only fraction.
#'
#' @param X Numeric matrix (observations x variables). It may hold missing values
#'   (`NA`; decision D26): the fractions of a variable are then computed on its observed
#'   rows (the tables restricted to those rows, `n` the number of them).
#' @param position The position basis: one matrix with one row per observation (for
#'   example a basis from [build_group_basis()]) for every variable, or a list with one
#'   such matrix per variable (named by the columns of `X`, or in their order). `NULL`
#'   or a matrix with no column for no position term.
#' @param drivers Numeric matrix or data frame of drivers, one row per observation of
#'   `X` (in [correct()], the temporal drivers repeated over the sites).
#' @return A data frame with one row per variable: `variable`, `k_position` and
#'   `k_drivers` (number of independent columns of each table), and the fractions
#'   `position_only`, `shared`, `drivers_only` and `residual`.
#' @references Blanchet, F. G., Legendre, P. and Borcard, D. (2008). Forward selection
#'   of explanatory variables. Ecology 89, 2623-2632.
#'
#'   Peres-Neto, P. R., Legendre, P., Dray, S. and Borcard, D. (2006). Variation
#'   partitioning of species data matrices: estimation and comparison of fractions.
#'   Ecology 87, 2614-2625.
#' @examples
#' set.seed(1)
#' t <- 1:120
#' position <- cbind(sin(2 * pi * t / 12), cos(2 * pi * t / 12))   # a seasonal basis
#' driver <- cbind(index = stats::rnorm(120))
#' x <- 2 * position[, 1] + driver[, 1] + stats::rnorm(120)
#' variation_partition(cbind(x = x), position, driver)
#' @export
variation_partition <- function(X, position, drivers) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  n <- nrow(X)
  D <- .driver_matrix(drivers, n)
  if (is.null(position) || !is.list(position) || is.data.frame(position)) {
    position <- rep(list(position), ncol(X))
  } else if (!is.null(names(position))) {
    if (!all(colnames(X) %in% names(position))) stop("`position` must name the basis of every column of `X`.", call. = FALSE)
    position <- position[colnames(X)]
  } else if (length(position) != ncol(X)) {
    stop("`position` must be one basis or one basis per column of `X`.", call. = FALSE)
  }
  QD <- .orthonormal(D)
  rows <- lapply(seq_len(ncol(X)), function(j) {
    P <- if (is.null(position[[j]])) matrix(0, n, 0) else as.matrix(unclass(position[[j]]))
    if (nrow(P) != n) stop("The position basis must have one row per observation of `X`.", call. = FALSE)
    QP <- if (ncol(P)) .orthonormal(P) else P
    QA <- if (ncol(P)) .orthonormal(cbind(P, D)) else QD
    x <- X[, j, drop = FALSE]
    a <- .adj_r2(x, QP); b <- .adj_r2(x, QD); ab <- .adj_r2(x, QA)
    data.frame(variable = colnames(X)[j], k_position = ncol(QP), k_drivers = ncol(QD),
               position_only = ab - b, shared = a + b - ab, drivers_only = ab - a, residual = 1 - ab)
  })
  do.call(rbind, rows)
}

# Adjusted R-squared of the column x on the orthonormal (centred) columns of Q. With
# missing values, on the observed rows, where Q is no longer orthonormal: a regression
# with the constant, k the rank of the table there.
.adj_r2 <- function(x, Q) {
  k <- ncol(Q); n <- nrow(Q)
  if (!k) return(0)
  if (anyNA(x)) {
    o <- !is.na(x); x <- as.vector(x)[o]; n <- sum(o)
    q <- qr(cbind(1, Q[o, , drop = FALSE]), tol = 1e-7)
    k <- q$rank - 1
    if (k >= n - 1) stop("A table of the partition has as many columns as observed values.", call. = FALSE)
    r2 <- 1 - sum(qr.resid(q, x)^2) / sum((x - mean(x))^2)
    return(1 - (1 - r2) * (n - 1) / (n - k - 1))
  }
  if (k >= n - 1) stop("A table of the partition has as many columns as observations.", call. = FALSE)
  x <- scale(x, scale = FALSE)
  r2 <- sum(crossprod(Q, x)^2) / sum(x^2)
  1 - (1 - r2) * (n - 1) / (n - k - 1)
}

# Drivers as a numeric matrix with `n` rows and no missing value.
.driver_matrix <- function(drivers, n) {
  D <- if (is.data.frame(drivers)) data.matrix(drivers) else as.matrix(drivers)
  if (!is.numeric(D) || !ncol(D)) stop("`drivers` must be a numeric matrix or data frame with at least one column.", call. = FALSE)
  if (nrow(D) != n) stop("`drivers` must have one row per observation.", call. = FALSE)
  if (any(!is.finite(D))) stop("`drivers` has missing or infinite values.", call. = FALSE)
  if (is.null(colnames(D))) colnames(D) <- paste0("driver", seq_len(ncol(D)))
  D
}
