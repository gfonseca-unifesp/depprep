#' Joint decorrelation by latent components (comparators for spatial data)
#'
#' Two baselines from the literature that decorrelate the variables jointly instead of
#' correcting them group by group: spatial blind source separation (SBSS,
#' `SpatialBSS::sbss`) and min/max autocorrelation factors (MAF). Both rotate the
#' standardised data into latent components ordered by spatial structure; components
#' with significant spatial autocorrelation (Moran's I, permutation test) are dropped
#' and the rest are rotated back to the variables. They target a different estimand
#' (uncorrelated latent sources), so the project compares them on the output methods
#' (protocol, comparators 6 and 7).
#'
#' * SBSS: ring kernel from 0 to `radius` (default: three times the median distance to
#'   the nearest neighbour).
#' * MAF: generalised eigenvectors of the covariance of the differences between
#'   k-nearest neighbours and the covariance of the data; the first factors are the
#'   most autocorrelated.
#'
#' @param X Numeric matrix (sites x variables).
#' @param coords Site coordinates.
#' @param method `"sbss"` or `"maf"`.
#' @param k Nearest neighbours (Moran's I and MAF).
#' @param radius Ring kernel radius for SBSS.
#' @param nsim Permutations of the Moran test.
#' @param alpha Significance level for dropping a component.
#' @param seed Seed of the permutations.
#' @return A list of class `depprep_decorrelated` with `X` (the decorrelated data), `R`
#'   (its correlation matrix), `removed` (number of components dropped) and `moran`
#'   (Moran's I and p-value of each component).
#' @examples
#' g <- sim_groups(sim_group("a", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(200, g, seed = 1)
#' decorrelate_latent(s$X, s$coords, method = "maf")$removed
#' @export
decorrelate_latent <- function(X, coords, method = c("sbss", "maf"), k = 8, radius = NULL, nsim = 199,
                               alpha = 0.05, seed = 1L) {
  method <- match.arg(method)
  Z <- scale(as.matrix(X))
  coords <- as.matrix(coords)
  lw <- .knn_weights(coords, k)
  if (method == "sbss") {
    if (is.null(radius)) {
      nn <- apply(as.matrix(stats::dist(coords)) + diag(Inf, nrow(coords)), 1, min)
      radius <- 3 * stats::median(nn)
    }
    fit <- SpatialBSS::sbss(Z, coords, kernel_type = "ring", kernel_parameters = c(0, radius))
    S <- as.matrix(fit$s)
    back <- function(keep) S[, keep, drop = FALSE] %*% t(fit$w_inv[, keep, drop = FALSE])
  } else {
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k))
    pairs <- do.call(rbind, lapply(seq_along(nb), function(i) cbind(i, nb[[i]])))
    D <- Z[pairs[, 1], , drop = FALSE] - Z[pairs[, 2], , drop = FALSE]
    C <- stats::cov(Z); Cd <- stats::cov(D)
    e <- eigen(solve(C, Cd))                               # small eigenvalues: most autocorrelated
    W <- Re(e$vectors[, order(Re(e$values)), drop = FALSE])
    S <- Z %*% W
    W_inv <- solve(W)
    back <- function(keep) S[, keep, drop = FALSE] %*% W_inv[keep, , drop = FALSE]
  }
  moran <- .with_stream(seed, .streams[["permutation"]],
                        t(vapply(seq_len(ncol(S)), function(j) .moran(S[, j], lw, nsim), numeric(2))))
  structured <- moran[, "moran_p"] < alpha & moran[, "moran_I"] > 0
  Xd <- back(!structured)
  dimnames(Xd) <- dimnames(Z)
  R <- suppressWarnings(stats::cor(Xd))
  R[!is.finite(R)] <- 0; diag(R) <- 1
  structure(list(X = Xd, R = R, removed = sum(structured), moran = moran, method = method),
            class = "depprep_decorrelated")
}
