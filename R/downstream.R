#' Run an exploratory method and score it against the truth
#'
#' Applies PCA, MDS, k-means or a self-organizing map (SOM) to a data set (raw or
#' corrected) and returns the metrics of the project's protocol, comparing the result
#' with the truth of a simulation:
#' * `"pca"` (on the correlation matrix `R`, or on `cor(X)`): distance between the
#'   estimated and the true leading subspace (eigenvectors of `R_true` with eigenvalue
#'   above 1, see [basis_distance()]); number of components above the broken-stick
#'   expectation; whether the first eigenvalue exceeds it (a spurious axis when
#'   `R_true` is the identity);
#' * `"mds"` (classical MDS of the observations, Euclidean or Gower distances):
#'   Procrustes correlation of the 2-D configuration with the configuration of the true
#'   signal; stress; Procrustes correlation with the positions (a configuration that
#'   only reproduces geography or time is spurious);
#' * `"kmeans"`: number of groups chosen by the gap statistic (`cluster::clusGap`), the
#'   adjusted Rand index (ARI) with the true classes for that number and for the true
#'   number of classes;
#' * `"som"` (5 x 5 hexagonal map; units clustered into the true number of classes, or
#'   3): ARI, quantisation error, topographic error and the canonical correlation
#'   between the position on the map and the position in space or time.
#'
#' @param X Numeric matrix (observations x variables).
#' @param method `"pca"`, `"mds"`, `"kmeans"` or `"som"`.
#' @param truth List with any of: `R_true` (p x p), `classes` (one label per row of
#'   `X`), `signal` (matrix, the true signal of each observation, for MDS) and
#'   `position` (coordinates or times of the rows, for the spurious-pattern metrics).
#' @param R Correlation matrix for PCA (`cor(X)` when `NULL`).
#' @param distance `"euclidean"` or `"gower"` (MDS).
#' @param k_max Largest number of groups for the gap statistic.
#' @param B_gap Reference data sets for the gap statistic.
#' @param seed Seed of k-means, the gap statistic and the SOM.
#' @return A data frame with `method`, `metric` and `value`.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(150, 0), 50), matrix(rnorm(150, 3), 50))
#' run_downstream(X, "kmeans", truth = list(classes = rep(1:2, each = 50)))
#' @export
run_downstream <- function(X, method = c("pca", "mds", "kmeans", "som"), truth = list(), R = NULL,
                           distance = c("euclidean", "gower"), k_max = 6, B_gap = 20, seed = 1L) {
  method <- match.arg(method)
  distance <- match.arg(distance)
  X <- as.matrix(X)
  out <- function(...) { v <- c(...); data.frame(method = if (method == "mds") paste0("mds_", distance) else method,
                                                  metric = names(v), value = unname(v), stringsAsFactors = FALSE) }
  if (method == "pca") {
    R <- if (is.null(R)) stats::cor(X) else R
    p <- ncol(R)
    ev <- eigen(R, symmetric = TRUE)
    bs <- rev(cumsum(1 / rev(seq_len(p))))                  # broken-stick: sum_{i=k}^{p} 1/i
    above <- ev$values > bs
    n_bs <- if (above[1]) which.min(c(above, FALSE)) - 1 else 0
    dist_true <- NA_real_
    if (!is.null(truth$R_true)) {
      et <- eigen(truth$R_true, symmetric = TRUE)
      q <- sum(et$values > 1 + 1e-8)
      if (q > 0) dist_true <- basis_distance(structure(ev$vectors[, seq_len(q), drop = FALSE], class = "depprep_basis"),
                                             structure(et$vectors[, seq_len(q), drop = FALSE], class = "depprep_basis"))
    }
    return(out(subspace_distance = dist_true, components_above_broken_stick = n_bs,
               first_above_broken_stick = as.numeric(above[1]), first_eigen_share = ev$values[1] / p))
  }
  Z <- scale(X)
  Z[!is.finite(Z)] <- 0
  .with_stream(seed, .streams[["downstream"]], {
    if (method == "mds") {
      D <- if (distance == "euclidean") stats::dist(Z) else cluster::daisy(as.data.frame(Z), metric = "gower")
      conf <- stats::cmdscale(D, k = 2)
      d_conf <- stats::dist(conf)
      stress <- sqrt(sum((as.vector(D) - as.vector(d_conf))^2) / sum(as.vector(D)^2))
      truth_r <- if (!is.null(truth$signal)) .procrustes_r(conf, stats::cmdscale(stats::dist(scale(truth$signal)), k = 2)) else NA_real_
      pos_r <- if (!is.null(truth$position)) .procrustes_r(conf, as.matrix(truth$position)) else NA_real_
      return(out(procrustes_truth = truth_r, stress = stress, procrustes_position = pos_r))
    }
    k_true <- if (!is.null(truth$classes)) length(unique(truth$classes)) else NA
    if (method == "kmeans") {
      km <- function(x, k) stats::kmeans(x, k, nstart = 10, iter.max = 50)
      gap <- cluster::clusGap(Z, km, K.max = k_max, B = B_gap, verbose = FALSE)
      k_gap <- cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
      ari <- function(k) if (is.na(k) || k < 1) NA_real_ else if (k == 1) mclust::adjustedRandIndex(rep(1, nrow(Z)), truth$classes) else
        mclust::adjustedRandIndex(km(Z, k)$cluster, truth$classes)
      return(out(k_gap = k_gap, k_gap_above_1 = as.numeric(k_gap > 1),
                 ari_gap = if (is.na(k_true)) NA_real_ else ari(k_gap), ari_true_k = if (is.na(k_true)) NA_real_ else ari(k_true)))
    }
    # SOM
    grid <- kohonen::somgrid(5, 5, "hexagonal")
    sm <- kohonen::som(Z, grid = grid, rlen = 100)
    codes <- sm$codes[[1]]
    d_all <- as.matrix(stats::dist(rbind(codes, Z)))[-seq_len(nrow(codes)), seq_len(nrow(codes)), drop = FALSE]
    o <- t(apply(d_all, 1, order))[, 1:2, drop = FALSE]
    qe <- mean(d_all[cbind(seq_len(nrow(Z)), o[, 1])])
    ud <- kohonen::unit.distances(grid)
    te <- mean(ud[o] > 1.0001)
    k_som <- if (is.na(k_true)) 3 else k_true
    unit_class <- stats::cutree(stats::hclust(stats::dist(codes), method = "ward.D2"), k_som)
    obs_class <- unit_class[sm$unit.classif]
    pos_cor <- if (!is.null(truth$position)) {
      pm <- grid$pts[sm$unit.classif, , drop = FALSE]
      tryCatch(stats::cancor(pm, as.matrix(truth$position))$cor[1], error = function(e) NA_real_)
    } else NA_real_
    out(ari = if (is.na(k_true)) NA_real_ else mclust::adjustedRandIndex(obs_class, truth$classes),
        quantisation_error = qe, topographic_error = te, map_position_cor = pos_cor)
  })
}

# Symmetric Procrustes correlation between two configurations (same rows): both are
# centred and scaled to unit size; the correlation is the sum of the singular values of
# their cross-product (1 = identical up to rotation, reflection and scale). A
# configuration with fewer columns is padded with zeros.
.procrustes_r <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  k <- max(ncol(A), ncol(B))
  pad <- function(M) cbind(M, matrix(0, nrow(M), k - ncol(M)))
  A <- pad(scale(A, scale = FALSE)); B <- pad(scale(B, scale = FALSE))
  A <- A / sqrt(sum(A^2)); B <- B / sqrt(sum(B^2))
  sum(svd(crossprod(A, B))$d)
}
