#' Run an exploratory method and score it against the truth
#'
#' Applies PCA, MDS, k-means, a self-organizing map (SOM) or contiguity-constrained
#' clustering to a data set (raw or corrected) and returns the metrics of the project's
#' protocol, comparing the result with the truth of a simulation:
#' * `"pca"` (on the correlation matrix `R`, or on `cor(X)`): distance between the
#'   estimated and the true leading subspace (eigenvectors of `R_true` with eigenvalue
#'   above 1, see [basis_distance()]); the same distance to the leading subspace of the
#'   correlation of the true signal (`truth$signal`: residual plus class effects, as many
#'   axes as pass the broken-stick; `subspace_distance_signal`), because the classes add
#'   structure shared by all variables that `R_true` does not hold; number of components above the broken-stick
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
#'   between the position on the map and the position in space or time;
#' * `"hdbscan"` (density-based clustering, `dbscan::hdbscan()`; decision D18): no
#'   distance threshold to tune, which would be unfair across corrections that change the
#'   scale of the distances. The minimum cluster size is `max(10, ceiling(0.02 n))`
#'   observations (a count, not a distance): micro-clusters below 2% of the data are not
#'   groups, while the simulated classes hold about a third each. Returns the number of
#'   clusters without the noise points (`n_clusters`), whether there is more than one
#'   (`clusters_above_1`, the spurious pattern of the null scenario, as `k > 1` by the gap
#'   statistic for k-means), the share of noise points (`noise_share`) and the ARI with
#'   the true classes, the noise points counted as a group of their own (`ari`); also the
#'   ARI on the clustered points only (`ari_clustered`, descriptive): HDBSCAN leaves many
#'   points as noise even around clean classes, so `ari` mixes purity and coverage.
#' * `"constrained"` (contiguity-constrained Ward clustering, `adespatial::constr.hclust`;
#'   Legendre & Legendre 2012, sections 12.6 and 13.3): only observations linked in
#'   space or time may merge. The links come from `truth$position` (see below); the tree
#'   is cut at the true number of classes, or 3. Returns that number `k`, the ARI with
#'   the true classes (`ari_true_k`, the same gap-free metric as for k-means) and the
#'   Rand index between the constrained partition and the k-means partition with the
#'   same `k` (`rand_kmeans`). It is the reference method when the classes sought are
#'   regions or periods, which any decorrelation removes.
#'
#' @param X Numeric matrix (observations x variables).
#' @param method `"pca"`, `"mds"`, `"kmeans"`, `"som"`, `"hdbscan"` or `"constrained"`.
#' @param truth List with any of: `R_true` (p x p), `classes` (one label per row of
#'   `X`), `signal` (matrix, the true signal of each observation, for MDS) and
#'   `position` (coordinates or times of the rows, for the spurious-pattern metrics and
#'   required by `"constrained"`). One column is time: consecutive times are linked.
#'   Two columns are space: each site is linked to its 4 nearest neighbours (more when
#'   needed to connect the graph). Three columns are space-time: the first two identify
#'   the site and the third the time; a row is linked to the rows of the neighbouring
#'   sites at the same time and to its own site at the next time.
#' @param R Correlation matrix for PCA (`cor(X)` when `NULL`).
#' @param max_n Largest number of observations used by MDS, k-means, the SOM and the
#'   constrained clustering. Larger data sets (space-time grids) are scored on a random
#'   subsample of `max_n` rows, drawn with `seed`; the truth is subset alike, and the
#'   metric `subsample` records the share used. PCA always uses every row.
#' @param distance `"euclidean"` or `"gower"` (MDS).
#' @param k_max Largest number of groups for the gap statistic.
#' @param B_gap Reference data sets for the gap statistic.
#' @param seed Seed of k-means, the gap statistic and the SOM.
#' @references Legendre, P. & Legendre, L. (2012). *Numerical Ecology*, 3rd ed.
#'   Elsevier, sections 12.6 and 13.3.
#' @return A data frame with `method`, `metric` and `value`.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(150, 0), 50), matrix(rnorm(150, 3), 50))
#' run_downstream(X, "kmeans", truth = list(classes = rep(1:2, each = 50)))
#' @export
run_downstream <- function(X, method = c("pca", "mds", "kmeans", "som", "hdbscan", "constrained"), truth = list(), R = NULL,
                           distance = c("euclidean", "gower"), k_max = 6, B_gap = 20, seed = 1L, max_n = 3000) {
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
    dist_to <- function(vectors, q) if (q > 0) basis_distance(
      structure(ev$vectors[, seq_len(q), drop = FALSE], class = "depprep_basis"),
      structure(vectors[, seq_len(q), drop = FALSE], class = "depprep_basis")) else NA_real_
    dist_true <- NA_real_
    if (!is.null(truth$R_true)) {
      et <- eigen(truth$R_true, symmetric = TRUE)
      dist_true <- dist_to(et$vectors, sum(et$values > 1 + 1e-8))
    }
    # against the correlation of the true signal (residual plus class effects): the
    # classes add structure shared by all variables that R_true does not hold
    dist_signal <- NA_real_
    if (!is.null(truth$signal)) {
      es <- eigen(suppressWarnings(stats::cor(as.matrix(truth$signal))), symmetric = TRUE)
      a_s <- es$values > bs
      q_s <- if (a_s[1]) which.min(c(a_s, FALSE)) - 1 else 0
      dist_signal <- dist_to(es$vectors, max(1, q_s))
    }
    return(out(subspace_distance = dist_true, subspace_distance_signal = dist_signal, components_above_broken_stick = n_bs,
               first_above_broken_stick = as.numeric(above[1]), first_eigen_share = ev$values[1] / p))
  }
  share <- 1
  if (nrow(X) > max_n) {                                  # space-time grids: score a seeded subsample
    keep <- sort(.with_stream(seed, .streams[["subsample"]], sample(nrow(X), max_n)))
    share <- max_n / nrow(X)
    X <- X[keep, , drop = FALSE]
    for (nm in intersect(c("classes", "signal", "position"), names(truth))) {
      truth[[nm]] <- if (is.null(dim(truth[[nm]]))) truth[[nm]][keep] else as.matrix(truth[[nm]])[keep, , drop = FALSE]
    }
  }
  out0 <- out
  out <- function(...) { r <- out0(...); if (share < 1) r <- rbind(r, out0(subsample = share)); r }
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
    if (method == "hdbscan") {
      min_pts <- max(10L, as.integer(ceiling(0.02 * nrow(Z))))
      cl <- dbscan::hdbscan(Z, minPts = min_pts)$cluster        # 0 = noise
      n_cl <- length(setdiff(unique(cl), 0))
      return(out(n_clusters = n_cl, clusters_above_1 = as.numeric(n_cl > 1), noise_share = mean(cl == 0),
                 ari = if (is.na(k_true)) NA_real_ else mclust::adjustedRandIndex(cl, truth$classes),
                 ari_clustered = if (is.na(k_true) || sum(cl > 0) < 2) NA_real_ else
                   mclust::adjustedRandIndex(cl[cl > 0], truth$classes[cl > 0]),
                 min_points = min_pts))
    }
    if (method == "constrained") {
      if (is.null(truth$position)) stop("`truth$position` is required for the constrained clustering.")
      k <- if (is.na(k_true)) 3 else k_true
      links <- .contiguity_links(as.matrix(truth$position))
      tree <- adespatial::constr.hclust(stats::dist(Z), method = "ward.D2", links = links)
      cl <- stats::cutree(tree, k)
      km_cl <- stats::kmeans(Z, k, nstart = 10, iter.max = 50)$cluster
      return(out(k = k, ari_true_k = if (is.na(k_true)) NA_real_ else mclust::adjustedRandIndex(cl, truth$classes),
                 rand_kmeans = .rand_index(cl, km_cl)))
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

# Edges (two-column matrix of row indices, i < j) of the contiguity graph of the rows,
# from their positions: 1 column = times (consecutive rows in time order linked), 2 columns =
# sites (k nearest neighbours, symmetrised, k raised until the graph is connected), 3
# columns = sites (first two) observed at times (third): neighbouring sites at the same
# time and the same site at consecutive times.
.contiguity_links <- function(position, k_nn = 4) {
  edges <- function(a, b) { e <- cbind(pmin(a, b), pmax(a, b)); unique(e[e[, 1] != e[, 2], , drop = FALSE]) }
  time_links <- function(t) {
    o <- order(t)
    edges(o[-length(o)], o[-1])                             # neighbours in time order
  }
  space_links <- function(xy) {
    n <- nrow(xy)
    if (n < 2) return(matrix(integer(0), 0, 2))
    D <- as.matrix(stats::dist(xy))
    diag(D) <- Inf
    nn <- t(apply(D, 1, order))
    for (k in seq(min(k_nn, n - 1), n - 1)) {
      e <- edges(rep(seq_len(n), each = k), as.vector(t(nn[, seq_len(k), drop = FALSE])))
      if (.n_components(n, e) == 1) return(e)
    }
    e
  }
  if (ncol(position) == 1) return(time_links(position[, 1]))
  if (ncol(position) == 2) return(space_links(position))
  if (ncol(position) != 3) stop("`truth$position` must have 1 (time), 2 (space) or 3 (space and time) columns.")
  site_key <- paste(position[, 1], position[, 2], sep = "\r")
  site <- match(site_key, unique(site_key))
  time <- match(position[, 3], sort(unique(position[, 3])))
  row_of <- matrix(NA_integer_, max(site), max(time))
  row_of[cbind(site, time)] <- seq_len(nrow(position))
  xy <- position[match(seq_len(max(site)), site), 1:2, drop = FALSE]
  sl <- space_links(xy)
  in_space <- do.call(rbind, lapply(seq_len(max(time)), function(t) cbind(row_of[sl[, 1], t], row_of[sl[, 2], t])))
  in_time <- if (max(time) > 1) cbind(as.vector(row_of[, -max(time)]), as.vector(row_of[, -1])) else matrix(integer(0), 0, 2)
  e <- rbind(in_space, in_time)
  edges(e[stats::complete.cases(e), 1], e[stats::complete.cases(e), 2])
}

# Number of connected components of a graph with n nodes and an edge matrix.
.n_components <- function(n, e) {
  comp <- seq_len(n)
  repeat {
    old <- comp
    m <- pmin(comp[e[, 1]], comp[e[, 2]])
    comp[e[, 1]] <- pmin(comp[e[, 1]], m)
    comp[e[, 2]] <- pmin(comp[e[, 2]], m)
    comp <- comp[comp]
    if (identical(comp, old)) break
  }
  length(unique(comp))
}

# Rand index between two partitions (the share of pairs of observations on which they
# agree: together in both or apart in both).
.rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  together_both <- sum(choose(tab, 2))
  together_a <- sum(choose(rowSums(tab), 2))
  together_b <- sum(choose(colSums(tab), 2))
  total <- choose(n, 2)
  (total + 2 * together_both - together_a - together_b) / total
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
