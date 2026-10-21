#' Build the shared basis of a group of variables
#'
#' Selects one spatial, temporal or spatio-temporal basis for a whole group, treating
#' the variables of the group as a multivariate response. Every variable of the group
#' is later corrected with this same basis (the core rule of the project).
#'
#' Candidates depend on the type of the group:
#' * spatial (`"isotropic"`, `"anisotropic"`): distance-based Moran eigenvector maps
#'   with positive autocorrelation (`adespatial::dbmem()`);
#' * temporal: `"seasonal"`, harmonics of `period`; `"trend"`, a cubic polynomial
#'   followed by temporal dbMEM; `"ar"`, temporal dbMEM (slow persistence);
#' * spatio-temporal: products of a spatial and a temporal basis, each first selected
#'   on the stacked data (all times as responses for the spatial part, all sites for
#'   the temporal part) and each including the constant, so that purely spatial and
#'   purely temporal terms are candidates too. `"temporal_common"` gets temporal terms
#'   only and `"spatial_only"` spatial terms only;
#' * `"none"`: an empty basis.
#'
#' Candidates are made orthonormal, so the share of variance a candidate explains does
#' not depend on the ones already selected. Forward selection adds candidates in order
#' of explained variance (summed over the variables of the group, which uses no
#' cross-products between variables) and stops at the first that is not significant
#' in a permutation test, when the adjusted R-squared of the model reaches that of all
#' candidates (the double stopping of Blanchet et al. 2008, here checked after adding
#' the term, so that one dominant pattern can enter), or when `k / n` would exceed
#' `kn_max`. A global permutation test of all candidates comes first; if it is not
#' significant the basis is empty.
#'
#' @param X_group Numeric matrix, one column per variable of the group.
#' @param coords Site coordinates (spatial and spatio-temporal data).
#' @param time Times of the rows (temporal data).
#' @param index Data frame with the `site` (row of `coords`) and `time` of each row
#'   (spatio-temporal data, complete grid).
#' @param type Regime type of the group, as returned by [classify_regime()].
#' @param period Seasonal period, required for the seasonal types.
#' @param kn_max Largest share of degrees of freedom the basis may use (`k / n`).
#' @param alpha Significance level of the permutation tests.
#' @param nperm Number of permutations per test.
#' @param seed Seed of the permutations.
#' @return A matrix with orthonormal columns (one row per observation) of class
#'   `depprep_basis`, with attributes `type`, `k`, `kn` (k / n), `r2` (share of the
#'   group's variance explained) and `terms` (names of the selected candidates).
#' @references Blanchet, F. G., Legendre, P. and Borcard, D. (2008). Forward selection
#'   of explanatory variables. Ecology 89, 2623-2632.
#' @examples
#' g <- sim_groups(sim_group("long", 4, "exponential", range = 0.5))
#' s <- sim_spatial(200, g, seed = 1)
#' b <- build_group_basis(s$X, coords = s$coords, type = "isotropic", seed = 1)
#' attr(b, "k")
#' @export
build_group_basis <- function(X_group, coords = NULL, time = NULL, index = NULL, type, period = NULL,
                              kn_max = 0.25, alpha = 0.05, nperm = 199, seed = 1L) {
  Y <- scale(as.matrix(X_group))
  Y[!is.finite(Y)] <- 0                                   # constant columns carry no information
  n <- nrow(Y)
  spatial_types <- c("isotropic", "anisotropic")
  temporal_types <- c("seasonal", "trend", "ar")
  st_types <- c("temporal_common", "spatial_only", "moving", "separable_seasonal", "separable_ar1")
  if (!type %in% c("none", spatial_types, temporal_types, st_types)) stop("Unknown `type`: ", type, call. = FALSE)
  if (type %in% c("seasonal", "separable_seasonal") && is.null(period)) {
    stop("`period` is required for seasonal groups.", call. = FALSE)
  }
  if (type == "none") return(.as_basis(matrix(0, n, 0), type, n, 0, character(0)))

  .with_stream(seed, .streams[["basis"]], {
    if (type %in% spatial_types) {
      if (is.null(coords)) stop("`coords` is required for spatial groups.", call. = FALSE)
      sel <- .forward_select(Y, .spatial_candidates(coords), kn_max, alpha, nperm)
    } else if (type %in% temporal_types) {
      if (is.null(time)) stop("`time` is required for temporal groups.", call. = FALSE)
      ord <- order(time)
      cand <- .temporal_candidates(time[ord], type, period)
      sel <- .forward_select(Y[ord, , drop = FALSE], cand, kn_max, alpha, nperm)
      sel$basis[ord, ] <- sel$basis                        # back to the order of the rows
    } else {
      if (is.null(coords) || is.null(index)) stop("`coords` and `index` are required for spatio-temporal groups.", call. = FALSE)
      sel <- .st_select(Y, coords, index, type, period, kn_max, alpha, nperm)
    }
    .as_basis(sel$basis, type, n, sel$r2, sel$terms)
  })
}

#' Union of bases
#'
#' The basis spanned by all the given bases together, used to estimate the association
#' between variables of different groups. Redundant directions are dropped.
#'
#' @param ... Bases from [build_group_basis()] (same observations, same row order).
#' @param kn_max Largest share of degrees of freedom; a warning is given above it. The
#'   default, 0.3, was not exceeded by any union of two groups in the simulation grid.
#' @return A matrix with orthonormal columns of class `depprep_basis` (type `"union"`).
#' @examples
#' g <- sim_groups(sim_group("long", 3, "exponential", range = 0.5),
#'                 sim_group("short", 3, "exponential", range = 0.1))
#' s <- sim_spatial(200, g, seed = 1)
#' b1 <- build_group_basis(s$X[, 1:3], coords = s$coords, type = "isotropic", seed = 1)
#' b2 <- build_group_basis(s$X[, 4:6], coords = s$coords, type = "isotropic", seed = 1)
#' attr(union_basis(b1, b2), "k")
#' @export
union_basis <- function(..., kn_max = 0.3) {
  bases <- list(...)
  n <- unique(vapply(bases, nrow, integer(1)))
  if (length(n) != 1) stop("All bases must have the same number of rows.", call. = FALSE)
  M <- do.call(cbind, lapply(bases, unclass))
  terms <- unlist(lapply(bases, function(b) attr(b, "terms")))
  if (!ncol(M)) return(.as_basis(matrix(0, n, 0), "union", n, NA_real_, character(0)))
  q <- qr(M, tol = 1e-7)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  if (q$rank / n > kn_max) {
    warning(sprintf("The union basis uses k/n = %.2f, above kn_max = %.2f.", q$rank / n, kn_max), call. = FALSE)
  }
  .as_basis(Q, "union", n, NA_real_, terms[q$pivot[seq_len(q$rank)]])
}

#' @export
print.depprep_basis <- function(x, ...) {
  cat(sprintf("<depprep_basis> type %s: k = %d of n = %d (k/n = %.3f)\n",
              attr(x, "type"), attr(x, "k"), nrow(x), attr(x, "kn")))
  invisible(x)
}

.as_basis <- function(B, type, n, r2, terms) {
  B <- as.matrix(B)
  dimnames(B) <- list(NULL, terms)
  structure(B, class = "depprep_basis", type = type, k = ncol(B), kn = ncol(B) / n, r2 = r2, terms = terms)
}

# Orthonormal columns, centred, without (near) null directions; keeps the column order
# so that earlier candidates (for example the trend) take precedence over later ones.
.orthonormal <- function(M) {
  M <- scale(as.matrix(M), scale = FALSE)
  keep <- sqrt(colSums(M^2)) > 1e-8 * sqrt(nrow(M))
  M <- M[, keep, drop = FALSE]
  if (!ncol(M)) return(M)
  q <- qr(M, tol = 1e-7)
  idx <- sort(q$pivot[seq_len(q$rank)])                  # independent columns, original order
  q <- qr(M[, idx, drop = FALSE])
  Q <- qr.Q(q)
  Q <- sweep(Q, 2, sign(colSums(Q * M[, idx, drop = FALSE])), "*")   # same sign as the candidate
  colnames(Q) <- colnames(M)[idx]
  Q
}

.spatial_candidates <- function(coords) {
  m <- as.matrix(adespatial::dbmem(as.matrix(coords), MEM.autocor = "positive", silent = TRUE))
  .orthonormal(m)
}

.temporal_candidates <- function(time, type, period = NULL) {
  t <- as.numeric(time)
  if (type %in% c("seasonal", "separable_seasonal")) {
    h <- seq_len(floor(period / 2))
    m <- do.call(cbind, lapply(h, function(i) cbind(sin(2 * pi * i * t / period), cos(2 * pi * i * t / period))))
    colnames(m) <- paste0(rep(c("sin", "cos"), length(h)), rep(h, each = 2), "_p", period)
    return(.orthonormal(m))
  }
  mem <- as.matrix(adespatial::dbmem(matrix(t), MEM.autocor = "positive", silent = TRUE))
  colnames(mem) <- paste0("tMEM", seq_len(ncol(mem)))
  if (type == "trend") {
    p <- stats::poly(t, 3)
    colnames(p) <- paste0("trend", 1:3)
    mem <- cbind(p, mem)
  }
  .orthonormal(mem)
}

# Forward selection over orthonormal candidates (see build_group_basis()). Uses the
# random stream of the caller.
.forward_select <- function(Y, cand, kn_max, alpha, nperm) {
  n <- nrow(Y)
  empty <- list(basis = matrix(0, n, 0), r2 = 0, terms = character(0))
  if (!ncol(cand)) return(empty)
  tss <- sum(Y^2)
  contrib <- rowSums(crossprod(cand, Y)^2)               # variance explained by each candidate
  m <- ncol(cand)
  # global test: all candidates together
  r2_all <- sum(contrib) / tss
  perm_all <- vapply(seq_len(nperm), function(b) sum(crossprod(cand, Y[sample.int(n), , drop = FALSE])^2), 0)
  if ((1 + sum(perm_all >= sum(contrib))) / (nperm + 1) > alpha) return(empty)
  adj_all <- if (m < n - 1) 1 - (1 - r2_all) * (n - 1) / (n - m - 1) else Inf
  ord <- order(contrib, decreasing = TRUE)
  chosen <- integer(0)
  explained <- 0
  for (j in ord) {
    k <- length(chosen) + 1
    if (k / n > kn_max) break
    res <- Y - cand[, chosen, drop = FALSE] %*% crossprod(cand[, chosen, drop = FALSE], Y)
    perm <- vapply(seq_len(nperm), function(b) sum(crossprod(cand[, j], res[sample.int(n), , drop = FALSE])^2), 0)
    if ((1 + sum(perm >= contrib[j])) / (nperm + 1) > alpha) break
    chosen <- c(chosen, j)
    explained <- explained + contrib[j]
    # second stopping rule, checked after adding the term: otherwise a single dominant
    # pattern (one harmonic explaining most of a seasonal group) could never enter
    if (1 - (1 - explained / tss) * (n - 1) / (n - k - 1) >= adj_all) break
  }
  list(basis = cand[, chosen, drop = FALSE], r2 = unname(explained / tss), terms = colnames(cand)[chosen])
}

# Spatio-temporal basis: marginal selections on the stacked data, then forward
# selection among their products.
.st_select <- function(Y, coords, index, type, period, kn_max, alpha, nperm) {
  sites <- index$site
  times <- match(index$time, sort(unique(index$time)))
  n_s <- nrow(as.matrix(coords))
  n_t <- max(times)
  if (nrow(Y) != n_s * n_t || anyDuplicated(cbind(sites, times))) {
    stop("Spatio-temporal bases need a complete grid (every site at every time).", call. = FALSE)
  }
  # one matrix per variable, sites x times; stacked by columns for each margin
  cube <- lapply(seq_len(ncol(Y)), function(j) { m <- matrix(0, n_s, n_t); m[cbind(sites, times)] <- Y[, j]; m })
  by_site <- do.call(cbind, cube)                          # sites x (times * variables)
  by_time <- do.call(cbind, lapply(cube, t))               # times x (sites * variables)
  s_terms <- if (type == "temporal_common") NULL else {
    .forward_select(scale(by_site, scale = FALSE), .spatial_candidates(coords), kn_max, alpha, nperm)$basis
  }
  t_terms <- if (type == "spatial_only") NULL else {
    t_type <- if (type == "separable_seasonal") "seasonal" else "ar"
    .forward_select(scale(by_time, scale = FALSE), .temporal_candidates(sort(unique(index$time)), t_type, period),
                    kn_max, alpha, nperm)$basis
  }
  S <- cbind(const = rep(1 / sqrt(n_s), n_s), s_terms)
  Tm <- cbind(const = rep(1 / sqrt(n_t), n_t), t_terms)
  pairs <- expand.grid(s = seq_len(ncol(S)), t = seq_len(ncol(Tm)))[-1, , drop = FALSE]   # drop constant x constant
  if (!nrow(pairs)) return(list(basis = matrix(0, nrow(Y), 0), r2 = 0, terms = character(0)))
  cand <- vapply(seq_len(nrow(pairs)), function(i) S[sites, pairs$s[i]] * Tm[times, pairs$t[i]], numeric(nrow(Y)))
  cand <- matrix(cand, nrow(Y))
  colnames(cand) <- paste(colnames(S)[pairs$s], colnames(Tm)[pairs$t], sep = ":")
  .forward_select(Y, cand, kn_max, alpha, nperm)
}

# Share of the true large-scale structure removed by a basis (capture) and share of the
# true residual removed beyond what k/n alone implies (absorption). Simulation only.
.capture_absorption <- function(basis, large, residual) {
  B <- unclass(basis)
  proj <- function(M) { M <- scale(as.matrix(M), scale = FALSE); if (!ncol(B)) return(0); sum(crossprod(B, M)^2) / sum(M^2) }
  c(capture = proj(large), absorption = proj(residual) - ncol(B) / nrow(B))
}
