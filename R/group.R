#' Group variables with similar dependence regimes
#'
#' Groups the variables by their dependence signatures, so that each group can later
#' share one basis ([build_group_basis()]). Only the signatures are used, never the
#' correlations between variables.
#'
#' Variables classified as `"none"` by [classify_regime()] form their own group, with
#' no basis. The others are grouped on the signature components of their family:
#' * spatial: log range, nugget ratio, Moran's I, log anisotropy ratio;
#' * temporal: AR(1) coefficient, seasonal strength, trend R-squared;
#' * spatio-temporal: log range, nugget ratio, Moran's I, AR(1) coefficient, seasonal
#'   strength, trend R-squared and shift indicator.
#'
#' Each component is standardised robustly (median and median absolute deviation).
#' With `method = "mclust"` (the default) a Gaussian mixture is fitted and the number of
#' groups is chosen by BIC. With `method = "hclust"` the variables are grouped by Ward's
#' method and the number of groups is the one with the largest mean silhouette (2 to
#' `max_groups`); if that silhouette is below `min_silhouette` they form a single group.
#' In the stage 3 simulation grid the mixture recovered the true groups better in
#' temporal and spatio-temporal data, while the silhouette rule almost always chose two
#' groups and merged different regimes. With `method = "type"` the first level is just
#' the regime type of [classify_regime()]; it is meant to be combined with
#' `split = "subspace"`.
#'
#' **Second level (subspace).** Variables with the same regime can still vary in
#' different places (independent large-scale fields with the same range); one basis
#' for all of them would have to cover all their patterns. With `split = "subspace"`
#' each dependent variable first gets its own basis ([build_group_basis()] on that
#' variable alone, so still marginal information), and within each regime group the
#' variables are split by the distance between their bases ([basis_distance()];
#' average linkage, cut at `split_threshold`). Subgroups are named `g1.1`, `g1.2`, ...
#'
#' **Stability.** With `B > 0` and `data`, the grouping is repeated on `B` subsamples
#' of the observations (a share `subsample` of the sites, drawn without replacement, and
#' a contiguous window of the same share of the times), with the signatures computed
#' again (without permutations) and the same number of groups. The co-association of two
#' variables is the share of subsamples in which they fall in the same group; the
#' stability of a variable is its mean co-association with the other members of its
#' group. Variables below `ambiguity` are flagged as ambiguous. Membership of the
#' `"none"` group is decided on the full data and is not resampled. Stability refers to
#' the first (regime) level.
#'
#' @param signature Output of [diagnose_variables()] (one family).
#' @param data The data the signature came from, for the stability: a `depprep_sim`
#'   object or a list with `X` and `coords`, `time` or `coords` and `index`.
#' @param method `"mclust"` (default), `"hclust"` or `"type"`.
#' @param n_groups Number of groups of dependent variables; chosen from the data when
#'   `NULL`.
#' @param max_groups Largest number of groups considered.
#' @param min_silhouette Mean silhouette below which the dependent variables form one
#'   group (`"hclust"` only).
#' @param B Number of subsamples for the stability (0 = none; 200 is a good choice for
#'   real data).
#' @param subsample Share of the observations kept in each subsample.
#' @param ambiguity Stability below which a variable is flagged as ambiguous.
#' @param thresholds Thresholds of [classify_regime()].
#' @param split `"none"` (one level) or `"subspace"` (split each regime group by the
#'   overlap of the variables' own bases; needs `data` and `split_threshold`).
#' @param split_threshold Largest basis distance (0 to 1) within a subgroup.
#' @param seed Seed of the subsamples and of the bases of the second level; required
#'   when `B > 0`.
#' @return A data frame of class `depprep_grouping` with one row per variable:
#'   `variable`, `type` (its regime type), `group`, `group_type` (most frequent type in
#'   the group), `stability` and `ambiguous`. Attributes: `method`, `n_groups`,
#'   `silhouette`, `coassociation` (matrix, dependent variables), `tree` (hclust),
#'   `period` and `ar_order` (median by group), `family`, and with
#'   `split = "subspace"` `basis_distance` (one matrix per regime group).
#' @examples
#' g <- sim_groups(sim_group("long", 4, "exponential", range = 0.5),
#'                 sim_group("short", 4, "exponential", range = 0.1),
#'                 sim_group("none", 3, "none"))
#' s <- sim_spatial(400, g, seed = 1)
#' group_regimes(diagnose_variables(s))
#' @importFrom mclust mclustBIC
#' @export
group_regimes <- function(signature, data = NULL, method = c("mclust", "hclust", "type"), n_groups = NULL,
                          max_groups = 8, min_silhouette = 0.25, B = 0, subsample = 0.8, ambiguity = 0.7,
                          thresholds = default_thresholds(), split = c("none", "subspace"), split_threshold = NULL,
                          seed = NULL) {
  method <- match.arg(method)
  split <- match.arg(split)
  if (split == "subspace" && (is.null(data) || is.null(split_threshold))) {
    stop("`data` and `split_threshold` are required when `split = \"subspace\"`.", call. = FALSE)
  }
  stopifnot(inherits(signature, "data.frame"), "family" %in% names(signature))
  family <- unique(signature$family)
  if (length(family) != 1) stop("`signature` must hold a single family.", call. = FALSE)
  if (B > 0 && (is.null(data) || is.null(seed))) stop("`data` and `seed` are required when `B > 0`.", call. = FALSE)

  type <- unname(classify_regime(signature, thresholds))
  dep <- type != "none"
  group <- rep("none", nrow(signature))
  fit <- list(cluster = integer(0), k = 0L, silhouette = NA_real_, tree = NULL)
  if (any(dep)) {
    Z <- .grouping_components(signature[dep, , drop = FALSE])
    fit <- if (method == "type") {
      list(cluster = as.integer(factor(type[dep])), k = length(unique(type[dep])), silhouette = NA_real_, tree = NULL)
    } else .cluster(Z, method, n_groups, max_groups, min_silhouette)
    # groups named g1, g2, ... by decreasing size
    lab <- names(sort(table(fit$cluster), decreasing = TRUE))
    group[dep] <- paste0("g", match(as.character(fit$cluster), lab))
  }
  group_type <- vapply(group, function(g) {
    t <- table(type[group == g])
    names(t)[which.max(t)]
  }, character(1), USE.NAMES = FALSE)

  stability <- rep(NA_real_, nrow(signature))
  co <- NULL
  if (B > 0 && sum(dep) > 1) {
    co <- .with_stream(seed, .streams[["subsample"]],
                       .coassociation(signature, data, dep, method, fit$k, B, subsample))
    gd <- group[dep]
    stability[dep] <- vapply(seq_len(sum(dep)), function(i) {
      same <- setdiff(which(gd == gd[i]), i)
      if (length(same)) mean(co[i, same]) else NA_real_
    }, numeric(1))
    dimnames(co) <- list(signature$variable[dep], signature$variable[dep])
  }
  distances <- NULL
  if (split == "subspace" && any(dep)) {
    sp <- .split_by_subspace(signature, data, group, group_type, split_threshold, if (is.null(seed)) 1L else seed)
    group <- sp$group
    distances <- sp$distances
  }
  by_group <- function(col) {
    v <- tapply(signature[[col]], group, stats::median, na.rm = TRUE)
    v[setdiff(names(v), "none")]
  }
  out <- data.frame(variable = signature$variable, type = type, group = group, group_type = group_type,
                    stability = stability, ambiguous = !is.na(stability) & stability < ambiguity,
                    stringsAsFactors = FALSE)
  structure(out, class = c("depprep_grouping", "data.frame"), method = method, n_groups = fit$k,
            silhouette = fit$silhouette, coassociation = co, tree = fit$tree,
            period = round(by_group("period")), ar_order = round(by_group("ar_order")), family = family,
            basis_distance = distances)
}

# The data as a list with X and the coordinates, times or index of its family.
.data_parts <- function(data) {
  if (inherits(data, "depprep_sim")) {
    data <- list(X = data$X, coords = data$coords, time = if (data$family == "temporal") data$time, index = data$index)
  }
  data$X <- as.matrix(data$X)
  data
}

# Second grouping level: within each regime group, split the variables by the distance
# between their own bases (each selected on the variable alone).
.split_by_subspace <- function(signature, data, group, group_type, threshold, seed) {
  d <- .data_parts(data)
  args <- list(coords = d$coords, time = d$time, index = d$index)
  args <- args[!vapply(args, is.null, TRUE)]
  if (signature$family[1] == "spatiotemporal") args$time <- NULL
  distances <- list()
  var_seed <- seed + rank(signature$variable, ties.method = "first")   # by name, not by column position
  for (g in setdiff(unique(group), "none")) {
    idx <- which(group == g)
    if (length(idx) < 2) next
    per <- stats::median(signature$period[idx], na.rm = TRUE)
    bases <- lapply(idx, function(i) do.call(build_group_basis, c(list(
      X_group = d$X[, i, drop = FALSE], type = group_type[i], period = if (is.finite(per)) round(per),
      seed = var_seed[i]), args)))
    D <- outer(seq_along(idx), seq_along(idx), Vectorize(function(a, b) basis_distance(bases[[a]], bases[[b]])))
    diag(D) <- 0
    dimnames(D) <- list(signature$variable[idx], signature$variable[idx])
    distances[[g]] <- D
    sub <- stats::cutree(stats::hclust(stats::as.dist(D), method = "average"), h = threshold)
    if (length(unique(sub)) > 1) {
      lab <- names(sort(table(sub), decreasing = TRUE))     # subgroups by decreasing size
      group[idx] <- paste0(g, ".", match(as.character(sub), lab))
    }
  }
  list(group = group, distances = distances)
}

# Signature components used for grouping (decision D5), robustly standardised.
.grouping_components <- function(sig) {
  fam <- sig$family[1]
  lg <- function(x) log(pmax(x, 1e-6))
  Z <- switch(fam,
    spatial = cbind(log_range = lg(sig$range), nugget_ratio = sig$nugget_ratio, moran_I = sig$moran_I,
                    log_aniso = lg(sig$aniso_ratio)),
    temporal = cbind(phi = sig$phi, seasonal_strength = sig$period_strength, trend_r2 = sig$trend_r2),
    spatiotemporal = cbind(log_range = lg(sig$range), nugget_ratio = sig$nugget_ratio, moran_I = sig$moran_I,
                           phi = sig$phi, seasonal_strength = sig$period_strength, trend_r2 = sig$trend_r2,
                           shift = sig$shift_indicator))
  apply(Z, 2, function(x) {
    x[!is.finite(x)] <- stats::median(x[is.finite(x)])    # a failed fit takes the typical value
    if (all(is.na(x))) return(rep(0, length(x)))
    s <- stats::mad(x)
    if (!is.finite(s) || s < 1e-8) s <- stats::sd(x)
    if (!is.finite(s) || s < 1e-8) return(rep(0, length(x)))
    (x - stats::median(x)) / s
  })
}

.cluster <- function(Z, method, n_groups = NULL, max_groups = 8, min_silhouette = 0.25) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  one <- list(cluster = rep(1L, n), k = 1L, silhouette = NA_real_, tree = NULL)
  if (n < 3) return(one)
  d <- stats::dist(Z)
  if (method == "hclust") {
    tree <- stats::hclust(d, method = "ward.D2")
    one$tree <- tree
    if (!is.null(n_groups)) {
      k <- min(n_groups, n)
      cl <- stats::cutree(tree, k)
      return(list(cluster = cl, k = as.integer(k), silhouette = .mean_silhouette(cl, d), tree = tree))
    }
    ks <- 2:min(max_groups, n - 1)
    sil <- vapply(ks, function(k) .mean_silhouette(stats::cutree(tree, k), d), 0)
    if (!length(ks) || max(sil) < min_silhouette) { one$silhouette <- if (length(sil)) max(sil) else NA; return(one) }
    k <- ks[which.max(sil)]
    return(list(cluster = stats::cutree(tree, k), k = as.integer(k), silhouette = max(sil), tree = tree))
  }
  G <- if (is.null(n_groups)) seq_len(min(max_groups, n - 1)) else n_groups
  m <- tryCatch(mclust::Mclust(Z, G = G, verbose = FALSE), error = function(e) NULL)
  if (is.null(m)) return(one)
  cl <- m$classification
  list(cluster = cl, k = as.integer(length(unique(cl))), silhouette = .mean_silhouette(cl, d), tree = NULL)
}

.mean_silhouette <- function(cl, d) {
  if (length(unique(cl)) < 2) return(NA_real_)
  mean(cluster::silhouette(cl, d)[, "sil_width"])
}

# Co-association of the dependent variables over B subsamples of the observations.
# Uses the random stream of the caller.
.coassociation <- function(signature, data, dep, method, k, B, subsample) {
  data <- .data_parts(data)
  X <- data$X[, dep, drop = FALSE]
  p <- ncol(X)
  co <- matrix(0, p, p)
  if (k <= 1) return(co + 1)
  family <- signature$family[1]
  window <- function(times) {
    u <- sort(unique(times))
    L <- max(3, round(subsample * length(u)))
    s <- sample.int(length(u) - L + 1, 1)
    u[s:(s + L - 1)]
  }
  for (b in seq_len(B)) {
    sig <- if (family == "spatial") {
      keep <- sort(sample.int(nrow(X), round(subsample * nrow(X))))
      diagnose_variables(X[keep, , drop = FALSE], coords = as.matrix(data$coords)[keep, , drop = FALSE], nsim = 0)
    } else if (family == "temporal") {
      w <- data$time %in% window(data$time)
      diagnose_variables(X[w, , drop = FALSE], time = data$time[w], nsim = 0)
    } else {
      coords <- as.matrix(data$coords)
      sites <- sort(sample.int(nrow(coords), round(subsample * nrow(coords))))
      tw <- window(data$index$time)
      w <- data$index$site %in% sites & data$index$time %in% tw
      idx <- data.frame(site = match(data$index$site[w], sites), time = match(data$index$time[w], tw))
      diagnose_variables(X[w, , drop = FALSE], coords = coords[sites, , drop = FALSE], index = idx, nsim = 0)
    }
    cl <- if (method == "type") as.integer(factor(classify_regime(sig))) else
      .cluster(.grouping_components(sig), method, n_groups = k)$cluster
    co <- co + outer(cl, cl, "==")
  }
  co / B
}
