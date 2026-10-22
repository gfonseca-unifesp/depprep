#' Correct a data set for the dependence of each group
#'
#' Removes the large-scale structure of each group of variables with the group's own
#' basis, so that what remains (the residual components) can be compared across
#' variables. Two outputs, as fixed by the project's decision P1:
#' * `R`: a correlation matrix assembled pair by pair — within a group with the group's
#'   basis, across groups with the union of the two bases — and projected to the
#'   nearest valid correlation matrix (`adjustment` records how much it moved); for
#'   PCA and other methods that take a correlation matrix;
#' * `X`: a data matrix of residuals on the union of all group bases; for k-means,
#'   SOM and MDS on observations.
#'
#' Modes:
#' * `"conditional"` (residuals on the basis): `(I - P) x`;
#' * `"whiten"`: the observations are whitened with the covariance estimated for each
#'   group — a variogram model pooled over the group's variables in space, an AR model
#'   in time (after removing the trend or seasonal terms of the group's basis), and both
#'   (separably) in space-time. Uses the whole covariance, not only a basis, which can
#'   matter for rough processes. Pairs across groups use each group's own whitening;
#' * `"scale"`: one correlation matrix per scale band (broad, medium, fine thirds of the
#'   Moran eigenvector maps, or of the temporal ones), from the projection of every
#'   variable on the band; spatial and temporal data only.
#'
#' Scope of the bases: `"groups"` (one basis per group, unions across groups) or
#' `"global"` (one basis for all dependent variables: the union of the bases selected on
#' all of them for each regime type present). Variables without dependence are never
#' corrected among themselves.
#'
#' Defaults (decision D7, calibrated on the stage 5 simulations, where all basis
#' approaches gave nearly the same error): spatial and temporal data use the conditional
#' mode with the global basis, the simplest; spatio-temporal data use whitening with the
#' groups, which erred about half as much there, because the shared process persists in
#' time and the residual does not. In space-time the result is then the correlation of
#' the innovations, which equals the residual correlation when the residual association
#' does not itself persist in time.
#'
#' @param X Numeric matrix (observations x variables), or a `depprep_sim` object.
#' @param grouping Output of [group_regimes()], or a named character vector giving the
#'   group of each variable (`"none"` for variables without dependence).
#' @param coords,time,index Coordinates, times or spatio-temporal index, as in
#'   [diagnose_variables()].
#' @param bases Optional named list with one basis per group (from
#'   [build_group_basis()]); built when `NULL`.
#' @param mode `"conditional"`, `"whiten"` or `"scale"`; `NULL` for the default of the
#'   data family (see above).
#' @param scope `"global"` or `"groups"`; `NULL` for the default (`"global"` for the
#'   conditional and scale modes in spatial and temporal data, `"groups"` otherwise and
#'   whenever `bases` is given).
#' @param kn_max Share of degrees of freedom above which a union of bases is flagged.
#' @param seed Seed of the basis selection (when `bases` is `NULL`).
#' @return An object of class `depprep_corrected`: a list with `R`, `R_pairwise`
#'   (before the projection), `adjustment`, `X` (residuals or whitened data; `NULL` in
#'   scale mode), `R_scale` (scale mode), `bases`, `kn` (by group), `union_kn` (by
#'   pair of groups), `warnings`, `mode`, `scope` and `group` (`"global"` for every
#'   dependent variable under the global scope).
#' @examples
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(300, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 1)
#' gr <- group_regimes(diagnose_variables(s, nsim = 99))
#' cr <- correct(s, gr)
#' round(cr$R, 2)
#' @export
correct <- function(X, grouping, coords = NULL, time = NULL, index = NULL, bases = NULL,
                    mode = NULL, scope = NULL, kn_max = 0.3, seed = 1L) {
  if (inherits(X, "depprep_sim")) {
    d <- .data_parts(X); X <- d$X; coords <- d$coords; time <- d$time; index <- d$index
  }
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else if (!is.null(coords)) "spatial" else
    if (!is.null(time)) "temporal" else stop("Give `coords`, `time` or `coords` and `index`.", call. = FALSE)
  if (is.null(mode)) mode <- if (family == "spatiotemporal") "whiten" else "conditional"
  mode <- match.arg(mode, c("conditional", "whiten", "scale"))
  if (is.null(scope)) scope <- if (is.null(bases) && mode != "whiten" && family != "spatiotemporal") "global" else "groups"
  scope <- match.arg(scope, c("groups", "global"))
  grp <- .group_vector(grouping, colnames(X))
  gtype <- .group_types(grouping, grp)
  args <- switch(family, spatial = list(coords = coords), temporal = list(time = time),
                 spatiotemporal = list(coords = coords, index = index))
  if (scope == "global" && any(grp != "none")) {
    dep <- grp != "none"
    vtype <- if (inherits(grouping, "data.frame") && "type" %in% names(grouping)) {
      grouping$type[match(colnames(X), grouping$variable)]
    } else unlist(gtype[grp])
    types <- unique(vtype[dep])
    if (is.null(bases)) {
      period <- if (inherits(grouping, "depprep_grouping")) attr(grouping, "period") else NULL
      per <- if (!is.null(period) && any(is.finite(period))) round(stats::median(period[is.finite(period)])) else 12
      bs <- lapply(types, function(t) do.call(build_group_basis, c(list(X_group = X[, dep, drop = FALSE], type = t,
                                                                        period = per, seed = seed), args)))
      bases <- list(global = suppressWarnings(do.call(union_basis, c(bs, kn_max = kn_max))))
    }
    grp[dep] <- "global"
    gtype <- list(global = if (any(types %in% c("seasonal", "trend"))) "seasonal" else types[1], none = "none")
  }
  groups <- setdiff(unique(grp), "none")
  if (is.null(bases)) {
    period <- if (inherits(grouping, "depprep_grouping")) attr(grouping, "period") else NULL
    bases <- stats::setNames(lapply(groups, function(g) {
      per <- if (!is.null(period) && g %in% names(period) && is.finite(period[[g]])) period[[g]] else 12
      do.call(build_group_basis, c(list(X_group = X[, grp == g, drop = FALSE], type = gtype[[g]], period = per,
                                        seed = seed), args))
    }), groups)
  }
  n <- nrow(X)
  empty <- .as_basis(matrix(0, n, 0), "none", n, 0, character(0))
  basis_of <- function(g) if (g == "none") empty else bases[[g]]
  kn <- vapply(groups, function(g) attr(bases[[g]], "kn"), 0)

  if (mode == "scale") {
    out <- .correct_scale(X, grp, family, coords, time, groups, bases, kn)
    out$scope <- scope
    return(out)
  }

  # transformed data per group: residuals on a basis, or whitened observations
  transform <- if (mode == "conditional") {
    function(cols, B) .resid_on(X[, cols, drop = FALSE], B)
  } else {
    W <- lapply(stats::setNames(c(groups, "none"), c(groups, "none")), function(g)
      .whitener(X[, grp == g, drop = FALSE], g, gtype[[g]], bases[[g]], family, coords, time, index))
    function(cols, B, g) W[[g]](X[, cols, drop = FALSE])
  }
  ug <- c(groups, if (any(grp == "none")) "none")
  R <- diag(ncol(X)); dimnames(R) <- list(colnames(X), colnames(X))
  union_kn <- NULL
  warn <- character(0)
  for (a in seq_along(ug)) for (b in a:length(ug)) {
    ga <- ug[a]; gb <- ug[b]
    ia <- which(grp == ga); ib <- which(grp == gb)
    if (mode == "conditional") {
      U <- if (ga == gb) basis_of(ga) else suppressWarnings(union_basis(basis_of(ga), basis_of(gb), kn_max = kn_max))
      Ea <- transform(ia, U); Eb <- transform(ib, U)
      if (ga != gb) {
        union_kn <- rbind(union_kn, data.frame(group_a = ga, group_b = gb, k = attr(U, "k"), kn = attr(U, "kn")))
        if (attr(U, "kn") > kn_max) warn <- c(warn, sprintf("union of %s and %s: k/n = %.2f > %.2f", ga, gb, attr(U, "kn"), kn_max))
      }
    } else {
      Ea <- transform(ia, NULL, ga); Eb <- transform(ib, NULL, gb)
      m <- min(nrow(Ea), nrow(Eb))                       # AR filters drop the first observations
      Ea <- utils::tail(Ea, m); Eb <- utils::tail(Eb, m)
    }
    R[ia, ib] <- stats::cor(Ea, Eb)
    R[ib, ia] <- t(R[ia, ib])
  }
  diag(R) <- 1
  near <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  dimnames(near) <- dimnames(R)
  Xout <- if (mode == "conditional") {
    U_all <- if (length(groups)) suppressWarnings(do.call(union_basis, c(unname(bases), kn_max = kn_max))) else empty
    E <- .resid_on(X, U_all); colnames(E) <- colnames(X); E
  } else {
    parts <- lapply(ug, function(g) { cols <- which(grp == g); E <- transform(cols, NULL, g); colnames(E) <- colnames(X)[cols]; E })
    m <- min(vapply(parts, nrow, 0L))
    do.call(cbind, lapply(parts, utils::tail, m))[, colnames(X), drop = FALSE]
  }
  structure(list(R = near, R_pairwise = R, adjustment = sqrt(sum((near - R)^2)), X = Xout, R_scale = NULL,
                 bases = bases, kn = kn, union_kn = union_kn, warnings = warn, mode = mode, scope = scope, group = grp),
            class = "depprep_corrected")
}

#' @export
print.depprep_corrected <- function(x, ...) {
  cat(sprintf("<depprep_corrected> mode %s, %s bases: %d variables in %d groups (+ %d without dependence)\n", x$mode, x$scope,
              length(x$group), length(setdiff(unique(x$group), "none")), sum(x$group == "none")))
  if (length(x$kn)) cat("k/n by group:", paste(sprintf("%s %.3f", names(x$kn), x$kn), collapse = "; "), "\n")
  if (length(x$warnings)) cat("Warnings:", paste(x$warnings, collapse = "; "), "\n")
  invisible(x)
}

# Group of each variable, from a grouping table or a named vector.
.group_vector <- function(grouping, vars) {
  g <- if (inherits(grouping, "data.frame")) stats::setNames(grouping$group, grouping$variable) else grouping
  if (is.null(names(g)) || !all(vars %in% names(g))) stop("`grouping` must name the group of every column of `X`.", call. = FALSE)
  unname(g[vars])
}

.group_types <- function(grouping, grp) {
  if (inherits(grouping, "data.frame")) {
    t <- tapply(grouping$group_type, grouping$group, `[`, 1)
    return(as.list(t))
  }
  stats::setNames(as.list(rep("isotropic", length(unique(grp)))), unique(grp))
}

# Centred residuals of the columns of M on the basis B.
.resid_on <- function(M, B) {
  M <- scale(as.matrix(M), scale = FALSE)
  B <- unclass(B)
  if (!ncol(B)) return(M)
  M - B %*% crossprod(B, M)
}

# A function that whitens the columns of a group. Spatial: pooled variogram model of
# the group (nugget + structured part) as the covariance. Temporal: AR filter fitted on
# the pooled series after removing the deterministic (trend, seasonal) terms of the
# group's basis. Space-time: spatial whitening of each time and AR filtering of each
# site. No dependence: centring only.
.whitener <- function(Xg, g, type, basis, family, coords, time, index) {
  if (g == "none" || !ncol(Xg)) {
    ord0 <- if (family == "temporal") order(time) else NULL
    return(function(M) { M <- as.matrix(M); if (!is.null(ord0)) M <- M[ord0, , drop = FALSE]; scale(M, scale = FALSE) })
  }
  Z <- scale(Xg)
  spatial_W <- function(cc, zz) {
    zz <- as.matrix(zz)
    fits <- lapply(seq_len(ncol(zz)), function(j) .fit_variogram(cc, zz[, j]))
    ok <- vapply(fits, function(f) !is.null(f$fit), TRUE)
    if (!any(ok)) return(diag(nrow(zz)))
    fit <- fits[[which(ok)[which.min(abs(vapply(fits[ok], `[[`, 0, "range") -
                                           stats::median(vapply(fits[ok], `[[`, 0, "range"))))]]]$fit
    d <- as.matrix(stats::dist(as.matrix(cc)))
    S <- gstat::variogramLine(fit, dist_vector = d, covariance = TRUE)
    solve(t(chol(S + diag(1e-8 * mean(diag(S)), nrow(S)))))
  }
  ar_coef <- function(series) {
    orders <- vapply(seq_len(ncol(series)), function(j) stats::ar(series[, j], aic = TRUE, order.max = 5)$order, 0)
    p <- max(1, round(stats::median(orders)))
    coefs <- vapply(seq_len(ncol(series)), function(j) stats::ar(series[, j], aic = FALSE, order.max = p)$ar, numeric(p))
    matrix(coefs, nrow = p) |> rowMeans()
  }
  ar_filter <- function(M, phi) {
    p <- length(phi)
    M <- as.matrix(M)
    out <- M[-(1:p), , drop = FALSE]
    for (i in seq_len(p)) out <- out - phi[i] * M[(p + 1 - i):(nrow(M) - i), , drop = FALSE]
    out
  }
  if (family == "spatial") {
    W <- spatial_W(coords, Z)
    return(function(M) W %*% scale(as.matrix(M), scale = FALSE))
  }
  if (family == "temporal") {
    ord <- order(time)
    det <- if (type %in% c("seasonal", "trend")) basis else NULL
    pre <- function(M) { M <- as.matrix(M)[ord, , drop = FALSE]; if (is.null(det)) scale(M, scale = FALSE) else .resid_on(M, unclass(det)[ord, , drop = FALSE]) }
    phi <- ar_coef(pre(Z))
    return(function(M) ar_filter(pre(M), phi))
  }
  # space-time: rows ordered by index; whiten sites at each time, then filter each site
  n_s <- nrow(as.matrix(coords)); tt <- match(index$time, sort(unique(index$time))); n_t <- max(tt)
  to_cube <- function(v) { m <- matrix(0, n_s, n_t); m[cbind(index$site, tt)] <- v; m }
  site_means <- do.call(cbind, lapply(seq_len(ncol(Z)), function(j) rowMeans(to_cube(Z[, j]))))
  W <- spatial_W(coords, site_means)
  time_means <- do.call(cbind, lapply(seq_len(ncol(Z)), function(j) colMeans(to_cube(Z[, j]))))
  phi <- ar_coef(scale(time_means, scale = FALSE))
  function(M) {
    M <- scale(as.matrix(M), scale = FALSE)
    do.call(cbind, lapply(seq_len(ncol(M)), function(j) {
      cube <- W %*% to_cube(M[, j])                       # spatial whitening of every time
      as.vector(t(ar_filter(t(cube), phi)))               # AR filtering of every site
    }))
  }
}

# Scale mode: projection of every variable on three bands of the candidate basis.
.correct_scale <- function(X, grp, family, coords, time, groups, bases, kn) {
  if (family == "spatiotemporal") stop("The scale mode is available for spatial and temporal data only.", call. = FALSE)
  cand <- if (family == "spatial") .spatial_candidates(coords) else {
    ord <- order(time); m <- .temporal_candidates(time[ord], "ar"); m[ord, ] <- m; m
  }
  bands <- split(seq_len(ncol(cand)), cut(seq_len(ncol(cand)), 3, labels = c("broad", "medium", "fine")))
  Xc <- scale(X, scale = FALSE)
  R_scale <- lapply(bands, function(b) {
    P <- cand[, b, drop = FALSE] %*% crossprod(cand[, b, drop = FALSE], Xc)
    r <- suppressWarnings(stats::cor(P)); dimnames(r) <- list(colnames(X), colnames(X)); r
  })
  structure(list(R = NULL, R_pairwise = NULL, adjustment = NA_real_, X = NULL, R_scale = R_scale, bases = bases,
                 kn = kn, union_kn = NULL, warnings = character(0), mode = "scale", group = grp),
            class = "depprep_corrected")
}
