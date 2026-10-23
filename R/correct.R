#' Correct a data set for the dependence of each group
#'
#' Removes the large-scale structure of each group of variables with the group's own
#' basis, so that what remains (the residual components) can be compared across
#' variables. Two outputs, as fixed by the project's decision P1:
#' * `R`: a correlation matrix assembled pair by pair — within a group with the group's
#'   basis, across groups with the union of the two bases — and projected to the
#'   nearest valid correlation matrix (`adjustment` records how much it moved); for
#'   PCA and other methods that take a correlation matrix;
#' * `X`: a data matrix for k-means, SOM and MDS on observations. By default
#'   (`data_basis = "own"`, decisions D13 and D17) each variable is corrected with a
#'   basis selected on that variable alone, in every mode and family: in the stage 6
#'   pilot this removed the spurious patterns almost as well as a shared basis and,
#'   unlike it, kept the classes of the observations, because a large shared basis
#'   removes part of every variable's class signal. The columns of `X` then carry the
#'   small between-basis error of their correlations (of order rho k/n); `R` is the
#'   estimate of the correlations. `data_basis = "union"` gives the residuals on the
#'   union of all the bases instead (one basis for every variable).
#'
#' Modes:
#' * `"conditional"` (residuals on the basis): `(I - P) x`;
#' * `"whiten"`: the observations are whitened with the covariance estimated for each
#'   group — a variogram model pooled over the group's variables in space, an AR model
#'   in time (after removing the trend or seasonal terms of the group's basis), and both
#'   (separably) in space-time. Uses the whole covariance, not only a basis, which can
#'   matter for rough processes. Pairs across groups use each group's own whitening;
#' * `"scale"`: one correlation matrix per scale band, from the projection of every
#'   variable on the band; spatial and temporal data only. The bands are made of the
#'   functions the bases selected (Borcard et al. 2004): the trend terms form their own
#'   band, and the other functions, ordered from broad to fine (the order of the Moran
#'   eigenvector maps in space; the wavelength in time, read from the lag-one
#'   autocorrelation of each function, which assumes a regular series), are split into
#'   up to three bands (broad, medium, fine) at gaps of the scalogram, a gap being a run
#'   of candidates that were not selected across which the scale changes by at least a
#'   factor 1.25. Among the gaps, the cuts are those that leave the least spread of
#'   log-scale within the bands, each function weighted by the variance it explains, so
#'   that a function with little signal does not make a band of its own. The
#'   share of variance of each variable in each band (`r2_scale`) is tested with an F
#'   test given all selected functions (`p_scale`); bands where no variable has signal
#'   (Holm adjustment over the variables) are flagged in `warnings`.
#'
#' Scope of the bases: `"groups"` (one basis per group, unions across groups),
#' `"global"` (one basis for all dependent variables: the union of the bases selected on
#' all of them for each regime type present) or `"variables"` (each dependent variable
#' selects its own basis, and the union of these bases, capped at `kn_max`, is the
#' shared basis of all of them; Peres-Neto and Legendre 2010). With `"variables"`, `R`
#' uses the union for every pair, so the bases of two variables never differ. Variables
#' without dependence are never corrected among themselves.
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
#'   [build_group_basis()]), or one basis per dependent variable (named by variable)
#'   when `scope = "variables"`; built when `NULL`.
#' @param mode `"conditional"`, `"whiten"` or `"scale"`; `NULL` for the default of the
#'   data family (see above).
#' @param scope `"global"`, `"groups"` or `"variables"`; `NULL` for the default
#'   (`"global"` for the conditional and scale modes in spatial and temporal data,
#'   `"groups"` otherwise and whenever `bases` is given).
#' @param data_basis Basis of each column of the data matrix `X`: `"own"` (default; a
#'   basis selected on each variable alone, its selection under `scope = "variables"`,
#'   in any mode; when `bases` are given, the basis given for its group) or `"union"`
#'   (the union of all bases, the same for every variable; conditional mode).
#' @param period Seasonal period for the seasonal groups: one number, or a vector named
#'   by group. When `NULL`, the periods recorded in `grouping` by [group_regimes()] are
#'   used; a seasonal group without any period gets 12, with a warning.
#' @param kn_max Share of degrees of freedom above which a union of bases is flagged;
#'   under `scope = "variables"`, the cap of the union of the variables' bases.
#' @param seed Seed of the basis selection (when `bases` is `NULL`).
#' @return An object of class `depprep_corrected`: a list with `R`, `R_pairwise`
#'   (before the projection), `adjustment`, `X` (residuals or whitened data; `NULL` in
#'   scale mode), `R_scale`, `r2_scale`, `p_scale`, `bands` and `scalogram` (scale
#'   mode), `bases`, `variable_bases` (scope `"variables"`), `kn` (by group),
#'   `union_kn` (by pair of groups), `warnings`, `mode`, `scope`, `data_basis` and
#'   `group` (`"global"` for every dependent variable under the global scope, `"union"`
#'   under the variables scope).
#' @references Borcard, D., Legendre, P., Avois-Jacquet, C. and Tuomisto, H. (2004).
#'   Dissecting the spatial structure of ecological data at multiple scales. Ecology 85,
#'   1826-1832.
#'
#'   Peres-Neto, P. R. and Legendre, P. (2010). Estimating and controlling for spatial
#'   structure in the study of ecological communities. Global Ecology and Biogeography
#'   19, 174-184.
#' @examples
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(300, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 1)
#' gr <- group_regimes(diagnose_variables(s, nsim = 99))
#' cr <- correct(s, gr)
#' round(cr$R, 2)
#' @export
correct <- function(X, grouping, coords = NULL, time = NULL, index = NULL, bases = NULL,
                    mode = NULL, scope = NULL, data_basis = c("own", "union"), period = NULL,
                    kn_max = 0.3, seed = 1L) {
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
  scope <- match.arg(scope, c("groups", "global", "variables"))
  data_basis <- match.arg(data_basis)
  grp <- .group_vector(grouping, colnames(X))
  gtype <- .group_types(grouping, grp)
  # regime type of each variable, before any union: from the grouping table, else from
  # the basis given for its group, else the group's type
  vtype0 <- if (inherits(grouping, "data.frame") && "type" %in% names(grouping)) {
    grouping$type[match(colnames(X), grouping$variable)]
  } else vapply(grp, function(g) {
    if (g == "none") return("none")
    tb <- if (!is.null(bases) && !is.null(bases[[g]])) attr(bases[[g]], "type") else NULL
    if (!is.null(tb)) tb else gtype[[g]]
  }, "")
  grp0 <- grp
  bases_given <- !is.null(bases)
  args <- switch(family, spatial = list(coords = coords), temporal = list(time = time),
                 spatiotemporal = list(coords = coords, index = index))
  n <- nrow(X)
  warn <- character(0)
  parts <- NULL                                            # bases before any union (for the scale mode)
  variable_bases <- NULL
  if (scope %in% c("global", "variables") && any(grp != "none")) {
    dep <- grp != "none"
    vtype <- if (inherits(grouping, "data.frame") && "type" %in% names(grouping)) {
      grouping$type[match(colnames(X), grouping$variable)]
    } else unlist(gtype[grp])
    types <- unique(vtype[dep])
    if (scope == "global") {
      if (is.null(bases)) {
        parts <- lapply(types, function(t) do.call(build_group_basis, c(list(X_group = X[, dep, drop = FALSE], type = t,
          period = .period_of(t, NULL, period, grouping), seed = seed), args)))
        bases <- list(global = suppressWarnings(do.call(union_basis, c(parts, kn_max = kn_max))))
      }
      label <- "global"
    } else {
      variable_bases <- if (is.null(bases)) {
        stats::setNames(lapply(which(dep), function(j) do.call(build_group_basis, c(list(X_group = X[, j, drop = FALSE],
          type = vtype[j], period = .period_of(vtype[j], grp[j], period, grouping), seed = seed), args))), colnames(X)[dep])
      } else {
        if (!all(colnames(X)[dep] %in% names(bases))) {
          stop("With `scope = \"variables\"`, `bases` must name one basis per dependent variable.", call. = FALSE)
        }
        bases[colnames(X)[dep]]
      }
      parts <- variable_bases
      u <- .capped_union(variable_bases, X[, dep, drop = FALSE], kn_max)
      if (isTRUE(attr(u, "capped"))) {
        warn <- c(warn, sprintf("union of the variable bases capped at k/n = %.2f (%d of %d directions kept)",
                                attr(u, "kn"), attr(u, "k"), attr(u, "rank")))
      }
      bases <- list(union = u)
      label <- "union"
    }
    grp[dep] <- label
    gtype <- stats::setNames(list(if (any(types %in% c("seasonal", "trend"))) "seasonal" else types[1], "none"),
                             c(label, "none"))
  }
  groups <- setdiff(unique(grp), "none")
  if (is.null(bases)) {
    bases <- stats::setNames(lapply(groups, function(g) {
      do.call(build_group_basis, c(list(X_group = X[, grp == g, drop = FALSE], type = gtype[[g]],
                                        period = .period_of(gtype[[g]], g, period, grouping), seed = seed), args))
    }), groups)
  }
  if (is.null(parts)) parts <- bases[groups]
  empty <- .as_basis(matrix(0, n, 0), "none", n, 0, character(0))
  basis_of <- function(g) if (g == "none") empty else bases[[g]]
  kn <- vapply(groups, function(g) attr(bases[[g]], "kn"), 0)

  if (mode == "scale") {
    out <- .correct_scale(X, grp, family, coords, time, bases, parts, kn)
    out$warnings <- c(warn, out$warnings)
    out$scope <- scope
    out$variable_bases <- variable_bases
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
  Xout <- if (data_basis == "own") {
    # each variable on a basis selected on that variable alone (decisions D13 and D17), in
    # any mode: R keeps the mode's estimate, the data matrix keeps each variable's classes
    if (is.null(variable_bases) && !bases_given) {
      variable_bases <- stats::setNames(lapply(which(grp0 != "none"), function(j) do.call(build_group_basis, c(list(
        X_group = X[, j, drop = FALSE], type = vtype0[j], period = .period_of(vtype0[j], grp0[j], period, grouping),
        seed = seed + j), args))), colnames(X)[grp0 != "none"])
    }
    own <- function(j) if (grp0[j] == "none") empty else
      if (!is.null(variable_bases)) variable_bases[[colnames(X)[j]]] else basis_of(grp[j])   # bases given: use them
    E <- vapply(seq_len(ncol(X)), function(j) .resid_on(X[, j, drop = FALSE], own(j))[, 1], numeric(n))
    E <- matrix(E, n); colnames(E) <- colnames(X); E
  } else if (mode == "conditional") {
    U_all <- if (length(groups)) suppressWarnings(do.call(union_basis, c(unname(bases), kn_max = kn_max))) else empty
    E <- .resid_on(X, U_all); colnames(E) <- colnames(X); E
  } else {
    pieces <- lapply(ug, function(g) { cols <- which(grp == g); E <- transform(cols, NULL, g); colnames(E) <- colnames(X)[cols]; E })
    m <- min(vapply(pieces, nrow, 0L))
    do.call(cbind, lapply(pieces, utils::tail, m))[, colnames(X), drop = FALSE]
  }
  structure(list(R = near, R_pairwise = R, adjustment = sqrt(sum((near - R)^2)), X = Xout, R_scale = NULL,
                 bases = bases, variable_bases = variable_bases, kn = kn, union_kn = union_kn, warnings = warn,
                 mode = mode, scope = scope, data_basis = data_basis, group = grp),
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

# Seasonal period of a group (`g`, or NULL for one basis of all groups): the `period`
# argument first, then the periods recorded by group_regimes(); 12 with a warning when
# none is known. NULL for the types that need no period.
.period_of <- function(type, g, period, grouping) {
  if (!type %in% c("seasonal", "separable_seasonal")) return(NULL)
  known <- if (!is.null(period)) period else if (inherits(grouping, "depprep_grouping")) attr(grouping, "period")
  p <- if (length(known) == 1 && is.null(names(known))) known[[1]] else
    if (!is.null(g) && g %in% names(known)) known[[g]] else
      if (is.null(g) && any(is.finite(known))) round(stats::median(known[is.finite(known)])) else NA
  if (!is.finite(p)) {
    warning(sprintf("No seasonal period known for %s; using 12. Give `period` to set it.",
                    if (is.null(g)) "the global basis" else paste("group", g)), call. = FALSE)
    p <- 12
  }
  p
}

# Union of the bases of the variables, capped at kn_max: the columns that explain most of
# the (standardised) variables come first, so the cap drops the least useful directions.
.capped_union <- function(bases, Y, kn_max) {
  n <- nrow(Y)
  M <- do.call(cbind, c(list(matrix(0, n, 0)), lapply(unname(bases), function(b) as.matrix(unclass(b)))))
  if (!ncol(M)) return(.as_basis(matrix(0, n, 0), "union", n, NA_real_, character(0)))
  terms <- if (is.null(colnames(M))) paste0("u", seq_len(ncol(M))) else colnames(M)
  Z <- scale(as.matrix(Y)); Z[!is.finite(Z)] <- 0
  o <- order(-rowSums(crossprod(M, Z)^2))
  M <- M[, o, drop = FALSE]; terms <- terms[o]
  keep <- .distinct_columns(M, terms)                     # many variables select the same functions
  M <- M[, keep, drop = FALSE]; terms <- terms[keep]
  if (ncol(M) >= n) { M <- M[, seq_len(n - 1), drop = FALSE]; terms <- terms[seq_len(n - 1)] }  # QR needs ncol < n
  q <- qr(M, tol = 1e-7)                                  # dependent directions are dropped here
  k <- min(q$rank, floor(kn_max * n))
  B <- .as_basis(qr.Q(q)[, seq_len(k), drop = FALSE], "union", n, NA_real_, terms[q$pivot[seq_len(k)]])
  attr(B, "rank") <- q$rank
  attr(B, "capped") <- q$rank > k
  B
}

# Columns to keep so that no column repeats an earlier one. Same-named columns are
# compared numerically: a name can stand for different columns (a temporal map is
# orthogonalised against the trend terms only in the "trend" type).
.distinct_columns <- function(M, terms) {
  first <- match(terms, terms)
  keep <- rep(TRUE, ncol(M))
  for (j in which(first != seq_along(terms))) {
    if (max(abs(M[, j] - M[, first[j]])) < 1e-8) keep[j] <- FALSE
  }
  keep
}

# Scale mode: projection of every variable on bands of the selected functions (see
# correct()).
.correct_scale <- function(X, grp, family, coords, time, bases, parts, kn, alpha = 0.05, min_ratio = 1.25) {
  if (family == "spatiotemporal") stop("The scale mode is available for spatial and temporal data only.", call. = FALSE)
  n <- nrow(X)
  warn <- character(0)
  S <- do.call(cbind, c(list(matrix(0, n, 0)), lapply(unname(parts), unclass)))   # selected functions, each once
  if (is.null(colnames(S))) colnames(S) <- paste0("u", seq_len(ncol(S)))
  S <- S[, .distinct_columns(S, colnames(S)), drop = FALSE]
  trend <- grepl("^trend", colnames(S))
  axis <- .scale_axis(colnames(S)[!trend], family, time, unique(unlist(lapply(parts, function(b) {
    p <- sub(".*_p", "", grep("_p[0-9.]+$", attr(b, "terms"), value = TRUE)); as.numeric(p)
  }))))
  if (family == "temporal" && length(unique(round(diff(sort(unique(time))), 8))) > 1) {
    warn <- c(warn, "the series is not regular: the wavelengths of the bands are approximate")
  }
  sel <- colnames(S)[!trend]
  pos <- match(sel, axis$term)
  if (anyNA(pos)) {
    warn <- c(warn, paste("functions without a known scale left out of the bands:", paste(sel[is.na(pos)], collapse = ", ")))
    sel <- sel[!is.na(pos)]; pos <- pos[!is.na(pos)]
  }
  o <- order(pos); sel <- sel[o]; pos <- pos[o]
  ls <- log(axis$scale[pos])
  gaps <- which(diff(pos) > 1 & abs(diff(ls)) >= log(min_ratio))  # unselected candidates between, and a change of scale
  Z <- scale(X[, grp != "none", drop = FALSE]); Z[!is.finite(Z)] <- 0
  w <- rowSums(crossprod(S[, sel, drop = FALSE], Z)^2) + 1e-12      # weight of each function: variance it explains
  spread <- function(cut) {                                # within-band spread of log-scale, weighted
    b <- findInterval(seq_along(ls), cut + 1)
    sum(vapply(split(seq_along(ls), b), function(i) sum(w[i] * (ls[i] - stats::weighted.mean(ls[i], w[i]))^2), 0))
  }
  options <- if (length(gaps) > 2) utils::combn(gaps, 2, simplify = FALSE) else list(gaps)
  cuts <- options[[which.min(vapply(options, spread, 0))]]
  band <- findInterval(seq_along(pos), cuts + 1) + 1
  labels <- list("all", c("broad", "fine"), c("broad", "medium", "fine"))[[max(1, length(cuts) + 1)]]
  members <- c(if (any(trend)) list(trend = colnames(S)[trend]),
               if (length(sel)) split(sel, factor(labels[band], labels)))
  if (!length(members)) warn <- c(warn, "no function was selected: there is no scale band")

  Xc <- scale(X, scale = FALSE)
  tss <- colSums(Xc^2)
  q_all <- qr(S[, unlist(members), drop = FALSE])
  k_all <- q_all$rank
  r2_all <- if (k_all) colSums(crossprod(qr.Q(q_all)[, seq_len(k_all), drop = FALSE], Xc)^2) / tss else 0 * tss
  fits <- lapply(members, function(cols) {
    q <- qr(S[, cols, drop = FALSE]); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    P <- Q %*% crossprod(Q, Xc)
    r2 <- colSums(P^2) / tss
    f <- (r2 / q$rank) / ((1 - r2_all) / (n - 1 - k_all))
    r <- suppressWarnings(stats::cor(P)); dimnames(r) <- list(colnames(X), colnames(X))
    list(R = r, r2 = r2, p = stats::pf(f, q$rank, n - 1 - k_all, lower.tail = FALSE))
  })
  by_band <- function(what) {
    m <- matrix(vapply(fits, `[[`, numeric(ncol(X)), what), ncol(X)); dimnames(m) <- list(colnames(X), names(members)); m
  }
  r2_scale <- by_band("r2"); p_scale <- by_band("p")
  for (b in names(members)) if (!any(stats::p.adjust(p_scale[, b], "holm") < alpha, na.rm = TRUE)) {
    warn <- c(warn, sprintf("band %s: no variable has signal (F test, Holm over the variables, p >= %.2f)", b, alpha))
  }
  wl <- function(cols, f) if (family == "temporal" && !any(grepl("^trend", cols))) f(axis$scale[match(cols, axis$term)]) else NA_real_
  bands <- data.frame(band = as.character(names(members)), k = vapply(members, length, 0L),
                      terms = vapply(members, paste, "", collapse = ", "),
                      wavelength_max = vapply(members, wl, 0, max), wavelength_min = vapply(members, wl, 0, min),
                      row.names = NULL)
  structure(list(R = NULL, R_pairwise = NULL, adjustment = NA_real_, X = NULL,
                 R_scale = lapply(fits, `[[`, "R"), r2_scale = r2_scale, p_scale = p_scale, bands = bands,
                 scalogram = lapply(parts, attr, "scalogram"), bases = bases, kn = kn, union_kn = NULL,
                 warnings = warn, mode = "scale", data_basis = "union", group = grp),
            class = "depprep_corrected")
}

# Candidates ordered from broad to fine, with their scale: for space the dbMEM, whose
# scale is taken as rank^(-1/2) (the wavelength of a planar map falls with the square
# root of its rank); for time the dbMEM and the harmonics of the given periods, with
# the wavelength read from their lag-one autocorrelation (2 pi dt / acos r1).
.scale_axis <- function(terms, family, time, periods) {
  if (family == "spatial") {
    k <- suppressWarnings(as.integer(sub("^MEM", "", grep("^MEM[0-9]+$", terms, value = TRUE))))
    m <- max(c(k, 0))
    return(data.frame(term = paste0("MEM", seq_len(m)), scale = seq_len(m)^-0.5))
  }
  o <- order(time)
  cand <- .temporal_candidates(time[o], "ar", trend = FALSE)
  for (p in periods[is.finite(periods)]) cand <- cbind(cand, .temporal_candidates(time[o], "seasonal", p, trend = FALSE))
  cand <- scale(cand, scale = FALSE)
  r1 <- colSums(cand[-1, , drop = FALSE] * cand[-nrow(cand), , drop = FALSE]) / colSums(cand^2)
  dt <- stats::median(diff(sort(unique(time))))
  wavelength <- 2 * pi * dt / acos(pmin(pmax(r1, -1), 1 - 1e-12))
  axis <- data.frame(term = colnames(cand), scale = wavelength)
  axis[order(-axis$scale), , drop = FALSE]
}
