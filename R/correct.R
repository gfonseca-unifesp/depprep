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
#' * `"conditional_ar"` (space-time; decision D29): the residuals of the conditional mode,
#'   pair by pair on the same bases, then an AR filter of each site with the coefficients
#'   of each variable's group, for `R`. It removes the persistence in time of what the basis
#'   leaves, which inflates the correlations of the conditional mode, and keeps the spatial
#'   part of the association, which the whitening removes. The filter is applied to `R`
#'   only: the data matrix `X` holds the residuals, without it;
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
#' mode with the global basis, the simplest; spatio-temporal data use the conditional
#' mode followed by the AR filter, with the groups (decision D29, which revises D7). When
#' the association between variables is itself autocorrelated and each variable has noise
#' of its own, as with satellite products, the whitening removed the association, the
#' conditional mode inflated it, and this mode erred least (`analysis/09f`).
#'
#' **Missing values** (`NA` in `X`; decision D26, no imputation). Each variable uses its
#' observed rows. The bases are selected with the missing values counting as the
#' variable's mean ([build_group_basis()]); each variable is then regressed on its basis
#' (and the constant) over its observed rows, so that two variables of a group keep the
#' same basis but each is projected on the rows where it was observed (when their gaps
#' differ, their projections differ a little: the rule "same basis, same projection" then
#' holds only approximately). In the whitening mode, each time (in space and space-time)
#' is whitened at its observed sites, with the covariance of those sites, and an AR filter
#' gives a missing value where the value or one of its lags is missing. The correlations of
#' `R` are computed on the rows where both variables have a value (pairwise), and the
#' corrected `X` keeps the missing values in place. The scale mode needs complete data.
#'
#' @param X Numeric matrix (observations x variables), or a `depprep_sim` object; may hold
#'   missing values (`NA`, see above).
#' @param grouping Output of [group_regimes()], or a named character vector giving the
#'   group of each variable (`"none"` for variables without dependence).
#' @param coords,time,index Coordinates, times or spatio-temporal index, as in
#'   [diagnose_variables()].
#' @param bases Optional named list with one basis per group (from
#'   [build_group_basis()]), or one basis per dependent variable (named by variable)
#'   when `scope = "variables"`; built when `NULL`.
#' @param mode `"conditional"`, `"conditional_ar"`, `"whiten"` or `"scale"`; `NULL` for the default of the
#'   data family (see above).
#' @param scope `"global"`, `"groups"` or `"variables"`; `NULL` for the default
#'   (`"global"` for the conditional and scale modes in spatial and temporal data,
#'   `"groups"` otherwise and whenever `bases` is given).
#' @param data_basis Basis of each column of the data matrix `X`: `"own"` (default; a
#'   basis selected on each variable alone, its selection under `scope = "variables"`,
#'   in any mode; when `bases` are given, the basis given for its group) or `"union"`
#'   (the union of all bases, the same for every variable; conditional modes).
#' @param period Seasonal period for the seasonal groups: one number, or a vector named
#'   by group. When `NULL`, the periods recorded in `grouping` by [group_regimes()] are
#'   used; a seasonal group without any period gets 12, with a warning.
#' @param kn_max Share of degrees of freedom above which a union of bases is flagged;
#'   under `scope = "variables"`, the cap of the union of the variables' bases.
#' @param seed Seed of the basis selection (when `bases` is `NULL`).
#' @param menu Optional menu of structures for the bases ([basis_menu()]): one menu for
#'   every group, or a list of menus named by group (groups without one keep the
#'   automatic choice). Under `scope = "global"` with a list, each group's basis is built
#'   with its menu and the global basis is their union. Not used by the scale mode.
#' @param cores Number of processes that build the bases at the same time (the bases of
#'   the groups, and the basis of each variable for the data matrix). Each basis has its
#'   own seed, so the result is the same with any number of cores. Forked processes on
#'   Linux and macOS, a socket cluster on Windows.
#' @param sample_sites Share of the sites on which the terms of every basis are selected
#'   (spatial and spatio-temporal data), passed to [build_group_basis()]; 1 (the
#'   default) uses all of them. Experimental and off by default: it saves time but the
#'   bases capture less of the large-scale structure (see [build_group_basis()]).
#' @param covariates Optional drivers to condition on (decision D15): a numeric matrix or
#'   data frame of drivers that vary in time only (river discharge, climate indices),
#'   with one row per time step (temporal and spatio-temporal data). The rows are matched
#'   to the times of the data by a column named `time` when there is one, and are
#'   otherwise taken in increasing order of the distinct times. Lags are columns: for the
#'   discharge of the same month and of the two months before, give three columns (say
#'   `discharge`, `discharge_lag1` and `discharge_lag2`), taken from the series before
#'   the first time where needed, since no value may be missing. The drivers are added
#'   to the basis of every group, and are the basis of the variables without dependence,
#'   so that the correction also removes what they explain from every variable; in
#'   space-time each driver takes the same value at every site. The bases are still
#'   selected on the data alone. In the whitening mode the drivers are removed from the
#'   data before the whitening. Not used by the scale mode.
#' @return An object of class `depprep_corrected`: a list with `R`, `R_pairwise`
#'   (before the projection), `adjustment`, `X` (residuals or whitened data; `NULL` in
#'   scale mode), `R_scale`, `r2_scale`, `p_scale`, `bands` and `scalogram` (scale
#'   mode), `bases`, `variable_bases` (scope `"variables"`), `kn` (by group),
#'   `union_kn` (by pair of groups), `warnings`, `mode`, `scope`, `data_basis` and
#'   `group` (`"global"` for every dependent variable under the global scope, `"union"`
#'   under the variables scope). With `covariates`, also `covariates` (the drivers as an
#'   orthonormal basis with one row per observation; `bases` then include them) and
#'   `partition`, the variation partition of each variable between its position (the
#'   basis of its group, without the drivers; none for the variables without dependence)
#'   and the drivers, from [variation_partition()].
#' @references Borcard, D., Legendre, P., Avois-Jacquet, C. and Tuomisto, H. (2004).
#'   Dissecting the spatial structure of ecological data at multiple scales. Ecology 85,
#'   1826-1832.
#'
#'   Peres-Neto, P. R. and Legendre, P. (2010). Estimating and controlling for spatial
#'   structure in the study of ecological communities. Global Ecology and Biogeography
#'   19, 174-184.
#' @examples
#' \donttest{
#' g <- sim_groups(sim_group("long", 4, "matern", range = 0.5), sim_group("none", 2, "none"))
#' s <- sim_spatial(300, g, R_true = make_R_true(g, "blocks", within = 0.5, between = 0), seed = 1)
#' gr <- group_regimes(diagnose_variables(s, nsim = 99))
#' cr <- correct(s, gr)
#' round(cr$R, 2)
#' }
#' @export
correct <- function(X, grouping, coords = NULL, time = NULL, index = NULL, bases = NULL,
                    mode = NULL, scope = NULL, data_basis = c("own", "union"), period = NULL,
                    kn_max = 0.3, seed = 1L, menu = NULL, cores = 1L, sample_sites = 1, covariates = NULL) {
  if (inherits(X, "depprep_sim")) {
    d <- .data_parts(X); X <- d$X; coords <- d$coords; time <- d$time; index <- d$index
  }
  X <- as.matrix(X)
  .check_values(X, coords, time)
  gaps <- anyNA(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  family <- if (!is.null(coords) && !is.null(index)) "spatiotemporal" else if (!is.null(coords)) "spatial" else
    if (!is.null(time)) "temporal" else stop("Give `coords`, `time` or `coords` and `index`.", call. = FALSE)
  if (is.null(mode)) mode <- if (family == "spatiotemporal") "conditional_ar" else "conditional"
  mode <- match.arg(mode, c("conditional", "conditional_ar", "whiten", "scale"))
  if (mode == "conditional_ar" && family != "spatiotemporal") stop("The conditional_ar mode is for space-time data.", call. = FALSE)
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
  if (!is.null(menu)) {
    if (mode == "scale") stop("`menu` is not used by the scale mode.", call. = FALSE)
    if (!inherits(menu, "depprep_menu") && !(is.list(menu) && !is.null(names(menu)) &&
                                              all(vapply(menu, inherits, TRUE, "depprep_menu")))) {
      stop("`menu` must be a basis_menu() or a list of them named by group.", call. = FALSE)
    }
  }
  if (!is.null(covariates) && mode == "scale") stop("`covariates` is not used by the scale mode.", call. = FALSE)
  if (gaps && mode == "scale") stop("The scale mode needs complete data: `X` has missing values.", call. = FALSE)
  C <- if (is.null(covariates)) NULL else .covariate_basis(covariates, family, time, index)
  menu_of <- function(g) if (is.null(menu)) NULL else if (inherits(menu, "depprep_menu")) menu else menu[[g]]
  args <- switch(family, spatial = list(coords = coords), temporal = list(time = time),
                 spatiotemporal = list(coords = coords, index = index))
  if (sample_sites < 1 && family != "temporal") args$sample_sites <- sample_sites
  force(seed); force(period)                               # the processes get values, not names to look up
  pool <- .pool(cores)                                     # processes made on first use, closed on exit
  on.exit(.pool_close(pool), add = TRUE)
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
        parts <- if (is.null(menu)) {
          .lapply_cores(types, function(t) do.call(build_group_basis, c(list(X_group = X[, dep, drop = FALSE], type = t,
            period = .period_of(t, NULL, period, grouping), seed = seed), args)), pool)
        } else if (inherits(menu, "depprep_menu")) {
          list(do.call(build_group_basis, c(list(X_group = X[, dep, drop = FALSE], menu = menu, seed = seed), args)))
        } else .lapply_cores(setdiff(unique(grp), "none"), function(g) do.call(build_group_basis, c(list(X_group = X[, grp == g, drop = FALSE],
          type = gtype[[g]], period = .period_of(gtype[[g]], g, period, grouping), menu = menu_of(g), seed = seed), args)), pool)
        bases <- list(global = suppressWarnings(do.call(union_basis, c(parts, kn_max = kn_max))))
      }
      label <- "global"
    } else {
      variable_bases <- if (is.null(bases)) {
        stats::setNames(.lapply_cores(which(dep), function(j) do.call(build_group_basis, c(list(X_group = X[, j, drop = FALSE],
          type = vtype[j], period = .period_of(vtype[j], grp[j], period, grouping), menu = menu_of(grp[j]), seed = seed), args)), pool),
          colnames(X)[dep])
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
    bases <- stats::setNames(.lapply_cores(groups, function(g) {
      do.call(build_group_basis, c(list(X_group = X[, grp == g, drop = FALSE], type = gtype[[g]],
                                        period = .period_of(gtype[[g]], g, period, grouping), menu = menu_of(g), seed = seed), args))
    }, pool), groups)
  }
  if (is.null(parts)) parts <- bases[groups]
  position <- bases                                        # the bases before the drivers (for the partition)
  if (!is.null(C)) {
    bases[groups] <- lapply(bases[groups], .with_covariates, C)
    for (g in groups) if (attr(bases[[g]], "kn") > kn_max) {
      warn <- c(warn, sprintf("basis of %s with the drivers: k/n = %.2f > %.2f", g, attr(bases[[g]], "kn"), kn_max))
    }
  }
  empty <- if (is.null(C)) .as_basis(matrix(0, n, 0), "none", n, 0, character(0)) else C
  basis_of <- function(g) if (g == "none") empty else bases[[g]]
  kn <- vapply(groups, function(g) attr(bases[[g]], "kn"), 0)

  if (mode == "scale") {
    out <- .correct_scale(X, grp, family, coords, time, bases, parts, kn)
    out$warnings <- c(warn, out$warnings)
    out$scope <- scope
    out$variable_bases <- variable_bases
    return(out)
  }

  # transformed data per group: residuals on a basis (then filtered in time in the
  # conditional_ar mode), or whitened observations
  transform <- if (mode != "whiten") {
    function(cols, B) .resid_on(X[, cols, drop = FALSE], B)
  } else {
    Xw <- if (is.null(C)) X else .resid_on(X, C)             # the drivers are removed before the whitening
    W <- lapply(stats::setNames(c(groups, "none"), c(groups, "none")), function(g)
      .whitener(Xw[, grp == g, drop = FALSE], g, gtype[[g]], bases[[g]], family, coords, time, index))
    function(cols, B, g) W[[g]](Xw[, cols, drop = FALSE])
  }
  if (mode == "conditional_ar") {
    n_s <- nrow(as.matrix(coords))
    Xd <- if (is.null(C)) X else .resid_on(X, C)           # the AR models do not see the drivers
    phi <- lapply(stats::setNames(groups, groups), function(g) .st_phi(Xd[, grp == g, drop = FALSE], index, n_s))
    in_time <- function(E, g) .st_ar_filter(E, phi[[g]], index, n_s)
  }
  ug <- c(groups, if (any(grp == "none")) "none")
  R <- diag(ncol(X)); dimnames(R) <- list(colnames(X), colnames(X))
  union_kn <- NULL
  for (a in seq_along(ug)) for (b in a:length(ug)) {
    ga <- ug[a]; gb <- ug[b]
    ia <- which(grp == ga); ib <- which(grp == gb)
    if (mode == "whiten") {
      Ea <- transform(ia, NULL, ga); Eb <- transform(ib, NULL, gb)
    } else {
      U <- if (ga == gb) basis_of(ga) else suppressWarnings(union_basis(basis_of(ga), basis_of(gb), kn_max = kn_max))
      Ea <- transform(ia, U); Eb <- transform(ib, U)
      if (ga != gb) {
        union_kn <- rbind(union_kn, data.frame(group_a = ga, group_b = gb, k = attr(U, "k"), kn = attr(U, "kn")))
        if (attr(U, "kn") > kn_max) warn <- c(warn, sprintf("union of %s and %s: k/n = %.2f > %.2f", ga, gb, attr(U, "kn"), kn_max))
      }
      if (mode == "conditional_ar") { Ea <- in_time(Ea, ga); Eb <- in_time(Eb, gb) }
    }
    if (mode != "conditional") {
      m <- min(nrow(Ea), nrow(Eb))                       # AR filters drop the first observations
      Ea <- utils::tail(Ea, m); Eb <- utils::tail(Eb, m)
    }
    R[ia, ib] <- if (gaps || mode == "conditional_ar") stats::cor(Ea, Eb, use = "pairwise.complete.obs") else stats::cor(Ea, Eb)
    R[ib, ia] <- t(R[ia, ib])
  }
  diag(R) <- 1
  if (anyNA(R)) {
    bad <- which(is.na(R) & upper.tri(R), arr.ind = TRUE)
    stop(sprintf("Too few rows where both variables have a value to estimate their correlation: %s.",
                 paste(colnames(X)[bad[, 1]], colnames(X)[bad[, 2]], sep = " and ", collapse = "; ")), call. = FALSE)
  }
  near <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  dimnames(near) <- dimnames(R)
  Xout <- if (data_basis == "own") {
    # each variable on a basis selected on that variable alone (decisions D13 and D17), in
    # any mode: R keeps the mode's estimate, the data matrix keeps each variable's classes
    if (is.null(variable_bases) && !bases_given) {
      variable_bases <- stats::setNames(.lapply_cores(which(grp0 != "none"), function(j) do.call(build_group_basis, c(list(
        X_group = X[, j, drop = FALSE], type = vtype0[j], period = .period_of(vtype0[j], grp0[j], period, grouping),
        menu = menu_of(grp0[j]), seed = seed + j), args)), pool), colnames(X)[grp0 != "none"])
    }
    own <- function(j) if (grp0[j] == "none") empty else
      if (!is.null(variable_bases)) .with_covariates(variable_bases[[colnames(X)[j]]], C) else basis_of(grp[j])   # bases given: use them
    E <- vapply(seq_len(ncol(X)), function(j) .resid_on(X[, j, drop = FALSE], own(j))[, 1], numeric(n))
    E <- matrix(E, n); colnames(E) <- colnames(X); E
  } else if (mode != "whiten") {
    U_all <- if (length(groups)) suppressWarnings(do.call(union_basis, c(unname(bases), kn_max = kn_max))) else empty
    E <- .resid_on(X, U_all); colnames(E) <- colnames(X); E
  } else {
    pieces <- lapply(ug, function(g) { cols <- which(grp == g); E <- transform(cols, NULL, g); colnames(E) <- colnames(X)[cols]; E })
    m <- min(vapply(pieces, nrow, 0L))
    do.call(cbind, lapply(pieces, utils::tail, m))[, colnames(X), drop = FALSE]
  }
  out <- list(R = near, R_pairwise = R, adjustment = sqrt(sum((near - R)^2)), X = Xout, R_scale = NULL,
              bases = bases, variable_bases = variable_bases, kn = kn, union_kn = union_kn, warnings = warn,
              mode = mode, scope = scope, data_basis = data_basis, group = grp)
  if (!is.null(C)) {
    out$covariates <- C
    out$partition <- variation_partition(X, stats::setNames(lapply(grp, function(g) if (g == "none") NULL else position[[g]]),
                                                            colnames(X)), C)
  }
  structure(out, class = "depprep_corrected")
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

# The drivers of correct() as an orthonormal basis with one row per observation: one row
# per time step, matched by a `time` column or in increasing order of the times, and
# repeated over the sites in space-time.
.covariate_basis <- function(covariates, family, time, index) {
  if (family == "spatial") stop("`covariates` are drivers that vary in time: give `time`, or `coords` and `index`.", call. = FALSE)
  tt <- if (family == "temporal") time else index$time
  steps <- sort(unique(tt))
  if (is.data.frame(covariates) && "time" %in% names(covariates)) {
    at <- match(steps, covariates$time)
    if (anyNA(at)) stop("`covariates` has no row for some times of the data.", call. = FALSE)
    covariates <- covariates[at, setdiff(names(covariates), "time"), drop = FALSE]
  } else if (NROW(covariates) != length(steps)) {
    stop(sprintf("`covariates` must have one row per time step (%d), or a `time` column.", length(steps)), call. = FALSE)
  }
  D <- .driver_matrix(covariates, length(steps))
  Q <- .orthonormal(D[match(tt, steps), , drop = FALSE])
  if (!ncol(Q)) stop("No column of `covariates` varies in time.", call. = FALSE)
  if (ncol(Q) >= length(steps) - 1) stop("`covariates` has as many columns as time steps.", call. = FALSE)
  .as_basis(Q, "covariates", length(tt), NA_real_, colnames(Q))
}

# A basis together with the drivers (none: the basis itself), keeping its type.
.with_covariates <- function(b, C) {
  if (is.null(C)) return(b)
  u <- suppressWarnings(union_basis(b, C, kn_max = 1))
  attr(u, "type") <- attr(b, "type")
  u
}

# Centred residuals of the columns of M on the basis B. A column with missing values is
# regressed on B and the constant over its observed rows (the columns with the same gaps
# together); its missing values stay missing.
.resid_on <- function(M, B) {
  M <- as.matrix(M)
  B <- unclass(B)
  if (!anyNA(M)) {
    M <- scale(M, scale = FALSE)
    if (!ncol(B)) return(M)
    return(M - B %*% crossprod(B, M))
  }
  out <- M
  for (s in .gap_sets(M)) {
    out[s$rows, s$cols] <- stats::lm.fit(cbind(1, B[s$rows, , drop = FALSE]), M[s$rows, s$cols, drop = FALSE])$residuals
  }
  out
}

# The columns of M grouped by their missing values: for each set, the columns and the rows
# they have observed.
.gap_sets <- function(M) {
  obs <- !is.na(M)
  key <- apply(obs, 2, function(o) paste(which(!o), collapse = " "))
  lapply(unique(key), function(k) { cols <- which(key == k); list(cols = cols, rows = obs[, cols[1]]) })
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
  gaps <- anyNA(Z)
  # spatial whitening: the pooled variogram model of the group as the covariance S of the
  # sites (each variable fitted at its observed sites), and W the inverse of its Cholesky
  # factor. A column with missing values is whitened at its observed sites, with the
  # covariance of those sites (one factor per set of sites, kept for the other columns).
  inv_chol <- function(S) solve(t(chol(S + diag(1e-8 * mean(diag(S)), nrow(S)))))
  spatial_W <- function(cc, zz) {
    zz <- as.matrix(zz)
    fits <- lapply(seq_len(ncol(zz)), function(j) {
      o <- !is.na(zz[, j])
      .fit_variogram(if (all(o)) cc else as.matrix(cc)[o, , drop = FALSE], zz[o, j])
    })
    ok <- vapply(fits, function(f) !is.null(f$fit), TRUE)
    if (!any(ok)) return(list(W = diag(nrow(zz)), S = NULL))
    fit <- fits[[which(ok)[which.min(abs(vapply(fits[ok], `[[`, 0, "range") -
                                           stats::median(vapply(fits[ok], `[[`, 0, "range"))))]]]$fit
    d <- as.matrix(stats::dist(as.matrix(cc)))
    S <- gstat::variogramLine(fit, dist_vector = d, covariance = TRUE)
    list(W = inv_chol(S), S = S, cache = new.env())
  }
  spatial_whiten <- function(sw, M) {
    if (!anyNA(M)) return(sw$W %*% M)
    if (is.null(sw$S)) return(M)
    for (s in .gap_sets(M)) {
      if (!any(s$rows)) next
      k <- paste("gaps", paste(which(!s$rows), collapse = " "))
      if (is.null(sw$cache[[k]])) sw$cache[[k]] <- inv_chol(sw$S[s$rows, s$rows, drop = FALSE])
      M[s$rows, s$cols] <- sw$cache[[k]] %*% M[s$rows, s$cols, drop = FALSE]
    }
    M
  }
  if (family == "spatial") {
    sw <- spatial_W(coords, Z)
    return(function(M) spatial_whiten(sw, scale(as.matrix(M), scale = FALSE)))
  }
  if (family == "temporal") {
    ord <- order(time)
    det <- if (type %in% c("seasonal", "trend")) basis else NULL
    pre <- function(M) { M <- as.matrix(M)[ord, , drop = FALSE]; if (is.null(det)) scale(M, scale = FALSE) else .resid_on(M, unclass(det)[ord, , drop = FALSE]) }
    phi <- .ar_coef(pre(Z))
    return(function(M) .ar_filter(pre(M), phi))
  }
  # space-time: rows ordered by index; whiten sites at each time, then filter each site
  n_s <- nrow(as.matrix(coords)); tt <- match(index$time, sort(unique(index$time))); n_t <- max(tt)
  to_cube <- function(v) { m <- matrix(0, n_s, n_t); m[cbind(index$site, tt)] <- v; m }
  site_means <- do.call(cbind, lapply(seq_len(ncol(Z)), function(j) rowMeans(to_cube(Z[, j]), na.rm = gaps)))
  site_means[is.nan(site_means)] <- NA                    # a site never observed in a variable
  sw <- spatial_W(coords, site_means)
  phi <- .st_phi(Xg, index, n_s)
  function(M) {
    M <- scale(as.matrix(M), scale = FALSE)
    do.call(cbind, lapply(seq_len(ncol(M)), function(j) {
      cube <- spatial_whiten(sw, to_cube(M[, j]))         # spatial whitening of every time
      as.vector(t(.ar_filter(t(cube), phi)))              # AR filtering of every site
    }))
  }
}

# AR coefficients pooled over the columns; with gaps the AR fit skips them (and the orders
# whose AIC cannot be computed are not chosen).
.ar_coef <- function(series) {
  na <- if (anyNA(series)) stats::na.pass else stats::na.fail
  fit <- function(x, ...) suppressWarnings(stats::ar(x, na.action = na, ...))
  orders <- vapply(seq_len(ncol(series)), function(j) fit(series[, j], aic = TRUE, order.max = 5)$order, 0)
  p <- max(1, round(stats::median(orders)))
  coefs <- vapply(seq_len(ncol(series)), function(j) fit(series[, j], aic = FALSE, order.max = p)$ar, numeric(p))
  matrix(coefs, nrow = p) |> rowMeans()
}

# Each column of M (rows in time order) filtered with the AR coefficients phi, without the
# first rows. A value is missing where it or one of its lags is missing.
.ar_filter <- function(M, phi) {
  p <- length(phi)
  M <- as.matrix(M)
  out <- M[-(1:p), , drop = FALSE]
  for (i in seq_len(p)) out <- out - phi[i] * M[(p + 1 - i):(nrow(M) - i), , drop = FALSE]
  out
}

# The values of one space-time variable as sites x times, NA where there is no row.
.to_cube <- function(v, index, n_s) {
  tt <- match(index$time, sort(unique(index$time)))
  m <- matrix(NA_real_, n_s, max(tt))
  m[cbind(index$site, tt)] <- v
  m
}

# AR coefficients of a space-time group, from the mean over the sites of each of its
# standardised columns.
.st_phi <- function(Xg, index, n_s) {
  Z <- scale(Xg)
  means <- do.call(cbind, lapply(seq_len(ncol(Z)), function(j) colMeans(.to_cube(Z[, j], index, n_s), na.rm = TRUE)))
  means[is.nan(means)] <- NA                               # a time never observed in a variable
  .ar_coef(scale(means, scale = FALSE))
}

# The columns of M filtered in time, site by site, with the coefficients phi (`NULL`: none).
# Rows sorted by time and by site within a time, without the first times.
.st_ar_filter <- function(M, phi, index, n_s) {
  M <- as.matrix(M)
  do.call(cbind, lapply(seq_len(ncol(M)), function(j) {
    cube <- .to_cube(M[, j], index, n_s)
    as.vector(if (length(phi)) t(.ar_filter(t(cube), phi)) else cube)
  }))
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

# A pool of `cores` processes for .lapply_cores(), made on first use and closed by
# .pool_close(): forked processes where they exist; on Windows one socket cluster per
# pool (its processes are slow to start) that loads this package (from its source folder
# when it was loaded with pkgload, as in development).
.pool <- function(cores = 1L) {
  p <- new.env(parent = emptyenv())
  p$cores <- max(1L, as.integer(cores)); p$cl <- NULL
  p
}
.pool_close <- function(p) {
  if (!is.null(p$cl)) parallel::stopCluster(p$cl)
  p$cl <- NULL
  invisible(NULL)
}

# lapply() over `x` with a pool (or a number of cores, for a pool used once). Every basis
# carries its own seed (through .with_stream()), so the result does not depend on the
# number of processes.
.lapply_cores <- function(x, FUN, pool = 1L) {
  if (!is.environment(pool)) { pool <- .pool(pool); on.exit(.pool_close(pool)) }
  cores <- min(pool$cores, length(x))
  if (cores <= 1L) return(lapply(x, FUN))
  if (.Platform$OS.type != "windows") {
    out <- parallel::mclapply(x, FUN, mc.cores = cores, mc.preschedule = FALSE)
    bad <- vapply(out, inherits, TRUE, "try-error")
    if (any(bad)) stop(conditionMessage(attr(out[[which(bad)[1]]], "condition")), call. = FALSE)
    return(out)
  }
  if (is.null(pool$cl)) {
    pool$cl <- parallel::makePSOCKcluster(pool$cores, setup_timeout = 600)
    ns_path <- getNamespaceInfo(asNamespace("depprep"), "path")
    dev <- isNamespaceLoaded("pkgload") && isTRUE(pkgload::is_dev_package("depprep"))
    parallel::clusterCall(pool$cl, function(path, dev) {
      if (dev) pkgload::load_all(path, quiet = TRUE, export_all = FALSE, helpers = FALSE) else loadNamespace("depprep")
      NULL
    }, ns_path, dev)
  }
  parallel::clusterApplyLB(pool$cl, x, FUN)
}
