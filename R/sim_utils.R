# Building blocks shared by sim_spatial(), sim_temporal() and sim_spatiotemporal().
#
# Generator model (simulation protocol, section 3): for variable j in group g,
#   x_j = sqrt(s_g) * L_g + sqrt(1 - s_g) * e_j + class effect,
# where L_g is the standardised large-scale component of the group, e has
# correlation R_true across variables and is independent across observations,
# and s_g is the share of variance that comes from the large scale.

# Independent random streams, one per component, so that changing one component
# (for example the residuals) leaves the others unchanged.
.streams <- c(coords = 1L, large = 2L, residual = 3L, classes = 4L, bootstrap = 5L, permutation = 6L,
              basis = 7L, subsample = 8L, downstream = 9L, counts = 10L)

# Evaluate `expr` with an L'Ecuyer-CMRG stream derived from `seed` and restore the
# caller's random number generator afterwards.
.with_stream <- function(seed, stream, expr) {
  env <- globalenv()
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env)
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- get(".Random.seed", envir = env)
  for (i in seq_len(stream)) s <- parallel::nextRNGStream(s)
  assign(".Random.seed", s, envir = env)
  force(expr)
}

# Exponential correlation with practical range `range` (correlation 0.05 at `range`).
.exp_cor <- function(h, range) exp(-3 * h / range)

# Matern correlation with smoothness nu in {0.5, 1.5, 2.5} and practical range `range`.
# nu = 0.5 is the exponential (rough fields); larger nu gives smoother fields. With
# s = sqrt(2 nu) h / l, the correlation is exp(-s), (1 + s) exp(-s) or
# (1 + s + s^2 / 3) exp(-s); l is chosen so that the correlation is 0.05 at `range`.
.matern_shape <- function(s, nu) {
  switch(as.character(nu), "0.5" = exp(-s), "1.5" = (1 + s) * exp(-s), "2.5" = (1 + s + s^2 / 3) * exp(-s),
         stop("`smoothness` must be 0.5, 1.5 or 2.5.", call. = FALSE))
}
.matern_s0 <- function(nu) {
  if (nu == 0.5) return(3)
  stats::uniroot(function(s) .matern_shape(s, nu) - 0.05, c(0.1, 20), tol = 1e-12)$root
}
.matern_cor <- function(h, range, nu) {
  if (nu == 0.5) return(.exp_cor(h, range))
  .matern_shape(.matern_s0(nu) * h / range, nu)
}

# Geometric anisotropy: rotate so that the major axis (at `angle` degrees from the
# x axis) lies on x, then stretch the minor axis by `ratio`. The practical range is
# `range` along the major axis and `range / ratio` along the minor axis.
.aniso_coords <- function(coords, ratio = 1, angle = 0) {
  th <- angle * pi / 180
  rot <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  xy <- as.matrix(coords) %*% rot
  xy[, 2] <- xy[, 2] * ratio
  xy
}

# Gaussian random field(s) with Matern covariance (exponential by default), by
# Cholesky factorisation. Returns an n x n_sim matrix with unit variance. Smooth
# covariances are ill-conditioned, so they get a larger jitter (a tiny nugget).
.grf_chol <- function(coords, range, ratio = 1, angle = 0, n_sim = 1, smoothness = 0.5) {
  xy <- .aniso_coords(coords, ratio, angle)
  cmat <- .matern_cor(as.matrix(stats::dist(xy)), range, smoothness)
  u <- chol(cmat + diag(if (smoothness == 0.5) 1e-10 else 1e-6, nrow(cmat)))
  t(u) %*% matrix(stats::rnorm(nrow(cmat) * n_sim), nrow(cmat))
}

# Random Fourier features for the Matern covariance. The spectral measure of a Matern
# with smoothness nu and length l (correlation M(sqrt(2 nu) |h| / l)) in two dimensions
# is a multivariate t with 2 nu degrees of freedom and scale 1 / l; for the exponential
# (nu = 0.5) it is a bivariate Cauchy with scale 3 / range. Returns a function that
# evaluates one field realisation at any set of coordinates, which is what a moving
# pattern needs.
.rff_sampler <- function(range, ratio = 1, angle = 0, n_features = 1000, smoothness = 0.5) {
  a <- range / 3
  w <- stats::rchisq(n_features, df = 2 * smoothness)
  z <- matrix(stats::rnorm(2 * n_features), n_features)
  omega <- if (smoothness == 0.5) z / (a * sqrt(w)) else {
    l <- sqrt(2 * smoothness) * range / .matern_s0(smoothness)
    z / (l * sqrt(w / (2 * smoothness)))
  }
  b <- stats::runif(n_features, 0, 2 * pi)
  function(coords) {
    xy <- .aniso_coords(coords, ratio, angle)
    phase <- sweep(xy %*% t(omega), 2, b, "+")
    sqrt(2 / n_features) * rowSums(cos(phase))
  }
}

# Stationary AR(1) series with unit variance; returns an n_times x n_series matrix.
.ar1 <- function(n_times, phi, n_series = 1) {
  stopifnot(abs(phi) < 1)
  eps <- matrix(stats::rnorm(n_times * n_series), n_times)
  x <- eps
  if (n_times > 1) {
    for (t in 2:n_times) x[t, ] <- phi * x[t - 1, ] + sqrt(1 - phi^2) * eps[t, ]
  }
  x
}

# Harmonic with unit variance.
.seasonal <- function(n_times, period, phase = 0) {
  sqrt(2) * sin(2 * pi * seq_len(n_times) / period + phase)
}

# Linear trend with a jump at `shift_at` (fraction of the series), standardised.
.trend_shift <- function(n_times, slope = 1, jump = 1, shift_at = 0.5) {
  t <- seq_len(n_times)
  f <- slope * t / n_times + jump * (t > shift_at * n_times)
  as.vector(scale(f))
}

# Harmonic plus linear trend, standardised. `strength` is the standard deviation of
# the trend relative to that of the harmonic (1: equal variances).
.seasonal_trend <- function(n_times, period, phase = 0, strength = 1) {
  t <- seq_len(n_times)
  as.vector(scale(.seasonal(n_times, period, phase) + strength * as.vector(scale(t))))
}

# Deterministic surface along `direction` (degrees from the x axis): the coordinate
# projected on that direction from the centre of the sites' bounding box ("linear") or
# its square ("quadratic", a trough across the direction), standardised over the
# sites. With a unit-variance `field`, the sum is rescaled to unit variance and
# `strength` is the standard deviation of the surface relative to that of the field.
.gradient <- function(coords, direction = 0, shape = "linear", strength = 1, field = NULL) {
  th <- direction * pi / 180
  ctr <- apply(coords, 2, function(v) mean(range(v)))
  d <- (coords[, 1] - ctr[1]) * cos(th) + (coords[, 2] - ctr[2]) * sin(th)
  s <- switch(shape, linear = d, quadratic = d^2,
              stop("`shape` must be \"linear\" or \"quadratic\".", call. = FALSE))
  s <- as.vector(scale(s))
  if (is.null(field)) return(s)
  (strength * s + field) / sqrt(strength^2 + 1)
}

# Residual components: n_obs x p, rows independent, columns with correlation R.
.residual_mvn <- function(n_obs, R) {
  u <- tryCatch(chol(R), error = function(e) {
    stop("`R_true` must be a positive definite correlation matrix.", call. = FALSE)
  })
  matrix(stats::rnorm(n_obs * ncol(R)), n_obs) %*% u
}

# Contiguous classes: k spatial regions (k-means of the coordinates) or k blocks of
# consecutive times.
.contiguous_classes <- function(k, coords = NULL, n_times = NULL) {
  if (!is.null(coords)) return(stats::kmeans(as.matrix(coords), centers = k, nstart = 10)$cluster)
  as.integer(ceiling(seq_len(n_times) * k / n_times))
}

# Balanced random classes.
.random_classes <- function(k, n_obs) sample(rep_len(seq_len(k), n_obs))

# Class means: one row per class, non-zero only in `vars`; centred over classes.
.class_effects <- function(k, p, effect, vars) {
  mu <- matrix(0, k, p)
  if (effect != 0 && length(vars)) {
    m <- matrix(stats::rnorm(k * length(vars)) * effect, k)
    mu[, vars] <- sweep(m, 2, colMeans(m))
  }
  mu
}

#' Define one dependence group for the simulators
#'
#' A group is a set of variables that share the same large-scale component (the
#' same dependence regime). Regime parameters are passed in `...` and depend on the
#' simulator: see [sim_spatial()], [sim_temporal()] and [sim_spatiotemporal()].
#'
#' @param name Group name (used to name the variables).
#' @param n_vars Number of variables in the group.
#' @param regime Regime name, for example `"exponential"`, `"ar1"` or `"none"`.
#' @param ... Regime parameters, for example `range = 0.25` or `phi = 0.9`.
#' @param large_scale_share Share of each variable's variance that comes from the
#'   large-scale component (between 0 and 1). Forced to 0 for regime `"none"`.
#' @param shared_with Name of another group whose large-scale component this group
#'   reuses. The regime of this group is then ignored.
#' @return An object of class `depprep_group`.
#' @examples
#' sim_group("long", 5, "exponential", range = 0.5)
#' @export
sim_group <- function(name, n_vars, regime, ..., large_scale_share = 0.6, shared_with = NULL) {
  stopifnot(is.character(name), length(name) == 1, n_vars >= 1, n_vars == round(n_vars),
            is.character(regime), length(regime) == 1,
            large_scale_share >= 0, large_scale_share <= 1)
  if (regime == "none") large_scale_share <- 0
  structure(list(name = name, n_vars = as.integer(n_vars), regime = regime,
                 params = list(...), large_scale_share = large_scale_share,
                 shared_with = shared_with),
            class = "depprep_group")
}

#' Combine dependence groups
#'
#' @param ... Groups created with [sim_group()].
#' @return An object of class `depprep_groups`: a list of groups with the variable
#'   names and the group of each variable as attributes.
#' @examples
#' g <- sim_groups(sim_group("long", 3, "exponential", range = 0.5),
#'                 sim_group("none", 2, "none"))
#' attr(g, "var_group")
#' @export
sim_groups <- function(...) {
  groups <- list(...)
  if (!length(groups) || !all(vapply(groups, inherits, logical(1), "depprep_group"))) {
    stop("Pass one or more groups created with sim_group().", call. = FALSE)
  }
  names(groups) <- vapply(groups, `[[`, "", "name")
  if (anyDuplicated(names(groups))) stop("Group names must be unique.", call. = FALSE)
  for (g in groups) {
    if (!is.null(g$shared_with)) {
      src <- groups[[g$shared_with]]
      if (is.null(src)) stop("Group '", g$name, "' shares with unknown group '", g$shared_with, "'.", call. = FALSE)
      if (!is.null(src$shared_with)) stop("Group '", g$shared_with, "' already shares another group's field.", call. = FALSE)
      if (src$regime == "none") stop("Group '", g$name, "' cannot share the field of a group without dependence.", call. = FALSE)
    }
  }
  var_group <- rep(names(groups), vapply(groups, `[[`, 0L, "n_vars"))
  var_names <- unlist(lapply(groups, function(g) paste0(g$name, "_", seq_len(g$n_vars))), use.names = FALSE)
  structure(groups, class = "depprep_groups", var_group = var_group, var_names = var_names)
}

#' Correlation matrix of the residual components
#'
#' @param groups Groups created with [sim_groups()].
#' @param type `"identity"` (no association between residual components, the null
#'   structure), `"blocks"` (constant correlation within and between groups) or
#'   `"crossed"` (constant correlation within and between blocks that cross the
#'   groups: the variables of every group are dealt in turn to the `n_blocks`
#'   blocks, so each block holds variables of every group).
#' @param within Correlation between variables of the same group (`"blocks"`) or
#'   block (`"crossed"`).
#' @param between Correlation between variables of different groups (`"blocks"`) or
#'   blocks (`"crossed"`).
#' @param n_blocks Number of blocks for `"crossed"`.
#' @return A p x p positive definite correlation matrix with variable names.
#' @examples
#' g <- sim_groups(sim_group("a", 2, "none"), sim_group("b", 2, "none"))
#' make_R_true(g, "blocks", within = 0.5, between = 0.2)
#' make_R_true(g, "crossed", within = 0.5, between = 0)
#' @export
make_R_true <- function(groups, type = c("identity", "blocks", "crossed"), within = 0.5, between = 0.2,
                        n_blocks = 2) {
  stopifnot(inherits(groups, "depprep_groups"))
  type <- match.arg(type)
  vg <- attr(groups, "var_group")
  p <- length(vg)
  r <- diag(p)
  if (type != "identity") {
    b <- if (type == "blocks") vg else .crossed_blocks(vg, n_blocks)
    r <- ifelse(outer(b, b, "=="), within, between)
    diag(r) <- 1
  }
  if (min(eigen(r, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
    stop("This combination of `within` and `between` is not positive definite.", call. = FALSE)
  }
  dimnames(r) <- list(attr(groups, "var_names"), attr(groups, "var_names"))
  r
}

# Blocks that cross the groups: the i-th variable of every group goes to block
# ((i - 1) mod n_blocks) + 1.
.crossed_blocks <- function(var_group, n_blocks) {
  stopifnot(n_blocks >= 2, n_blocks == round(n_blocks))
  pos <- stats::ave(seq_along(var_group), var_group, FUN = seq_along)
  as.integer((pos - 1) %% n_blocks + 1)
}

.check_R <- function(R_true, groups) {
  p <- length(attr(groups, "var_group"))
  if (is.null(R_true)) R_true <- diag(p)
  R_true <- as.matrix(R_true)
  if (!all(dim(R_true) == p)) stop("`R_true` must be ", p, " x ", p, " (one row per variable).", call. = FALSE)
  if (max(abs(R_true - t(R_true))) > 1e-10 || max(abs(diag(R_true) - 1)) > 1e-10) {
    stop("`R_true` must be a symmetric correlation matrix.", call. = FALSE)
  }
  dimnames(R_true) <- list(attr(groups, "var_names"), attr(groups, "var_names"))
  R_true
}

.check_regime <- function(groups, allowed) {
  for (g in groups) {
    if (is.null(g$shared_with) && !g$regime %in% allowed) {
      stop("Regime '", g$regime, "' of group '", g$name, "' is not available here. Use one of: ",
           paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
  }
}

# Spread each group's large-scale component (one column per group) over its
# variables, following `shared_with`.
.large_by_variable <- function(groups, fields, n_obs) {
  cols <- lapply(groups, function(g) {
    src <- if (is.null(g$shared_with)) g$name else g$shared_with
    f <- fields[[src]]
    if (is.null(f)) f <- numeric(n_obs)
    matrix(f, n_obs, g$n_vars)
  })
  do.call(cbind, cols)
}

# Class labels and class effects. `spec` is NULL or a list with k, type
# ("random" or "contiguous"), effect and vars (indices of affected variables).
.make_classes <- function(spec, n_obs, p, seed, coords = NULL, n_times = NULL, expand = 1L) {
  if (is.null(spec)) return(list(labels = rep(1L, n_obs), effect = matrix(0, n_obs, p), spec = NULL))
  spec <- utils::modifyList(list(k = 3L, type = "random", effect = 1, vars = seq_len(p)), spec)
  stopifnot(spec$type %in% c("random", "contiguous"))
  .with_stream(seed, .streams[["classes"]], {
    # contiguous classes of sites are repeated at every time (`expand` times)
    labels <- if (spec$type == "contiguous") {
      rep(.contiguous_classes(spec$k, coords = coords, n_times = n_times), times = expand)
    } else {
      .random_classes(spec$k, n_obs)
    }
    mu <- .class_effects(spec$k, p, spec$effect, spec$vars)
    list(labels = as.integer(labels), effect = mu[labels, , drop = FALSE], spec = spec)
  })
}

# Assemble X and the simulation object.
.assemble <- function(family, groups, large, R_true, classes, seed, extra, call, counts = 0) {
  share <- rep(vapply(groups, `[[`, 0, "large_scale_share"), vapply(groups, `[[`, 0L, "n_vars"))
  n_obs <- nrow(large)
  residual <- .with_stream(seed, .streams[["residual"]], .residual_mvn(n_obs, R_true))
  x <- sweep(large, 2, sqrt(share), "*") + sweep(residual, 2, sqrt(1 - share), "*") + classes$effect
  vn <- attr(groups, "var_names")
  colnames(x) <- colnames(large) <- colnames(residual) <- colnames(classes$effect) <- vn
  count_vars <- .count_vars(ncol(x), counts)
  latent <- NULL
  if (length(count_vars)) {                               # perturbation: some variables observed as counts
    latent <- x
    x[, count_vars] <- .with_stream(seed, .streams[["counts"]], .to_counts(x[, count_vars, drop = FALSE]))
  }
  structure(c(list(
    family = family,
    X = x,
    group = stats::setNames(attr(groups, "var_group"), vn),
    components = list(large = large, residual = residual, class = classes$effect),
    large_scale_share = stats::setNames(share, vn),
    R_true = R_true,
    classes = classes$labels,
    class_spec = classes$spec,
    groups = groups,
    count_vars = vn[count_vars],
    X_latent = latent,
    seed = seed,
    call = call), extra), class = "depprep_sim")
}

#' @export
print.depprep_sim <- function(x, ...) {
  cat("<depprep simulation: ", x$family, ">\n", sep = "")
  cat(nrow(x$X), " observations x ", ncol(x$X), " variables; seed ", x$seed, "\n", sep = "")
  for (g in x$groups) {
    reg <- if (is.null(g$shared_with)) g$regime else paste0("shares field of '", g$shared_with, "'")
    par <- if (length(g$params)) paste(names(g$params), unlist(lapply(g$params, paste, collapse = "/")), sep = " = ", collapse = ", ") else ""
    cat(sprintf("  %-12s %2d vars  %-34s %s  (large-scale share %.2f)\n",
                g$name, g$n_vars, reg, par, g$large_scale_share))
  }
  invisible(x)
}

# Perturbation "counts" (protocol §4): the share `counts` of the variables, every k-th
# one (k = round(1 / counts), starting at the k-th), is observed as counts. Each chosen
# variable z (standardised) becomes Poisson with log-mean a + b z, b = 0.5 and
# a = log(5) - b^2 / 2, so the mean count is about 5 whatever the variable. The latent
# variables, the residual correlation R_true and the classes stay on the latent scale.
.count_vars <- function(p, counts) {
  if (is.null(counts) || counts <= 0) return(integer(0))
  if (counts > 1) stop("`counts` is a share of the variables, between 0 and 1.", call. = FALSE)
  k <- max(1, round(1 / counts))
  seq(k, p, by = k)
}
.to_counts <- function(Z, b = 0.5, mean_count = 5) {
  Z <- scale(as.matrix(Z))
  a <- log(mean_count) - b^2 / 2
  matrix(stats::rpois(length(Z), exp(a + b * Z)), nrow(Z), dimnames = dimnames(Z))
}

# Perturbation "non-stationary range" (protocol §4): a field whose practical range grows
# along x, from `ratio[1]` to `ratio[2]` times the nominal one. Two independent stationary
# fields with those ranges are mixed with weights sqrt(w) and sqrt(1 - w), w falling
# linearly from 1 at the smallest x to 0 at the largest: the variance stays 1 at every
# site and the local range moves smoothly from the short to the long one. The nominal
# range is the one recorded as truth (the middle of the domain is about halfway).
.grf_nonstationary <- function(coords, range, ratio = c(0.5, 1.5), anisotropy = 1, angle = 0, smoothness = 0.5) {
  f1 <- as.vector(.grf_chol(coords, range * ratio[1], anisotropy, angle, smoothness = smoothness))
  f2 <- as.vector(.grf_chol(coords, range * ratio[2], anisotropy, angle, smoothness = smoothness))
  x <- coords[, 1]
  w <- (max(x) - x) / (max(x) - min(x))
  sqrt(w) * f1 + sqrt(1 - w) * f2
}
