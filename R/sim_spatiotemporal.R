#' Simulate a spatio-temporal data set with known truth
#'
#' `n` sites in the unit square (or `coords`) observed at `n_times` regular times.
#' Rows of `X` are ordered time by time (all sites at time 1, then all sites at
#' time 2, ...); `index` gives the site and time of each row. The residual
#' components have correlation `R_true` and are independent across rows.
#'
#' Regimes for [sim_group()]:
#' * `"separable_ar1"`: exponential spatial field evolving as AR(1) in time,
#'   `Z_t = phi Z_(t-1) + sqrt(1 - phi^2) E_t`; parameters `range`, `phi`,
#'   `anisotropy`, `angle`.
#' * `"separable_seasonal"`: exponential spatial field times a harmonic (an amplitude
#'   map with a seasonal cycle); parameters `range`, `period`, `phase`.
#' * `"moving"`: non-separable pattern that moves with constant `velocity` (a
#'   two-element vector, in units of the side of the domain per time step); the
#'   field is drawn with random Fourier features of the exponential covariance;
#'   parameters `range`, `velocity`, `anisotropy`, `angle`, `n_features`.
#' * `"temporal_common"`: one AR(1) series shared by all sites; parameter `phi`.
#'
#' The spatial fields of the first three regimes are exponential by default; the
#' optional parameter `smoothness` (0.5, 1.5 or 2.5) makes them Matern, smoother for
#' larger values (see [sim_spatial()]).
#' * `"none"`: no large-scale component.
#'
#' @inheritParams sim_spatial
#' @param n_times Number of time steps.
#' @param classes Optional list as in [sim_spatial()]; contiguous classes are
#'   spatial regions that keep their class at every time.
#' @return An object of class `depprep_sim` with `X`, `coords`, `time`, `index`,
#'   `group`, `components`, `R_true`, `classes`, `groups` and `seed`.
#' @examples
#' g <- sim_groups(sim_group("sep", 3, "separable_ar1", range = 0.5, phi = 0.9),
#'                 sim_group("move", 3, "moving", range = 0.25, velocity = c(0.02, 0)),
#'                 sim_group("none", 2, "none"))
#' sim_spatiotemporal(50, 60, g, seed = 1)
#' @export
sim_spatiotemporal <- function(n, n_times, groups, R_true = NULL, coords = NULL, classes = NULL, seed) {
  stopifnot(inherits(groups, "depprep_groups"), !missing(seed), n_times >= 2)
  .check_regime(groups, c("separable_ar1", "separable_seasonal", "moving", "temporal_common", "none"))
  R_true <- .check_R(R_true, groups)
  if (is.null(coords)) {
    coords <- .with_stream(seed, .streams[["coords"]], matrix(stats::runif(2 * n), ncol = 2))
  }
  coords <- as.matrix(coords)
  colnames(coords) <- c("x", "y")
  n <- nrow(coords)
  need <- function(p, what, g) {
    miss <- setdiff(what, names(p))
    if (length(miss)) stop("Group '", g, "' needs ", paste0("`", miss, "`", collapse = ", "), ".", call. = FALSE)
  }

  fields <- .with_stream(seed, .streams[["large"]], {
    out <- list()
    for (g in groups) {
      if (!is.null(g$shared_with) || g$regime == "none") next
      p <- utils::modifyList(list(anisotropy = 1, angle = 0, smoothness = 0.5), g$params)
      z <- switch(g$regime,
        separable_ar1 = {
          need(p, c("range", "phi"), g$name)
          e <- .grf_chol(coords, p$range, p$anisotropy, p$angle, n_sim = n_times, smoothness = p$smoothness)
          for (t in seq_len(n_times)[-1]) e[, t] <- p$phi * e[, t - 1] + sqrt(1 - p$phi^2) * e[, t]
          e
        },
        separable_seasonal = {
          need(p, c("range", "period"), g$name)
          phase <- if (is.null(p$phase)) stats::runif(1, 0, 2 * pi) else p$phase
          f <- .grf_chol(coords, p$range, p$anisotropy, p$angle, smoothness = p$smoothness)
          as.vector(f) %o% .seasonal(n_times, p$period, phase)
        },
        moving = {
          need(p, c("range", "velocity"), g$name)
          f <- .rff_sampler(p$range, p$anisotropy, p$angle,
                            n_features = if (is.null(p$n_features)) 1000 else p$n_features, smoothness = p$smoothness)
          vapply(seq_len(n_times), function(t) {
            f(coords - matrix(p$velocity * t, n, 2, byrow = TRUE))
          }, numeric(n))
        },
        temporal_common = {
          need(p, "phi", g$name)
          matrix(rep(.ar1(n_times, p$phi), each = n), n, n_times)
        })
      out[[g$name]] <- as.vector(z)  # time-major: sites vary fastest
    }
    out
  })
  n_obs <- n * n_times
  large <- .large_by_variable(groups, fields, n_obs)
  cls <- .make_classes(classes, n_obs, ncol(large), seed, coords = coords, expand = n_times)
  index <- data.frame(site = rep(seq_len(n), times = n_times), time = rep(seq_len(n_times), each = n))
  .assemble("spatiotemporal", groups, large, R_true, cls, seed,
            list(coords = coords, time = seq_len(n_times), index = index), match.call())
}
