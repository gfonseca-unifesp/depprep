#' Simulate a temporal data set with known truth
#'
#' One site observed at `n_times` regular times (the `T` of the protocol). Each
#' variable is the large-scale component of its group plus a residual component; the
#' residual components have correlation `R_true` and are independent across times.
#'
#' Regimes for [sim_group()]:
#' * `"ar1"`: stationary AR(1) with unit variance; parameter `phi`.
#' * `"seasonal"`: harmonic with unit variance; parameters `period` and `phase`
#'   (random when omitted).
#' * `"trend_shift"`: linear trend with a jump, standardised; parameters `slope`
#'   (default 1), `jump` (default 1) and `shift_at` (fraction of the series where
#'   the jump happens, default 0.5).
#' * `"seasonal_trend"`: harmonic plus linear trend, standardised; parameters
#'   `period`, `phase` (random when omitted) and `strength` (standard deviation of
#'   the trend relative to that of the harmonic, default 1).
#' * `"none"`: no large-scale component.
#'
#' @inheritParams sim_spatial
#' @param n_times Number of time steps.
#' @param classes Optional list as in [sim_spatial()]; contiguous classes are blocks
#'   of consecutive times.
#' @return An object of class `depprep_sim` with `X`, `time`, `group`, `components`,
#'   `R_true`, `classes`, `groups` and `seed`.
#' @examples
#' g <- sim_groups(sim_group("strong", 3, "ar1", phi = 0.9),
#'                 sim_group("season", 3, "seasonal", period = 12),
#'                 sim_group("none", 2, "none"))
#' sim_temporal(120, g, seed = 1)
#' @export
sim_temporal <- function(n_times, groups, R_true = NULL, classes = NULL, seed, counts = 0) {
  stopifnot(inherits(groups, "depprep_groups"), !missing(seed), n_times >= 2)
  .check_regime(groups, c("ar1", "seasonal", "trend_shift", "seasonal_trend", "none"))
  R_true <- .check_R(R_true, groups)

  fields <- .with_stream(seed, .streams[["large"]], {
    out <- list()
    for (g in groups) {
      if (!is.null(g$shared_with) || g$regime == "none") next
      p <- g$params
      out[[g$name]] <- switch(g$regime,
        ar1 = {
          if (is.null(p$phi)) stop("Group '", g$name, "' needs `phi`.", call. = FALSE)
          as.vector(.ar1(n_times, p$phi))
        },
        seasonal = {
          if (is.null(p$period)) stop("Group '", g$name, "' needs `period`.", call. = FALSE)
          phase <- if (is.null(p$phase)) stats::runif(1, 0, 2 * pi) else p$phase
          .seasonal(n_times, p$period, phase)
        },
        trend_shift = {
          p <- utils::modifyList(list(slope = 1, jump = 1, shift_at = 0.5), p)
          .trend_shift(n_times, p$slope, p$jump, p$shift_at)
        },
        seasonal_trend = {
          if (is.null(p$period)) stop("Group '", g$name, "' needs `period`.", call. = FALSE)
          phase <- if (is.null(p$phase)) stats::runif(1, 0, 2 * pi) else p$phase
          .seasonal_trend(n_times, p$period, phase, if (is.null(p$strength)) 1 else p$strength)
        })
    }
    out
  })
  large <- .large_by_variable(groups, fields, n_times)
  cls <- .make_classes(classes, n_times, ncol(large), seed, n_times = n_times)
  .assemble("temporal", groups, large, R_true, cls, seed, list(time = seq_len(n_times)), match.call(), counts = counts)
}
