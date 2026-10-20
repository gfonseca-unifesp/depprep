#' Simulate a spatial data set with known truth
#'
#' Generates `n` sites in the unit square (or uses `coords`) and one time. Each
#' variable is the large-scale field of its group plus a residual component; the
#' residual components have correlation `R_true` and are independent across sites.
#'
#' Regimes for [sim_group()]:
#' * `"exponential"`: Gaussian random field with exponential covariance. Parameters
#'   `range` (practical range, in units of the side of the domain; correlation 0.05
#'   at this distance), `anisotropy` (ratio between the ranges along the major and
#'   minor axes, default 1) and `angle` (direction of the major axis in degrees,
#'   default 0).
#' * `"none"`: no large-scale component.
#'
#' @param n Number of sites (ignored when `coords` is given).
#' @param groups Groups created with [sim_groups()].
#' @param R_true Correlation of the residual components (p x p); identity if `NULL`.
#'   See [make_R_true()].
#' @param coords Optional two-column matrix of site coordinates.
#' @param classes Optional list describing observation classes: `k` (number of
#'   classes), `type` (`"random"` or `"contiguous"`), `effect` (standard deviation of
#'   the class means) and `vars` (indices of the variables that carry the classes).
#' @param seed Integer seed; recorded in the output.
#' @return An object of class `depprep_sim` with `X`, `coords`, `group`,
#'   `components` (`large`, `residual`, `class`), `R_true`, `classes`, `groups` and
#'   `seed`.
#' @examples
#' g <- sim_groups(sim_group("long", 3, "exponential", range = 0.5),
#'                 sim_group("short", 3, "exponential", range = 0.1),
#'                 sim_group("none", 2, "none"))
#' s <- sim_spatial(100, g, seed = 1)
#' s
#' @export
sim_spatial <- function(n, groups, R_true = NULL, coords = NULL, classes = NULL, seed) {
  stopifnot(inherits(groups, "depprep_groups"), !missing(seed))
  .check_regime(groups, c("exponential", "none"))
  R_true <- .check_R(R_true, groups)
  if (is.null(coords)) {
    coords <- .with_stream(seed, .streams[["coords"]], matrix(stats::runif(2 * n), ncol = 2))
  }
  coords <- as.matrix(coords)
  colnames(coords) <- c("x", "y")
  n <- nrow(coords)

  fields <- .with_stream(seed, .streams[["large"]], {
    out <- list()
    for (g in groups) {
      if (!is.null(g$shared_with) || g$regime == "none") next
      p <- utils::modifyList(list(anisotropy = 1, angle = 0), g$params)
      if (is.null(p$range)) stop("Group '", g$name, "' needs `range`.", call. = FALSE)
      out[[g$name]] <- as.vector(.grf_chol(coords, p$range, p$anisotropy, p$angle))
    }
    out
  })
  large <- .large_by_variable(groups, fields, n)
  cls <- .make_classes(classes, n, ncol(large), seed, coords = coords)
  .assemble("spatial", groups, large, R_true, cls, seed, list(coords = coords), match.call())
}
