#' Default thresholds for classify_regime()
#'
#' Values calibrated on the stage 3 simulations (60,000 variables with known regime):
#' the p-value threshold keeps false positives near 5%; the anisotropy and seasonal
#' thresholds are the ones with the best accuracy on the threshold curves.
#'
#' @return A list with:
#' * `p_dependence`: p-value below which a variable is declared dependent: the
#'   combined nearest-neighbour Moran and Moran-correlogram test in space
#'   (`dependence_p`; `moran_p` when absent), Ljung-Box in time;
#' * `aniso_ratio`: ratio between directional ranges above which a spatial variable is
#'   declared anisotropic;
#' * `seasonal_strength`: share of variance at the dominant seasonal period above which
#'   a series is declared seasonal;
#' * `trend_r2`: R-squared of a linear trend above which a series is declared trending;
#' * `shift`: shift indicator above which a spatio-temporal pattern is declared moving.
#' @examples
#' default_thresholds()
#' @export
default_thresholds <- function() {
  list(p_dependence = 0.05, aniso_ratio = 2.5, seasonal_strength = 0.2, trend_r2 = 0.4, shift = 0.08)
}

#' Classify the dependence regime of each variable
#'
#' Turns the signatures from [diagnose_variables()] into a regime type with simple,
#' explicit rules. The scale (range, phi) stays in the signature; the type says what
#' kind of dependence the variable has:
#' * spatial: `"none"`, `"isotropic"`, `"anisotropic"`;
#' * temporal: `"none"`, `"seasonal"`, `"trend"`, `"ar"`;
#' * spatio-temporal: `"none"`, `"temporal_common"` (no spatial structure),
#'   `"moving"`, `"separable_seasonal"`, `"separable_ar1"`, `"spatial_only"`.
#'
#' @param signature Output of [diagnose_variables()].
#' @param thresholds A list as returned by [default_thresholds()].
#' @return A character vector of regime types, named by variable.
#' @examples
#' g <- sim_groups(sim_group("long", 2, "exponential", range = 0.5), sim_group("none", 2, "none"))
#' classify_regime(diagnose_variables(sim_spatial(200, g, seed = 1)))
#' @export
classify_regime <- function(signature, thresholds = default_thresholds()) {
  stopifnot(inherits(signature, "data.frame"), "family" %in% names(signature))
  th <- utils::modifyList(default_thresholds(), thresholds)
  s <- signature
  ge <- function(x, t) !is.na(x) & x > t
  sp_p <- if ("dependence_p" %in% names(s)) ifelse(is.na(s$dependence_p), s$moran_p, s$dependence_p) else s$moran_p
  sp_dep <- !is.na(sp_p) & sp_p < th$p_dependence
  tm_dep <- !is.na(s$ljung_p) & s$ljung_p < th$p_dependence
  out <- ifelse(s$family == "spatial",
           ifelse(!sp_dep, "none", ifelse(ge(s$aniso_ratio, th$aniso_ratio), "anisotropic", "isotropic")),
         ifelse(s$family == "temporal",
           ifelse(!tm_dep, "none", ifelse(ge(s$period_strength, th$seasonal_strength), "seasonal",
                  ifelse(ge(s$trend_r2, th$trend_r2), "trend", "ar"))),
           ifelse(!sp_dep & !tm_dep, "none",
             ifelse(!sp_dep, "temporal_common",
               ifelse(!tm_dep, "spatial_only",
                 ifelse(ge(s$shift_indicator, th$shift), "moving",
                   ifelse(ge(s$period_strength, th$seasonal_strength), "separable_seasonal", "separable_ar1")))))))
  stats::setNames(out, s$variable)
}
