# Simulation grid of decision D4 (with the D10 revision), shared by the stage scripts
# from stage 3 on: the groups of each family, the scenarios, the truth of each variable
# and the simulation of one replicate. Sourced by analysis/03_*.R and 04_*.R.

spatial_groups <- function(short_range) sim_groups(
  sim_group("long", 8, "exponential", range = 0.5),
  sim_group("short", 8, "exponential", range = short_range),
  sim_group("aniso", 7, "exponential", range = 0.25, anisotropy = 4, angle = 30),
  sim_group("none", 7, "none"))
temporal_groups <- function(weak_phi) sim_groups(
  sim_group("strong", 6, "ar1", phi = 0.9),
  sim_group("weak", 6, "ar1", phi = weak_phi),
  sim_group("seasonal", 6, "seasonal", period = 12),
  sim_group("trend", 6, "trend_shift"),
  sim_group("none", 6, "none"))
st_groups <- function() sim_groups(
  sim_group("sep_long", 6, "separable_ar1", range = 0.5, phi = 0.9),
  sim_group("sep_short", 6, "separable_seasonal", range = 0.25, period = 12),   # 0.1 was not identifiable (D10)
  sim_group("moving", 6, "moving", range = 0.25, velocity = c(0.02, 0)),
  sim_group("common", 6, "temporal_common", phi = 0.9),
  sim_group("none", 6, "none"))

scenarios <- rbind(
  data.frame(family = "spatial", n = c(100, 400, 100, 400), n_times = NA,
             heterogeneity = c("high", "high", "low", "low")),
  data.frame(family = "temporal", n = NA, n_times = c(120, 480, 120, 480),
             heterogeneity = c("high", "high", "low", "low")),
  data.frame(family = "spatiotemporal", n = c(50, 100), n_times = c(60, 120), heterogeneity = "D4"))
scenarios$id <- seq_len(nrow(scenarios))
scenarios$label <- with(scenarios, ifelse(family == "spatial", sprintf("S n=%d %s", n, heterogeneity),
                    ifelse(family == "temporal", sprintf("T T=%d %s", n_times, heterogeneity),
                           sprintf("ST %dx%d", n, n_times))))

# truth of each variable: regime type (as classify_regime() names it) and scale
truth_table <- function(groups) {
  do.call(rbind, lapply(groups, function(g) {
    src <- if (is.null(g$shared_with)) g else groups[[g$shared_with]]
    p <- src$params
    type <- switch(src$regime,
      exponential = if (!is.null(p$anisotropy) && p$anisotropy > 1) "anisotropic" else "isotropic",
      ar1 = "ar", trend_shift = "trend", seasonal = "seasonal", src$regime)
    data.frame(group = g$name, true_type = type,
               true_range = if (is.null(p$range)) NA_real_ else p$range,
               true_phi = if (is.null(p$phi)) NA_real_ else p$phi,
               true_period = if (is.null(p$period)) NA_real_ else p$period,
               share = g$large_scale_share, stringsAsFactors = FALSE)[rep(1, g$n_vars), ]
  }))
}

simulate <- function(sc, seed) {
  switch(sc$family,
    spatial = sim_spatial(sc$n, spatial_groups(if (sc$heterogeneity == "high") 0.1 else 0.25), seed = seed),
    temporal = sim_temporal(sc$n_times, temporal_groups(if (sc$heterogeneity == "high") 0.3 else 0.6), seed = seed),
    spatiotemporal = sim_spatiotemporal(sc$n, sc$n_times, st_groups(), seed = seed))
}
