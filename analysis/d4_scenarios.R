# Simulation grid of decision D4 (with the D10 revision), shared by the stage scripts
# from stage 3 on: the groups of each family, the scenarios, the truth of each variable
# and the simulation of one replicate. Sourced by analysis/03_*.R, 04_*.R and 05_*.R.
# `scenarios` is the D4 grid (ids 1-10); `extra_scenarios` (ids 11 on, stage 5) adds
# smooth large-scale fields (decision D12) and the "same regime, different fields" case
# for the second grouping level (D11). Smoothness 0.5 reproduces the D4 data exactly.

spatial_field <- function(name, n_vars, smoothness, ...) {
  if (smoothness == 0.5) sim_group(name, n_vars, "exponential", ...) else
    sim_group(name, n_vars, "matern", smoothness = smoothness, ...)
}
spatial_groups <- function(short_range, smoothness = 0.5) sim_groups(
  spatial_field("long", 8, smoothness, range = 0.5),
  spatial_field("short", 8, smoothness, range = short_range),
  spatial_field("aniso", 7, smoothness, range = 0.25, anisotropy = 4, angle = 30),
  sim_group("none", 7, "none"))
fields_groups <- function(smoothness) sim_groups(
  spatial_field("long_a", 6, smoothness, range = 0.5),
  spatial_field("long_b", 6, smoothness, range = 0.5),          # same regime, independent field
  spatial_field("short", 6, smoothness, range = 0.15),
  sim_group("none", 6, "none"))
temporal_groups <- function(weak_phi) sim_groups(
  sim_group("strong", 6, "ar1", phi = 0.9),
  sim_group("weak", 6, "ar1", phi = weak_phi),
  sim_group("seasonal", 6, "seasonal", period = 12),
  sim_group("trend", 6, "trend_shift"),
  sim_group("none", 6, "none"))
st_groups <- function(smoothness = 0.5) {
  sm <- if (smoothness == 0.5) list() else list(smoothness = smoothness)
  sim_groups(
  do.call(sim_group, c(list("sep_long", 6, "separable_ar1", range = 0.5, phi = 0.9), sm)),
  do.call(sim_group, c(list("sep_short", 6, "separable_seasonal", range = 0.25, period = 12), sm)),   # 0.1 was not identifiable (D10)
  do.call(sim_group, c(list("moving", 6, "moving", range = 0.25, velocity = c(0.02, 0)), sm)),
  sim_group("common", 6, "temporal_common", phi = 0.9),
  sim_group("none", 6, "none"))
}

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
scenarios$smoothness <- 0.5

extra_scenarios <- rbind(
  transform(scenarios[scenarios$family != "temporal", ], smoothness = 2.5),
  data.frame(family = "fields", n = 400, n_times = NA, heterogeneity = "same range", id = NA, label = NA,
             smoothness = c(0.5, 2.5)))
extra_scenarios$id <- 10L + seq_len(nrow(extra_scenarios))
extra_scenarios$label <- ifelse(extra_scenarios$family == "fields",
                                sprintf("F n=400 %s", ifelse(extra_scenarios$smoothness == 0.5, "rough", "smooth")),
                                paste(extra_scenarios$label, "smooth"))

# truth of each variable: regime type (as classify_regime() names it) and scale
truth_table <- function(groups) {
  do.call(rbind, lapply(groups, function(g) {
    src <- if (is.null(g$shared_with)) g else groups[[g$shared_with]]
    p <- src$params
    type <- switch(src$regime,
      exponential = , matern = if (!is.null(p$anisotropy) && p$anisotropy > 1) "anisotropic" else "isotropic",
      ar1 = "ar", trend_shift = "trend", seasonal = "seasonal", src$regime)
    data.frame(group = g$name, true_type = type,
               true_range = if (is.null(p$range)) NA_real_ else p$range,
               true_phi = if (is.null(p$phi)) NA_real_ else p$phi,
               true_period = if (is.null(p$period)) NA_real_ else p$period,
               share = g$large_scale_share, stringsAsFactors = FALSE)[rep(1, g$n_vars), ]
  }))
}

simulate <- function(sc, seed) {
  sm <- if (is.null(sc$smoothness)) 0.5 else sc$smoothness
  switch(sc$family,
    spatial = sim_spatial(sc$n, spatial_groups(if (sc$heterogeneity == "high") 0.1 else 0.25, sm), seed = seed),
    fields = sim_spatial(sc$n, fields_groups(sm), seed = seed),
    temporal = sim_temporal(sc$n_times, temporal_groups(if (sc$heterogeneity == "high") 0.3 else 0.6), seed = seed),
    spatiotemporal = sim_spatiotemporal(sc$n, sc$n_times, st_groups(sm), seed = seed))
}
