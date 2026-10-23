# Simulation grid of decision D4 (with the D10 revision), shared by the stage scripts
# from stage 3 on: the groups of each family, the scenarios, the truth of each variable
# and the simulation of one replicate. Sourced by analysis/03_*.R, 04_*.R and 05_*.R.
# `scenarios` is the D4 grid (ids 1-10); `extra_scenarios` (ids 11 on, stage 5) adds
# smooth large-scale fields (decision D12) and the "same regime, different fields" case
# for the second grouping level (D11). Smoothness 0.5 reproduces the D4 data exactly.
# `d14_scenarios` (ids 19 on, decision D14) adds a spatial gradient, a seasonal cycle
# with a trend and a period off the Fourier grid; `pilot_scenarios` gains a truth with
# blocks that cross the regime groups and one with correlation between groups. The
# columns `variant` and `between` default to "base" and 0, which give the earlier data.

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
# D14: every D4 regime of the family plus the new one, 6 variables per group
gradient_groups <- function() sim_groups(
  sim_group("long", 6, "exponential", range = 0.5),
  sim_group("short", 6, "exponential", range = 0.1),
  sim_group("aniso", 6, "exponential", range = 0.25, anisotropy = 4, angle = 30),
  sim_group("gradient", 6, "gradient", direction = 60, strength = 1, range = 0.25),   # plane plus field
  sim_group("none", 6, "none"))
seasonal_trend_groups <- function() sim_groups(
  sim_group("strong", 5, "ar1", phi = 0.9),
  sim_group("weak", 5, "ar1", phi = 0.3),
  sim_group("seasonal", 5, "seasonal", period = 12),
  sim_group("trend", 5, "trend_shift"),
  sim_group("seas_trend", 5, "seasonal_trend", period = 12, strength = 1),
  sim_group("none", 5, "none"))
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
scenarios$variant <- "base"
scenarios$between <- 0

extra_scenarios <- rbind(
  transform(scenarios[scenarios$family != "temporal", ], smoothness = 2.5),
  data.frame(family = "fields", n = 400, n_times = NA, heterogeneity = "same range", id = NA, label = NA,
             smoothness = c(0.5, 2.5), variant = "base", between = 0))
extra_scenarios$id <- 10L + seq_len(nrow(extra_scenarios))
extra_scenarios$label <- ifelse(extra_scenarios$family == "fields",
                                sprintf("F n=400 %s", ifelse(extra_scenarios$smoothness == 0.5, "rough", "smooth")),
                                paste(extra_scenarios$label, "smooth"))

# D14 (2026-10-22): a gradient plus a field in space; a seasonal cycle plus a trend;
# the D4 temporal groups with T = 100, where the period 12 falls between Fourier
# frequencies (100 / 12 is not an integer).
d14_scenarios <- data.frame(family = c("spatial", "temporal", "temporal"), n = c(400, NA, NA),
                            n_times = c(NA, 480, 100), heterogeneity = "high", id = NA,
                            label = c("S n=400 high gradient", "T T=480 high seasonal+trend", "T T=100 high"),
                            smoothness = 0.5, variant = c("gradient", "seasonal_trend", "base"), between = 0)
d14_scenarios$id <- max(extra_scenarios$id) + seq_len(nrow(d14_scenarios))

# truth of each variable: regime type (as classify_regime() names it), scale and
# whether the large-scale component holds a trend
truth_table <- function(groups) {
  do.call(rbind, lapply(groups, function(g) {
    src <- if (is.null(g$shared_with)) g else groups[[g$shared_with]]
    p <- src$params
    iso <- if (!is.null(p$anisotropy) && p$anisotropy > 1) "anisotropic" else "isotropic"
    type <- switch(src$regime,
      exponential = , matern = iso, gradient = if (is.null(p$range)) "trend" else iso,
      ar1 = "ar", trend_shift = "trend", seasonal = , seasonal_trend = "seasonal", src$regime)
    data.frame(group = g$name, true_type = type,
               true_range = if (is.null(p$range)) NA_real_ else p$range,
               true_phi = if (is.null(p$phi)) NA_real_ else p$phi,
               true_period = if (is.null(p$period)) NA_real_ else p$period,
               true_trend = src$regime %in% c("trend_shift", "gradient", "seasonal_trend"),
               share = g$large_scale_share, stringsAsFactors = FALSE)[rep(1, g$n_vars), ]
  }))
}

# `truth` (stage 6 pilot): "null" (default: identity R_true, no classes, the D4 data),
# "structured" (R_true in blocks, 0.5 within groups and `between` across them; 3 random
# classes), "contiguous" (the same with 3 classes contiguous in space or time) or
# "crossed" (D14: 0.5 within two blocks that cross the regime groups, each holding
# every other variable of every group, and `between` across them; 3 random classes).
# Stage 7 perturbations (protocol-v1, section 4): `counts` turns that share of the
# variables into counts; `nonstationary` (e.g. c(0.5, 1.5)) makes the range of every
# spatial field grow along x from 0.5 to 1.5 times its nominal value. Both default to none,
# which gives the earlier data.
simulate <- function(sc, seed, truth = "null", between = if (is.null(sc$between)) 0 else sc$between,
                     counts = 0, nonstationary = NULL) {
  sm <- if (is.null(sc$smoothness)) 0.5 else sc$smoothness
  variant <- if (is.null(sc$variant)) "base" else sc$variant
  g <- switch(variant,
    gradient = gradient_groups(),
    seasonal_trend = seasonal_trend_groups(),
    switch(sc$family,
      spatial = spatial_groups(if (sc$heterogeneity == "high") 0.1 else 0.25, sm),
      fields = fields_groups(sm),
      temporal = temporal_groups(if (sc$heterogeneity == "high") 0.3 else 0.6),
      spatiotemporal = st_groups(sm)))
  if (!is.null(nonstationary)) {
    attrs <- attributes(g)
    g <- lapply(g, function(x) { if (x$regime %in% c("exponential", "matern") && is.null(x$shared_with))
      x$params$nonstationary <- nonstationary; x })
    attributes(g) <- attrs
  }
  R_true <- switch(truth, null = NULL,
                   crossed = make_R_true(g, "crossed", within = 0.5, between = between, n_blocks = 2),
                   make_R_true(g, "blocks", within = 0.5, between = between))
  classes <- switch(truth, null = NULL, structured = , crossed = list(k = 3, type = "random", effect = 1),
                    contiguous = list(k = 3, type = "contiguous", effect = 1))
  switch(sc$family,
    spatial = , fields = sim_spatial(sc$n, g, R_true = R_true, classes = classes, seed = seed, counts = counts),
    temporal = sim_temporal(sc$n_times, g, R_true = R_true, classes = classes, seed = seed, counts = counts),
    spatiotemporal = sim_spatiotemporal(sc$n, sc$n_times, g, R_true = R_true, classes = classes, seed = seed,
                                        counts = counts))
}

# Stage 6 pilot: the largest affordable size of each family, in three truths.
pilot_scenarios <- merge(rbind(scenarios[scenarios$label %in% c("S n=400 high", "T T=480 high", "ST 50x60"), ],
                               extra_scenarios[extra_scenarios$label == "S n=400 high smooth", ]),
                         data.frame(truth = c("null", "structured", "contiguous"), stringsAsFactors = FALSE))
pilot_scenarios <- pilot_scenarios[order(match(pilot_scenarios$family, c("spatial", "temporal", "spatiotemporal")),
                                         pilot_scenarios$smoothness, match(pilot_scenarios$truth, c("null", "structured", "contiguous"))), ]
pilot_scenarios$base_id <- pilot_scenarios$id
pilot_scenarios$id <- seq_len(nrow(pilot_scenarios))
pilot_scenarios$label <- paste(pilot_scenarios$label, pilot_scenarios$truth)

# D14 (2026-10-22): the same four scenarios with blocks that cross the regime groups
# (the pilot's own blocks coincide with the groups, which favours raw data in PCA) and
# with a between-group correlation of 0.3; ids 13 on.
d14_truths <- data.frame(truth = c("crossed", "structured"), between = c(0, 0.3), stringsAsFactors = FALSE)
d14_base <- pilot_scenarios[pilot_scenarios$truth == "null", setdiff(names(pilot_scenarios), c("truth", "between"))]
d14_pilot <- merge(d14_base, d14_truths)
d14_pilot$label <- paste(sub(" null$", "", d14_pilot$label), d14_pilot$truth,
                         ifelse(d14_pilot$between > 0, paste("between", d14_pilot$between), ""))
d14_pilot$label <- trimws(d14_pilot$label)
d14_pilot <- d14_pilot[order(match(d14_pilot$family, c("spatial", "temporal", "spatiotemporal")),
                             d14_pilot$smoothness, match(d14_pilot$truth, d14_truths$truth)), ]
d14_pilot$id <- max(pilot_scenarios$id) + seq_len(nrow(d14_pilot))
pilot_scenarios <- rbind(pilot_scenarios, d14_pilot[names(pilot_scenarios)])
