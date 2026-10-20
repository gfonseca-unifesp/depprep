# Stage 2 — Examples of the three simulation families, to support decision D4
# (final grid values). Every figure uses a fixed seed, recorded in
# results/sim_examples/seeds.csv. Candidate values shown here are not decisions.
#
# Outputs (results/figs/):
#   sim_scheme_{spatial,temporal,spatiotemporal}.png  generator schemes
#   sim_spatial_gallery.png, sim_spatial_variograms.png
#   sim_temporal_gallery.png, sim_temporal_acf.png
#   sim_st_maps.png, sim_st_hovmoller.png
#   sim_apparent_association.png
#   sim_classes.png
#   sim_cost.png
# and results/sim_examples/{seeds,cost,association}.csv

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})

fig_dir <- file.path("results", "figs")
out_dir <- file.path("results", "sim_examples")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidates <- list(            # candidate values for D4 (protocol section 4)
  range = c(0.05, 0.1, 0.25, 0.5),
  anisotropy = c(1, 2, 4),
  phi = c(0.3, 0.6, 0.9),
  period = 12,
  velocity = c(0.02, 0),
  large_scale_share = 0.6
)
seeds <- c(gallery = 2101, variogram = 2102, temporal = 2103, st = 2104,
           association = 2105, classes = 2106, cost = 2107)
write.csv(data.frame(figure = names(seeds), seed = seeds), file.path(out_dir, "seeds.csv"), row.names = FALSE)
grid_xy <- function(m) as.matrix(expand.grid(x = (seq_len(m) - 0.5) / m, y = (seq_len(m) - 0.5) / m))
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
theme_set(theme_minimal(base_size = 10))

# --- 1. Generator schemes -------------------------------------------------------------
scheme <- function(title, groups_txt, shared_from = NULL) {
  ng <- length(groups_txt)
  boxes <- data.frame(id = paste0("g", seq_len(ng)), x = 0, y = rev(seq_len(ng)) * 1.2, label = groups_txt)
  extra <- data.frame(id = c("res", "cls", "x"), x = c(0, 0, 3.2),
                      y = c(-0.4, -1.6, mean(range(c(boxes$y, -1.6)))),
                      label = c("Residual e ~ N(0, R_true)\nindependent across observations",
                                "Classes (optional)\nmean shift per class, random or contiguous",
                                "x_j = sqrt(s_g) L_g\n + sqrt(1 - s_g) e_j\n + class effect"))
  all <- rbind(boxes, extra)
  arrows <- data.frame(x = 1.55, y = all$y[all$id != "x"], xend = 2.35, yend = all$y[all$id == "x"])
  p <- ggplot() +
    geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), colour = "grey55",
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_label(data = all, aes(x, y, label = label), size = 3, linewidth = 0.3, hjust = 0.5,
               fill = ifelse(all$id == "x", "#e0eff0", ifelse(all$id %in% c("res", "cls"), "#f8ecd9", "white"))) +
    labs(title = title, subtitle = "Large-scale component of each group (left) + residual + classes = variable") +
    xlim(-1.6, 4.6) + theme_void(base_size = 10)
  if (!is.null(shared_from)) {
    p <- p + geom_curve(data = data.frame(x = -1.55, y = boxes$y[shared_from[1]], xend = -1.55, yend = boxes$y[shared_from[2]]),
                        aes(x, y, xend = xend, yend = yend), curvature = 0.4, linetype = "dashed", colour = "#9a5806",
                        arrow = arrow(length = unit(2, "mm"))) +
      annotate("text", x = -1.5, y = mean(boxes$y[shared_from]), label = "same field", colour = "#9a5806", size = 3, hjust = 0)
  }
  p
}
save_fig(scheme("Spatial family (S): n sites, one time",
                c("long_a: exponential field, long range",
                  "long_b: shares the field of long_a",
                  "short: exponential field, short range",
                  "aniso: exponential field, anisotropic",
                  "none: no large-scale component"), shared_from = c(1, 2)),
         "sim_scheme_spatial.png", 8, 4.5)
save_fig(scheme("Temporal family (T): one site, T times",
                c("strong: AR(1), phi high", "weak: AR(1), phi low", "seasonal: harmonic, period 12",
                  "trend: linear trend with a jump", "none: no large-scale component")),
         "sim_scheme_temporal.png", 8, 4.5)
save_fig(scheme("Spatio-temporal family (ST): n sites x T times",
                c("sep_long: long-range field x AR(1)", "sep_short: short-range map x seasonal cycle",
                  "moving: pattern displaced each step (non-separable)", "common: one series at all sites",
                  "none: no large-scale component")),
         "sim_scheme_spatiotemporal.png", 8, 4.5)

# --- 2. Spatial gallery and variograms ------------------------------------------------
xy <- grid_xy(40)
gal <- do.call(rbind, c(
  lapply(candidates$range, function(r) {
    s <- sim_spatial(groups = sim_groups(sim_group("f", 1, "exponential", range = r, large_scale_share = 1)),
                     coords = xy, seed = seeds[["gallery"]])
    data.frame(xy, z = s$components$large[, 1], panel = sprintf("range %.2f", r), row = "Range (isotropic)")
  }),
  lapply(candidates$anisotropy, function(a) {
    s <- sim_spatial(groups = sim_groups(sim_group("f", 1, "exponential", range = 0.25, anisotropy = a,
                                                   angle = 30, large_scale_share = 1)),
                     coords = xy, seed = seeds[["gallery"]])
    data.frame(xy, z = s$components$large[, 1], panel = sprintf("anisotropy %g (range 0.25, 30 deg)", a), row = "Anisotropy")
  }),
  list(data.frame(xy, z = sim_spatial(groups = sim_groups(sim_group("f", 1, "none")), coords = xy,
                                      seed = seeds[["gallery"]])$X[, 1], panel = "no dependence", row = "Anisotropy"))))
gal$panel <- factor(gal$panel, unique(gal$panel))
p <- ggplot(gal, aes(x, y, fill = z)) + geom_raster() + coord_equal() +
  facet_wrap(~ panel, nrow = 2) + scale_fill_viridis_c(name = "value") +
  labs(title = "Spatial regimes: candidate ranges and anisotropy (unit square)",
       subtitle = "Same random numbers in every panel; only the parameter changes", x = NULL, y = NULL) +
  theme(axis.text = element_blank())
save_fig(p, "sim_spatial_gallery.png", 10, 5.5)

vg <- do.call(rbind, lapply(candidates$range, function(r) {
  g <- sim_groups(sim_group("f", 1, "exponential", range = r, large_scale_share = 1))
  do.call(rbind, lapply(1:10, function(i) {
    s <- sim_spatial(400, g, seed = seeds[["variogram"]] * 100 + i)
    v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(s$coords, z = s$components$large[, 1]),
                          cutoff = 0.7, width = 0.035)
    data.frame(dist = v$dist, gamma = v$gamma, rep = i, panel = sprintf("range %.2f", r), range = r)
  }))
}))
truth <- do.call(rbind, lapply(candidates$range, function(r) {
  h <- seq(0, 0.7, length.out = 100)
  data.frame(dist = h, gamma = 1 - exp(-3 * h / r), panel = sprintf("range %.2f", r))
}))
p <- ggplot(vg, aes(dist, gamma)) +
  geom_line(aes(group = rep), colour = "grey60", linewidth = 0.3) +
  geom_line(data = truth, colour = "#0d6b72", linewidth = 1) +
  geom_vline(data = unique(vg[, c("panel", "range")]), aes(xintercept = range), linetype = "dashed", colour = "#9a5806") +
  facet_wrap(~ panel, nrow = 1) +
  labs(title = "Empirical variograms (10 replicates, n = 400) against the true curve",
       subtitle = "Dashed line: simulated practical range. Long ranges reach the size of the domain and are hard to estimate",
       x = "distance (side of the domain = 1)", y = "semivariance")
save_fig(p, "sim_spatial_variograms.png", 10, 3.4)

# --- 3. Temporal gallery and ACF --------------------------------------------------------
tg <- sim_groups(sim_group("phi 0.3", 1, "ar1", phi = 0.3), sim_group("phi 0.6", 1, "ar1", phi = 0.6),
                 sim_group("phi 0.9", 1, "ar1", phi = 0.9), sim_group("seasonal", 1, "seasonal", period = 12),
                 sim_group("trend + jump", 1, "trend_shift"), sim_group("none", 1, "none"))
st <- sim_temporal(120, tg, seed = seeds[["temporal"]])
ts_df <- do.call(rbind, lapply(seq_len(ncol(st$X)), function(j) {
  data.frame(time = st$time, large = unname(st$components$large[, j]), x = unname(st$X[, j]), panel = unname(st$group[j]))
}))
ts_df$panel <- factor(ts_df$panel, unique(ts_df$panel))
p <- ggplot(ts_df, aes(time)) +
  geom_line(aes(y = x), colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = large), colour = "#0d6b72", linewidth = 0.8) +
  facet_wrap(~ panel, ncol = 2) +
  labs(title = "Temporal regimes (T = 120 months)",
       subtitle = "Blue: large-scale component. Grey: observed variable (large scale + residual, share 0.6)",
       x = "time step", y = NULL)
save_fig(p, "sim_temporal_gallery.png", 9, 6)
acf_df <- do.call(rbind, lapply(seq_len(ncol(st$X)), function(j) {
  a <- stats::acf(st$X[, j], lag.max = 24, plot = FALSE)$acf[-1]
  data.frame(lag = 1:24, acf = a, panel = unname(st$group[j]))
}))
acf_df$panel <- factor(acf_df$panel, levels(ts_df$panel))
p <- ggplot(acf_df, aes(lag, acf)) + geom_col(fill = "#0d6b72", width = 0.7) +
  geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(120), linetype = "dashed", colour = "#9a5806") +
  facet_wrap(~ panel, ncol = 3) +
  labs(title = "Autocorrelation function of the observed variables", subtitle = "Dashed: 95% band for white noise",
       x = "lag (time steps)", y = "ACF")
save_fig(p, "sim_temporal_acf.png", 9, 4.5)

# --- 4. Spatio-temporal maps and Hovmoller --------------------------------------------
xy <- grid_xy(30)
stg <- sim_groups(
  sim_group("separable: long range x AR(1) 0.9", 1, "separable_ar1", range = 0.5, phi = 0.9),
  sim_group("separable: short range x seasonal", 1, "separable_seasonal", range = 0.1, period = 12),
  sim_group("non-separable: moving pattern", 1, "moving", range = 0.25, velocity = candidates$velocity),
  sim_group("temporal only, common to all sites", 1, "temporal_common", phi = 0.9))
sst <- sim_spatiotemporal(groups = stg, n_times = 24, coords = xy, seed = seeds[["st"]])
maps <- do.call(rbind, lapply(seq_len(ncol(sst$X)), function(j) {
  d <- data.frame(sst$index, z = unname(sst$components$large[, j]), regime = unname(sst$group[j]))
  d <- cbind(d, xy[d$site, ])
  d[d$time %in% c(1, 4, 7, 10), ]
}))
maps$regime <- factor(maps$regime, names(stg))
p <- ggplot(maps, aes(x, y, fill = z)) + geom_raster() + coord_equal() +
  facet_grid(regime ~ factor(paste("t =", time), paste("t =", sort(unique(time))))) +
  scale_fill_viridis_c(name = "value") +
  labs(title = "Spatio-temporal regimes: the same grid at four times", x = NULL, y = NULL) +
  theme(axis.text = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0))
save_fig(p, "sim_st_maps.png", 10, 7)
hov <- do.call(rbind, lapply(seq_len(ncol(sst$X)), function(j) {
  d <- data.frame(sst$index, z = unname(sst$components$large[, j]), regime = unname(sst$group[j]))
  d <- cbind(d, xy[d$site, ])
  d[abs(d$y - xy[which.min(abs(xy[, 2] - 0.5)), 2]) < 1e-9, ]
}))
hov$regime <- factor(hov$regime, names(stg))
p <- ggplot(hov, aes(x, time, fill = z)) + geom_raster() + facet_wrap(~ regime, nrow = 1) +
  scale_fill_viridis_c(name = "value") +
  labs(title = "Hovmoller diagrams: one transect (y = 0.5) through time",
       subtitle = "Tilted stripes = a pattern that moves (non-separable). Vertical stripes = a pattern that stays in place",
       x = "position along the transect", y = "time step")
save_fig(p, "sim_st_hovmoller.png", 10, 3.6)

# --- 5. Apparent association ------------------------------------------------------------
ag <- sim_groups(sim_group("long_a", 6, "exponential", range = 0.5),
                 sim_group("long_b", 6, "exponential", shared_with = "long_a"),
                 sim_group("short", 6, "exponential", range = 0.1),
                 sim_group("aniso", 6, "exponential", range = 0.25, anisotropy = 4, angle = 30),
                 sim_group("none", 6, "none"))
assoc <- list()
for (structure_type in c("identity", "blocks")) {
  r <- make_R_true(ag, structure_type, within = 0.5, between = 0)
  s <- sim_spatial(400, ag, R_true = r, seed = seeds[["association"]])
  mats <- list("R_true (imposed)" = r, "Raw correlation of X" = stats::cor(s$X),
               "Correlation of the true residuals" = stats::cor(s$components$residual))
  for (m in names(mats)) {
    v <- mats[[m]]
    assoc[[length(assoc) + 1]] <- data.frame(expand.grid(i = seq_len(30), j = seq_len(30)), r = as.vector(v),
                                             matrix = m, truth = ifelse(structure_type == "identity", "null structure", "block structure"))
  }
}
assoc <- do.call(rbind, assoc)
assoc$matrix <- factor(assoc$matrix, c("R_true (imposed)", "Raw correlation of X", "Correlation of the true residuals"))
write.csv(aggregate(abs(r) ~ matrix + truth, data = assoc[assoc$i != assoc$j, ], FUN = mean),
          file.path(out_dir, "association.csv"), row.names = FALSE)
cuts <- cumsum(table(attr(ag, "var_group"))[names(ag)]) + 0.5
p <- ggplot(assoc, aes(i, j, fill = r)) + geom_raster() +
  geom_hline(yintercept = cuts[-length(cuts)], colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = cuts[-length(cuts)], colour = "white", linewidth = 0.3) +
  facet_grid(truth ~ matrix) + coord_equal() + scale_y_reverse() +
  scale_fill_gradient2(low = "#9a5806", mid = "white", high = "#0d6b72", limits = c(-1, 1), name = "correlation") +
  labs(title = "Apparent association: what the large-scale fields add to the correlation between variables",
       subtitle = "Groups (lines): long_a, long_b (shares long_a's field), short, aniso, none. n = 400, large-scale share 0.6",
       x = "variable", y = "variable")
save_fig(p, "sim_apparent_association.png", 10, 6.5)

# --- 6. Classes ---------------------------------------------------------------------------
cg <- sim_groups(sim_group("f", 2, "exponential", range = 0.25))
cls <- do.call(rbind, lapply(c("random", "contiguous"), function(ty) {
  s <- sim_spatial(400, cg, classes = list(k = 3, type = ty, effect = 1), seed = seeds[["classes"]])
  data.frame(s$coords, class = factor(s$classes), type = paste(ty, "classes"))
}))
p <- ggplot(cls, aes(x, y, colour = class)) + geom_point(size = 1.4) + coord_equal() +
  facet_wrap(~ type) + scale_colour_manual(values = c("#0d6b72", "#9a5806", "#586972")) +
  labs(title = "Observation classes (k = 3) for k-means and SOM",
       subtitle = "Contiguous classes coincide with large-scale structure: correcting it can erase them (Q3)",
       x = NULL, y = NULL)
save_fig(p, "sim_classes.png", 8, 4.2)

# --- 7. Generation cost -------------------------------------------------------------------
protocol_groups <- list(
  spatial = sim_groups(sim_group("long", 8, "exponential", range = 0.5), sim_group("short", 8, "exponential", range = 0.1),
                       sim_group("aniso", 7, "exponential", range = 0.25, anisotropy = 4), sim_group("none", 7, "none")),
  temporal = sim_groups(sim_group("strong", 6, "ar1", phi = 0.9), sim_group("weak", 6, "ar1", phi = 0.3),
                        sim_group("seasonal", 6, "seasonal", period = 12), sim_group("trend", 6, "trend_shift"),
                        sim_group("none", 6, "none")),
  spatiotemporal = sim_groups(sim_group("sep_long", 6, "separable_ar1", range = 0.5, phi = 0.9),
                              sim_group("sep_short", 6, "separable_seasonal", range = 0.1, period = 12),
                              sim_group("move", 6, "moving", range = 0.25, velocity = candidates$velocity),
                              sim_group("common", 6, "temporal_common", phi = 0.9), sim_group("none", 6, "none")))
sizes <- rbind(data.frame(family = "spatial", n = c(100, 400), n_times = 1),
               data.frame(family = "temporal", n = 1, n_times = c(120, 480)),
               data.frame(family = "spatiotemporal", n = c(50, 100), n_times = c(60, 120)))
cost <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
  sz <- sizes[i, ]
  secs <- vapply(1:5, function(k) {
    system.time(switch(sz$family,
      spatial = sim_spatial(sz$n, protocol_groups$spatial, seed = seeds[["cost"]] + k),
      temporal = sim_temporal(sz$n_times, protocol_groups$temporal, seed = seeds[["cost"]] + k),
      spatiotemporal = sim_spatiotemporal(sz$n, sz$n_times, protocol_groups$spatiotemporal, seed = seeds[["cost"]] + k)))[["elapsed"]]
  }, numeric(1))
  cbind(sz, seconds_per_replicate = stats::median(secs))
}))
cost$label <- ifelse(cost$family == "spatial", paste("n =", cost$n),
                     ifelse(cost$family == "temporal", paste("T =", cost$n_times), paste(cost$n, "x", cost$n_times)))
cost$minutes_200_reps_8_combinations <- cost$seconds_per_replicate * 200 * 8 / 60
write.csv(cost, file.path(out_dir, "cost.csv"), row.names = FALSE)
p <- ggplot(cost, aes(label, seconds_per_replicate)) + geom_col(fill = "#0d6b72", width = 0.6) +
  geom_text(aes(label = ifelse(seconds_per_replicate < 0.01, "< 0.01 s", sprintf("%.2f s", seconds_per_replicate))),
            vjust = -0.4, size = 3) +
  facet_wrap(~ family, scales = "free_x") +
  labs(title = "Time to generate one replicate (30 variables, protocol regimes)",
       subtitle = "Generation only; the correction and the four output methods (stage 6 pilot) will dominate the cost",
       x = NULL, y = "seconds per replicate (median of 5)")
save_fig(p, "sim_cost.png", 8, 3.6)

message("Figures written to ", fig_dir)
