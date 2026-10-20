# Stage 2 — Recovery of the simulated parameters (protocol section 7), for the
# candidate grid values of decision D4. 50 replicates per scenario.
#
# Seeds follow the scheme proposed for the protocol (section 9):
#   seed = base_seed + 1000 * scenario_index + replicate
#
# For each scenario the parameter is estimated both on the large-scale component
# (does the generator produce what it promises?) and on the observed variable
# (large scale + residual, share 0.6: how hard is it for stage 3?).
#
# Outputs:
#   results/sim_recovery.rds            all estimates, with seeds
#   results/sim_recovery_summary.csv    one line per scenario
#   results/figs/sim_recovery.png

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
base_seed <- 2200L
n_rep <- 50L
share <- 0.6
candidates <- list(range = c(0.05, 0.1, 0.25, 0.5), n = c(100, 400),
                   phi = c(0.3, 0.6, 0.9), n_times = c(120, 480))

# Exponential variogram fit; the observed variable needs a nugget (its residual).
# Singular or non-positive fits count as failures (NA).
fit_range <- function(coords, z, start, nugget = FALSE) {
  v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(coords, z = z), cutoff = 0.7)
  vz <- stats::var(z)
  model <- if (nugget) gstat::vgm(0.6 * vz, "Exp", start / 3, nugget = 0.4 * vz) else gstat::vgm(vz, "Exp", start / 3)
  fit <- tryCatch(suppressWarnings(gstat::fit.variogram(v, model)), error = function(e) NULL)
  if (is.null(fit) || isTRUE(attr(fit, "singular"))) return(NA_real_)
  a <- fit$range[fit$model == "Exp"]
  if (!length(a) || !is.finite(a) || a <= 0) return(NA_real_)
  3 * a
}

# --- spatial range -------------------------------------------------------------------
sc_s <- expand.grid(range = candidates$range, n = candidates$n)
spatial <- do.call(rbind, lapply(seq_len(nrow(sc_s)), function(k) {
  sc <- sc_s[k, ]
  g <- sim_groups(sim_group("f", 1, "exponential", range = sc$range, large_scale_share = share))
  do.call(rbind, lapply(seq_len(n_rep), function(r) {
    seed <- base_seed + 1000L * k + r
    s <- sim_spatial(sc$n, g, seed = seed)
    data.frame(parameter = "range", true = sc$range, size = sc$n, replicate = r, seed = seed,
               component = fit_range(s$coords, s$components$large[, 1], sc$range),
               observed = fit_range(s$coords, s$X[, 1], sc$range, nugget = TRUE))
  }))
}))

# --- temporal phi ----------------------------------------------------------------------
sc_t <- expand.grid(phi = candidates$phi, n_times = candidates$n_times)
offset <- nrow(sc_s)
temporal <- do.call(rbind, lapply(seq_len(nrow(sc_t)), function(k) {
  sc <- sc_t[k, ]
  g <- sim_groups(sim_group("f", 1, "ar1", phi = sc$phi, large_scale_share = share))
  do.call(rbind, lapply(seq_len(n_rep), function(r) {
    seed <- base_seed + 1000L * (offset + k) + r
    s <- sim_temporal(sc$n_times, g, seed = seed)
    lag1 <- function(x) stats::acf(x, 1, plot = FALSE)$acf[2]
    data.frame(parameter = "phi", true = sc$phi, size = sc$n_times, replicate = r, seed = seed,
               component = lag1(s$components$large[, 1]), observed = lag1(s$X[, 1]))
  }))
}))

# --- residual correlation ------------------------------------------------------------
ag <- sim_groups(sim_group("long_a", 6, "exponential", range = 0.5),
                 sim_group("long_b", 6, "exponential", shared_with = "long_a"),
                 sim_group("short", 6, "exponential", range = 0.1),
                 sim_group("aniso", 6, "exponential", range = 0.25, anisotropy = 4),
                 sim_group("none", 6, "none"))
offset <- offset + nrow(sc_t)
sc_r <- expand.grid(structure = c("identity", "blocks"), n = candidates$n, stringsAsFactors = FALSE)
corr <- do.call(rbind, lapply(seq_len(nrow(sc_r)), function(k) {
  sc <- sc_r[k, ]
  rt <- make_R_true(ag, sc$structure, within = 0.5, between = 0)
  expected <- sqrt(sum(((1 - rt^2)^2)[upper.tri(rt)]) * 2 / sc$n)
  do.call(rbind, lapply(seq_len(n_rep), function(r) {
    seed <- base_seed + 1000L * (offset + k) + r
    s <- sim_spatial(sc$n, ag, R_true = rt, seed = seed)
    data.frame(structure = sc$structure, size = sc$n, replicate = r, seed = seed,
               residual_error = norm(stats::cor(s$components$residual) - rt, "F"),
               raw_error = norm(stats::cor(s$X) - rt, "F"),
               expected_sampling_error = expected)
  }))
}))

est <- rbind(spatial, temporal)
# Range: relative bias (mean / true - 1). Phi: difference from the expected value;
# for the observed variable the expected lag-1 autocorrelation is share * phi.
summary_tab <- do.call(rbind, lapply(split(est, list(est$parameter, est$true, est$size), drop = TRUE), function(d) {
  is_range <- d$parameter[1] == "range"
  target_obs <- if (is_range) d$true[1] else share * d$true[1]
  bias <- function(x, target) if (is_range) mean(x, na.rm = TRUE) / target - 1 else mean(x, na.rm = TRUE) - target
  data.frame(parameter = d$parameter[1], true = d$true[1], size = d$size[1], n_rep = nrow(d),
             component_mean = mean(d$component, na.rm = TRUE), component_median = stats::median(d$component, na.rm = TRUE),
             component_bias = bias(d$component, d$true[1]), component_failed = sum(is.na(d$component)),
             observed_target = target_obs,
             observed_mean = mean(d$observed, na.rm = TRUE), observed_median = stats::median(d$observed, na.rm = TRUE),
             observed_bias = bias(d$observed, target_obs), observed_failed = sum(is.na(d$observed)),
             tolerance = if (is_range) 0.25 else 0.05)
}))
summary_tab$component_ok <- abs(summary_tab$component_bias) <= summary_tab$tolerance
summary_tab <- summary_tab[order(summary_tab$parameter, summary_tab$size, summary_tab$true), ]
corr_sum <- aggregate(cbind(residual_error, raw_error, expected_sampling_error) ~ structure + size, data = corr, FUN = mean)

saveRDS(list(estimates = est, correlation = corr, summary = summary_tab, correlation_summary = corr_sum,
             base_seed = base_seed, n_rep = n_rep, large_scale_share = share),
        file.path("results", "sim_recovery.rds"))
write.csv(summary_tab, file.path("results", "sim_recovery_summary.csv"), row.names = FALSE)
write.csv(corr_sum, file.path("results", "sim_recovery_correlation.csv"), row.names = FALSE)

long <- rbind(data.frame(est[, c("parameter", "true", "size")], what = "large-scale component", value = est$component),
              data.frame(est[, c("parameter", "true", "size")], what = "observed variable (share 0.6)", value = est$observed))
target <- ifelse(long$parameter == "phi" & long$what != "large-scale component", share * long$true, long$true)
long$ratio <- ifelse(long$parameter == "range", long$value / target, long$value - target)
# display only: range ratios are shown between 0 and 3 (full values in the .rds)
long$ratio <- ifelse(long$parameter == "range", pmin(pmax(long$ratio, 0), 3), long$ratio)
long$label <- ifelse(long$parameter == "range", sprintf("%.2f", long$true), sprintf("%.1f", long$true))
long$size_label <- ifelse(long$parameter == "range", paste("range, n =", long$size), paste("phi, T =", long$size))
band <- unique(long[, c("parameter", "size_label")])
band$lo <- ifelse(band$parameter == "range", 0.75, -0.05)
band$hi <- ifelse(band$parameter == "range", 1.25, 0.05)
band$mid <- ifelse(band$parameter == "range", 1, 0)
long$size_label <- factor(long$size_label, c("range, n = 100", "range, n = 400", "phi, T = 120", "phi, T = 480"))
band$size_label <- factor(band$size_label, levels(long$size_label))
p <- ggplot(long, aes(label, ratio, fill = what)) +
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), inherit.aes = FALSE,
            fill = "#e2f1e5", alpha = 0.8) +
  geom_hline(data = band, aes(yintercept = mid), colour = "#2c7340") +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(0.8)) +
  facet_wrap(~ size_label, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("#0d6b72", "#f8ecd9"), name = NULL) +
  coord_cartesian(ylim = NULL) +
  labs(title = "Recovery of the simulated parameters (50 replicates per scenario)",
       subtitle = "Range: estimate / true (band: within 25%). Phi: estimate - expected (band: within 0.05); expected lag 1 of the observed variable = 0.6 x phi.
Range ratios above 3 are drawn at 3; failed variogram fits are left out (see the summary table)",
       x = "simulated value", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path("results", "figs", "sim_recovery.png"), p, width = 10, height = 6, dpi = 150, bg = "white")

print(summary_tab, digits = 2)
print(corr_sum, digits = 3)
