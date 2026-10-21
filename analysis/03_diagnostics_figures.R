# Stage 3 — Figures and schemes of the diagnosis, to support decision D5. Reads the
# output of analysis/03_diagnostics.R. Figures stay local (results/ is not tracked).
#
# Outputs (results/figs/):
#   diag_scheme.png                   how a signature is built
#   diag_variogram_example.png        nugget, sill and range on two simulated variables
#   diag_signature_by_regime.png      what distinguishes each regime
#   diag_confusion.png                true vs classified regime, by scenario
#   diag_scale_range.png, diag_scale_phi.png   range and phi estimates against the truth
#   diag_thresholds_dependence.png    false positives and negatives of "no dependence"
#   diag_thresholds_other.png         accuracy against the other D5 thresholds
# and results/diagnostics_thresholds.csv, results/diagnostics_bootstrap_coverage.csv

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
fig_dir <- file.path("results", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
theme_set(theme_minimal(base_size = 10))
pal <- c(teal = "#0d6b72", brown = "#9a5806", grey = "#586972", green = "#2c7340")

all <- readRDS(file.path("results", "diagnostics_variables.rds"))
boot <- readRDS(file.path("results", "diagnostics_bootstrap.rds"))
scen <- read.csv(file.path("results", "diagnostics", "scenarios.csv"), stringsAsFactors = FALSE)
all$label <- factor(all$label, scen$label)

# --- 1. Scheme ---------------------------------------------------------------------------
boxes <- data.frame(
  x = c(0, 2.6, 2.6, 2.6, 5.2, 7.6, 7.6),
  y = c(2, 3.4, 2, 0.6, 2, 2.7, 1.3),
  fill = c("in", "sp", "tm", "st", "sig", "out", "out"),
  label = c("One variable\n(its own values only)",
            "Space\nMoran's I; variogram:\nnugget, sill, range;\nanisotropy",
            "Time\nACF, AR(1) phi;\nperiodogram: period,\nstrength; trend",
            "Space-time\ntime slices, sites;\nshift indicator\n(moving pattern)",
            "Signature\none row of numbers\nper variable",
            "classify_regime()\nregime type",
            "Stage 4\ngroup similar\nsignatures"))
arrows <- data.frame(x = c(0.95, 0.95, 0.95, 3.65, 3.65, 3.65, 6.1, 6.1),
                     y = c(2, 2, 2, 3.4, 2, 0.6, 2, 2),
                     xend = c(1.55, 1.55, 1.55, 4.35, 4.35, 4.35, 6.75, 6.75),
                     yend = c(3.4, 2, 0.6, 2, 2, 2, 2.7, 1.3))
p <- ggplot() +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), colour = "grey55", arrow = arrow(length = unit(2, "mm"))) +
  geom_label(data = boxes, aes(x, y, label = label, fill = fill), size = 3, linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c(`in` = "white", sp = "#e0eff0", tm = "#e0eff0", st = "#e0eff0", sig = "#f8ecd9", out = "#e2f1e5")) +
  annotate("text", x = 3.8, y = -0.5, size = 3.2, colour = pal[["brown"]],
           label = "Never uses correlations between variables: each column is diagnosed alone") +
  coord_cartesian(xlim = c(-0.8, 8.5), ylim = c(-0.8, 4.2)) + theme_void(base_size = 10) +
  labs(title = "How diagnose_variables() builds a dependence signature")
save_fig(p, "diag_scheme.png", 9, 4.2)

# --- 2. Didactic variogram ---------------------------------------------------------------
ex <- do.call(rbind, lapply(c(short = 0.1, long = 0.5), function(r) {
  s <- sim_spatial(400, sim_groups(sim_group("f", 1, "exponential", range = r)), seed = 3901)
  f <- depprep:::.fit_variogram(s$coords, s$X[, 1])
  v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(s$coords, z = s$X[, 1]), cutoff = 0.7)
  line <- gstat::variogramLine(f$fit, maxdist = 0.7)
  rbind(data.frame(kind = "points", dist = v$dist, gamma = v$gamma),
        data.frame(kind = "fit", dist = line$dist, gamma = line$gamma)) |>
    transform(panel = sprintf("simulated range %.1f (estimated %.2f)", r, f$range), nugget = f$nugget, sill = f$sill, range = f$range)
}))
marks <- unique(ex[, c("panel", "nugget", "sill", "range")])
p <- ggplot(ex, aes(dist, gamma)) +
  geom_point(data = ex[ex$kind == "points", ], colour = pal[["grey"]]) +
  geom_line(data = ex[ex$kind == "fit", ], colour = pal[["teal"]], linewidth = 1) +
  geom_hline(data = marks, aes(yintercept = sill), linetype = "dashed", colour = pal[["brown"]]) +
  geom_hline(data = marks, aes(yintercept = nugget), linetype = "dotted", colour = pal[["brown"]]) +
  geom_vline(data = marks, aes(xintercept = range), linetype = "dashed", colour = pal[["green"]]) +
  geom_text(data = marks, aes(x = 0.68, y = sill, label = "sill"), vjust = -0.5, hjust = 1, size = 3.2, colour = pal[["brown"]]) +
  geom_text(data = marks, aes(x = 0.68, y = nugget, label = "nugget"), vjust = -0.5, hjust = 1, size = 3.2, colour = pal[["brown"]]) +
  geom_text(data = marks, aes(x = range, y = 0.05, label = "range"), hjust = -0.1, size = 3.2, colour = pal[["green"]]) +
  facet_wrap(~ panel) + expand_limits(y = 0) +
  labs(title = "Reading a variogram: nugget, sill and range",
       subtitle = "Points: empirical semivariance by distance. Line: fitted model. Both variables have 40% residual (nugget)",
       x = "distance (side of the domain = 1)", y = "semivariance")
save_fig(p, "diag_variogram_example.png", 10, 3.8)

# --- 3. Signature by regime ---------------------------------------------------------------
comp <- list(spatial = c("moran_I", "range", "nugget_ratio", "aniso_ratio"),
             temporal = c("acf1", "phi", "period_strength", "trend_r2", "ljung_p"),
             spatiotemporal = c("moran_I", "range", "acf1", "period_strength", "shift_indicator"))
pick <- c(spatial = "S n=400 high", temporal = "T T=480 high", spatiotemporal = "ST 100x120")
sig <- do.call(rbind, lapply(names(comp), function(fam) {
  d <- all[all$label == pick[[fam]], ]
  do.call(rbind, lapply(comp[[fam]], function(cc) {
    m <- tapply(d[[cc]], d$true_type, stats::median, na.rm = TRUE)
    data.frame(family = fam, component = cc, true_type = names(m), value = as.numeric(m),
               scaled = as.numeric((m - min(m, na.rm = TRUE)) / diff(range(m, na.rm = TRUE))))
  }))
}))
sig$family <- factor(sig$family, names(comp))
p <- ggplot(sig, aes(component, true_type, fill = scaled)) + geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = formatC(value, format = "g", digits = 2)), size = 3) +
  facet_wrap(~ family, scales = "free", nrow = 1) +
  scale_fill_gradient(low = "white", high = pal[["teal"]], name = "relative\nvalue", labels = NULL) +
  labs(title = "Typical signature of each regime (median over all replicates)",
       subtitle = sprintf("Scenarios: %s. Colour: position of the value within each column", paste(pick, collapse = "; ")),
       x = NULL, y = "true regime") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p, "diag_signature_by_regime.png", 11, 4)

# --- 4. Confusion matrices ---------------------------------------------------------------
conf <- as.data.frame(table(label = all$label, true = all$true_type, classified = all$type))
conf <- conf[conf$Freq > 0 | conf$true == conf$classified, ]
tot <- aggregate(Freq ~ label + true, data = conf, FUN = sum)
conf <- merge(conf, setNames(tot, c("label", "true", "total")))
conf <- conf[conf$total > 0, ]
conf$pct <- 100 * conf$Freq / conf$total
p <- ggplot(conf, aes(classified, true, fill = pct)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct >= 1, sprintf("%.0f", pct), "")), size = 2.6) +
  facet_wrap(~ label, scales = "free", ncol = 4) +
  scale_fill_gradient(low = "white", high = pal[["teal"]], limits = c(0, 100), name = "% of the\ntrue regime") +
  labs(title = "Regime classified by the diagnosis vs true regime (200 replicates per scenario)",
       subtitle = "Each row sums to 100%. The diagonal is the share classified correctly", x = "classified as", y = "true regime") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7), axis.text.y = element_text(size = 7))
save_fig(p, "diag_confusion.png", 11, 7.5)

# --- 5. Scale uncertainty -------------------------------------------------------------------
sp <- all[!is.na(all$true_range) & all$family != "temporal", ]
sp$ratio <- pmin(pmax(sp$range / sp$true_range, 0), 4)
sp$size <- ifelse(sp$family == "spatial", paste0("S n=", sub("S n=([0-9]+).*", "\\1", sp$label)), as.character(sp$label))
tm <- all[!is.na(all$true_phi) & all$family == "temporal", ]
tm$diff <- tm$phi - tm$share * tm$true_phi
tm$size <- sub(" (high|low)$", "", tm$label)
p1 <- ggplot(sp, aes(factor(true_range), ratio, fill = size)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 1.25, fill = "#e2f1e5") +
  geom_hline(yintercept = 1, colour = pal[["green"]]) +
  geom_boxplot(outlier.size = 0.4, width = 0.7, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("#f8ecd9", "#e0eff0", "#cfe3e5", "#9fc9cc"), name = NULL) +
  labs(title = "Range estimated / simulated (observed variables)", subtitle = "Band: within 25%. Values above 4 drawn at 4",
       x = "simulated range", y = NULL)
p2 <- ggplot(tm, aes(factor(true_phi), diff, fill = size)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.05, ymax = 0.05, fill = "#e2f1e5") +
  geom_hline(yintercept = 0, colour = pal[["green"]]) +
  geom_boxplot(outlier.size = 0.4, width = 0.7, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("#f8ecd9", "#9fc9cc"), name = NULL) +
  labs(title = "phi estimated - expected (observed variables)", subtitle = "Expected = 0.6 x phi. Band: within 0.05",
       x = "simulated phi", y = NULL)
ggsave(file.path(fig_dir, "diag_scale_range.png"), p1, width = 8, height = 4, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "diag_scale_phi.png"), p2, width = 6, height = 4, dpi = 150, bg = "white")

# --- 6. Thresholds (D5) ---------------------------------------------------------------------
pvals <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2)
dep <- do.call(rbind, lapply(pvals, function(pt) {
  do.call(rbind, lapply(split(all, all$family), function(d) {
    lab <- classify_regime(structure(d, class = c("depprep_signature", "data.frame")), list(p_dependence = pt))
    none_true <- d$true_type == "none"
    data.frame(family = d$family[1], threshold = pt,
               false_positive = mean(lab[none_true] != "none"),
               false_negative = mean(lab[!none_true] == "none"))
  }))
}))
write.csv(dep, file.path("results", "diagnostics_thresholds.csv"), row.names = FALSE)
dl <- rbind(data.frame(dep[, c("family", "threshold")], error = "dependence missed (false negative)", rate = dep$false_negative),
            data.frame(dep[, c("family", "threshold")], error = "dependence invented (false positive)", rate = dep$false_positive))
p <- ggplot(dl, aes(threshold, rate, colour = error)) + geom_line(linewidth = 0.9) + geom_point() +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = pal[["grey"]]) +
  scale_x_log10(breaks = pvals) + scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ family) + scale_colour_manual(values = c(pal[["brown"]], pal[["teal"]]), name = NULL) +
  labs(title = "Declaring 'no dependence': the p-value threshold (all scenarios together)",
       subtitle = "Dashed line: candidate threshold 0.05", x = "p-value threshold (log scale)", y = "share of variables") +
  theme(legend.position = "bottom")
save_fig(p, "diag_thresholds_dependence.png", 10, 4)

grid_acc <- function(fam, types, col, values, key) {
  d <- all[all$family == fam & all$true_type %in% types, ]
  do.call(rbind, lapply(values, function(v) {
    th <- default_thresholds(); th[[key]] <- v
    lab <- classify_regime(structure(d, class = c("depprep_signature", "data.frame")), th)
    data.frame(threshold_name = key, value = v, accuracy = mean(lab == d$true_type))
  }))
}
other <- rbind(
  grid_acc("spatial", c("isotropic", "anisotropic"), "aniso_ratio", c(1.25, 1.5, 1.75, 2, 2.5, 3, 4), "aniso_ratio"),
  grid_acc("temporal", c("seasonal", "trend", "ar"), "period_strength", c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4), "seasonal_strength"),
  grid_acc("temporal", c("seasonal", "trend", "ar"), "trend_r2", c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), "trend_r2"),
  grid_acc("spatiotemporal", c("moving", "separable_ar1", "separable_seasonal"), "shift_indicator",
           c(0.02, 0.04, 0.06, 0.08, 0.1, 0.15), "shift"))
cand <- data.frame(threshold_name = names(default_thresholds())[-1], value = unlist(default_thresholds()[-1]))
write.csv(other, file.path("results", "diagnostics_thresholds_other.csv"), row.names = FALSE)
p <- ggplot(other, aes(value, accuracy)) + geom_line(colour = pal[["teal"]], linewidth = 0.9) + geom_point(colour = pal[["teal"]]) +
  geom_vline(data = cand, aes(xintercept = value), linetype = "dashed", colour = pal[["grey"]]) +
  facet_wrap(~ threshold_name, scales = "free_x", nrow = 1) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Accuracy of the other rules against their threshold (dependent variables only)",
       subtitle = "Dashed line: candidate value", x = "threshold", y = "classified correctly")
save_fig(p, "diag_thresholds_other.png", 11, 3.6)

cov <- rbind(data.frame(parameter = "range", aggregate(range_cover ~ label, data = boot, FUN = mean)) |> setNames(c("parameter", "label", "coverage")),
             data.frame(parameter = "phi", aggregate(phi_cover ~ label, data = boot, FUN = mean)) |> setNames(c("parameter", "label", "coverage")))
write.csv(cov, file.path("results", "diagnostics_bootstrap_coverage.csv"), row.names = FALSE)
message("Figures written to ", fig_dir)
