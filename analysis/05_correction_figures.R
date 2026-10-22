# Stage 5 — Figures of the mechanism (H2a, H2b), of the second grouping level (D11) and
# of the correction on the grid (D7). Reads the outputs of analysis/05a_h2_mechanism.R
# and analysis/05b_correction.R (the grid figures are skipped while 05b has not run).
#
# Outputs (results/figs/):
#   corr_scheme.png        didactic: what a large-scale basis leaves of a rough or smooth shared pattern
#   corr_heatmaps.png      R_true, raw R, R with independent bases, R with shared bases
#   h2a_distance.png       between-basis error against the distance between the bases
#   h2b_group_size.png     error and capture by group size, roughness and residual
#   d11_distance.png       distance between individual bases: same field vs independent fields
#   corr_grid_error.png    error of the residual correlations by approach and scenario (05b)
#   corr_grid_pairs.png    the same by kind of pair (05b)
#   corr_grid_capture.png  capture against degrees of freedom by approach (05b)
#   corr_grid_vs_independent.png  paired difference of each approach from "one basis per
#                          variable", replicate by replicate, with +- 2 standard errors (05b)
# and results/correction_vs_independent.csv

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
fig_dir <- file.path("results", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
theme_set(theme_minimal(base_size = 10))
pal <- c(teal = "#0d6b72", brown = "#9a5806", grey = "#586972", green = "#2c7340", purple = "#7a4f9a")

# --- 1. Didactic scheme on a transect ------------------------------------------------------
# A pattern shared by two variables, rough (exponential) or smooth (Matern 2.5), along a
# transect. A large-scale basis (slow waves) removes the smooth part; what it leaves is
# still common to both variables and looks like association.
set.seed(5901)
t <- seq(0, 1, length.out = 300)
dd <- as.matrix(stats::dist(t))
shared <- list(rough = depprep:::.matern_cor(dd, 0.6, 0.5), smooth = depprep:::.matern_cor(dd, 0.6, 2.5))
B <- qr.Q(qr(cbind(1, sapply(1:8, function(i) cos(pi * i * t)))))
did1 <- do.call(rbind, lapply(names(shared), function(nm) {
  L <- as.numeric(t(chol(shared[[nm]] + diag(1e-6, 300))) %*% stats::rnorm(300)); L <- L / stats::sd(L)
  e1 <- stats::rnorm(300); e2 <- stats::rnorm(300)
  x <- sqrt(0.6) * L + sqrt(0.4) * e1; y <- sqrt(0.6) * L + sqrt(0.4) * e2
  captured <- as.numeric(B %*% crossprod(B, L)); left <- L - captured
  rx <- x - as.numeric(B %*% crossprod(B, x)); ry <- y - as.numeric(B %*% crossprod(B, y))
  lab <- sprintf("%s shared pattern: basis captures %.0f%%; corrected x and y still correlate %.2f (truth 0)",
                 nm, 100 * (1 - sum(left^2) / sum((L - mean(L))^2)), stats::cor(rx, ry))
  rbind(data.frame(panel = lab, t = t, value = L, what = "shared pattern"),
        data.frame(panel = lab, t = t, value = captured, what = "removed by the basis"),
        data.frame(panel = lab, t = t, value = left - 3, what = "left behind (common to x and y)"))
}))
did1$panel <- factor(did1$panel, unique(did1$panel))
did1$what <- factor(did1$what, c("shared pattern", "removed by the basis", "left behind (common to x and y)"))
p <- ggplot(did1, aes(t, value, colour = what)) + geom_line(linewidth = 0.7) + facet_wrap(~ panel, ncol = 1) +
  scale_colour_manual(values = c("grey60", pal[["teal"]], pal[["brown"]]), name = NULL) +
  labs(title = "Why part of the shared structure stays: a large-scale basis removes only the smooth part of a pattern",
       subtitle = "Two variables share one pattern along a transect (plus independent noise). The basis: 8 slow waves. The leftover (shifted down) is the same in both",
       x = "position", y = NULL) + theme(axis.text.y = element_blank(), legend.position = "bottom")
save_fig(p, "corr_scheme.png", 10, 5.5)

# --- 2. Didactic heatmaps ------------------------------------------------------------------
did <- readRDS(file.path("results", "h2", "didactic.rds"))
v <- colnames(did$R_raw)
hm <- do.call(rbind, lapply(c(R_true = "R_true", raw = "R_raw", `independent bases` = "R_independent",
                              `shared bases (groups)` = "R_shared"), function(nm) {
  m <- did[[nm]]; dimnames(m) <- list(v, v)
  data.frame(a = factor(rep(v, length(v)), v), b = factor(rep(v, each = length(v)), rev(v)), r = as.vector(m))
}))
hm$matrix <- factor(rep(c("1. truth (R_true)", "2. raw data", "3. independent bases", "4. shared bases (groups)"),
                        each = length(v)^2))
p <- ggplot(hm, aes(a, b, fill = r)) + geom_tile() + facet_wrap(~ matrix, nrow = 1) +
  scale_fill_gradient2(low = pal[["brown"]], mid = "white", high = pal[["teal"]], limits = c(-1, 1), name = "correlation") +
  labs(title = "Groups A and B share one smooth large-scale field; C has its own. True residual correlation: 0.5 within, 0 between",
       subtitle = "The raw data show A and B as one block (apparent association); both corrections remove most of it", x = NULL, y = NULL) +
  theme(axis.text = element_text(size = 5), axis.text.x = element_text(angle = 90))
save_fig(p, "corr_heatmaps.png", 13, 4)

# --- 3. H2a ----------------------------------------------------------------------------------
h2a <- readRDS(file.path("results", "h2", "h2a.rds"))
s2a <- do.call(rbind, lapply(split(h2a, list(h2a$j, h2a$rho)), function(x) rbind(
  data.frame(rho = x$rho[1], distance = mean(x$distance), approach = "different bases", m = mean(x$error_different),
             se = stats::sd(x$error_different) / sqrt(nrow(x))),
  data.frame(rho = x$rho[1], distance = mean(x$distance), approach = "one shared basis", m = mean(x$error_shared),
             se = stats::sd(x$error_shared) / sqrt(nrow(x))))))
s2a$rho_lab <- sprintf("true residual correlation %.1f", s2a$rho)
th <- data.frame(distance = seq(0, 0.75, by = 0.01))
th <- rbind(transform(th, rho = 0), transform(th, rho = 0.5))
th$theory <- -th$rho * th$distance * 40 / 360
th$rho_lab <- sprintf("true residual correlation %.1f", th$rho)
p <- ggplot(s2a, aes(distance, m, colour = approach)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(data = th, aes(distance, theory), inherit.aes = FALSE, linetype = "dashed", colour = pal[["grey"]]) +
  geom_pointrange(aes(ymin = m - 2 * se, ymax = m + 2 * se), size = 0.3) + facet_wrap(~ rho_lab) +
  scale_colour_manual(values = c(`different bases` = pal[["brown"]], `one shared basis` = pal[["teal"]]), name = NULL) +
  labs(title = "H2a: the error that comes only from using different bases",
       subtitle = "The large scale lies in both bases; 500 replicates per point (mean +- 2 SE). Dashed: theory, -rho x (k - sum cos^2) / (n - k)",
       x = "distance between the two bases (0 = same, 1 = orthogonal)", y = "error of the residual correlation") +
  theme(legend.position = "bottom")
save_fig(p, "h2a_distance.png", 10, 4)

# --- 4. H2b ----------------------------------------------------------------------------------
h2b <- read.csv(file.path("results", "h2b_summary.csv"))
l2b <- rbind(
  data.frame(h2b[c("smoothness", "m", "rho", "residual")], approach = "one basis per variable", error = h2b$error_independent_mean),
  data.frame(h2b[c("smoothness", "m", "rho", "residual")], approach = "shared basis of the group", error = h2b$error_shared_mean),
  data.frame(h2b[c("smoothness", "m", "rho", "residual")], approach = "whitening", error = h2b$error_whiten_mean))
l2b$field <- ifelse(l2b$smoothness == 0.5, "rough field (exponential)", "smooth field (Matern 2.5)")
l2b$rho_lab <- sprintf("rho = %.1f", l2b$rho)
p1 <- ggplot(l2b, aes(m, error, colour = approach, linetype = residual)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) + facet_grid(rho_lab ~ field) + scale_x_continuous(breaks = c(2, 4, 8, 16), trans = "log2") +
  scale_colour_manual(values = c(pal[["brown"]], pal[["teal"]], pal[["purple"]]), name = NULL) +
  labs(title = "H2b: the shared structure left behind dominates the error, and depends on roughness",
       subtitle = "Mean error of the residual correlation of a pair (200 replicates). A larger group gives the shared basis more information",
       x = "variables in the group", y = "error (estimated - true)", linetype = "residual") + theme(legend.position = "bottom")
save_fig(p1, "h2b_group_size.png", 10, 6)

# --- 5. D11 ----------------------------------------------------------------------------------
d11 <- readRDS(file.path("results", "h2", "d11.rds"))
thr <- read.csv(file.path("results", "d11_threshold.csv"))
best <- thr$threshold[which.max(thr$youden)]
d11$pair <- ifelse(d11$same_field, "same field", "independent fields")
d11$cond <- sprintf("n = %d, %s", d11$n, ifelse(d11$smoothness == 0.5, "rough", "smooth"))
p <- ggplot(d11, aes(pair, distance, fill = pair)) + geom_hline(yintercept = best, linetype = "dashed", colour = pal[["brown"]]) +
  geom_boxplot(outlier.size = 0.4) + facet_grid(rho ~ cond, labeller = labeller(rho = function(x) paste("rho =", x))) +
  scale_fill_manual(values = c(`same field` = "#e0eff0", `independent fields` = "#f8ecd9"), guide = "none") +
  labs(title = "D11: distance between the bases selected for each variable alone",
       subtitle = sprintf("Pairs of variables with the same large-scale field or independent fields of the same range. Dashed: threshold %.2f (largest Youden index)", best),
       x = NULL, y = "distance between the two bases")
save_fig(p, "d11_distance.png", 11, 5)

# --- 6. The grid (05b) -----------------------------------------------------------------------
f_grid <- file.path("results", "correction_summary.csv")
if (file.exists(f_grid)) {
  g <- read.csv(f_grid)
  lev_app <- c("none", "independent", "global", "groups", "groups_split", "type_split", "oracle_groups", "whiten")
  g$approach <- factor(g$approach, lev_app)
  g$label <- factor(g$label, unique(g$label))
  p <- ggplot(g, aes(approach, rms, fill = approach == "groups")) + geom_col() + facet_wrap(~ label, ncol = 6) +
    scale_fill_manual(values = c(`TRUE` = pal[["teal"]], `FALSE` = "grey75"), guide = "none") +
    labs(title = "Error of the residual correlations by approach (RMS over all pairs, mean over replicates)",
         subtitle = "Teal: the proposal adopted at stage 4 (Gaussian-mixture groups, one basis per group)", x = NULL, y = "RMS error") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
  save_fig(p, "corr_grid_error.png", 13, 7)
  long <- rbind(data.frame(g[c("label", "approach")], pair = "within a true group", error = g$within_group),
                data.frame(g[c("label", "approach")], pair = "across groups", error = g$across_groups),
                data.frame(g[c("label", "approach")], pair = "involving no dependence", error = g$involving_none))
  p <- ggplot(long[long$approach != "none", ], aes(approach, error, colour = pair, group = pair)) + geom_point() + geom_line() +
    facet_wrap(~ label, ncol = 6) + scale_colour_manual(values = c(pal[["teal"]], pal[["brown"]], pal[["grey"]]), name = NULL) +
    labs(title = "Error of the residual correlations by kind of pair", x = NULL, y = "RMS error") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), legend.position = "bottom")
  save_fig(p, "corr_grid_pairs.png", 13, 7.5)
  p <- ggplot(g[is.finite(g$capture), ], aes(kn_mean, capture, colour = approach)) + geom_point(size = 2) +
    facet_wrap(~ family) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
    labs(title = "What each approach removes (capture of the true large scale) and what it costs (k/n)",
         x = "mean k/n of the bases applied", y = "capture") + theme(legend.position = "bottom")
  save_fig(p, "corr_grid_capture.png", 11, 4.5)

  # paired differences from "one basis per variable" (same data in every replicate)
  reps <- readRDS(file.path("results", "correction_replicates.rds"))
  w <- reshape(reps[, c("scenario", "label", "family", "replicate", "approach", "rms")],
               idvar = c("scenario", "label", "family", "replicate"), timevar = "approach", direction = "wide")
  names(w) <- sub("rms.", "", names(w), fixed = TRUE)
  others <- c("global", "groups", "groups_split", "type_split", "oracle_groups", "whiten")
  vs <- do.call(rbind, lapply(split(w, w$scenario), function(d) do.call(rbind, lapply(others, function(ap) {
    x <- d[[ap]] - d$independent
    data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], approach = ap, n_rep = nrow(d),
               diff = mean(x), se = stats::sd(x) / sqrt(length(x)), rel = mean(x) / mean(d$independent))
  }))))
  write.csv(vs, file.path("results", "correction_vs_independent.csv"), row.names = FALSE)
  vs$label <- factor(vs$label, rev(levels(g$label)))
  vs$approach <- factor(vs$approach, others)
  p <- ggplot(vs[vs$approach != "whiten", ], aes(diff, label)) + geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbarh(aes(xmin = diff - 2 * se, xmax = diff + 2 * se), height = 0, colour = pal[["grey"]]) +
    geom_point(aes(colour = diff < 0), size = 1.8) + facet_wrap(~ approach, nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = pal[["teal"]], `FALSE` = pal[["brown"]]), labels = c(`TRUE` = "smaller error", `FALSE` = "larger error"),
                        name = "than one basis per variable") +
    labs(title = "Each approach minus 'one basis per variable' (error of the residual correlations; same data, replicate by replicate)",
         subtitle = "Points: mean difference; bars: +- 2 standard errors. Whitening is left out (differences up to -0.07 in space-time, +0.03 elsewhere)",
         x = "difference in RMS error (negative = better)", y = NULL) + theme(legend.position = "bottom")
  save_fig(p, "corr_grid_vs_independent.png", 13, 5.5)
}
message("Figures written to ", fig_dir)
