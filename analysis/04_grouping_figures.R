# Stage 4 — Figures and schemes of the grouping and of the shared bases, to support
# decision D6. Reads the output of analysis/04_grouping.R. Figures stay local.
#
# Outputs (results/figs/):
#   grp_scheme.png                 from signatures to one basis per group
#   grp_basis_space.png            three Moran eigenvector maps: broad, medium, fine
#   grp_basis_time.png             temporal candidates: harmonics and a slow dbMEM
#   grp_dendrogram_example.png     one replicate: tree, cut and true groups
#   grp_coassociation_example.png  the same replicate: co-association over subsamples
#   grp_ari.png                    agreement with the true groups, hclust vs mclust
#   grp_ngroups.png                number of groups chosen vs true
#   grp_stability.png              stability of correctly and wrongly assigned variables
#   grp_kn.png                     degrees of freedom used by group bases and unions
#   grp_capture_absorption.png     what the bases remove: structure vs residual
#   grp_kn_threshold.png           share of bases and unions above each k/n threshold
# and results/kn_threshold.csv

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
source(file.path("analysis", "d4_scenarios.R"))
fig_dir <- file.path("results", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
theme_set(theme_minimal(base_size = 10))
pal <- c(teal = "#0d6b72", brown = "#9a5806", grey = "#586972", green = "#2c7340")
grp_cols <- c("#0d6b72", "#9a5806", "#2c7340", "#7a4f9a", "#b23a48", "#586972")

grp <- readRDS(file.path("results", "grouping_replicates.rds"))
stab <- readRDS(file.path("results", "stability_replicates.rds"))
bas <- readRDS(file.path("results", "bases_replicates.rds"))
lev <- scenarios$label
grp$label <- factor(grp$label, lev)

# --- 1. Scheme ---------------------------------------------------------------------------
boxes <- data.frame(
  x = c(0, 2.4, 4.8, 7.2, 7.2, 9.6),
  y = c(2, 2, 2, 2.9, 1.1, 2),
  fill = c("in", "step", "step", "out", "none", "out"),
  label = c("Signatures\n(stage 3, one row\nper variable)",
            "Standardised\ncomponents\n(median, MAD)",
            "Grouping\nhclust or mclust;\nstability over\nsubsamples",
            "One basis per group\nMEM, harmonics,\ntrend, products",
            "'No dependence'\ngroup: empty basis",
            "Pairs across groups:\nunion of the two bases;\nwarning if k/n is large"))
arrows <- data.frame(x = c(0.95, 3.35, 5.75, 5.75, 8.35), y = c(2, 2, 2, 2, 2.9),
                     xend = c(1.45, 3.85, 6.2, 6.2, 8.65), yend = c(2, 2, 2.9, 1.1, 2.1))
p <- ggplot() +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), colour = "grey55", arrow = arrow(length = unit(2, "mm"))) +
  geom_label(data = boxes, aes(x, y, label = label, fill = fill), size = 3, linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c(`in` = "white", step = "#e0eff0", out = "#e2f1e5", none = "#f1f1f1")) +
  annotate("text", x = 4.8, y = 0, size = 3.2, colour = pal[["brown"]],
           label = "Same group -> same basis. Only the signatures and each group's own variance are used, never correlations between variables") +
  coord_cartesian(xlim = c(-0.8, 10.6), ylim = c(-0.4, 3.8)) + theme_void(base_size = 10) +
  labs(title = "Stage 4: from signatures to one shared basis per group")
save_fig(p, "grp_scheme.png", 10, 3.8)

# --- 2. What a basis is ---------------------------------------------------------------------
# a regular 20 x 20 grid shows the scales clearly; on scattered sites the maps look patchier
grid_xy <- as.matrix(expand.grid(x = seq(0.025, 0.975, by = 0.05), y = seq(0.025, 0.975, by = 0.05)))
mem <- depprep:::.spatial_candidates(grid_xy)
pick <- c(1, 12, 80)
mp <- do.call(rbind, lapply(pick, function(j) data.frame(grid_xy, value = mem[, j],
  panel = sprintf("MEM %d: %s scale", j, c("broad", "medium", "fine")[match(j, pick)]))))
mp$panel <- factor(mp$panel, unique(mp$panel))
p <- ggplot(mp, aes(x, y, fill = value)) + geom_tile() + coord_equal() + facet_wrap(~ panel) +
  scale_fill_gradient2(low = pal[["brown"]], mid = "grey95", high = pal[["teal"]], guide = "none") +
  labs(title = "A spatial basis: Moran eigenvector maps (dbMEM) on a regular grid of 400 sites",
       subtitle = "Each map is one column of the candidate basis. A group's basis is the few maps that describe its large-scale pattern",
       x = NULL, y = NULL) + theme(axis.text = element_blank())
save_fig(p, "grp_basis_space.png", 10, 3.8)
tt <- 1:120
tc <- cbind(depprep:::.temporal_candidates(tt, "seasonal", 12)[, c("sin1_p12", "cos1_p12")],
            depprep:::.temporal_candidates(tt, "ar")[, c("tMEM1", "tMEM5")])
tp <- data.frame(time = rep(tt, ncol(tc)), value = as.vector(tc),
                 term = factor(rep(c("harmonic sin, period 12", "harmonic cos, period 12", "slow dbMEM 1", "slow dbMEM 5"), each = 120),
                               c("harmonic sin, period 12", "harmonic cos, period 12", "slow dbMEM 1", "slow dbMEM 5")))
p <- ggplot(tp, aes(time, value)) + geom_line(colour = pal[["teal"]]) + facet_wrap(~ term, nrow = 1) +
  labs(title = "Temporal candidates: harmonics for seasonal groups, slow dbMEM for persistent ones (and trend)",
       x = "time step", y = NULL) + theme(axis.text.y = element_blank())
save_fig(p, "grp_basis_time.png", 10, 2.6)

# --- 3. One replicate step by step --------------------------------------------------------
dendro_segments <- function(h) {
  n <- length(h$order)
  x <- numeric(n); x[h$order] <- seq_len(n)
  pos <- numeric(n - 1)
  getx <- function(j) if (j < 0) x[-j] else pos[j]
  geth <- function(j) if (j < 0) 0 else h$height[j]
  do.call(rbind, lapply(seq_len(n - 1), function(i) {
    a <- h$merge[i, 1]; b <- h$merge[i, 2]
    xa <- getx(a); xb <- getx(b); pos[i] <<- (xa + xb) / 2
    data.frame(x = c(xa, xb, xa), y = c(geth(a), geth(b), h$height[i]),
               xend = c(xa, xb, xb), yend = rep(h$height[i], 3))
  }))
}
ex_k <- 2L; ex_r <- 1L                                      # S n=400 high, replicate 1
d3 <- readRDS(file.path("results", "diagnostics", sprintf("scenario_%02d.rds", ex_k)))
sig <- d3[d3$replicate == ex_r, ]
g <- group_regimes(sig)
tree <- attr(g, "tree")
dep <- g$group != "none"
seg <- dendro_segments(tree)
leaves <- data.frame(x = match(seq_along(tree$order), tree$order), variable = g$variable[dep], true = sig$group[dep], group = g$group[dep])
k <- attr(g, "n_groups")
cut_h <- if (k > 1) mean(sort(tree$height, decreasing = TRUE)[c(k - 1, k)]) else NA
p <- ggplot() + geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), colour = "grey40") +
  geom_text(data = leaves, aes(x, -0.15, label = variable, colour = true), angle = 90, hjust = 1, size = 2.8,
            key_glyph = "point") +
  geom_hline(yintercept = cut_h, linetype = "dashed", colour = pal[["brown"]]) +
  annotate("text", x = 0.5, y = cut_h, label = sprintf("cut: %d groups (largest mean silhouette, %.2f)", k, attr(g, "silhouette")),
           hjust = 0, vjust = -0.4, size = 3, colour = pal[["brown"]]) +
  scale_colour_manual(values = grp_cols, name = "true group") +
  scale_y_continuous(expand = expansion(mult = c(0.35, 0.05))) +
  labs(title = sprintf("One replicate (%s): Ward tree of the dependent variables", lev[ex_k]),
       subtitle = sprintf("%d variables classified 'no dependence' are kept apart, with no basis. Leaf colour: true group",
                          sum(!dep)), x = NULL, y = "merge height") +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom")
save_fig(p, "grp_dendrogram_example.png", 10, 5)

st_tab <- stab$table
first <- which(st_tab$scenario == ex_k & st_tab$replicate == ex_r & st_tab$method == "hclust")
keys <- unique(paste(st_tab$scenario, st_tab$replicate, st_tab$method))   # same order as the matrices
co <- stab$coassociation[[match(paste(ex_k, ex_r, "hclust"), keys)]]
tab_ex <- st_tab[first, ]
tab_ex <- tab_ex[tab_ex$group != "none", ]
ord <- order(tab_ex$true_group, tab_ex$group)
v <- tab_ex$variable[ord]
cd <- expand.grid(a = v, b = v, stringsAsFactors = FALSE)
cd$co <- co[cbind(match(cd$a, rownames(co)), match(cd$b, colnames(co)))]
cd$a <- factor(cd$a, v); cd$b <- factor(cd$b, rev(v))
amb <- tab_ex$variable[tab_ex$ambiguous]
p <- ggplot(cd, aes(a, b, fill = co)) + geom_tile() +
  scale_fill_gradient(low = "white", high = pal[["teal"]], limits = c(0, 1), name = "together in\nthis share of\nsubsamples") +
  labs(title = "Co-association over 50 subsamples of 80% of the sites (same replicate, hclust)",
       subtitle = sprintf("Variables ordered by true group. Ambiguous (stability < 0.7): %d of %d", length(amb), nrow(tab_ex)),
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 7))
save_fig(p, "grp_coassociation_example.png", 8, 7)

# --- 4. Agreement with the true groups ------------------------------------------------------
p <- ggplot(grp, aes(label, ari, fill = method)) + geom_boxplot(outlier.size = 0.4, width = 0.7, position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, colour = pal[["green"]]) +
  scale_fill_manual(values = c(hclust = "#e0eff0", mclust = "#f8ecd9"), name = NULL) +
  labs(title = "Agreement between estimated and true groups (adjusted Rand index, 200 replicates per scenario)",
       subtitle = "1 = identical partitions, 0 = no better than chance. All variables, 'no dependence' group included",
       x = NULL, y = "ARI") + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_fig(p, "grp_ari.png", 10, 4.2)
ng <- as.data.frame(table(label = grp$label, method = grp$method, diff = pmax(pmin(grp$n_groups - grp$true_groups, 3), -3)))
ng <- ng[ng$Freq > 0, ]
ng$diff <- factor(ng$diff, as.character(-3:3))
p <- ggplot(ng, aes(diff, label, fill = Freq / 200)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f", 100 * Freq / 200)), size = 2.8) + facet_wrap(~ method) +
  scale_fill_gradient(low = "white", high = pal[["teal"]], name = "% of\nreplicates", labels = scales::percent) +
  labs(title = "Number of groups chosen minus true number of dependent groups",
       subtitle = "0 = right number; negative = groups merged; positive = groups split (clipped at +-3)",
       x = "chosen - true", y = NULL)
save_fig(p, "grp_ngroups.png", 9, 4.2)

# --- 5. Stability and wrong assignments ---------------------------------------------------
sd <- st_tab[st_tab$group != "none" & !is.na(st_tab$stability), ]
sd$label <- factor(sd$label, lev)
sd$assigned <- ifelse(sd$wrong, "wrong group", "right group")
p <- ggplot(sd, aes(label, stability, fill = assigned)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = pal[["brown"]]) +
  geom_boxplot(outlier.size = 0.4, width = 0.7, position = position_dodge(0.8)) + facet_wrap(~ method, ncol = 1) +
  scale_fill_manual(values = c(`right group` = "#e0eff0", `wrong group` = "#f8ecd9"), name = NULL) +
  labs(title = "Stability of each variable (mean co-association with its group, 50 subsamples)",
       subtitle = "20 replicates per scenario. Dashed line: 0.7, below which a variable is flagged as ambiguous",
       x = NULL, y = "stability") + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_fig(p, "grp_stability.png", 10, 6.5)

# --- 6. Degrees of freedom -------------------------------------------------------------------
b <- bas$bases[!is.na(bas$bases$kn), ]
u <- bas$unions
kd <- rbind(data.frame(label = b$label, family = b$family, grouping = b$grouping, what = paste("group:", b$group_type), kn = b$kn),
            data.frame(label = u$label, family = u$family, grouping = u$grouping, what = "union of two groups", kn = u$kn))
kd$label <- factor(kd$label, lev)
p <- ggplot(kd[kd$grouping == "estimated", ], aes(what, kn)) +
  geom_hline(yintercept = c(0.1, 0.2, 0.3), linetype = "dashed", colour = "grey60") +
  geom_boxplot(fill = "#e0eff0", outlier.size = 0.4) + coord_flip() + facet_wrap(~ label, ncol = 5, scales = "free_y") +
  labs(title = "Share of degrees of freedom used by each basis (k/n), estimated grouping, 50 replicates per scenario",
       subtitle = "Dashed lines: candidate thresholds 0.1, 0.2, 0.3. Selection stops at 0.25 per group",
       x = NULL, y = "k / n") + theme(axis.text.y = element_text(size = 7))
save_fig(p, "grp_kn.png", 12, 6.5)

# --- 7. Capture and absorption ---------------------------------------------------------------
ca <- b[is.finite(b$capture), ]
ca$label <- factor(ca$label, lev)
p <- ggplot(ca, aes(capture, absorption, colour = group_type, shape = grouping)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_point(alpha = 0.5, size = 1.3) +
  facet_wrap(~ family) + scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_shape_manual(values = c(estimated = 16, oracle = 4)) +
  labs(title = "What each basis removes: the true large-scale structure (capture) vs the true residual (absorption)",
       subtitle = "Capture: share of the large-scale component removed. Absorption: share of the residual removed beyond k/n. Ideal: right, near 0",
       x = "capture", y = "absorption", colour = "group type", shape = "grouping")
save_fig(p, "grp_capture_absorption.png", 11, 4.4)

# --- 8. Support for D6: the k/n threshold ---------------------------------------------------
ths <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)
kt <- do.call(rbind, lapply(ths, function(t) {
  e <- kd[kd$grouping == "estimated", ]
  e$kind <- ifelse(e$what == "union of two groups", "unions of two groups", "group bases")
  aggregate(kn ~ family + kind, data = e, FUN = function(x) mean(x > t)) |> transform(threshold = t)
}))
names(kt)[names(kt) == "kn"] <- "share_above"
write.csv(kt, file.path("results", "kn_threshold.csv"), row.names = FALSE)
p <- ggplot(kt, aes(threshold, share_above, colour = kind)) + geom_line(linewidth = 0.9) + geom_point() +
  facet_wrap(~ family) + scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = c(pal[["teal"]], pal[["brown"]]), name = NULL) +
  labs(title = "Share of bases that would get a warning, by k/n threshold (estimated grouping, all scenarios)",
       x = "k/n threshold", y = "share above the threshold") + theme(legend.position = "bottom")
save_fig(p, "grp_kn_threshold.png", 10, 3.8)
message("Figures written to ", fig_dir)
