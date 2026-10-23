# Stage 6 — Figures of the pilot, to support decision D8 and the review of D7 in time.
# Reads the output of analysis/06_pilot.R. Figures stay local.
#
# Outputs (results/figs/):
#   pilot_kmeans_example.png  a null spatial data set: k-means groups before and after correction
#   pilot_h1.png              spurious patterns in the null scenarios (H1)
#   pilot_recovery.png        recovery of the true structure, structured scenarios (H3/H4)
#   pilot_q3.png              contiguous classes: erased by the correction, kept by the broad band
#   pilot_d7_temporal.png     temporal data: global basis vs groups vs others (review of D7)
#   pilot_d8.png              replicates needed to detect each effect at stage 7 (D8)
#   pilot_crossed.png         D14: recovery when the true blocks cross the regime groups
#   pilot_between.png         D14: recovery with correlation 0.3 between groups

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
source(file.path("analysis", "d4_scenarios.R"))
fig_dir <- file.path("results", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
theme_set(theme_minimal(base_size = 10))
pal <- c(teal = "#0d6b72", brown = "#9a5806", grey = "#586972", green = "#2c7340", purple = "#7a4f9a")
app_lev <- c("none", "independent", "default", "groups", "variables", "whiten", "broad", "sbss", "maf")
app_cols <- c(none = "grey60", independent = pal[["brown"]], default = pal[["teal"]], groups = pal[["green"]],
              variables = "#d07a2b", whiten = pal[["purple"]], broad = "#c9a227", sbss = "#b23a48", maf = "#5b7fb8")

summ <- read.csv(file.path("results", "pilot_summary.csv"), stringsAsFactors = FALSE)
summ$approach <- factor(summ$approach, app_lev)
summ$fam <- sub(" (null|structured|contiguous|crossed).*$", "", summ$label)
summ_all <- summ
summ <- summ[summ$between == 0 & summ$truth != "crossed", ]       # the original pilot truths
summ$fam <- factor(summ$fam, unique(summ$fam[order(summ$scenario)]))

# --- 1. Didactic example: k-means on a null spatial data set ---------------------------------
sc <- pilot_scenarios[pilot_scenarios$label == "S n=400 high null", ]
seed <- 6600L + 1000L * sc$id + 1L
s <- simulate(sc, seed, "null")
gr <- group_regimes(diagnose_variables(s))
cr <- correct(s, gr, seed = seed, data_basis = "union")   # as in the pilot, before D17
km <- function(Z) { set.seed(1); stats::kmeans(scale(Z), 3, nstart = 20)$cluster }
nb <- spdep::knn2nb(spdep::knearneigh(s$coords, 8))
same_nb <- function(cl) mean(unlist(lapply(seq_along(nb), function(i) cl[nb[[i]]] == cl[i])))
d_ex <- rbind(data.frame(s$coords, group = factor(km(s$X)), panel = "raw data"),
              data.frame(s$coords, group = factor(km(cr$X)), panel = "after correction (default)"))
lab <- c(sprintf("raw data: %.0f%% of the neighbours in the same group", 100 * same_nb(km(s$X))),
         sprintf("after correction: %.0f%% (chance: about 33%%)", 100 * same_nb(km(cr$X))))
d_ex$panel <- factor(ifelse(d_ex$panel == "raw data", lab[1], lab[2]), lab)
p <- ggplot(d_ex, aes(x, y, colour = group)) + geom_point(size = 1.6) + coord_equal() + facet_wrap(~ panel) +
  scale_colour_manual(values = c(pal[["teal"]], pal[["brown"]], pal[["purple"]]), guide = "none") +
  labs(title = "A data set with no real association and no classes: k-means forced to find 3 groups",
       subtitle = "Raw data: the groups are regions of the map (the shared large-scale fields). After correction: no regional pattern left",
       x = NULL, y = NULL) + theme(axis.text = element_blank())
save_fig(p, "pilot_kmeans_example.png", 10, 4.5)

# --- 2. H1: spurious patterns in the null scenarios ------------------------------------------
h1_metrics <- c(pca = "first_above_broken_stick", mds_euclidean = "procrustes_position", kmeans = "k_gap_above_1",
                som = "map_position_cor")
h1_lab <- c(pca = "PCA: 1st axis above broken-stick (share)", mds_euclidean = "MDS: Procrustes with the positions",
            kmeans = "k-means: gap chooses k > 1 (share)", som = "SOM: map vs position correlation")
h1 <- summ[summ$truth == "null" & paste(summ$method, summ$metric) %in% paste(names(h1_metrics), h1_metrics), ]
h1$panel <- factor(h1_lab[h1$method], h1_lab)
p <- ggplot(h1, aes(approach, mean, fill = approach)) + geom_col() + facet_grid(panel ~ fam, scales = "free_y", labeller = label_wrap_gen(22)) +
  scale_fill_manual(values = app_cols, guide = "none") +
  labs(title = "H1: patterns found where there is nothing (null scenarios: no association, no classes; 20 replicates)",
       subtitle = "Any value above 0 is spurious. Mean over replicates", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), strip.text.y = element_text(size = 7))
save_fig(p, "pilot_h1.png", 12, 8)

# --- 3. Recovery in the structured scenarios ----------------------------------------------------
rec_metrics <- c(pca = "subspace_distance", mds_euclidean = "procrustes_truth", kmeans = "ari_true_k", som = "ari")
rec_lab <- c(pca = "PCA: distance to the true axes (lower = better)", mds_euclidean = "MDS: Procrustes with the truth",
             kmeans = "k-means: ARI with the true classes", som = "SOM: ARI with the true classes")
rec <- summ[summ$truth == "structured" & paste(summ$method, summ$metric) %in% paste(names(rec_metrics), rec_metrics), ]
rec$panel <- factor(rec_lab[rec$method], rec_lab)
p <- ggplot(rec, aes(approach, median, colour = approach)) + geom_pointrange(aes(ymin = q05, ymax = q95), size = 0.3) +
  facet_grid(panel ~ fam, scales = "free_y", labeller = label_wrap_gen(22)) + scale_colour_manual(values = app_cols, guide = "none") +
  labs(title = "H3/H4: recovery of the true structure (structured scenarios: correlations in blocks + 3 random classes)",
       subtitle = "Median and 5-95% interval of 20 replicates", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), strip.text.y = element_text(size = 7))
save_fig(p, "pilot_recovery.png", 12, 8)

# --- 4. Q3: contiguous classes -----------------------------------------------------------------
q3 <- summ[summ$truth == "contiguous" & ((summ$method %in% c("kmeans", "constrained") & summ$metric == "ari_true_k") |
                                           (summ$method == "som" & summ$metric == "ari")), ]
q3$panel <- factor(c(kmeans = "k-means", som = "SOM", constrained = "constrained clustering (D14)")[q3$method],
                   c("k-means", "SOM", "constrained clustering (D14)"))
p <- ggplot(q3, aes(approach, median, colour = approach)) + geom_pointrange(aes(ymin = q05, ymax = q95), size = 0.3) +
  facet_grid(panel ~ fam) + scale_colour_manual(values = app_cols, guide = "none") +
  labs(title = "Q3: classes contiguous in space or time (ARI with the true classes)",
       subtitle = "The correction removes the large scale, and here the classes are the large scale; the broad band of scales keeps them",
       x = NULL, y = "ARI") + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_fig(p, "pilot_q3.png", 11, 5)

# --- 5. Review of D7 in time ---------------------------------------------------------------------
tm <- summ[summ$family == "temporal" & summ$truth != "null" & paste(summ$method, summ$metric) %in% paste(names(rec_metrics), rec_metrics) &
             summ$approach %in% c("none", "independent", "default", "groups", "whiten"), ]
tm$panel <- factor(rec_lab[tm$method], rec_lab)
p <- ggplot(tm, aes(approach, median, colour = approach)) + geom_pointrange(aes(ymin = q05, ymax = q95), size = 0.3) +
  facet_grid(panel ~ truth, scales = "free_y", labeller = label_wrap_gen(24)) + scale_colour_manual(values = app_cols, guide = "none") +
  labs(title = "Review of D7 in time (T T=480): the global basis (default) against the others",
       subtitle = "Median and 5-95% interval of 20 replicates", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(size = 7))
save_fig(p, "pilot_d7_temporal.png", 9, 8)

# --- 6. Support for D8: replicates needed ----------------------------------------------------------
ct <- read.csv(file.path("results", "pilot_contrasts.csv"), stringsAsFactors = FALSE)
key <- rbind(data.frame(method = names(h1_metrics), metric = h1_metrics, truth = "null"),
             data.frame(method = names(rec_metrics), metric = rec_metrics, truth = "structured"))
ct <- merge(ct, key)
ct <- ct[!grepl("between|crossed", ct$label), ]
ct <- ct[ct$approach %in% c("independent", "default", "groups") & is.finite(ct$replicates_needed), ]
ct$fam <- sub(" (null|structured|contiguous|crossed).*$", "", ct$label)
ct$what <- paste0(ct$method, ": ", ct$metric)
p <- ggplot(ct, aes(pmin(replicates_needed, 1000), what, colour = approach)) + geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = c(20, 100, 200), linetype = "dashed", colour = "grey60") + facet_wrap(~ fam, nrow = 1) +
  scale_x_log10() + scale_colour_manual(values = app_cols[c("independent", "default", "groups")], name = "vs no correction") +
  labs(title = "Replicates needed at stage 7 to detect each effect seen in the pilot (paired, 5% two-sided, power 80%)",
       subtitle = "Capped at 1000. Dashed lines: 20, 100 and 200 replicates", x = "replicates (log scale)", y = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
save_fig(p, "pilot_d8.png", 13, 5.5)
# --- 7. D14: blocks that cross the groups, and correlation between groups ---------------------
d14_fig <- function(truth, between, file, title) {
  d <- summ_all[summ_all$truth == truth & summ_all$between == between &
                  paste(summ_all$method, summ_all$metric) %in% paste(names(rec_metrics), rec_metrics), ]
  d$fam <- factor(d$fam, levels(summ$fam))
  d$approach <- factor(d$approach, app_lev)
  d$panel <- factor(rec_lab[d$method], rec_lab)
  p <- ggplot(d, aes(approach, median, colour = approach)) + geom_pointrange(aes(ymin = q05, ymax = q95), size = 0.3) +
    facet_grid(panel ~ fam, scales = "free_y", labeller = label_wrap_gen(22)) + scale_colour_manual(values = app_cols, guide = "none") +
    labs(title = title, subtitle = "Median and 5-95% interval of 20 replicates", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), strip.text.y = element_text(size = 7))
  save_fig(p, file, 12, 8)
}
d14_fig("crossed", 0, "pilot_crossed.png",
        "D14: true correlation blocks that cross the regime groups (the apparent association no longer has the shape of the true one)")
d14_fig("structured", 0.3, "pilot_between.png", "D14: true correlation 0.5 within groups and 0.3 between groups")
message("Figures written to ", fig_dir)
