# Stage 7 — the pre-registered tests of protocol-v1 (interno/docs/protocol.md, tag
# protocol-v1, section 2.3) and the rule of decision D9 (section 11), on the output of
# analysis/07_full_grid.R and analysis/07b_h2a.R.
#
# Tests: paired differences (the same data in every replicate), oriented so that a
# positive value means the predicted winner wins. The estimate of a test is the mean over
# the scenarios of its data family (spatial, including the "fields" scenarios; temporal;
# space-time, both sizes), each scenario weighted equally; its standard error is the
# stratified one, sqrt(sum s_c^2 / n_c) / C; two-sided p from the normal reference (every
# scenario has at least 10 replicates). Holm within each hypothesis. A predicted effect is
# confirmed when the estimate reaches the minimum relevant effect (D8) and the
# Holm-adjusted p is below 0.05; a non-inferiority test when the Holm-adjusted one-sided p
# of "effect <= -margin" is below 0.05. Exploratory tests are reported, not counted.
# Truths: H1 on the null truth; H3 and H4 on the structured, crossed and between-0.3
# truths pooled (also reported apart). The perturbation scenarios are reported apart as
# robustness and do not enter the primary tests.
#
# Outputs (not tracked): results/stage7_tests.csv, results/stage7_hypotheses.csv,
#   results/stage7_d9.csv, results/stage7_q3.csv, results/stage7_summary.csv,
#   results/stage7_perturbations.csv, results/stage7_density.csv

suppressMessages(devtools::load_all(quiet = TRUE))
grid <- read.csv(file.path("results", "stage7", "grid.csv"), stringsAsFactors = FALSE)
grid$data_family <- ifelse(grid$family == "fields", "spatial", grid$family)
all <- do.call(rbind, lapply(grid$combo, function(k) readRDS(file.path("results", "stage7", sprintf("combo_%03d.rds", k)))))
all <- merge(all, grid[, c("combo", "base_id", "label", "data_family", "smoothness", "truth_label", "perturbation")], by = "combo")
main <- all[all$perturbation == "none", ]
structured_truths <- c("structured", "crossed", "between03")

summ <- stats::aggregate(value ~ combo + label + data_family + truth_label + perturbation + approach + method + metric,
                         data = all, FUN = function(v) mean(v[is.finite(v)]))
write.csv(summ, file.path("results", "stage7_summary.csv"), row.names = FALSE)

# paired differences of approach a minus approach b, by combination and replicate
paired <- function(d, a, b, method, metric) {
  x <- d[d$approach == a & d$method == method & d$metric == metric, c("combo", "replicate", "value")]
  y <- d[d$approach == b & d$method == method & d$metric == metric, c("combo", "replicate", "value")]
  m <- merge(x, y, by = c("combo", "replicate"))
  m$diff <- m$value.x - m$value.y
  m[is.finite(m$diff), ]
}
stratified <- function(m, sign) {
  if (!nrow(m)) return(c(estimate = NA, se = NA, scenarios = 0, replicates = 0))
  by <- split(sign * m$diff, m$combo)
  mu <- vapply(by, mean, 0); v <- vapply(by, function(x) if (length(x) > 1) stats::var(x) / length(x) else NA_real_, 0)
  c(estimate = mean(mu), se = sqrt(sum(v, na.rm = TRUE)) / length(mu), scenarios = length(mu), replicates = nrow(m))
}
families <- c("spatial", "temporal", "spatiotemporal")
test_row <- function(hyp, technique, fam, contrast, metric, delta, prediction, est, kind = "superiority") {
  data.frame(hypothesis = hyp, technique = technique, data_family = fam, contrast = contrast, metric = metric, delta = delta,
             prediction = prediction, kind = kind, estimate = est[["estimate"]], se = est[["se"]],
             scenarios = est[["scenarios"]], replicates = est[["replicates"]], stringsAsFactors = FALSE)
}
rows <- list()
add <- function(r) rows[[length(rows) + 1]] <<- r
pred <- function(set) function(tech, fam) if (paste(tech, fam) %in% set) "+" else "exploratory"

# --- H1: null truth, the proposal against no correction ---------------------------------------
h1_pred <- pred(c(paste("correlation", families), paste("distances", families), paste("prototypes", families),
                  "density spatial", "density temporal"))
h1 <- data.frame(technique = c("correlation", "distances", "prototypes", "density"),
                 method = c("pca", "mds_euclidean", "som", "dbscan"),
                 metric = c("first_above_broken_stick", "procrustes_position", "map_position_cor", "clusters_above_1"),
                 proposal = c("default", "independent", "independent", "independent"), delta = c(0.10, 0.05, 0.05, 0.10))
for (i in seq_len(nrow(h1))) for (f in families) {
  d <- main[main$truth_label == "null" & main$data_family == f, ]
  add(test_row("H1", h1$technique[i], f, paste("none -", h1$proposal[i]), h1$metric[i], h1$delta[i], h1_pred(h1$technique[i], f),
               stratified(paired(d, "none", h1$proposal[i], h1$method[i], h1$metric[i]), 1)))
}
# --- H4: structured truths, the proposal of each family's output against its comparator ----
h4_pred <- pred(c("distances temporal", "prototypes spatial", "prototypes temporal", "density spatial", "density temporal",
                  "density spatiotemporal"))
h4 <- data.frame(technique = c("correlation", "distances", "prototypes", "density"),
                 method = c("pca", "mds_euclidean", "kmeans", "hdbscan"),
                 metric = c("subspace_distance_signal", "procrustes_truth", "ari_true_k", "ari"),
                 proposal = c("default", "independent", "independent", "independent"),
                 comparator = c("independent", "default", "default", "default"), sign = c(-1, 1, 1, 1))
for (i in seq_len(nrow(h4))) for (f in families) {
  d <- main[main$truth_label %in% structured_truths & main$data_family == f, ]
  add(test_row("H4", h4$technique[i], f, paste(h4$proposal[i], "-", h4$comparator[i]), h4$metric[i], 0.05, h4_pred(h4$technique[i], f),
               stratified(paired(d, h4$proposal[i], h4$comparator[i], h4$method[i], h4$metric[i]), h4$sign[i])))
}
# --- H3-X: k-means ARI, one basis per variable against none (non-inferiority), the shared
# basis and the groups (superiority) --------------------------------------------------------------
for (f in families) for (comp in c("none", "default", "groups")) {
  d <- main[main$truth_label %in% structured_truths & main$data_family == f, ]
  ni <- comp == "none"
  add(test_row("H3-X", "prototypes", f, paste("independent -", comp), "ari_true_k", 0.05,
               if (ni) "non-inferior" else if (f == "spatiotemporal") "exploratory" else "+",
               stratified(paired(d, "independent", comp, "kmeans", "ari_true_k"), 1), if (ni) "non-inferiority" else "superiority"))
}
# --- H3-R: RMS error of the correlations, the proposal for R against one basis per variable;
# minimum effect 10% of the error of one basis per variable ---------------------------------------
for (f in families) {
  d <- main[main$truth_label %in% structured_truths & main$data_family == f, ]
  rms_ind <- mean(d$value[d$approach == "independent" & d$method == "correlation" & d$metric == "rms"], na.rm = TRUE)
  add(test_row("H3-R", "correlation", f, "independent - default", "rms", 0.10 * rms_ind,
               if (f == "spatiotemporal") "+" else "exploratory",
               stratified(paired(d, "independent", "default", "correlation", "rms"), 1)))
}
# --- H2b: rough minus smooth error of the proposal's R (not paired: different data) ----------
rough_smooth <- list(spatial = list(rough = c(1, 2, 3, 4, 17), smooth = c(11, 12, 13, 14, 18)),
                     spatiotemporal = list(rough = c(9, 10), smooth = c(15, 16)))
for (f in names(rough_smooth)) {
  d <- main[main$truth_label %in% structured_truths & main$approach == "default" & main$method == "correlation" &
              main$metric == "rms", ]
  r <- d[d$base_id %in% rough_smooth[[f]]$rough, ]; s <- d[d$base_id %in% rough_smooth[[f]]$smooth, ]
  mr <- vapply(split(r$value, r$combo), mean, 0); ms <- vapply(split(s$value, s$combo), mean, 0)
  vr <- vapply(split(r$value, r$combo), function(x) stats::var(x) / length(x), 0)
  vs <- vapply(split(s$value, s$combo), function(x) stats::var(x) / length(x), 0)
  est <- c(estimate = mean(mr) - mean(ms), se = sqrt(sum(vr) / length(mr)^2 + sum(vs) / length(ms)^2),
           scenarios = length(mr) + length(ms), replicates = nrow(r) + nrow(s))
  add(test_row("H2b", "correlation", f, "rough - smooth", "rms", 0.05, "+", est))
}
tests <- do.call(rbind, rows)

# p-values, Holm within each hypothesis (H1: over its predicted tests, as registered; the
# others over every listed test), and confirmation
tests$p <- ifelse(tests$kind == "non-inferiority",
                  stats::pnorm((tests$estimate + tests$delta) / tests$se, lower.tail = FALSE),
                  2 * stats::pnorm(abs(tests$estimate) / tests$se, lower.tail = FALSE))
tests$p[!is.finite(tests$p)] <- ifelse(tests$estimate[!is.finite(tests$p)] != 0, 0, 1)   # se = 0 (e.g. all replicates equal)
tests$holm_set <- tests$hypothesis != "H1" | tests$prediction != "exploratory"
tests$p_holm <- NA_real_
for (h in unique(tests$hypothesis)) {
  i <- which(tests$hypothesis == h & tests$holm_set)
  tests$p_holm[i] <- stats::p.adjust(tests$p[i], "holm")
}
tests$confirmed <- ifelse(tests$prediction == "exploratory", NA,
  ifelse(tests$kind == "non-inferiority", tests$p_holm < 0.05,
         tests$estimate >= tests$delta & tests$p_holm < 0.05))
write.csv(tests, file.path("results", "stage7_tests.csv"), row.names = FALSE)

# --- Hypotheses, H2a and the D9 rule ---------------------------------------------------------------
h2a <- read.csv(file.path("results", "stage7_h2a_summary.csv"))
predicted <- function(h) tests[tests$hypothesis == h & tests$prediction != "exploratory", ]
hyp <- data.frame(hypothesis = c("H1", "H2a", "H2b", "H4"),
                  predicted_tests = c(nrow(predicted("H1")), nrow(h2a), nrow(predicted("H2b")), nrow(predicted("H4"))),
                  confirmed_tests = c(sum(predicted("H1")$confirmed), sum(h2a$within_0.005 & h2a$shared_below_0.002),
                                      sum(predicted("H2b")$confirmed), sum(predicted("H4")$confirmed)))
hyp$confirmed <- hyp$predicted_tests == hyp$confirmed_tests
fam_ok <- function(h) vapply(families, function(f) {
  t <- tests[tests$hypothesis == h & tests$data_family == f, ]
  if (h == "H3-X") all(t$confirmed %in% TRUE) else any(t$prediction != "exploratory") && all(t$confirmed[t$prediction != "exploratory"] %in% TRUE)
}, TRUE)
d9 <- data.frame(data_family = families, h3r_confirmed = fam_ok("H3-R"), h3x_confirmed = fam_ok("H3-X"))
d9_method <- sum(d9$h3r_confirmed) >= 2 || sum(d9$h3x_confirmed) >= 2
hyp <- rbind(hyp, data.frame(hypothesis = c("H3-R", "H3-X"), predicted_tests = c(1, 3), confirmed_tests = c(sum(d9$h3r_confirmed), sum(d9$h3x_confirmed)),
                             confirmed = c(sum(d9$h3r_confirmed) >= 2, sum(d9$h3x_confirmed) >= 2)))
write.csv(hyp, file.path("results", "stage7_hypotheses.csv"), row.names = FALSE)
d9$outcome <- if (d9_method) "method" else "diagnostic study"
write.csv(d9, file.path("results", "stage7_d9.csv"), row.names = FALSE)

# --- Q3 warning: the registered criterion ---------------------------------------------------------
q <- main[main$method == "verify" & main$metric == "q3_warning", ]
q$gradient <- q$base_id == 19
q3 <- stats::aggregate(value ~ data_family + truth_label + gradient, data = q, FUN = mean)
names(q3)[names(q3) == "value"] <- "warning_rate"
q3$criterion <- ifelse(q3$gradient, "reported apart", ifelse(q3$truth_label == "contiguous", ">= 0.90", "<= 0.05"))
q3$met <- ifelse(q3$gradient, NA, ifelse(q3$truth_label == "contiguous", q3$warning_rate >= 0.90, q3$warning_rate <= 0.05))
write.csv(q3, file.path("results", "stage7_q3.csv"), row.names = FALSE)

# --- Robustness: the perturbations, and the density comparison (D18) ---------------------------
pert <- all[all$perturbation != "none" & all$approach %in% c("none", "independent", "default") &
              paste(all$method, all$metric) %in% c("kmeans ari_true_k", "pca subspace_distance_signal", "correlation rms"), ]
pt <- stats::aggregate(value ~ label + perturbation + approach + method + metric, data = pert, FUN = function(v) mean(v[is.finite(v)]))
base <- summ[summ$perturbation == "none" & summ$truth_label == "structured" & summ$label %in% unique(pert$label) &
               summ$approach %in% c("none", "independent", "default") &
               paste(summ$method, summ$metric) %in% c("kmeans ari_true_k", "pca subspace_distance_signal", "correlation rms"), ]
base$perturbation <- "none (structured)"
write.csv(rbind(pt, base[, names(pt)]), file.path("results", "stage7_perturbations.csv"), row.names = FALSE)
dens <- summ[summ$perturbation == "none" & summ$method %in% c("hdbscan", "dbscan") &
               summ$metric %in% c("ari", "ari_clustered", "noise_share", "n_clusters", "clusters_above_1"), ]
write.csv(dens, file.path("results", "stage7_density.csv"), row.names = FALSE)

print(tests[, c("hypothesis", "technique", "data_family", "contrast", "prediction", "estimate", "se", "delta", "p_holm", "confirmed")],
      digits = 3, row.names = FALSE)
print(hyp, row.names = FALSE)
print(d9, row.names = FALSE)
print(q3, digits = 3, row.names = FALSE)
