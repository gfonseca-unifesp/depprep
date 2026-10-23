# Before the protocol is frozen — the primary tests of H1, H3 and H4, their predicted
# directions and the replicates they need at stage 7, with Holm's correction.
#
# Families of output methods (by what they receive and how dependence distorts them):
#   correlation  PCA      spurious axis above the broken-stick | distance to the true axes
#   distances    MDS      configuration that maps the positions | Procrustes with the truth
#   prototypes   k-means and SOM: map that reproduces the positions (SOM) | ARI (k-means)
#   density      DBSCAN: more than one cluster in null data | HDBSCAN: ARI (noise as a group)
# For each family and data family (spatial, temporal, space-time): one primary metric in
# the null truth (H1) and one in the structured truth (H3, H4). The proposal is the
# default of correct(): R (D7) for PCA, the data matrix with each variable on its own
# basis (D17; the "independent" approach of the pilot) for the others.
#
# Effect sizes and spreads come from the pilot (paired differences, the same data in
# every replicate): results/pilot_replicates.rds, results/pilot_hdbscan/ (HDBSCAN and the
# signal-based PCA metric) and, for H3-R, results/correction_replicates.rds (stage 5).
# Replicates needed for 80% power at the minimum relevant effect (D8: 0.05 on 0-1
# metrics, 0.10 on shares of replicates; H3-R: 10% of the error of one basis per
# variable), two-sided, at the Holm level of the first step (alpha / number of tests of
# the hypothesis), the most demanding one.
#
# Outputs (not tracked): results/protocol_tests.csv, results/protocol_cost.csv

suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("analysis", "d4_scenarios.R"))

pr <- readRDS(file.path("results", "pilot_replicates.rds"))
hd <- do.call(rbind, lapply(list.files(file.path("results", "pilot_hdbscan"), "^scenario_.*rds$", full.names = TRUE), readRDS))
hd <- merge(hd, pilot_scenarios[, c("id", "label", "family", "truth")], by.x = "scenario", by.y = "id")
keep_cols <- c("scenario", "replicate", "label", "family", "truth", "approach", "method", "metric", "value")
long <- rbind(pr[pr$scenario <= 12, keep_cols],
              hd[hd$method %in% c("hdbscan", "dbscan") | hd$metric == "subspace_distance_signal", keep_cols])

fam <- data.frame(
  technique = c("correlation", "distances", "prototypes", "density"),
  h1_method = c("pca", "mds_euclidean", "som", "dbscan"),     # density, H1: DBSCAN (D18, amended 2026-10-23)
  h1_metric = c("first_above_broken_stick", "procrustes_position", "map_position_cor", "clusters_above_1"),
  h1_delta = c(0.10, 0.05, 0.05, 0.10),
  h3_method = c("pca", "mds_euclidean", "kmeans", "hdbscan"),
  h3_metric = c("subspace_distance_signal", "procrustes_truth", "ari_true_k", "ari"),
  h3_good = c(-1, 1, 1, 1),                                  # -1: lower is better
  proposal = c("default", "independent", "independent", "independent"),
  stringsAsFactors = FALSE)

z_need <- function(alpha) stats::qnorm(1 - alpha / 2) + stats::qnorm(0.8)
paired <- function(d, a, b, method, metric) {
  x <- d[d$approach == a & d$method == method & d$metric == metric, c("scenario", "replicate", "value")]
  y <- d[d$approach == b & d$method == method & d$metric == metric, c("scenario", "replicate", "value")]
  m <- merge(x, y, by = c("scenario", "replicate"))
  m$diff <- m$value.x - m$value.y
  m[is.finite(m$diff), ]
}
summarise_test <- function(m, delta, sign, alpha_holm) {
  if (!nrow(m)) return(data.frame(effect = NA, sd_max = NA, n_holm = NA, n_nominal = NA))
  by <- split(m$diff, m$scenario)
  sds <- vapply(by, function(v) if (length(v) > 1) stats::sd(v) else NA_real_, 0)
  sd_max <- max(sds, na.rm = TRUE)
  data.frame(effect = sign * mean(m$diff), sd_max = sd_max,
             n_holm = ceiling((z_need(alpha_holm) * sd_max / delta)^2),
             n_nominal = ceiling((z_need(0.05) * sd_max / delta)^2))
}
fams <- c(spatial = "spatial", temporal = "temporal", spatiotemporal = "spatiotemporal")
rows <- list()
# H1: in the null truth, the proposal against no correction (12 tests)
for (i in seq_len(nrow(fam))) for (f in fams) {
  d <- long[long$truth == "null" & long$family == f, ]
  m <- paired(d, "none", fam$proposal[i], fam$h1_method[i], fam$h1_metric[i])
  rows[[length(rows) + 1]] <- cbind(hypothesis = "H1", technique = fam$technique[i], family = f,
    contrast = paste("none -", fam$proposal[i]), metric = fam$h1_metric[i], delta = fam$h1_delta[i],
    summarise_test(m, fam$h1_delta[i], 1, 0.05 / (nrow(fam) * length(fams))))
}
# H4 (and H3-X for the data-matrix families): in the structured truth, the proposal
# against the shared basis of D7 ("default"; whitening in space-time)
for (i in seq_len(nrow(fam))) for (f in fams) {
  comp <- if (fam$technique[i] == "correlation") "independent" else "default"
  d <- long[long$truth == "structured" & long$family == f, ]
  m <- paired(d, fam$proposal[i], comp, fam$h3_method[i], fam$h3_metric[i])
  rows[[length(rows) + 1]] <- cbind(hypothesis = "H4", technique = fam$technique[i], family = f,
    contrast = paste(fam$proposal[i], "-", comp), metric = fam$h3_metric[i], delta = 0.05,
    summarise_test(m, 0.05, fam$h3_good[i], 0.05 / (nrow(fam) * length(fams))))
}
# H3-X: data matrix (k-means ARI), the proposal against no correction (non-inferiority,
# margin 0.05), the shared basis and the groups (superiority), by data family (9 tests)
for (f in fams) for (comp in c("none", "default", "groups")) {
  d <- long[long$truth == "structured" & long$family == f, ]
  m <- paired(d, "independent", comp, "kmeans", "ari_true_k")
  rows[[length(rows) + 1]] <- cbind(hypothesis = "H3-X", technique = "prototypes", family = f,
    contrast = paste("independent -", comp), metric = "ari_true_k", delta = 0.05,
    summarise_test(m, 0.05, 1, 0.05 / (3 * length(fams))))
}
# H3-R: RMS error of the correlations (stage 5 grid), the proposal (D7: global in space
# and time, whitening in space-time) against one basis per variable; minimum effect 10%
cr <- readRDS(file.path("results", "correction_replicates.rds"))
for (f in fams) {
  d <- cr[cr$family == f & cr$scenario <= 18, ]
  prop <- if (f == "spatiotemporal") "whiten" else "global"
  x <- d[d$approach == prop, c("scenario", "replicate", "rms")]; y <- d[d$approach == "independent", c("scenario", "replicate", "rms")]
  m <- merge(x, y, by = c("scenario", "replicate")); m$diff <- m$rms.x - m$rms.y
  delta <- 0.10 * mean(m$rms.y)
  rows[[length(rows) + 1]] <- cbind(hypothesis = "H3-R", technique = "correlation", family = f,
    contrast = paste(prop, "- independent"), metric = "rms (relative 10%)", delta = delta,
    summarise_test(m, delta, -1, 0.05 / length(fams)))
}
tests <- do.call(rbind, rows)
# predicted direction, fixed mechanically from the pilot: a cell gets a prediction ("+",
# the proposal wins) only when the pilot effect reached the minimum relevant effect;
# otherwise it is exploratory (two-sided, reported, not counted in H1 or H4)
tests$prediction <- ifelse(is.finite(tests$effect) & abs(tests$effect) >= tests$delta, ifelse(tests$effect > 0, "+", "-"),
                           "exploratory")
# H3-X against no correction is a non-inferiority test: the proposal may lose at most the
# margin (0.05 of ARI)
ni <- tests$hypothesis == "H3-X" & grepl("- none$", tests$contrast)
tests$prediction[ni] <- ifelse(tests$effect[ni] > -tests$delta[ni], "non-inferior", "inferior")
# stage 7 tests pool the paired differences of all the scenarios of a data family
# (stratified mean, each scenario weighted equally): the replicates needed per scenario
# are the total divided by the number of scenarios with that truth (spatial 11,
# temporal 6, space-time 2 at 50 x 60; the 100 x 120 ones add 10 replicates each)
n_scen <- c(spatial = 11, temporal = 6, spatiotemporal = 2)
tests$scenarios <- n_scen[tests$family]
tests$n_per_scenario <- ceiling(tests$n_holm / tests$scenarios)
tests$replicates_planned <- ifelse(tests$family == "spatiotemporal", 50L, 100L)
tests$enough <- tests$n_per_scenario <= tests$replicates_planned |
  (tests$family == "spatiotemporal" & tests$n_holm <= 2 * 50 + 2 * 10)
write.csv(tests, file.path("results", "protocol_tests.csv"), row.names = FALSE)

# cost of stage 7 by family: seconds per replicate (pilot, one core, every approach and
# method), scaled for the smaller sizes from the stage 5 timings
pt <- read.csv(file.path("results", "pilot_times.csv"))
sec <- function(pattern) mean(pt$seconds_replicate[grepl(pattern, pt$label)])
cost <- data.frame(
  family = c("spatial n=400", "spatial n=100", "temporal T=480", "temporal T=100-120", "space-time 50x60"),
  seconds = c(sec("^S n=400"), sec("^S n=400") * 2.5 / 7, sec("^T T=480"), sec("^T T=480") * 2.9 / 6.8, sec("^ST 50x60")),
  combinations = c(7 * 5, 4 * 5, 3 * 5, 3 * 5, 2 * 5),     # scenarios x truths (ids: 2,4,12,14,17,18,19 | 1,3,11,13 | 6,8,20 | 5,7,21 | 9,15)
  replicates = c(100, 100, 100, 100, 50))
cost$core_hours <- cost$seconds * cost$combinations * cost$replicates / 3600
write.csv(cost, file.path("results", "protocol_cost.csv"), row.names = FALSE)

# D18: HDBSCAN against DBSCAN (radius by the silhouette), descriptive, for discussion
dens <- hd[hd$method %in% c("hdbscan", "dbscan") & hd$metric %in% c("ari", "ari_clustered", "noise_share", "n_clusters", "clusters_above_1"), ]
dc <- stats::aggregate(value ~ label + truth + approach + method + metric, data = dens, FUN = function(v) mean(v[is.finite(v)]))
write.csv(dc, file.path("results", "density_comparison.csv"), row.names = FALSE)
print(tests, digits = 3, row.names = FALSE)
print(cost, digits = 3, row.names = FALSE)
