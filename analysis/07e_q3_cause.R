# Stage 7 — why the Q3 statistic depends on the data family (exploratory, decision D22).
#
# The Q3 statistic of verify() is the mean over the pairs of standardised variables of
# |cov_large|, the covariance of their projected (large-scale) parts on the union of the
# bases. In stage 7 it was higher without contiguous classes in time than in space, which
# is why one threshold failed (amendment protocol-v1.2). This script describes, in the
# truths without contiguous classes of the validation run of analysis/07d_q3_recalibration.R
# (new seeds, k/n of the union recorded with the statistic):
# * how the statistic grows with k/n, by family;
# * how far it is above white noise. For two independent white-noise variables projected on
#   a basis of k orthonormal columns, |cov_large| is about sqrt(2 / pi) * sqrt(k) / n. Above
#   it, the projected parts of independent autocorrelated variables are correlated by
#   chance: a large-scale field has few effective degrees of freedom, so two independent
#   ones look alike (Yule's nonsense correlation).
# No test and no decision: the output feeds the option of changing the statistic (D22).
#
# Outputs (not tracked): results/stage7_q3_cause_scenarios.csv, results/stage7_q3_cause_fit.csv

grid <- read.csv(file.path("results", "stage7", "grid.csv"), stringsAsFactors = FALSE)
grid <- grid[grid$perturbation == "none" & grid$base_id != 19 & grid$truth_label != "contiguous", ]
grid$data_family <- ifelse(grid$family == "fields", "spatial", grid$family)
grid$n_obs <- ifelse(grid$data_family == "spatial", grid$n, ifelse(grid$data_family == "temporal", grid$n_times, grid$n * grid$n_times))
q <- do.call(rbind, lapply(grid$combo, function(k) readRDS(file.path("results", "stage7_q3v", sprintf("combo_%03d.rds", k)))))
q <- merge(q, grid[, c("combo", "base_id", "label", "data_family", "truth_label", "n_obs")], by = "combo")
q <- q[is.finite(q$kn) & q$kn > 0, ]
q$white <- sqrt(2 / pi) * sqrt(q$kn * q$n_obs) / q$n_obs

scen <- do.call(rbind, lapply(split(q, q$base_id), function(d) data.frame(
  base_id = d$base_id[1], label = d$label[1], data_family = d$data_family[1], n_obs = d$n_obs[1], replicates = nrow(d),
  kn_mean = mean(d$kn), statistic_mean = mean(d$statistic), statistic_q95 = unname(stats::quantile(d$statistic, 0.95)),
  white_mean = mean(d$white), ratio_to_white = mean(d$statistic / d$white))))
scen <- scen[order(match(scen$data_family, c("spatial", "temporal", "spatiotemporal")), scen$n_obs, scen$base_id), ]

fit <- do.call(rbind, lapply(split(q, q$data_family), function(d) {
  m <- stats::lm(statistic ~ kn, data = d)
  data.frame(data_family = d$data_family[1], replicates = nrow(d), cor_kn = stats::cor(d$statistic, d$kn),
             intercept = unname(stats::coef(m)[1]), slope_kn = unname(stats::coef(m)[2]), r_squared = summary(m)$r.squared)
}))
write.csv(scen, file.path("results", "stage7_q3_cause_scenarios.csv"), row.names = FALSE)
write.csv(fit, file.path("results", "stage7_q3_cause_fit.csv"), row.names = FALSE)
print(fit, digits = 3, row.names = FALSE)
print(scen, digits = 3, row.names = FALSE)
