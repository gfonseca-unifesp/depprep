# Stage 4 — Grouping of regimes and shared bases, on the D4 grid.
#
# Three parts, all on the very data of stage 3 (same seeds, regenerated on the fly):
#   1. grouping: group_regimes() on the stage 3 signatures of every replicate (200 per
#      scenario), with hclust and mclust; agreement with the true groups (ARI) and the
#      number of groups chosen;
#   2. stability: 20 replicates per scenario, B = 50 subsamples, with hclust and mclust;
#      co-association and the link between instability and wrong assignments;
#   3. bases: 50 replicates per scenario; one basis per group with the estimated and with
#      the true (oracle) grouping; k, k/n, capture and absorption of each basis, and k/n
#      of the union of every pair of groups.
#
# Seeds: the data use the stage 3 seeds (3300 + 1000 * scenario + replicate); the
# subsamples and the permutations of the basis selection use 4400 + 1000 * scenario +
# replicate. Runs in parallel and can be resumed (one file per scenario and part).
#
# Outputs (not tracked by git):
#   results/grouping/{grouping,stability,bases}_NN.rds
#   results/grouping_summary.csv   ARI and number of groups by scenario and method
#   results/stability_summary.csv  stability, ambiguity and wrong assignments
#   results/bases_summary.csv      k, k/n, capture and absorption by group type
#   results/union_summary.csv      k/n of the unions of pairs of groups

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

data_seed <- 3300L
stage_seed <- 4400L
n_stab_rep <- 20L
B_stab <- 50L
n_basis_rep <- 50L
in_dir <- file.path("results", "diagnostics")
out_dir <- file.path("results", "grouping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log")

signatures <- function(k) {
  d <- readRDS(file.path(in_dir, sprintf("scenario_%02d.rds", k)))
  names(d)[names(d) == "group"] <- "true_group"
  split(d, d$replicate)
}

# --- 1. Grouping of every replicate --------------------------------------------------------
one_grouping <- function(sig) {
  do.call(rbind, lapply(c("hclust", "mclust"), function(m) {
    t0 <- Sys.time()
    g <- group_regimes(sig, method = m)
    dep <- sig$true_type != "none"
    data.frame(scenario = sig$scenario[1], replicate = sig$replicate[1], method = m,
               ari = mclust::adjustedRandIndex(g$group, sig$true_group),
               ari_dependent = if (sum(dep) > 1) mclust::adjustedRandIndex(g$group[dep], sig$true_group[dep]) else NA,
               n_groups = attr(g, "n_groups"), true_groups = length(unique(sig$true_group[dep])),
               none_recall = mean(g$group[!dep] == "none"), none_precision = mean(!dep[g$group == "none"]),
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }))
}

# --- 2. Stability on a subsample of replicates ---------------------------------------------
one_stability <- function(sc, sig, method = "hclust") {
  r <- sig$replicate[1]
  s <- simulate(sc, data_seed + 1000L * sc$id + r)
  g <- group_regimes(sig, data = s, method = method, B = B_stab, seed = stage_seed + 1000L * sc$id + r)
  # a variable is wrongly assigned when its true group is not the most common true group
  # of the estimated group it was put in
  major <- tapply(sig$true_group, g$group, function(x) names(which.max(table(x))))
  list(table = data.frame(scenario = sc$id, replicate = r, method = method, variable = g$variable, true_group = sig$true_group,
                          true_type = sig$true_type, group = g$group, group_type = g$group_type,
                          stability = g$stability, ambiguous = g$ambiguous,
                          wrong = sig$true_group != major[g$group], stringsAsFactors = FALSE),
       coassociation = attr(g, "coassociation"))
}

# --- 3. Bases with the estimated and the oracle grouping -----------------------------------
group_bases <- function(s, members, types, periods, seed) {
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time),
                 spatiotemporal = list(coords = s$coords, index = s$index))
  bases <- list()
  rows <- list()
  for (g in setdiff(names(members), "none")) {
    v <- members[[g]]
    t0 <- Sys.time()
    per <- periods[[g]]
    if (!is.null(per) && is.na(per)) per <- NULL
    b <- tryCatch(do.call(build_group_basis, c(list(X_group = s$X[, v, drop = FALSE], type = types[[g]],
                                                    period = per, seed = seed), args)),
                  error = function(e) conditionMessage(e))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.character(b)) {                                  # recorded, not fatal
      rows[[g]] <- data.frame(group = g, group_type = types[[g]], n_vars = length(v), k = NA, kn = NA, r2 = NA,
                              capture = NA, absorption = NA, seconds = secs, error = b)
      next
    }
    ca <- .capture_absorption(b, s$components$large[, v, drop = FALSE], s$components$residual[, v, drop = FALSE])
    bases[[g]] <- b
    rows[[g]] <- data.frame(group = g, group_type = types[[g]], n_vars = length(v), k = attr(b, "k"),
                            kn = attr(b, "kn"), r2 = attr(b, "r2"), capture = ca[["capture"]],
                            absorption = ca[["absorption"]], seconds = secs, error = NA_character_)
  }
  pairs <- if (length(bases) > 1) utils::combn(names(bases), 2, simplify = FALSE) else list()
  unions <- do.call(rbind, lapply(pairs, function(p) {
    u <- suppressWarnings(union_basis(bases[[p[1]]], bases[[p[2]]]))
    data.frame(group_a = p[1], group_b = p[2], type_a = types[[p[1]]], type_b = types[[p[2]]],
               k = attr(u, "k"), kn = attr(u, "kn"))
  }))
  list(bases = do.call(rbind, rows), unions = unions)
}

one_bases <- function(sc, sig) {
  r <- sig$replicate[1]
  s <- simulate(sc, data_seed + 1000L * sc$id + r)
  seed <- stage_seed + 1000L * sc$id + r
  g <- group_regimes(sig)
  est <- group_bases(s, split(g$variable, g$group),
                     types = tapply(g$group_type, g$group, `[`, 1),
                     periods = as.list(c(attr(g, "period"), none = NA)), seed = seed)
  # oracle: the true groups, their true types and periods
  truth <- unique(sig[, c("true_group", "true_type", "true_period")])
  types <- stats::setNames(as.list(truth$true_type), truth$true_group)
  names(types)[truth$true_type == "none"] <- "none"
  periods <- stats::setNames(as.list(truth$true_period), names(types))
  members <- split(sig$variable, ifelse(sig$true_type == "none", "none", sig$true_group))
  ora <- group_bases(s, members, types, periods, seed)
  tag <- function(x, grouping) if (!is.null(x)) cbind(scenario = sc$id, replicate = r, grouping = grouping, x)
  list(bases = rbind(tag(est$bases, "estimated"), tag(ora$bases, "oracle")),
       unions = rbind(tag(est$unions, "estimated"), tag(ora$unions, "oracle")))
}

# --- Run -------------------------------------------------------------------------------
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE)); NULL }, getwd())
clusterExport(cl, c("data_seed", "stage_seed", "B_stab", "spatial_groups", "temporal_groups", "st_groups",
                    "truth_table", "simulate", "one_grouping", "one_stability", "group_bases", "one_bases"))
cat(format(Sys.time()), "start\n", file = log_file, append = TRUE)
run_part <- function(part, k, fun, reps) {
  f <- file.path(out_dir, sprintf("%s_%02d.rds", part, k))
  if (file.exists(f)) return(invisible())
  t0 <- Sys.time()
  sc <- scenarios[k, ]
  sigs <- signatures(k)[seq_len(reps)]
  res <- parLapply(cl, sigs, fun, sc = sc)
  saveRDS(res, f)
  cat(format(Sys.time()), part, sc$label, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      "min\n", file = log_file, append = TRUE)
}
for (k in scenarios$id) run_part("grouping", k, function(sig, sc) one_grouping(sig), 200L)
for (k in scenarios$id) run_part("bases", k, function(sig, sc) one_bases(sc, sig), n_basis_rep)
for (k in scenarios$id) run_part("stability", k, function(sig, sc) one_stability(sc, sig), n_stab_rep)
for (k in scenarios$id) run_part("stability_mclust", k, function(sig, sc) one_stability(sc, sig, "mclust"), n_stab_rep)
stopCluster(cl)

# --- Summaries -----------------------------------------------------------------------------
read_part <- function(part) lapply(scenarios$id, function(k) readRDS(file.path(out_dir, sprintf("%s_%02d.rds", part, k))))
lab <- function(d) { d$label <- scenarios$label[d$scenario]; d$family <- scenarios$family[d$scenario]; d }

grp <- lab(do.call(rbind, unlist(read_part("grouping"), recursive = FALSE)))
saveRDS(grp, file.path("results", "grouping_replicates.rds"))
grp_sum <- do.call(rbind, lapply(split(grp, list(grp$scenario, grp$method), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], method = d$method[1],
             ari_median = stats::median(d$ari), ari_q05 = stats::quantile(d$ari, 0.05), ari_q95 = stats::quantile(d$ari, 0.95),
             ari_dependent_median = stats::median(d$ari_dependent, na.rm = TRUE),
             true_groups = d$true_groups[1], n_groups_median = stats::median(d$n_groups),
             n_groups_exact = mean(d$n_groups == d$true_groups), none_recall = mean(d$none_recall, na.rm = TRUE),
             none_precision = mean(d$none_precision, na.rm = TRUE), seconds = mean(d$seconds))))
grp_sum <- grp_sum[order(grp_sum$scenario, grp_sum$method), ]
write.csv(grp_sum, file.path("results", "grouping_summary.csv"), row.names = FALSE)

stab_raw <- c(unlist(read_part("stability"), recursive = FALSE), unlist(read_part("stability_mclust"), recursive = FALSE))
stab <- lab(do.call(rbind, lapply(stab_raw, `[[`, "table")))
saveRDS(list(table = stab, coassociation = lapply(stab_raw, `[[`, "coassociation")), file.path("results", "stability_replicates.rds"))
dep <- stab[stab$group != "none" & !is.na(stab$stability), ]
stab_sum <- do.call(rbind, lapply(split(dep, list(dep$scenario, dep$method), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], method = d$method[1],
             stability_median = stats::median(d$stability), ambiguous = mean(d$ambiguous), wrong = mean(d$wrong),
             wrong_if_ambiguous = if (any(d$ambiguous)) mean(d$wrong[d$ambiguous]) else NA,
             wrong_if_stable = if (any(!d$ambiguous)) mean(d$wrong[!d$ambiguous]) else NA)))
stab_sum <- stab_sum[order(stab_sum$scenario, stab_sum$method), ]
write.csv(stab_sum, file.path("results", "stability_summary.csv"), row.names = FALSE)

bas_raw <- unlist(read_part("bases"), recursive = FALSE)
bas <- lab(do.call(rbind, lapply(bas_raw, `[[`, "bases")))
uni <- lab(do.call(rbind, lapply(bas_raw, `[[`, "unions")))
saveRDS(list(bases = bas, unions = uni), file.path("results", "bases_replicates.rds"))
bas_sum <- do.call(rbind, lapply(split(bas, list(bas$scenario, bas$grouping, bas$group_type), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], grouping = d$grouping[1],
             group_type = d$group_type[1], n_bases = nrow(d), k_median = stats::median(d$k),
             kn_median = stats::median(d$kn), kn_q95 = stats::quantile(d$kn, 0.95),
             capture_median = stats::median(d$capture, na.rm = TRUE), absorption_median = stats::median(d$absorption, na.rm = TRUE),
             seconds_median = stats::median(d$seconds))))
bas_sum <- bas_sum[order(bas_sum$scenario, bas_sum$grouping, bas_sum$group_type), ]
write.csv(bas_sum, file.path("results", "bases_summary.csv"), row.names = FALSE)
uni_sum <- do.call(rbind, lapply(split(uni, list(uni$scenario, uni$grouping), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], label = d$label[1], family = d$family[1], grouping = d$grouping[1],
             n_unions = nrow(d), kn_median = stats::median(d$kn), kn_q95 = stats::quantile(d$kn, 0.95),
             above_0.1 = mean(d$kn > 0.1), above_0.2 = mean(d$kn > 0.2), above_0.3 = mean(d$kn > 0.3))))
write.csv(uni_sum, file.path("results", "union_summary.csv"), row.names = FALSE)
cat(format(Sys.time()), "summaries written\n", file = log_file, append = TRUE)
print(grp_sum, digits = 2)
print(stab_sum, digits = 2)
print(bas_sum, digits = 2)
