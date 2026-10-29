# Stage 8 — descriptive check of the basis menu (decision D16). Not a pre-registered test.
#
# Question. When the user knows the structure (the period, the direction of a current),
# does a menu with it build a better basis than the automatic choice, and does the choice
# of the neighbourhood graph (term_mem(W = "auto")) help against the default dbMEM?
# "Better" is read on the true components of the simulation: capture (share of the true
# large-scale structure removed) and absorption (share of the true residual removed beyond
# k/n), per variable, by true group (.capture_absorption()); and on the outputs: the
# error of R against the correlation of the true signal, and the k-means ARI on the data
# matrix (one basis per variable, D17).
#
# Scenarios (analysis/d4_scenarios.R), structured truth, 20 replicates:
# * S n=400 high (id 2): automatic; term_trend() + term_mem("auto"); term_trend() +
#   term_aem(30), the direction of the anisotropic group;
# * T T=480 high seasonal+trend (id 20): automatic; term_trend() + term_cycle(12), the
#   true period; term_trend() + term_cycle(), the period searched;
# * ST 50x60 (id 9): automatic; term_trend() + term_aem(0) + term_mem(margin = "time"),
#   the direction of the moving group (velocity along x).
# The menus go to the groups of the matching diagnosed type (see `variants`), in the
# conditional mode with one basis per group (scope "groups") for R; the data matrix keeps
# one basis per variable, built with the menu of its group.
# Seeds: 13200 + 1000 * scenario id + replicate.
#
# Outputs (not tracked): results/menu_check.csv, results/menu_check_summary.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

base_seed <- 13200L
n_rep <- 20L
bases_all <- rbind(scenarios, extra_scenarios, d14_scenarios)
# A menu with a known structure goes only to the groups of that diagnosed type (a cycle to
# the seasonal groups, a direction to the anisotropic or moving ones); the other groups keep
# the automatic menu. term_mem(W = "auto") goes to every group.
for_types <- function(types, m) function(gr) {
  g <- unique(gr$group[gr$group_type %in% types])
  if (!length(g)) NULL else stats::setNames(rep(list(m), length(g)), g)
}
seasonal_types <- c("seasonal", "separable_seasonal")
variants <- list(
  "2" = list(automatic = function(gr) NULL, mem_auto = function(gr) basis_menu(term_trend(), term_mem(W = "auto")),
             aem_30 = for_types("anisotropic", basis_menu(term_trend(), term_aem(30)))),
  "20" = list(automatic = function(gr) NULL, cycle_12 = for_types(seasonal_types, basis_menu(term_trend(), term_cycle(12))),
              cycle_searched = for_types(seasonal_types, basis_menu(term_trend(), term_cycle()))),
  "9" = list(automatic = function(gr) NULL,
             aem_0 = for_types("moving", basis_menu(term_trend(), term_aem(0), term_mem(margin = "time")))))

one_check <- function(id, r) {
  sc <- bases_all[bases_all$id == as.integer(id), ]
  seed <- base_seed + 1000L * as.integer(id) + r
  s <- simulate(sc, seed, "structured")
  gr <- group_regimes(diagnose_variables(s))
  n <- nrow(s$X)
  latent <- if (!is.null(s$X_latent)) s$X_latent else s$X
  share <- truth_table(s$groups)$share
  signal <- latent - sweep(s$components$large, 2, sqrt(share), "*")
  position <- switch(s$family, spatial = s$coords, temporal = matrix(s$time),
                     spatiotemporal = cbind(s$coords[s$index$site, ], s$index$time / max(s$index$time)))
  truth <- list(R_true = s$R_true, classes = s$classes, signal = signal, position = position)
  R_sig <- stats::cor(signal); ut <- upper.tri(R_sig)
  dep <- s$group != "none"
  do.call(rbind, lapply(names(variants[[id]]), function(v) {
    m <- variants[[id]][[v]](gr)
    targeted <- if (is.null(m)) 0L else if (inherits(m, "depprep_menu")) sum(unique(gr$group) != "none") else length(m)
    cr <- tryCatch(correct(s, gr, mode = "conditional", scope = "groups", seed = seed, menu = m), error = function(e) e)
    if (inherits(cr, "error")) return(data.frame(scenario = as.integer(id), label = sc$label, variant = v, replicate = r, seed = seed,
      true_group = NA, capture = NA, absorption = NA, kn = NA, own_capture = NA, own_absorption = NA, rms_signal = NA,
      kmeans_ari = NA, targeted = targeted, error = conditionMessage(cr)))
    ca <- t(vapply(which(dep), function(j) {
      g <- cr$group[j]
      b <- if (g == "none") matrix(0, n, 0) else cr$bases[[g]]
      ob <- cr$variable_bases[[colnames(s$X)[j]]]
      if (is.null(ob)) ob <- matrix(0, n, 0)
      c(.capture_absorption(b, s$components$large[, j], s$components$residual[, j]), kn = ncol(b) / n,
        .capture_absorption(ob, s$components$large[, j], s$components$residual[, j]))
    }, numeric(5)))
    colnames(ca) <- c("capture", "absorption", "kn", "own_capture", "own_absorption")
    km <- run_downstream(cr$X, "kmeans", truth = truth, seed = seed)
    by_g <- stats::aggregate(as.data.frame(ca), list(true_group = s$group[dep]), mean)
    cbind(data.frame(scenario = as.integer(id), label = sc$label, variant = v, replicate = r, seed = seed), by_g,
          rms_signal = sqrt(mean((cr$R - R_sig)[ut]^2)), kmeans_ari = km$value[km$metric == "ari_true_k"], targeted = targeted,
          error = NA_character_)
  }))
}

if (!identical(Sys.getenv("DEPPREP_SOURCE_ONLY"), "1")) {
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
  cl <- makeCluster(max(1, min(24, detectCores() - 4)), setup_timeout = 900)
  clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
    source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
  clusterExport(cl, c("base_seed", "bases_all", "variants", "for_types", "seasonal_types", "one_check"))
  jobs <- expand.grid(id = names(variants), r = seq_len(n_rep), stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$id != "9"), ]                                            # the slow ones first
  res <- parLapply(cl, seq_len(nrow(jobs)), function(k, jobs) one_check(jobs$id[k], jobs$r[k]), jobs = jobs)
  stopCluster(cl)
  res <- do.call(rbind, res)
  dir.create("results", showWarnings = FALSE)
  write.csv(res, file.path("results", "menu_check.csv"), row.names = FALSE)
  ok <- res[is.na(res$error), ]
  summ <- stats::aggregate(cbind(capture, absorption, kn, own_capture, own_absorption, rms_signal, kmeans_ari) ~ label + variant + true_group,
                           data = ok, FUN = mean)
  summ <- merge(summ, stats::aggregate(targeted ~ label + variant, data = res, FUN = mean))
  write.csv(summ, file.path("results", "menu_check_summary.csv"), row.names = FALSE)
  print(summ, digits = 3, row.names = FALSE)
  cat("errors:", sum(!is.na(res$error)), "\n")
}
