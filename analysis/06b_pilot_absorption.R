# Stage 6 — Why the shared bases lower class recovery in the pilot.
#
# In the structured pilot scenarios (3 classes at random positions) the global basis and
# the union of the group bases lowered the ARI of k-means and SOM, and one basis per
# variable did not. Random classes add the same pattern to every variable, so a basis
# selected on all variables at once may pick the eigenfunctions that best follow the
# class pattern. This script measures, for each basis, the share of the class space
# (the centred class indicators) that lies in the basis, against k/n, the share expected
# for a basis unrelated to the classes.
#
# Seeds: the pilot's (6600 + 1000 * scenario + replicate), replicates 1-10.
# Output (not tracked): results/pilot_absorption.csv

suppressMessages(devtools::load_all(quiet = TRUE))
library(parallel)
source(file.path("analysis", "d4_scenarios.R"))

one_absorption <- function(sc, r) {
  seed <- 6600L + 1000L * sc$id + r
  s <- simulate(sc, seed, sc$truth)
  n <- nrow(s$X)
  C <- scale(stats::model.matrix(~ factor(s$classes) - 1), scale = FALSE)
  Qc <- qr.Q(qr(C))[, seq_len(length(unique(s$classes)) - 1), drop = FALSE]
  class_share <- function(B) if (ncol(B) == 0) 0 else sum(crossprod(qr.Q(qr(B)), Qc)^2) / ncol(Qc)
  union_of <- function(cobj) {
    b <- Filter(function(x) !is.null(x) && ncol(x) > 0, cobj$bases)
    if (!length(b)) matrix(0, n, 0) else do.call(union_basis, c(b, kn_max = 1))
  }
  sig <- diagnose_variables(s)
  type <- unname(classify_regime(sig))
  gr <- group_regimes(sig)
  args <- switch(s$family, spatial = list(coords = s$coords), temporal = list(time = s$time))
  Bd <- union_of(correct(s, gr, seed = seed))
  Bg <- union_of(correct(s, gr, mode = "conditional", scope = "groups", seed = seed))
  ind <- t(vapply(seq_len(ncol(s$X)), function(j) {
    if (type[j] == "none") return(c(0, 0))
    B <- do.call(build_group_basis, c(list(X_group = s$X[, j, drop = FALSE], type = type[j],
      period = if (is.finite(sig$period[j])) round(sig$period[j]) else 12, seed = seed + j), args))
    c(ncol(B) / n, class_share(B))
  }, numeric(2)))
  data.frame(scenario = sc$id, label = sc$label, family = s$family, replicate = r, seed = seed,
             approach = c("default", "groups", "independent"),
             kn = c(ncol(Bd) / n, ncol(Bg) / n, mean(ind[, 1])),
             class_share = c(class_share(Bd), class_share(Bg), mean(ind[, 2])))
}

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
cl <- makeCluster(max(1, min(20, detectCores() - 4)), setup_timeout = 900)
clusterCall(cl, function(wd) { setwd(wd); suppressMessages(devtools::load_all(quiet = TRUE))
  source(file.path("analysis", "d4_scenarios.R")); NULL }, getwd())
clusterExport(cl, "one_absorption")
todo <- pilot_scenarios[pilot_scenarios$truth == "structured" & pilot_scenarios$family != "spatiotemporal", ]
jobs <- expand.grid(i = seq_len(nrow(todo)), r = 1:10)
res <- parLapply(cl, seq_len(nrow(jobs)), function(j, jobs, todo) one_absorption(todo[jobs$i[j], ], jobs$r[j]),
                 jobs = jobs, todo = todo)
stopCluster(cl)
res <- do.call(rbind, res)
res$excess <- res$class_share / res$kn
write.csv(res, file.path("results", "pilot_absorption.csv"), row.names = FALSE)
print(aggregate(cbind(kn, class_share, excess) ~ label + approach, data = res, FUN = mean), digits = 3)
