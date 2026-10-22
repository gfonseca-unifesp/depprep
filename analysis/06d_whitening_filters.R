# D14, block D2 — Do different whitening filters attenuate correlations across groups?
#
# Legendre & Legendre (2012, 12.2): every filter distorts the series. In whiten mode
# each group gets its own filter, and a derivation suggests that two different AR
# filters shrink the correlation of a pair across groups. Truth with correlation 0.3
# between groups (pilot scenarios added by D14); for every pair across groups, the
# correlation with each group's own filter (as correct() does) is compared with a
# common filter (both members filtered with one group's filter, averaged over the two).
# Pairs within groups are the reference for the bias that whitening has anyway (its
# estimand is the correlation of the innovations).
#
# Seeds: 7800 + 1000 * scenario + replicate (recorded in the output).
# Output (not tracked): results/whitening_filters.csv

suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("analysis", "d4_scenarios.R"))
one <- function(lab, r) {
  sc <- pilot_scenarios[pilot_scenarios$label == lab, ]
  seed <- 7800L + 1000L * sc$id + r
  s <- simulate(sc, seed, sc$truth)
  gr <- group_regimes(diagnose_variables(s, nsim = 99))
  cr <- correct(s, gr, mode = "whiten", seed = seed)
  grp <- cr$group; X <- s$X; colnames(X) <- colnames(cr$R)
  fam <- s$family
  gt <- vapply(setNames(unique(grp), unique(grp)), function(g) if (g == "none") "none" else gr$group_type[match(g, gr$group)], "")
  W <- lapply(setNames(unique(grp), unique(grp)), function(g) depprep:::.whitener(X[, grp == g, drop = FALSE], g,
          if (g == "none") "none" else gt[[g]], cr$bases[[g]], fam, s$coords, s$time, s$index))
  R_common <- cr$R_pairwise
  ug <- unique(grp)
  for (a in ug) for (b in ug) if (a < b) {
    ia <- which(grp == a); ib <- which(grp == b)
    f <- function(g) { Ea <- W[[g]](X[, ia, drop = FALSE]); Eb <- W[[g]](X[, ib, drop = FALSE]); m <- min(nrow(Ea), nrow(Eb)); cor(tail(Ea, m), tail(Eb, m)) }
    rc <- (f(a) + f(b)) / 2
    R_common[ia, ib] <- rc; R_common[ib, ia] <- t(rc)
  }
  tg <- s$group; dep_true <- tg != "none"
  across <- outer(tg, tg, "!=") & outer(dep_true, dep_true) & upper.tri(diag(length(tg)))
  within <- outer(tg, tg, "==") & outer(dep_true, dep_true) & upper.tri(diag(length(tg)))
  rt <- s$R_true
  data.frame(label = lab, r = r, across_own = mean(cr$R_pairwise[across] - rt[across]), across_common = mean(R_common[across] - rt[across]),
             within_own = mean(cr$R_pairwise[within] - rt[within]), truth_across = mean(rt[across]))
}
labs <- c("T T=480 high structured between 0.3", "ST 50x60 structured between 0.3")
res <- do.call(rbind, c(lapply(1:20, function(r) one(labs[1], r)), lapply(1:5, function(r) one(labs[2], r))))
res$seed <- 7800L + 1000L * pilot_scenarios$id[match(res$label, pilot_scenarios$label)] + res$r
write.csv(res, file.path("results", "whitening_filters.csv"), row.names = FALSE)
print(aggregate(cbind(across_own, across_common, within_own, truth_across) ~ label, res, mean), digits = 3)
