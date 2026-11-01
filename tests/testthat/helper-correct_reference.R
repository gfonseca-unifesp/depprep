# Small temporal and spatio-temporal corrections compared, in test-correct.R, with the
# outputs that the correct() of commit 1fe741c gave (fixtures/correct_reference.rds), so
# that new options left at their defaults keep every output the same.
correct_reference_cases <- function(...) {
  gt <- sim_groups(sim_group("ar", 3, "ar1", phi = 0.9), sim_group("seas", 3, "seasonal", period = 12),
                   sim_group("none", 1, "none"))
  st <- sim_temporal(120, gt, seed = 5)
  grt <- data.frame(variable = colnames(st$X), group = rep(c("g1", "g2", "none"), c(3, 3, 1)),
                    group_type = rep(c("ar", "seasonal", "none"), c(3, 3, 1)))
  class(grt) <- c("depprep_grouping", "data.frame"); attr(grt, "period") <- c(g1 = NA, g2 = 12)
  gs <- sim_groups(sim_group("sep", 3, "separable_ar1", range = 0.5, phi = 0.9), sim_group("none", 2, "none"))
  ss <- sim_spatiotemporal(20, 16, gs, seed = 8)
  grs <- data.frame(variable = colnames(ss$X), group = rep(c("g1", "none"), c(3, 2)),
                    group_type = rep(c("separable_ar1", "none"), c(3, 2)))
  terms <- function(b) lapply(b, attr, "terms")
  keep <- function(cr) list(R = cr$R, X = cr$X, bases = terms(cr$bases), variable_bases = terms(cr$variable_bases),
                            kn = cr$kn, union_kn = cr$union_kn, warnings = cr$warnings)
  list(temporal_conditional = keep(correct(st, grt, ...)),
       temporal_whiten = keep(correct(st, grt, mode = "whiten", ...)),
       st_whiten = keep(correct(ss, grs, ...)),
       st_conditional = keep(correct(ss, grs, mode = "conditional", data_basis = "union", ...)))
}
