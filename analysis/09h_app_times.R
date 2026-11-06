# Stage 9 — times of what the app runs on real data, for the messages of steps 4 and 5
# (decision D30; protocol-v1.3.1, section 17.5). Measures, on the shelf and on the slope of
# the Foz: correct() with the default of space-time, the other mode with the bases already
# built, and verify() as the app calls it. Also checks that R is the one of the registered
# runs. No random numbers beyond the seed of 09d.
#
# Requires results/fza/ocean_cube.csv and the corrections of 09g.
# Output (not tracked): results/fza/app_times.csv

run_region <- function(region) {
  Sys.setenv(DEPPREP_OCEAN_REGION = region)
  on.exit(Sys.unsetenv("DEPPREP_OCEAN_REGION"))
  env <- new.env(parent = globalenv())
  # the data and the grouping of 09d, the part of that script before the corrections
  src <- parse("analysis/09d_fza_ocean_analysis.R", keep.source = TRUE)
  stop_at <- grep("^# --- Corrections and verification", readLines("analysis/09d_fza_ocean_analysis.R", warn = FALSE))
  for (e in src[utils::getSrcLocation(src, "line") < stop_at]) eval(e, envir = env)
  X <- env$data_of[["0.25"]]; gr <- env$groupings[["0.25"]]
  pos <- list(coords = env$coords, index = env$index)
  cores <- env$params$cores
  out <- list()
  timed <- function(step, expr) {
    t0 <- Sys.time(); value <- expr
    out[[step]] <<- data.frame(design = sub(" .*", "", region), rows = nrow(X), variables = ncol(X), groups = length(unique(gr$group[gr$group != "none"])),
                               cores = cores, step = step, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    value
  }
  first <- timed("correct_default", do.call(correct, c(list(X, gr, mode = "conditional_ar", scope = "groups", seed = env$seed, cores = cores), pos)))
  other <- timed("correct_other_mode", do.call(correct, c(list(X, gr, bases = first$bases, mode = "whiten", scope = "groups", data_basis = "union",
                                                               seed = env$seed, cores = cores), pos)))
  saved <- function(run) readRDS(env$out_file(paste0("correction_", run, ".rds")))$correction$R
  stopifnot(max(abs(first$R - saved("conditional_ar"))) < 1e-10, max(abs(other$R - saved("main"))) < 1e-10)
  timed("verify_default", do.call(verify, c(list(first, X, nperm = 99), pos)))
  timed("verify_whiten", do.call(verify, c(list(other, X, nperm = 99), pos)))
  do.call(rbind, out)
}

times <- do.call(rbind, lapply(c("shelf", "slope and deep"), run_region))
utils::write.csv(times, file.path("results", "fza", "app_times.csv"), row.names = FALSE)
print(times, digits = 4, row.names = FALSE)
