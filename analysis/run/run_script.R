# Run one analysis script and record its provenance: the commit of the repository and the
# tracked files changed since it, the seeds written in the script, the start and end
# times, the exit status, the outputs the run wrote in results/ and the R session (R
# version, platform, the hash of renv.lock and the versions of the packages the package
# imports). One record per script, results/provenance/<script>.json (not tracked, like
# every output; a later run of the same script replaces it).
#
# Usage, from the root of the repository:
#   Rscript analysis/run/run_script.R analysis/NN_name.R [arguments of the script]
# The script runs in its own R process, exactly as `Rscript analysis/NN_name.R`.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !file.exists(args[1])) stop("Give the path of an analysis script.", call. = FALSE)
script <- args[1]
extra <- args[-1]
git <- function(...) {
  out <- tryCatch(suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (!is.null(attr(out, "status"))) character() else out
}
sha256 <- function(f) if (file.exists(f)) digest::digest(file = f, algo = "sha256") else NA_character_

src <- readLines(script, warn = FALSE)
seeds <- unique(unlist(regmatches(src, gregexpr("\\b(base_seed|seed)\\s*(<-|=)\\s*[0-9]+L?|set\\.seed\\([^)]*\\)", src))))
stamp <- function() list.files("results", recursive = TRUE, full.names = TRUE)
mtimes <- function(f) stats::setNames(file.info(f)$mtime, f)
before <- mtimes(stamp())
commit <- git("rev-parse", "HEAD")                       # at the start: commits made during a long run do not count
changed_tracked <- git("status", "--porcelain", "--untracked-files=no")
start <- Sys.time()
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(extra)))
end <- Sys.time()
after <- mtimes(stamp())
changed <- names(after)[is.na(before[names(after)]) | after > before[names(after)]]
outputs <- sort(changed[!grepl("^results/(provenance|logs)/", changed)])

desc <- read.dcf("DESCRIPTION", fields = "Imports")[1, 1]
imports <- trimws(sub("\\(.*", "", strsplit(desc, ",")[[1]]))
imports <- setdiff(imports[nzchar(imports)], c("graphics", "grDevices", "parallel", "stats", "utils"))
version_of <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)

record <- list(
  script = script, arguments = I(extra), script_sha256 = sha256(script),
  commit = commit,
  changed_tracked_files = I(changed_tracked),
  seeds_in_script = I(seeds),
  start = format(start, "%Y-%m-%d %H:%M:%S %z"), end = format(end, "%Y-%m-%d %H:%M:%S %z"),
  seconds = round(as.numeric(difftime(end, start, units = "secs"))), exit_status = status,
  outputs = I(outputs),
  environment = list(R = R.version.string, platform = R.version$platform, os = utils::sessionInfo()$running,
                     renv_lock_sha256 = sha256("renv.lock"),
                     packages = as.list(stats::setNames(vapply(c(imports, "devtools", "pkgload"), version_of, ""),
                                                        c(imports, "devtools", "pkgload")))))
dir.create(file.path("results", "provenance"), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(record, file.path("results", "provenance", paste0(sub("\\.R$", "", basename(script)), ".json")),
                     auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("%s: exit %d, %.1f min, %d outputs\n", basename(script), status, record$seconds / 60, length(outputs)))
quit(status = status)
