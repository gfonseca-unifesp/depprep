# analysis/

Scripts numbered by project stage (`01_...`, `02_...`). Each script calls `set.seed()`
explicitly, records the seed in its output and writes its results to `results/`
(figures to `results/figs/`).

## Running

- `analysis/run/run_all.sh` runs the scripts in order, from the table
  `analysis/run/steps.tsv` (step, script, arguments, the steps it needs, the network or
  login it needs, the time it took). `bash analysis/run/run_all.sh 05a 06g` runs a range,
  `DRY_RUN=1` lists the steps, `SKIP_NETWORK=1` leaves out the downloads, and `RSCRIPT`
  sets the R to use.
- Each script runs through `analysis/run/run_script.R`, which records its provenance in
  `results/provenance/<script>.json`: the commit and the tracked files changed since it,
  the seeds written in the script, the time, the outputs and the R session (R, platform,
  the hash of `renv.lock`, the versions of the imported packages). A script run directly
  with `Rscript` records nothing.
- `analysis/run/history/`: the shell scripts with which earlier stages were rerun, kept as
  a record.

## Raw data

`data/raw/` is not tracked. `analysis/00_raw_manifest.R` writes `data/raw_manifest.csv`
(tracked): size, SHA-256, date written, source, product and version of every file.
`Rscript analysis/00_raw_manifest.R check` compares the files present with it.
