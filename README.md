# depprep

**Dependence-aware preprocessing for multivariate exploratory analysis.**

> **Status:** stages 0 to 8 are complete: simulation generators, the diagnosis of each
> variable, grouping and shared bases, the correction, a pilot, the full pre-registered
> simulation study (protocol-v1, with two labelled amendments) and the package with its
> menu of structures, report and app. The real-data application (stage 9, analysis plan
> registered as protocol-v1.3) is in progress. Missing values are accepted, without
> imputation: each variable uses the rows where it was observed (section 4).

---

## 1. The problem, in plain words

Exploratory methods such as **PCA, MDS, k-means and self-organizing maps (SOM)** assume
that observations are independent. Environmental data rarely are: samples that are close
in space or in time tend to be similar (*autocorrelation*).

The trouble gets worse when a data set has **dozens of variables with different
dependence regimes**: one variable varies over long distances, another over short ones,
another follows a seasonal cycle, another has no structure at all. Two things then go
wrong:

1. **Apparent association.** Variables that follow the same large-scale structure (for
   example, the same coast-to-offshore gradient) look correlated even when they have no
   relation beyond that shared structure. An exploratory method will happily "discover"
   that association.
2. **Distortion by different corrections.** The usual fix is to correct each variable on
   its own (its own spatial filter, its own AR model). But correcting two variables with
   *different* bases changes the correlation between them in a way that corresponds to no
   real structure.

## 2. The idea

Each variable is seen as two parts:

```
 variable  =  large-scale component of its group  +  residual component
              (what the correction removes)          (what the correction keeps:
                                                      the correlations between residuals
                                                      are the target)
```

The rule that makes the correction coherent:

- **Same group, same basis.** Variables with a similar dependence regime form a group and
  are all corrected with the *same* spatial/temporal basis.
- **Pairs across groups** use the union of the two groups' bases.
- **Variables without dependence** are not corrected at all.
- Groups are formed from each variable's **own** dependence signature (variogram, ACF),
  never from the correlations that are estimated afterwards.

The planned workflow:

```
 diagnose each variable  ->  group similar regimes  ->  build one basis per group
        ->  correct  ->  verify (what was removed, what was kept)  ->  corrected matrix
                                                                     for PCA, MDS,
                                                                     k-means, SOM
```

## 3. What is in this repository

This repository holds **code only**: the R package, its tests and the script of each
project stage. Figures and tables are *not* versioned; every script regenerates them
from recorded seeds into `results/` on your machine.

| Path | Content |
| --- | --- |
| `R/` | Package functions: simulation generators, diagnosis, grouping, bases and their menu, correction, verification, report and app |
| `tests/testthat/` | Tests with known truth for every function |
| `man/` | Function documentation (generated with roxygen2) |
| `inst/extdata/` | Example tables of the Foz do Amazonas case study (section 4) |
| `analysis/` | One script per project stage, numbered in the order they are run |
| `results/` | Created by the scripts; not tracked by git |
| `data/raw/` | Downloaded raw data; not tracked by git, never modified |
| `renv.lock` | Exact versions of every R package used |

## 4. The app: the workflow without code

`run_app()` opens a point-and-click interface (Shiny) in the browser. It runs the same
functions as the code in sections 9 to 12 and, at the end, gives you the R code that
reproduces every choice you made.

**Start it from R or RStudio:**

```r
install.packages(c("remotes", "shiny"))
remotes::install_github("gfonseca-unifesp/depprep")   # once
depprep::run_app()
```

**Or from the command line** (a terminal, with R on the PATH):

```bash
Rscript -e "depprep::run_app(launch.browser = TRUE)"
```

The app opens in your default browser; close the terminal (or press Ctrl+C) to stop it.
To try it without your own data, click **Use example data** in step 1 (a spatial
simulation with three groups of variables and three hidden classes).

### The input table

One CSV file, **one row per observation**, with a header row. Some columns give the
position of the observation; every other numeric column is a variable. Column names are
free (you mark the position columns in step 1). The three layouts (the values below are
illustrative):

**Spatial** (one map): one row per site, with its longitude and latitude in degrees.

| lon | lat | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- | --- |
| -48.125 | 1.375 | 1.84 | 28.31 | 0.21 | 210.4 |
| -47.875 | 1.375 | 1.12 | 28.45 | 0.15 | 168.9 |
| -48.125 | 1.625 | 2.07 | 28.22 | 0.24 | 231.0 |

**Temporal** (one series): one row per time, with a numeric time.

| time | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- |
| 1 | 1.52 | 28.10 | 0.19 | 195.2 |
| 2 | 1.61 | 28.34 | 0.20 | 201.7 |
| 3 | 1.95 | 28.62 | 0.23 | 224.8 |

**Space-time** (a series at every site): one row per site and time, with a site
identifier, the site's coordinates (repeated on each of its rows) and the time.

| site | lon | lat | time | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | -48.125 | 1.375 | 1 | 1.84 | 28.31 | 0.21 | 210.4 |
| 1 | -48.125 | 1.375 | 2 | 1.90 | 28.40 | 0.22 | 214.1 |
| 2 | -47.875 | 1.375 | 1 | 1.12 | 28.45 | 0.15 | 168.9 |
| 2 | -47.875 | 1.375 | 2 | 1.20 | 28.51 | 0.16 | 172.3 |

Rules:

- **At least two numeric variables**, each with at least three observed values.
- **Missing values are allowed in the variables** (empty cells, or `NA` in code), not in
  the positions: the coordinates, the time and the site must be complete. Nothing is
  imputed; each variable uses the rows where it was observed (see "Missing values"
  below). Infinite values stop with an error.
- **Time is a number**: a month counter (1, 2, 3, ...) or a decimal year (2015.25), not a
  date text such as "2015-04". Equal steps are expected.
- **Coordinates as longitude and latitude in degrees.** They treat a degree of
  longitude as a degree of latitude, while it is shorter by a factor of cos(latitude):
  only 1.2% at 9 degrees (the north of the Foz do Amazonas basin). The app shows that
  difference; beyond about 11 degrees of latitude (more than 2%) project the coordinates
  first (for example to UTM, in km) and use them in place of lon and lat. The app leaves
  columns named like coordinates (`lon`, `lat`, `x`, `y`, ...) out of the variables it
  proposes.
- **Space-time needs a complete grid**: every site has a row at every time, once. A
  value that was not observed stays as an empty cell in its row; the row itself is kept.
- Text columns (names, dates as text, labels) are ignored as variables; you can keep them
  in the file.

The same layouts work in code: `diagnose_variables()` takes the matrix of variables with
`coords` (spatial), `time` (temporal), or `coords` of the sites plus `index` (a data
frame with the `site` and `time` of each row, space-time).

### Missing values

Satellite fields have gaps (clouds, sun glint, the edge of a swath), and filling them
would add a dependence that belongs to the interpolation. So nothing is imputed: each
variable uses the rows where it was observed (decision D26).

| Step | What happens with the missing values of a variable |
| --- | --- |
| Diagnosis (`diagnose_variables()`) | In space, the statistics use the sites where the variable was observed (the nearest neighbours among them). In time, the series keeps its time steps with the gaps in place: the trend, the cycle and the periodogram are fitted on the observed times, and the autocorrelation and the AR fit skip the missing ones. In space-time, the time slices and the sites summarised are those with at least half of their values observed |
| Basis (`build_group_basis()`) | The terms of a group are selected with each missing value counting as the variable's mean, so it adds nothing to what a candidate explains, and a variable with many gaps weighs less in the group. The selected terms are defined at every row |
| Correction (`correct()`) | Each variable is regressed on its group's basis (and the constant) over its own observed rows. In the whitening mode, each time is whitened at its observed sites, and the AR filter gives a missing value where the value or one of its lags is missing. The corrected data keep the missing values in place; the correlations of `R` are computed on the rows where both variables have a value |
| Verification (`verify()`) | Each variable is described on its observed rows: the autocorrelation, the Moran correlograms and the local Moran's I at its observed sites (in space-time, the site means over the observed times), and the Ljung-Box test with the gaps in place. The covariances behind the Q3 warning use the rows where both variables have a value; the Mantel correlogram uses the rows where every variable has a value |
| Variation partition (`variation_partition()`) | Each variable on its observed rows |
| App (`run_app()`) | Empty cells are accepted. Step 1 notes the share of missing values and the variables that have them. The previews of the techniques (PCA, k-means) use only the rows where every variable has a value, raw and corrected, and the periodogram plot is left out for a series with gaps (its period is in the table of the diagnosis) |

What to keep in mind:

- **The shared basis holds only approximately.** The variables of a group share the same
  basis, but each is projected on its own observed rows. When their gaps differ, their
  projections differ a little, and the rule "same basis, same projection" is then only
  approximate. With complete data every output is the same as before.
- **Errors you may see.** A variable with fewer than three observed values; a pair of
  variables with too few rows observed together to estimate their correlation; missing
  values in the positions; and the scale mode (`mode = "scale"`), which still needs
  complete data. The drivers of `correct(covariates = )` cannot have missing values
  either: take the lags from the series before the first time.
- **What the gaps cost** was measured on simulated data with 10% to 30% of the values
  missing, at random and in patches like clouds (in space-time, a region during one
  season): `analysis/09c_missing_check.R`.

### Example with real data: the Foz do Amazonas ocean

Two tables of the project's case study come with the package, in the layouts above.
They are built from monthly Level-3 satellite fields of the marine Foz do Amazonas basin
(0.25-degree cells, April 2015 to October 2021) by `analysis/08d_fza_example_tables.R`:

| File | Layout | Rows | Columns |
| --- | --- | --- | --- |
| `foz_ocean_map.csv` | spatial | 544 cells | `lon`, `lat` (cell centres, degrees), then the mean of each variable over its cloud-free months |
| `foz_ocean_series.csv` | temporal | 79 months | `month` (text, ignored as a variable), `time` (1 to 79), then the basin mean of each variable in that month |
| `foz_ocean_sources.csv` | — | 5 variables | name, unit, ERDDAP dataset, provider and licence of each variable |

The variables: `chl` (chlorophyll-a, mg m-3), `par` (photosynthetically available
radiation, einstein m-2 day-1), `poc` (particulate organic carbon, mg m-3), `sst` (sea
surface temperature by day, degree C) from Aqua MODIS (NASA GSFC OBPG), and `wind_x`
(eastward wind, m s-1) from Metop ASCAT (NOAA ERD), all read from the NOAA CoastWatch
ERDDAP server. They are free to use and redistribute; please cite the providers.

They have no gaps because they are averages: cells enter the map when every variable is
valid in at least half of the months, and the basin mean of each month uses the valid
cells. They are examples for learning the tool, not the data of the case study (stage 9),
which will work on every cell and month.

In the app, download a file (for example
[foz_ocean_map.csv](https://raw.githubusercontent.com/gfonseca-unifesp/depprep/main/inst/extdata/foz_ocean_map.csv)
or [foz_ocean_series.csv](https://raw.githubusercontent.com/gfonseca-unifesp/depprep/main/inst/extdata/foz_ocean_series.csv))
and load it in step 1, marking `lon` and `lat` (spatial; the cells lie between the
equator and 9 N) or `time` (temporal). In code:

```r
library(depprep)
map <- read.csv(system.file("extdata", "foz_ocean_map.csv", package = "depprep"))
X   <- as.matrix(map[, c("chl", "par", "poc", "sst", "wind_x")])
xy  <- as.matrix(map[, c("lon", "lat")])
gr  <- group_regimes(diagnose_variables(X, coords = xy))
cr  <- correct(X, gr, coords = xy)                 # about half a minute
v   <- verify(cr, X, coords = xy)

series <- read.csv(system.file("extdata", "foz_ocean_series.csv", package = "depprep"))
Xt <- as.matrix(series[, c("chl", "par", "poc", "sst", "wind_x")])
crt <- correct(Xt, group_regimes(diagnose_variables(Xt, time = series$time)), time = series$time)
```

### The six steps

| Step | What you do | What you get |
| --- | --- | --- |
| 1. Data | Load a CSV with one row per observation and mark the position columns (x and y, time, site); choose spatial, temporal or space-time | A check of the columns and the first rows |
| 2. Context | Say which technique the data are for (PCA, or MDS, k-means, SOM and HDBSCAN) and whether the classes you seek could be regions or periods | Which corrected matrix the technique reads, and a warning when correcting would remove the classes you are looking for |
| 3. Diagnosis | Look at each variable's dependence and its plots (variogram, ACF, periodogram); change a regime by hand, or build the menu of structures of a group (trend, cycle with its period, Moran maps with a chosen graph, directional maps) | The groups of variables that will share a basis |
| 4. Correction | Accept the default (the data matrix with one basis per variable, the correlation matrix with the shared basis) or choose another mode (conditional, whitening) and other bases for the correlation matrix | The bases used, with their size (k, k/n) and first terms |
| 5. Verification | Check what was removed and what was kept | The autocorrelation of each variable before and after, the change in the correlations, the warning of signal at the scale of the basis, and a preview of the technique (PCA or k-means) before and after |
| 6. Export | Download the results | The corrected data (CSV), the corrected correlation matrix (CSV), the HTML report of `report_correction()` and the R code |

Notes and limits shown in the app: the share of missing values, when there are any
(empty cells are accepted; missing positions are not); coordinates in degrees should be
projected beyond about 11 degrees of latitude; and data with more than a few thousand
observations are slow (the correction of a 628-site x 93-month grid takes hours).

## 5. Known limits

### Why purely spatial regimes are hard to tell apart

In the simulations the diagnosis and the grouping work well in time and in space-time,
but much less well for purely spatial data with up to a few hundred sites. The reasons
are informative about real data, not only about the method:

1. **Different processes can look alike on average.** An anisotropic field (patches
   elongated in one direction, as a river plume along a coast) has a long range along
   its axis and a short one across it. Seen in all directions at once, its average scale
   is close to that of a short-range isotropic field. In the simulation grid the two
   regimes had nearly the same estimated range, nugget and Moran's I; only the
   anisotropy ratio, the noisiest part of the signature, set them apart.
2. **A map is a single realisation.** A time series of a few hundred steps repeats the
   same process many times, so its persistence or its seasonal cycle is estimated
   precisely. A map of a few hundred sites is one draw of a spatial field: neighbouring
   sites repeat each other's information, and the estimated range varies a lot from one
   map to another.
3. **Small-scale noise hides short-range structure.** Variation without spatial
   structure (measurement error, processes finer than the sampling) adds to the
   variogram at every distance. When the range is short, the structure lives at the
   few short distances where there are fewest pairs of sites, and it blends into the
   nugget.
4. **Sampling density sets what can be seen.** When the typical distance between
   neighbouring sites is close to the range of a process, the variogram has almost no
   points before it levels off, and the structure is read as "no dependence".

What this means in practice: a grouping error between two spatial regimes of similar
scale costs little, because they receive similar bases; the tool flags these variables
as ambiguous (unstable over subsamples) instead of hiding the problem; and data observed
repeatedly in time (space-time) carry much more information about each spatial regime
than a single map. Results: `analysis/03_*` and `analysis/04_*`.

### What the correction can and cannot do for the exploratory methods (stage 6 pilot)

- **It removes most spurious patterns.** Without correction, MDS configurations and
  SOM maps reproduce the positions of the samples, and k-means finds groups in data
  without classes (mostly in time). After correction these patterns almost disappear.
- **It cannot remove shared structure finer than its basis.** With rough spatial fields
  the first PCA axis stays above the broken-stick expectation after any basis
  correction: the fine part of the shared field is common to all variables, and no
  method that uses only the positions of the samples can tell it from a real
  association.
- **A large shared basis also removes part of the real signal.** When true classes are
  scattered among the samples, a basis selected on all variables at once takes away a
  share of every variable's class signal, roughly in proportion to its size. A small
  basis per variable, or whitening in time, keeps the classes.
- **Classes that are contiguous in space or time are removed by design.** There the
  classes are the large-scale structure. Keeping only the broad scales preserves them.
  So when the classes sought are regions or periods, the data should not be
  decorrelated: the reference tool is clustering constrained by spatial or temporal
  contiguity (Legendre & Legendre 2012, sections 12.6 and 13.3), available as
  `run_downstream(X, "constrained", truth = list(position = ...))`.

### Not supported yet

- **Missing positions and drivers.** Missing values in the variables are accepted
  (section 4, "Missing values"), but the coordinates, the times and the drivers must be
  complete. Nothing is imputed, by design: filling gaps would add a dependence that
  belongs to the interpolation. The scale mode of `correct()` still needs complete
  data.
- **Filters other than the autoregressive ones** (SAR/CAR whitening in space, STAR for
  moving patterns) are not implemented.
- **Coordinates in degrees** are used as they are, with a degree of longitude taken as a
  degree of latitude: fine near the equator, but beyond about 11 degrees of latitude
  they should be projected first (section 4).

## 6. Reproducibility

- **Seeds.** Every simulation uses an explicit seed, stored in its output. The generators
  use the L'Ecuyer-CMRG generator with one independent stream per component (coordinates,
  large-scale fields, residuals, classes), so changing one component leaves the others
  unchanged, and they restore your random number generator afterwards. In the stage
  scripts, replicate `r` of scenario `k` uses `base_seed + 1000 * k + r`.
- **Tests with known truth.** Each function is checked against the value it was asked to
  produce: variogram range within 25% (mean of 50 replicates), AR(1) coefficient within
  0.05, residual correlations within their sampling error, exact displacement of the
  moving pattern, and more.
- **Package versions.** `renv.lock` pins every package.
- **One table of every script.** Section 7 lists each analysis script with the steps whose
  results it reads, the network or login it needs and the time it took, so any script
  can be found, run and audited on its own.
- **Provenance of each run.** Run through `analysis/run/` (section 7), a script records in
  `results/provenance/<script>.json` the commit, the tracked files changed since it, the
  seeds written in the script, the time, its outputs and the R session (R, platform, the
  hash of `renv.lock`, the versions of the imported packages).
- **Raw data.** `data/raw/` is not tracked; `data/raw_manifest.csv` (tracked, written by
  `analysis/00_raw_manifest.R`) gives the size, SHA-256, date, source, product and
  version of every downloaded file. `Rscript analysis/00_raw_manifest.R check` compares
  your copy with it.

## 7. Project stages and their scripts

Stages 0 to 8 are done; stage 9 (the real-data application) is in progress. Every
analysis script is listed below, in the order it is run, so that each one can be found,
run and audited on its own. **Needs** gives the steps whose outputs in `results/` the
script reads: run those first. The same table, in tab-separated form, drives the runner
(`analysis/run/steps.tsv`).

```bash
bash analysis/run/run_all.sh 05a                # one step, recording its provenance
bash analysis/run/run_all.sh 05a 06g            # a range of steps, in the order of the table
DRY_RUN=1 bash analysis/run/run_all.sh          # list every step without running it
SKIP_NETWORK=1 bash analysis/run/run_all.sh     # leave out the steps that download
Rscript analysis/05a_h2_mechanism.R             # a script run directly records no provenance
```

Run from any folder; `RSCRIPT` sets the R to use. Each step writes its log to
`results/logs/<step>.log`; a failed step stops the run.

| Stage (package code) | Step | Script | What it produces | Needs | Access | Time |
| --- | --- | --- | --- | --- | --- | --- |
| **0. Environment**<br>`renv.lock` | 00 | `00_raw_manifest.R` | Manifest of `data/raw/` (size, SHA-256, date, source, product, version); `check` compares your copy with it | — | — | minutes |
| **1. Real data: is there a gap?** | 01a | `01a_fza_study_area.R` | Study area from the ANP basin limits and ETOPO1 topography | — | ANP WFS, NOAA ERDDAP | minutes |
| | 01b | `01b_fza_data_feasibility.R` | Catalogue of candidate L3 products; monthly cloud-free coverage of the products served anonymously by ERDDAP | 01a | NOAA ERDDAP | minutes |
| | 01c | `01c_fza_nasa_feasibility.R` | The same for NASA products (SMAP salinity and soil moisture, MODIS land and ocean colour, GPM IMERG), subset with OPeNDAP | 01a | Earthdata login (note 1) | under 1 h |
| **2. Simulation with known truth**<br>`R/sim_*.R` | 02 | `02_simulation_examples.R` | Schemes and galleries of the generators, variograms, space-time maps, apparent-association matrices, observation classes, time per replicate | — | — | minutes |
| | 02b | `02b_generator_recovery.R` | Recovery of the simulated range and AR coefficient (50 replicates per scenario); error of the residual and raw correlations | — | — | — |
| **3. Diagnosis of each variable**<br>`R/diagnose*.R`, `R/classify.R` | 03 | `03_diagnostics.R` | Signatures of every variable of the grid (200 replicates per scenario, parallel, resumable): accuracy, scale errors, bootstrap coverage | — | — | — |
| | 03f | `03_diagnostics_figures.R` | Figures: scheme, didactic variogram, typical signatures, confusion matrices, threshold curves | 03 | — | minutes |
| **4. Groups and shared bases**<br>`R/group.R`, `R/basis.R` | 04 | `04_grouping.R` | Grouping of the signatures (hclust, mclust), stability over subsamples, one basis per group: k/n, capture of the true large scale, absorption of the true residual (resumable) | 03 | — | — |
| | 04f | `04_grouping_figures.R` | Figures: scheme, Moran eigenvector maps, one replicate step by step, agreement, stability, k/n curves | 03 04 | — | minutes |
| **5. Correction and the mechanism**<br>`R/correct.R`, `R/verify.R` | 05a | `05a_h2_mechanism.R` | H2a: the error of different bases against their distance; H2b: the shared structure a basis leaves; calibration of the second grouping level | — | — | 30 min |
| | 05b | `05b_correction.R` | The correction on the grid: eight approaches, from none to the true grouping; error of the residual correlations, capture, degrees of freedom, time | 03 05a | — | 5 h |
| | 05c | `05c_whitening_st.R` | Why whitening helps in space-time, by kind of shared process | — | — | 10 min |
| | 05f | `05_correction_figures.R` | Figures of stage 5 | 05a 05b | — | minutes |
| **6. Pilot with the exploratory methods**<br>`R/downstream.R`, `R/comparators.R` | 06 | `06_pilot.R` | Pilot scenarios, 20 replicates each: every correction followed by PCA, MDS (Euclidean and Gower), k-means and SOM, scored with `run_downstream()`; SBSS and MAF (`decorrelate_latent()`) as comparators; paired contrasts and the replicates needed for stage 7 | — | — | 2.5 h |
| | 06f | `06_pilot_figures.R` | Figures of stage 6 | 06 | — | minutes |
| | 06b | `06b_pilot_absorption.R` | How much class signal each basis removes, against what its size alone implies | — | — | 5 min |
| | 06c | `06c_q3_warning.R` | Calibration of the first threshold of the Q3 warning | — | — | 5 min |
| | 06d | `06d_whitening_filters.R` | Whitening filters compared | — | — | 5 min |
| | 06h | `06f_hdbscan_pca.R` | HDBSCAN as a fifth method, and the metric of the PCA | — | — | — |
| | 06g | `06g_protocol_tests.R` | Tests on the pilot that fixed the thresholds of the protocol | 05b 06 | — | — |
| **7. Full pre-registered simulations** | 07 | `07_full_grid.R` | Grid of protocol-v1: 21 base scenarios x five truths (plus perturbations), 100 replicates in space and time, 50 in space-time; every approach and method (resumable, 16 processes) | — | — | 9 h |
| | 07b | `07b_h2a.R` | The between-basis error (H2a) on new seeds | 05a | — | — |
| | 07a | `07_analysis.R` | Pre-registered paired tests with Holm's correction, and the decision rule | 07 07b | — | minutes |
| | 07c | `07c_amendment.R` | Amendment protocol-v1.1: the correlation error against the true signal | 07 | — | — |
| | 07d1 | `07d_q3_recalibration.R calibrate` | Amendment protocol-v1.2: one Q3 threshold per data family | 07 | — | minutes |
| | 07d2 | `07d_q3_recalibration.R validate` | Validation of those thresholds on new seeds | 07d1 | — | — |
| | 07e | `07e_q3_cause.R` | Why the Q3 statistic depends on the family (exploratory) | 07 07d2 | — | — |
| **8. Package: menu, report, app**<br>`R/menu.R`, `R/report.R`, `R/app.R` | 08a | `08a_q3_reference_pilot.R` | Pilot of a Q3 statistic normalised by an independence reference that keeps each variable's autocorrelation | 06 07 | — | — |
| | 08b | `08b_menu_check.R` | Descriptive check of the menu of structures: known period, chosen graph, directional maps | — | — | — |
| | 08c | `08c_q3_size_check.R` | The Q3 thresholds at the size of the real case (93 months, 628 cells, 628 x 93) | 07d1 | — | — |
| | 08d | `08d_fza_example_tables.R` | The example tables of the Foz do Amazonas ocean in `inst/extdata/` | 01a 01b | — | minutes |
| | 08e | `08e_sample_sites_check.R` | Bases chosen on a sample of sites, against all sites | — | — | — |
| | 08f | `08f_correct_time.R` | Time of the correction on one and on several processes | — | — | — |
| **9. Real data: the Foz do Amazonas ocean** (in progress)<br>`R/partition.R` | 09a | `09a_fza_ocean_cube.R` | The ocean cube (cell x month) and the drivers; resumes from the files present | 01a | Earthdata, NOAA ERDDAP, ANA, NOAA PSL | hours (downloads), minutes (cached) |
| | 09b | `09b_fza_gap_census.R` | Census of the gaps of the cube | 01a 01b 09a | — | minutes |
| | 09c | `09c_missing_check.R` | What gaps cost the correction: 10% to 30% missing, at random and in cloud-like patches | — | — | 1-2 h |
| | 09d | `09d_fza_ocean_analysis.R` | The ocean case of protocol-v1.3: five runs and the output methods, raw against corrected (resumable) | 09a | — | 2-3 h |

Notes:

1. **Earthdata login.** A free [NASA Earthdata](https://urs.earthdata.nasa.gov) account,
   stored in `~/_netrc` (`machine urs.earthdata.nasa.gov login <user> password
   <password>`); the scripts never print it.
2. **Why Level-3 data.** The application uses L3 satellite time series (gridded, not
   gap-filled) of the Foz do Amazonas basin, Brazilian Equatorial Margin, in two case
   studies sharing one rectangle: the ocean (the marine part of the basin, official ANP
   limit) and the land (within 200 km of that coast). Gap-filled (L4) products add
   spatial autocorrelation that belongs to the interpolation, not to the process.
3. **Amendments.** In stage 7 the pre-registered results stay the main ones; the amended
   results are reported beside them and labelled.
4. **Shared scenarios.** `analysis/d4_scenarios.R` defines the simulation grid (groups,
   scenarios, truth) and is read by the stage scripts; it is not run on its own.
5. **Earlier reruns.** `analysis/run/history/` keeps the shell scripts with which earlier
   stages were rerun, as a record.

## 8. Getting started

```r
# R >= 4.3 (developed with R 4.6.1)
install.packages("renv")
renv::restore()          # installs the exact package versions in renv.lock
devtools::load_all()     # loads the package from source
devtools::test()         # runs the tests (a few take about a minute)
```

Run the analysis scripts from the repository root, alone or in order, as shown in
section 7. Outputs go to `results/`.

## 9. Quick start: simulating data with known truth

Every variable is simulated as

```
x_j = sqrt(s) * L_g  +  sqrt(1 - s) * e_j  +  class effect
```

where `L_g` is the large-scale component of the variable's group (unit variance), `e` is
the residual with a known correlation `R_true` across variables, and `s` is the share of
variance that comes from the large scale.

```r
library(depprep)

groups <- sim_groups(
  sim_group("long_a", 6, "exponential", range = 0.5),      # long-range field
  sim_group("long_b", 6, "exponential", shared_with = "long_a"),  # same field as long_a
  sim_group("short",  6, "exponential", range = 0.1),      # short-range field
  sim_group("none",   6, "none")                           # no spatial structure
)
R_true <- make_R_true(groups, "blocks", within = 0.5, between = 0)

sim <- sim_spatial(n = 400, groups = groups, R_true = R_true, seed = 1)
sim                                  # summary of the groups
dim(sim$X)                           # 400 sites x 24 variables
cor(sim$components$residual)         # close to R_true
cor(sim$X)                           # inflated by the shared large-scale fields
```

Regimes available in each simulator (`range` is the practical range, in units of the side
of the unit square):

| Simulator | Regimes and parameters |
| --- | --- |
| `sim_spatial(n, groups, ...)` | `"exponential"` (`range`, `anisotropy`, `angle`), `"matern"` (the same plus `smoothness` = 0.5, 1.5 or 2.5; smoother fields with the same practical range), `"none"` |
| `sim_temporal(n_times, groups, ...)` | `"ar1"` (`phi`), `"seasonal"` (`period`, `phase`), `"trend_shift"` (`slope`, `jump`, `shift_at`), `"none"` |
| `sim_spatiotemporal(n, n_times, groups, ...)` | `"separable_ar1"` (`range`, `phi`), `"separable_seasonal"` (`range`, `period`), `"moving"` (`range`, `velocity`), `"temporal_common"` (`phi`), `"none"`; the spatial fields accept `smoothness` |

Every simulator returns the data `X`, the coordinates and/or times, the true group of each
variable, the separate components (`large`, `residual`, `class`), `R_true`, the true
class of each observation, and the seed. Optional observation classes
(`classes = list(k = 3, type = "random" | "contiguous", effect = 1)`) provide the truth
for k-means and SOM.

## 10. Quick start: diagnosing each variable

```r
sig <- diagnose_variables(sim)            # one dependence signature per variable
sig[, c("variable", "moran_p", "range", "nugget_ratio", "aniso_ratio")]
classify_regime(sig)                      # "none", "isotropic", "anisotropic", ...
```

`diagnose_variables()` looks at **one variable at a time** — its variogram, its
autocorrelation, its periodogram — and never at the correlations between variables, so
that the later choice of a correction does not reuse the associations it is meant to
estimate. It works on spatial (`coords`), temporal (`time`) and spatio-temporal
(`coords` + `index`) data, or directly on a simulation object:

| Family | Signature components |
| --- | --- |
| Spatial | Moran's I and p-value; nugget, sill and practical range of the best of an exponential and a spherical variogram; anisotropy ratio and direction |
| Temporal | lag-1 autocorrelation, AR(1) coefficient and AIC order, Ljung-Box p-value, dominant seasonal period and its strength, trend R-squared |
| Spatio-temporal | the spatial statistics on time slices, the temporal ones on sites, and a shift indicator that detects a pattern moving in one direction |

With `B > 0` it adds a parametric bootstrap of the range and of phi (standard error and
90% interval). `classify_regime()` turns a signature into a regime type with explicit
thresholds (`default_thresholds()`), calibrated on the stage 3 simulations: a variable is
dependent when Moran's or the Ljung-Box p-value is below 0.05, anisotropic when the
anisotropy ratio is above 2.5, seasonal when the seasonal strength is above 0.2.

## 11. Quick start: groups and shared bases

```r
sig <- diagnose_variables(sim)
gr  <- group_regimes(sig, data = sim, B = 50, seed = 1)   # groups + stability
gr[, c("variable", "group", "group_type", "stability", "ambiguous")]

members <- gr$variable[gr$group == "g1"]
b <- build_group_basis(sim$X[, members], coords = sim$coords,
                       type = gr$group_type[gr$group == "g1"][1], seed = 1)
b                                         # one basis for every variable of the group
```

`group_regimes()` keeps the variables without dependence in their own group (no basis)
and groups the others on their signature components, using only the signatures. With
`B > 0` it repeats the grouping on subsamples and flags unstable (ambiguous) variables.
`build_group_basis()` selects one basis for the whole group by forward selection among
orthonormal candidates — Moran eigenvector maps in space, harmonics, trend and slow
temporal eigenvectors in time, and their products in space-time — with permutation tests
and a cap on the share of degrees of freedom (`kn_max`). `union_basis()` joins the bases
of two groups for the pairs of variables across them.

## 12. Quick start: correcting and verifying

```r
sig <- diagnose_variables(sim)
gr <- group_regimes(sig)
diagnosis_table(sig, gr)                   # what the diagnosis found, variable by variable
cr <- correct(sim, gr)                     # the defaults below
cr$R                                       # corrected correlation matrix (for PCA)
head(cr$X)                                 # corrected data (for k-means, SOM, MDS, HDBSCAN)
v <- verify(cr, sim)                       # what was removed and what was kept
report_correction(cr, sim, v, file = "report.html",    # one HTML page for the user
                  signature = sig, grouping = gr)
correct(sim, gr, scope = "groups")         # one basis per group, unions across groups
run_app()                                  # the same workflow, point and click
```

Defaults, chosen from the simulation studies:

- **The correlation matrix `R`** (for PCA): in space and in time, one global basis of the
  dependent variables, conditional mode; in space-time, whitening group by group.
- **The data matrix `X`** (for MDS, k-means, SOM, HDBSCAN): each variable corrected with
  a basis selected on that variable alone (`data_basis = "own"`). A shared basis takes
  part of every variable's class signal away; a basis per variable keeps the classes.
- Variables without dependence are never corrected among themselves.

### Choosing the structures yourself: the menu

When you know the system (the period of a cycle, the direction of a current), say so
with a menu; the parameters are still searched within it, with the same tests and cap
on k/n:

```r
m <- basis_menu(term_trend(), term_cycle(period = 12))      # time: trend and an annual cycle
m <- basis_menu(term_trend(), term_mem(W = "auto"))         # space: choose the neighbourhood graph
m <- basis_menu(term_trend(), term_aem(direction = 45))     # space: a directional process
cr <- correct(sim, gr, menu = m)                            # one menu, or a list named by group
```

Every term is a function of the positions only, so it can be shared by a group. Terms
that are functions of the variable itself (its lags, the mean of its neighbours) are
filters: they change what is estimated (the correlation of the innovations) and are not
in the menu.

### What `verify()` reports

The residual structure of each variable (tested against white noise with the same basis
removed), Moran correlograms and a multivariate Mantel correlogram, how each
correlation changed (Delta R), the part of each correlation carried by the large scale,
the local Moran's I of the corrected variables (where structure is left, in space), and
the warning of signal at the scale of the basis (question Q3): when much of the
association between variables lived at the scale of the basis, classes that are regions
or periods would be removed by the correction. The same Moran correlograms, local
Moran's I and (in time) Ljung-Box tests are given for the raw data
(`correlograms_before`, `lisa_before`), so that what the diagnosis saw can be compared
with what the correction left.

### Exploring the diagnosis before the correction

Everything the diagnosis computed is kept for the user to explore:
`diagnosis_table(sig, gr)` puts each variable's signature (Moran's I and its p-value,
the variogram and its range, nugget and anisotropy in space; the lag-1 autocorrelation,
the AR coefficient, the Ljung-Box test, the period and its strength and the trend in
time) beside its regime, group and stability. The report of `report_correction()` opens
with this table, and the app shows it in step 3 and lets you download it, together with
the local Moran's I before and after.

### The warning of signal at the scale of the basis (Q3)

The correction removes what varies at the scale of its basis. When the classes you look
for are regions or periods, they live at that scale too and are removed with the
dependence. `verify()` measures how much of the association between variables was
carried by the removed part:

```
Q3 statistic = mean over all pairs of variables of |cov_large| / (sd_x sd_y)
```

where `cov_large` is the covariance of the two variables' parts on the basis. The warning
(`v$q3$warning`) is raised when the statistic exceeds the threshold of the data family:

| Data | Threshold | Validated on new seeds (warnings without / with contiguous classes) | At the size of the case study |
| --- | --- | --- | --- |
| Spatial | 0.16 | at most 1.2% / at least 98.8% (stage 7 grid) | n = 628 sites: 0% / 97-99% |
| Temporal | 0.19 | at most 0.3% / 99.7% (stage 7 grid) | T = 93 months: 0-2% / 99-100% |
| Space-time | 0.13 | 0% / 99% (stage 7 grid) | 628 sites x 93 months: 0% / 100% |

The thresholds differ by family because the statistic grows with the size of the basis
even without classes: the large-scale parts of independent autocorrelated variables are
correlated by chance, more so in short series. Each threshold is the midpoint between the
97.5th percentile of the truths without contiguous classes and the 5th percentile of the
contiguous truth (`analysis/07d_q3_recalibration.R`), checked at the size of the real
case by `analysis/08c_q3_size_check.R` (criterion: at most 5% of warnings without
contiguous classes, at least 90% with them). A statistic normalised by a reference of
independent surrogates (`analysis/08a_q3_reference_pilot.R`) was tried and not adopted:
it separated the truths worse, because contiguous classes also make each variable more
autocorrelated, so the reference rises with them and hides the signal. The threshold can be set with `q3_threshold`.

When the warning fires, the correlations are not wrong: much of the association lived at
the scale of the basis. If the classes you seek are regions or periods, analyse the raw
data, or use clustering constrained by contiguity (section 5).

`correct()` has three modes: `"conditional"` (residuals on the group bases),
`"whiten"` (the observations are whitened with each group's estimated covariance) and
`"scale"` (one correlation matrix per band of scales). What the simulations show so far:
the error of using different bases for two variables is real but small (of order
rho k/n); the larger error is the part of a shared large-scale process that no basis
captures, which grows with the roughness of that process. A correction cannot always
remove it: at scales finer than the basis, a shared process and a genuine association
leave the same trace. In space-time, observation over time separates them better, and
whitening reduces the error substantially.

### The basis is feature engineering: three rules for your own features

The basis of the correction is a set of positional features (trend surfaces, Moran
eigenvector maps, harmonics, autoregressive lags). In predictive regression, features
that explain more are better. Here they are not: features that explain more of the
large scale also remove more of the association you want to keep, which is why the
bases are selected with tests and a cap on k/n. Three consequences for the features
you build yourself:

1. **Build derived features before the diagnosis.** Anomalies from a climatology,
   indices (such as NDVI), ratios and transformations change a variable's dependence
   regime: an anomaly is already a seasonal correction, and a ratio mixes the regimes
   of its two variables. Diagnose and correct the variables you will actually analyse.
2. **Do not build features from correlations between variables before correcting**
   (for example, PCA axes): the diagnosis uses each variable alone, so that the
   associations it is meant to estimate are not used twice.
3. **Reuse the basis in a downstream model.** If a regression on the corrected data
   uses positional features (smooths of the coordinates, Moran maps, lags), use the
   basis that `correct()` returns (`cr$bases`). Different positional features bring
   back the error of correcting two variables with different bases.

With autocorrelated data, positional features also inflate the apparent explanatory
power of predictive models under random validation; blocked validation in space or
time removes much of it (Roberts et al. 2017, Ecography 40:913-929; Ploton et al. 2020,
Nature Communications 11:4540).

## 13. License

MIT — see `LICENSE.md`.
