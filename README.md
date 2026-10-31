# depprep

**Dependence-aware preprocessing for multivariate exploratory analysis.**

> **Status:** stages 0 to 7 are complete: simulation generators, the diagnosis of each
> variable, grouping and shared bases, the correction, a pilot and the full
> pre-registered simulation study (protocol-v1, with two labelled amendments). Stage 8
> (package, menu of structures, report and app) is in progress. The real-data application
> (stage 9) has not been done yet, and missing values are not supported.

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
| `analysis/` | One script per project stage, numbered in the order they are run |
| `results/` | Created by the scripts; not tracked by git |
| `data/raw/` | Downloaded raw data; not tracked by git, never modified |
| `renv.lock` | Exact versions of every R package used |

## 4. Project stages

| Stage | Question | Code | Status |
| --- | --- | --- | --- |
| 0 | Reproducible environment | package skeleton, `renv.lock` | done |
| 1 | Is there a gap? Which real data? | `analysis/01a_*`, `01b_*`, `01c_*` | done |
| 2 | Can we simulate data with known truth? | `R/sim_*.R`, `analysis/02_*`, `02b_*` | done |
| 3 | Diagnose each variable's dependence | `R/diagnose*.R`, `R/classify.R`, `analysis/03_*` | done |
| 4 | Group regimes and build shared bases | `R/group.R`, `R/basis.R`, `analysis/04_*` | done |
| 5 | Correct, and measure the mechanism (H2a, H2b) | `R/correct.R`, `R/verify.R`, `analysis/05*` | done |
| 6 | Pilot with PCA, MDS, k-means and SOM | `R/downstream.R`, `R/comparators.R`, `analysis/06*` | done |
| 7 | Full pre-registered simulations | `analysis/07*` | done |
| 8 | Package: menu of structures, report, app | `R/menu.R`, `R/report.R`, `R/app.R`, `analysis/08*` | in progress |
| 9 | Real-data application | — | planned |

### Stage 1 — real-data case study (Foz do Amazonas)

The application uses **Level-3 satellite time series** (gridded, not gap-filled) of the
Foz do Amazonas basin, Brazilian Equatorial Margin, in two case studies sharing one
rectangle: the **ocean** (the marine part of the basin, official ANP limit) and the
**land** (within 200 km of that coast). L3 is required because gap-filled (L4) products
add spatial autocorrelation that belongs to the interpolation, not to the process.

| Script | What it does | Needs |
| --- | --- | --- |
| `analysis/01a_fza_study_area.R` | Builds the study area from the ANP basin limits and ETOPO1 topography | internet |
| `analysis/01b_fza_data_feasibility.R` | Catalogue of candidate L3 products; monthly cloud-free coverage of the products served anonymously by NOAA ERDDAP | internet |
| `analysis/01c_fza_nasa_feasibility.R` | Same for NASA products (SMAP salinity and soil moisture, MODIS land and ocean colour, GPM IMERG), subset with OPeNDAP | a free [NASA Earthdata](https://urs.earthdata.nasa.gov) login stored in `~/_netrc` (`machine urs.earthdata.nasa.gov login <user> password <password>`); the script never prints it |

### Stage 2 — simulation generators

| Script | What it produces |
| --- | --- |
| `analysis/02_simulation_examples.R` | Schemes of the generator, galleries of the candidate regimes, variograms, space-time maps and Hovmöller diagrams, the apparent-association matrices, observation classes and the time to generate one replicate |
| `analysis/02b_generator_recovery.R` | Recovery of the simulated range and AR coefficient in 50 replicates per scenario, and the error of the residual and raw correlations |

### Stage 3 — diagnosis of each variable

| Script | What it produces |
| --- | --- |
| `analysis/03_diagnostics.R` | Signatures of every variable in the simulation grid (200 replicates per scenario, run in parallel and resumable), classification accuracy, scale errors and bootstrap coverage |
| `analysis/03_diagnostics_figures.R` | Scheme of the diagnosis, a didactic variogram, the typical signature of each regime, confusion matrices, scale uncertainty and the threshold calibration curves |

### Stage 4 — grouping and shared bases

| Script | What it produces |
| --- | --- |
| `analysis/d4_scenarios.R` | The simulation grid (groups, scenarios, truth, one replicate), shared by the stage scripts |
| `analysis/04_grouping.R` | Grouping of the stage 3 signatures (hclust and mclust, 200 replicates per scenario), stability over subsamples, and one basis per group with the estimated and the true grouping: k/n, capture of the true large-scale structure and absorption of the true residual |
| `analysis/04_grouping_figures.R` | Scheme, didactic Moran eigenvector maps, one replicate step by step, agreement with the true groups, stability, degrees of freedom and the k/n threshold curves |

### Stage 5 — correction and the mechanism

| Script | What it produces |
| --- | --- |
| `analysis/05a_h2_mechanism.R` | H2a: the error of using different bases against the distance between them (and its theory); H2b: the shared structure a basis leaves behind, by roughness of the shared field and group size; the calibration of the second grouping level |
| `analysis/05b_correction.R` | The correction on the grid (D4 plus smooth fields and "same regime, different fields"): eight approaches, from no correction to the true grouping, with the error of the residual correlations, capture, degrees of freedom and time |
| `analysis/05c_whitening_st.R` | Why whitening helps in space-time: its spatial and temporal parts, by kind of shared process |
| `analysis/05_correction_figures.R` | The didactic scheme, example correlation matrices, and the figures of the three scripts |

### Stage 6 — pilot with the four exploratory methods

| Script | What it produces |
| --- | --- |
| `analysis/06_pilot.R` | Twelve pilot scenarios (spatial rough and smooth, temporal, space-time; with no association, with true blocks and random classes, and with contiguous classes), 20 replicates each: every correction approach followed by PCA, MDS (Euclidean and Gower), k-means and SOM, scored with `run_downstream()`; SBSS and MAF (`decorrelate_latent()`) as comparators in space; paired contrasts and the replicates needed for stage 7 |
| `analysis/06_pilot_figures.R` | A didactic k-means example, spurious patterns, recovery of the truth, contiguous classes, the temporal default and the replicates needed |
| `analysis/06b_pilot_absorption.R` | How much of the class signal each basis removes, against what its size alone implies |

### Stage 7 — full simulations (pre-registered)

| Script | What it produces |
| --- | --- |
| `analysis/07_full_grid.R` | The grid of protocol-v1: 21 base scenarios x five truths (plus perturbations), 100 replicates in space and time and 50 in space-time; every correction approach and every output method |
| `analysis/07b_h2a.R` | The between-basis error (H2a) on new seeds |
| `analysis/07_analysis.R` | The pre-registered paired tests with Holm's correction, and the decision rule of the project |
| `analysis/07c_amendment.R` | Amendment protocol-v1.1: the correlation error against the correlation of the true signal |
| `analysis/07d_q3_recalibration.R` | Amendment protocol-v1.2: one threshold of the Q3 warning per data family, calibrated on stage 7 and validated on new seeds |
| `analysis/07e_q3_cause.R` | Why the Q3 statistic depends on the family (exploratory) |

The pre-registered results stay the main ones; the amended results are reported beside
them and labelled.

### Stage 8 — package

| Script | What it produces |
| --- | --- |
| `analysis/08a_q3_reference_pilot.R` | Pilot of a Q3 statistic normalised by an independence reference that keeps each variable's autocorrelation |
| `analysis/08b_menu_check.R` | Descriptive check of the menu of structures: known period, chosen neighbourhood graph, directional maps |

## 5. Getting started

```r
# R >= 4.3 (developed with R 4.6.1)
install.packages("renv")
renv::restore()          # installs the exact package versions in renv.lock
devtools::load_all()     # loads the package from source
devtools::test()         # runs the tests (a few take about a minute)
```

Run the scripts from the repository root, in order (for example
`Rscript analysis/02_simulation_examples.R`). Outputs go to `results/`.

## 6. Quick start: simulating data with known truth

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

## 7. Quick start: diagnosing each variable

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

## 8. Quick start: groups and shared bases

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

## 9. Quick start: correcting and verifying

```r
gr <- group_regimes(diagnose_variables(sim))
cr <- correct(sim, gr)                     # the defaults below
cr$R                                       # corrected correlation matrix (for PCA)
head(cr$X)                                 # corrected data (for k-means, SOM, MDS, HDBSCAN)
v <- verify(cr, sim)                       # what was removed and what was kept
report_correction(cr, sim, v, file = "report.html")   # one HTML page for the user
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
or periods would be removed by the correction.

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

## 10. The app: the same workflow without code

`run_app()` opens a point-and-click interface (Shiny) in the browser. It runs the same
functions as the code above and, at the end, gives you the R code that reproduces every
choice you made.

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

**Spatial** (one map): one row per site, with its two coordinates.

| x_km | y_km | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- | --- |
| 412.5 | 9876.0 | 1.84 | 28.31 | 0.21 | 210.4 |
| 440.3 | 9876.0 | 1.12 | 28.45 | 0.15 | 168.9 |
| 412.5 | 9903.8 | 2.07 | 28.22 | 0.24 | 231.0 |

**Temporal** (one series): one row per time, with a numeric time.

| time | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- |
| 1 | 1.52 | 28.10 | 0.19 | 195.2 |
| 2 | 1.61 | 28.34 | 0.20 | 201.7 |
| 3 | 1.95 | 28.62 | 0.23 | 224.8 |

**Space-time** (a series at every site): one row per site and time, with a site
identifier, the site's coordinates (repeated on each of its rows) and the time.

| site | x_km | y_km | time | chl | sst | kd490 | poc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 412.5 | 9876.0 | 1 | 1.84 | 28.31 | 0.21 | 210.4 |
| 1 | 412.5 | 9876.0 | 2 | 1.90 | 28.40 | 0.22 | 214.1 |
| 2 | 440.3 | 9876.0 | 1 | 1.12 | 28.45 | 0.15 | 168.9 |
| 2 | 440.3 | 9876.0 | 2 | 1.20 | 28.51 | 0.16 | 172.3 |

Rules:

- **At least two numeric variables**, and **no missing values** in the variables, the
  coordinates or the time (fill or drop them before; gaps are not supported yet).
- **Time is a number**: a month counter (1, 2, 3, ...) or a decimal year (2015.25), not a
  date text such as "2015-04". Equal steps are expected.
- **Coordinates in a projected system** (for example UTM, in km or m). Longitude and
  latitude in degrees are accepted with a warning, since distances in degrees are
  distorted away from the equator.
- **Space-time needs a complete grid**: every site has a row at every time, once.
- Text columns (names, dates as text, labels) are ignored as variables; you can keep them
  in the file.

The same layouts work in code: `diagnose_variables()` takes the matrix of variables with
`coords` (spatial), `time` (temporal), or `coords` of the sites plus `index` (a data
frame with the `site` and `time` of each row, space-time).

### The six steps

| Step | What you do | What you get |
| --- | --- | --- |
| 1. Data | Load a CSV with one row per observation and mark the position columns (x and y, time, site); choose spatial, temporal or space-time | A check of the columns and the first rows |
| 2. Context | Say which technique the data are for (PCA, or MDS, k-means, SOM and HDBSCAN) and whether the classes you seek could be regions or periods | Which corrected matrix the technique reads, and a warning when correcting would remove the classes you are looking for |
| 3. Diagnosis | Look at each variable's dependence and its plots (variogram, ACF, periodogram); change a regime by hand, or build the menu of structures of a group (trend, cycle with its period, Moran maps with a chosen graph, directional maps) | The groups of variables that will share a basis |
| 4. Correction | Accept the default (the data matrix with one basis per variable, the correlation matrix with the shared basis) or choose another mode (conditional, whitening) and other bases for the correlation matrix | The bases used, with their size (k, k/n) and first terms |
| 5. Verification | Check what was removed and what was kept | The autocorrelation of each variable before and after, the change in the correlations, the warning of signal at the scale of the basis, and a preview of the technique (PCA or k-means) before and after |
| 6. Export | Download the results | The corrected data (CSV), the corrected correlation matrix (CSV), the HTML report of `report_correction()` and the R code |

Limits shown in the app: missing values are not supported yet, coordinates in degrees
should be projected first, and data with more than a few thousand observations are slow
(the correction of a 628-site x 93-month grid takes hours).

## 11. Known limits

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

- **Missing values.** The diagnosis, the bases and the correction stop on missing values.
- **Filters other than the autoregressive ones** (SAR/CAR whitening in space, STAR for
  moving patterns) are not implemented.
- **Coordinates in degrees** should be projected first; distances in degrees are
  distorted away from the equator.

## 12. Reproducibility

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

## 13. License

MIT — see `LICENSE.md`.
