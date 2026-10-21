# depprep

**Dependence-aware preprocessing for multivariate exploratory analysis.**

> **Status:** early development — stages 0 to 4 are complete (simulation generators, the
> diagnosis of each variable, grouping and shared bases); stage 5 (the correction) is in
> progress. Do not use the package on real data yet.

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
| `R/` | Package functions: simulation generators and the diagnosis of each variable |
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
| 5 | Correct, and show that shared bases avoid the distortion | — | in progress |
| 6 | Pilot with PCA, MDS, k-means and SOM | — | planned |
| 7 | Full simulations | — | planned |

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

## 9. Known limits: why purely spatial regimes are hard to tell apart

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

## 10. Reproducibility

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

## 11. License

MIT — see `LICENSE.md`.
