# depprep

**Dependence-aware preprocessing for multivariate exploratory analysis.**

> **Status:** early development — stage 2 of 7 (simulation generators) is complete.
> The correction itself is not implemented yet; do not use the package on real data.

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
| `R/` | Package functions (currently the simulation generators) |
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
| 3 | Diagnose each variable's dependence | — | next |
| 4 | Group regimes and build shared bases | — | planned |
| 5 | Correct, and show that shared bases avoid the distortion | — | planned |
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
| `sim_spatial(n, groups, ...)` | `"exponential"` (`range`, `anisotropy`, `angle`), `"none"` |
| `sim_temporal(n_times, groups, ...)` | `"ar1"` (`phi`), `"seasonal"` (`period`, `phase`), `"trend_shift"` (`slope`, `jump`, `shift_at`), `"none"` |
| `sim_spatiotemporal(n, n_times, groups, ...)` | `"separable_ar1"` (`range`, `phi`), `"separable_seasonal"` (`range`, `period`), `"moving"` (`range`, `velocity`), `"temporal_common"` (`phi`), `"none"` |

Every simulator returns the data `X`, the coordinates and/or times, the true group of each
variable, the separate components (`large`, `residual`, `class`), `R_true`, the true
class of each observation, and the seed. Optional observation classes
(`classes = list(k = 3, type = "random" | "contiguous", effect = 1)`) provide the truth
for k-means and SOM.

## 7. Reproducibility

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

## 8. License

MIT — see `LICENSE.md`.
