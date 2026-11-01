# depprep 0.1.0

## Drivers (stage 9)

* `variation_partition()` splits the variance of each variable into the parts explained
  by its position only, by the drivers only, by both, and the residual, with adjusted
  R-squared (Peres-Neto et al. 2006).

## Package, menu and interface (stage 8)

* `basis_menu()` with `term_trend()`, `term_cycle()`, `term_mem()` and `term_aem()`:
  the user chooses the structures a basis is built from, and `build_group_basis()` and
  `correct()` search the parameters within the menu. Without a menu the
  bases are the same as before.
* `verify()` gives the local Moran's I of the corrected variables (`lisa`), to see
  where structure is left, and the same descriptions of the raw data before the
  correction (`correlograms_before`, `lisa_before`, and the Ljung-Box test in time).
* `diagnosis_table()` gathers what the diagnosis found for each variable (signature,
  regime, group, stability) in one table to read or save.
* `report_correction()` writes one self-contained HTML page that explains a correction,
  opening with the diagnosis and showing every test before and after.
* `run_app()` opens a point-and-click interface to the whole workflow, with longitude
  and latitude as coordinates. The diagnosis shows each variable's variogram, Moran
  correlogram and map, or its series with the linear and cubic trends, ACF and
  periodogram; the menu of structures explains each term and draws the neighbourhood
  graph and the direction of the directional maps. The diagnosis and the local Moran's
  I can be downloaded.
* `term_aem(direction)` takes the direction counter-clockwise from the x axis, as
  documented; it was taken as a compass bearing before.
* Example tables of the Foz do Amazonas ocean case in `inst/extdata/`.
* The diagnosis, the bases and the correction stop with a clear message on missing
  values instead of filling them.
* The Q3 warning of `verify()` uses one threshold per data family (amendment
  protocol-v1.2 of the simulation study).

## Method and simulation study (stages 0 to 7)

* Simulation generators with known truth (`sim_spatial()`, `sim_temporal()`,
  `sim_spatiotemporal()`).
* The diagnosis of each variable (`diagnose_variables()`, `classify_regime()`), the
  grouping by regime (`group_regimes()`) and the shared bases (`build_group_basis()`,
  `union_basis()`).
* The correction (`correct()`, with the conditional, whitening and scale modes) and its
  verification (`verify()`), with the defaults chosen in the simulation studies.
