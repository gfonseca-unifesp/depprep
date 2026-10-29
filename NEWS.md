# depprep 0.1.0

## Package, menu and interface (stage 8)

* `basis_menu()` with `term_trend()`, `term_cycle()`, `term_mem()` and `term_aem()`:
  the user chooses the structures a basis is built from, and `build_group_basis()` and
  `correct()` search the parameters within the menu. Without a menu the
  bases are the same as before.
* `verify()` gives the local Moran's I of the corrected variables (`lisa`), to see
  where structure is left.
* `report_correction()` writes one self-contained HTML page that explains a correction.
* `run_app()` opens a point-and-click interface to the whole workflow.
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
