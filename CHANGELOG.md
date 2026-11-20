# Changelog

## 1.0.0 (first public release)

Numerics
- Caputo solver: Anderson-accelerated fixed-point iteration (depth 10) with a relative stopping test
  replaces plain Picard iteration; it converges for small ν (0.05) and amplitudes of order one, and a
  non-converging step is reported as an error.
- Caputo solver: L1 weights evaluated as w^{1−ν} expm1((1−ν) log1p(τ_k/w)), which removes the catastrophic
  cancellation on strongly graded meshes (r up to 10).
- Caputo model with ν = 1 runs the conservative split-step solver with dt = T/N, with a warning;
  ν ≥ 0.95 carries a dissipation warning.
- Space-fractional model with α ≤ 1 carries a warning; exact mode is refused there.
- Convergence studies save the same number of snapshots at every level.
- Exact mode, classical soliton: the reference is the soliton periodized over the box, Σₘ p(x + mL, τ), so a
  soliton may cross x = ±L/2. Exact mode is refused when soliton and images overlap (2 sech(ηL) > 10⁻³, field
  `L`) or when the carrier exceeds the grid (2|ξ| ≥ πNx/L, field `Nx`).
- Exact mode: a reference floor with four parts (`exact.refFloorParts`: periodic images, spatial truncation of
  the soliton spectrum, profile residual, round-off) is reported with each run; observed orders of errors below
  ten times the floor are omitted, and a warning names the dominant part and the L or Nx that removes it. The
  image and spatial parts are estimates, not bounds.
- Convergence studies test their own errors a posteriori: from three consecutive levels they fit
  e = C h^p + S (p solves (e_a − e_b)/(e_b − e_c) = (h_a^p − h_b^p)/(h_b^p − h_c^p), so the step ratios need not be
  equal; `Simulation::fitPowerLaw`, `fit_power_law`) and withhold an order when p ≤ 0, when the step-independent part S exceeds 10 % of the error
  (`convergence.saturated`, with a warning), or when the error exceeds 10 % of the peak modulus
  (`convergence.preasymptotic`); the step-independent parts are reported (`convergence.dtIndependent`, the largest S found). A study with only two levels cannot separate S and
  publishes no order (`convergence.untested`, with a warning that names the dt or N giving three levels).
  `src/php/tests/sweep.php` runs this test over the 636 admitted configurations of two grids of classic-exact
  and writes every configuration, level, error and published or withheld order (`results/sweep.txt`, `.json`).
- Exact mode, split-step models: when the output times under-sample β, γ or g (a coefficient departs from its
  linear interpolation between outputs by more than 5 % of max(1, max|·|)), the exact error is measured at
  every step (`exact.errorSampling = "every step"`), with a warning.
- Exact mode, fractional soliton: the stationary solution computed on the box is described (`meta.profile`:
  localized, periodic or homogeneous branch; Petviashvili iterations or non-convergence), with warnings; the
  residual of the profile actually used is reported (`meta.profileResidual`); an off-node centre whose
  Fourier-shifted profile is not stationary on the grid (relative residual > 10⁻⁹) is moved to the nearest node,
  and a residual that remains above 10⁻⁹ enters the reference floor.
- Reference integrals G = ∫g and τ = ∫β use an adaptive Simpson table (panels doubled from 4096 until two
  tables agree to 10⁻⁸), so rapidly varying coefficients do not corrupt the exact solutions or the gain diagnostic.
- Split-step runs compare the quadrature that the solver applies (Simpson on half steps for β and g; one
  gain-weighted Simpson panel per step for the nonlinear phase γ; one Simpson panel per step for a
  time-dependent V(x,t), checked at 16 grid points and at the maxima of |V| at t = 0 and T) with the reference integrals
  (`meta.quadrature`) and warn when they differ by more than 10⁻⁶; the adaptive reference tables are also
  cross-checked against a composite 3-point Gauss–Legendre rule, which detects coefficients aliased on the
  Simpson nodes (e.g. g = cos 32768πt on [0, 1]). These checks are estimates, not a certification; convergence levels that fail this test are
  listed in `convergence.unresolved` and give no observed order.
- The step actually taken, T/N with N = round(T/dt), is reported (`meta.dtEffective`) and is the abscissa of
  convergence studies. For ν = 1 the Caputo model always uses dt = T/N, in both implementations; a dt that
  differs is ignored with a warning.
- The mathematical domain of every model (0 < ν ≤ 1, 0 < α ≤ 2, 0 < dt ≤ T, N ≥ 1, r ≥ 1, Nx a power of two,
  η > 0) is enforced in the numerical core of both implementations, not only in the web API.
- `SplitStepSolver` and `CaputoSolver` validate their own contract (L > 0, Nx a power of two, 0 < α ≤ 2,
  0 < ν ≤ 1, T > 0, dt > 0, N ≥ 1, r ≥ 1) when used directly as library classes. `CaputoSolver` accepts
  ν = 1, where the L1 weights reduce to backward Euler on the graded mesh (first order, dissipative); the web
  interface and `Simulation` run ν = 1 with the split-step solver.

Interface and API
- Home page (`index.php`): the model in brief, the three variants, what the laboratory does, how it is validated,
  how to cite and a static preview of the preset classic-collision (`assets/img/collision-preview.svg`, built by
  `tests/make_home_preview.php`). The laboratory moves to `lab.php`, with unchanged behaviour. The Zenodo DOI is set
  with `MANAKOV_DOI`; links that still carry the release placeholder are shown as pending on the home page.
- Input errors return the offending field; warnings for rounded Nx and capped grading exponent.
- Formula preview sent by POST (`api/latex.php` also accepts GET); expressions limited to 2000 characters.
- Plotly bundled locally; `.htaccess` rules for `src/`, `tests/` and `config.php`.
- Version, licence, source-code and citation links in the footer; keyboard navigation of the model selector.
- Run report shows the iteration statistics of the nonlinear solver; `meta.picardSteps` lists the iterations
  of every step.
- The scheme badge of the equation panel follows ν and r (ν = 1: Strang split-step with uniform Δt = T/N;
  automatic grading capped at r = 10, with the expected order min(rν, 2 − ν)) and, after a run, the server's
  `meta.r` and `meta.expectedOrder`.
- Convergence table names the dominant part of the reference floor, the unresolved, saturated and
  pre-asymptotic levels, and whether the error was measured at every step; the run summary and the data panel
  show the mesh actually used (uniform Δt = T/N, uniform Strang steps for ν = 1, or the graded Caputo mesh); the cost estimate distinguishes work (∝ N²Nx) from memory (∝ N·Nx); the
  classic-trap description states β ∈ [0.2, 1.8]; the theory panel describes the step quadrature of the gain
  and the periodized reference.

Repository
- Source code under `src/` (`src/php/`, `src/python/`); licence in `LICENSE.txt`.
- Reproducibility material only: `results/` holds the archived outputs and `reproduce/make_figures.py` redraws the
  figures of the article into `reproduce/figures/`. The sources of the article are not part of the repository.

Testing and reproducibility
- `src/php/tests/run_tests.php` (30 tests) and `src/php/tests/reproduce.php` (every number of Section 3 except
  run times and the N = 2400, 4800 levels); section names follow the table numbers of the article (Tables 3-5),
  with sections `mittag` (analytic Mittag-Leffler test of the Caputo solver) and `spatial` (Nx refinement at
  fixed L). The Mittag-Leffler row takes the maximum over all output times (attained at t = T). Table 5 was measured
  again with this version; `meta.runtime_total` now starts when the run is set up, so it includes the
  reference tables, every check, the diagnostics and the convergence runs. The README maps every result of the article to its command and output.
- New Python port of the numerical core (`src/python/`) with pytest suite, PHP-Python cross-check and `caputo_refinement.py`.
- The PHP-Python cross-check compares P1 and P2 separately, the final complex field (phases included), the
  Hamiltonian, the exact-error history and the iteration counts of every step, and a case passes only if all agree.
