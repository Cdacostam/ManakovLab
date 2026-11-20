# Manakov Lab — Python implementation

A NumPy port of the numerical core of Manakov Lab, written line by line from the PHP classes in
`../php/src/Numerics` (`src/php/` in the repository). It solves the same three models with the same algorithms,
settings and conventions, needs no PHP at run time, and is used to cross-check the PHP results (to about 1e-12),
to produce the figures of the accompanying article and to run refinements beyond the web work budget
(`caputo_refinement.py`). Because it is a port, the cross-check detects translation, language and FFT errors,
not errors in the specification of the schemes; those are tested by the exact and manufactured solutions.

```
i q_t + β(t) q_xx + 2γ(t)|q|² q + V(x,t) q = i g(t) q                 classical
i q_t − β(t)(−∂_x²)^{α/2} q + 2γ(t)|q|² q + V q = i g(t) q              Riesz space-fractional; tested at α = 1.2, 1.5, 2, recommended 1 < α ≤ 2
i ᶜD_t^ν q + β(t) q_xx + 2γ(t)|q|² q + V q = i g(t) q                   Caputo time-fractional, 0 < ν < 1 (α = 2)
```

with q = (q₁, q₂) and |q|² = |q₁|² + |q₂|².

## Installation

Python ≥ 3.10 and NumPy are required.

```bash
cd src/python
pip install -e ".[test,figures]"   # figures: matplotlib, for reproduce/make_figures.py
```

## Usage

```python
from manakov_lab import run_preset, problem, Simulation

r = run_preset("classic-exact")                       # one of the eight presets of the web interface
print(r["convergence"]["order"])                      # observed orders of the split-step scheme

r = run_preset("caputo-decay", params=dict(nu=0.7, N=1200))
d = r["diagnostics"]
print((d["P1"][-1] + d["P2"][-1]) / (d["P1"][0] + d["P2"][0]))

# custom coefficients are Python callables
import math
P = problem("classic-collision", coeffs=dict(beta=lambda t: 1 + 0.3 * math.sin(t)))
r = Simulation(P).run()
```

The result dictionary has the keys `x`, `t`, `initial`, `final`, `diagnostics` (`P1`, `P2`, `peak`, `H`,
`gain`), `warnings`, `meta`, and, when requested, `exact` and `convergence`, with the same meaning as the
JSON returned by `src/php/api/run.php`. Invalid input raises `SimulationError`, whose `field` attribute names
the offending parameter. The core enforces the mathematical domain of each model, as the PHP core does:
0 < ν ≤ 1, 0 < α ≤ 2, 0 < dt ≤ T, N ≥ 1, r ≥ 1, Nx a power of two, η > 0; the work-budget limits of the web
API do not apply.

In exact mode the reference must solve the same discrete periodic problem as the integrator. The classical
reference is the soliton periodized over the box, Σₘ p(x + mL, τ), so the soliton may cross x = ±L/2; exact
mode is refused only when soliton and images overlap (2 sech(ηL) > 10⁻³). What the reference does not resolve
is bounded a priori as a reference floor, `exact['refFloorParts']` (periodic images, spatial truncation of the
soliton spectrum, residual of the fractional profile, round-off), summed in `convergence['refFloor']`; observed
orders are omitted where the error is below ten times it, and a warning names the dominant part. For the
fractional soliton, `meta['profile']` describes the stationary solution computed on the box (localized,
periodic or homogeneous branch; Petviashvili iterations, `None` if not converged), an off-node centre whose
Fourier-shifted profile is not stationary on the grid (relative residual above 10⁻⁹) is moved to the nearest
node, and `meta['profileResidual']` gives the residual of the profile used. The integrals G and τ use an
adaptive Simpson table (`core.antiderivative_adaptive`). Split-step runs compare their step quadrature of β, γ
and g with these integrals (`meta['quadrature']`), warn above 10⁻⁶ and list unresolved convergence levels in
`convergence['unresolved']`; `meta['dtEffective']` and `convergence['h']` give the step actually taken, T/N.
The solver classes `SplitStepSolver` and `CaputoSolver` also validate their own parameters.

## Methods

| Model | Scheme |
|---|---|
| classical, space-fractional | Strang split-step Fourier, exact linear flow with ∫β, exact pointwise nonlinear flow with gain, Simpson-integrated coefficients |
| Caputo | L1 scheme on the graded mesh t_n = T(n/N)^r, r = min(10, (2−ν)/ν); implicit dispersion and gain; Anderson-accelerated fixed-point iteration (depth 10, relative tolerance 10⁻¹²) for the nonlinear, potential and forcing terms; L1 weights in the cancellation-free form w^{1−ν} expm1((1−ν) log1p(τ_k/w)) |
| fractional solitons | Petviashvili iteration for (μ + (−∂²)^{α/2})φ = 2φ³ |
| exact solutions | gauge reduction q = e^{G}p(x,τ), G = ∫g, τ = ∫β, γ = βe^{−2G} (V = 0); manufactured solution for Caputo |

For the Caputo model with ν = 1 the run uses the split-step solver with dt = T/N, and a warning is issued; a
`dt` passed together with `N` is ignored, with a warning. Used directly as a library class, `CaputoSolver`
also accepts ν = 1: the L1 weights then reduce to backward Euler on the graded mesh (first order and
dissipative), which is useful as a limit check but is not what the web interface or `Simulation` run.

Exact mode reports a reference floor (`exact.refFloorParts`: periodic images, spatial truncation, profile
residual, round-off). The image and spatial parts are estimates, not bounds, so the convergence study also
tests its own errors a posteriori: an observed order is withheld (`convergence.saturated`) when three
consecutive levels (step ratios need not be equal) do not fit e = C h^p + S with p > 0 and a step-independent
part S below 10 % of the error, and (`convergence.preasymptotic`) when the error exceeds 10 % of the peak
modulus; a study with only two levels publishes no order (`convergence.untested`). When the output times
under-sample β, γ or g (a coefficient differs by more than 5 % from its linear interpolation between outputs),
the exact-solution error is measured at every step (`exact.errorSampling = "every step"`). Split-step runs
compare the quadrature that the solver applies (Simpson on half steps for β and g; one gain-weighted Simpson
panel per step for the nonlinear phase; one Simpson panel per step for a time-dependent V, at 16 grid points
and the maxima of |V|) with the reference integrals (`meta.quadrature`); the reference
tables are themselves cross-checked against a 3-point Gauss–Legendre rule. These checks are estimates: a
coefficient that oscillates in phase with every sampling grid can still escape them.

## Tests and cross-check

```bash
pytest                                        # unit and convergence tests
python crosscheck.py ../../results/php_reference.json
```

`crosscheck.py` reruns the ten cases of `../php/tests/crosscheck_cases.json` and compares them with the
values exported by `php ../php/tests/export_reference.php ../../results/php_reference.json`. Each case
passes only if the component powers P₁ and P₂ (relative to the initial total power), the peak intensity, the
Hamiltonian, the final complex field (real and imaginary parts, hence also phases), the exact-solution error
history and the nonlinear iteration counts (identical at every step) all agree within the tolerances printed by the script.

## Licence

MIT, see `../../LICENSE.txt`.
