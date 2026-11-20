# Manakov Lab — PHP web application

Interactive web laboratory for the **variable-coefficient Manakov system** and its
**Riesz space-fractional** and **Caputo time-fractional** variants — written in **pure PHP**.

- Numerical core in PHP 8: radix-2 FFT, Strang split-step Fourier, L1 scheme on graded mesh,
  Petviashvili iteration, Lanczos Γ, Simpson integration — no Python, no extensions, no shell calls
- Formula input parsed by a PHP recursive-descent parser → validated AST → LaTeX (KaTeX) and compiled closures
- KaTeX and Plotly (both bundled, no CDN needed) for the mathematics and the interactive 3D and 2D charts
- Exact-solution mode and automatic convergence studies for every model

```
i ∂t q + β(t) ∂x² q + 2γ(t)(|q₁|²+|q₂|²) q + V(x,t) q = i g(t) q          classical
i ∂t q − β(t)(−∂x²)^{α/2} q + 2γ(t)(|q₁|²+|q₂|²) q + V q = i g(t) q        space-fractional, tested at α = 1.2, 1.5, 2; recommended 1 < α ≤ 2
i ᶜDₜ^ν q + β(t) ∂x² q + 2γ(t)(|q₁|²+|q₂|²) q + V q = i g(t) q            Caputo, 0 < ν < 1 (α = 2)
```

## Requirements

- PHP **8.1+** (only the standard `json` extension; `zlib` recommended for gzip responses)
- Recommended: OPcache JIT (≈ 4× faster solvers)
- Any current browser (Plotly 2.35.2 is bundled in `assets/vendor/plotly/`)

## Quick start

The web application is the folder `src/php/` of the repository; `index.php` in it is the home page and `lab.php`
the laboratory.

**XAMPP / WAMP / MAMP / any Apache + PHP:** copy the contents of `src/php/` into a folder
`manakov-lab` inside `htdocs` (or `www`) and open **http://localhost/manakov-lab/**

**PHP built-in server:**

```bash
cd src/php
php -S localhost:8000
# open http://localhost:8000
```

With the JIT enabled (recommended, about 3x faster):

```bash
php -d opcache.enable=1 -d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=128M -S localhost:8000
```

## Production

Apache: the included `.htaccess` files block `config.php`, `src/` and `tests/` and disable directory listings
(needs `AllowOverride All`). Nginx: deny the same paths:

```nginx
server {
    root /var/www/manakov-lab;
    index index.php;
    location ^~ /src/ { deny all; }
    location ^~ /tests/ { deny all; }
    location = /config.php { deny all; }
    location ~ \.php$ {
        include fastcgi_params;
        fastcgi_param SCRIPT_FILENAME $document_root$fastcgi_script_name;
        fastcgi_pass unix:/run/php/php8.3-fpm.sock;
        fastcgi_read_timeout 200s;
    }
}
```

Recommended `php.ini`:

```ini
opcache.enable=1
opcache.jit=tracing
opcache.jit_buffer_size=128M
max_execution_time=180
memory_limit=512M
```

**Docker** (JIT enabled, run from `src/php/`): `docker build -t manakov-lab . && docker run -p 8080:80 manakov-lab`

## Configuration (`config.php`, each key overridable by environment variable)

| Variable | Default | Meaning |
|---|---|---|
| `MANAKOV_TIMEOUT` | `150` | wall-clock limit per simulation (s), checked inside the time loops |
| `MANAKOV_MAX_CONCURRENT` | `2` | simultaneous simulations (file-lock slots) |
| `MANAKOV_MAX_WORK_SPLIT` | `4e8` | budget `Nx·log2(Nx)·T/dt` (≈ 3.5·10⁶/s without JIT) |
| `MANAKOV_MAX_WORK_CAPUTO` | `8e8` | budget `N²·Nx` (≈ 10⁷/s without JIT) |
| `MANAKOV_MEMORY_LIMIT` | `512M` | memory limit for simulation requests |

Raise the budgets when the JIT is on. `MANAKOV_REPOSITORY_URL` and `MANAKOV_CITATION_URL` set the links to the
source code and the citation file shown in the page footers, and `MANAKOV_DOI` the Zenodo DOI shown on the home
page (links that still contain the placeholder `xxxx` are shown there as pending). The bundled `assets/vendor/plotly/plotly.min.js`
is used when present; otherwise the page falls back to the CDN.

## Using the app

The home page links to the laboratory (**Open the laboratory**, or `lab.php` directly).

1. Pick a **preset** (8 ready experiments) or set up your own problem.
2. Choose the **model** and type the coefficients as formulas, e.g. `1 + 0.5cos(2t)`,
   `0.05 sin(t) exp(-t/10)`, `-0.01x^2`. The PHP parser shows them as LaTeX while you type.
   Allowed: `t` (and `x` in V), `pi`, `e`, `+ - * / ^`,
   `sin cos tan sinh cosh tanh sech exp log sqrt abs atan sign heaviside`, implicit products.
3. Set the grid, time stepping and up to three bright vector solitons
   (amplitude 2η, velocity parameter ξ, centre x₀, polarization c ∈ ℂ²).
4. Optional: **Exact-solution mode** and **Convergence study**.
5. **Run** (or Ctrl/⌘ + Enter).

Views: rotatable 3D surface and space–time maps of |q₁|², |q₂|², total; final-time moduli with exact
overlay, log-scale intensity, real/imaginary parts; power, power balance, peak intensity, Hamiltonian;
error versus time; log–log convergence with reference slope and observed orders; run report.
PNG export on every chart, CSV/JSON data export, printable report, dark mode.

## Numerical methods (all in `src/Numerics/`)

| File | Content |
|---|---|
| `FFT.php` | iterative radix-2 Cooley–Tukey FFT, cached plans, NumPy conventions |
| `SplitStepSolver.php` | Strang splitting: exact linear flow `exp(−i|k|^α∫β)` in Fourier space, exact pointwise nonlinear flow with gain, Simpson-integrated coefficients, fused half-steps |
| `CaputoSolver.php` | L1 scheme on graded mesh `tₙ = T(n/N)^r`, implicit dispersion/gain, Anderson-accelerated fixed-point iteration (depth 10, relative tolerance 1e-12) for nonlinear + potential + forcing, cancellation-free L1 weights, full memory |
| `Solitons.php` | Manakov bright solitons, fractional solitons by Petviashvili iteration, Fourier shifts |
| `ExprCompiler.php` | vectorised closures from the validated AST (constant folding, IEEE `fdiv`) |
| `Integral.php`, `Gamma.php` | Simpson rule, tabulated antiderivatives, Lanczos Γ |
| `Simulation.php` | initial data, exact solutions, diagnostics, convergence study, JSON result |

### Validation (pure PHP, reproduced by the presets)

| Preset | max error | observed order | expected |
|---|---|---|---|
| Classical, β = 1+0.5cos2t, g = 0.05 sin t | 8.08·10⁻⁴ (Δt = 0.005) | 1.99, 2.00, 2.00 | 2 |
| Space-fractional α = 1.5 | 5.08·10⁻⁴ | 1.99, 2.00, 2.00 | 2 |
| Caputo ν = 0.6, all coefficients varying | 4.35·10⁻⁵ (N = 512) | 1.39, 1.39, 1.40 | 2 − ν = 1.4 |

Exact solutions: gauge reduction `q = e^{G} p(x, τ)`, `G = ∫g`, `τ = ∫β`, γ = β·e^{−2G} (classical and
space-fractional; V = 0), with p a Manakov soliton (α = 2) or a stationary Petviashvili fractional soliton
(α < 2); for Caputo a manufactured solution `q* = c(1 + t^ν + i t²)ψ(x)` with automatically computed forcing.

## API

`POST api/latex.php` with form fields `expr=1+0.5cos(2t)`, `vars=t` (`vars=xt` for V; GET with the same query parameters also works) →
`{"ok":true,"latex":"1 + 0.5\\,\\cos\\left(2\\,t\\right)","canonical":"(1+(0.5*cos((2*t))))","uses":["t"]}`

`POST api/run.php`:

```json
{
  "model": "caputo",
  "params": {"L": 160, "Nx": 512, "T": 15, "N": 600, "nu": 0.8, "r": "auto"},
  "coeffs": {"beta": "1 + 0.5cos(2t)", "gamma": "1", "g": "0", "V": "0"},
  "solitons": [{"eta": 0.6, "xi": 0.1, "x0": 5, "p1re": 1, "p1im": 0, "p2re": 0, "p2im": 0.7}],
  "exact": false, "convergence": false
}
```

Response: `result.surface`, `result.final`, `result.diagnostics`, `result.exact`, `result.convergence`,
`result.meta` (including `version`, `picardMean`, `picardMax`), `result.warnings`; errors as
`{"ok": false, "error": "...", "field": "..."}` with HTTP 4xx/5xx.

## Project layout

```
config.php
index.php          home page: model, variants, validation, citation (KaTeX)
lab.php            the laboratory (KaTeX + Plotly)
api/               run.php (simulation endpoint), latex.php (formula -> LaTeX)
assets/            css, js, img (home-page preview), bundled KaTeX and Plotly
src/               ExprParser.php, Presets.php, Http.php, ConcurrencyGate.php, bootstrap.php
src/Numerics/      FFT, solvers, solitons, expression compiler, integration, Γ, Simulation
tests/             run_tests.php, reproduce.php, export_reference.php, make_home_preview.php (command line only)
```

## Tests and reproduction of the article

```bash
J="-d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=128M"
php $J tests/run_tests.php                                        # 30 tests, about 20 s with the JIT
php $J tests/reproduce.php all --json=../../results/reproduce.json > ../../results/reproduce.txt   # Sec. 3
php $J tests/reproduce.php timing                                 # Table 5 (median of three runs)
php $J tests/export_reference.php ../../results/php_reference.json   # input of src/python/crosscheck.py
```

`reproduce.php` accepts section names instead of `all`: `table3`, `mittag`, `general`, `fractional`, `alpha`,
`table4`, `spatial`, `nu1`, `solver`, `collision`, `timing`.

## Notes

- Periodic boundaries: choose L large enough; for α < 2 tails decay algebraically.
- Exact mode reports a reference floor (`exact.refFloorParts`: periodic images, spatial truncation, profile
  residual, round-off); observed orders below ten times it are omitted. The image and spatial parts are
  estimates, not bounds, so the convergence study also fits e = C h^p + S to three consecutive levels (step
  ratios need not be equal) and withholds orders whose errors stop decreasing as a power of the step
  (`convergence.saturated`) or exceed 10 % of the peak modulus (`convergence.preasymptotic`). A study with only
  two levels (T/40 < dt <= T/20, or N < 34 for Caputo) cannot separate S and publishes no order
  (`convergence.untested`). `tests/sweep.php` runs this test over 636 admitted configurations (`results/sweep.txt`).
  The classical reference is the soliton periodized over the box, so a soliton may cross x = ±L/2. When the
  output times under-sample β, γ or g, the error is measured at every step (`exact.errorSampling`).
  Split-step runs warn when dt does not resolve β, g (Simpson on half steps), the nonlinear phase (one
  gain-weighted Simpson panel per step) or a time-dependent V(x,t) (one Simpson panel per step, checked at 16
  grid points and at the maxima of |V| at t = 0 and T) (`meta.quadrature`); `meta.dtEffective` is the step actually taken, T/N.
  These checks are estimates, not a certification: a coefficient that oscillates in phase with every
  sampling grid can escape them.
- Caputo keeps the full memory (work ∝ N²·Nx, memory ∝ N·Nx). For ν = 1 the run uses the split-step solver with dt = T/N
  and a warning; for ν ≥ 0.95 the L1 scheme is close to backward Euler and a dissipation warning is shown.
  `CaputoSolver` used directly as a library class also accepts ν = 1, where the L1 scheme is backward Euler.
- For α ≤ 1 the cubic fractional NLS is L²-critical or supercritical: runs carry a warning and exact mode is refused.
- The Hamiltonian is invariant only for constant β and γ, g = 0, time-independent V and ν = 1.
