# Manakov Lab

Interactive laboratory for vector solitons of the **variable-coefficient Manakov system** and its **Riesz
space-fractional** and **Caputo time-fractional** variants. An experiment is set up in the browser, run on the
server by a numerical core written in plain PHP, and checked against exact solutions and automatic convergence
studies. A NumPy port of the same core cross-checks the PHP results and runs refinements beyond the web budget.

| Channel | Link |
|---|---|
| Online instance | https://fizmako.com/manakov/ |
| Source code (GitHub) | https://github.com/Cdacostam/ManakovLab |
| Archived release (Zenodo) | [doi:10.5281/zenodo.24406892](https://doi.org/10.5281/zenodo.24406892) |

Manakov Lab is described in an article submitted to *SoftwareX*.

## The model

The unknown is a two-component field $\mathbf q=(q_1,q_2)$, for example two polarizations of light in an
optical fibre or two species in a Bose–Einstein condensate. Both components share one cubic nonlinearity:

$$
i\,\partial_t q_j+\beta(t)\,\partial_x^2 q_j+2\gamma(t)\big(|q_1|^2+|q_2|^2\big)q_j+V(x,t)\,q_j=i\,g(t)\,q_j,\qquad j=1,2 .
$$

The coefficients are typed as formulas: $\beta(t)$ is the dispersion, $\gamma(t)$ the strength of the
nonlinearity, $g(t)$ a gain or loss and $V(x,t)$ an external potential. With constant coefficients and
$g=V=0$ the system is integrable.

| Variant | Change to the equation | Numerical method |
|---|---|---|
| Classical | the equation above | Strang split-step Fourier, second order in time |
| Space-fractional (Riesz) | $\partial_x^2 \to -(-\partial_x^2)^{\alpha/2}$, validated for $1<\alpha\le 2$ | split-step with symbol $\lvert k\rvert^\alpha$; solitons by Petviashvili iteration |
| Time-fractional (Caputo) | $i\,\partial_t \to i\,{}^{C}\!D_t^{\nu}$, $0<\nu\le 1$ | L1 scheme on a graded mesh, order $2-\nu$, Anderson-accelerated nonlinear solves |

What the laboratory offers:

- coefficients typed as formulas, validated and typeset as they are written; one to three initial solitons with
  their own amplitude, velocity, position and polarization; eight ready-made experiments;
- exact references: nonautonomous solitons by a gauge reduction, fractional solitons, manufactured Caputo
  solutions and a Mittag-Leffler mode of the linear Caputo problem;
- convergence studies that test their own errors and withhold an order, with the reason, when the error has
  reached the accuracy of the reference or has stopped decreasing;
- CSV and JSON export of profiles, diagnostics, the space–time field, the full result and the problem setup,
  and a printable run report.

## Repository layout

```
src/php/      web application: pure-PHP numerical core, API, home page and laboratory (deploy this folder)
src/python/   NumPy port of the same core (cross-check, figures, large refinements)
results/      archived outputs of the reproduction commands below
reproduce/    make_figures.py, which redraws the figures of the article (written to reproduce/figures/)
```

The two parts are self-contained: `src/php/` needs only PHP ≥ 8.1, `src/python/` needs only
Python ≥ 3.10 and NumPy (matplotlib for the figures). See `src/php/README.md` and `src/python/README.md`.

## Quick start

```bash
# web application: home page at http://localhost:8000, laboratory at http://localhost:8000/lab.php
cd src/php && php -S localhost:8000

# Python package
cd src/python && pip install -e ".[test,figures]" && pytest
```

## Reproducing the reported results

All commands run from the repository root. `J` enables the OPcache JIT; it only changes run times.

```bash
J="-d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=128M"
pip install -e "src/python[test,figures]"                                  # NumPy, pytest, matplotlib
php $J src/php/tests/run_tests.php                                         # 30 tests
php $J src/php/tests/reproduce.php all --json=results/reproduce.json > results/reproduce.txt
php $J src/php/tests/sweep.php --json=results/sweep.json > results/sweep.txt
php $J src/php/tests/reproduce.php timing > results/timing_jit.txt      # machine-dependent
php src/php/tests/reproduce.php timing > results/timing_nojit.txt        # same, without JIT
php $J src/php/tests/export_reference.php results/php_reference.json
python src/python/crosscheck.py results/php_reference.json > results/crosscheck.txt
python src/python/caputo_refinement.py > results/caputo_refinement.txt
python reproduce/make_figures.py results/reproduce.json                 # Figures 1-4 -> reproduce/figures/
```

What each command produces (table and section numbers refer to the article):

| Result | Command | Output |
|---|---|---|
| Convergence orders of the exact and manufactured references, Caputo memory study (N = 600, 1200) and every other number of Sections 3.1-3.3 (Tables 3-4); data of Figure 2 | `reproduce.php all` | `results/reproduce.txt`, `results/reproduce.json` |
| Sweep of the a posteriori checks of exact mode: configurations, levels, errors, published and withheld orders (Sec. 3.1) | `sweep.php` | `results/sweep.txt`, `results/sweep.json` |
| Levels N = 2400 and 4800 of caputo-decay for ν = 0.6-0.9, observed orders and Richardson limits (Sec. 3.3) | `caputo_refinement.py` | `results/caputo_refinement.txt` |
| Agreement between the PHP and Python implementations (Sec. 3.1) | `export_reference.php`, then `crosscheck.py` | `results/php_reference.json`, `results/crosscheck.txt` |
| Run times per preset, with and without JIT (Table 5; machine-dependent, archived values measured on the machine of Sec. 3.4) | `reproduce.php timing` | `results/timing_jit.txt`, `results/timing_nojit.txt` |
| Figures 1-4: architecture, convergence (from `reproduce.json`), collision and power under Caputo memory (computed with the Python port) | `make_figures.py` | `reproduce/figures/Figure_*.pdf`, `.png` |

## Authors

Alvaro H. Salas (corresponding author, ahsalass@unal.edu.co), Carlos D. Acosta M. and Elizabeth Solorzano Tovar,
Department of Mathematics and Statistics, Universidad Nacional de Colombia, Sede Manizales, Colombia.

## Citation

Please cite the archived software version, [doi:10.5281/zenodo.24406892](https://doi.org/10.5281/zenodo.24406892);
see `CITATION.cff`.

## Licence

MIT (see `LICENSE.txt`). Bundled third-party libraries are listed in `THIRD_PARTY_NOTICES.md`.
