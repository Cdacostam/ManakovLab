"""Experiment driver: a port of ``src/php/src/Numerics/Simulation.php``.

A *problem* is a dict with keys
    model ('classic' | 'space' | 'caputo'), L, Nx, T, dt (split-step) or N and r (Caputo), alpha, nu,
    beta, gamma, g (callables of t), V (callable of (x, t) or None),
    solitons (list of dicts eta, xi, x0, p1re, p1im, p2re, p2im), exact (bool), convergence (bool),
    warnings (optional list of notes).
The result mirrors the JSON returned by ``src/php/api/run.php`` (without rounding).
"""
from __future__ import annotations

import math
import time
from typing import Optional

import numpy as np

from . import __version__
import sys

from .core import (CaputoSolver, SimulationError, SplitStepSolver, antiderivative_adaptive, fourier_shift,
                   fractional_soliton, grid, manakov, manakov_periodic, nonlinear_quadrature_error, normalise,
                   panel_quadrature_error_sampled, php_round, relative_residual, step_quadrature_error)

N_SAVE = 400
REF_REFUSE = 1e-3        # relative overlap with the periodic images above which the classical reference is refused
FLOOR_REPORT = 1e-10     # absolute reference defect below which nothing is reported
PROFILE_TOL = 1e-9       # relative residual above which a fractional profile is not treated as stationary
LOCALIZED = 1e-2         # edge/peak ratio above which a stationary solution is reported as not localized
QUAD_TOL = 1e-6          # relative discrepancy of the step quadrature above which dt does not resolve a coefficient
ROUNDOFF = 1e-12         # rounding level of an exact-mode comparison, relative to the peak of the initial field
SATURATION = 0.1         # step-independent fraction of the error above which an observed order is withheld
PREASYMPTOTIC = 0.1      # error, relative to the peak modulus, above which an order is not asymptotic
OUTPUT_ALIAS = 0.05      # relative interpolation defect of a coefficient between outputs above which they alias it


def _cosh(v: float) -> float:
    return math.inf if abs(v) > 700.0 else math.cosh(v)


def php_round_half(v: float) -> float:
    """PHP round(): half away from zero."""
    return float(math.floor(v + 0.5)) if v >= 0 else -float(math.floor(-v + 0.5))


def _finite_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)



def fit_power_law(ha: float, hb: float, hc: float, ea: float, eb: float, ec: float):
    """Fit e = C h^p + S to three levels ha > hb > hc with unequal ratios (see Simulation.php::fitPowerLaw).

    Returns (p, S), or None when the errors do not decrease or no p > 0 fits them.
    """
    d1, d2 = ea - eb, eb - ec
    if not (d1 > 0.0 and d2 > 0.0 and ha > hb > hc > 0.0):
        return None
    lA, lB = math.log(ha / hb), math.log(hb / hc)
    target = math.log(d1 / d2)

    def lnf(p):
        return p * lB + math.log(math.expm1(p * lA)) - math.log(math.expm1(p * lB))
    if target <= math.log(lA / lB):            # f(0+) = lA/lB: only p <= 0 fits
        return None
    lo, hi = 0.0, 30.0
    if lnf(hi) <= target:
        lo = hi
    else:
        k = 0
        while k < 200 and hi - lo > 1e-13 * hi:
            mid = 0.5 * (lo + hi)
            if lnf(mid) < target:
                lo = mid
            else:
                hi = mid
            k += 1
    p = 0.5 * (lo + hi)
    C = d2 / (hb ** p - hc ** p)
    return p, ec - C * hc ** p

class Simulation:
    def __init__(self, P: dict):
        self._t_created = time.perf_counter()     # runtime_total includes the reference tables built here
        self.P = dict(P)
        self.P.setdefault('alpha', 2.0)
        self.P.setdefault('nu', 1.0)
        self.P.setdefault('V', None)
        self.P.setdefault('exact', False)
        self.P.setdefault('convergence', False)
        self.warnings = list(P.get('warnings', []))
        self.floor = {'images': 0.0, 'spatial': 0.0, 'profile': 0.0, 'roundoff': 0.0}
        self.profile_residual = None
        self.profile_info = None
        self._G_gamma = None
        self._quad_cache: dict = {}
        self._domain_checks()
        self.x, self.k, self.dx = grid(self.P['L'], self.P['Nx'])
        self.beta, self.gamma, self.g, self.V = self.P['beta'], self.P['gamma'], self.P['g'], self.P['V']
        self._check_coefficients()
        self._V_sample, self._V_points = self._potential_samples()
        self.G, self.G_info = antiderivative_adaptive(self.g, self.P['T'])
        self.tau, self.tau_info = antiderivative_adaptive(self.beta, self.P['T'])
        for name, info in (('g', self.G_info), ('beta', self.tau_info)):
            if not info['resolved']:
                self.warnings.append(f"{name}(t) varies too fast on [0, T] for the reference quadrature ({info['m']} "
                                     'Simpson panels): its integral, used by the exact solutions and the gain '
                                     'diagnostic, may be inaccurate.')

    @property
    def ref_floor(self) -> float:
        """Total reference floor (absolute)."""
        return sum(self.floor.values())

    # ------------------------------------------------------------------ checks
    def _mathematical_domain(self) -> None:
        """Mathematical domain of every model, as in the PHP core; the web API adds work-budget limits."""
        P = self.P
        pos = lambda v: _finite_number(v) and v > 0
        if P.get('model') not in ('classic', 'space', 'caputo'):
            raise SimulationError("model must be 'classic', 'space' or 'caputo'.", 'model')
        if not pos(P.get('L')):
            raise SimulationError('L must be a positive finite number.', 'L')
        Nx = P.get('Nx')
        if not isinstance(Nx, int) or isinstance(Nx, bool) or Nx < 4 or Nx & (Nx - 1):
            raise SimulationError('Nx must be a power of two, at least 4 (radix-2 FFT).', 'Nx')
        if not pos(P.get('T')):
            raise SimulationError('T must be a positive finite number.', 'T')
        if P['model'] == 'caputo':
            nu = P.get('nu')
            if not _finite_number(nu) or nu <= 0.0 or nu > 1.0:
                raise SimulationError('The Caputo order must satisfy 0 < nu <= 1.', 'nu')
            N = P.get('N')
            if not isinstance(N, int) or isinstance(N, bool) or N < 1:
                raise SimulationError('N (number of time steps) must be a positive integer.', 'N')
            if nu < 1.0 and (not pos(P.get('r')) or P['r'] < 1.0):
                raise SimulationError('The grading exponent must satisfy r >= 1.', 'r')
        else:
            dt = P.get('dt')
            if not pos(dt) or dt > P['T']:
                raise SimulationError('dt must satisfy 0 < dt <= T.', 'dt')
        if P['model'] == 'space':
            a = P.get('alpha')
            if not _finite_number(a) or a <= 0.0 or a > 2.0:
                raise SimulationError('The fractional order must satisfy 0 < alpha <= 2.', 'alpha')
        if not P.get('solitons'):
            raise SimulationError('At least one soliton is required.', 'solitons')
        for s in P['solitons']:
            if not pos(s.get('eta')):
                raise SimulationError('Every soliton needs eta > 0.', 'eta')
            for f in ('xi', 'x0', 'p1re', 'p1im', 'p2re', 'p2im'):
                if not _finite_number(s.get(f)):
                    raise SimulationError(f'Soliton parameter {f} must be a finite number.', f)

    def _domain_checks(self) -> None:
        P = self.P
        self._mathematical_domain()
        if P['model'] == 'space' and P['alpha'] <= 1.0:
            if P['exact']:
                raise SimulationError('Exact mode is available only for 1 < alpha <= 2: for alpha <= 1 the cubic '
                                      'fractional NLS is L2-critical (alpha = 1) or supercritical, and fractional '
                                      'solitons may be unstable.', 'alpha')
            self.warnings.append('alpha <= 1: the cubic fractional NLS is L2-critical (alpha = 1) or supercritical '
                                 '(alpha < 1); solitons may be unstable and solutions may collapse. Results are '
                                 'validated only for 1 < alpha <= 2.')
        if P['model'] == 'caputo':
            if P['nu'] >= 1.0:
                if P['exact']:
                    raise SimulationError('Caputo exact mode requires nu < 1. For nu = 1 use the Classical model.', 'nu')
                self.warnings.append('nu = 1: the Caputo derivative is the ordinary derivative, so the run uses the '
                                     'conservative split-step solver with dt = T/N (the L1 scheme would reduce to the '
                                     'dissipative backward Euler method).')
                dt = P.get('dt')
                if dt is not None and abs(dt - P['T'] / P['N']) > 1e-12 * P['T']:
                    self.warnings.append(f"nu = 1: the step is set by N; dt = {dt:.6g} was ignored and dt = T/N = "
                                         f"{P['T'] / P['N']:.6g} used.")
            elif P['nu'] >= 0.95:
                self.warnings.append('nu close to 1: the L1 scheme is close to backward Euler and adds numerical '
                                     'dissipation; check P(T)/P(0) under refinement of N before interpreting the power loss.')

    def _check_coefficients(self) -> None:
        T = self.P['T']
        for name, f in (('beta', self.beta), ('gamma', self.gamma), ('g', self.g)):
            for i in range(41):
                if not math.isfinite(f(T * i / 40)):
                    raise SimulationError(f'{name}(t) is not finite on [0, T].', name)
        if self.V is not None:
            for i in range(11):
                if not np.all(np.isfinite(self.V(self.x, T * i / 10))):
                    raise SimulationError('V(x,t) is not finite on the domain.', 'V')

    # ------------------------------------------------------------------ setup
    def _fractional_profile(self, eta: float, x0: float, snap: bool = False) -> np.ndarray:
        """Fractional profile centred at x0 (see Simulation.php): described (branch, convergence), residual of
        the profile actually used, off-node centres moved to a node in exact mode, residual floor."""
        alpha, mu = self.P['alpha'], 4.0 * eta * eta
        maxit = 3000
        phi, _, iters = fractional_soliton(self.x, self.k, alpha, mu, 1e-13, maxit)
        self._describe_profile(phi, None if iters > maxit else iters, eta)
        prof = fourier_shift(phi, self.k, x0)
        rel = relative_residual(prof, self.k, alpha, mu)
        m = php_round_half(x0 / self.dx)
        if rel > PROFILE_TOL and abs(x0 - m * self.dx) > 1e-12 * max(1.0, abs(x0)):
            if snap:
                x0n = m * self.dx
                self.warnings.append(f'x0 = {x0:.6g} lies between grid nodes: the Fourier-shifted fractional profile is not '
                                     f'a stationary solution on this grid (relative residual {rel:.1e}), and exact mode '
                                     f'would report that defect as time-integration error. The soliton was moved to the '
                                     f'nearest node, x0 = {x0n:.6g}.')
                prof = fourier_shift(phi, self.k, x0n)
                rel = relative_residual(prof, self.k, alpha, mu)
            else:
                self.warnings.append(f'The fractional soliton at x0 = {x0:.6g} lies between grid nodes and its shifted '
                                     f'profile is not stationary on this grid (relative residual {rel:.1e}). Increase Nx '
                                     f'or place x0 on a node (a multiple of dx = {self.dx:.6g}).')
        if snap and rel > PROFILE_TOL:
            # the defect r of the stationary equation acts as a forcing: |error| <~ |r| |tau| e^G
            c3 = float(np.max(np.abs(2.0 * prof ** 3)))
            self.floor['profile'] = rel * c3 * max(1.0, self.tau_info['maxAbs']) * math.exp(self.G_info['max'])
            how = 'did not converge' if self.profile_info['iterations'] is None else 'converged to its tolerance'
            self.warnings.append(f"Exact mode: the fractional profile satisfies its stationary equation only to a relative "
                                 f"residual {rel:.1e} (Petviashvili {how}). This defect alone can produce errors of about "
                                 f"{self.floor['profile']:.1e}, which are not attributed to the time integrator: the "
                                 f"corresponding observed orders are omitted.")
        self.profile_residual = max(self.profile_residual or 0.0, rel)
        return prof

    def _describe_profile(self, phi: np.ndarray, iters: Optional[int], eta: float) -> None:
        """Identifies the stationary solution computed on the box and warns if it is not a localized soliton."""
        a = np.abs(phi)
        mx, mn = float(np.max(a)), float(np.min(a))
        edge = abs(float(phi[0])) / mx if mx > 0.0 else 1.0          # phi[0] is at x = -L/2
        branch = 'homogeneous' if mx > 0.0 and mn / mx > 1.0 - 1e-6 else ('periodic' if edge > LOCALIZED else 'localized')
        self.profile_info = {'branch': branch, 'edge': edge, 'iterations': iters, 'amplitude': mx}
        if iters is None:
            self.warnings.append('Petviashvili iteration did not converge in 3000 iterations: the fractional profile is '
                                 'not a converged stationary solution on this grid.')
        if branch == 'homogeneous':
            self.warnings.append(f'The stationary solution found on this box is spatially constant, phi = sqrt(mu/2) = '
                                 f'{math.sqrt(2.0 * eta * eta):.6g} (the homogeneous branch), not a soliton: the box is '
                                 f'too small for a localized state with eta = {eta:.4g}. Exact mode then tests the time '
                                 f'integration of a plane wave. Increase L.')
        elif branch == 'periodic':
            self.warnings.append(f'The stationary solution found on this box is not localized: its value at x = +-L/2 is '
                                 f'{edge:.2g} of the peak, so it is a periodic stationary state of the box rather than an '
                                 f'isolated soliton. Increase L for a localized soliton.')

    def _check_classical_reference(self, eta: float, xi: float) -> None:
        """Periodized classical reference: refusal when the soliton overlaps its images, images and spatial floors
        (see Simulation.php::checkClassicalReference)."""
        L = self.P['L']
        eG = math.exp(self.G_info['max'])
        overlap = 2.0 / _cosh(eta * L)
        if overlap > REF_REFUSE:
            raise SimulationError(f'Exact mode: the soliton (eta = {eta:.4g}) is too wide for the periodic box: midway to '
                                  f'its periodic image its tail is {overlap:.1e} of the peak, so it interacts with its '
                                  f'images and is no longer the reference solution. Use L >= '
                                  f'{math.ceil(math.log(4.0 / REF_REFUSE) / eta):.4g}.', 'L')
        growth = 96.0 * eta ** 3 * max(1.0, self.tau_info['maxAbs']) * eG
        self.floor['images'] = growth / _cosh(2.0 * eta * L)
        kN = math.pi * self.P['Nx'] / L
        d = kN - 2.0 * abs(xi)
        if d <= 0.0:
            raise SimulationError(f'Exact mode: the carrier wavenumber 2|xi| = {2.0 * abs(xi):.4g} exceeds the largest '
                                  f'grid wavenumber pi Nx/L = {kN:.4g}. Increase Nx.', 'Nx')
        self.floor['spatial'] = 4.0 * eta / math.pi * eG * math.exp(-math.pi * d / (4.0 * eta))
        if self.floor['images'] > FLOOR_REPORT:
            Lmin = math.ceil(math.log(2.0 * growth / FLOOR_REPORT) / (2.0 * eta))
            self.warnings.append(f"Exact mode: the soliton interacts weakly with its periodic images (tail {overlap:.1e} of "
                                 f"the peak midway); the periodized reference is estimated to be accurate to about "
                                 f"{self.floor['images']:.1e} (L >= {Lmin:.4g} makes this negligible). Observed orders of "
                                 f"errors below ten times this value are omitted.")
        if self.floor['spatial'] > FLOOR_REPORT:
            need = 2.0 * abs(xi) + 4.0 * eta / math.pi * math.log(4.0 * eta / math.pi * eG / FLOOR_REPORT)
            nx_min = 2 ** int(math.ceil(math.log2(need * L / math.pi)))
            self.warnings.append(f"Exact mode: the grid does not resolve the soliton spectrum (dx = {self.dx:.4g} for a "
                                 f"soliton of width 1/(2 eta) = {0.5 / eta:.4g}). The truncation of the initial soliton "
                                 f"alone is estimated at {self.floor['spatial']:.1e}; the spatial error can grow with "
                                 f"time beyond this estimate. Observed orders of errors below ten times it are omitted, "
                                 f"and the convergence study also omits orders whose errors stop decreasing with dt. "
                                 f"Nx >= {nx_min} makes it negligible.")

    def initial_condition(self) -> np.ndarray:
        Q = np.zeros((2, self.P['Nx']), dtype=complex)
        for s in self.P['solitons']:
            c = (s['p1re'] + 1j * s['p1im'], s['p2re'] + 1j * s['p2im'])
            if self.P['model'] == 'space' and self.P['alpha'] < 2.0:
                prof = self._fractional_profile(s['eta'], s['x0'])
                Q += normalise(c)[:, None] * (prof * np.exp(-2j * s['xi'] * self.x))[None, :]
            else:
                Q += manakov(self.x, 0.0, s['eta'], s['xi'], c, s['x0'])
        return Q

    def _setup_exact(self):
        P = self.P
        if len(P['solitons']) > 1:
            self.warnings.append('Exact mode uses only the first soliton.')
        s = P['solitons'][0]
        eta, xi, x0 = s['eta'], s['xi'], s['x0']
        c = normalise((s['p1re'] + 1j * s['p1im'], s['p2re'] + 1j * s['p2im']))
        x, G, tau, beta = self.x, self.G, self.tau, self.beta

        if P['model'] != 'caputo':
            if self.V is not None:
                self.warnings.append('Exact mode requires V = 0; the potential was switched off.')
                self.V = None
            gamma_ex = lambda t: beta(t) * math.exp(-2.0 * G(t))
            if P['model'] == 'space' and P['alpha'] < 2.0:
                if xi != 0.0:
                    self.warnings.append('For alpha < 2 the exact soliton is stationary: xi was set to 0.')
                mu = 4.0 * eta * eta
                prof = self._fractional_profile(eta, x0, snap=True)
                exact = lambda t: c[:, None] * (math.exp(G(t)) * prof * np.exp(1j * mu * tau(t)))[None, :]
            else:
                self._check_classical_reference(eta, xi)
                L = P['L']
                exact = lambda t: manakov_periodic(x, L, tau(t), eta, xi, c, x0, math.exp(G(t)))
            return exact, gamma_ex, None

        # Caputo: manufactured solution q* = c (1 + t^nu + i t^2) psi(x)
        nu = P['nu']
        psi = 2.0 * eta / np.cosh(2.0 * eta * (x - x0)) * np.exp(-2j * xi * x)
        psixx = np.fft.ifft(-self.k ** 2 * np.fft.fft(psi))
        g1, g3 = math.gamma(1.0 + nu), math.gamma(3.0 - nu)
        A = lambda t: (1.0 + t ** nu) + 1j * t * t
        dA = lambda t: g1 + 1j * 2.0 * t ** (2.0 - nu) / g3
        exact = lambda t: c[:, None] * (A(t) * psi)[None, :]
        gam, gf, V = self.gamma, self.g, self.V

        def forcing(t: float) -> np.ndarray:
            q = exact(t)
            w = 2.0 * gam(t) * np.sum(np.abs(q) ** 2, axis=0)
            if V is not None:
                w = w + V(x, t)
            return (1j * c[:, None] * (dA(t) * psi)[None, :] + beta(t) * c[:, None] * (A(t) * psixx)[None, :]
                    + w[None, :] * q - 1j * gf(t) * q)

        self.warnings.append('Caputo exact mode uses a manufactured solution q* = c(1 + t^nu + i t^2) psi(x) '
                             'with a computed forcing term.')
        return exact, None, forcing

    def _make_solver(self, gamma_override=None, forcing=None):
        P = self.P
        gam = gamma_override or self.gamma
        if P['model'] == 'caputo' and P['nu'] >= 1.0:
            return SplitStepSolver(P['L'], P['Nx'], 2.0, self.beta, gam, self.g, self.V)
        if P['model'] == 'caputo':
            return CaputoSolver(P['L'], P['Nx'], P['nu'], self.beta, gam, self.g, self.V, forcing,
                                anderson=int(P.get('anderson', 10)))
        alpha = P['alpha'] if P['model'] == 'space' else 2.0
        return SplitStepSolver(P['L'], P['Nx'], alpha, self.beta, gam, self.g, self.V)

    def _run_solver(self, s, Q0, nSave, on_save, dt=None, N=None):
        P = self.P
        if isinstance(s, CaputoSolver):
            return s.solve(Q0, P['T'], N or P['N'], P['r'], nSave, on_save)
        if P['model'] == 'caputo':                  # nu = 1: the step is always T/N
            return s.solve(Q0, P['T'], P['T'] / (N or P['N']), nSave, on_save)
        return s.solve(Q0, P['T'], dt if dt is not None else P['dt'], nSave, on_save)

    def _steps_for(self, dt: float) -> int:
        """Number of steps actually taken by the split-step solver for a requested dt (it uses T/N)."""
        return max(1, php_round(self.P['T'] / dt))

    @staticmethod
    def expected_order(P: dict) -> float:
        """Expected order of the time discretization: 2 for Strang (also nu = 1), min(r nu, 2 - nu) for L1."""
        return min(P['r'] * P['nu'], 2.0 - P['nu']) if P['model'] == 'caputo' and P['nu'] < 1.0 else 2.0

    def _potential_samples(self):
        """Sample points of the quadrature diagnostic of a time-dependent V (see Simulation.php constructor).

        A V whose values at the sample points do not change at five probe times is treated as time-independent,
        for which the Simpson panel is exact.
        """
        V, x, T = self.V, self.x, self.P['T']
        if V is None:
            return None, 0
        Nx = self.P['Nx']
        n16 = min(16, Nx)
        idx = [j * Nx // n16 for j in range(n16)]
        for t in (0.0, float(T)):
            idx.append(int(np.argmax(np.abs(np.broadcast_to(V(x, t), x.shape)))))
        idx = list(dict.fromkeys(idx))
        xs = x[idx]

        def Vs(t):
            return np.broadcast_to(np.asarray(V(xs, t), dtype=float), xs.shape)
        v0 = Vs(0.0)
        if all(np.array_equal(Vs(f * T), v0) for f in (0.1234, 0.381966, 0.618034, 0.8765, 1.0)):
            return None, 0
        return Vs, len(idx)

    def _step_quadrature(self, N: int, gam_eff) -> dict:
        """Relative discrepancy between the split-step quadrature and the reference integral of beta, gamma, g."""
        if N in self._quad_cache:
            return self._quad_cache[N]
        if self._G_gamma is None:
            G = self.G
            self._G_gamma, _ = antiderivative_adaptive(lambda t: gam_eff(t) * math.exp(2.0 * G(t)), self.P['T'])
        T = self.P['T']
        q = {'beta': step_quadrature_error(self.beta, self.tau, T, N),
             'gamma': nonlinear_quadrature_error(gam_eff, self.g, self._G_gamma, self.G, T, N),
             'g': step_quadrature_error(self.g, self.G, T, N)}
        if self.V is not None and self._V_sample is not None:   # a time-independent V is integrated exactly
            q['V'] = panel_quadrature_error_sampled(self._V_sample, T, N)
        self._quad_cache[N] = q
        return q

    @staticmethod
    def max_diff(A: np.ndarray, B: np.ndarray) -> float:
        d = np.abs(A - B)
        return math.inf if np.any(np.isnan(d)) else float(np.max(d))

    # ------------------------------------------------------------------ run
    def run(self) -> dict:
        start = time.perf_counter()
        P = self.P
        model = P['model']
        exact = gamma_ex = forcing = None
        if P['exact']:
            exact, gamma_ex, forcing = self._setup_exact()
            Q0 = exact(0.0)
            # rounding level of the comparison: errors near it carry no information on the order
            self.floor['roundoff'] = ROUNDOFF * float(np.max(np.abs(Q0)))
        else:
            Q0 = self.initial_condition()
        solver = self._make_solver(gamma_ex, forcing)
        alpha = P['alpha'] if model == 'space' else 2.0
        Dk = np.abs(self.k) ** alpha
        Nx, dx = P['Nx'], self.dx
        beta, G, V = self.beta, self.G, self.V
        gam_eff = gamma_ex or self.gamma
        rec = {'t': [], 'P1': [], 'P2': [], 'peak': [], 'H': [], 'gain': [], 'err': []}
        state = {'initial': None, 'final': None}

        def on_save(t: float, Q: np.ndarray) -> None:
            a = np.abs(Q) ** 2
            I = a[0] + a[1]
            Qh = np.fft.fft(Q, axis=-1)
            kin = float(np.sum(Dk * (np.abs(Qh[0]) ** 2 + np.abs(Qh[1]) ** 2)))
            VI = float(np.sum(V(self.x, t) * I)) if V is not None else 0.0
            rec['t'].append(t)
            rec['P1'].append(float(np.sum(a[0])) * dx)
            rec['P2'].append(float(np.sum(a[1])) * dx)
            rec['peak'].append(float(np.max(I)))
            rec['H'].append(beta(t) * kin * dx / Nx - gam_eff(t) * float(np.sum(I * I)) * dx - VI * dx)
            rec['gain'].append(math.exp(2.0 * G(t)))
            if exact is not None:
                rec['err'].append(self.max_diff(Q, exact(t)))
            if state['initial'] is None:
                state['initial'] = Q.copy()
            state['final'] = Q.copy()

        quad = None
        N_main = 0
        if isinstance(solver, SplitStepSolver):
            N_main = P['N'] if model == 'caputo' else self._steps_for(P['dt'])
            quad = self._step_quadrature(N_main, gam_eff)
            for name, d in quad.items():
                if d > QUAD_TOL:
                    extra = ('; in particular the power follows exp(2 x step sum of g), not exp(2 int g)'
                             if name == 'g' else '')
                    how = {'gamma': 'one gain-weighted Simpson panel per step for the nonlinear phase',
                           'V': f'one Simpson panel per step, checked at {self._V_points} grid points'}.get(
                        name, 'Simpson on half steps')
                    label = 'V(x,t)' if name == 'V' else f'{name}(t)'
                    self.warnings.append(f"dt = {P['T'] / N_main:.4g} does not resolve {label}: the quadrature that the "
                                         f"solver applies ({how}) differs from the integral by {d:.1e} (relative). The run "
                                         f"then follows the discretized coefficient{extra}. Reduce dt.")
        self._run_solver(solver, Q0, N_SAVE, on_save)
        if not np.all(np.isfinite(state['final'])):
            raise SimulationError('The solution blew up (non-finite values). Reduce dt / increase N, or check the '
                                  'coefficients.', 'N' if model == 'caputo' else 'dt')
        t_main = time.perf_counter() - start
        res = {
            'model': model, 'x': self.x, 't': np.array(rec['t']), 'initial': state['initial'], 'final': state['final'],
            'diagnostics': {k: np.array(rec[k]) for k in ('P1', 'P2', 'peak', 'H', 'gain')},
            'warnings': self.warnings,
        }
        sampling = None
        if exact is not None:
            err_max = max(rec['err'])
            sampling = 'outputs'
            aliased = self._output_aliasing(N_main, gam_eff) if isinstance(solver, SplitStepSolver) else []
            if aliased:
                # the output times can fall where the error happens to vanish (e.g. at the periods of beta)
                sampling = 'every step'
                err_max = max(err_max, self._error_at_every_step(exact, gamma_ex, forcing, Q0))
                self.warnings.append(f"The output times (every {max(1, N_main // N_SAVE)} steps) under-sample "
                                     f"{', '.join(n + '(t)' for n in aliased)}, so the error between them is not "
                                     f"represented by the stored outputs. The exact error (and the convergence study) is "
                                     f"therefore measured at every time step.")
            res['exact'] = {'final': exact(rec['t'][-1]), 'err_t': np.array(rec['err']),
                            'err_final': rec['err'][-1], 'err_max': err_max, 'errorSampling': sampling,
                            'refFloor': self.ref_floor, 'refFloorParts': dict(self.floor)}
        if P['convergence'] and exact is not None:
            peak = max(float(np.max(np.abs(Q0))), float(np.max(np.abs(exact(rec['t'][-1])))))
            res['convergence'] = self._convergence(exact, gamma_ex, forcing, Q0, err_max, sampling == 'every step', peak)
        caputo = isinstance(solver, CaputoSolver)
        res['meta'] = {
            'version': __version__, 'Nx': Nx, 'L': P['L'], 'dx': dx, 'T': P['T'],
            'dt': P.get('dt'), 'N': P.get('N'), 'r': P.get('r'), 'alpha': alpha, 'nu': P['nu'],
            'dtEffective': (P['T'] / (P['N'] if model == 'caputo' else self._steps_for(P['dt']))
                            if isinstance(solver, SplitStepSolver) else None),
            'profileResidual': self.profile_residual,
            'profile': self.profile_info,
            'quadrature': quad,
            'picardMax': solver.max_picard if caputo else None,
            'picardMean': round(solver.sum_picard / solver.steps_done, 2) if caputo and solver.steps_done else None,
            'picardSteps': list(solver.picard_steps) if caputo else None,
            'runtime_main': t_main, 'runtime_total': time.perf_counter() - self._t_created,
            'expectedOrder': self.expected_order(P),
        }
        return res

    def _output_aliasing(self, N: int, gam_eff) -> list:
        """Coefficients that the output times (every k steps) under-sample (see Simulation.php::outputAliasing)."""
        k = max(1, N // N_SAVE)
        if k == 1:
            return []
        T = self.P['T']
        D = k * T / N
        th = 0.5 * (3.0 - math.sqrt(5.0))
        out = []
        for name, f in (('beta', self.beta), ('gamma', gam_eff), ('g', self.g)):
            d = 0.0
            scale = 1.0
            f0 = f(0.0)
            j = 0
            while (j + 1) * D <= T * (1.0 + 1e-12):
                f1 = f((j + 1) * D)
                d = max(d, abs(f(j * D + th * D) - ((1.0 - th) * f0 + th * f1)))
                scale = max(scale, abs(f0), abs(f1))
                f0 = f1
                j += 1
            if d > OUTPUT_ALIAS * scale:
                out.append(name)
        return out

    def _error_at_every_step(self, exact, gamma_ex, forcing, Q0, dt=None, N=None) -> float:
        s = self._make_solver(gamma_ex, forcing)
        e = [0.0]

        def cb(t, Q):
            e[0] = max(e[0], self.max_diff(Q, exact(t)))
        self._run_solver(s, Q0, sys.maxsize, cb, dt, N)
        return e[0]

    def _check_asymptotic(self, conv: dict, caputo: bool, peak: float) -> None:
        """A posteriori test of the observed orders (see Simulation.php::checkAsymptotic)."""
        e = conv['err']
        n = len(e)
        hs = [1.0 / h if caputo else float(h) for h in conv['h']]
        S_max = None
        for i in range(1, n):
            if conv['order'][i] is None:
                continue
            if e[i] > PREASYMPTOTIC * peak:
                conv['order'][i] = None
                conv['preasymptotic'].append(conv['h'][i])
                continue
            if n < 3:
                if e[i] >= e[i - 1]:              # the error does not decrease: nothing to attribute
                    conv['order'][i] = None
                    conv['saturated'].append(conv['h'][i])
                    S_max = max(S_max or 0.0, e[i])
                else:                             # decreasing, but S cannot be separated from C h^p
                    conv['order'][i] = None
                    conv['untested'].append(conv['h'][i])
                continue
            a, b, c = (i - 1, i, i + 1) if i + 1 < n else (i - 2, i - 1, i)
            sat = True
            if None not in (e[a], e[b], e[c]):
                fit = fit_power_law(hs[a], hs[b], hs[c], e[a], e[b], e[c])
                if fit is not None:
                    sat = abs(fit[1]) > SATURATION * e[i]
                    if sat:
                        S_max = max(S_max or 0.0, abs(fit[1]))
                else:
                    S_max = max(S_max or 0.0, e[c])
            if sat:
                conv['order'][i] = None
                conv['saturated'].append(conv['h'][i])
        conv['dtIndependent'] = S_max
        if conv['saturated']:
            about = f', about {S_max:.1e},' if S_max is not None else ''
            self.warnings.append(f"Convergence study: the errors do not decrease as a power of the time step at "
                                 f"{len(conv['saturated'])} level(s); a component that does not depend on the step{about} "
                                 f"dominates there (spatial resolution, the reference or round-off), so those observed "
                                 f"orders are omitted. Increase Nx (or L) and check again.")
        if conv['untested']:
            self.warnings.append("Convergence study: the study has only two levels, so the a posteriori saturation test, "
                                 "which needs three, cannot be applied and the observed order is omitted. "
                                 + ("Use N >= 34 for three levels or N >= 68 for four." if caputo else
                                    "Use dt <= T/40 for three levels or dt <= T/80 for four (levels need dt <= T/10)."))
        if conv['preasymptotic']:
            self.warnings.append(f"Convergence study: at {len(conv['preasymptotic'])} level(s) the error exceeds "
                                 f"{round(100 * PREASYMPTOTIC)}% of the peak amplitude, so the run is not in the "
                                 f"asymptotic regime and those observed orders are omitted. Reduce the time step.")

    def _convergence(self, exact, gamma_ex, forcing, Q0, err_main: float, dense: bool = False,
                     peak: float = math.inf) -> dict:
        P = self.P
        conv = {'h': [], 'err': [], 'order': [], 'refFloor': self.ref_floor, 'unresolved': [], 'saturated': [],
                'preasymptotic': [], 'untested': [], 'dtIndependent': None, 'expected': self.expected_order(P)}
        floor = 10.0 * self.ref_floor     # errors below this are not attributable to the integrator
        gam_eff = gamma_ex or self.gamma

        def error_of(dt=None, N=None) -> Optional[float]:
            s = self._make_solver(gamma_ex, forcing)
            e = [0.0]

            def cb(t, Q):
                e[0] = max(e[0], self.max_diff(Q, exact(t)))
            try:
                self._run_solver(s, Q0, sys.maxsize if dense else N_SAVE, cb, dt, N)
            except SimulationError:
                return None
            return e[0] if math.isfinite(e[0]) else None

        if P['model'] == 'caputo':
            N = P['N']
            Ns = sorted(set([max(16, N // 8), max(16, N // 4), max(16, N // 2), N]))
            for n in Ns:
                conv['h'].append(n)
                conv['err'].append(err_main if n == N else error_of(N=n))
            for i, h in enumerate(conv['h']):
                e0, e1 = (conv['err'][i - 1] if i else None), conv['err'][i]
                conv['order'].append(math.log(e0 / e1) / math.log(h / conv['h'][i - 1])
                                     if i and e0 and e1 and e1 > floor else None)
        else:
            dt = P['dt']
            dts = [d for d in (8 * dt, 4 * dt, 2 * dt, dt) if d <= P['T'] / 10] or [dt]
            bad = []
            for d in dts:
                n = self._steps_for(d)
                conv['h'].append(P['T'] / n)              # the step actually taken (T/N), not the requested one
                conv['err'].append(err_main if d == dt else error_of(dt=d))
                bad.append(max(self._step_quadrature(n, gam_eff).values()) > QUAD_TOL)
                if bad[-1]:
                    conv['unresolved'].append(P['T'] / n)
            for i, h in enumerate(conv['h']):
                e0, e1 = (conv['err'][i - 1] if i else None), conv['err'][i]
                ok = i and e0 and e1 and e1 > floor and not bad[i - 1] and not bad[i]
                conv['order'].append(math.log(e0 / e1) / math.log(conv['h'][i - 1] / h) if ok else None)
        self._check_asymptotic(conv, P['model'] == 'caputo', peak)
        return conv
