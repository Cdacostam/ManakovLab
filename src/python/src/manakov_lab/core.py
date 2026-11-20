"""Numerical core of Manakov Lab in Python: a line-by-line port of ``src/php/src/Numerics``.

Field layout: complex array of shape (2, Nx), q[j] = q_{j+1}(x).
Grid: x_i = (i - Nx/2) dx on the periodic domain [-L/2, L/2), angular wavenumbers in NumPy fftfreq order.
"""
from __future__ import annotations

import math
from typing import Callable, Optional

import numpy as np

TimeFn = Callable[[float], float]
FieldFn = Callable[[np.ndarray, float], np.ndarray]


class SimulationError(RuntimeError):
    """Raised when a run cannot be completed; ``field`` names the offending input, if any."""

    def __init__(self, message: str, field: Optional[str] = None):
        super().__init__(message)
        self.field = field


# --------------------------------------------------------------------------- grid, integrals
def grid(L: float, Nx: int) -> tuple[np.ndarray, np.ndarray, float]:
    if L <= 0 or Nx < 4 or Nx % 2:
        raise ValueError("L must be positive and Nx an even integer >= 4")
    dx = L / Nx
    x = (np.arange(Nx) - Nx // 2) * dx
    k = 2.0 * np.pi * np.fft.fftfreq(Nx, d=dx)
    return x, k, dx


def simpson(f: TimeFn, a: float, b: float) -> float:
    return (b - a) / 6.0 * (f(a) + 4.0 * f(0.5 * (a + b)) + f(b))


def _table(f: TimeFn, T: float, m: int) -> list:
    """Composite Simpson table F[i] = int_0^{i h} f, h = T/m."""
    h = T / m
    F = [0.0]
    for i in range(m):
        F.append(F[i] + simpson(f, i * h, (i + 1) * h))
    return F


def _closure(f: TimeFn, F: list, h: float, m: int) -> TimeFn:
    def G(t: float) -> float:
        if t <= 0.0:
            return 0.0
        i = min(m - 1, int(math.floor(t / h)))
        return F[i] + simpson(f, i * h, t)

    return G


def antiderivative(f: TimeFn, T: float, m: int = 4096) -> TimeFn:
    """F(t) = int_0^t f, tabulated by composite Simpson on m panels and completed by one local panel."""
    return _closure(f, _table(f, T, m), T / m, m)


def _agrees_with_gauss(f: TimeFn, T: float, m: int, F: list, tol: float) -> bool:
    h = T / m
    r = 0.5 * math.sqrt(0.6)
    S = 0.0
    d = 0.0
    scale = 1.0
    for i in range(m):
        c = (i + 0.5) * h
        S += h / 18.0 * (5.0 * f(c - r * h) + 8.0 * f(c) + 5.0 * f(c + r * h))
        d = max(d, abs(S - F[i + 1]))
        scale = max(scale, abs(S))
    return d <= tol * scale


def antiderivative_adaptive(f: TimeFn, T: float, tol: float = 1e-8, m0: int = 4096,
                            m_max: int = 1048576) -> tuple[TimeFn, dict]:
    """Antiderivative with a controlled table (port of Integral::antiderivativeAdaptive).

    The number of panels is doubled from m0 until two successive tables agree at their common nodes to
    tol * max(1, max|F|), or until m_max; the coarser of the two agreeing tables is used, and it must also agree
    with a composite 3-point Gauss-Legendre table on the same panels (nodes not commensurate with Simpson's).
    Both tests estimate the quadrature error; they do not certify every formula.
    """
    m = m0
    F = _table(f, T, m)
    resolved = False
    while 2 * m <= m_max:
        F2 = _table(f, T, 2 * m)
        d = 0.0
        scale = 1.0
        for i, v in enumerate(F):
            d = max(d, abs(v - F2[2 * i]))
            scale = max(scale, abs(F2[2 * i]))
        if d <= tol * scale and _agrees_with_gauss(f, T, m, F, tol):
            resolved = True
            break
        m *= 2
        F = F2
    info = {"m": m, "resolved": resolved, "max": max(F), "min": min(F), "maxAbs": max(abs(max(F)), abs(min(F)))}
    return _closure(f, F, T / m, m), info


def nonlinear_quadrature_error(gamma: TimeFn, g: TimeFn, Gam: TimeFn, G: TimeFn, T: float, N: int) -> float:
    """Discrepancy of the nonlinear phase quadrature of a split-step run (port of
    Integral::nonlinearQuadratureError): the gain-weighted Simpson step of SplitStepSolver._nonlinear, weighted by
    e^{2G(t0)}, against Gam(t) = int_0^t gamma e^{2G}, relative to max(1, max|Gam|)."""
    dt = T / N
    S = 0.0
    d = 0.0
    scale = 1.0
    for n in range(1, N + 1):
        t0 = (n - 1) * dt
        tm = t0 + 0.5 * dt
        t1 = n * dt
        Gm = simpson(g, t0, tm)
        G1 = Gm + simpson(g, tm, t1)
        step = dt / 6.0 * (gamma(t0) + 4.0 * gamma(tm) * math.exp(2.0 * Gm) + gamma(t1) * math.exp(2.0 * G1))
        S += math.exp(2.0 * G(t0)) * step
        ref = Gam(t1)
        d = max(d, abs(S - ref))
        scale = max(scale, abs(ref))
    return d / scale


def step_quadrature_error(f: TimeFn, F: TimeFn, T: float, N: int) -> float:
    """Largest discrepancy, relative to max(1, max|F|), between the split-step quadrature (Simpson on half
    steps) and the reference antiderivative F at the step nodes t_n = n T/N."""
    dt = T / N
    S = 0.0
    d = 0.0
    scale = 1.0
    f0 = f(0.0)
    for n in range(1, N + 1):
        t0 = (n - 1) * dt
        tm = t0 + 0.5 * dt
        t1 = n * dt
        # Simpson on [t0, tm] and [tm, t1], as in SplitStepSolver, with the shared endpoint values reused
        fm = f(tm)
        f1 = f(t1)
        S += dt / 12.0 * (f0 + 4.0 * f(0.5 * (t0 + tm)) + 2.0 * fm + 4.0 * f(0.5 * (tm + t1)) + f1)
        f0 = f1
        ref = F(t1)
        d = max(d, abs(S - ref))
        scale = max(scale, abs(ref))
    return d / scale


def panel_quadrature_error_sampled(Vs, T: float, N: int, tol: float = 1e-8, m0: int = 4096,
                                   m_max: int = 1048576) -> float:
    """Discrepancy of the potential quadrature at sample points (see Integral::panelQuadratureErrorSampled).

    Vs(t) returns the values of V at the K sample points. Each step adds h/6 [V(t0) + 4 V(tm) + V(t1)]; the
    reference integrates every step with M Simpson panels, M doubled from ceil(m0/N) until it agrees with 2M
    panels and with three-point Gauss-Legendre on M panels to tol * scale. Returns max |sum - ref| / scale.
    """
    dt = T / N
    M = max(2, math.ceil(m0 / N))
    r = 0.5 * math.sqrt(0.6)
    while True:
        hq, hp = dt / (4 * M), dt / M
        prev = np.asarray(Vs(0.0), dtype=float)
        AM = np.zeros_like(prev); A2 = np.zeros_like(prev); AG = np.zeros_like(prev); S = np.zeros_like(prev)
        d_ref = d_gl = d_sol = d_sol2 = 0.0
        scale = 1.0
        for n in range(1, N + 1):
            t0 = (n - 1) * dt
            v = [prev] + [np.asarray(Vs(t0 + j * hq), dtype=float) for j in range(1, 4 * M + 1)]
            prev = v[4 * M]
            iM = sum(v[4 * p] + 4.0 * v[4 * p + 2] + v[4 * p + 4] for p in range(M))
            i2 = sum(v[2 * p] + 4.0 * v[2 * p + 1] + v[2 * p + 2] for p in range(2 * M))
            AM = AM + hp / 6.0 * iM
            A2 = A2 + hp / 12.0 * i2
            S = S + dt / 6.0 * (v[0] + 4.0 * v[2 * M] + v[4 * M])
            for p in range(M):
                c = t0 + (p + 0.5) * hp
                AG = AG + hp / 18.0 * (5.0 * np.asarray(Vs(c - r * hp), dtype=float) + 8.0 * v[4 * p + 2]
                                       + 5.0 * np.asarray(Vs(c + r * hp), dtype=float))
            d_ref = max(d_ref, float(np.max(np.abs(AM - A2))))
            d_gl = max(d_gl, float(np.max(np.abs(AM - AG))))
            d_sol = max(d_sol, float(np.max(np.abs(S - AM))))
            d_sol2 = max(d_sol2, float(np.max(np.abs(S - A2))))
            scale = max(scale, float(np.max(np.abs(A2))))
        if d_ref <= tol * scale and d_gl <= tol * scale:
            return d_sol / scale
        if 2 * M * N > m_max:
            return d_sol2 / scale
        M *= 2


def php_round(v: float) -> int:
    """Round half away from zero, as PHP's round()."""
    return int(math.floor(v + 0.5)) if v >= 0 else -int(math.floor(-v + 0.5))


# --------------------------------------------------------------------------- solitons
def normalise(c) -> np.ndarray:
    c = np.asarray(c, dtype=complex)
    n = np.linalg.norm(c)
    if c.shape != (2,) or n == 0.0:
        raise ValueError("polarization vector cannot be zero")
    return c / n


def manakov(x: np.ndarray, tau: float, eta: float, xi: float, c, x0: float, scale: float = 1.0) -> np.ndarray:
    """Bright Manakov soliton of i p_tau + p_xx + 2|p|^2 p = 0 (convention of the PHP core):

    p = c 2 eta sech(2 eta (x - x0 + 4 xi tau)) exp(-2 i xi x - 4 i (xi^2 - eta^2) tau),

    i.e. carrier exp(-2 i xi x) and velocity -4 xi.
    """
    c = normalise(c)
    w = 2.0 * eta
    env = scale * w / np.cosh(w * (x - x0 + 4.0 * xi * tau))
    ph = -2.0 * xi * x - 4.0 * (xi * xi - eta * eta) * tau
    return c[:, None] * (env * np.exp(1j * ph))[None, :]


def manakov_periodic(x: np.ndarray, L: float, tau: float, eta: float, xi: float, c, x0: float,
                     scale: float = 1.0) -> np.ndarray:
    """Soliton periodized over the box, sum_m p(x + mL, tau), with the images whose tail reaches the box above
    ~1e-17 of the peak (port of Solitons::manakovPeriodic)."""
    R = math.log(2.0e17) / (2.0 * eta)
    cc = x0 - 4.0 * xi * tau
    lo = int(math.ceil((cc - 0.5 * L - R) / L))
    hi = int(math.floor((cc + 0.5 * L + R) / L))
    Q = manakov(x, tau, eta, xi, c, x0, scale)
    for m in range(lo, hi + 1):
        if m == 0:
            continue
        sel = np.abs(x + m * L - cc) < R          # an image contributes only within R of its centre
        if np.any(sel):
            Q[:, sel] += manakov(x[sel] + m * L, tau, eta, xi, c, x0, scale)
    return Q


def fractional_soliton(x: np.ndarray, k: np.ndarray, alpha: float, mu: float,
                       tol: float = 1e-13, maxit: int = 3000) -> tuple[np.ndarray, float, int]:
    """Real, even solution of (mu + (-d_xx)^{alpha/2}) phi = 2 phi^3 by Petviashvili iteration.

    Returns (phi, residual, iterations); iterations = maxit + 1 means that the iteration did not converge.
    """
    Lk = mu + np.abs(k) ** alpha
    sq = math.sqrt(mu)
    phi = sq / np.cosh(sq * x)
    it = 0
    converged = False
    for it in range(maxit):
        ph = np.fft.fft(phi)
        N = np.fft.fft(2.0 * phi ** 3)
        M = np.sum(Lk * np.abs(ph) ** 2) / np.real(np.sum(np.conj(ph) * N))
        new = np.real(np.fft.ifft(M ** 1.5 * N / Lk))
        diff = np.max(np.abs(new - phi))
        phi = new
        if diff < tol:
            converged = True
            break
    res = np.max(np.abs(np.real(np.fft.ifft(Lk * np.fft.fft(phi))) - 2.0 * phi ** 3))
    return phi, float(res), it + 1 if converged else maxit + 1      # > maxit: not converged (as in PHP)


def relative_residual(phi: np.ndarray, k: np.ndarray, alpha: float, mu: float) -> float:
    """Relative residual max|(mu + |D|^alpha) phi - 2 phi^3| / max|2 phi^3| of a real profile on the grid."""
    Lk = mu + np.abs(k) ** alpha
    c = 2.0 * phi ** 3
    scale = float(np.max(np.abs(c)))
    res = float(np.max(np.abs(np.real(np.fft.ifft(Lk * np.fft.fft(phi))) - c)))
    return res / scale if scale > 0.0 else math.inf


def fourier_shift(f: np.ndarray, k: np.ndarray, x0: float) -> np.ndarray:
    """Periodic shift f(x - x0) by a Fourier phase."""
    if x0 == 0.0:
        return f
    return np.real(np.fft.ifft(np.fft.fft(f) * np.exp(-1j * k * x0)))


# --------------------------------------------------------------------------- split-step
class SplitStepSolver:
    """Strang split-step Fourier solver for

    i q_t - beta(t) (-d_xx)^{alpha/2} q + 2 gamma(t)(|q1|^2+|q2|^2) q + V(x,t) q = i g(t) q,

    with exact linear and nonlinear sub-flows, Simpson-integrated coefficients and fused half-steps.
    """

    def __init__(self, L: float, Nx: int, alpha: float, beta: TimeFn, gamma: TimeFn, g: TimeFn,
                 V: Optional[FieldFn] = None):
        if not math.isfinite(L) or L <= 0.0:
            raise SimulationError("SplitStepSolver: L must be positive.", "L")
        if Nx < 4 or Nx & (Nx - 1):
            raise SimulationError("SplitStepSolver: Nx must be a power of two, at least 4.", "Nx")
        if not math.isfinite(alpha) or alpha <= 0.0 or alpha > 2.0:
            raise SimulationError("SplitStepSolver: the order must satisfy 0 < alpha <= 2.", "alpha")
        self.L, self.n, self.alpha = L, Nx, alpha
        self.beta, self.gamma, self.g, self.V = beta, gamma, g, V
        self.x, self.k, self.dx = grid(L, Nx)
        self.Dk = np.abs(self.k) ** alpha

    def _linear(self, Q: np.ndarray, B: float) -> np.ndarray:
        return np.fft.ifft(np.fft.fft(Q, axis=-1) * np.exp(-1j * self.Dk * B)[None, :], axis=-1)

    def _nonlinear(self, Q: np.ndarray, t0: float, t1: float) -> np.ndarray:
        g, gam = self.g, self.gamma
        tm = 0.5 * (t0 + t1)
        h = t1 - t0
        Gm = simpson(g, t0, tm)
        G1 = Gm + simpson(g, tm, t1)
        cNl = h / 6.0 * 2.0 * (gam(t0) + 4.0 * gam(tm) * math.exp(2.0 * Gm) + gam(t1) * math.exp(2.0 * G1))
        ph = cNl * np.sum(np.abs(Q) ** 2, axis=0)
        if self.V is not None:
            ph = ph + h / 6.0 * (self.V(self.x, t0) + 4.0 * self.V(self.x, tm) + self.V(self.x, t1))
        return Q * (math.exp(G1) * np.exp(1j * ph))[None, :]

    def solve(self, Q: np.ndarray, T: float, dt: float, nSave: int, on_save) -> np.ndarray:
        if not math.isfinite(T) or T <= 0.0 or not math.isfinite(dt) or dt <= 0.0:
            raise SimulationError("SplitStepSolver: T and dt must be positive.", "dt")
        beta = self.beta
        N = max(1, php_round(T / dt))
        dt = T / N
        every = max(1, N // max(1, nSave))
        on_save(0.0, Q)
        Q = self._linear(Q, simpson(beta, 0.0, 0.5 * dt))
        for n in range(1, N + 1):
            t0 = (n - 1) * dt
            t1 = n * dt
            Q = self._nonlinear(Q, t0, t1)
            Bsecond = simpson(beta, t0 + 0.5 * dt, t1)
            if n % every == 0 or n == N:
                Q = self._linear(Q, Bsecond)
                on_save(t1, Q)
                if n < N:
                    Q = self._linear(Q, simpson(beta, t1, t1 + 0.5 * dt))
            else:
                Q = self._linear(Q, Bsecond + simpson(beta, t1, t1 + 0.5 * dt))
        return Q


# --------------------------------------------------------------------------- Caputo L1
def l1_weights(t: np.ndarray, n: int, e: float, G2: float) -> np.ndarray:
    """L1 weights a_{n,k}, k = 1..n, computed without cancellation.

    a_{n,k} = [(t_n - t_{k-1})^e - (t_n - t_k)^e] / (G2 tau_k),  e = 1 - nu,  tau_k = t_k - t_{k-1}.
    With w = t_n - t_k > 0 the difference is w^e expm1(e log1p(tau_k / w)); the direct formula loses
    all digits when tau_k << w, which happens on strongly graded meshes (first steps ~ T N^{-r}).
    """
    tau = t[1:n + 1] - t[:n]
    w = t[n] - t[1:n + 1]
    out = np.empty(n)
    pos = w > 0
    out[pos] = w[pos] ** e * np.expm1(e * np.log1p(tau[pos] / w[pos]))
    out[~pos] = tau[~pos] ** e
    return out / (G2 * tau)


class CaputoSolver:
    """L1 scheme on the graded mesh t_n = T (n/N)^r for

    i D_t^nu q + beta(t) q_xx + 2 gamma(t)(|q1|^2+|q2|^2) q + V q = i g(t) q + f(x,t),   0 < nu <= 1.

    Dispersion and gain implicit (diagonal in Fourier space); nonlinear, potential and forcing terms
    implicit, solved by Anderson-accelerated fixed-point iteration. Full memory is kept.
    With nu = 1 the L1 weights reduce to backward Euler on the graded mesh (first order, dissipative); this
    is accepted for library use and limit checks, while Simulation runs nu = 1 with the split-step solver.
    """

    def __init__(self, L: float, Nx: int, nu: float, beta: TimeFn, gamma: TimeFn, g: TimeFn,
                 V: Optional[FieldFn] = None, forcing: Optional[Callable[[float], np.ndarray]] = None,
                 anderson: int = 10):
        # contract of the L1 scheme (also enforced by Simulation): 0 < nu <= 1; nu = 1 is backward Euler
        if not math.isfinite(nu) or nu <= 0.0 or nu > 1.0:
            raise SimulationError("CaputoSolver: the Caputo order must satisfy 0 < nu <= 1.", "nu")
        if not math.isfinite(L) or L <= 0.0:
            raise SimulationError("CaputoSolver: L must be positive.", "L")
        if Nx < 4 or Nx & (Nx - 1):
            raise SimulationError("CaputoSolver: Nx must be a power of two, at least 4.", "Nx")
        self.L, self.n, self.nu = L, Nx, nu
        self.beta, self.gamma, self.g, self.V, self.forcing = beta, gamma, g, V, forcing
        self.anderson = anderson
        self.x, self.k, self.dx = grid(L, Nx)
        self.k2 = self.k ** 2
        self.max_picard = 0
        self.sum_picard = 0
        self.steps_done = 0
        self.picard_steps: list = []

    @staticmethod
    def mesh(T: float, N: int, r: float) -> np.ndarray:
        return np.array([T * (i / N) ** r for i in range(N + 1)])

    def _rhs(self, Q: np.ndarray, t: float) -> np.ndarray:
        """Nonlinear + potential + forcing terms divided by i:  i(2 gamma |q|^2 + V) q - i f."""
        w = 2.0 * self.gamma(t) * np.sum(np.abs(Q) ** 2, axis=0)
        if self.V is not None:
            w = w + self.V(self.x, t)
        out = 1j * w[None, :] * Q
        if self.forcing is not None:
            out = out - 1j * self.forcing(t)
        return out

    @staticmethod
    def _flat(Q: np.ndarray) -> np.ndarray:
        return np.concatenate([Q[0].real, Q[0].imag, Q[1].real, Q[1].imag])

    @staticmethod
    def _unflat(v: np.ndarray, Nx: int) -> np.ndarray:
        return np.stack([v[:Nx] + 1j * v[Nx:2 * Nx], v[2 * Nx:3 * Nx] + 1j * v[3 * Nx:]])

    @staticmethod
    def _anderson_step(X: list, F: list) -> np.ndarray:
        k = len(X) - 1
        dF = np.array([F[j + 1] - F[j] for j in range(k)])
        dX = np.array([X[j + 1] - X[j] for j in range(k)])
        A = dF @ dF.T
        tr = np.trace(A)
        A = A + 1e-12 * (tr / k if tr > 0.0 else 1.0) * np.eye(k)
        gam = np.linalg.solve(A, dF @ F[k])
        return X[k] + F[k] - (dX + dF).T @ gam

    def solve(self, Q: np.ndarray, T: float, N: int, r: float, nSave: int, on_save,
              tol: float = 1e-12, maxit: int = 100) -> np.ndarray:
        if not math.isfinite(T) or T <= 0.0 or N < 1 or not math.isfinite(r) or r < 1.0:
            raise SimulationError("CaputoSolver: requires T > 0, N >= 1 and r >= 1.",
                                  "N" if N < 1 else ("r" if r < 1.0 else "T"))
        nu = self.nu
        t = self.mesh(T, N, r)
        G2 = math.gamma(2.0 - nu)
        e = 1.0 - nu
        Nx = self.n

        save = set()
        m = min(nSave, N)
        j = 1
        for s in range(m + 1):
            target = T * s / m
            while j < N and t[j] < target:
                j += 1
            save.add(max(1, j))
        save.add(N)

        on_save(0.0, Q)
        Qh = np.fft.fft(Q, axis=-1)
        dU = np.zeros((N + 1, 2, Nx), dtype=complex)
        self.max_picard = self.sum_picard = self.steps_done = 0
        self.picard_steps = []

        for n in range(1, N + 1):
            tn = t[n]
            a = np.concatenate([[0.0], l1_weights(t, n, e, G2)])
            ann = a[n]
            H = np.tensordot(a[1:n], dU[1:n], axes=(0, 0)) if n > 1 else np.zeros_like(Qh)
            known = ann * Qh - H
            den = (ann - self.g(tn)) + 1j * self.beta(tn) * self.k2

            Qn = Q
            QhNew = None
            converged = False
            diff = math.inf
            hX: list = []
            hF: list = []
            it = 0
            for it in range(1, maxit + 1):
                GQh = (known + np.fft.fft(self._rhs(Qn, tn), axis=-1)) / den[None, :]
                GQ = np.fft.ifft(GQh, axis=-1)
                fv = self._flat(GQ) - self._flat(Qn)
                diff = float(np.max(np.abs(fv)))
                scale = max(1.0, float(np.max(np.abs(self._flat(GQ)))))
                if not math.isfinite(diff):
                    break
                if diff < tol * scale:
                    Qn, QhNew, converged = GQ, GQh, True
                    break
                if self.anderson < 1:
                    Qn = GQ
                    continue
                hX.append(self._flat(Qn))
                hF.append(fv)
                if len(hX) > self.anderson + 1:
                    hX.pop(0)
                    hF.pop(0)
                Qn = GQ if len(hX) == 1 else self._unflat(self._anderson_step(hX, hF), Nx)
            if not converged:
                raise SimulationError(
                    f"The implicit nonlinear solve did not converge at step {n} (t = {tn:.4g}, last change {diff:.2e}). "
                    "The nonlinearity is too strong for this step: increase N, reduce the amplitude, or increase nu.", "N")
            self.max_picard = max(self.max_picard, it)
            self.sum_picard += it
            self.picard_steps.append(it)
            dU[n] = QhNew - Qh
            self.steps_done = n
            Qh = QhNew
            Q = Qn
            if n in save:
                on_save(tn, Q)
        return Q
