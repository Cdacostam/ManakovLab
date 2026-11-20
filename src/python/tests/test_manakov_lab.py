"""Tests of the Python implementation of Manakov Lab (mirror of src/php/tests/run_tests.php)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from manakov_lab import SimulationError, Simulation, fractional_soliton, grid, problem, run_preset
from manakov_lab.core import l1_weights


def test_soliton_convention():
    """Carrier exp(-2 i xi x): the peak travels with velocity -4 xi (beta = gamma = 1, g = 0)."""
    xi, T = 0.25, 4.0
    r = run_preset('classic-exact', params=dict(T=T, dt=0.01),
                   coeffs=dict(beta=lambda t: 1.0, g=lambda t: 0.0),
                   solitons=[dict(eta=0.6, xi=xi, x0=0, p1re=1, p1im=0, p2re=0, p2im=0)],
                   exact=False, convergence=False)
    I = np.abs(r['final'][0]) ** 2 + np.abs(r['final'][1]) ** 2
    assert abs(r['x'][np.argmax(I)] - (-4 * xi * T)) < 0.2


def test_petviashvili_residual():
    x, k, _ = grid(80.0, 512)
    alpha, mu = 1.5, 1.0
    phi, res, _ = fractional_soliton(x, k, alpha, mu)
    Lphi = np.real(np.fft.ifft((mu + np.abs(k) ** alpha) * np.fft.fft(phi)))
    assert res < 1e-12 and np.max(np.abs(Lphi - 2 * phi ** 3)) < 1e-12


def test_power_conservation_general_coefficients():
    r = run_preset('classic-collision', params=dict(T=5.0))
    d = r['diagnostics']
    P = (d['P1'] + d['P2']) / d['gain']
    assert np.max(np.abs(P / P[0] - 1)) < 1e-12


@pytest.mark.parametrize('preset', ['classic-exact', 'space-exact'])
def test_split_step_order_two(preset):
    c = run_preset(preset)['convergence']
    assert all(abs(o - 2) < 0.02 for o in c['order'][1:])


def test_caputo_orders():
    graded = run_preset('caputo-mms', params=dict(N=256))['convergence']
    uniform = run_preset('caputo-mms', params=dict(N=256, r=1))['convergence']
    assert abs(graded['order'][-1] - 1.4) < 0.02
    assert abs(uniform['order'][-1] - 0.6) < 0.02


def test_l1_weights_stable():
    """Weights a_{n,k} (stored at index k - 1) agree with a 60-digit evaluation on a mesh with r = 10."""
    T, N, r, nu = 1.0, 64, 10.0, 0.1
    t = T * (np.arange(N + 1) / N) ** r
    n = N
    a = l1_weights(t, n, 1 - nu, math.gamma(2 - nu))
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    e, G2 = Decimal(1 - nu), Decimal(math.gamma(2 - nu))
    tn = Decimal(t[n])
    for kk in (1, 2, 5, N):
        a0, a1 = Decimal(t[kk - 1]), Decimal(t[kk])
        tau = a1 - a0
        ref = ((tn - a0) ** e - ((tn - a1) ** e if tn > a1 else Decimal(0))) / (G2 * tau)
        assert abs(a[kk - 1] - float(ref)) <= 1e-12 * abs(float(ref))


@pytest.mark.parametrize('nu', [0.1, 0.05])
def test_caputo_small_nu_converges(nu):
    r = run_preset('caputo-decay', params=dict(nu=nu, T=2.0, N=128, L=80.0, Nx=128))
    assert r['meta']['picardMax'] < 100
    assert np.all(np.isfinite(r['final']))


def test_caputo_amplitude_one_mms():
    sol = [dict(eta=0.5, xi=-0.25, x0=0, p1re=1, p1im=0, p2re=0.5, p2im=-0.5)]
    r = run_preset('caputo-mms', params=dict(N=256), solitons=sol, convergence=False)
    assert r['exact']['err_max'] < 5e-2


def test_caputo_nu_one_uses_split_step():
    r = run_preset('caputo-decay', params=dict(nu=1.0, T=2.0, N=200))
    assert r['meta']['picardMax'] is None
    assert any('nu = 1' in w for w in r['warnings'])
    d = r['diagnostics']
    assert abs((d['P1'][-1] + d['P2'][-1]) / (d['P1'][0] + d['P2'][0]) - 1) < 1e-12
    with pytest.raises(SimulationError) as e:
        run_preset('caputo-mms', params=dict(nu=1.0))
    assert e.value.field == 'nu'


def test_caputo_close_to_one_warning():
    r = run_preset('caputo-decay', params=dict(nu=0.97, T=1.0, N=50, L=80.0, Nx=128))
    assert any('close to 1' in w for w in r['warnings'])


def test_alpha_le_one():
    r = run_preset('space-collision', params=dict(alpha=1.0, T=1.0))
    assert any('alpha <= 1' in w for w in r['warnings'])
    with pytest.raises(SimulationError) as e:
        run_preset('space-exact', params=dict(alpha=0.9))
    assert e.value.field == 'alpha'


def test_non_finite_coefficient():
    with pytest.raises(SimulationError) as e:
        Simulation(problem('classic-exact', coeffs=dict(beta=lambda t: math.inf)))
    assert e.value.field == 'beta'


def test_unknown_override():
    with pytest.raises(TypeError):
        problem('classic-exact', foo=1)


def test_crosscheck_against_php():
    """Agreement with the PHP engine on the cases of src/php/tests/crosscheck_cases.json (if exported)."""
    import importlib.util
    from pathlib import Path
    root = Path(__file__).resolve().parents[3]
    ref = root / 'results' / 'php_reference.json'
    if not ref.exists():
        pytest.skip('results/php_reference.json not found (run src/php/tests/export_reference.php)')
    spec = importlib.util.spec_from_file_location('crosscheck', root / 'src' / 'python' / 'crosscheck.py')
    cc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cc)
    assert cc.main(str(ref)) == 0


# ------------------------------------------------------------------ validity of the references and domain
def _boundary_soliton(L, Nx, T):
    return problem('classic-exact', params=dict(L=L, Nx=Nx, T=T, dt=0.01),
                   coeffs=dict(beta=lambda t: 1.0, gamma=lambda t: 1.0, g=lambda t: 0.0, V=None),
                   solitons=[dict(eta=0.5, xi=-0.5, x0=0, p1re=1, p1im=0, p2re=0, p2im=0)], convergence=True)


@pytest.mark.parametrize('T', [7.0, 14.0])
def test_exact_crossing_boundary_periodized(T):
    """Centre at x = 14 and x = 28 > L/2 = 20: the periodized reference keeps order 2."""
    r = Simulation(_boundary_soliton(40.0, 512, T)).run()
    assert abs(r['convergence']['order'][-1] - 2.0) < 0.02
    assert r['exact']['refFloorParts']['images'] < 1e-14


def test_exact_refused_when_overlapping_images():
    with pytest.raises(SimulationError) as e:
        Simulation(_boundary_soliton(12.0, 128, 1.0)).run()       # 2 sech(eta L) = 9.9e-3 > 1e-3
    assert e.value.field == 'L'


def test_exact_spatial_floor():
    r = run_preset('classic-exact', params=dict(T=1.0, Nx=128, dt=0.0025))
    C = r['convergence']
    assert any('does not resolve the soliton spectrum' in w for w in r['warnings']) and C['refFloor'] > 1e-4
    assert all(o is None or e > 10 * C['refFloor'] for o, e in zip(C['order'], C['err']))
    assert C['order'][-1] is None
    r = run_preset('classic-exact', params=dict(T=1.0, Nx=512, dt=0.0025))
    assert abs(r['convergence']['order'][-1] - 2.0) < 0.01


def test_adaptive_reference_quadrature():
    from manakov_lab.core import antiderivative, antiderivative_adaptive
    f = lambda t: 300.0 * math.sin(40.0 * math.pi * t)
    F = lambda t: 300.0 * (1.0 - math.cos(40.0 * math.pi * t)) / (40.0 * math.pi)
    fixed = antiderivative(f, 100.0)
    adapt, info = antiderivative_adaptive(f, 100.0)
    ts = [99.9 * i / 4000 + 0.0123 for i in range(4001)]
    e0 = max(abs(fixed(t) - F(t)) for t in ts)
    e1 = max(abs(adapt(t) - F(t)) for t in ts)
    assert info['resolved'] and info['m'] > 4096
    assert e0 > 1e-2 and e1 < 1e-7
    _, info = antiderivative_adaptive(lambda t: 1.0 + 0.5 * math.sin(t), 5.0)
    assert info['resolved'] and info['m'] == 4096


def test_fractional_profile_branches():
    sol = lambda eta, x0: [dict(eta=eta, xi=0, x0=x0, p1re=1, p1im=0, p2re=0.6, p2im=0.4)]
    r = run_preset('space-exact', params=dict(L=5.0, Nx=128, T=0.5, dt=0.0025), solitons=sol(0.25, 0.0))
    assert r['meta']['profile']['branch'] == 'homogeneous' and any('spatially constant' in w for w in r['warnings'])
    assert all(o is None for o in r['convergence']['order'])
    r = run_preset('space-exact', params=dict(L=10.0, Nx=128, T=1.0, dt=0.005), solitons=sol(0.2495, 0.0390625))
    assert r['meta']['profile']['iterations'] is None and any('did not converge' in w for w in r['warnings'])
    assert r['meta']['profile']['branch'] == 'periodic'
    assert r['exact']['refFloorParts']['profile'] > 1e-7
    assert all(o is None for o in r['convergence']['order'])


def test_step_quadrature_and_effective_dt():
    c = dict(beta=lambda t: 0.0, gamma=lambda t: 0.0, g=lambda t: math.cos(8 * math.pi * t), V=None)
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.5), coeffs=c)
    assert any('does not resolve g(t)' in w for w in r['warnings']) and r['meta']['quadrature']['g'] > 1e-2
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.25), coeffs=c)
    assert not any('does not resolve' in w for w in r['warnings'])
    r = run_preset('classic-exact', params=dict(T=1.0, dt=0.0119))
    assert all(abs(1 / h - round(1 / h)) < 1e-9 for h in r['convergence']['h'])
    assert abs(r['meta']['dtEffective'] - 1 / 84) < 1e-15


def test_caputo_nu_one_step_is_T_over_N():
    a = run_preset('caputo-decay', params=dict(nu=1.0, T=2.0, N=200, dt=0.5))
    b = run_preset('caputo-decay', params=dict(nu=1.0, T=2.0, N=200))
    assert np.array_equal(a['final'], b['final']) and any('was ignored' in w for w in a['warnings'])
    assert a['meta']['dtEffective'] == 2.0 / 200


@pytest.mark.parametrize('make,field', [
    (lambda c: c.CaputoSolver(2 * math.pi, 64, 1.1, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0), 'nu'),
    (lambda c: c.CaputoSolver(2 * math.pi, 64, 0.0, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0), 'nu'),
    (lambda c: c.CaputoSolver(2 * math.pi, 60, 0.5, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0), 'Nx'),
    (lambda c: c.SplitStepSolver(2 * math.pi, 64, 2.5, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0), 'alpha'),
    (lambda c: c.SplitStepSolver(-1.0, 64, 2.0, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0), 'L'),
])
def test_solvers_validate_contract(make, field):
    from manakov_lab import core
    with pytest.raises(SimulationError) as e:
        make(core)
    assert e.value.field == field


def test_fractional_off_node_centre_snapped():
    sol = [dict(eta=0.5, xi=0, x0=0.15625, p1re=1, p1im=0, p2re=0.6, p2im=0.4)]
    r = run_preset('space-exact', params=dict(Nx=256, T=1.0), solitons=sol)
    assert any('nearest node' in w for w in r['warnings'])
    assert r['meta']['profileResidual'] < 1e-11
    assert abs(r['convergence']['order'][-1] - 2.0) < 0.02


@pytest.mark.parametrize('preset,params,field', [
    ('caputo-decay', dict(nu=1.1, N=50, T=1.0, Nx=128), 'nu'),
    ('caputo-decay', dict(nu=0.0, N=50, T=1.0, Nx=128), 'nu'),
    ('space-collision', dict(alpha=2.5, T=1.0), 'alpha'),
    ('classic-collision', dict(dt=-0.01, T=1.0), 'dt'),
    ('classic-collision', dict(Nx=500, T=1.0), 'Nx'),
])
def test_domain_rejected(preset, params, field):
    with pytest.raises(SimulationError) as e:
        run_preset(preset, params=params)
    assert e.value.field == field


def test_caputo_mittag_leffler():
    """i D^nu q + q_xx = 0, q = c e^{ix}: q(t) = E_{1/2}(-i sqrt t) q(0) = e^{-t}(1 - i erfi(sqrt t)) q(0)."""
    from manakov_lab.core import CaputoSolver

    def erfi(x):
        return 2 / math.sqrt(math.pi) * sum(x ** (2 * n + 1) / (math.factorial(n) * (2 * n + 1)) for n in range(80))
    T = 0.25
    ref = math.exp(-T) * (1 - 1j * erfi(math.sqrt(T)))
    errs = []
    for N in (32, 64, 128, 256):
        s = CaputoSolver(2 * math.pi, 64, 0.5, lambda t: 1.0, lambda t: 0.0, lambda t: 0.0, None, None)
        Q0 = np.vstack([np.exp(1j * s.x), 0.5 * np.exp(1j * s.x)])
        out = {}
        s.solve(Q0, T, N, 3.0, 400, lambda t, Q: out.__setitem__('Q', Q.copy()))
        errs.append(float(np.max(np.abs(out['Q'] - ref * Q0))))
    assert abs(errs[-1] - 8.058e-5) < 1e-7
    assert abs(math.log2(errs[-2] / errs[-1]) - 1.5) < 0.05


def test_saturated_orders_withheld():
    """D2-H2: orders whose errors stop decreasing with dt are withheld (a posteriori test)."""
    r = run_preset('classic-exact', params=dict(Nx=128))
    C = r['convergence']
    assert len(C['saturated']) >= 2 and all(o is None for o in C['order'])
    assert any('do not decrease as a power' in w for w in r['warnings'])
    r = run_preset('classic-exact', params=dict(T=4.0, Nx=128, dt=0.0025))
    assert all(o is None for o in r['convergence']['order'])
    r = run_preset('classic-exact', params=dict(Nx=256))
    C = r['convergence']
    assert C['saturated'] == [] and C['preasymptotic'] == [] and abs(C['order'][-1] - 2.0) < 0.01


def test_error_sampled_at_every_step_when_outputs_alias():
    """N1: beta = gamma = 300 sin(40 pi t) with outputs every 10 steps (one period) is stroboscopic."""
    osc = lambda t: 300.0 * math.sin(40.0 * math.pi * t)
    c = dict(beta=osc, gamma=osc, g=lambda t: 0.0, V=None)
    sol = [dict(eta=0.5, xi=2.5, x0=0, p1re=1, p1im=0, p2re=0, p2im=0)]
    r = run_preset('classic-exact', params=dict(L=80.0, Nx=256, T=10.0, dt=0.0025), coeffs=c, solitons=sol,
                   convergence=False)
    assert r['exact']['errorSampling'] == 'every step' and float(r['exact']['err_max']) > 0.1
    assert any('under-sample' in w for w in r['warnings'])
    assert run_preset('classic-exact')['exact']['errorSampling'] == 'outputs'


def test_nonlinear_phase_quadrature_matches_solver():
    c = dict(beta=lambda t: 0.0, gamma=lambda t: math.cos(4 * math.pi * t), g=lambda t: 0.0, V=None)
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.5), coeffs=c)
    assert abs(r['meta']['quadrature']['gamma'] - 1 / 3) < 1e-6
    assert any('gain-weighted Simpson' in w for w in r['warnings'])
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.25), coeffs=c)
    assert r['meta']['quadrature']['gamma'] < 1e-12 and not any('does not resolve' in w for w in r['warnings'])


def test_reference_table_gauss_legendre_check():
    from manakov_lab.core import antiderivative_adaptive
    w = 32768 * math.pi
    F, info = antiderivative_adaptive(lambda t: math.cos(w * t), 1.0)
    assert info['resolved'] and info['m'] > 4096
    assert abs(F(1.0)) < 1e-6 and abs(F(0.3) - math.sin(w * 0.3) / w) < 1e-6
    c = dict(beta=lambda t: 0.0, gamma=lambda t: 0.0, g=lambda t: math.cos(w * t), V=None)
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.5), coeffs=c)
    assert r['meta']['quadrature']['g'] > 0.5 and any('does not resolve g(t)' in w for w in r['warnings'])


def test_saturation_fit_unequal_ratios_and_two_levels():
    from manakov_lab.simulation import fit_power_law
    h = [1 / 10, 1 / 20, 1 / 41, 1 / 81]
    p, S = fit_power_law(h[1], h[2], h[3], h[1] ** 2, h[2] ** 2, h[3] ** 2)
    assert abs(p - 2.0) < 1e-9 and abs(S) < 1e-15
    C = run_preset('caputo-mms', params=dict(N=36))['convergence']          # levels N = 16, 18, 36
    assert C['saturated'] == [] and abs(C['order'][2] - 1.3336) < 1e-3
    C = run_preset('classic-exact', params=dict(T=0.25, dt=0.25 / 81))['convergence']
    assert C['saturated'] == [] and abs(C['order'][-1] - 2.0) < 0.01
    c0 = dict(beta=lambda t: 1.0, gamma=lambda t: 1.0, g=lambda t: 0.0, V=None)
    sol = lambda eta, xi: [dict(eta=eta, xi=xi, x0=0, p1re=1, p1im=0, p2re=0.6, p2im=0.4)]
    C = run_preset('classic-exact', params=dict(L=80.0, Nx=64, T=4.0, dt=0.16), coeffs=c0,
                   solitons=sol(0.35, 0))['convergence']
    assert len(C['h']) == 2 and C['order'] == [None, None] and C['saturated'] == [0.16]
    r = run_preset('classic-exact', params=dict(L=80.0, Nx=128, T=4.0, dt=0.16), coeffs=c0, solitons=sol(0.5, -1))
    assert r['convergence']['order'] == [None, None] and r['convergence']['untested'] == [0.16]
    assert any('only two levels' in w for w in r['warnings'])


def test_potential_quadrature():
    c = lambda V: dict(beta=lambda t: 0.0, gamma=lambda t: 0.0, g=lambda t: 0.0, V=V)
    V8 = lambda x, t: math.cos(8 * math.pi * t) + 0.0 * x
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.5), coeffs=c(V8))
    assert abs(r['meta']['quadrature']['V'] - 1.0) < 1e-9 and any('does not resolve V(x,t)' in w for w in r['warnings'])
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=1 / 64), coeffs=c(V8))
    assert r['meta']['quadrature']['V'] < 1e-6 and not any('does not resolve' in w for w in r['warnings'])
    r = run_preset('classic-collision', params=dict(T=1.0, Nx=128, dt=0.5), coeffs=c(lambda x, t: 0.01 * x ** 2))
    assert 'V' not in r['meta']['quadrature']


def test_expected_order_metadata():
    r = run_preset('caputo-decay', params=dict(nu=0.05, r=10.0, N=32, T=0.5, Nx=128))
    assert abs(r['meta']['expectedOrder'] - 0.5) < 1e-12
    r = run_preset('caputo-decay', params=dict(nu=1.0, N=20, T=1.0, Nx=128))
    assert r['meta']['expectedOrder'] == 2.0 and abs(r['meta']['dtEffective'] - 0.05) < 1e-15
    assert r['meta']['runtime_total'] >= r['meta']['runtime_main']
