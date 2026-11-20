"""The eight presets of the web interface (``src/php/src/Presets.php``), with coefficients as Python callables."""
from __future__ import annotations

import copy
import math

import numpy as np


def _sol(eta, xi, x0, p1re, p1im, p2re, p2im):
    return dict(eta=eta, xi=xi, x0=x0, p1re=p1re, p1im=p1im, p2re=p2re, p2im=p2im)


ONE = lambda t: 1.0
ZERO = lambda t: 0.0

PRESETS = {
    'classic-exact': dict(
        title='Classical: nonautonomous soliton, exact + convergence', model='classic',
        params=dict(L=80.0, Nx=512, T=10.0, dt=0.005),
        coeffs=dict(beta=lambda t: 1 + 0.5 * math.cos(2 * t), gamma=ONE, g=lambda t: 0.05 * math.sin(t), V=None),
        solitons=[_sol(0.6, 0.15, 10, 1, 0, 0.6, 0.4)], exact=True, convergence=True),
    'classic-collision': dict(
        title='Classical: shape-changing collision', model='classic',
        params=dict(L=100.0, Nx=512, T=30.0, dt=0.01),
        coeffs=dict(beta=lambda t: 1 + 0.5 * math.cos(2 * t),
                    gamma=lambda t: (1 + 0.5 * math.cos(2 * t)) * math.exp(-0.1 * (1 - math.cos(t))),
                    g=lambda t: 0.05 * math.sin(t), V=None),
        solitons=[_sol(0.5, -0.25, -15, 1, 0, 1, 0), _sol(0.6, 0.25, 15, 1, 0, 0, 0)], exact=False, convergence=False),
    'classic-trap': dict(
        title='Classical: dispersion management in a harmonic trap', model='classic',
        params=dict(L=120.0, Nx=512, T=40.0, dt=0.01),
        coeffs=dict(beta=lambda t: 1 + 0.8 * math.sin(4 * t), gamma=ONE, g=ZERO, V=lambda x, t: -0.01 * x ** 2),
        solitons=[_sol(0.5, 0, -12, 1, 0, 0, 0.5)], exact=False, convergence=False),
    'space-exact': dict(
        title='Space-fractional: soliton α = 1.5, exact + convergence', model='space',
        params=dict(L=80.0, Nx=512, T=10.0, dt=0.005, alpha=1.5),
        coeffs=dict(beta=lambda t: 1 + 0.5 * math.cos(2 * t), gamma=ONE, g=lambda t: 0.05 * math.sin(t), V=None),
        solitons=[_sol(0.5, 0, 0, 1, 0, 0.6, 0.4)], exact=True, convergence=True),
    'space-collision': dict(
        title='Space-fractional: inelastic collision, α = 1.5', model='space',
        params=dict(L=160.0, Nx=1024, T=25.0, dt=0.01, alpha=1.5),
        coeffs=dict(beta=ONE, gamma=ONE, g=ZERO, V=None),
        solitons=[_sol(0.5, -0.5, -15, 1, 0, 1, 0), _sol(0.5, 0.5, 15, 1, 0, 0, 0)], exact=False, convergence=False),
    'caputo-mms': dict(
        title='Caputo: manufactured solution ν = 0.6, exact + convergence', model='caputo',
        params=dict(L=60.0, Nx=256, T=1.0, N=512, nu=0.6, r='auto'),
        coeffs=dict(beta=lambda t: 1 + 0.5 * math.cos(2 * t), gamma=lambda t: 1 + 0.3 * math.sin(t),
                    g=lambda t: 0.1 * math.cos(t), V=lambda x, t: 0.05 * np.cos(0.5 * x) * (1 + t)),
        solitons=[_sol(0.25, -0.25, 0, 1, 0, 0.5, -0.5)], exact=True, convergence=True),
    'caputo-decay': dict(
        title='Caputo: soliton decay under memory, ν = 0.8', model='caputo',
        params=dict(L=160.0, Nx=512, T=15.0, N=600, nu=0.8, r='auto'),
        coeffs=dict(beta=lambda t: 1 + 0.5 * math.cos(2 * t), gamma=ONE, g=ZERO, V=None),
        solitons=[_sol(0.6, 0.1, 5, 1, 0, 0, 0.7)], exact=False, convergence=False),
    'caputo-collision': dict(
        title='Caputo: two-soliton interaction, ν = 0.9', model='caputo',
        params=dict(L=100.0, Nx=256, T=20.0, N=700, nu=0.9, r='auto'),
        coeffs=dict(beta=ONE, gamma=ONE, g=ZERO, V=None),
        solitons=[_sol(0.6, -0.4, -12, 1, 0, 1, 0), _sol(0.6, 0.4, 12, 1, 0, 0, 0)], exact=False, convergence=False),
}


def problem(preset_id: str, params: dict | None = None, coeffs: dict | None = None, **over) -> dict:
    """Problem dict for :class:`~manakov_lab.simulation.Simulation` from a preset and overrides.

    ``params`` and ``coeffs`` are merged into the preset's; ``model``, ``solitons``, ``exact`` and
    ``convergence`` replace it. ``r='auto'`` becomes min(10, (2 - nu)/nu), as in the web API.
    """
    p = copy.deepcopy(PRESETS[preset_id])
    p['params'].update(params or {})
    p['coeffs'].update(coeffs or {})
    for key in ('model', 'solitons', 'exact', 'convergence', 'anderson'):
        if key in over:
            p[key] = over.pop(key)
    if over:
        raise TypeError(f'unknown overrides: {sorted(over)}')
    q = p['params']
    P = dict(model=p['model'], L=float(q['L']), Nx=int(q['Nx']), T=float(q['T']), dt=q.get('dt'), N=q.get('N'),
             alpha=float(q.get('alpha', 2.0)), nu=float(q.get('nu', 1.0)), solitons=p['solitons'],
             exact=p['exact'], convergence=p['convergence'], anderson=int(p.get('anderson', 10)), **p['coeffs'])
    if P['model'] == 'caputo':
        r = q.get('r', 'auto')
        if r != 'auto':
            P['r'] = float(r)
        else:   # out-of-range nu is reported by Simulation, not by a division error here
            P['r'] = min(10.0, (2 - P['nu']) / P['nu']) if 0.0 < P['nu'] < 1.0 else 1.0
    return P
