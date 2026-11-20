"""Figures 1-4 of the article describing Manakov Lab (submitted to SoftwareX).

    python reproduce/make_figures.py [results/reproduce.json] [--out DIR]

Fig. 2 uses the convergence data written by ``src/php/tests/reproduce.php``; Figs. 3 and 4 are computed
with the Python port of the numerical core, which agrees with the PHP engine to about 1e-12
(see ``src/python/crosscheck.py``). Output, by default in ``reproduce/figures/``: ``Figure_1.pdf`` (architecture of
the web application), ``Figure_2.pdf`` (convergence), ``Figure_3.pdf`` (collision) and ``Figure_4.pdf`` (power under
Caputo memory), with PNG copies alongside.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python" / "src"))
from manakov_lab import Simulation, problem  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "savefig.bbox": "tight", "font.family": "serif", "mathtext.fontset": "cm"})
OUT = ROOT / "reproduce" / "figures"


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)


def power_history(preset, **over):
    r = Simulation(problem(preset, **over)).run()
    d = r["diagnostics"]
    P = d["P1"] + d["P2"]
    return r["t"], P / P[0]


def fig_validation(rep: dict):
    t2 = rep["table3"]
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.1))

    a = ax[0]
    for key, lab, mk in (("classic_exact", "classical, exact", "o"), ("space_exact", r"fractional $\alpha=1.5$, exact", "s")):
        a.loglog(t2[key]["h"], t2[key]["err"], mk + "-", label=lab, ms=4)
    if "alpha" in rep:
        r12 = rep["alpha"][0]
        a.loglog(r12["h"], r12["err"], "^-", label=r"fractional $\alpha=1.2$, exact", ms=4)
    h = np.array(t2["classic_exact"]["h"])
    a.loglog(h, 0.6 * t2["classic_exact"]["err"][0] * (h / h[0]) ** 2, "k--", lw=0.8, label="slope 2")
    a.set_xlabel(r"$\Delta t$")
    a.set_xticks(t2["classic_exact"]["h"])
    a.set_xticklabels([f"{v:g}" for v in t2["classic_exact"]["h"]])
    a.minorticks_off()
    a.set_ylabel(r"max error")
    a.set_title("(a) split-step", loc="left")
    a.legend(frameon=False)

    a = ax[1]
    g, u = t2["caputo_graded"], t2["caputo_uniform"]
    N = np.array(g["h"], float)
    a.loglog(N, g["err"], "o-", ms=4, label=r"graded $r=7/3$")
    a.loglog(N, u["err"], "s-", ms=4, label=r"uniform $r=1$")
    a.loglog(N, 0.6 * g["err"][0] * (N / N[0]) ** -1.4, "k--", lw=0.8, label=r"slope $-(2-\nu)$")
    a.loglog(N, 1.4 * u["err"][0] * (N / N[0]) ** -0.6, "k:", lw=0.8, label=r"slope $-\nu$")
    a.set_xlabel(r"$N$")
    a.set_xticks(N)
    a.set_xticklabels([f"{int(v)}" for v in N])
    a.minorticks_off()
    a.set_title(r"(b) Caputo L1, $\nu=0.6$", loc="left")
    a.legend(frameon=False)

    fig.tight_layout()
    save(fig, "Figure_2")


def fig_power():
    fig, a = plt.subplots(figsize=(4.6, 3.1))
    t, P = power_history("caputo-decay", model="classic", params=dict(dt=0.01))
    a.plot(t, P, "k-", lw=1.2)
    a.annotate(r"$\nu=1$", (t[-1], P[-1]), xytext=(-4, -9), textcoords="offset points", ha="right")
    colors = {0.9: "C0", 0.8: "C1", 0.7: "C2", 0.6: "C3"}
    ends = {}
    for nu, c in colors.items():
        t, P = power_history("caputo-decay", params=dict(nu=nu, N=600))
        a.plot(t, P, color=c, lw=1.2)
        ends[nu] = P[-1]
    # label positions separated by at least 0.06 (nu = 0.8 and 0.9 end 0.012 apart)
    ypos, last = {}, -1.0
    for nu in sorted(ends, key=ends.get):
        ypos[nu] = max(ends[nu], last + 0.06)
        last = ypos[nu]
    for nu, c in colors.items():
        a.annotate(rf"$\nu={nu}$", (t[-1], ends[nu]), xytext=(t[-1] + 0.4, ypos[nu]), textcoords="data",
                   va="center", color=c, fontsize=8)
    a.set_xlim(0, 18.5)
    a.set_ylim(0, 1.05)
    a.set_xlabel(r"$t$")
    a.set_ylabel(r"$P(t)/P(0)$")
    fig.tight_layout()
    save(fig, "Figure_4")


def fig_architecture():
    """Fig. 1: architecture of the web application (boxes and arrows)."""
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    def layer(y0, y1, label):
        ax.add_patch(FancyBboxPatch((1, y0), 98, y1 - y0, boxstyle="round,pad=0,rounding_size=1.5",
                                    fc="0.97", ec="0.7", lw=0.8))
        ax.text(2.2, y1 - 1.2, label, ha="left", va="top", fontsize=8, fontweight="bold", color="0.4",
                family="sans-serif")

    def box(xc, yc, w, title, sub, key=False):
        ax.add_patch(FancyBboxPatch((xc - w / 2, yc - 3.4), w, 6.8, boxstyle="round,pad=0,rounding_size=1",
                                    fc="#eef2fb" if key else "white", ec="#1f3f99" if key else "0.45",
                                    lw=1.2 if key else 0.8))
        ax.text(xc, yc + 1.3, title, ha="center", va="center", fontsize=9, fontweight="bold", family="sans-serif")
        ax.text(xc, yc - 1.5, sub, ha="center", va="center", fontsize=7.5, family="sans-serif")
        return (xc, yc, w)

    def arrow(p, q, both=False, label=None, side="right", ly=None):
        ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="<->" if both else "->", color="0.35", lw=0.9))
        if label:
            my = ly if ly is not None else (p[1] + q[1]) / 2
            mx = p[0] + (q[0] - p[0]) * (my - p[1]) / (q[1] - p[1])
            ax.text(mx + (1 if side == "right" else -1), my, label, ha="left" if side == "right" else "right",
                    va="center", fontsize=7.5, color="0.35", family="sans-serif")

    layer(48, 61, "Browser")
    idx = box(14, 53, 22, "lab.php", "page, forms, theory")
    app = box(38.5, 53, 21, "app.js", "controller, export", key=True)
    ktx = box(62, 53, 20, "KaTeX", "renders all mathematics")
    plt_ = box(85, 53, 22, "Plotly.js", "3D surfaces, 2D charts")
    layer(28, 44, "API layer (PHP)")
    ltx = box(24, 34, 30, "api/latex.php", "formula → AST → LaTeX")
    run = box(66, 34, 34, "api/run.php", "validate, check budget, run")
    layer(1, 24, "Numerical core (src/Numerics, pure PHP)")
    ssf = box(17, 15, 27, "SplitStepSolver", "Strang, |k|$^\\alpha$, gain")
    sim = box(50, 15, 27, "Simulation", "exact solutions, diagnostics")
    cap = box(83, 15, 27, "CaputoSolver", "L1 graded mesh, Anderson")
    cmp_ = box(17, 5.5, 27, "ExprCompiler", "AST → closures")
    fft = box(50, 5.5, 27, "FFT", "radix-2, cached plans")
    sol = box(83, 5.5, 27, "Solitons", "Manakov, Petviashvili")
    arrow((33, 49.6), (24, 37.4), both=True, label="formula ↔ LaTeX", side="left", ly=46)
    arrow((44, 49.6), (62, 37.4), both=True, label="problem JSON ↔ results", ly=46)
    arrow((60, 30.6), (52, 18.4), label="in-process call", ly=26)
    arrow((36.5, 15), (30.5, 15))
    arrow((63.5, 15), (69.5, 15))
    arrow((50, 11.6), (50, 8.9))
    arrow((40, 11.6), (24, 8.9))
    arrow((60, 11.6), (76, 8.9))
    save(fig, "Figure_1")


def fig_collision():
    P = problem("classic-collision")
    sim = Simulation(P)
    Q0 = sim.initial_condition()
    solver = sim._make_solver()
    ts, I1, I2 = [], [], []

    def cb(t, Q):
        ts.append(t)
        I1.append(np.abs(Q[0]) ** 2)
        I2.append(np.abs(Q[1]) ** 2)

    sim._run_solver(solver, Q0, 300, cb)
    x = sim.x
    sel = np.abs(x) <= 30
    T, X = np.array(ts), x[sel]
    I1, I2 = np.array(I1)[:, sel], np.array(I2)[:, sel]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.0), sharey=True)
    for a, data, title in zip(ax, (I1 + I2, I1, I2),
                              (r"(a) $|q_1|^2+|q_2|^2$", r"(b) $|q_1|^2$", r"(c) $|q_2|^2$")):
        m = a.pcolormesh(X, T, data, shading="auto", cmap="magma", rasterized=True,
                         norm=matplotlib.colors.PowerNorm(0.5))
        a.set_xlabel(r"$x$")
        a.set_title(title, loc="left")
        a.set_xlim(-30, 30)
        fig.colorbar(m, ax=a, pad=0.02)
    ax[0].set_ylabel(r"$t$")
    fig.tight_layout()
    save(fig, "Figure_3")


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--out" in args:
        i = args.index("--out")
        OUT = Path(args[i + 1]).resolve()
        del args[i:i + 2]
    OUT.mkdir(parents=True, exist_ok=True)
    path = Path(args[0]) if args else ROOT / "results" / "reproduce.json"
    rep = json.loads(path.read_text())
    fig_architecture()
    fig_validation(rep)
    fig_collision()
    fig_power()
    print("written:", ", ".join(str(p) for p in sorted(OUT.glob("Figure_*.pdf"))))
