"""Time refinement of the caputo-decay preset beyond the web work budget (Table 4 of the article).

    python caputo_refinement.py > ../../results/caputo_refinement.txt     # from src/python

Runs caputo-decay (L = 160, Nx = 512, T = 15, automatic r) for nu = 0.9, 0.8, 0.7, 0.6 with
N = 600, 1200, 2400, 4800 and prints P(T)/P(0), the peak intensity at T, the observed order of the
successive differences and a Richardson estimate of the limit. N = 600 and 1200 coincide with the PHP
values of results/reproduce.txt; N = 2400 and 4800 exceed the default work budget of the web instance.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from manakov_lab import __version__, run_preset  # noqa: E402

NS = (600, 1200, 2400, 4800)


def main() -> None:
    print(f"caputo-decay time refinement, Manakov Lab {__version__} (Python)")
    for nu in (0.9, 0.8, 0.7, 0.6):
        P, pk = [], []
        for N in NS:
            r = run_preset("caputo-decay", params=dict(nu=nu, N=N))
            d = r["diagnostics"]
            tot = d["P1"] + d["P2"]
            P.append(float(tot[-1] / tot[0]))
            pk.append(float(d["peak"][-1]))
        for name, v in (("P(T)/P(0)", P), ("peak", pk)):
            dif = [v[i] - v[i + 1] for i in range(3)]
            p = math.log2(dif[1] / dif[2]) if dif[1] * dif[2] > 0 else float("nan")
            lim = v[-1] - dif[2] / (2 ** p - 1) if p == p and p > 0 else float("nan")
            print(f"nu={nu}  {name:10s} " + "  ".join(f"N={N}: {x:.6f}" for N, x in zip(NS, v))
                  + f"  | order {math.log2(dif[0] / dif[1]):.2f}, {p:.2f}  | Richardson limit {lim:.6f}"
                  + f"  | rel. diff N=1200 vs limit {abs(v[1] - lim) / abs(lim):.2%}")


if __name__ == "__main__":
    main()
