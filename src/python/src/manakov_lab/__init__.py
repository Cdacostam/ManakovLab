"""Manakov Lab (Python reference implementation).

A line-by-line port of the PHP numerical core of Manakov Lab: Strang split-step Fourier for the
variable-coefficient and Riesz space-fractional Manakov systems, the L1 scheme on a graded mesh for the
Caputo time-fractional system, Petviashvili iteration, exact and manufactured solutions and convergence
studies. It is used to cross-check the PHP results (see ``crosscheck.py``).
"""
__version__ = "1.0.0"

from .core import (CaputoSolver, SimulationError, SplitStepSolver, fractional_soliton, grid,  # noqa: E402
                   manakov)
from .presets import PRESETS, problem  # noqa: E402
from .simulation import Simulation  # noqa: E402


def run_preset(preset_id: str, **overrides) -> dict:
    """Run a preset of the web interface, with optional overrides (see :func:`problem`)."""
    return Simulation(problem(preset_id, **overrides)).run()


__all__ = ["CaputoSolver", "SimulationError", "SplitStepSolver", "Simulation", "PRESETS", "problem",
           "run_preset", "fractional_soliton", "grid", "manakov", "__version__"]
