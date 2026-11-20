"""Cross-check of the Python implementation against the PHP reference results.

    php ../php/tests/export_reference.php ../../results/php_reference.json   # PHP side (from src/python)
    python crosscheck.py ../../results/php_reference.json                   # comparison

For every case of src/php/tests/crosscheck_cases.json the script reruns the same problem in Python and
compares, with the stated normalisation, every quantity it reports:

  P1, P2  component powers, max over the output times of |P_j^py - P_j^php| / (P1 + P2)(0)
  peak    peak intensity, absolute difference (PHP rounds it to 1e-8)
  H       Hamiltonian, max |H^py - H^php| / max(1, max |H^php|)
  q       final complex field (real and imaginary parts, so phases are compared), max |q^py - q^php|
          (PHP rounds the exported parts to 1e-7)
  err     exact-solution error history, relative to its maximum (exact and manufactured cases)
  iter    nonlinear iterations of every step (Caputo cases, L1 scheme), required to be identical step by step;
          their mean and maximum are printed

A case passes only if every one of these is within its tolerance. The PHP values are rounded in the JSON
(powers and H to 1e-14, peak to 1e-8, field to 1e-7), which sets the attainable agreement.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from manakov_lab import run_preset  # noqa: E402

TOL = {"P1": 1e-9, "P2": 1e-9, "peak": 1e-6, "H": 1e-8, "q": 1e-6, "err": 1e-8}


def compare(case: dict) -> dict:
    r = run_preset(case["id"], **case["over"])
    d, ref = r["diagnostics"], case["diagnostics"]
    P0 = ref["P1"][0] + ref["P2"][0]
    q_php = np.array([np.array(re) + 1j * np.array(im) for re, im in case["final"]])
    out = {
        "P1": float(np.max(np.abs(d["P1"] - np.array(ref["P1"]))) / P0),
        "P2": float(np.max(np.abs(d["P2"] - np.array(ref["P2"]))) / P0),
        "peak": float(np.max(np.abs(d["peak"] - np.array(ref["peak"])))),
        "H": float(np.max(np.abs(d["H"] - np.array(ref["H"]))) / max(1.0, float(np.max(np.abs(ref["H"]))))),
        "q": float(np.max(np.abs(r["final"] - q_php))),
    }
    if case["err_t"] is not None:
        e_py, e_php = r["exact"]["err_t"], np.array(case["err_t"])
        out["err"] = float(np.max(np.abs(e_py - e_php) / np.max(e_php)))
    out["iter_py"] = [r["meta"]["picardMean"], r["meta"]["picardMax"]]
    out["iter_php"] = case["picard"]
    out["iter_steps"] = len(case["picardSteps"] or [])
    out["iter_ok"] = out["iter_py"] == out["iter_php"] and (r["meta"]["picardSteps"] or []) == (case["picardSteps"] or [])
    out["ok"] = out["iter_ok"] and all(out[k] <= TOL[k] for k in TOL if k in out)
    return out


def main(path: str) -> int:
    ref = json.loads(Path(path).read_text())
    bad = 0
    print(f"PHP reference: {ref['engine']}, Manakov Lab {ref['version']}")
    for case in ref["cases"]:
        c = compare(case)
        bad += not c["ok"]
        print(f"{'ok  ' if c['ok'] else 'FAIL'} {case['id']:<18} {json.dumps(case['over']):<62} "
              f"P1 {c['P1']:.1e}  P2 {c['P2']:.1e}  peak {c['peak']:.1e}  H {c['H']:.1e}  q {c['q']:.1e}"
              + (f"  err {c['err']:.1e}" if "err" in c else "")
              + (f"  iterations {'equal at all ' + str(c['iter_steps']) + ' steps' if c['iter_ok'] else 'DIFFER'} "
                 f"(mean, max {c['iter_py']})" if c["iter_php"][0] is not None else ""))
    print(f"{len(ref['cases']) - bad}/{len(ref['cases'])} cases agree within the tolerances {TOL}"
          " and with identical iteration counts at every step")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[2] / "results/php_reference.json")))
