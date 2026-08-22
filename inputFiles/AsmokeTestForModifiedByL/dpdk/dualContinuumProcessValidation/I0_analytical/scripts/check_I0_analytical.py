#!/usr/bin/env python3
# Purpose: verify the analytical I0 matrix water-saturation relaxation rate and closed-system conservation.
"""Check the analytical I0 capillary-exchange exponential solution."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


SEQ = 0.6
MW = 0.99 / 1.0e-3
MG = 0.99 / 2.0e-5
BG = MW / (MW + MG)
TB = 1.0e-12 * 4.0 * (1.0 / (0.3 ** 2) * 3.0)
SLOPE = 5.0e5 / SEQ
PHI = 0.2
LAMBDA = TB * MW * (1.0 - BG) * SLOPE / PHI


def load(path: Path, field: str) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5[f"{field} Time"])[:, 0], np.asarray(h5[field])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "runs" / "20260822_analytic")
    parser.add_argument("--max-fit-error", type=float, default=2.0e-3)
    parser.add_argument("--max-mass-drift", type=float, default=1.0e-8)
    args = parser.parse_args()

    t, matrix_phase = load(args.run_dir / "I0a_matrix_phase_history.hdf5", "phaseVolumeFraction")
    tf, fracture_phase = load(args.run_dir / "I0a_fracture_phase_history.hdf5", "phaseVolumeFraction")
    if not np.allclose(t, tf):
        raise AssertionError("matrix and fracture histories use different time grids")

    swm = matrix_phase[:, 0, 1]
    swf = fracture_phase[:, 0, 1]
    sum_sw = swm[0] + swf[0]
    eq_f = sum_sw - SEQ

    if not np.isfinite(swm).all() or not np.isfinite(swf).all():
        raise AssertionError("phase history contains NaN or Inf")

    # water mass conservation: psi*rho_w*sum(Sw) is constant (psi=phi, both cells equal volume)
    water_total = swm + swf
    mass_drift = abs(water_total[-1] - water_total[0]) / water_total[0]

    # matrix water saturation approaches capillary equilibrium as a single exponential
    pred_m = SEQ - (SEQ - swm[0]) * np.exp(-LAMBDA * t)
    fit_error_m = float(np.max(np.abs(pred_m - swm) / np.maximum(np.abs(swm), 1.0e-9)))
    # fracture follows from total water conservation
    pred_f = sum_sw - pred_m
    fit_error_f = float(np.max(np.abs(pred_f - swf) / np.maximum(np.abs(swf), 1.0e-9)))
    # capillary-equilibrium limit: Swm -> Seq, Swf -> sum_sw - Seq
    eq_error_m = abs(swm[-1] - SEQ) / SEQ
    eq_error_f = abs(swf[-1] - eq_f) / eq_f

    if not np.all(swm[1:] >= swm[:-1]):
        raise AssertionError("matrix water saturation is not monotonically increasing")
    if not np.all(swf[1:] <= swf[:-1]):
        raise AssertionError("fracture water saturation is not monotonically decreasing")
    if fit_error_m >= args.max_fit_error:
        raise AssertionError(f"matrix fit criteria failed: fit_error_m={fit_error_m:.3e}")
    if fit_error_f >= args.max_fit_error:
        raise AssertionError(f"fracture fit criteria failed: fit_error_f={fit_error_f:.3e}")
    if mass_drift >= args.max_mass_drift:
        raise AssertionError(f"conservation criteria failed: mass_drift={mass_drift:.3e}")

    print(f"I0a lambda_ref={LAMBDA:.9f} s^-1, tau={1.0 / LAMBDA:.6f} s")
    print(f"I0a matrix_exp_fit_error={fit_error_m:.3e}")
    print(f"I0a fracture_exp_fit_error={fit_error_f:.3e}")
    print(f"I0a water_mass_relative_drift={mass_drift:.3e}")
    print(f"I0a final_matrix_Sw={swm[-1]:.6f}, equilibrium_Seq={SEQ:.6f}")
    print(f"I0a final_fracture_Sw={swf[-1]:.6f}, equilibrium_Swf={eq_f:.6f}")
    print(f"I0a matrix_equilibrium_error={eq_error_m:.3e}, fracture_equilibrium_error={eq_error_f:.3e}")


if __name__ == "__main__":
    main()
