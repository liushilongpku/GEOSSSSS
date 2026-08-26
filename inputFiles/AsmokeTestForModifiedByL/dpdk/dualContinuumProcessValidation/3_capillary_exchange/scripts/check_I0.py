#!/usr/bin/env python3
# Purpose: check I0 capillary-driven phase response and closed-system inventory conservation.
"""Independent quantitative checks for the I0 archived history."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def load(path: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5[name]), np.asarray(h5[f"{name} Time"])[:, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parents[1] / "runs" / "20260821")
    parser.add_argument("--max-mass-drift", type=float, default=1.0e-6)
    args = parser.parse_args()
    matrix_amount, time = load(args.run_dir / "I0_matrix_amount_history.hdf5", "compAmount")
    fracture_amount, fracture_time = load(args.run_dir / "I0_fracture_amount_history.hdf5", "compAmount")
    matrix_phase, _ = load(args.run_dir / "I0_matrix_phase_history.hdf5", "phaseVolumeFraction")
    fracture_phase, _ = load(args.run_dir / "I0_fracture_phase_history.hdf5", "phaseVolumeFraction")
    matrix_cap, _ = load(args.run_dir / "I0_matrix_cap_history.hdf5", "matrixCapPressure_phaseCapPressure")
    if not np.allclose(time, fracture_time):
        raise AssertionError("I0 matrix and fracture histories use different time grids")
    if not np.isfinite(matrix_amount).all() or not np.isfinite(fracture_amount).all():
        raise AssertionError("I0 component history contains NaN or Inf")
    total = matrix_amount[:, 0, :] + fracture_amount[:, 0, :]
    drift = (total[-1] - total[0]) / total[0]
    phase_change = max(float(np.max(np.abs(matrix_phase[-1, 0] - matrix_phase[0, 0]))),
                       float(np.max(np.abs(fracture_phase[-1, 0] - fracture_phase[0, 0]))))
    if np.max(np.abs(drift)) > args.max_mass_drift or phase_change <= 0.0 or matrix_cap[0, 0, 1] <= 0.0:
        raise AssertionError(f"I0 criteria failed: drift={drift}, phase_change={phase_change:g}")
    print(f"I0 final_time={time[-1]:.12g} s")
    print(f"I0 maximum_relative_inventory_drift={np.max(np.abs(drift)):.12e}")
    print(f"I0 maximum_phase_volume_fraction_change={phase_change:.12e}")
    print(f"I0 initial_matrix_water_capillary_pressure={matrix_cap[0, 0, 1]:.12e} Pa")


if __name__ == "__main__":
    main()
