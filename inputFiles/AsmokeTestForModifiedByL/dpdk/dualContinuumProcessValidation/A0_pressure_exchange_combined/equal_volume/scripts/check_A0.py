#!/usr/bin/env python3
# Purpose: check A0 pressure-exchange decay, conservation, and analytical agreement.
"""Independent quantitative checks for the A0 archived history."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5["pressure Time"])[:, 0], np.asarray(h5["pressure"])[:, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parents[1] / "runs" / "20260821_vf05_flow")
    args = parser.parse_args()
    tm, pm = load(args.run_dir / "A0_matrix_pressure_history.hdf5")
    tf, pf = load(args.run_dir / "A0_fracture_pressure_history.hdf5")
    if not np.allclose(tm, tf):
        raise AssertionError("matrix and fracture time grids differ")
    delta = pm - pf
    dt, storage, exchange = 0.02, 4.0e-10, 1.2e-9
    q = 1.0 / (1.0 + 2.0 * exchange * dt / storage)
    q_num = delta[1] / delta[0]
    q_error = abs(q_num - q) / q
    delta_be = delta[0] * q ** (np.arange(delta.size))
    curve_error = np.max(np.abs(delta - delta_be) / np.maximum(np.abs(delta_be), 1.0))
    mean = 0.5 * (pm + pf)
    mean_drift = np.max(np.abs(mean - 1.5e6)) / 1.5e6
    if not np.all(delta[1:] < delta[:-1]):
        raise AssertionError("pressure contrast is not monotonically decreasing")
    if q_error >= 1.0e-5 or curve_error >= 1.0e-5 or mean_drift >= 1.0e-6:
        raise AssertionError(f"A0 criteria failed: q_error={q_error:g}, curve_error={curve_error:g}, "
                             f"mean_drift={mean_drift:g}")
    print(f"A0 q_num={q_num:.10f}, q_ref={q:.10f}, relative_q_error={q_error:.3e}")
    print(f"A0 relative_BE_curve_error={curve_error:.3e}")
    print(f"A0 relative_mean_pressure_drift={mean_drift:.3e}, final_contrast={delta[-1]:.6g} Pa")


if __name__ == "__main__":
    main()
