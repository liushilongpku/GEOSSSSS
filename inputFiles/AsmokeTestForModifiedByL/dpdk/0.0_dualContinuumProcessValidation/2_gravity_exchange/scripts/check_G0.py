#!/usr/bin/env python3
# Purpose: check the G0 signed GDP equilibrium from the archived pressure histories.
"""Independent quantitative checks for the G0 archived history."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def load(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5["pressure"])[:, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parents[1] / "runs" / "20260821")
    args = parser.parse_args()
    matrix = load(args.run_dir / "G0_matrix_pressure_history.hdf5")
    fracture = load(args.run_dir / "G0_fracture_pressure_history.hdf5")
    gdp = (-9.81) * (800.0 - 1100.0) * 100.0 / 2.0
    final_contrast = matrix[-1] - fracture[-1]
    absolute_error = abs(final_contrast + gdp)
    relative_error = absolute_error / gdp
    if final_contrast >= 0.0 or absolute_error >= 1.0 or relative_error >= 1.0e-5:
        raise AssertionError(f"G0 criteria failed: contrast={final_contrast:g}, error={absolute_error:g}")
    print(f"G0 GDP={gdp:.6f} Pa, final_contrast={final_contrast:.6f} Pa")
    print(f"G0 absolute_error={absolute_error:.6e} Pa, relative_error={relative_error:.3e}")


if __name__ == "__main__":
    main()
