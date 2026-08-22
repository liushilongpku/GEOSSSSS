#!/usr/bin/env python3
# Purpose: check the unequal-volume A0 pressure-exchange balance and decay factor.
"""Validate the 0.95/0.05 REV-volume pressure-exchange case."""

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
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parents[1] / "runs" / "latest")
    args = parser.parse_args()
    tm, pm = load(args.run_dir / "A0_unequal_matrix_pressure_history.hdf5")
    tf, pf = load(args.run_dir / "A0_unequal_fracture_pressure_history.hdf5")
    if not np.allclose(tm, tf):
        raise AssertionError("matrix and fracture time grids differ")

    vm, vf = 0.95, 0.05
    phi_m, phi_f = 0.2, 1.0
    fluid_c, porosity_c = 1.0e-12, 1.0e-12
    rho, mobility = 1000.0, 1.0e6
    dt, permeability, shape_factor = 0.02, 1.0e-12, 1.2e-3
    c0 = rho * (fluid_c + porosity_c)
    cm, cf = vm * phi_m * c0, vf * phi_f * c0
    equilibrium = (cm * 2.0e6 + cf * 1.0e6) / (cm + cf)
    transfer = vm * permeability * shape_factor * mobility
    lam = transfer * (1.0 / cm + 1.0 / cf)
    q = 1.0 / (1.0 + lam * dt)
    delta = pm - pf
    delta_be = delta[0] * q ** np.arange(delta.size)
    curve_error = np.max(np.abs(delta - delta_be) / np.maximum(np.abs(delta_be), 1.0))
    weighted = (cm * pm + cf * pf) / (cm + cf)
    weighted_drift = np.max(np.abs(weighted - equilibrium)) / equilibrium
    if not np.all(delta[1:] < delta[:-1]):
        raise AssertionError("pressure contrast is not monotonically decreasing")
    if curve_error >= 1.0e-5 or weighted_drift >= 1.0e-6:
        raise AssertionError(f"criteria failed: curve_error={curve_error:g}, weighted_drift={weighted_drift:g}")
    print(f"A0 unequal q_ref={q:.10f}, curve_error={curve_error:.3e}")
    print(f"A0 unequal equilibrium={equilibrium:.6g} Pa, weighted_drift={weighted_drift:.3e}")
    print(f"A0 unequal final_contrast={delta[-1]:.6g} Pa")


if __name__ == "__main__":
    main()
