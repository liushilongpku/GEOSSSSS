#!/usr/bin/env python3
"""Compare the GEOS single-porosity case with the classical Mandel solution."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "geos_mandel_mplconfig"))

import h5py
import matplotlib
import numpy as np
from scipy.optimize import brentq

matplotlib.use("Agg")
import matplotlib.pyplot as plt


K = 1.1e9
G = 7.5738e8
NU = 0.22
KS = 27.0e9
KF = 1.744e9
PHI = 0.14
PERMEABILITY = 4.9346e-21
VISCOSITY = 1.0e-3
A = 0.03
PC = 1.0e6

VALIDATION_DIR = Path(__file__).resolve().parent.parent
DEFAULT_RESULT_DIR = VALIDATION_DIR / "result"
DEFAULT_GEOS_OUTPUT_DIR = DEFAULT_RESULT_DIR / "single_porosity_geos"

ALPHA = 1.0 - K / KS
SKEMPTON_B = (1.0 / K - 1.0 / KS) / (
    1.0 / K - 1.0 / KS + PHI * (1.0 / KF - 1.0 / KS)
)
NU_U = (3.0 * NU + SKEMPTON_B * ALPHA * (1.0 - 2.0 * NU)) / (
    3.0 - SKEMPTON_B * ALPHA * (1.0 - 2.0 * NU)
)
C_DIFF = (
    2.0 * PERMEABILITY * G * (1.0 - NU) * (NU_U - NU)
    / (VISCOSITY * ALPHA**2 * (1.0 - 2.0 * NU) ** 2 * (1.0 - NU_U))
)
P0 = PC * SKEMPTON_B * (1.0 + NU_U) / 3.0
T0 = A**2 / C_DIFF


def mandel_roots(count: int = 80) -> np.ndarray:
  coefficient = (1.0 - NU) / (NU_U - NU)
  roots = []
  for index in range(count):
    left = index * np.pi + 1.0e-10
    right = (index + 0.5) * np.pi - 1.0e-10
    roots.append(brentq(lambda value: np.tan(value) - coefficient * value, left, right))
  return np.asarray(roots)


ROOTS = mandel_roots()
DENOMINATOR = ROOTS - np.sin(ROOTS) * np.cos(ROOTS)


def analytical_pressure(time: np.ndarray, x: float) -> np.ndarray:
  decay = np.exp(-np.outer(time / T0, ROOTS**2))
  spatial = np.cos(ROOTS * x / A) - np.cos(ROOTS)
  coefficients = np.sin(ROOTS) * spatial / DENOMINATOR
  return 2.0 * P0 * decay @ coefficients


def load_history(path: Path, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  with h5py.File(path, "r") as h5:
    time = np.asarray(h5[f"{field} Time"]).squeeze()
    values = np.asarray(h5[field])
    centers = np.asarray(h5[f"{field} elementCenter"])[0]
  return time, values, centers


def nearest_center(centers: np.ndarray) -> int:
  target = np.asarray([A / 40.0, A / 2.0, A / 2.0])
  return int(np.argmin(np.linalg.norm(centers - target, axis=1)))


def interpolate_reference(time: np.ndarray, analytical_time: np.ndarray, values: np.ndarray) -> np.ndarray:
  return np.interp(np.log10(time), np.log10(analytical_time), values)


def metrics(numerical: np.ndarray, reference: np.ndarray, scale: float) -> tuple[float, float]:
  error = np.abs(numerical - reference) / scale
  return float(np.sqrt(np.mean(error**2))), float(np.max(error))


def main() -> None:
  parser = argparse.ArgumentParser(description="Compare GEOS single-porosity pressure with the Mandel solution.")
  parser.add_argument("output_dir", nargs="?", type=Path, default=DEFAULT_GEOS_OUTPUT_DIR,
                      help=f"GEOS output directory (default: {DEFAULT_GEOS_OUTPUT_DIR})")
  parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULT_DIR,
                      help=f"comparison output directory (default: {DEFAULT_RESULT_DIR})")
  args = parser.parse_args()

  output_dir = args.output_dir.resolve()
  out_dir = args.out_dir.resolve()
  out_dir.mkdir(parents=True, exist_ok=True)

  tp, pressure, pressure_centers = load_history(output_dir / "pressure_history.hdf5", "pressure")

  pressure_index = nearest_center(pressure_centers)
  x_pressure = float(pressure_centers[pressure_index, 0])

  p_geos = np.asarray(pressure[:, pressure_index]).squeeze()

  positive_times = tp[tp > 0.0]
  analytical_time = np.logspace(np.log10(positive_times.min()), np.log10(positive_times.max()), 600)
  p_analytical = analytical_pressure(analytical_time, x_pressure)

  pressure_mask = tp > 0.0
  p_reference = interpolate_reference(tp[pressure_mask], analytical_time, p_analytical)
  p_rmse, p_max = metrics(p_geos[pressure_mask], p_reference, P0)

  plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
  fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
  tau_geos = tp[pressure_mask] / T0
  ax.semilogx(analytical_time / T0, p_analytical / P0, color="#202020", linewidth=1.7,
              label="Mandel analytical")
  marker_step = max(1, len(tau_geos) // 35)
  ax.semilogx(tau_geos, p_geos[pressure_mask] / P0, color="#B33A3A", linewidth=1.0, marker="o",
              markevery=marker_step, markersize=3.0, markerfacecolor="white", label="GEOS single medium")
  ax.set_xlabel(r"Dimensionless time $\tau=c t/a^2$")
  ax.set_ylabel(r"Center pressure $p/P_0^+$")
  ax.set_title("Single-porosity GEOS versus Mandel analytical solution")
  ax.legend(frameon=False)

  png_path = out_dir / "single_porosity_vs_mandel.png"
  fig.savefig(png_path, dpi=300)
  plt.close(fig)

  csv_path = out_dir / "single_porosity_vs_mandel_metrics.csv"
  with csv_path.open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("quantity", "normalized_rmse", "normalized_max_abs_error"))
    writer.writerow(("pressure", p_rmse, p_max))

  summary_path = out_dir / "comparison_summary.txt"
  summary_path.write_text(
      "Single-porosity GEOS versus classical Mandel solution\n"
      f"GEOS output: {output_dir}\n"
      f"K={K:.6e} Pa, G={G:.6e} Pa, nu={NU:.6f}, alpha={ALPHA:.8f}\n"
      f"B={SKEMPTON_B:.8f}, nu_u={NU_U:.8f}, c={C_DIFF:.8e} m^2/s, t0={T0:.8e} s\n"
      f"P0+={P0:.8e} Pa, Pc={PC:.8e} Pa\n"
      f"Pressure normalized RMSE={p_rmse:.8e}, max error={p_max:.8e}\n"
  )

  print(f"PNG: {png_path}")
  print(f"Metrics: {csv_path}")
  print(f"Summary: {summary_path}")


if __name__ == "__main__":
  main()
