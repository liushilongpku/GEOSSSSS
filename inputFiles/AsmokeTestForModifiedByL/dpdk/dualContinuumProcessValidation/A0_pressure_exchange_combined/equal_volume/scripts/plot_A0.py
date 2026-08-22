#!/usr/bin/env python3
# Purpose: plot the A0 pressure-exchange history against its two analytical references.
"""Generate the isolated A0 pressure-exchange figure."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dual_continuum_validation_mplconfig"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = CASE_ROOT / "runs" / "20260821_vf05_flow"
FIGURE_ROOT = CASE_ROOT / "figures"
BLUE, ORANGE, BLACK, GRAY = "#0072B2", "#D55E00", "#222222", "#666666"


def load_history(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5["pressure Time"])[:, 0], np.asarray(h5["pressure"])[:, 0]


def style_axes(ax: plt.Axes) -> None:
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.8)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#222222")
        spine.set_linewidth(0.9)


def main() -> None:
    time_m, matrix = load_history(RUN_ROOT / "A0_matrix_pressure_history.hdf5")
    time_f, fracture = load_history(RUN_ROOT / "A0_fracture_pressure_history.hdf5")
    if not np.allclose(time_m, time_f):
        raise ValueError("A0 matrix and fracture histories use different time grids")

    dt, storage, exchange = 0.02, 4.0e-10, 1.2e-9
    q = 1.0 / (1.0 + 2.0 * exchange * dt / storage)
    tau = storage / (2.0 * exchange)
    raw_mean, raw_delta = 1.5e6, 1.0e6
    reference_time = time_m + dt
    delta_be = raw_delta * q ** (reference_time / dt)
    delta_exp = raw_delta * np.exp(-reference_time / tau)
    matrix_be, fracture_be = raw_mean + 0.5 * delta_be, raw_mean - 0.5 * delta_be
    matrix_exp, fracture_exp = raw_mean + 0.5 * delta_exp, raw_mean - 0.5 * delta_exp
    q_num = float((matrix[1] - fracture[1]) / (matrix[0] - fracture[0]))
    q_error = abs(q_num - q) / q
    resolved = np.abs(delta_be) > 100.0 * np.finfo(float).eps * raw_mean
    curve_error = float(np.max(np.abs((matrix - fracture)[resolved] - delta_be[resolved]) / delta_be[resolved]))

    # Near-square 6.2:6.0 inch canvas: 1860x1800 px at 300 dpi.
    fig, ax = plt.subplots(figsize=(6.2, 6.0), constrained_layout=True)
    ax.plot(time_m, matrix / 1e6, color=BLUE, marker="o", markevery=10, ms=3, label="Matrix GEOS")
    ax.plot(time_f, fracture / 1e6, color=ORANGE, marker="s", markevery=10, ms=3, label="Fracture GEOS")
    ax.plot(time_m, matrix_be / 1e6, color=BLUE, ls="--", lw=1, label="Matrix backward-Euler reference")
    ax.plot(time_f, fracture_be / 1e6, color=ORANGE, ls="--", lw=1, label="Fracture backward-Euler reference")
    ax.plot(time_m, matrix_exp / 1e6, color=BLUE, ls=":", lw=1, label="Matrix continuous reference")
    ax.plot(time_f, fracture_exp / 1e6, color=ORANGE, ls=":", lw=1, label="Fracture continuous reference")
    ax.set(xlabel="Time (s)", ylabel="Pressure (MPa)", title="A0 pressure-driven exchange")
    ax.legend(frameon=False, ncols=1, fontsize=8, loc="upper right")
    style_axes(ax)
    ax.text(0.02, 0.04, r"Reference: $p_m=1.5+0.5e^{-6t}$ MPa; $p_f=1.5-0.5e^{-6t}$ MPa",
            transform=ax.transAxes, fontsize=8, color=BLACK,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.9})
    FIGURE_ROOT.mkdir(exist_ok=True)
    metadata = {"Title": "pressure exchange", "Subject": "GEOS and analytical references", "Creator": "plot_A0.py"}
    fig.savefig(FIGURE_ROOT / "A0_pressure_exchange.png", dpi=300, metadata=metadata)
    fig.savefig(FIGURE_ROOT / "A0_pressure_exchange.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {FIGURE_ROOT / 'A0_pressure_exchange.png'}")


if __name__ == "__main__":
    main()
