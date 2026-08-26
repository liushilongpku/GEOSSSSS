#!/usr/bin/env python3
# Purpose: plot the G0 gravity-exchange pressure contrast against GDP equilibrium.
"""Generate the isolated G0 gravity-drainage figure."""

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
RUN_ROOT = CASE_ROOT / "runs" / "20260822_dt001"
FIGURE_ROOT = CASE_ROOT / "figures"
BLUE, RED, BLACK = "#0072B2", "#C62828", "#222222"


def load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5["pressure Time"])[:, 0], np.asarray(h5["pressure"])[:, 0]


def style(ax: plt.Axes) -> None:
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.8)


def main() -> None:
    time_m, matrix = load(RUN_ROOT / "G0_matrix_pressure_history.hdf5")
    time_f, fracture = load(RUN_ROOT / "G0_fracture_pressure_history.hdf5")
    if not np.allclose(time_m, time_f):
        raise ValueError("G0 matrix and fracture histories use different time grids")
    gdp = 147150.0
    contrast = matrix - fracture
    dt = 0.01
    lam = 7.125
    q = 1.0 / (1.0 + lam * dt)
    step = np.arange(len(time_m), dtype=float) + 1.0
    contrast_be = -gdp * (1.0 - q ** step)
    fig, ax = plt.subplots(figsize=(6.4, 6.2), constrained_layout=True)
    ax.plot(time_m, contrast / 1e3, color=BLUE, marker="o", markevery=10, ms=3,
            lw=1.8, zorder=3, label="GEOS result: $p_m-p_f$")
    ax.plot(time_m, contrast_be / 1e3, color=RED, ls=(0, (7, 3)), lw=3.2, zorder=2,
            label="Backward-Euler transient benchmark")
    ax.axhline(-gdp / 1e3, color=BLACK, ls=(0, (2, 2)), lw=2.2, zorder=1,
               label="Final GDP equilibrium: $-147.15$ kPa")
    ax.set(xlabel="Time (s)", ylabel="Pressure contrast (kPa)", title="G0 gravity-driven exchange")
    ax.set_xlim(left=0.0)
    ax.set_ylim(-152.0, 0.0)
    ax.axhline(0.0, color="#AAAAAA", linewidth=0.8, alpha=0.55)
    ax.axvline(0.0, color="#AAAAAA", linewidth=0.8, alpha=0.55)
    t_target = 0.32
    xy = (t_target, contrast_be[np.argmin(np.abs(time_m - t_target))] / 1e3)
    ax.annotate("Transient benchmark\n$\\Delta p_n=-P_{GDP}(1-q^n)$\n$q=0.9334889$",
                xy=xy, xycoords="data",
                xytext=(0.60, -42.0), textcoords="data",
                color=RED, fontsize=9, ha="left", va="top",
                arrowprops={"arrowstyle": "-", "color": RED, "lw": 1.2},
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                      "edgecolor": RED, "alpha": 0.95})
    ax.annotate("Final GDP equilibrium:\n$\\Delta p_{eq}=-P_{GDP}=-147.15$ kPa",
                xy=(1.55, -gdp / 1e3), xycoords="data",
                xytext=(1.18, -140.5), textcoords="data",
                color=BLACK, fontsize=9, ha="left", va="bottom",
                arrowprops={"arrowstyle": "-", "color": BLACK, "lw": 1.2},
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                      "edgecolor": BLACK, "alpha": 0.95})
    ax.legend(frameon=True, facecolor="white", edgecolor="#D0D0D0", framealpha=0.92,
              loc="upper right", fontsize=9)
    style(ax)
    FIGURE_ROOT.mkdir(exist_ok=True)
    metadata = {"Title": "G0 gravity exchange", "Subject": "GEOS versus transient GDP benchmark", "Creator": "plot_G0.py"}
    fig.savefig(FIGURE_ROOT / "G0_gravity_exchange.png", dpi=300, metadata=metadata)
    fig.savefig(FIGURE_ROOT / "G0_gravity_exchange.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {FIGURE_ROOT / 'G0_gravity_exchange.png'}")


if __name__ == "__main__":
    main()
