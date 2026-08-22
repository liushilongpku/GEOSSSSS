#!/usr/bin/env python3
# Purpose: plot I0 capillary state, imbibition response, and component conservation.
"""Generate the isolated I0 capillary/imbibition figure."""

from __future__ import annotations

import os
import tempfile
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dual_continuum_validation_mplconfig"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = CASE_ROOT / "runs" / "20260821"
FIGURE_ROOT = CASE_ROOT / "figures"
BLUE, ORANGE, PURPLE, GREEN, BLACK = "#0072B2", "#D55E00", "#7B2CBF", "#009E73", "#222222"


def load(path: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5[name]), np.asarray(h5[f"{name} Time"])[:, 0]


def style(ax: plt.Axes) -> None:
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    matrix_pc, time = load(RUN_ROOT / "I0_matrix_cap_history.hdf5", "matrixCapPressure_phaseCapPressure")
    fracture_pc, time_f = load(RUN_ROOT / "I0_fracture_cap_history.hdf5", "fractureCapPressure_phaseCapPressure")
    matrix_phase, time_pm = load(RUN_ROOT / "I0_matrix_phase_history.hdf5", "phaseVolumeFraction")
    fracture_phase, time_pf = load(RUN_ROOT / "I0_fracture_phase_history.hdf5", "phaseVolumeFraction")
    matrix_amount, time_am = load(RUN_ROOT / "I0_matrix_amount_history.hdf5", "compAmount")
    fracture_amount, time_af = load(RUN_ROOT / "I0_fracture_amount_history.hdf5", "compAmount")
    if not all(np.allclose(time, other) for other in (time_f, time_pm, time_pf, time_am, time_af)):
        raise ValueError("I0 histories use different time grids")
    inventory_drift = (matrix_amount + fracture_amount - matrix_amount[0] - fracture_amount[0]) / (matrix_amount[0] + fracture_amount[0])
    fig = plt.figure(figsize=(10.0, 8.35), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.46])
    axes = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
            [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]
    ax = axes[0][0]
    ax.plot(time, matrix_pc[:, 0, 1] / 1e3, color=BLUE, label="Matrix water Pc")
    ax.plot(time_f, fracture_pc[:, 0, 1] / 1e3, color=ORANGE, ls="--", label="Fracture water Pc")
    ax.set(xlabel="Time (s)", ylabel="Water capillary pressure (kPa)", title="Capillary state")
    ax.legend(frameon=False)
    style(ax)
    ax = axes[0][1]
    ax.plot(time_pm, matrix_phase[:, 0, 0], color=BLUE, label="Matrix gas volume fraction")
    ax.plot(time_pf, fracture_phase[:, 0, 0], color=ORANGE, ls="--", label="Fracture gas volume fraction")
    ax.set(xlabel="Time (s)", ylabel="Gas volume fraction (-)", title="Imbibition response")
    ax.legend(frameon=False)
    style(ax)
    ax = axes[1][0]
    ax.plot(time_pm, matrix_phase[:, 0, 1], color=BLUE, label="Matrix water volume fraction")
    ax.plot(time_pf, fracture_phase[:, 0, 1], color=ORANGE, ls="--", label="Fracture water volume fraction")
    ax.set(xlabel="Time (s)", ylabel="Water volume fraction (-)", title="Complementary phase response")
    ax.legend(frameon=False)
    style(ax)
    ax = axes[1][1]
    ax.plot(time_am, inventory_drift[:, 0, 0], color=PURPLE, label="CO2 inventory drift")
    ax.plot(time_am, inventory_drift[:, 0, 1], color=GREEN, ls="--", label="Water inventory drift")
    ax.axhline(0.0, color=BLACK, lw=0.8)
    ax.set(xlabel="Time (s)", ylabel="Relative total inventory drift (-)", title="Closed-system conservation")
    ax.legend(frameon=False)
    style(ax)
    info = fig.add_subplot(grid[2, :])
    info.axis("off")
    values = [
        ("Reference", "Closed compositional system; gravity and external mass boundaries are off.\n"
                      "Matrix water Pc is nonzero; fracture water Pc is the zero reference."),
        ("This run", f"Matrix water Pc: {matrix_pc[0, 0, 1]/1e3:.2f} -> {matrix_pc[-1, 0, 1]/1e3:.2f} kPa.\n"
                      f"Gas fraction: matrix {matrix_phase[0, 0, 0]:.5f}->{matrix_phase[-1, 0, 0]:.5f}; "
                      f"fracture {fracture_phase[0, 0, 0]:.5f}->{fracture_phase[-1, 0, 0]:.5f}."),
        ("Physical picture", "Capillary suction draws water from the fracture into the matrix.\n"
                            "Matrix gas decreases, fracture gas increases, and total component\n"
                            "inventories remain constant."),
    ]
    for index, (heading, body) in enumerate(values):
        left = index / len(values) + 0.012
        info.text(left, 0.95, heading, transform=info.transAxes, va="top", fontsize=9, fontweight="bold", color=BLACK)
        wrapped = "\n".join(textwrap.fill(line, width=27, break_long_words=False) for line in body.splitlines())
        info.text(left, 0.72, wrapped, transform=info.transAxes, va="top", fontsize=7.2, linespacing=1.25,
                  bbox={"boxstyle": "round,pad=0.35", "facecolor": "#F5F5F5", "edgecolor": "#D0D0D0"})
    FIGURE_ROOT.mkdir(exist_ok=True)
    metadata = {"Title": "I0 capillary exchange", "Subject": "GEOS capillary imbibition", "Creator": "plot_I0.py"}
    fig.savefig(FIGURE_ROOT / "I0_capillary_exchange.png", dpi=300, metadata=metadata)
    fig.savefig(FIGURE_ROOT / "I0_capillary_exchange.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {FIGURE_ROOT / 'I0_capillary_exchange.png'}")


if __name__ == "__main__":
    main()
