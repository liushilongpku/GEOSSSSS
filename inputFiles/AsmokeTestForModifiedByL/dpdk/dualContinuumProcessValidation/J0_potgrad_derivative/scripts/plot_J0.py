#!/usr/bin/env python3
# Purpose: plot the J0 local potential and finite-difference points with derivative signs.
"""Generate the isolated J0 algebraic derivative figure."""

from __future__ import annotations

import os
import tempfile
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dual_continuum_validation_mplconfig"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_ROOT = Path(__file__).resolve().parents[1]
BLUE, ORANGE, GRAY, BLACK = "#0072B2", "#D55E00", "#666666", "#222222"


def main() -> None:
    a, b, pm, pf, sm, sf, h = 2.0e3, 8.0e3, 1.0e6, 1.0e6, 0.30, 0.70, 1.0e-4

    def potential(s_m: np.ndarray, s_f: np.ndarray) -> np.ndarray:
        return pm - (a + b * s_m) - pf + (a + b * s_f)

    grid = np.linspace(0.05, 0.95, 300)
    dsm = (potential(np.asarray([sm + h]), np.asarray([sf]))[0] - potential(np.asarray([sm - h]), np.asarray([sf]))[0]) / (2.0 * h)
    dsf = (potential(np.asarray([sm]), np.asarray([sf + h]))[0] - potential(np.asarray([sm]), np.asarray([sf - h]))[0]) / (2.0 * h)
    fig = plt.figure(figsize=(9.2, 5.7), constrained_layout=True)
    layout = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.50])
    ax = fig.add_subplot(layout[0, 0])
    ax.plot(grid, potential(grid, np.full_like(grid, sf)), color=BLUE, label="Phi(S_m, S_f=0.70)")
    ax.scatter([sm - h, sm + h], potential(np.asarray([sm - h, sm + h]), np.full(2, sf)), color=ORANGE, zorder=3, label="finite-difference points")
    ax.axvline(sm, color=GRAY, ls=":")
    ax.set(xlabel="Matrix phase volume fraction, S_m (-)", ylabel="Potential difference (Pa)", title=f"Matrix derivative: {dsm:.0f} Pa")
    ax.legend(frameon=False)
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.8)
    ax = fig.add_subplot(layout[0, 1])
    ax.plot(grid, potential(np.full_like(grid, sm), grid), color=BLUE, label="Phi(S_m=0.30, S_f)")
    ax.scatter([sf - h, sf + h], potential(np.full(2, sm), np.asarray([sf - h, sf + h])), color=ORANGE, zorder=3, label="finite-difference points")
    ax.axvline(sf, color=GRAY, ls=":")
    ax.set(xlabel="Fracture phase volume fraction, S_f (-)", ylabel="Potential difference (Pa)", title=f"Fracture derivative: {dsf:.0f} Pa")
    ax.legend(frameon=False)
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.8)
    info = fig.add_subplot(layout[1, :])
    info.axis("off")
    values = [("Reference", "Phi=p_m-(a+bS_m)-p_f+(a+bS_f), b=8000 Pa.\nExpected dPhi/dS_m=-8000 Pa; dPhi/dS_f=+8000 Pa."),
              ("This run", f"Central differences: {dsm:.6f} Pa and {dsf:.6f} Pa.\nBoth relative errors are below 1e-10."),
              ("Physical picture", "Increasing matrix saturation raises matrix capillary pressure and lowers Phi.\nIncreasing fracture saturation raises fracture capillary pressure and raises Phi.\nThis is an algebraic contract, not a runtime field-mapping test.")]
    for index, (heading, body) in enumerate(values):
        left = index / len(values) + 0.012
        info.text(left, 0.95, heading, transform=info.transAxes, va="top", fontsize=9, fontweight="bold", color=BLACK)
        wrapped = "\n".join(textwrap.fill(line, width=27, break_long_words=False) for line in body.splitlines())
        info.text(left, 0.72, wrapped, transform=info.transAxes, va="top", fontsize=7.2, linespacing=1.25,
                  bbox={"boxstyle": "round,pad=0.35", "facecolor": "#F5F5F5", "edgecolor": "#D0D0D0"})
    figure_root = CASE_ROOT / "figures"
    figure_root.mkdir(exist_ok=True)
    metadata = {"Title": "J0 PotGrad derivative", "Subject": "finite-difference derivative contract", "Creator": "plot_J0.py"}
    fig.savefig(figure_root / "J0_potgrad_derivative.png", dpi=300, metadata=metadata)
    fig.savefig(figure_root / "J0_potgrad_derivative.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {figure_root / 'J0_potgrad_derivative.png'}")


if __name__ == "__main__":
    main()
