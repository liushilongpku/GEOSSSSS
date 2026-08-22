#!/usr/bin/env python3
# Purpose: plot the analytical I0 water-saturation relaxation against the exponential reference.
"""Generate the analytical I0 capillary-imbibition figure."""

from __future__ import annotations

import os
import argparse
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dual_continuum_validation_mplconfig"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SEQ = 0.6
MW = 0.99 / 1.0e-3
MG = 0.99 / 2.0e-5
BG = MW / (MW + MG)
TB = 1.0e-12 * 4.0 * (1.0 / (0.3 ** 2) * 3.0)
SLOPE = 5.0e5 / SEQ
PHI = 0.2
LAMBDA = TB * MW * (1.0 - BG) * SLOPE / PHI

CASE_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = CASE_ROOT / "runs" / "20260822_analytic"
FIGURE_ROOT = CASE_ROOT / "figures"
BLUE, ORANGE, PURPLE, BLACK = "#0072B2", "#D55E00", "#7B2CBF", "#222222"
GREEN = "#009E73"

# Wetting (water) and non-wetting (gas) relative-permeability sample points (S in [0,1]).
RL_S = np.array([0.0, 0.2, 0.4, 1.0])
RL_KR = np.array([0.0, 0.9, 0.99, 1.0])
# Matrix capillary pressure: linear from 5e5 Pa at S=0 to 0 at S=0.6. Fracture is 0 everywhere.
CAP_S = np.array([0.0, 0.6])
CAP_MATRIX = np.array([5.0e5, 0.0])
CAP_FRACTURE = np.array([0.0, 0.0])


def load(path: Path, field: str) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5[f"{field} Time"])[:, 0], np.asarray(h5[field])


def plot_constitutive(figure_dir: Path) -> None:
    """Render the relative-permeability and capillary-pressure curves used by the case."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)

    ax = axes[0]
    ax.plot(RL_S, RL_KR, "o-", color=BLUE, lw=1.8, label="$k_{rw}(S_w)$")
    ax.plot(RL_S, RL_KR, "s--", color=ORANGE, lw=1.8, label="$k_{rg}(S_g)$")
    ax.axvspan(0.4, 0.6, color="#E8E8E8", alpha=0.7, linewidth=0)
    ax.set(xlabel="Phase saturation (-)", ylabel="Relative permeability (-)",
           title="Relative permeability")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8)
    ax.annotate("operative $S$ window", xy=(0.5, 0.03), ha="center", fontsize=7,
                xycoords="data", color=BLACK)

    ax = axes[1]
    ax.plot(CAP_S, CAP_MATRIX / 1.0e3, "o-", color=BLUE, lw=1.8, label="Matrix water $P_c$")
    ax.plot([0.0, 1.0], [0.0, 0.0], "s--", color=ORANGE, lw=1.8, label="Fracture water $P_c$")
    ax.axvline(SEQ, color=BLACK, ls=":", lw=0.9)
    ax.set(xlabel="Water saturation (-)", ylabel="Water capillary pressure (kPa)",
           title="Capillary pressure")
    ax.set_xlim(0, 1)
    ax.set_ylim(-20, 520)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.annotate("$S_{eq}=0.6$", xy=(SEQ, 480), ha="right", va="top", fontsize=7,
                color=BLACK)

    for axis in axes:
        axis.grid(True, color="#D9D9D9", linewidth=0.6)
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_color(BLACK)

    figure_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"Title": "I0 constitutive curves", "Creator": "plot_I0_analytical.py"}
    fig.savefig(figure_dir / "I0_constitutive_curves.png", dpi=300, metadata=metadata)
    fig.savefig(figure_dir / "I0_constitutive_curves.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {figure_dir / 'I0_constitutive_curves.png'}")


def main(args: argparse.Namespace) -> None:
    t, matrix_phase = load(args.run_dir / "I0a_matrix_phase_history.hdf5", "phaseVolumeFraction")
    tf, fracture_phase = load(args.run_dir / "I0a_fracture_phase_history.hdf5", "phaseVolumeFraction")
    if not np.allclose(t, tf):
        raise ValueError("matrix and fracture histories use different time grids")

    swm = matrix_phase[:, 0, 1]
    swf = fracture_phase[:, 0, 1]
    sum_sw = swm[0] + swf[0]
    # Analytical references: matrix relaxes as a single exponential; fracture follows from
    # total water conservation (Swm + Swf = constant). At capillary equilibrium Pc_m -> 0,
    # which fixes Swm -> Seq and Swf -> (Swm0 + Swf0) - Seq.
    pred_m = SEQ - (SEQ - swm[0]) * np.exp(-LAMBDA * t)
    pred_f = sum_sw - pred_m
    eq_f = sum_sw - SEQ

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.5), constrained_layout=True)

    ax1.plot(t, swm, "o-", ms=3, markevery=5, label="Matrix Sw GEOS", color=BLUE)
    ax1.plot(t, pred_m, "--", label="Matrix Sw analytical", color=PURPLE)
    ax1.axhline(SEQ, color=GREEN, lw=1.2, ls=":", label="Capillary equilibrium")
    ax1.set(xlabel="Time (s)", ylabel="Water saturation (-)",
            title="Matrix water saturation")
    ax1.legend(frameon=False, fontsize=9, loc="lower right")
    ax1.text(0.03, 0.90,
             rf"$S_w^m(t)=S_{{eq}}-(S_{{eq}}-S_w^m(0))e^{{-\lambda t}}$" + "\n" +
             rf"$\lambda={LAMBDA:.6f}\ \mathrm{{s^{{-1}}}},\ S_{{eq}}={SEQ:.3f}$",
             transform=ax1.transAxes, fontsize=9, va="top",
             bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D0D0"})

    ax2.plot(t, swf, "s-", ms=3, markevery=5, label="Fracture Sw GEOS", color=ORANGE)
    ax2.plot(t, pred_f, "--", label="Fracture Sw analytical", color=PURPLE)
    ax2.axhline(eq_f, color=GREEN, lw=1.0, ls=":", label="Capillary equilibrium")
    ax2.set(xlabel="Time (s)", title="Fracture water saturation")
    ax2.legend(frameon=False, fontsize=9, loc="upper right")
    ax2.text(0.04, 0.08,
             rf"$S_w^f(t)=S_{{tot}}-S_w^m(t)$" + "\n" +
             rf"$S_{{tot}}={sum_sw:.6f},\ S_{{eq}}^f={eq_f:.6f}$",
             transform=ax2.transAxes, fontsize=9, va="bottom",
             bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D0D0"})

    for ax in (ax1, ax2):
        ax.grid(True, color="#D9D9D9", linewidth=0.6)
        ax.set_xlim(0, t[-1])
        ax.margins(y=0.06)
        ax.tick_params(labelsize=9)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#222222")
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"Title": "I0 analytical capillary imbibition", "Creator": "plot_I0_analytical.py"}
    fig.savefig(args.figure_dir / "I0_analytical.png", dpi=300, metadata=metadata)
    fig.savefig(args.figure_dir / "I0_analytical.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {args.figure_dir / 'I0_analytical.png'}")
    plot_constitutive(args.figure_dir)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--run-dir", type=Path, default=RUN_ROOT)
    p.add_argument("--figure-dir", type=Path, default=FIGURE_ROOT)
    args = p.parse_args()
    main(args)
