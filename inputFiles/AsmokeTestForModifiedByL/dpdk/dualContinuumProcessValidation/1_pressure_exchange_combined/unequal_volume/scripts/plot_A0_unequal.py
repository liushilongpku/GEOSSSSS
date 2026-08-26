#!/usr/bin/env python3
# Purpose: plot unequal-volume A0 pressure histories against discrete and continuous references.
"""Generate the unequal-volume A0 pressure-exchange figure."""

from __future__ import annotations

import os
import tempfile
import argparse
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dual_continuum_validation_mplconfig"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = CASE_ROOT / "runs" / "latest"
FIGURE_ROOT = CASE_ROOT / "figures"


def load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        return np.asarray(h5["pressure Time"])[:, 0], np.asarray(h5["pressure"])[:, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=RUN_ROOT)
    parser.add_argument("--figure-dir", type=Path, default=FIGURE_ROOT)
    args = parser.parse_args()
    t, matrix = load(args.run_dir / "A0_unequal_matrix_pressure_history.hdf5")
    tf, fracture = load(args.run_dir / "A0_unequal_fracture_pressure_history.hdf5")
    if not np.allclose(t, tf):
        raise ValueError("matrix and fracture histories use different time grids")
    vm, vf = 0.95, 0.05
    cm, cf = vm * 0.2 * 1000.0 * 2.0e-12, vf * 1.0 * 1000.0 * 2.0e-12
    pbar = (cm * 2.0e6 + cf * 1.0e6) / (cm + cf)
    delta0 = 1.0e6
    lam = (vm * 1.0e-12 * 1.2e-3 * 1.0e6) * (1.0 / cm + 1.0 / cf)
    q = 1.0 / (1.0 + lam * 0.02)
    tau = 1.0 / lam
    reference_time = t + 0.02
    delta_be = delta0 * q ** (reference_time / 0.02)
    delta_exp = delta0 * np.exp(-reference_time / tau)
    matrix_be = pbar + cf / (cm + cf) * delta_be
    fracture_be = pbar - cm / (cm + cf) * delta_be
    matrix_exp = pbar + cf / (cm + cf) * delta_exp
    fracture_exp = pbar - cm / (cm + cf) * delta_exp
    fig, ax = plt.subplots(figsize=(6.2, 6.0), constrained_layout=True)
    ax.plot(t, matrix / 1.0e6, "o-", ms=3, markevery=10, label="Matrix GEOS", color="#0072B2")
    ax.plot(t, fracture / 1.0e6, "s-", ms=3, markevery=10, label="Fracture GEOS", color="#D55E00")
    ax.plot(t, matrix_be / 1.0e6, "--", label="Matrix backward-Euler reference", color="#0072B2")
    ax.plot(t, fracture_be / 1.0e6, "--", label="Fracture backward-Euler reference", color="#D55E00")
    ax.plot(t, matrix_exp / 1.0e6, ":", label="Matrix continuous reference", color="#0072B2")
    ax.plot(t, fracture_exp / 1.0e6, ":", label="Fracture continuous reference", color="#D55E00")
    ax.set(xlabel="Time (s)", ylabel="Pressure (MPa)", title="A0 unequal-volume pressure exchange")
    ax.grid(True, color="#D9D9D9", linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#222222")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.text(0.02, 0.04, r"$v_m=0.95$, $v_f=0.05$, $\phi_m=0.2$, $\phi_f=1.0$; "
            rf"$p_{{eq}}={pbar / 1.0e6:.6f}$ MPa, $q_{{BE}}={q:.6f}$",
            transform=ax.transAxes, fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D0D0"})
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"Title": "A0 unequal-volume pressure exchange", "Creator": "plot_A0_unequal.py"}
    fig.savefig(args.figure_dir / "A0_unequal_volume_pressure_exchange.png", dpi=300, metadata=metadata)
    fig.savefig(args.figure_dir / "A0_unequal_volume_pressure_exchange.pdf", metadata=metadata)
    plt.close(fig)
    print(f"saved {args.figure_dir / 'A0_unequal_volume_pressure_exchange.png'}")


if __name__ == "__main__":
    main()
