#!/usr/bin/env python3
"""Create the Chinese single-axis DPDP Mandel pressure-validation figure."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "geos_mandel_mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import h5py
import numpy as np
from scipy.interpolate import PchipInterpolator


VALIDATION_DIR = Path(__file__).resolve().parent.parent
RUNS_DIR = VALIDATION_DIR / "runs"
FIM_DIR = RUNS_DIR / "fim_eff_mesh10_20260728_dense"
SEQ_DIR = RUNS_DIR / "seq_eff_mesh10_20260728_dense"
OUTPUT_STEM = VALIDATION_DIR / "analitical_result" / "DPDP_Mandel_压力验证_解析_FIM_迭代耦合"
CHINESE_FONT = Path("/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc")
PM0 = 4.55e5
PF0 = 4.88e5
FIGURE_WIDTH_IN = 8.02 / 2.54
FIGURE_HEIGHT_IN = 5.20 / 2.54


def center_index(element_centers: np.ndarray) -> int:
  """Return the element nearest the no-flow axis at mid-height z=0.015 m."""
  x = element_centers[:, 0]
  z = element_centers[:, 2]
  return int(np.argmin(np.abs(x - x.min()) + np.abs(z - 0.015)))


def load_pressure_history(output_dir: Path, medium: str, t0: float) -> tuple[np.ndarray, np.ndarray]:
  path = output_dir / f"pressure_{medium}_history.hdf5"
  with h5py.File(path, "r") as h5:
    time = np.asarray(h5["pressure Time"])[:, 0]
    pressure = np.asarray(h5["pressure"])
    centers = np.asarray(h5["pressure elementCenter"])[0]

  index = center_index(centers)
  mask = time > 0.0
  return time[mask] / t0, pressure[mask, index]


def densify_log_curve(tau: np.ndarray, value: np.ndarray, points: int = 600) -> tuple[np.ndarray, np.ndarray]:
  """Smooth the analytical curve for display on the logarithmic time axis."""
  log_tau = np.log10(tau)
  dense_log_tau = np.linspace(log_tau.min(), log_tau.max(), points)
  return 10.0**dense_log_tau, PchipInterpolator(log_tau, value)(dense_log_tau)


def log_spaced_marker_indices(tau: np.ndarray, count: int = 24) -> np.ndarray:
  """Select actual history samples with visually uniform log-time spacing."""
  targets = np.linspace(np.log10(tau[0]), np.log10(tau[-1]), count)
  indices = np.searchsorted(np.log10(tau), targets)
  return np.unique(np.clip(indices, 0, len(tau) - 1))


def main() -> None:
  font_manager.fontManager.addfont(str(CHINESE_FONT))
  plt.rcParams.update({
      "font.family": font_manager.FontProperties(fname=CHINESE_FONT).get_name(),
      "font.size": 7,
      "axes.unicode_minus": False,
      "pdf.fonttype": 42,
      "ps.fonttype": 42,
      "axes.linewidth": 0.6,
  })

  import dpdp_mandel_analytical as analytical

  fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), constrained_layout=True)
  media = (
      ("matrix", "基质", "#1B6A8F"),
      ("fracture", "裂缝", "#C84A44"),
  )
  cases = (
      ("全耦合", FIM_DIR, "--", "o"),
      ("迭代耦合", SEQ_DIR, (0, (5, 1.2, 1.2, 1.2)), "s"),
  )
  solution = analytical.solve(tau_min=1e-3, tau_max=3e4, n=25)

  for key, medium_name, color in media:
    analytical_tau, analytical_pressure = densify_log_curve(
        solution["tau"], solution["pm" if key == "matrix" else "pf"]
    )
    ax.semilogx(analytical_tau, analytical_pressure, color=color, linewidth=1.4,
                label=f"{medium_name}：解析解", zorder=2)
    p0 = PM0 if key == "matrix" else PF0
    for method_name, output_dir, linestyle, marker in cases:
      tau, pressure = load_pressure_history(output_dir, key, analytical.t0)
      marker_indices = log_spaced_marker_indices(tau)
      ax.semilogx(tau, pressure / p0, color=color, linestyle=linestyle, linewidth=1.1,
                  marker=marker, markevery=marker_indices, markersize=2.8, markerfacecolor="white",
                  markeredgewidth=0.65, label=f"{medium_name}：{method_name}", zorder=3)

  ax.set_xlim(1e-3, 3e4)
  ax.set_ylim(-0.06, 1.16)
  ax.set_xlabel(r"无量纲时间 $\tau=t/t_0$", fontsize=7.5, labelpad=4)
  ax.set_ylabel(r"无量纲压力 $p_i/p_{i0}^{+}$", fontsize=7.5, labelpad=5)
  ax.tick_params(axis="both", which="major", labelsize=6, length=3.0, width=0.6)
  ax.tick_params(axis="both", which="minor", length=1.5, width=0.45)
  ax.grid(True, which="major", color="#AEB8C2", linewidth=0.28, alpha=0.70)
  ax.grid(True, which="minor", color="#DCE2E7", linewidth=0.18, alpha=0.75)
  legend = ax.legend(ncol=1, loc="lower left", bbox_to_anchor=(0.012, 0.012), frameon=True, framealpha=0.97,
                     edgecolor="#97A4AF", handlelength=2.1, fontsize=5.5, markerscale=1.0, borderaxespad=0.0,
                     columnspacing=1.0, labelspacing=0.35)
  legend.get_frame().set_linewidth(0.35)

  for suffix in (".png", ".pdf"):
    fig.savefig(OUTPUT_STEM.with_suffix(suffix), dpi=600)
  plt.close(fig)


if __name__ == "__main__":
  main()
