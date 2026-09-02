#!/usr/bin/env python3
"""Plot the archived Thomas fine-grid recovery against the Fig. 4 reference."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT.parent / "single_block/reference/thomas_fig4_3d_model.csv"


def read_csv(path: Path, value_column: str) -> tuple[np.ndarray, np.ndarray]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    rows = list(csv.DictReader(lines))
    return (
        np.asarray([float(row["time_years"]) for row in rows]),
        np.asarray([float(row[value_column]) for row in rows]),
    )


def configure_plotting() -> None:
    font_path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if font_path.is_file():
        font_manager.fontManager.addfont(font_path)
        family = font_manager.FontProperties(fname=font_path).get_name()
    else:
        family = "DejaVu Sans"
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [family, "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 11,
        }
    )


def main() -> None:
    geos_time, geos_recovery = read_csv(
        ROOT / "analysis/recovery_results.csv", "oil_component_mass_recovery_pct"
    )
    reference_time, reference_recovery = read_csv(REFERENCE, "oil_recovery_pct")

    configure_plotting()
    figure, axis = plt.subplots(figsize=(7.4, 5.5), dpi=260, layout="constrained")
    axis.fill_between(
        reference_time,
        reference_recovery - 1.0,
        reference_recovery + 1.0,
        color="#BDBDBD",
        alpha=0.35,
        label="Thomas 图像读取不确定度（约 ±1 个百分点）",
    )
    axis.plot(
        reference_time,
        reference_recovery,
        color="#222222",
        linestyle=(0, (5, 2)),
        linewidth=2.1,
        marker="s",
        markersize=4.3,
        label=f"Thomas Fig. 4（2.5 年 {np.interp(2.5, reference_time, reference_recovery):.1f}%）",
    )
    axis.plot(
        geos_time,
        geos_recovery,
        color="#D55E00",
        linewidth=2.2,
        marker="o",
        markersize=4.2,
        label=f"GEOS 7×7×8 细网格（2.5 年 {geos_recovery[-1]:.2f}%）",
    )
    axis.set_xlim(0.0, 2.55)
    axis.set_ylim(0.0, 50.0)
    axis.set_xlabel("时间 / 年")
    axis.set_ylabel("基质原油质量采收率 / %")
    axis.set_title("Thomas 10 ft 细网格：气油重力排驱")
    axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    axis.legend(loc="lower right", frameon=True, framealpha=0.96)
    for spine in axis.spines.values():
        spine.set_linewidth(0.9)

    output = ROOT / "figures/recovery_comparison.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
