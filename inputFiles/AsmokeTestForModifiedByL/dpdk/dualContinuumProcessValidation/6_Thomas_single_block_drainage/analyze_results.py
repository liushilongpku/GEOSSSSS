#!/usr/bin/env python3
"""Aggregate the Thomas 1983 dual-continuum single-block recovery and plot it.

Key information: the matrix oil-component amount history is read from the GEOS
TimeHistory output and compared against the digitized Thomas 1983 Fig. 4 curve.
Both saturation-inventory and oil-component-mass recovery are reported, because
SPE6/Thomas text only states an oil-recovery percentage without a normalization
formula.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


ROOT = Path(__file__).resolve().parents[0]
DEFAULT_OUTPUT = Path("/tmp/thomas_single_block_dual_continuum")
YEAR = 31_536_000.0


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


def read_history(output_dir: Path, stem: str) -> list[dict[str, float]]:
    roots = sorted(output_dir.glob(f"*{stem}.hdf5"))
    if not roots:
        raise FileNotFoundError(f"no {stem} history in {output_dir}")
    path = roots[0]
    rows: list[dict[str, float]] = []
    with h5py.File(path, "r") as f:
        phase = f["phaseVolumeFraction"]
        time = np.ravel(f["phaseVolumeFraction Time"])
        for index, t in enumerate(time):
            values = np.ravel(phase[index])
            rows.append(
                {
                    "time_seconds": float(t),
                    "oil_saturation": float(values[0]),
                    "gas_saturation": float(values[1]),
                    "water_saturation": float(values[2]),
                }
            )
    return rows


def read_amount_history(output_dir: Path, stem: str) -> list[dict[str, float]]:
    roots = sorted(output_dir.glob(f"*{stem}.hdf5"))
    if not roots:
        raise FileNotFoundError(f"no {stem} history in {output_dir}")
    path = roots[0]
    rows: list[dict[str, float]] = []
    with h5py.File(path, "r") as f:
        amount = f["compAmount"]
        time = np.ravel(f["compAmount Time"])
        for i, t in enumerate(time):
            vals = np.asarray(amount[i])
            rows.append(
                {"time_seconds": float(t), "oil_amount": float(np.ravel(vals)[0])}
            )
    return rows


def reference_curve() -> tuple[np.ndarray, np.ndarray]:
    lines = [
        line
        for line in (ROOT / "reference/thomas_fig4_3d_model.csv").read_text(
            encoding="utf-8"
        ).splitlines()
        if not line.startswith("#")
    ]
    rows = list(csv.DictReader(lines))
    return (
        np.asarray([float(row["time_years"]) for row in rows]),
        np.asarray([float(row["oil_recovery_pct"]) for row in rows]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_root = args.output.resolve()

    amount_rows = read_amount_history(output_root, "thomas_matrix_amount_history")
    phase_rows = read_history(output_root, "thomas_matrix_phase_history")
    if not amount_rows or not phase_rows:
        raise RuntimeError("no history samples found")

    initial_oil_amount = float(amount_rows[0]["oil_amount"])
    initial_phase = phase_rows[0]
    reference_time, reference_recovery = reference_curve()

    rows: list[dict[str, float]] = []
    for index, (amount, phase) in enumerate(zip(amount_rows, phase_rows)):
        time_years = float(amount["time_seconds"]) / YEAR
        mass_recovery = 100.0 * (1.0 - float(amount["oil_amount"]) / initial_oil_amount)
        benchmark = float(np.interp(time_years, reference_time, reference_recovery))
        rows.append(
            {
                "time_years": time_years,
                "matrix_oil_component_mass_recovery_pct": mass_recovery,
                "oil_component_mass_minus_thomas_percentage_points": mass_recovery - benchmark,
                "thomas_fig4_recovery_pct": benchmark,
                "matrix_oil_saturation": float(phase["oil_saturation"]),
                "matrix_gas_saturation": float(phase["gas_saturation"]),
                "matrix_water_saturation": float(phase["water_saturation"]),
            }
        )

    analysis_dir = ROOT / "analysis"
    figures_dir = ROOT / "figures"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    with (analysis_dir / "recovery_results.csv").open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write(
            "# Purpose: compare dual-continuum single-block matrix oil recovery with the digitized Thomas 1983 Fig. 4 curve.\n"
        )
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    configure_plotting()
    times = np.asarray([row["time_years"] for row in rows])
    mass = np.asarray([row["matrix_oil_component_mass_recovery_pct"] for row in rows])
    fig, ax = plt.subplots(figsize=(7.4, 5.5), dpi=260, layout="constrained")
    ax.fill_between(
        reference_time,
        reference_recovery - 1.0,
        reference_recovery + 1.0,
        color="#BDBDBD",
        alpha=0.35,
        label="Thomas 图像读取不确定度（约 ±1 个百分点）",
    )
    ax.plot(
        reference_time,
        reference_recovery,
        color="#222222",
        linestyle=(0, (5, 2)),
        linewidth=2.1,
        marker="s",
        markersize=4.3,
        label=f"Thomas Fig. 4（2.5 年 {reference_recovery[-1]:.1f}%）",
    )
    ax.plot(
        times,
        mass,
        color="#0072B2",
        linewidth=2.2,
        marker="o",
        markersize=3.8,
        label=f"GEOS 油组分质量（{times[-1]:.1f} 年 {mass[-1]:.2f}%）",
    )
    ax.set_xlim(0.0, 2.55)
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel("基质原油采收率 / %")
    ax.set_title("Thomas 10 ft 双重介质单块：气油重力排驱")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.96)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    fig.savefig(figures_dir / "recovery_comparison.png", facecolor="white")
    fig.savefig(figures_dir / "recovery_comparison.pdf", facecolor="white")
    plt.close(fig)

    for row in rows:
        if abs(row["time_years"] - 0.5) < 1.0e-6 or abs(row["time_years"] - 2.5) < 1.0e-6:
            print(
                f"{row['time_years']:.1f} yr: GEOS={row['matrix_oil_component_mass_recovery_pct']:.4f}%, "
                f"Thomas={row['thomas_fig4_recovery_pct']:.4f}%, "
                f"difference={row['oil_component_mass_minus_thomas_percentage_points']:+.4f} pp, "
                f"So={row['matrix_oil_saturation']:.5f}, Sg={row['matrix_gas_saturation']:.5f}, "
                f"Sw={row['matrix_water_saturation']:.5f}"
            )


if __name__ == "__main__":
    main()
