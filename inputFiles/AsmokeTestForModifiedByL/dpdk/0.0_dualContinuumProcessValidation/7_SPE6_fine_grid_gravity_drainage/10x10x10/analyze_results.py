#!/usr/bin/env python3
"""Aggregate SPE6 10x10x10 fine-grid recovery and create comparison plots.

Key information: only the 512 interior matrix cells are counted; saturation-
  inventory and oil-component-mass recovery are reported separately. The run
  contains VTK output only; no restart or HDF5 output is required.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path("/tmp/fine_grid_spe6_10x10x10_vtk_only_unlimited")
COARSE_RESULTS = ROOT.parent / "7x7x8/analysis/recovery_results.csv"
YEAR = 31_536_000.0


def vtk_checkpoint_files(output_dir: Path) -> list[Path]:
    """Resolve ordered matrix-region VTK files from one VTK-only run."""

    files = list(output_dir.glob("vtkOutput/*/mesh1/Level0/matrixRegion/rank_0.vtu"))
    return sorted(files, key=lambda path: int(path.parents[3].name))


def read_result_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read time and volume-inventory recovery from an analysis CSV file."""

    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")
    ]
    rows = list(csv.DictReader(lines))
    return (
        np.asarray([float(row["time_years"]) for row in rows]),
        np.asarray([float(row["oil_saturation_inventory_recovery_pct"]) for row in rows]),
    )


def matrix_state(path: Path) -> dict[str, float | np.ndarray]:
    """Aggregate oil amount and volume-weighted phase fractions from one VTK file."""

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    cell_data = grid.GetCellData()

    def values(name: str) -> np.ndarray:
        array = cell_data.GetArray(name)
        if array is None:
            raise KeyError(f"VTK cell array not found: {name}")
        return np.asarray(vtk_to_numpy(array), dtype=float)

    volume = values("elementVolume")
    porosity = values("matrixPorosity_porosity")
    phase = values("phaseVolumeFraction")
    comp_amount = values("compAmount")
    pressure = values("pressure")
    centers = values("elementCenter")
    time_array = grid.GetFieldData().GetArray("TIME")
    if time_array is None:
        raise KeyError(f"VTK TIME field not found: {path}")
    time_seconds = float(vtk_to_numpy(time_array)[0])
    pore_volume = volume * porosity
    return {
        "cell_count": len(volume),
        "time_years": time_seconds / YEAR,
        "oil_amount": float(np.sum(comp_amount[:, 0])),
        "oil_inventory": float(np.sum(pore_volume * phase[:, 0])),
        "mean_phase": np.average(phase, axis=0, weights=pore_volume),
        "mean_pressure_pa": float(np.average(pressure, weights=pore_volume)),
        "volume": volume,
        "pore_volume": pore_volume,
        "phase": phase,
        "center": np.asarray(centers),
    }


def reference_curve() -> tuple[np.ndarray, np.ndarray]:
    """Read the digitized SPE6 zero-fracture-Pc engineering reference."""

    lines = [
        line
        for line in (ROOT / "reference/spe_fig1_gas_oil_zeroPcf.csv").read_text(
            encoding="utf-8"
        ).splitlines()
        if not line.startswith("#")
    ]
    rows = list(csv.DictReader(lines))
    return (
        np.asarray([float(row["time_years"]) for row in rows]),
        100.0 * np.asarray([float(row["recovery"]) for row in rows]),
    )


def configure_plotting() -> None:
    """Configure stable Chinese scientific plotting defaults."""

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
    """Write recovery/profile CSV files and the SPE6 comparison figure."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--artifact-dir", type=Path, default=ROOT)
    args = parser.parse_args()
    output_root = args.output.resolve()
    vtk_files = vtk_checkpoint_files(output_root)
    if len(vtk_files) < 2:
        raise RuntimeError(f"expected initial and final VTK files, found {len(vtk_files)}")
    states = [matrix_state(path) for path in vtk_files]
    initial = states[0]
    if int(initial["cell_count"]) != 512:
        raise RuntimeError(f"expected 512 matrix cells, found {initial['cell_count']}")
    initial_oil_amount = float(initial["oil_amount"])
    initial_oil_inventory = float(initial["oil_inventory"])
    reference_time, reference_recovery = reference_curve()

    rows: list[dict[str, float]] = []
    rows.append(
        {
            "time_years": 0.0,
            "oil_saturation_inventory_recovery_pct": 0.0,
            "oil_component_mass_recovery_pct": 0.0,
            "spe6_digitized_recovery_pct": 0.0,
            "volume_minus_spe6_percentage_points": 0.0,
            "matrix_oil_saturation": float(np.asarray(initial["mean_phase"])[0]),
            "matrix_gas_saturation": float(np.asarray(initial["mean_phase"])[1]),
            "matrix_water_saturation": float(np.asarray(initial["mean_phase"])[2]),
            "matrix_mean_pressure_pa": float(initial["mean_pressure_pa"]),
        }
    )
    for state in states[1:]:
        time_years = float(state["time_years"])
        volume_recovery = 100.0 * (1.0 - float(state["oil_inventory"]) / initial_oil_inventory)
        mass_recovery = 100.0 * (1.0 - float(state["oil_amount"]) / initial_oil_amount)
        benchmark = float(np.interp(time_years, reference_time, reference_recovery))
        mean_phase = np.asarray(state["mean_phase"])
        rows.append(
            {
                "time_years": time_years,
                "oil_saturation_inventory_recovery_pct": volume_recovery,
                "oil_component_mass_recovery_pct": mass_recovery,
                "spe6_digitized_recovery_pct": benchmark,
                "volume_minus_spe6_percentage_points": volume_recovery - benchmark,
                "matrix_oil_saturation": float(mean_phase[0]),
                "matrix_gas_saturation": float(mean_phase[1]),
                "matrix_water_saturation": float(mean_phase[2]),
                "matrix_mean_pressure_pa": float(state["mean_pressure_pa"]),
            }
        )

    artifact_root = args.artifact_dir.resolve()
    analysis_dir = artifact_root / "analysis"
    figures_dir = artifact_root / "figures"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    with (analysis_dir / "recovery_results.csv").open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write("# Purpose: compare 512-cell matrix recovery with the digitized SPE6 Fig. 1 reference.\n")
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    final = states[-1]
    centers = np.asarray(final["center"])
    pore_volume = np.asarray(final["pore_volume"])
    phase = np.asarray(final["phase"])
    profile_rows: list[dict[str, float]] = []
    for z_value in sorted(set(np.round(centers[:, 2], decimals=9))):
        mask = np.isclose(centers[:, 2], z_value, atol=1.0e-8)
        averaged = np.average(phase[mask], axis=0, weights=pore_volume[mask])
        profile_rows.append(
            {
                "height_m": float(z_value),
                "oil_saturation": float(averaged[0]),
                "gas_saturation": float(averaged[1]),
                "water_saturation": float(averaged[2]),
            }
        )
    with (analysis_dir / "final_vertical_profile.csv").open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write(
        f"# Purpose: record the {float(states[-1]['time_years']):g}-year "
            "volume-weighted phase saturations by matrix layer.\n"
        )
        writer = csv.DictWriter(csv_file, fieldnames=profile_rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(profile_rows)

    configure_plotting()
    time = np.asarray([row["time_years"] for row in rows])
    volume_recovery = np.asarray([row["oil_saturation_inventory_recovery_pct"] for row in rows])
    mass_recovery = np.asarray([row["oil_component_mass_recovery_pct"] for row in rows])
    final_years = float(time[-1])
    fig, ax = plt.subplots(figsize=(7.4, 5.5), dpi=260, layout="constrained")
    ax.plot(
        reference_time,
        reference_recovery,
        color="#222222",
        linestyle=(0, (5, 2)),
        linewidth=2.1,
        marker="s",
        markersize=4.3,
        label=f"SPE6 Fig. 1 数字化参考（5 年 {reference_recovery[-1]:.1f}%）",
    )
    ax.plot(
        time,
        volume_recovery,
        color="#0072B2",
        linewidth=2.2,
        marker="o",
        markersize=3.8,
        label=f"GEOS 油相饱和度库存（{final_years:g} 年 {volume_recovery[-1]:.2f}%）",
    )
    ax.plot(
        time,
        mass_recovery,
        color="#D55E00",
        linewidth=2.0,
        marker="^",
        markersize=3.6,
        label=f"GEOS 油组分质量（{final_years:g} 年 {mass_recovery[-1]:.2f}%）",
    )
    ax.set_xlim(0.0, 5.05)
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel("基质原油采收率 / %")
    ax.set_title("SPE6 单块物理条件下的 10×10×10 细网格试验")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.96)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    fig.savefig(figures_dir / "recovery_comparison.png", facecolor="white")
    plt.close(fig)

    if COARSE_RESULTS.is_file():
        coarse_time, coarse_recovery = read_result_columns(COARSE_RESULTS)
        fig, ax = plt.subplots(figsize=(7.4, 5.5), dpi=260, layout="constrained")
        ax.plot(
            reference_time,
            reference_recovery,
            color="#222222",
            linestyle=(0, (5, 2)),
            linewidth=2.2,
            marker="s",
            markersize=4.0,
            label=f"SPE6 Fig. 1 参考（5 年 {reference_recovery[-1]:.1f}%）",
        )
        ax.plot(
            coarse_time,
            coarse_recovery,
            color="#D55E00",
            linewidth=2.0,
            marker="^",
            markersize=3.8,
            label=f"7×7×8（5 年 {coarse_recovery[-1]:.2f}%）",
        )
        ax.plot(
            time,
            volume_recovery,
            color="#0072B2",
            linewidth=2.2,
            marker="o",
            markersize=3.8,
            label=f"10×10×10（5 年 {volume_recovery[-1]:.2f}%）",
        )
        ax.set_xlim(0.0, 5.05)
        ax.set_xlabel("时间 / 年")
        ax.set_ylabel("基质油相体积采收率 / %")
        ax.set_title("细网格分辨率敏感性")
        ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
        ax.legend(loc="lower right", frameon=True, framealpha=0.96)
        for spine in ax.spines.values():
            spine.set_linewidth(0.9)
        fig.savefig(figures_dir / "grid_sensitivity_comparison.png", facecolor="white")
        plt.close(fig)

    profile_height = np.asarray([row["height_m"] for row in profile_rows])
    profile_oil = np.asarray([row["oil_saturation"] for row in profile_rows])
    profile_gas = np.asarray([row["gas_saturation"] for row in profile_rows])
    profile_water = np.asarray([row["water_saturation"] for row in profile_rows])
    fig, ax = plt.subplots(figsize=(6.0, 5.5), dpi=260, layout="constrained")
    ax.plot(profile_oil, profile_height, color="#0072B2", marker="o", linewidth=2.2, label="油相")
    ax.plot(profile_gas, profile_height, color="#D55E00", marker="^", linewidth=2.2, label="气相")
    ax.plot(
        profile_water,
        profile_height,
        color="#009E73",
        linestyle=(0, (5, 2)),
        marker="s",
        linewidth=2.0,
        label="水相",
    )
    ax.set_xlim(0.0, 0.85)
    ax.set_ylim(0.0, 3.05)
    ax.set_xlabel("相饱和度")
    ax.set_ylabel("矩阵单元中心高度 / m")
    ax.set_title(f"{final_years:g} 年矩阵垂向饱和度剖面")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="center right", frameon=True, framealpha=0.96)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    fig.savefig(figures_dir / "final_vertical_profile.png", facecolor="white")
    plt.close(fig)
    print(
        f"{final_years:g}-year recovery: saturation inventory={volume_recovery[-1]:.6f}%, "
        f"oil-component mass={mass_recovery[-1]:.6f}%"
    )


if __name__ == "__main__":
    main()
