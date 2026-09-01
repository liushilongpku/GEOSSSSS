#!/usr/bin/env python3
"""Aggregate the Thomas 1983 dual-continuum single-block recovery and plot it.

Key information: matrix oil-component recovery is the Thomas-compatible result.
Normalized oil-saturation reduction is retained only as a diagnostic. GEOS does
not execute periodic time-history events at maxTime, so the final sample is read
from the VTK output when it is absent from HDF5.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import vtk


ROOT = Path(__file__).resolve().parent.parent
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


def read_history(output_dir: Path, stem: str) -> tuple[np.ndarray, np.ndarray]:
    roots = sorted(output_dir.glob(f"*{stem}.hdf5"))
    if not roots:
        raise FileNotFoundError(f"no {stem} history in {output_dir}")
    path = roots[0]
    with h5py.File(path, "r") as f:
        phase = np.asarray(f["phaseVolumeFraction"])
        time = np.ravel(np.asarray(f["phaseVolumeFraction Time"]))
    if phase.shape[0] != time.size or phase.shape[-1] != 3:
        raise RuntimeError(f"unexpected phase history shape in {path}: {phase.shape}")
    return time, phase.reshape(time.size, -1, 3).mean(axis=1)


def read_amount_history(output_dir: Path, stem: str) -> tuple[np.ndarray, np.ndarray]:
    roots = sorted(output_dir.glob(f"*{stem}.hdf5"))
    if not roots:
        raise FileNotFoundError(f"no {stem} history in {output_dir}")
    path = roots[0]
    with h5py.File(path, "r") as f:
        amount = np.asarray(f["compAmount"])
        time = np.ravel(np.asarray(f["compAmount Time"]))
    if amount.shape[0] != time.size or amount.shape[-1] != 3:
        raise RuntimeError(f"unexpected component history shape in {path}: {amount.shape}")
    return time, amount.reshape(time.size, -1, 3).sum(axis=1)


def read_final_vtk(output_dir: Path) -> tuple[float, np.ndarray, np.ndarray]:
    pvd = output_dir / "vtkFinalOutput.pvd"
    datasets = ET.parse(pvd).getroot().findall("./Collection/DataSet")
    if not datasets:
        raise RuntimeError(f"no VTK timesteps found in {pvd}")
    final = max(datasets, key=lambda node: float(node.attrib["timestep"]))
    vtm = output_dir / final.attrib["file"]
    matrix_files = [
        node.attrib["file"]
        for node in ET.parse(vtm).getroot().findall(".//DataSet")
        if "matrixRegion" in node.attrib.get("file", "")
    ]
    if len(matrix_files) != 1:
        raise RuntimeError(f"expected one matrixRegion dataset in {vtm}, found {matrix_files}")
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtm.parent / matrix_files[0]))
    reader.Update()
    data = reader.GetOutput().GetCellData()
    amount_array = data.GetArray("compAmount")
    phase_array = data.GetArray("phaseVolumeFraction")
    if amount_array is None or phase_array is None:
        raise RuntimeError(f"matrix compAmount or phaseVolumeFraction missing from {matrix_files[0]}")
    amounts = np.asarray(
        [[amount_array.GetComponent(i, j) for j in range(3)] for i in range(amount_array.GetNumberOfTuples())]
    ).sum(axis=0)
    phases = np.asarray(
        [[phase_array.GetComponent(i, j) for j in range(3)] for i in range(phase_array.GetNumberOfTuples())]
    ).mean(axis=0)
    return float(final.attrib["timestep"]), amounts, phases


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
    parser.add_argument("--analysis-dir", type=Path, default=ROOT / "analysis")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    output_root = args.output.resolve()

    amount_time, amounts = read_amount_history(output_root, "thomas_matrix_amount_history")
    phase_time, phases = read_history(output_root, "thomas_matrix_phase_history")
    if not amount_time.size or not phase_time.size:
        raise RuntimeError("no history samples found")
    if amount_time.shape != phase_time.shape or not np.allclose(amount_time, phase_time, rtol=0.0, atol=1.0e-6):
        raise RuntimeError("component and phase histories use different sample times")

    final_time, final_amounts, final_phases = read_final_vtk(output_root)
    if final_time > amount_time[-1] + 1.0e-6:
        amount_time = np.append(amount_time, final_time)
        amounts = np.vstack((amounts, final_amounts))
        phases = np.vstack((phases, final_phases))

    initial_oil_amount = float(amounts[0, 0])
    initial_oil_saturation = float(phases[0, 0])
    reference_time, reference_recovery = reference_curve()

    rows: list[dict[str, float]] = []
    for index, time_seconds in enumerate(amount_time):
        time_years = float(time_seconds) / YEAR
        mass_recovery = 100.0 * (1.0 - float(amounts[index, 0]) / initial_oil_amount)
        saturation_reduction = 100.0 * (1.0 - float(phases[index, 0]) / initial_oil_saturation)
        benchmark = float(np.interp(time_years, reference_time, reference_recovery))
        rows.append(
            {
                "time_years": time_years,
                "matrix_oil_component_mass_recovery_pct": mass_recovery,
                "normalized_oil_saturation_reduction_pct": saturation_reduction,
                "oil_component_mass_minus_thomas_percentage_points": mass_recovery - benchmark,
                "thomas_fig4_recovery_pct": benchmark,
                "matrix_oil_saturation": float(phases[index, 0]),
                "matrix_gas_saturation": float(phases[index, 1]),
                "matrix_water_saturation": float(phases[index, 2]),
            }
        )

    analysis_dir = args.analysis_dir.resolve()
    figures_dir = args.figures_dir.resolve()
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
