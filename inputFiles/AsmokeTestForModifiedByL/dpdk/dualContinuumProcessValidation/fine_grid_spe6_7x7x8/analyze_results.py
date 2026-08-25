#!/usr/bin/env python3
"""Aggregate SPE6 fine-grid matrix recovery and create the comparison plot.

Key information: only the 150 interior matrix cells are counted; saturation-
inventory and oil-component-mass recovery are reported separately, and partial
restart chains can be analyzed without overwriting final artifacts.
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


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path("/tmp/fine_grid_spe6_7x7x8_saturated_tol5e2_max12")
HDF5_PREFIX = (
    "Problem/domain/MeshBodies/mesh1/meshLevels/Level0/"
    "ElementRegions/elementRegionsGroup/matrixRegion/elementSubRegions"
)


def wrapper_values(group: h5py.Group, name: str) -> np.ndarray:
    """Read a restart wrapper in direct or grouped form."""

    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[...], dtype=float)
    return np.asarray(obj["__values__"][...], dtype=float)


def read_manifest() -> list[dict[str, str]]:
    """Read the generated segment manifest."""

    lines = [
        line
        for line in (ROOT / "segment_manifest.csv").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    return list(csv.DictReader(lines))


def restart_rank_file(output_dir: Path, stem: str, first: bool) -> Path:
    """Resolve the first or final rank file in one completed segment."""

    roots = sorted(output_dir.glob(f"{stem}_restart_*.root"))
    if not roots:
        raise FileNotFoundError(f"no restart roots in {output_dir}")
    root = roots[0] if first else roots[-1]
    rank_file = root.with_suffix("") / "rank_0000000.hdf5"
    if not rank_file.is_file():
        raise FileNotFoundError(rank_file)
    return rank_file


def matrix_state(path: Path) -> dict[str, float | np.ndarray]:
    """Aggregate oil amount and volume-weighted phase fractions over matrix cells."""

    oil_amount = 0.0
    volumes: list[float] = []
    porosities: list[float] = []
    phases: list[np.ndarray] = []
    pressures: list[float] = []
    centers: list[np.ndarray] = []
    with h5py.File(path, "r") as h5_file:
        subregions = h5_file[HDF5_PREFIX]
        for name, subregion in subregions.items():
            if name == "__size__":
                continue
            comp_amount = np.asarray(wrapper_values(subregion, "compAmount"))
            oil_amount += float(np.sum(comp_amount[..., 0]))
            volumes.append(float(np.ravel(wrapper_values(subregion, "elementVolume"))[0]))
            porosities.append(
                float(np.ravel(wrapper_values(subregion["matrixPorosity"], "porosity"))[0])
            )
            phases.append(np.ravel(wrapper_values(subregion, "phaseVolumeFraction"))[:3])
            pressures.append(float(np.ravel(wrapper_values(subregion, "pressure"))[0]))
            centers.append(np.ravel(wrapper_values(subregion, "elementCenter"))[:3])
    volume = np.asarray(volumes)
    porosity = np.asarray(porosities)
    pore_volume = volume * porosity
    phase = np.asarray(phases)
    pressure = np.asarray(pressures)
    return {
        "cell_count": len(volumes),
        "oil_amount": oil_amount,
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
    parser.add_argument("--segments", type=int, default=None)
    parser.add_argument("--artifact-dir", type=Path, default=ROOT)
    args = parser.parse_args()
    output_root = args.output.resolve()
    records = read_manifest()
    if args.segments is not None:
        if not 1 <= args.segments <= len(records):
            raise ValueError(f"segments must be in [1, {len(records)}]")
        records = records[: args.segments]
    first_deck = ROOT / records[0]["deck"]
    first_output = output_root / first_deck.stem
    initial = matrix_state(restart_rank_file(first_output, first_deck.stem, first=True))
    if int(initial["cell_count"]) != 150:
        raise RuntimeError(f"expected 150 matrix cells, found {initial['cell_count']}")
    initial_oil_amount = float(initial["oil_amount"])
    initial_oil_inventory = float(initial["oil_inventory"])
    reference_time, reference_recovery = reference_curve()

    rows: list[dict[str, float]] = []
    states: list[dict[str, float | np.ndarray]] = [initial]
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
    for record in records:
        deck = ROOT / record["deck"]
        output_dir = output_root / deck.stem
        state = matrix_state(restart_rank_file(output_dir, deck.stem, first=False))
        states.append(state)
        time_years = float(record["end_years"])
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
        csv_file.write("# Purpose: compare 150-cell matrix recovery with the digitized SPE6 Fig. 1 reference.\n")
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
            f"# Purpose: record the {float(records[-1]['end_years']):g}-year "
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
    ax.set_title("SPE6 单块物理条件下的 7×7×8 细网格近似")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.96)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    fig.savefig(figures_dir / "recovery_comparison.png", facecolor="white")
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
