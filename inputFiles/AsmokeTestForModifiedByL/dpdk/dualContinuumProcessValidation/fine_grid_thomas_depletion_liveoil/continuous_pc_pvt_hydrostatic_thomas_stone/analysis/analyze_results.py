#!/usr/bin/env python3
"""Analyze the Thomas Eq. (25) deck-approximation fine-grid chain.

Key information: recovery is based on matrix oil-component amount, while the
pressure audit checks the time-varying PVT gas gradient at every fracture cell.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
HDF5_PREFIX = (
    "Problem/domain/MeshBodies/mesh1/meshLevels/Level0/"
    "ElementRegions/elementRegionsGroup"
)
INITIAL_PRESSURE_PA = 38_332_780.5
PSI_TO_PA = 6_894.757293168
PRESSURE_DECLINE_PSI_PER_DAY = 0.75
TOP_DATUM_M = 3.0486096
GRAVITY_M_PER_S2 = 9.81


def values(group: h5py.Group, name: str) -> np.ndarray:
    """Read a GEOS restart wrapper in either direct or grouped form."""

    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[...], dtype=float)
    return np.asarray(obj["__values__"][...], dtype=float)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV that starts with traceability comments."""

    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines))


def write_csv(path: Path, purpose: str, rows: list[dict[str, float | str]]) -> None:
    """Write one traceable analysis table."""

    with path.open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write(f"# Purpose: {purpose}\n")
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def restart_rank_file(output_dir: Path, stem: str, first: bool) -> Path:
    """Resolve the initial or final rank file from a completed segment."""

    roots = sorted(output_dir.glob(f"{stem}_restart_*.root"))
    if not roots:
        raise FileNotFoundError(f"no restart root in {output_dir}")
    root = roots[0] if first else roots[-1]
    rank_file = root.with_suffix("") / "rank_0000000.hdf5"
    if not rank_file.is_file():
        raise FileNotFoundError(rank_file)
    return rank_file


def read_region(path: Path, region: str) -> dict[str, np.ndarray | float]:
    """Aggregate all single-cell subregions belonging to one GEOS region."""

    oil_amount = 0.0
    volume: list[float] = []
    pressure: list[float] = []
    phase_fraction: list[np.ndarray] = []
    centers: list[np.ndarray] = []
    with h5py.File(path, "r") as h5_file:
        subregions = h5_file[f"{HDF5_PREFIX}/{region}/elementSubRegions"]
        for name, subregion in subregions.items():
            if name == "__size__":
                continue
            comp_amount = np.asarray(values(subregion, "compAmount"))
            if comp_amount.ndim == 1:
                oil_amount += float(comp_amount[0])
            else:
                oil_amount += float(np.sum(comp_amount[..., 0]))
            volume.append(float(np.ravel(values(subregion, "elementVolume"))[0]))
            pressure.append(float(np.ravel(values(subregion, "pressure"))[0]))
            phase_fraction.append(np.ravel(values(subregion, "phaseVolumeFraction"))[:3])
            centers.append(np.ravel(values(subregion, "elementCenter"))[:3])

    volume_array = np.asarray(volume)
    phase_array = np.asarray(phase_fraction)
    return {
        "oil_amount": oil_amount,
        "volume": volume_array,
        "pressure": np.asarray(pressure),
        "phase_fraction": phase_array,
        "center": np.asarray(centers),
        "mean_phase_fraction": np.average(phase_array, axis=0, weights=volume_array),
    }


def available_variants(manifest: list[dict[str, str]]) -> list[str]:
    """Return only chains that have every final restart available."""

    variants: list[str] = []
    for variant in ("left", "right"):
        records = [row for row in manifest if row["variant"] == variant]
        if records and all(
            list((ROOT / "runs" / variant / Path(row["deck"]).stem).glob(f"{Path(row['deck']).stem}_restart_*.root"))
            for row in records
        ):
            variants.append(variant)
    if not variants:
        raise RuntimeError("no complete hydrostatic restart chain found")
    return variants


def configure_plotting() -> None:
    """Set stable Chinese scientific-figure defaults."""

    font_path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    font_manager.fontManager.addfont(font_path)
    font_name = font_manager.FontProperties(fname=font_path).get_name()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font_name, "Arial", "DejaVu Sans"],
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "axes.unicode_minus": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def main() -> None:
    """Extract metrics, compare with Thomas, and make validation figures."""

    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    manifest = read_csv(ROOT / "segment_manifest.csv")
    reference = read_csv(ROOT / "reference/thomas_fig4_3d_model.csv")
    variants = available_variants(manifest)
    reference_time = np.asarray([float(row["time_years"]) for row in reference])
    reference_recovery = np.asarray([float(row["oil_recovery_pct"]) for row in reference])

    recovery_rows: list[dict[str, float | str]] = []
    pressure_rows: list[dict[str, float | str]] = []
    profile_rows: list[dict[str, float | str]] = []
    final_regions: dict[str, dict[str, np.ndarray | float]] = {}

    for variant in variants:
        records = [row for row in manifest if row["variant"] == variant]
        first_deck = ROOT / records[0]["deck"]
        first_output = ROOT / "runs" / variant / first_deck.stem
        initial_path = restart_rank_file(first_output, first_deck.stem, first=True)
        initial_matrix = read_region(initial_path, "matrixRegion")
        initial_oil = float(initial_matrix["oil_amount"])
        recovery_rows.append(
            {
                "variant": variant,
                "time_years": 0.0,
                "matrix_oil_component_amount": initial_oil,
                "matrix_oil_recovery_pct": 0.0,
                "matrix_pressure_mean_pa": float(np.average(
                    np.asarray(initial_matrix["pressure"]), weights=np.asarray(initial_matrix["volume"])
                )),
                "matrix_oil_saturation": float(np.asarray(initial_matrix["mean_phase_fraction"])[0]),
                "matrix_gas_saturation": float(np.asarray(initial_matrix["mean_phase_fraction"])[1]),
                "matrix_water_saturation": float(np.asarray(initial_matrix["mean_phase_fraction"])[2]),
                "thomas_fig4_recovery_pct": 0.0,
                "difference_percentage_points": 0.0,
            }
        )

        for record in records:
            deck = ROOT / record["deck"]
            output_dir = ROOT / "runs" / variant / deck.stem
            rank_file = restart_rank_file(output_dir, deck.stem, first=False)
            matrix = read_region(rank_file, "matrixRegion")
            fracture = read_region(rank_file, "fractureRegion")
            time_years = float(record["end_years"])
            recovery = 100.0 * (1.0 - float(matrix["oil_amount"]) / initial_oil)
            benchmark = float(np.interp(time_years, reference_time, reference_recovery))
            mean_phase = np.asarray(matrix["mean_phase_fraction"])
            recovery_rows.append(
                {
                    "variant": variant,
                    "time_years": time_years,
                    "matrix_oil_component_amount": float(matrix["oil_amount"]),
                    "matrix_oil_recovery_pct": recovery,
                    "matrix_pressure_mean_pa": float(np.average(
                        np.asarray(matrix["pressure"]), weights=np.asarray(matrix["volume"])
                    )),
                    "matrix_oil_saturation": float(mean_phase[0]),
                    "matrix_gas_saturation": float(mean_phase[1]),
                    "matrix_water_saturation": float(mean_phase[2]),
                    "thomas_fig4_recovery_pct": benchmark,
                    "difference_percentage_points": recovery - benchmark,
                }
            )

            center_z = np.asarray(fracture["center"])[:, 2]
            fracture_pressure = np.asarray(fracture["pressure"])
            datum_pressure = INITIAL_PRESSURE_PA - (
                PRESSURE_DECLINE_PSI_PER_DAY * PSI_TO_PA * 365.0 * time_years
            )
            gas_density = float(record["gas_density_end_kg_per_m3"])
            expected_dp_dz = -gas_density * GRAVITY_M_PER_S2
            expected_pressure = datum_pressure + gas_density * GRAVITY_M_PER_S2 * (
                TOP_DATUM_M - center_z
            )
            fitted_slope, fitted_intercept = np.polyfit(center_z, fracture_pressure, 1)
            error = fracture_pressure - expected_pressure
            pressure_rows.append(
                {
                    "variant": variant,
                    "time_years": time_years,
                    "datum_pressure_pa": datum_pressure,
                    "gas_density_kg_per_m3": gas_density,
                    "expected_dp_dz_pa_per_m": expected_dp_dz,
                    "fitted_dp_dz_pa_per_m": float(fitted_slope),
                    "slope_error_pa_per_m": float(fitted_slope - expected_dp_dz),
                    "fitted_intercept_pa": float(fitted_intercept),
                    "max_abs_pressure_error_pa": float(np.max(np.abs(error))),
                    "rms_pressure_error_pa": float(np.sqrt(np.mean(error**2))),
                    "pressure_spread_pa": float(np.ptp(fracture_pressure)),
                }
            )
            final_regions[variant] = matrix

    for variant, matrix in final_regions.items():
        centers = np.asarray(matrix["center"])
        volume = np.asarray(matrix["volume"])
        phase = np.asarray(matrix["phase_fraction"])
        for z_value in sorted(set(np.round(centers[:, 2], decimals=9))):
            mask = np.isclose(centers[:, 2], z_value, atol=1.0e-8)
            averaged = np.average(phase[mask], axis=0, weights=volume[mask])
            profile_rows.append(
                {
                    "variant": variant,
                    "height_m": float(z_value),
                    "oil_saturation": float(averaged[0]),
                    "gas_saturation": float(averaged[1]),
                    "water_saturation": float(averaged[2]),
                }
            )

    write_csv(
        OUT / "recovery_results.csv",
        "Compare matrix oil-component recovery with the digitized Thomas 1983 Fig. 4 curve.",
        recovery_rows,
    )
    write_csv(
        OUT / "fracture_hydrostatic_diagnostics.csv",
        "Verify the imposed fixed-density fracture-gas pressure gradient and absolute pressure.",
        pressure_rows,
    )
    write_csv(
        OUT / "final_vertical_profile.csv",
        "Record volume-weighted matrix phase saturation by height at 2.5 years.",
        profile_rows,
    )

    configure_plotting()
    colors = {"left": "#0072B2", "right": "#D55E00"}
    labels = {"left": "左端点 Pc", "right": "右端点 Pc"}

    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=320, layout="constrained")
    ax.fill_between(
        reference_time,
        reference_recovery - 1.0,
        reference_recovery + 1.0,
        color="#BDBDBD",
        alpha=0.35,
        label="Thomas 图像读取不确定度（约 ±1 个百分点）",
    )
    ax.plot(reference_time, reference_recovery, color="#111111", linewidth=2.0, label="Thomas Fig. 4：3D model")
    if set(variants) == {"left", "right"}:
        left_rows = [row for row in recovery_rows if row["variant"] == "left"]
        right_rows = [row for row in recovery_rows if row["variant"] == "right"]
        envelope_time = np.asarray([float(row["time_years"]) for row in left_rows])
        left_recovery = np.asarray([float(row["matrix_oil_recovery_pct"]) for row in left_rows])
        right_recovery = np.asarray([float(row["matrix_oil_recovery_pct"]) for row in right_rows])
        ax.fill_between(
            envelope_time,
            np.minimum(left_recovery, right_recovery),
            np.maximum(left_recovery, right_recovery),
            color="#56B4E9",
            alpha=0.22,
            label="GEOS Pc 分段端点包络",
        )
    for variant in variants:
        rows = [row for row in recovery_rows if row["variant"] == variant]
        ax.plot(
            [float(row["time_years"]) for row in rows],
            [float(row["matrix_oil_recovery_pct"]) for row in rows],
            color=colors[variant],
            marker="o",
            markersize=3.2,
            linewidth=1.5,
            label=labels[variant],
        )
    ax.set(xlabel="时间（年）", ylabel="基质原油组分采收率（%）", title="Thomas 10 ft 细网格：PVT 动态静水裂缝压力")
    ax.set_xlim(0.0, 2.55)
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.legend(loc="lower right", frameon=True, framealpha=0.95)
    for spine in ax.spines.values():
        spine.set_visible(True)
    fig.savefig(FIGURES / "recovery_comparison.png", facecolor="white")
    fig.savefig(FIGURES / "recovery_comparison.svg", facecolor="white")
    fig.savefig(FIGURES / "recovery_comparison.pdf", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), dpi=320, layout="constrained")
    for variant in variants:
        rows = [row for row in pressure_rows if row["variant"] == variant]
        times = [float(row["time_years"]) for row in rows]
        axes[0].plot(
            times,
            [float(row["fitted_dp_dz_pa_per_m"]) for row in rows],
            color=colors[variant],
            marker="o",
            markersize=3,
            label=labels[variant],
        )
        axes[0].plot(
            times,
            [float(row["expected_dp_dz_pa_per_m"]) for row in rows],
            color="#111111",
            linestyle="--",
            linewidth=1.2,
            label="PVT 设定斜率" if variant == variants[0] else None,
        )
        axes[1].plot(
            times,
            [float(row["max_abs_pressure_error_pa"]) for row in rows],
            color=colors[variant],
            marker="o",
            markersize=3,
            label=labels[variant],
        )
    axes[0].set(xlabel="时间（年）", ylabel=r"$dp_f/dz$（Pa/m）", title="裂缝气体静水压力梯度")
    axes[1].set(xlabel="时间（年）", ylabel="最大绝对误差（Pa）", title="裂缝压力对设定函数的误差")
    for axis in axes:
        axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.7)
        axis.legend(frameon=True)
        for spine in axis.spines.values():
            spine.set_visible(True)
    fig.savefig(FIGURES / "fracture_hydrostatic_check.png", facecolor="white")
    fig.savefig(FIGURES / "fracture_hydrostatic_check.svg", facecolor="white")
    fig.savefig(FIGURES / "fracture_hydrostatic_check.pdf", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), dpi=320, layout="constrained")
    for variant in variants:
        rows = [row for row in profile_rows if row["variant"] == variant]
        height = [float(row["height_m"]) for row in rows]
        axes[0].plot(
            [float(row["oil_saturation"]) for row in rows],
            height,
            color=colors[variant],
            marker="o",
            label=labels[variant],
        )
        axes[1].plot(
            [float(row["gas_saturation"]) for row in rows],
            height,
            color=colors[variant],
            marker="o",
            label=labels[variant],
        )
    axes[0].set(xlabel="油相饱和度", ylabel="高度（m）", title="2.5 年基质油相饱和度")
    axes[1].set(xlabel="气相饱和度", ylabel="高度（m）", title="2.5 年基质气相饱和度")
    for axis in axes:
        axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.7)
        axis.legend(frameon=True)
        for spine in axis.spines.values():
            spine.set_visible(True)
    fig.savefig(FIGURES / "final_vertical_saturation.png", facecolor="white")
    fig.savefig(FIGURES / "final_vertical_saturation.svg", facecolor="white")
    fig.savefig(FIGURES / "final_vertical_saturation.pdf", facecolor="white")
    plt.close(fig)

    for variant in variants:
        selected = [row for row in recovery_rows if row["variant"] == variant]
        checkpoints = {round(float(row["time_years"]), 1): row for row in selected}
        print(f"{variant}:")
        for time_years in (0.1, 0.5, 1.0, 2.5):
            row = checkpoints[time_years]
            print(
                f"  {time_years:.1f} yr: GEOS={float(row['matrix_oil_recovery_pct']):.3f}%, "
                f"Thomas={float(row['thomas_fig4_recovery_pct']):.3f}%, "
                f"difference={float(row['difference_percentage_points']):+.3f} pp"
            )
        pressure = [row for row in pressure_rows if row["variant"] == variant][-1]
        print(
            f"  final dp/dz={float(pressure['fitted_dp_dz_pa_per_m']):.6f} Pa/m, "
            f"max pressure error={float(pressure['max_abs_pressure_error_pa']):.6e} Pa"
        )


if __name__ == "__main__":
    main()
