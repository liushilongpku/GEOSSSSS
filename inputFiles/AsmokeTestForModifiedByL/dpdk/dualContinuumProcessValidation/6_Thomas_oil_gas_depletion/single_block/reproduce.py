#!/usr/bin/env python3
"""Check or reproduce the canonical Thomas single-block validation results."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
EXPECTED_RECOVERY = {
    "gdp_off": {0.5: 13.2473, 2.5: 18.7311},
    "gdp_on": {0.5: 27.3123, 2.5: 45.5305},
}
DECKS = {
    "gdp_off": ROOT / "thomas_singleblock_gas_oil_gravity_drainage.xml",
    "gdp_on": ROOT / "thomas_singleblock_gas_oil_gravity_drainage_gdp_on.xml",
}
INPUT_FILES = [
    ROOT / "tables/pvto_bo.txt",
    ROOT / "tables/pvtg_norv_bo.txt",
    ROOT / "tables/pvtw_bo.txt",
    ROOT / "tables/ve_pseudo_pc_sg_axis.txt",
    ROOT / "tables/ve_pseudo_pc_pressure_axis.txt",
    ROOT / "tables/ve_pseudo_pc_values.txt",
    ROOT / "reference/thomas_fig4_3d_model.csv",
]
SOURCE_FILES = [
    "src/coreComponents/constitutive/relativePermeability/ThreePhaseTableRelativePermeability.cpp",
    "src/coreComponents/constitutive/gravityDrainagePressure/ThomasGasOilGravityDrainagePressure.cpp",
    "src/coreComponents/constitutive/capillaryPressure/PressureScaledTableCapillaryPressure.cpp",
    "src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumCrossFlow.cpp",
    "src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/PPUPhaseFlux.hpp",
]


def find_repository() -> Path:
    for parent in ROOT.parents:
        if (parent / "src/main/main.cpp").is_file():
            return parent
    raise RuntimeError("could not locate the GEOS repository above this validation directory")


def run(command: list[str], cwd: Path = ROOT) -> None:
    print("+ " + " ".join(command), flush=True)
    environment = os.environ | {"PYTHONDONTWRITEBYTECODE": "1"}
    subprocess.run(command, cwd=cwd, env=environment, check=True)


def require_files(paths: list[Path], description: str) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing {description}:\n" + "\n".join(missing))


def check_python_dependencies() -> None:
    missing = [name for name in ("numpy", "scipy", "h5py", "matplotlib", "vtk")
               if importlib.util.find_spec(name) is None]
    if missing:
        raise RuntimeError("missing Python packages: " + ", ".join(missing))


def check_generated_tables() -> None:
    names = (
        "ve_pseudo_pc_sg_axis.txt",
        "ve_pseudo_pc_pressure_axis.txt",
        "ve_pseudo_pc_values.txt",
    )
    with tempfile.TemporaryDirectory(prefix="thomas_ve_check_") as temporary:
        output = Path(temporary)
        run([sys.executable, "generate_vertical_equilibrium_pseudocapillary.py",
             "--output-dir", str(output)])
        different = [name for name in names
                     if (output / name).read_bytes() != (ROOT / "tables" / name).read_bytes()]
        if different:
            raise RuntimeError("checked-in VE tables differ from regenerated files: " + ", ".join(different))


def read_recovery(path: Path) -> dict[float, float]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines()
             if not line.startswith("#")]
    rows = csv.DictReader(lines)
    return {
        float(row["time_years"]): float(row["matrix_oil_component_mass_recovery_pct"])
        for row in rows
    }


def verify_recovery(case: str, path: Path) -> None:
    recovery = read_recovery(path)
    for time, expected in EXPECTED_RECOVERY[case].items():
        actual = recovery[time]
        if abs(actual - expected) > 0.001:
            raise RuntimeError(
                f"{case} recovery mismatch at {time} years: expected {expected:.4f}%, got {actual:.4f}%"
            )
        print(f"{case} {time:.1f} years: {actual:.4f}% (expected {expected:.4f}%)")


def main() -> None:
    repository = find_repository()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geosx", type=Path, default=repository / "build/bin/geosx")
    parser.add_argument("--output-root", type=Path, default=Path("/tmp/thomas_single_block_reproduction"))
    parser.add_argument("--run", action="store_true", help="run both 2.5-year simulations after dependency checks")
    args = parser.parse_args()

    require_files(list(DECKS.values()) + INPUT_FILES, "validation inputs")
    require_files([repository / path for path in SOURCE_FILES], "required GEOS implementation files")
    require_files([args.geosx], "GEOS executable")
    check_python_dependencies()
    check_generated_tables()
    run([sys.executable, "prepare_thomas_single_block.py"])

    validation_root = args.output_root / "input_validation"
    for case, deck in DECKS.items():
        run([str(args.geosx), "-i", deck.name, "-o", str(validation_root / case), "-v"])

    if not args.run:
        print("All inputs, source extensions, generated tables, Python dependencies, and XML decks are available.")
        print("Run again with --run to execute and verify both canonical simulations.")
        return

    for case, deck in DECKS.items():
        simulation_dir = args.output_root / case
        analysis_dir = args.output_root / "analysis" / case
        figures_dir = args.output_root / "figures" / case
        run([str(args.geosx), "-i", deck.name, "-o", str(simulation_dir)])
        run([
            sys.executable,
            "analyze_results.py",
            "--output", str(simulation_dir),
            "--analysis-dir", str(analysis_dir),
            "--figures-dir", str(figures_dir),
        ])
        verify_recovery(case, analysis_dir / "recovery_results.csv")

    oracle = args.output_root / "analysis/thomas_single_cell_oracle_ve.csv"
    run([sys.executable, "thomas_single_cell_oracle.py", "--output", str(oracle)])
    print(f"Reproduction completed successfully under {args.output_root}")


if __name__ == "__main__":
    main()
