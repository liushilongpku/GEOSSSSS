#!/usr/bin/env python3
"""Audit the Thomas 1983 dual-continuum single-block initialization and depletion.

Key information: 5545-psig live-oil matrix at the Thomas Table 4 state, a gas-filled
fracture continuum, and a 0.75-psi/day fracture depletion. The reference is the
Thomas 1983 Fig. 4 ultimate recovery of about 46% at 2.5 years for a 10-ft block.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import xml.etree.ElementTree as ET


INITIAL_PRESSURE_PA = 38_332_780.5
PSI_TO_PA = 6_894.757293168
DEPLETION_PSI_PER_DAY = 0.75
YEARS_TO_SECONDS = 31_536_000.0
MAX_TIME_YEARS = 2.5
INITIAL_COMPONENTS = (0.53678013, 0.16595388, 0.29726599)
INITIAL_GAS_OIL_PC_PA = -0.74 * PSI_TO_PA
ROOT = Path(__file__).resolve().parent.parent
DECKS = {
    "gdp_off": ROOT / "thomas_singleblock_gas_oil_gravity_drainage.xml",
    "gdp_on": ROOT / "thomas_singleblock_gas_oil_gravity_drainage_gdp_on.xml",
}
EXPECTED_CROSSFLOW = {
    "matrixControlledExchangeUpwinding": "1",
    "matrixControlledReverseExchangeRelPerm": "{ -1, 0.42, 0.03 }",
    "shapeFactorType": "direct",
    "shapeFactorValue": "0.215278208",
}


def audit_deck(label: str, deck: Path) -> tuple[float, float, float]:
    root = ET.parse(deck).getroot()
    fields = root.find("FieldSpecifications")
    if fields is None:
        raise RuntimeError("FieldSpecifications is missing")
    values = []
    for name in ("matrixOilFrac", "matrixGasFrac", "matrixWaterFrac"):
        field = fields.find(f"FieldSpecification[@name='{name}']")
        if field is None:
            raise RuntimeError(f"missing {name}")
        values.append(float(field.attrib["scale"]))
    components = tuple(values)
    crossflow = root.find("./Solvers/CompositionalMultiPhaseDualContinuumFVM/DualContinuumCrossFlow")
    if crossflow is None:
        raise RuntimeError(f"{deck}: DualContinuumCrossFlow is missing")
    expected = EXPECTED_CROSSFLOW | {"gravityDrainageFlag": "1" if label == "gdp_on" else "0"}
    for name, value in expected.items():
        if crossflow.attrib.get(name) != value:
            raise RuntimeError(f"{deck}: expected {name}={value!r}, got {crossflow.attrib.get(name)!r}")
    fluid = root.find("./Constitutive/BlackOilFluid")
    relperm = root.find("./Constitutive/TableRelativePermeability[@name='relperm']")
    if fluid is None or fluid.attrib.get("phaseNames") != "{ oil, gas, water }":
        raise RuntimeError(f"{deck}: BlackOilFluid phase order must be oil, gas, water")
    if relperm is None or relperm.attrib.get("phaseNames") != "{ oil, gas, water }":
        raise RuntimeError(f"{deck}: matrix relative-permeability phase order must be oil, gas, water")
    if relperm.attrib.get("threePhaseInterpolator") != "STONEII":
        raise RuntimeError(f"{deck}: matrix relative permeability must use STONEII")
    return components


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic-source", choices=DECKS)
    parser.add_argument("--diagnostic-output", type=Path)
    parser.add_argument("--vtk-frequency", type=float, default=788_400.0)
    parser.add_argument("--vertical-cells", type=int, default=1)
    parser.add_argument("--fracture-pseudocapillary", action="store_true")
    parser.add_argument("--free-fracture-composition", action="store_true")
    parser.add_argument("--initial-fracture-pc-equilibrium", action="store_true")
    parser.add_argument("--direct-shape-factor-no-axis-average", action="store_true")
    parser.add_argument("--no-matrix-pseudocapillary-pressure-scaling", action="store_true")
    parser.add_argument("--vertical-equilibrium-pseudocapillary", action="store_true")
    args = parser.parse_args()
    if bool(args.diagnostic_source) != bool(args.diagnostic_output):
        parser.error("--diagnostic-source and --diagnostic-output must be provided together")
    if args.vertical_cells < 1:
        parser.error("--vertical-cells must be positive")
    if args.vertical_cells != 1 and not args.diagnostic_source:
        parser.error("--vertical-cells requires --diagnostic-source and --diagnostic-output")

    for label, deck in DECKS.items():
        components = audit_deck(label, deck)
        assert abs(sum(components) - 1.0) < 1.0e-9
        assert all(abs(a - b) < 1.0e-9 for a, b in zip(components, INITIAL_COMPONENTS))
        print(f"{label} deck = {deck.name}")
        print("  matrix global component fractions = " + ", ".join(f"{value:.8f}" for value in components))

    end_time_s = MAX_TIME_YEARS * YEARS_TO_SECONDS
    end_pressure_pa = INITIAL_PRESSURE_PA - DEPLETION_PSI_PER_DAY * PSI_TO_PA * 365.0 * MAX_TIME_YEARS
    pressure_units_file = ROOT / "reference" / "thomas_fig4_3d_model.csv"
    print(f"initial pressure = {INITIAL_PRESSURE_PA:.4f} Pa")
    print(f"  5545 psig -> Pa = {5545.0 * PSI_TO_PA:.4f} Pa (atm offset not included)")
    print(f"0.75 psi/day over {MAX_TIME_YEARS} years = {DEPLETION_PSI_PER_DAY * PSI_TO_PA * 365.0 * MAX_TIME_YEARS:.2f} Pa")
    print(f"end pressure = {end_pressure_pa:.2f} Pa at {end_time_s:.1f} s")
    print(f"Thomas Fig. 4 reference file = {pressure_units_file} (exists={pressure_units_file.is_file()})")

    if args.diagnostic_source:
        source = DECKS[args.diagnostic_source]
        tree = ET.parse(source)
        vtk_event = tree.getroot().find("./Events/PeriodicEvent[@name='vtkOutput']")
        if vtk_event is None:
            raise RuntimeError(f"{source}: vtkOutput event is missing")
        vtk_event.attrib["timeFrequency"] = f"{args.vtk_frequency:.16g}"
        fluid = tree.getroot().find("./Constitutive/BlackOilFluid")
        if fluid is None:
            raise RuntimeError(f"{source}: BlackOilFluid is missing")
        table_files = [item.strip() for item in fluid.attrib["tableFiles"].strip("{} ").split(",")]
        fluid.attrib["tableFiles"] = "{ " + ", ".join(str((source.parent / item).resolve()) for item in table_files) + " }"
        for table in tree.getroot().findall("./Functions/TableFunction[@coordinateFiles]"):
            files = [item.strip() for item in table.attrib["coordinateFiles"].strip("{} ").split(",")]
            table.attrib["coordinateFiles"] = "{ " + ", ".join(
                str((source.parent / item).resolve()) if not Path(item).is_absolute() else item for item in files
            ) + " }"
            voxel_file = Path(table.attrib["voxelFile"])
            if not voxel_file.is_absolute():
                table.attrib["voxelFile"] = str((source.parent / voxel_file).resolve())
        for mesh in tree.getroot().findall("./Mesh/InternalMesh"):
            mesh.attrib["nz"] = f"{{ {args.vertical_cells} }}"
        if args.fracture_pseudocapillary:
            fracture_capillary = tree.getroot().find("./Constitutive/TableCapillaryPressure[@name='fractureCapPressure']")
            if fracture_capillary is None:
                raise RuntimeError(f"{source}: fractureCapPressure is missing or has an unexpected type")
            fracture_capillary.tag = "PressureScaledTableCapillaryPressure"
            fracture_capillary.attrib["nonWettingIntermediateCapPressureTableName"] = "gasCapPresTable"
            fracture_capillary.attrib["pressureScalingTableName"] = "gasOilSurfaceTensionScaling"
        if args.free_fracture_composition:
            fields = tree.getroot().find("./FieldSpecifications")
            if fields is None:
                raise RuntimeError(f"{source}: FieldSpecifications is missing")
            for name in (
                "fractureBoundaryGasFrac",
                "fractureBoundaryOilFrac",
                "fractureBoundaryWaterFrac",
            ):
                field = fields.find(f"FieldSpecification[@name='{name}']")
                if field is None:
                    raise RuntimeError(f"{source}: {name} is missing")
                fields.remove(field)
        if args.initial_fracture_pc_equilibrium:
            fields = tree.getroot().find("./FieldSpecifications")
            functions = tree.getroot().find("./Functions")
            if fields is None or functions is None:
                raise RuntimeError(f"{source}: FieldSpecifications or Functions is missing")
            fracture_pressure = fields.find("FieldSpecification[@name='fracturePressure']")
            if fracture_pressure is None:
                raise RuntimeError(f"{source}: fracturePressure is missing")
            fracture_pressure.attrib["scale"] = f"{INITIAL_PRESSURE_PA + INITIAL_GAS_OIL_PC_PA:.16g}"
            pressure_decline = functions.find("TableFunction[@name='pressureDecline']")
            if pressure_decline is None:
                raise RuntimeError(f"{source}: pressureDecline is missing")
            values = [float(value.strip()) for value in pressure_decline.attrib["values"].strip("{} ").split(",")]
            pressure_decline.attrib["values"] = "{ " + ", ".join(
                f"{value + INITIAL_GAS_OIL_PC_PA:.16g}" for value in values
            ) + " }"
        if args.direct_shape_factor_no_axis_average:
            crossflow = tree.getroot().find("./Solvers/CompositionalMultiPhaseDualContinuumFVM/DualContinuumCrossFlow")
            if crossflow is None:
                raise RuntimeError(f"{source}: DualContinuumCrossFlow is missing")
            shape_factor = float(crossflow.attrib["shapeFactorValue"])
            crossflow.attrib["shapeFactorValue"] = f"{3.0 * shape_factor:.16g}"
        if args.no_matrix_pseudocapillary_pressure_scaling:
            capillary = tree.getroot().find(
                "./Constitutive/PressureScaledTableCapillaryPressure[@name='matrixCapPressure']"
            )
            if capillary is None:
                raise RuntimeError(f"{source}: matrixCapPressure is missing or has an unexpected type")
            capillary.tag = "TableCapillaryPressure"
            capillary.attrib.pop("pressureScalingTableName", None)
            capillary.attrib.pop("pressureDependentTableName", None)
        if args.vertical_equilibrium_pseudocapillary:
            if args.no_matrix_pseudocapillary_pressure_scaling:
                parser.error("--vertical-equilibrium-pseudocapillary conflicts with "
                             "--no-matrix-pseudocapillary-pressure-scaling")
            root = tree.getroot()
            capillary = root.find(
                "./Constitutive/PressureScaledTableCapillaryPressure[@name='matrixCapPressure']"
            )
            functions = root.find("./Functions")
            if capillary is None or functions is None:
                raise RuntimeError(f"{source}: matrixCapPressure or Functions is missing")
            capillary.attrib.pop("pressureScalingTableName", None)
            capillary.attrib["pressureDependentTableName"] = "gasOilVerticalEquilibriumPseudoPc"
            table_dir = source.parent / "tables"
            table = functions.find("TableFunction[@name='gasOilVerticalEquilibriumPseudoPc']")
            if table is None:
                table = ET.SubElement(functions, "TableFunction", {"name": "gasOilVerticalEquilibriumPseudoPc"})
            table.attrib.update({
                "coordinateFiles": "{ " + ", ".join(
                    str((table_dir / name).resolve())
                    for name in ("ve_pseudo_pc_sg_axis.txt", "ve_pseudo_pc_pressure_axis.txt")
                ) + " }",
                "voxelFile": str((table_dir / "ve_pseudo_pc_values.txt").resolve()),
                "interpolation": "linear",
            })
        args.diagnostic_output.parent.mkdir(parents=True, exist_ok=True)
        tree.write(args.diagnostic_output, encoding="unicode")
        if args.vertical_equilibrium_pseudocapillary:
            matrix_pc_mode = "vertical-equilibrium table"
        elif args.no_matrix_pseudocapillary_pressure_scaling:
            matrix_pc_mode = "unscaled"
        else:
            matrix_pc_mode = "scalar pressure scaling"
        print(
            f"diagnostic deck = {args.diagnostic_output} "
            f"(1x1x{args.vertical_cells}, fracture pseudo Pc={args.fracture_pseudocapillary}, "
            f"free fracture composition={args.free_fracture_composition}, "
            f"initial fracture Pc equilibrium={args.initial_fracture_pc_equilibrium}, "
            f"direct shape factor no axis average={args.direct_shape_factor_no_axis_average}, "
            f"matrix pseudo Pc mode={matrix_pc_mode}, "
            f"VE pseudo Pc={args.vertical_equilibrium_pseudocapillary}, "
            f"VTK every {args.vtk_frequency:g} s)"
        )


if __name__ == "__main__":
    main()
