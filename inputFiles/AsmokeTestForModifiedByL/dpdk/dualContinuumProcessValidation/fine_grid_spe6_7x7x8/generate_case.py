#!/usr/bin/env python3
"""Generate the SPE6-constrained 10-ft, 7x7x8 fine-grid case.

Key information: SPE6 supplies the 4500-psig/5-year physical target and
water/oil Pc, while Thomas Table 1/3 supplies the fine grid and rock curves.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
VALIDATION_ROOT = ROOT.parent
THOMAS_SOURCE = VALIDATION_ROOT / "fine_grid_thomas_depletion_liveoil/thomas_10ft_7x7x8_liveoil_depletion.xml"
P0_ROOT = VALIDATION_ROOT / "P0_thomas_single_block"
THOMAS_FINAL_ROOT = (
    VALIDATION_ROOT
    / "fine_grid_thomas_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone"
)
MAIN_DECK = ROOT / "spe6_10ft_7x7x8_zeroPcf.xml"
YEAR = 31_536_000.0
REPORT_INTERVAL_YEARS = 0.25
END_YEARS = 5.0
PRESSURE_PA = 31_127_732.819256
PC_SCALE = 5.2442396313
NEWTON_TOL = "1.0e-3"
NEWTON_MAX_ITER = 80
BLOCK_CENTER_Z_M = 1.5243048
INITIAL_OIL_DENSITY_KG_PER_M3 = 607.43596104
INITIAL_GAS_DENSITY_KG_PER_M3 = 207.52880055
INITIAL_WATER_DENSITY_KG_PER_M3 = 978.6686918484461
GRAVITY_M_PER_S2 = 9.81
FRACTURE_OIL_SATURATION = 0.005
FRACTURE_GAS_SATURATION = 0.990
FRACTURE_WATER_SATURATION = 0.005
FRACTURE_DISSOLVED_GAS_TO_OIL_MASS_RATIO = 0.228027729306525
FRACTURE_OIL_PHASE_MASS = FRACTURE_OIL_SATURATION * INITIAL_OIL_DENSITY_KG_PER_M3
FRACTURE_GAS_PHASE_MASS = FRACTURE_GAS_SATURATION * INITIAL_GAS_DENSITY_KG_PER_M3
FRACTURE_WATER_PHASE_MASS = FRACTURE_WATER_SATURATION * INITIAL_WATER_DENSITY_KG_PER_M3
FRACTURE_TOTAL_MASS = FRACTURE_OIL_PHASE_MASS + FRACTURE_GAS_PHASE_MASS + FRACTURE_WATER_PHASE_MASS
FRACTURE_OIL_COMPONENT = (
    FRACTURE_OIL_PHASE_MASS
    / (1.0 + FRACTURE_DISSOLVED_GAS_TO_OIL_MASS_RATIO)
    / FRACTURE_TOTAL_MASS
)
FRACTURE_GAS_COMPONENT = (
    FRACTURE_GAS_PHASE_MASS
    + FRACTURE_OIL_PHASE_MASS
    * FRACTURE_DISSOLVED_GAS_TO_OIL_MASS_RATIO
    / (1.0 + FRACTURE_DISSOLVED_GAS_TO_OIL_MASS_RATIO)
) / FRACTURE_TOTAL_MASS
FRACTURE_WATER_COMPONENT = FRACTURE_WATER_PHASE_MASS / FRACTURE_TOTAL_MASS
SATURATED_WATER_MASS_FRACTION = "0.287145405913946-7.95155908893894e-06*z"
SATURATED_DISSOLVED_GAS_RATIO = "0.228027729306525-5.75194103403429e-05*z"
SATURATED_OIL_MASS_FRACTION = (
    f"(1-({SATURATED_WATER_MASS_FRACTION}))/(1+({SATURATED_DISSOLVED_GAS_RATIO}))"
)
SATURATED_GAS_MASS_FRACTION = (
    f"(1-({SATURATED_WATER_MASS_FRACTION}))*({SATURATED_DISSOLVED_GAS_RATIO})"
    f"/(1+({SATURATED_DISSOLVED_GAS_RATIO}))"
)
SATURATED_OIL_COMPONENT = SATURATED_OIL_MASS_FRACTION
SATURATED_GAS_COMPONENT = SATURATED_GAS_MASS_FRACTION
SATURATED_WATER_COMPONENT = SATURATED_WATER_MASS_FRACTION

BASE_GAS_PC_PA = (517.1, 586.1, 655.0, 792.9, 999.7, 1758.2, 2661.4)
ORIGINAL_KROW = (
    '    <TableFunction name="oilRelPermTableForOW" coordinates="{ 0.25, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80 }"\n'
    '      values="{ 0.0, 0.042, 0.154, 0.304, 0.392, 0.492, 0.600, 0.723, 0.860, 1.000 }" interpolation="linear"/>'
)
THOMAS_DRAINAGE_KROW = (
    '    <!-- Thomas Eq. (25): at Sw approximately Swc, krow remains at its connate-water endpoint. -->\n'
    '    <TableFunction name="oilRelPermTableForOW" coordinates="{ 0.0, 0.001, 1.0 }"\n'
    '      values="{ 0.0, 0.998, 1.0 }" interpolation="linear"/>'
)


def replace_once(text: str, old: str, new: str, label: str) -> str:
    """Replace one authoritative fragment and fail if its source changed."""

    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"expected one {label}, found {count}")
    return text.replace(old, new)


def copy_support_files() -> None:
    """Copy the audited PVT tables and digitized SPE6 comparison curve."""

    tables = ROOT / "tables"
    reference = ROOT / "reference"
    tables.mkdir(parents=True, exist_ok=True)
    reference.mkdir(parents=True, exist_ok=True)
    sources = {
        "pvto_bo.txt": THOMAS_FINAL_ROOT / "tables/pvto_surface_condition.txt",
        "pvtg_norv_bo.txt": THOMAS_FINAL_ROOT / "tables/pvtg_absolute_pressure.txt",
        "pvtw_bo.txt": THOMAS_FINAL_ROOT / "tables/pvtw_absolute_pressure.txt",
    }
    for name, source in sources.items():
        shutil.copy2(source, tables / name)
    shutil.copy2(
        P0_ROOT / "reference/spe_fig1_gas_oil_zeroPcf.csv",
        reference / "spe_fig1_gas_oil_zeroPcf.csv",
    )


def build_main_deck() -> str:
    """Transform the audited Thomas fine-grid deck into the SPE6 case."""

    text = THOMAS_SOURCE.read_text(encoding="utf-8")
    text = replace_once(
        text,
        "<!-- Purpose: implement the Thomas 1983 10-ft, 7x7x8 live-oil gravity-drainage fine-grid case. Key information: 5545-psig initial state, 0.75-psi/day pressure depletion, gas-filled outer fracture-valued cells. -->",
        "<!-- Purpose: approximate the SPE6 zero-fracture-Pc single-block problem with a 7x7x8 fine grid. Key information: 4500-psig live oil, fixed fracture-reservoir shell, and 5-year recovery. -->",
        "deck purpose header",
    )
    text = text.replace("pvto_bo_surface_condition.txt", "tables/pvto_bo.txt")
    text = text.replace("../tables/pvtg_norv_bo.txt", "tables/pvtg_norv_bo.txt")
    text = text.replace("../tables/pvtw_bo.txt", "tables/pvtw_bo.txt")
    text = replace_once(
        text,
        '    <PressurePorosity name="matrixPorosity" defaultReferencePorosity="0.30"\n'
        '      referencePressure="3.83327805e7" compressibility="5.076e-10"/>',
        '    <PressurePorosity name="matrixPorosity" defaultReferencePorosity="0.29"\n'
        f'      referencePressure="{PRESSURE_PA:.8g}" compressibility="5.076e-10"/>',
        "matrix porosity",
    )
    text = replace_once(
        text,
        '    <PressurePorosity name="fracturePorosity" defaultReferencePorosity="1.0"\n'
        '      referencePressure="3.83327805e7" compressibility="5.076e-10"/>',
        '    <PressurePorosity name="fracturePorosity" defaultReferencePorosity="1.0"\n'
        f'      referencePressure="{PRESSURE_PA:.8g}" compressibility="5.076e-10"/>',
        "shell porosity",
    )
    text = replace_once(
        text,
        ORIGINAL_KROW,
        THOMAS_DRAINAGE_KROW,
        "connate-water oil relative-permeability branch",
    )
    text = replace_once(
        text,
        'newtonTol="1.0e-3" newtonMaxIter="80"',
        f'newtonTol="{NEWTON_TOL}" newtonMaxIter="{NEWTON_MAX_ITER}"',
        "Newton tolerance",
    )
    old_initial = '''    <!-- Initial matrix: Thomas Table 4 state at 5545 psig; gas component is dissolved gas. -->
    <FieldSpecification name="initialMatrixPressure" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="pressure"
      setNames="{ all }" scale="3.83327805e7"/>
    <FieldSpecification name="initialMatrixOilFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="0" setNames="{ all }" scale="0.53683289"/>
    <FieldSpecification name="initialMatrixGasFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="1" setNames="{ all }" scale="0.16590112"/>
    <FieldSpecification name="initialMatrixWaterFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="2" setNames="{ all }" scale="0.29726599"/>

    <!-- Initial fracture shell: gas-filled and at the same pressure, matching the external gas boundary. -->
    <FieldSpecification name="initialFracturePressure" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="pressure"
      setNames="{ all }" scale="3.83327805e7"/>'''
    new_initial = f'''    <!-- Initial matrix: 4500-psig live oil with dissolved gas but no intended free-gas phase. -->
    <FieldSpecification name="initialMatrixPressure" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="pressure"
      setNames="{{ all }}" functionName="initialMatrixHydrostatic" scale="1.0"/>
    <FieldSpecification name="initialMatrixOilFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="0" setNames="{{ all }}" functionName="initialSaturatedOilComponent" scale="1.0"/>
    <FieldSpecification name="initialMatrixGasFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="1" setNames="{{ all }}" functionName="initialSaturatedGasComponent" scale="1.0"/>
    <FieldSpecification name="initialMatrixWaterFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/matrixRegion" fieldName="globalCompFraction"
      component="2" setNames="{{ all }}" functionName="initialSaturatedWaterComponent" scale="1.0"/>

    <!-- The 0.001-ft shell represents fractures initialized at
         Sg=0.990, So=0.005, and Sw=0.005. The global component fields below
         are the equivalent mass fractions at the 4500-psig center datum. -->
    <FieldSpecification name="initialFracturePressure" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="pressure"
      setNames="{{ all }}" functionName="fractureHydrostatic" scale="1.0"/>'''
    text = replace_once(text, old_initial, new_initial, "initial matrix and shell state")
    text = replace_once(
        text,
        '''    <FieldSpecification name="initialFractureOilFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="0" setNames="{ all }" scale="0.0"/>
    <FieldSpecification name="initialFractureGasFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="1" setNames="{ all }" scale="1.0"/>
    <FieldSpecification name="initialFractureWaterFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="2" setNames="{ all }" scale="0.0"/>''',
        f'''    <FieldSpecification name="initialFractureOilFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="0" setNames="{{ all }}" scale="{FRACTURE_OIL_COMPONENT:.16g}"/>
    <FieldSpecification name="initialFractureGasFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="1" setNames="{{ all }}" scale="{FRACTURE_GAS_COMPONENT:.16g}"/>
    <FieldSpecification name="initialFractureWaterFrac" initialCondition="1" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="2" setNames="{{ all }}" scale="{FRACTURE_WATER_COMPONENT:.16g}"/>''',
        "fracture saturation-equivalent component fractions",
    )

    old_boundary = '''    <!-- Thomas depletion: all six external faces follow 0.75 psi/day pressure decline. -->
    <FieldSpecification name="depletionPressure" targetMesh="mesh1" objectPath="faceManager"
      fieldName="pressure" setNames="{ xneg, xpos, yneg, ypos, zneg, zpos }"
      functionName="pressureDecline" scale="1.0"/>
    <FieldSpecification name="boundaryOil" targetMesh="mesh1" objectPath="faceManager"
      fieldName="globalCompFraction" component="0" setNames="{ xneg, xpos, yneg, ypos, zneg, zpos }" scale="0.0"/>
    <FieldSpecification name="boundaryGas" targetMesh="mesh1" objectPath="faceManager"
      fieldName="globalCompFraction" component="1" setNames="{ xneg, xpos, yneg, ypos, zneg, zpos }" scale="1.0"/>
    <FieldSpecification name="boundaryWater" targetMesh="mesh1" objectPath="faceManager"
      fieldName="globalCompFraction" component="2" setNames="{ xneg, xpos, yneg, ypos, zneg, zpos }" scale="0.0"/>
    <FieldSpecification name="boundaryTemperature" targetMesh="mesh1" objectPath="faceManager"
      fieldName="temperature" setNames="{ xneg, xpos, yneg, ypos, zneg, zpos }" scale="366.0"/>'''
    new_boundary = f'''    <!-- Treat the complete thin shell as an infinite fracture reservoir.
         GEOS constrains pressure and global component mass fractions; the
         resulting phase state is approximately Sg=0.990, So=Sw=0.005. -->
    <FieldSpecification name="fractureBoundaryPressure" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="pressure"
      setNames="{{ all }}" functionName="fractureHydrostatic" scale="1.0"/>
    <FieldSpecification name="fractureBoundaryOil" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="0" setNames="{{ all }}" scale="{FRACTURE_OIL_COMPONENT:.16g}"/>
    <FieldSpecification name="fractureBoundaryGas" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="1" setNames="{{ all }}" scale="{FRACTURE_GAS_COMPONENT:.16g}"/>
    <FieldSpecification name="fractureBoundaryWater" targetMesh="mesh1"
      objectPath="ElementRegions/fractureRegion" fieldName="globalCompFraction"
      component="2" setNames="{{ all }}" scale="{FRACTURE_WATER_COMPONENT:.16g}"/>'''
    text = replace_once(text, old_boundary, new_boundary, "constant fracture-shell state")

    hydrostatic_functions = f'''  <Functions>
    <!-- 4500 psig is the block-center datum; each phase starts with its PVT density gradient. -->
    <SymbolicFunction name="initialMatrixHydrostatic" inputVarNames="{{ elementCenter }}"
      variableNames="{{ x, y, z }}"
      expression="{PRESSURE_PA:.12g}+{INITIAL_OIL_DENSITY_KG_PER_M3 * GRAVITY_M_PER_S2:.12g}*({BLOCK_CENTER_Z_M:.10g}-z)"/>
    <!-- These mass fractions give So=0.8 and Sw=0.2 at the local saturated Rs,
         with no deliberate free-gas seed. -->
    <SymbolicFunction name="initialSaturatedOilComponent" inputVarNames="{{ elementCenter }}"
      variableNames="{{ x, y, z }}" expression="{SATURATED_OIL_COMPONENT}"/>
    <SymbolicFunction name="initialSaturatedGasComponent" inputVarNames="{{ elementCenter }}"
      variableNames="{{ x, y, z }}" expression="{SATURATED_GAS_COMPONENT}"/>
    <SymbolicFunction name="initialSaturatedWaterComponent" inputVarNames="{{ elementCenter }}"
      variableNames="{{ x, y, z }}" expression="{SATURATED_WATER_COMPONENT}"/>
    <SymbolicFunction name="fractureHydrostatic" inputVarNames="{{ elementCenter }}"
      variableNames="{{ x, y, z }}"
      expression="{PRESSURE_PA:.12g}+{INITIAL_GAS_DENSITY_KG_PER_M3 * GRAVITY_M_PER_S2:.12g}*({BLOCK_CENTER_Z_M:.10g}-z)"/>
'''
    text = replace_once(text, "  <Functions>", hydrostatic_functions, "hydrostatic pressure functions")

    old_water_pc = '''    <TableFunction name="waterCapPresTable" coordinates="{ 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.75 }"
      values="{ 344738, 62052.84, 13789.52, 3447.38, 0, -2757.90, -8273.71, -27579.04, -68947.6, -275790.4 }" interpolation="linear"/>'''
    new_water_pc = '''    <!-- SPE6 1990 Table 1 replaces the Thomas water/oil capillary-pressure table. -->
    <TableFunction name="waterCapPresTable" coordinates="{ 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.75 }"
      values="{ 6894.7573, 3447.3786, 2068.4272, 1034.2136, 0, -1378.9515, -8273.7088, -27579.0292, -68947.5729, -275790.2917 }" interpolation="linear"/>'''
    text = replace_once(text, old_water_pc, new_water_pc, "SPE6 water/oil capillary pressure")

    scaled_pc = ", ".join(f"{value * PC_SCALE:.10g}" for value in BASE_GAS_PC_PA)
    text = replace_once(
        text,
        '    <TableFunction name="gasCapPresTable" coordinates="{ 0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55 }"\n'
        '      values="{ 517.1, 586.1, 655.0, 792.9, 999.7, 1758.2, 2661.4 }" interpolation="linear"/>',
        '    <!-- Thomas Table 3 rock Pc scaled by sigma(4500 psig)/sigma(5545 psig). -->\n'
        '    <TableFunction name="gasCapPresTable" coordinates="{ 0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55 }"\n'
        f'      values="{{ {scaled_pc} }}" interpolation="linear"/>',
        "pressure-scaled matrix gas/oil capillary pressure",
    )
    text = replace_once(
        text,
        '    <TableFunction name="pressureDecline" inputVarNames="{time}"\n'
        '      coordinates="{ 0.0, 157680000.0 }"\n'
        '      values="{ 38332780.5, 28895581.5 }" interpolation="linear"/>\n',
        "",
        "Thomas pressure-decline function",
    )
    text = replace_once(
        text,
        '''    <TableFunction name="fractureWaterRelPermTable" coordinates="{ 0.0, 1.0 }"
      values="{ 0.0, 1.0 }" interpolation="linear"/>
    <TableFunction name="fractureOilRelPermTableForOW" coordinates="{ 0.0, 1.0 }"
      values="{ 0.0, 1.0 }" interpolation="linear"/>
    <TableFunction name="fractureGasRelPermTable" coordinates="{ 0.0, 1.0 }"
      values="{ 0.0, 1.0 }" interpolation="linear"/>
    <TableFunction name="fractureOilRelPermTableForOG" coordinates="{ 0.0, 1.0 }"
      values="{ 0.0, 1.0 }" interpolation="linear"/>''',
        f'''    <!-- The 0.005 oil/water phase seeds are residual and cannot feed the matrix. -->
    <TableFunction name="fractureWaterRelPermTable" coordinates="{{ {FRACTURE_WATER_SATURATION}, 1.0 }}"
      values="{{ 0.0, 1.0 }}" interpolation="linear"/>
    <TableFunction name="fractureOilRelPermTableForOW" coordinates="{{ {FRACTURE_OIL_SATURATION}, 1.0 }}"
      values="{{ 0.0, 1.0 }}" interpolation="linear"/>
    <TableFunction name="fractureGasRelPermTable" coordinates="{{ 0.0, 1.0 }}"
      values="{{ 0.0, 1.0 }}" interpolation="linear"/>
    <TableFunction name="fractureOilRelPermTableForOG" coordinates="{{ {FRACTURE_OIL_SATURATION}, 1.0 }}"
      values="{{ 0.0, 1.0 }}" interpolation="linear"/>''',
        "fracture residual-phase relative permeability",
    )
    text = replace_once(text, '  <Events maxTime="1.5768e8">', '  <Events maxTime="157680000.0">', "end time")
    text = replace_once(text, 'maxEventDt="2.0e5"', 'maxEventDt="5.0e4"', "maximum timestep")
    text = replace_once(
        text,
        "  </Events>",
        '    <PeriodicEvent name="restartOutput" timeFrequency="7884000.0" targetExactTimestep="1"\n'
        '      target="/Outputs/restartOutput"/>\n'
        "  </Events>",
        "restart event",
    )
    text = replace_once(
        text,
        '    <TimeHistory name="fracturePhaseHistory" sources="{ /Tasks/fracturePhaseCollection }"\n'
        '      filename="partitioned_fracture_phase_history"/>\n'
        "  </Outputs>",
        '    <TimeHistory name="fracturePhaseHistory" sources="{ /Tasks/fracturePhaseCollection }"\n'
        '      filename="partitioned_fracture_phase_history"/>\n'
        '    <Restart name="restartOutput"/>\n'
        "  </Outputs>",
        "restart output",
    )
    return text


def generate_segments(main_text: str) -> None:
    """Generate 0.25-year restart segments without changing the physics."""

    decks = ROOT / "decks"
    decks.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int | float]] = []
    segment_count = round(END_YEARS / REPORT_INTERVAL_YEARS)
    for index in range(segment_count):
        start_years = index * REPORT_INTERVAL_YEARS
        end_years = (index + 1) * REPORT_INTERVAL_YEARS
        end_time = end_years * YEAR
        stem = f"segment_{index:02d}_{start_years:.2f}_to_{end_years:.2f}yr".replace(".", "p")
        text = main_text
        text = text.replace("tables/pvto_bo.txt", "../tables/pvto_bo.txt")
        text = text.replace("tables/pvtg_norv_bo.txt", "../tables/pvtg_norv_bo.txt")
        text = text.replace("tables/pvtw_bo.txt", "../tables/pvtw_bo.txt")
        text = replace_once(text, 'maxTime="157680000.0"', f'maxTime="{end_time:.1f}"', "segment end time")
        text = replace_once(
            text,
            'name="restartOutput" timeFrequency="7884000.0"',
            f'name="restartOutput" timeFrequency="{end_time:.1f}"',
            "segment restart time",
        )
        header = (
            "<!-- Purpose: provide one reproducible SPE6 fine-grid restart segment. "
            f"Key information: interval={start_years:.2f}-{end_years:.2f} years; physical parameters are unchanged. -->\n"
        )
        relative_path = Path("decks") / f"{stem}.xml"
        (ROOT / relative_path).write_text(header + text, encoding="utf-8")
        rows.append(
            {
                "segment": index,
                "start_years": f"{start_years:.2f}",
                "end_years": f"{end_years:.2f}",
                "deck": str(relative_path),
            }
        )

    with (ROOT / "segment_manifest.csv").open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write("# Purpose: map each 0.25-year SPE6 reporting interval to its restart deck.\n")
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Generate the self-contained main deck, support files, and segments."""

    copy_support_files()
    text = build_main_deck()
    MAIN_DECK.write_text(text, encoding="utf-8")
    generate_segments(text)
    print(f"generated {MAIN_DECK}")


if __name__ == "__main__":
    main()
