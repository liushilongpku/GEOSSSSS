#!/usr/bin/env python3
"""Generate the SPE6-constrained 10-ft, 7x7x8 fine-grid case.

Key information: SPE6 supplies the 4500-psig/5-year physical target and
water/oil Pc, while Thomas Table 1/3 supplies the rock curves and thin-shell
construction. The 7x7x8 grid is a new spatial-resolution experiment.
"""

from __future__ import annotations

import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
VALIDATION_ROOT = ROOT.parent.parent
THOMAS_SOURCE = VALIDATION_ROOT / "fine_grid_thomas_depletion_liveoil/thomas_10ft_7x7x8_liveoil_depletion.xml"
P0_ROOT = VALIDATION_ROOT / "P0_thomas_single_block"
THOMAS_FINAL_ROOT = (
    VALIDATION_ROOT
    / "fine_grid_thomas_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone"
)
MAIN_DECK = ROOT / "spe6_10ft_7x7x8_zeroPcf.xml"
GRID_NX = 7
GRID_NY = 7
GRID_NZ = 8
SHELL_THICKNESS_FT = 0.001
M_TO_FT = 3.280839895013123
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


def replace_xml_section(text: str, start: str, end: str, replacement: str, label: str) -> str:
    """Replace one complete XML section used to define the spatial grid."""

    start_index = text.find(start)
    end_index = text.find(end, start_index)
    if start_index < 0 or end_index < 0:
        raise RuntimeError(f"could not find {label}")
    end_index += len(end)
    return text[:start_index] + replacement + text[end_index:]


def cell_block_text(names: list[str], indent: str = "        ") -> str:
    """Format generated cell-block names for a GEOS XML list."""

    rows = [names[index : index + 8] for index in range(0, len(names), 8)]
    return "\n".join(indent + ", ".join(row) + ("," if index < len(rows) - 1 else "")
                     for index, row in enumerate(rows))


def axis_nodes(count: int) -> str:
    """Return comma-separated node coordinates along one axis with `count` cells."""
    shell_m = SHELL_THICKNESS_FT / M_TO_FT
    interior_m = 10.0 / M_TO_FT
    interior_cells = count - 2
    nodes = [0.0, shell_m]
    nodes.extend(
        shell_m + interior_m * index / interior_cells for index in range(1, interior_cells)
    )
    nodes.extend((shell_m + interior_m, 2.0 * shell_m + interior_m))
    return ", ".join(f"{value:.10g}" for value in nodes)


def grid_sections() -> tuple[str, str]:
    """Build a 7x7x8 mesh with an internal-matrix block of 5x5x6 cells."""

    all_cells: list[str] = []
    matrix_cells: list[str] = []
    fracture_cells: list[str] = []
    for k in range(GRID_NZ):
        for j in range(GRID_NY):
            for i in range(GRID_NX):
                is_fracture = (
                    i in (0, GRID_NX - 1)
                    or j in (0, GRID_NY - 1)
                    or k in (0, GRID_NZ - 1)
                )
                name = f"{'f' if is_fracture else 'm'}_{i}_{j}_{k}"
                all_cells.append(name)
                (fracture_cells if is_fracture else matrix_cells).append(name)

    interior_count = (GRID_NX - 2) * (GRID_NY - 2) * (GRID_NZ - 2)
    shell_count = GRID_NX * GRID_NY * GRID_NZ - interior_count
    mesh = f'''  <!-- Outer 0.001-ft fracture shell plus an interior matrix block. -->
  <Mesh>
    <InternalMesh name="mesh1" elementTypes="{{ C3D8 }}"
      xCoords="{{ {axis_nodes(GRID_NX)} }}"
      yCoords="{{ {axis_nodes(GRID_NY)} }}"
      zCoords="{{ {axis_nodes(GRID_NZ)} }}"
      nx="{{ {', '.join(['1'] * GRID_NX)} }}"
      ny="{{ {', '.join(['1'] * GRID_NY)} }}"
      nz="{{ {', '.join(['1'] * GRID_NZ)} }}"
      cellBlockNames="{{
{cell_block_text(all_cells)}
      }}"/>
  </Mesh>'''
    regions = f'''  <!-- The outer shell has {shell_count} fracture cells; the {interior_count}-cell interior is the matrix. -->
  <ElementRegions>
    <CellElementRegion name="matrixRegion" meshBody="mesh1"
      cellBlocks="{{
{cell_block_text(matrix_cells)}
      }}"
      materialList="{{ matrixRock, fluid, matrixRelperm, matrixCapPressure }}"/>
    <CellElementRegion name="fractureRegion" meshBody="mesh1"
      cellBlocks="{{
{cell_block_text(fracture_cells)}
      }}"
      materialList="{{ fractureRock, fluid, fractureRelperm, fractureCapPressure }}"/>
  </ElementRegions>'''
    return mesh, regions


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
    mesh, regions = grid_sections()
    text = replace_xml_section(text, "  <Mesh>", "  </Mesh>", mesh, "7x7x8 mesh")
    text = replace_once(
        text,
        "  <!-- Table 1: x/y = 0.001, 1, 2, 4, 2, 1, 0.001 ft; z = 0.001, 1, 2, 2, 2, 2, 1, 0.001 ft. -->\n",
        "",
        "obsolete 7x7x8 mesh comment",
    )
    text = replace_xml_section(
        text, "  <ElementRegions>", "  </ElementRegions>", regions, "7x7x8 element regions"
    )
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
    text = replace_once(text, ' maxEventDt="2.0e5"', "", "maximum timestep")
    text = replace_once(
        text,
        '    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver"/>\n'
        '    <PeriodicEvent name="vtkOutput" timeFrequency="7.884e6" target="/Outputs/vtkOutput"/>',
        '    <!-- Write the current state before advancing it so VTK TIME and fields remain aligned. -->\n'
        '    <PeriodicEvent name="vtkOutput" timeFrequency="7.884e6" target="/Outputs/vtkOutput"/>\n'
        '    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver"/>',
        "VTK/solver event order",
    )
    text = replace_once(
        text,
        '    <!-- Update representative-cell collections before writing their time histories. -->\n'
        '    <PeriodicEvent name="matrixPhaseCollection" timeFrequency="7.884e6"\n'
        '      target="/Tasks/matrixPhaseCollection"/>\n'
        '    <PeriodicEvent name="matrixPhaseHistory" timeFrequency="7.884e6" targetExactTimestep="0"\n'
        '      target="/Outputs/matrixPhaseHistory"/>\n'
        '    <PeriodicEvent name="fracturePhaseCollection" timeFrequency="7.884e6"\n'
        '      target="/Tasks/fracturePhaseCollection"/>\n'
        '    <PeriodicEvent name="fracturePhaseHistory" timeFrequency="7.884e6" targetExactTimestep="0"\n'
        '      target="/Outputs/fracturePhaseHistory"/>\n',
        "",
        "non-VTK output events",
    )
    text = replace_once(
        text,
        '  <Tasks>\n'
        '    <!-- Representative cell histories supplement the all-cell VTK output. -->\n'
        '    <PackCollection name="matrixPhaseCollection"\n'
        '      objectPath="mesh1/ElementRegions/matrixRegion/m_1_1_1"\n'
        '      fieldName="phaseVolumeFraction"/>\n'
        '    <PackCollection name="fracturePhaseCollection"\n'
        '      objectPath="mesh1/ElementRegions/fractureRegion/f_0_0_0"\n'
        '      fieldName="phaseVolumeFraction"/>\n'
        '  </Tasks>\n\n',
        "",
        "non-VTK tasks",
    )
    text = replace_once(
        text,
        '    <TimeHistory name="matrixPhaseHistory" sources="{ /Tasks/matrixPhaseCollection }"\n'
        '      filename="partitioned_matrix_phase_history"/>\n'
        '    <TimeHistory name="fracturePhaseHistory" sources="{ /Tasks/fracturePhaseCollection }"\n'
        '      filename="partitioned_fracture_phase_history"/>\n',
        "",
        "non-VTK outputs",
    )
    return text


def main() -> None:
    """Generate the self-contained VTK-only main deck and support files."""

    copy_support_files()
    text = build_main_deck()
    MAIN_DECK.write_text(text, encoding="utf-8")
    print(f"generated {MAIN_DECK}")


if __name__ == "__main__":
    main()
