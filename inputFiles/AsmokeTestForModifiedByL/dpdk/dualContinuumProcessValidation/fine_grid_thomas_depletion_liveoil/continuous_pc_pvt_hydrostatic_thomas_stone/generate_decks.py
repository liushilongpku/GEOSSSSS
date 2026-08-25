#!/usr/bin/env python3
"""Generate the final self-contained Thomas 10 ft fine-grid restart chain.

Key information: build directly from the authoritative base XML, use
right-endpoint Pc sampling, PVT-consistent gas hydrostatics, strict numerics,
and the Thomas Eq. (25)-consistent water/oil relative-permeability branch.
"""

from __future__ import annotations

import csv
from bisect import bisect_right
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "thomas_10ft_7x7x8_liveoil_depletion.xml"
PVTG = ROOT / "tables/pvtg_absolute_pressure.txt"
YEAR = 31_536_000.0
STEP_YEARS = 0.1
END_YEARS = 2.5
INITIAL_PRESSURE_PA = 38_332_780.5
INITIAL_PRESSURE_PSIG = 5545.0
PRESSURE_DECLINE_PSI_PER_DAY = 0.75
PSI_TO_PA = 6_894.757293168
TOP_DATUM_M = 3.0486096
GRAVITY_M_PER_S2 = 9.81
INITIAL_OIL_DENSITY_KG_PER_M3 = 578.39
GAS_SURFACE_DENSITY_KG_PER_M3 = 0.929
BASE_PC_PA = (517.1, 586.1, 655.0, 792.9, 999.7, 1758.2, 2661.4)
NEWTON_TOL = "1.0e-6"
MAX_EVENT_DT = "5.0e4"

ORIGINAL_KROW = (
    '    <TableFunction name="oilRelPermTableForOW" coordinates="{ 0.25, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80 }"\n'
    '      values="{ 0.0, 0.042, 0.154, 0.304, 0.392, 0.492, 0.600, 0.723, 0.860, 1.000 }" interpolation="linear"/>'
)
THOMAS_DRAINAGE_KROW = (
    '    <!-- Thomas Eq. (25): krow is a function of Sw, which remains approximately Swc=0.20 here. -->\n'
    '    <TableFunction name="oilRelPermTableForOW" coordinates="{ 0.0, 0.001, 1.0 }"\n'
    '      values="{ 0.0, 0.998, 1.0 }" interpolation="linear"/>'
)

# Thomas 1983 Table 2, ordered by increasing gauge pressure.
PRESSURE_PSIG = (1674.0, 2031.0, 2530.0, 2991.0, 3553.0, 4110.0, 4544.0, 4935.0, 5255.0, 5545.0)
SURFACE_TENSION_DYN_PER_CM = (6.0, 4.7, 3.3, 2.2, 1.28, 0.72, 0.444, 0.255, 0.155, 0.090)
INITIAL_SURFACE_TENSION = SURFACE_TENSION_DYN_PER_CM[-1]


def replace_once(text: str, old: str, new: str, label: str) -> str:
    """Replace one source fragment and fail if the source deck changed."""

    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"expected one {label}, found {count}")
    return text.replace(old, new)


def linear_interpolate(x: float, x_values: tuple[float, ...], y_values: tuple[float, ...]) -> float:
    """Linearly interpolate within a strictly increasing tabulation."""

    if not x_values[0] <= x <= x_values[-1]:
        raise ValueError(f"value {x} is outside [{x_values[0]}, {x_values[-1]}]")
    upper = min(bisect_right(x_values, x), len(x_values) - 1)
    lower = max(upper - 1, 0)
    if upper == lower:
        return y_values[lower]
    weight = (x - x_values[lower]) / (x_values[upper] - x_values[lower])
    return y_values[lower] + weight * (y_values[upper] - y_values[lower])


def read_pvdg() -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Read absolute pressure and Bg from the local Thomas gas table."""

    rows = []
    for line in PVTG.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        rows.append((float(fields[0]), float(fields[1])))
    return tuple(row[0] for row in rows), tuple(row[1] for row in rows)


PVDG_PRESSURE_PA, PVDG_BG = read_pvdg()


def pressure_psig_at(time_years: float) -> float:
    """Return the prescribed top-datum gauge pressure."""

    return INITIAL_PRESSURE_PSIG - PRESSURE_DECLINE_PSI_PER_DAY * 365.0 * time_years


def pressure_pa_at(time_years: float) -> float:
    """Return the prescribed top-datum absolute pressure."""

    return INITIAL_PRESSURE_PA - PRESSURE_DECLINE_PSI_PER_DAY * PSI_TO_PA * 365.0 * time_years


def gas_density_at(time_years: float) -> float:
    """Compute gas density from GEOS' rho_g=sc-density/Bg convention."""

    bg = linear_interpolate(pressure_pa_at(time_years), PVDG_PRESSURE_PA, PVDG_BG)
    return GAS_SURFACE_DENSITY_KG_PER_M3 / bg


def surface_tension_at(time_years: float) -> float:
    """Interpolate the published pressure-dependent surface tension."""

    return linear_interpolate(pressure_psig_at(time_years), PRESSURE_PSIG, SURFACE_TENSION_DYN_PER_CM)


def segment_name(index: int, start_years: float, end_years: float) -> str:
    """Build a sortable restart-segment name."""

    start_tag = f"{start_years:.1f}".replace(".", "p")
    end_tag = f"{end_years:.1f}".replace(".", "p")
    return f"segment_{index:03d}_{start_tag}_to_{end_tag}yr"


def hydrostatic_functions(start_years: float, end_years: float) -> str:
    """Return initial and time-varying PVT gas hydrostatic functions."""

    start_time = start_years * YEAR
    duration = (end_years - start_years) * YEAR
    decline_pa_per_s = PRESSURE_DECLINE_PSI_PER_DAY * PSI_TO_PA / 86_400.0
    matrix_gradient = INITIAL_OIL_DENSITY_KG_PER_M3 * GRAVITY_M_PER_S2
    initial_gas_density = gas_density_at(0.0)
    rho_start = gas_density_at(start_years)
    rho_end = gas_density_at(end_years)
    rho_slope = (rho_end - rho_start) / duration
    density_expression = f"({rho_start:.12g}+{rho_slope:.12g}*(t-{start_time:.1f}))"
    return (
        '    <SymbolicFunction name="initialMatrixHydrostatic" inputVarNames="{ elementCenter }"\n'
        '      variableNames="{ x, y, z }"\n'
        f'      expression="{INITIAL_PRESSURE_PA:.10g}+{matrix_gradient:.10g}*({TOP_DATUM_M:.10g}-z)"/>\n'
        '    <SymbolicFunction name="initialFractureHydrostatic" inputVarNames="{ elementCenter }"\n'
        '      variableNames="{ x, y, z }"\n'
        f'      expression="{INITIAL_PRESSURE_PA:.10g}+{initial_gas_density * GRAVITY_M_PER_S2:.12g}'
        f'*({TOP_DATUM_M:.10g}-z)"/>\n'
        '    <SymbolicFunction name="fractureBoundaryHydrostatic" inputVarNames="{ faceCenter, time }"\n'
        '      variableNames="{ x, y, z, t }"\n'
        f'      expression="{INITIAL_PRESSURE_PA:.10g}-{decline_pa_per_s:.12g}*t+'
        f'{GRAVITY_M_PER_S2:.10g}*{density_expression}*({TOP_DATUM_M:.10g}-z)"/>'
    )


def make_deck(name: str, scale: float, start_years: float, end_years: float, variant: str) -> None:
    """Create one final Thomas Eq. (25)-consistent restart deck."""

    end_time = end_years * YEAR
    text = SOURCE.read_text(encoding="utf-8")
    text = text.replace("pvto_bo_surface_condition.txt", "../../tables/pvto_surface_condition.txt")
    text = text.replace("../tables/pvtg_norv_bo.txt", "../../tables/pvtg_absolute_pressure.txt")
    text = text.replace("../tables/pvtw_bo.txt", "../../tables/pvtw_absolute_pressure.txt")
    text = replace_once(
        text,
        'values="{ 517.1, 586.1, 655.0, 792.9, 999.7, 1758.2, 2661.4 }"',
        'values="{ ' + ", ".join(f"{value * scale:.10g}" for value in BASE_PC_PA) + ' }"',
        "matrix gas/oil capillary-pressure values",
    )
    text = replace_once(
        text,
        '    <FieldSpecification name="initialMatrixPressure" initialCondition="1" targetMesh="mesh1"\n'
        '      objectPath="ElementRegions/matrixRegion" fieldName="pressure"\n'
        '      setNames="{ all }" scale="3.83327805e7"/>',
        '    <FieldSpecification name="initialMatrixPressure" initialCondition="1" targetMesh="mesh1"\n'
        '      objectPath="ElementRegions/matrixRegion" fieldName="pressure"\n'
        '      setNames="{ all }" functionName="initialMatrixHydrostatic" scale="1.0"/>',
        "initial matrix pressure",
    )
    text = replace_once(
        text,
        '    <FieldSpecification name="initialFracturePressure" initialCondition="1" targetMesh="mesh1"\n'
        '      objectPath="ElementRegions/fractureRegion" fieldName="pressure"\n'
        '      setNames="{ all }" scale="3.83327805e7"/>',
        '    <FieldSpecification name="initialFracturePressure" initialCondition="1" targetMesh="mesh1"\n'
        '      objectPath="ElementRegions/fractureRegion" fieldName="pressure"\n'
        '      setNames="{ all }" functionName="initialFractureHydrostatic" scale="1.0"/>',
        "initial fracture pressure",
    )
    text = replace_once(
        text,
        '      functionName="pressureDecline" scale="1.0"/>',
        '      functionName="fractureBoundaryHydrostatic" scale="1.0"/>',
        "depletion pressure function",
    )
    text = replace_once(
        text,
        '    <TableFunction name="pressureDecline" inputVarNames="{time}"\n'
        '      coordinates="{ 0.0, 157680000.0 }"\n'
        '      values="{ 38332780.5, 28895581.5 }" interpolation="linear"/>',
        hydrostatic_functions(start_years, end_years),
        "pressure functions",
    )
    text = replace_once(text, 'newtonTol="1.0e-3"', f'newtonTol="{NEWTON_TOL}"', "Newton tolerance")
    text = replace_once(
        text,
        'maxEventDt="2.0e5"',
        f'maxEventDt="{MAX_EVENT_DT}"',
        "maximum event time step",
    )
    text = replace_once(text, ORIGINAL_KROW, THOMAS_DRAINAGE_KROW, "water/oil relative-permeability table")
    text = replace_once(text, '  <Events maxTime="1.5768e8">', f'  <Events maxTime="{end_time:.1f}">', "end time")
    text = text.replace('timeFrequency="7.884e6"', f'timeFrequency="{STEP_YEARS * YEAR:.1f}"')
    restart_event = (
        f'    <PeriodicEvent name="restartOutput" timeFrequency="{end_time:.1f}" '
        'targetExactTimestep="1" target="/Outputs/restartOutput"/>\n'
    )
    text = replace_once(text, "  </Events>", restart_event + "  </Events>", "restart event")
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
    header = (
        "<!-- Purpose: approximate Thomas Eq. (25) for near-connate-water gas/oil drainage. "
        "Key information: krow is held at 0.998-1.000 over the active oil-saturation range. -->\n"
        "<!-- Purpose: test Thomas fine-grid numerical sensitivity. "
        f"Key information: right-endpoint Pc, newtonTol={NEWTON_TOL}, maxEventDt={MAX_EVENT_DT} s. -->\n"
        "<!-- Purpose: validate Thomas 1983 with PVT-consistent hydrostatic fracture gas. "
        f"Key information: variant={variant}, Pc scale={scale:.8f}, interval={start_years:.1f}-{end_years:.1f} yr, "
        f"gas density={gas_density_at(start_years):.5f}-{gas_density_at(end_years):.5f} kg/m3. -->\n"
    )
    output_dir = ROOT / "decks" / variant
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{name}.xml").write_text(header + text, encoding="utf-8")


def main() -> None:
    """Generate the final right-endpoint chain and traceability manifest."""

    rows: list[dict[str, float | int | str]] = []
    segment_count = round(END_YEARS / STEP_YEARS)
    for index in range(segment_count):
        start_years = index * STEP_YEARS
        end_years = (index + 1) * STEP_YEARS
        name = segment_name(index, start_years, end_years)
        variant = "right"
        sample_years = end_years
        sigma = surface_tension_at(sample_years)
        scale = sigma / INITIAL_SURFACE_TENSION
        make_deck(name, scale, start_years, end_years, variant)
        rows.append(
            {
                "variant": variant,
                "segment": index,
                "start_years": start_years,
                "end_years": end_years,
                "pc_sample_years": sample_years,
                "pressure_psig": pressure_psig_at(sample_years),
                "surface_tension_dyn_per_cm": sigma,
                "pc_scale": scale,
                "gas_density_start_kg_per_m3": gas_density_at(start_years),
                "gas_density_end_kg_per_m3": gas_density_at(end_years),
                "deck": f"decks/{variant}/{name}.xml",
            }
        )

    with (ROOT / "segment_manifest.csv").open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write("# Purpose: trace the final Thomas Eq. (25) right-endpoint restart chain.\n")
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
