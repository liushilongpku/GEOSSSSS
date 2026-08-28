#!/usr/bin/env python3
"""Audit the Thomas 1983 dual-continuum single-block initialization and depletion.

Key information: 5545-psig live-oil matrix at the Thomas Table 4 state, a gas-filled
fracture continuum, and a 0.75-psi/day fracture depletion. The reference is the
Thomas 1983 Fig. 4 ultimate recovery of about 46% at 2.5 years for a 10-ft block.
"""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET


INITIAL_PRESSURE_PA = 38_332_780.5
PSI_TO_PA = 6_894.757293168
DEPLETION_PSI_PER_DAY = 0.75
YEARS_TO_SECONDS = 31_536_000.0
MAX_TIME_YEARS = 2.5
INITIAL_COMPONENTS = (0.53683289, 0.16590112, 0.29726599)
DECK = Path(__file__).with_name("thomas_singleblock_gas_oil_gravity_drainage.xml")


def read_matrix_components() -> tuple[float, float, float]:
    root = ET.parse(DECK).getroot()
    fields = root.find("FieldSpecifications")
    if fields is None:
        raise RuntimeError("FieldSpecifications is missing")
    values = []
    for name in ("matrixOilFrac", "matrixGasFrac", "matrixWaterFrac"):
        field = fields.find(f"FieldSpecification[@name='{name}']")
        if field is None:
            raise RuntimeError(f"missing {name}")
        values.append(float(field.attrib["scale"]))
    return tuple(values)


def main() -> None:
    components = read_matrix_components()
    assert abs(sum(components) - 1.0) < 1.0e-9
    assert all(abs(a - b) < 1.0e-9 for a, b in zip(components, INITIAL_COMPONENTS))

    end_time_s = MAX_TIME_YEARS * YEARS_TO_SECONDS
    end_pressure_pa = INITIAL_PRESSURE_PA - DEPLETION_PSI_PER_DAY * PSI_TO_PA * 365.0 * MAX_TIME_YEARS
    pressure_units_file = Path(__file__).parent / "reference" / "thomas_fig4_3d_model.csv"
    print(f"initial pressure = {INITIAL_PRESSURE_PA:.4f} Pa")
    print(f"  5545 psig -> Pa = {5545.0 * PSI_TO_PA:.4f} Pa (atm offset not included)")
    print(f"0.75 psi/day over {MAX_TIME_YEARS} years = {DEPLETION_PSI_PER_DAY * PSI_TO_PA * 365.0 * MAX_TIME_YEARS:.2f} Pa")
    print(f"end pressure = {end_pressure_pa:.2f} Pa at {end_time_s:.1f} s")
    print("matrix global component fractions = " + ", ".join(f"{value:.8f}" for value in components))
    print("sum = " + f"{sum(components):.8f}")
    print(f"Thomas Fig. 4 reference file = {pressure_units_file} (exists={pressure_units_file.is_file()})")


if __name__ == "__main__":
    main()
