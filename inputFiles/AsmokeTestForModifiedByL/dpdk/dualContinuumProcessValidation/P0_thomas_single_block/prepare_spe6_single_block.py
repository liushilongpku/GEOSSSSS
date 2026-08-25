#!/usr/bin/env python3
"""Audit the 4500-psig SPE6 single-block initialization and Pc scaling."""

# Purpose: audit the P0 SPE6 live-oil composition and reproduce the Thomas IFT-based gas/oil Pc scale.
# Key information: 4500 psig is linearly interpolated from Thomas Table 2, while oil viscosity remains PVTO-derived.

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET


PRESSURE_PSIG = (1674.0, 2031.0, 2530.0, 2991.0, 3553.0, 4110.0, 4544.0, 4935.0, 5255.0, 5545.0)
IFT_DYN_PER_CM = (6.0, 4.7, 3.3, 2.2, 1.28, 0.72, 0.444, 0.255, 0.155, 0.090)
TARGET_PRESSURE_PSIG = 4500.0
INITIAL_COMPONENTS = (0.58112639, 0.13187000, 0.28700361)
BASE_PC_PA = (-5102.12, -3240.54, -2275.27, -413.69, 1447.90, 3309.48, 5171.07, 8756.34)
DECK = Path(__file__).with_name("P0_spe_singleblock_gas_oil_zeroPcf.xml")


def linear_interpolate(x: float, xs: tuple[float, ...], ys: tuple[float, ...]) -> float:
    for lower, upper in zip(range(len(xs) - 1), range(1, len(xs))):
        if xs[lower] <= x <= xs[upper]:
            weight = (x - xs[lower]) / (xs[upper] - xs[lower])
            return ys[lower] + weight * (ys[upper] - ys[lower])
    raise ValueError(f"pressure {x} psig is outside the Thomas table")


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
    sigma_4500 = linear_interpolate(TARGET_PRESSURE_PSIG, PRESSURE_PSIG, IFT_DYN_PER_CM)
    sigma_5545 = IFT_DYN_PER_CM[-1]
    scale = sigma_4500 / sigma_5545
    scaled_pc = tuple(value * scale for value in BASE_PC_PA)
    components = read_matrix_components()
    assert abs(sum(components) - 1.0) < 1.0e-10
    assert all(abs(a - b) < 1.0e-10 for a, b in zip(components, INITIAL_COMPONENTS))
    print(f"sigma(4500 psig) = {sigma_4500:.10f} dyn/cm")
    print(f"sigma(5545 psig) = {sigma_5545:.10f} dyn/cm")
    print(f"Pc scale = {scale:.10f}")
    print("scaled Thomas gas/oil Pc [Pa] = " + ", ".join(f"{value:.10f}" for value in scaled_pc))
    print("matrix global component fractions = " + ", ".join(f"{value:.8f}" for value in components))
    print("sum = " + f"{sum(components):.8f}")


if __name__ == "__main__":
    main()
