#!/usr/bin/env python3
"""Generate Thomas pressure-dependent gas/oil VE pseudo-capillary pressure.

Thomas Eq. (28) scales the local Table 3 rock capillary curve. The 10-ft
vertical-equilibrium integration is then repeated at each pressure; gravity is
not scaled with surface tension. The printed Table 4 values anchor the initial
pressure, while pressure changes come only from the independent VE calculation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


PSI_TO_PA = 6894.757293168
ATM_PA = 101325.0
G = 9.80665
HEIGHT_M = 3.048
OIL_SURFACE_DENSITY = 819.18
GAS_SURFACE_DENSITY = 0.929

PRESSURE_PSIG = np.asarray([1674, 2031, 2530, 2991, 3553, 4110, 4544, 4935, 5255, 5545, 7000], float)
BO = np.asarray([1.3001, 1.3359, 1.3891, 1.4425, 1.5141, 1.5938, 1.6630, 1.7315, 1.7953, 1.8540, 2.1978])
BG_RB_SCF = np.asarray([0.00198, 0.00162, 0.00130, 0.00111, 0.000959, 0.000855, 0.000795, 0.000751, 0.000720, 0.000696, 0.000600])
RS_SCF_STB = np.asarray([367, 447, 564, 679, 832, 1000, 1143, 1285, 1413, 1530, 2259], float)
SIGMA = np.asarray([6.0, 4.7, 3.3, 2.2, 1.28, 0.72, 0.444, 0.255, 0.155, 0.090, 0.050])

LOCAL_SG = np.asarray([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55])
LOCAL_PC_PSI = np.asarray([0.075, 0.085, 0.095, 0.115, 0.145, 0.255, 0.386])
PSEUDO_SG = np.asarray([0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55])
TABLE4_PC_PSI = np.asarray([-0.74, -0.47, -0.33, -0.06, 0.21, 0.48, 0.75, 1.27])
REFERENCE_PRESSURE_PSIG = 5545.0


def interp(values: np.ndarray, pressure_psig: float) -> float:
    return float(np.interp(pressure_psig, PRESSURE_PSIG, values))


def phase_densities(pressure_psig: float) -> tuple[float, float]:
    bo = interp(BO, pressure_psig)
    bg_m3_sm3 = interp(BG_RB_SCF, pressure_psig) * 5.614583333333333
    rs_sm3_sm3 = interp(RS_SCF_STB, pressure_psig) * 0.17810760667903525
    oil = (OIL_SURFACE_DENSITY + rs_sm3_sm3 * GAS_SURFACE_DENSITY) / bo
    gas = GAS_SURFACE_DENSITY / bg_m3_sm3
    return oil, gas


def ve_pc(target_sg: float, pressure_psig: float, integration_points: int = 20001) -> float:
    oil_density, gas_density = phase_densities(pressure_psig)
    head_psi = (oil_density - gas_density) * G * HEIGHT_M / (2.0 * PSI_TO_PA)
    scale = interp(SIGMA, pressure_psig) / interp(SIGMA, REFERENCE_PRESSURE_PSIG)
    zeta = np.linspace(-1.0, 1.0, integration_points)

    def average_sg(reference_pc: float) -> float:
        local_pc = (reference_pc + head_psi * zeta) / scale
        saturation = np.interp(local_pc, LOCAL_PC_PSI, LOCAL_SG)
        return float(np.trapezoid(saturation, zeta) / 2.0)

    lower = scale * LOCAL_PC_PSI[0] - head_psi
    upper = scale * LOCAL_PC_PSI[-1] + head_psi
    if target_sg <= LOCAL_SG[0]:
        return lower
    if target_sg >= LOCAL_SG[-1]:
        return upper
    for _ in range(80):
        middle = 0.5 * (lower + upper)
        if average_sg(middle) < target_sg:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def corrected_curve(pressure_psig: float, reference_ve: np.ndarray) -> np.ndarray:
    current_ve = np.asarray([ve_pc(sg, pressure_psig) for sg in PSEUDO_SG])
    return TABLE4_PC_PSI + current_ve - reference_ve


def write_values(path: Path, values: np.ndarray) -> None:
    path.write_text("\n".join(f"{value:.12g}" for value in values) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("tables"))
    args = parser.parse_args()

    pressure_axis_psig = PRESSURE_PSIG[(PRESSURE_PSIG >= 3553.0) & (PRESSURE_PSIG <= REFERENCE_PRESSURE_PSIG)]
    pressure_axis_pa = pressure_axis_psig * PSI_TO_PA + ATM_PA
    reference_ve = np.asarray([ve_pc(sg, REFERENCE_PRESSURE_PSIG) for sg in PSEUDO_SG])
    interior_error = np.max(np.abs(reference_ve[1:-1] - TABLE4_PC_PSI[1:-1]))
    if interior_error > 0.01:
        raise RuntimeError(f"VE reconstruction does not reproduce Table 4: max error={interior_error:.6f} psi")

    # TableFunction stores axis 0 (Sg) as the fastest-varying dimension.
    values_pa = np.concatenate(
        [corrected_curve(pressure, reference_ve) * PSI_TO_PA for pressure in pressure_axis_psig]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_values(args.output_dir / "ve_pseudo_pc_sg_axis.txt", PSEUDO_SG)
    write_values(args.output_dir / "ve_pseudo_pc_pressure_axis.txt", pressure_axis_pa)
    write_values(args.output_dir / "ve_pseudo_pc_values.txt", values_pa)
    print(f"Table 4 interior VE error = {interior_error:.6f} psi")
    print(f"wrote {len(pressure_axis_pa)}x{len(PSEUDO_SG)} Pc(Sg,P) table to {args.output_dir}")


if __name__ == "__main__":
    main()
