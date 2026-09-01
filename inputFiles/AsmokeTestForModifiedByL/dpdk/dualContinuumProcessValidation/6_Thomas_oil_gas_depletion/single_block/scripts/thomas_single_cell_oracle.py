#!/usr/bin/env python3
"""Independent black-oil oracle for the Thomas 1983 10-ft single-cell case.

The implementation follows Thomas Eqs. (4)-(7), (10), (23)-(26), and (28).
It intentionally does not import GEOS output or reuse GEOS constitutive helpers.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from generate_vertical_equilibrium_pseudocapillary import (
    REFERENCE_PRESSURE_PSIG,
    corrected_curve,
    phase_densities,
    ve_pc,
)


PSI_TO_PA = 6894.757293168
FT3_TO_M3 = 0.028316846592
DAY = 86400.0
YEAR = 365.0 * DAY

INITIAL_PRESSURE_PSIG = 5545.0
PRESSURE_DECLINE_PSI_PER_DAY = 0.75
BLOCK_VOLUME_FT3 = 10.0**3
POROSITY = 0.30
SHAPE_FACTOR_FT2 = 0.02
PERMEABILITY_MD = 1.0
SW = 0.20

PRESSURE_PSIG = np.asarray([1674, 2031, 2530, 2991, 3553, 4110, 4544, 4935, 5255, 5545, 7000], dtype=float)
BO = np.asarray([1.3001, 1.3359, 1.3891, 1.4425, 1.5141, 1.5938, 1.6630, 1.7315, 1.7953, 1.8540, 2.1978])
BG = np.asarray([0.00198, 0.00162, 0.00130, 0.00111, 0.000959, 0.000855, 0.000795, 0.000751, 0.000720, 0.000696, 0.000600])
RS = np.asarray([367, 447, 564, 679, 832, 1000, 1143, 1285, 1413, 1530, 2259], dtype=float)
MU_O = np.asarray([0.529, 0.487, 0.436, 0.397, 0.351, 0.310, 0.278, 0.246, 0.229, 0.210, 0.109])
MU_G = np.asarray([0.0162, 0.0171, 0.0184, 0.0197, 0.0213, 0.0230, 0.0244, 0.0255, 0.0265, 0.0274, 0.0330])
SIGMA = np.asarray([6.0, 4.7, 3.3, 2.2, 1.28, 0.72, 0.444, 0.255, 0.155, 0.090, 0.050])

SG_TABLE = np.asarray([0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55])
PCGO_PSI = np.asarray([-0.74, -0.47, -0.33, -0.06, 0.21, 0.48, 0.75, 1.27])
KRG = np.asarray([0.0, 0.026, 0.062, 0.14, 0.21, 0.29, 0.37, 0.42])
KRO = np.asarray([1.0, 0.89, 0.80, 0.62, 0.44, 0.25, 0.070, 0.0])

VE_PRESSURE_PSIG = PRESSURE_PSIG[(PRESSURE_PSIG >= 3553.0) & (PRESSURE_PSIG <= REFERENCE_PRESSURE_PSIG)]
VE_REFERENCE = np.asarray([ve_pc(value, REFERENCE_PRESSURE_PSIG) for value in SG_TABLE])
VE_PC_TABLE = np.asarray([corrected_curve(pressure, VE_REFERENCE) for pressure in VE_PRESSURE_PSIG])

# 0.001127 converts md*ft*psi/cP to reservoir bbl/day. Here sigma*V has units ft.
TRANSMISSIBILITY = 0.001127 * SHAPE_FACTOR_FT2 * BLOCK_VOLUME_FT3 * PERMEABILITY_MD
PV_RB = POROSITY * BLOCK_VOLUME_FT3 / 5.614583333333333


def interp(values: np.ndarray, pressure_psig: float) -> float:
    return float(np.interp(pressure_psig, PRESSURE_PSIG, values))


def pseudo_pcgo(sg: float, pressure_psig: float) -> float:
    values_at_sg = np.asarray([np.interp(sg, SG_TABLE, curve) for curve in VE_PC_TABLE])
    return float(np.interp(pressure_psig, VE_PRESSURE_PSIG, values_at_sg))


def gravity_head_psi(pressure_psig: float) -> float:
    oil_density, gas_density = phase_densities(pressure_psig)
    return (oil_density - gas_density) * 9.80665 * 3.048 / (2.0 * PSI_TO_PA)


def matrix_inventories(po_psig: float, sg: float) -> np.ndarray:
    so = 1.0 - SW - sg
    pg_psig = po_psig + pseudo_pcgo(sg, po_psig)
    oil_stb = PV_RB * so / interp(BO, po_psig)
    gas_scf = PV_RB * sg / interp(BG, pg_psig) + interp(RS, po_psig) * oil_stb
    return np.asarray([oil_stb, gas_scf])


def inventory_jacobian(po_psig: float, sg: float) -> np.ndarray:
    dp = 0.01
    ds = 1.0e-7
    return np.column_stack(
        (
            (matrix_inventories(po_psig + dp, sg) - matrix_inventories(po_psig - dp, sg)) / (2.0 * dp),
            (matrix_inventories(po_psig, sg + ds) - matrix_inventories(po_psig, sg - ds)) / (2.0 * ds),
        )
    )


def fracture_gas_pressure(time_days: float) -> float:
    return INITIAL_PRESSURE_PSIG - PRESSURE_DECLINE_PSI_PER_DAY * time_days


def rates(time_days: float, po_psig: float, sg: float, fracture_pc_mode: str) -> tuple[float, float, float, float]:
    pg_matrix = po_psig + pseudo_pcgo(sg, po_psig)
    pg_fracture = fracture_gas_pressure(time_days)
    if fracture_pc_mode == "zero":
        po_fracture = pg_fracture
    elif fracture_pc_mode == "endpoint":
        po_fracture = pg_fracture - pseudo_pcgo(SG_TABLE[-1], pg_fracture)
    else:
        raise ValueError(fracture_pc_mode)

    delta_po = po_psig - po_fracture + gravity_head_psi(po_psig)
    delta_pg = pg_matrix - pg_fracture
    kro = float(np.interp(sg, SG_TABLE, KRO))
    krg_matrix = float(np.interp(sg, SG_TABLE, KRG))

    # Positive rates are matrix-to-fracture. A pure-gas fracture has zero oil coverage.
    oil_rate_stb_day = 0.0
    if delta_po > 0.0:
        oil_rate_stb_day = TRANSMISSIBILITY * kro / interp(MU_O, po_psig) * delta_po / interp(BO, po_psig)

    if delta_pg >= 0.0:
        gas_rate_scf_day = (
            TRANSMISSIBILITY * krg_matrix / interp(MU_G, pg_matrix) * delta_pg / interp(BG, pg_matrix)
        )
    else:
        # Thomas Eq. (26): endpoint matrix krg times fracture gas coverage (one here).
        gas_rate_scf_day = (
            TRANSMISSIBILITY * KRG[-1] / interp(MU_G, pg_fracture) * delta_pg / interp(BG, pg_fracture)
        )

    total_gas_rate_scf_day = gas_rate_scf_day + interp(RS, po_psig) * oil_rate_stb_day
    return oil_rate_stb_day, total_gas_rate_scf_day, delta_po, delta_pg


def solve(fracture_pc_mode: str) -> list[dict[str, float | str]]:
    initial_state = np.asarray([INITIAL_PRESSURE_PSIG, 0.0])

    def rhs(time_days: float, state: np.ndarray) -> np.ndarray:
        po_psig, sg = state
        oil_rate, gas_rate, _, _ = rates(time_days, po_psig, sg, fracture_pc_mode)
        inventory_rate = np.asarray([-oil_rate, -gas_rate])
        return np.linalg.solve(inventory_jacobian(po_psig, sg), inventory_rate)

    sample_days = np.linspace(0.0, 2.5 * 365.0, 11)
    result = solve_ivp(
        rhs,
        (sample_days[0], sample_days[-1]),
        initial_state,
        t_eval=sample_days,
        method="BDF",
        rtol=1.0e-9,
        atol=(1.0e-5, 1.0e-10),
        max_step=1.0,
    )
    if not result.success:
        raise RuntimeError(result.message)

    initial_oil = matrix_inventories(*initial_state)[0]
    rows: list[dict[str, float | str]] = []
    for time_days, (po_psig, sg) in zip(result.t, result.y.T):
        oil, gas = matrix_inventories(po_psig, sg)
        oil_rate, gas_rate, delta_po, delta_pg = rates(time_days, po_psig, sg, fracture_pc_mode)
        rows.append(
            {
                "fracture_pc_mode": fracture_pc_mode,
                "time_years": time_days / 365.0,
                "matrix_oil_pressure_psig": po_psig,
                "matrix_gas_saturation": sg,
                "matrix_oil_inventory_stb": oil,
                "matrix_total_gas_inventory_scf": gas,
                "oil_recovery_pct": 100.0 * (1.0 - oil / initial_oil),
                "oil_potential_difference_psi": delta_po,
                "gas_potential_difference_psi": delta_pg,
                "oil_rate_stb_day": oil_rate,
                "total_gas_rate_scf_day": gas_rate,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/thomas_single_cell_oracle.csv"))
    args = parser.parse_args()

    rows = solve("zero") + solve("endpoint")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for mode in ("zero", "endpoint"):
        final = [row for row in rows if row["fracture_pc_mode"] == mode][-1]
        print(
            f"fracture Pc={mode}: recovery={final['oil_recovery_pct']:.6f}%, "
            f"Sg={final['matrix_gas_saturation']:.6f}, "
            f"deltaPo={final['oil_potential_difference_psi']:.6f} psi"
        )


if __name__ == "__main__":
    main()
