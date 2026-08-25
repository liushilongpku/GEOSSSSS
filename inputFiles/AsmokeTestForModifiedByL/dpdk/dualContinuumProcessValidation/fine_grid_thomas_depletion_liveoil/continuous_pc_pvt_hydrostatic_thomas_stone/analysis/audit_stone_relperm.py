#!/usr/bin/env python3
"""Audit GEOS oil relative permeability against Thomas Eq. (25).

Key information: read the final right-endpoint restart without changing it,
then compare stored GEOS kro with direct Sw/Sg lookup of Thomas Table 3.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HDF5_PREFIX = (
    "Problem/domain/MeshBodies/mesh1/meshLevels/Level0/"
    "ElementRegions/elementRegionsGroup/matrixRegion/elementSubRegions"
)
FINAL_RUN = ROOT / "runs/right/segment_024_2p4_to_2p5yr"
FINAL_STEM = "segment_024_2p4_to_2p5yr"

SW = np.asarray([0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.75])
KRW = np.asarray([0.0, 0.005, 0.010, 0.020, 0.030, 0.045, 0.060, 0.110, 0.180, 0.230])
KROW = np.asarray([1.000, 0.860, 0.723, 0.600, 0.492, 0.392, 0.304, 0.154, 0.042, 0.000])
SG = np.asarray([0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55])
KRG = np.asarray([0.0, 0.015, 0.050, 0.103, 0.190, 0.310, 0.420])
KROG = np.asarray([1.000, 0.700, 0.450, 0.250, 0.110, 0.028, 0.000])


def values(group: h5py.Group, name: str) -> np.ndarray:
    """Read a GEOS wrapper stored directly or below __values__."""

    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[...], dtype=float)
    return np.asarray(obj["__values__"][...], dtype=float)


def final_rank_file() -> Path:
    """Resolve the final rank-zero restart file."""

    roots = sorted(FINAL_RUN.glob(f"{FINAL_STEM}_restart_*.root"))
    if not roots:
        raise FileNotFoundError(FINAL_RUN)
    return roots[-1].with_suffix("") / "rank_0000000.hdf5"


def interp(value: float, coordinates: np.ndarray, table: np.ndarray) -> float:
    """Linearly interpolate one Thomas table with endpoint clamping."""

    return float(np.interp(value, coordinates, table))


def thomas_kro(sw: float, sg: float) -> float:
    """Evaluate Thomas Eq. (25) using Table 3's native Sw and Sg axes."""

    krw = interp(sw, SW, KRW)
    krow = interp(sw, SW, KROW)
    krg = interp(sg, SG, KRG)
    krog = interp(sg, SG, KROG)
    return max(0.0, (krw + krow) * (krg + krog) - (krw + krg))


def main() -> None:
    """Write per-cell and height-averaged oil-relperm discrepancies."""

    rows: list[dict[str, float | str]] = []
    with h5py.File(final_rank_file(), "r") as h5_file:
        subregions = h5_file[HDF5_PREFIX]
        for name, subregion in subregions.items():
            if name == "__size__":
                continue
            phase = np.ravel(values(subregion, "phaseVolumeFraction"))[:3]
            mobility = np.ravel(values(subregion, "phaseMobility"))[:3]
            viscosity = np.ravel(values(subregion, "fluid/phaseViscosity"))[:3]
            mass_density = np.ravel(values(subregion, "fluid/phaseMassDensity"))[:3]
            relperm = mobility * viscosity / mass_density
            center = np.ravel(values(subregion, "elementCenter"))[:3]
            volume = float(np.ravel(values(subregion, "elementVolume"))[0])
            so, sg, sw = (float(value) for value in phase)
            geos_kro = float(relperm[0])
            expected_kro = thomas_kro(sw, sg)
            rows.append(
                {
                    "subregion": name,
                    "height_m": float(center[2]),
                    "volume_m3": volume,
                    "oil_saturation": so,
                    "gas_saturation": sg,
                    "water_saturation": sw,
                    "geos_oil_relperm": geos_kro,
                    "thomas_eq25_oil_relperm": expected_kro,
                    "geos_minus_thomas": geos_kro - expected_kro,
                    "geos_over_thomas": geos_kro / expected_kro if expected_kro > 0.0 else float("nan"),
                }
            )

    output = Path(__file__).resolve().parent / "stone_relperm_audit.csv"
    with output.open("w", newline="", encoding="utf-8") as csv_file:
        csv_file.write("# Purpose: compare stored GEOS kro with Thomas Table 3 and Eq. (25) at 2.5 years.\n")
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    geos = np.asarray([float(row["geos_oil_relperm"]) for row in rows])
    thomas = np.asarray([float(row["thomas_eq25_oil_relperm"]) for row in rows])
    volume = np.asarray([float(row["volume_m3"]) for row in rows])
    positive = thomas > 1.0e-12
    print(f"cells={len(rows)}")
    print(f"volume-weighted GEOS kro={np.average(geos, weights=volume):.8f}")
    print(f"volume-weighted Thomas Eq25 kro={np.average(thomas, weights=volume):.8f}")
    print(f"weighted GEOS/Thomas={np.average(geos[positive], weights=volume[positive]) / np.average(thomas[positive], weights=volume[positive]):.6f}")
    print(f"minimum cell GEOS/Thomas={np.min(geos[positive] / thomas[positive]):.6f}")


if __name__ == "__main__":
    main()
