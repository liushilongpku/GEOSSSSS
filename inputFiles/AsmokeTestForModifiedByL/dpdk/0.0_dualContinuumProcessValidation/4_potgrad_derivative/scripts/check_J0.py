#!/usr/bin/env python3
# Purpose: check the local PotGrad capillary-potential derivative signs by finite difference.
"""Independent algebraic contract check for the J0 derivative signs."""

from __future__ import annotations


def main() -> None:
    a, b = 2.0e3, 8.0e3
    pm = pf = 1.0e6
    sm, sf, h = 0.30, 0.70, 1.0e-4

    def potential(s_m: float, s_f: float) -> float:
        return pm - (a + b * s_m) - pf + (a + b * s_f)

    dsm = (potential(sm + h, sf) - potential(sm - h, sf)) / (2.0 * h)
    dsf = (potential(sm, sf + h) - potential(sm, sf - h)) / (2.0 * h)
    error_m, error_f = abs(dsm + b) / b, abs(dsf - b) / b
    if error_m >= 1.0e-10 or error_f >= 1.0e-10:
        raise AssertionError(f"J0 criteria failed: dsm={dsm:g}, dsf={dsf:g}")
    print(f"J0 dPhi/dSm={dsm:.9e} Pa, relative_error={error_m:.3e}")
    print(f"J0 dPhi/dSf={dsf:.9e} Pa, relative_error={error_f:.3e}")


if __name__ == "__main__":
    main()
