# Thomas Single-Block Reproducibility

This directory reproduces the Thomas (1983) case 6 gas/oil gravity-drainage
calculation with one matrix cell and one colocated fracture cell. No grid
refinement, shape-factor adjustment, endpoint adjustment, or recovery-target
fitting is used.

## Canonical Cases

The canonical decks differ only in the explicit Thomas gravity-drainage pressure
(GDP) switch:

- `thomas_singleblock_gas_oil_gravity_drainage.xml`: GDP off.
- `thomas_singleblock_gas_oil_gravity_drainage_gdp_on.xml`: GDP on.

Both decks use the Thomas shape factor `0.215278208 m^-2`, the existing GEOS
modified Stone-II table interpolator, matrix-controlled PPU exchange, phase order
`{ oil, gas, water }`, and reverse-exchange relative
permeabilities `{ -1, 0.42, 0.03 }`. A negative oil entry selects the current
matrix oil relative permeability. Gas and water use the Thomas endpoint values.

Run from this directory:

```bash
python3 scripts/reproduce.py
python3 scripts/reproduce.py --run
```

The first command checks every local data dependency, regenerates and compares
the VE tables in a temporary directory, verifies the required GEOS source
extensions and Python packages, and asks GEOS to load both XML decks. The second
command additionally runs both 2.5-year cases, performs postprocessing, runs the
independent oracle, and asserts the 0.5- and 2.5-year recoveries. Outputs are
written under `/tmp/thomas_single_block_reproduction` by default.

The validation directory contains all case-specific inputs and postprocessing
code. It is not intended to duplicate GEOS source code. Reproduction therefore
requires the surrounding GEOS checkout to contain the implementations used by
the decks: `ThomasGasOilGravityDrainagePressure`, the `pressureDependentTableName` path in
`PressureScaledTableCapillaryPressure`, and matrix-controlled dual-continuum PPU
exchange. `scripts/reproduce.py` checks these source and executable dependencies.

With GEOS commit `3de23377b` and the current worktree:

| Case | 0.5-year recovery | 2.5-year recovery | Final `So/Sg/Sw` | Steps | Cuts |
| --- | ---: | ---: | --- | ---: | ---: |
| GDP off | 13.2576% | 18.7370% | 0.60449/0.19537/0.20014 | 3963 | 0 |
| GDP on | 27.3364% | 45.5423% | 0.40507/0.39479/0.20014 | 3963 | 0 |
| Thomas Fig. 4 | 27.0% | 46.0% | - | - | - |

Both GEOS runs have zero discarded nonlinear iterations. The GDP-on differences
from Thomas are `+0.3364 pp` at 0.5 years and `-0.4577 pp` at 2.5 years.

## Recovery Definition

Thomas reports stock-tank oil. The matching GEOS quantity is matrix oil-component
mass recovery:

```text
R = 1 - matrix_oil_component_mass(t) / matrix_oil_component_mass(0)
```

This is the black-oil inventory `phi * b_o * S_o = phi * S_o / B_o` up to a
constant surface density. `1-So/So0` is only a saturation diagnostic and gives
`49.4061%`, not the 2.5-year stock-tank recovery. A vertically resolved reference
case independently gives `46.0030%` by the mass definition.

## Pressure-Dependent Pseudocapillary Pressure

Thomas Eq. (28) scales the local rock curve from Table 3:

```text
Pc_go,local(Sg,P) = sigma(P) / sigma_I * Pc_go,local(Sg,P_I)
```

Table 4 is not that local curve. It is a block-average pseudofunction that already
contains the 10-ft vertical gravity integration. Multiplying the entire Table 4
curve by `sigma(P)/sigma_I` incorrectly scales its embedded gravity offset.

`scripts/generate_vertical_equilibrium_pseudocapillary.py` therefore scales Table 3 at
each Table 2 pressure, recomputes oil and gas density, and repeats the full-height
vertical-equilibrium integration. It writes the ordinary two-dimensional
`TableFunction` files:

- `tables/ve_pseudo_pc_sg_axis.txt`
- `tables/ve_pseudo_pc_pressure_axis.txt`
- `tables/ve_pseudo_pc_values.txt`

At 5545 psig the independent integration reproduces the six internal Table 4
points with a maximum error of `0.006699 psi`. The printed Table 4 values anchor
the reference-pressure curve; pressure changes come from the independent
vertical-equilibrium calculation, not recovery fitting.

`PressureScaledTableCapillaryPressure` reads this closure through
`pressureDependentTableName`. It bilinearly evaluates `Pc_go(Sg,P)` and
`dPc_go/dSg`, retains the water/oil capillary table in three-phase flow, and is
mutually exclusive with the older scalar `pressureScalingTableName` option. The
old scalar path remains unchanged when no two-dimensional table is configured.

## GDP and Oracle

The Thomas GDP oil head uses half the block height:

```text
GDP_o = (rho_o - rho_g) * abs(g_z) * Lz / 2
GDP_g = GDP_w = 0
```

The pressure-dependent Table 4 closure and this half-height exchange head are
both needed in the GEOS mapping. GDP-off reaches only `18.7370%`; using the full
height overdrives the model to `54.6939%`.

`scripts/thomas_single_cell_oracle.py` independently integrates the single-cell black-oil
inventory and exchange equations. With zero fracture capillary pressure it gives
`45.4514%` at 2.5 years. Across all 11 quarter-year output points, the maximum
absolute recovery difference between GEOS and the oracle is `0.0909 pp`. The
alternative endpoint fracture-capillary diagnostic gives `65.3323%` and is not
the canonical boundary condition.

## Remaining Transient Difference

The canonical GDP-on history does not match every digitized Fig. 4 point:

| Time / years | GEOS | Thomas | Difference / pp |
| ---: | ---: | ---: | ---: |
| 0.25 | 17.4678% | 12.5% | +4.9678 |
| 0.50 | 27.3364% | 27.0% | +0.3364 |
| 0.75 | 33.5632% | 35.25% | -1.6868 |
| 1.00 | 37.7124% | 40.0% | -2.2876 |
| 1.50 | 42.5042% | 44.15% | -1.6458 |
| 2.00 | 44.7305% | 46.0% | -1.2695 |
| 2.50 | 45.5423% | 46.0% | -0.4577 |

The digitized figure has about `+/-1 pp` reading uncertainty. The 0.25-year and
middle-time differences remain unexplained and must not be described as strict
full-curve agreement.

## Diagnostic Limits

Splitting both continua into `1x1x2`, `1x1x4`, or `1x1x8` gives 2.5-year
recoveries of `48.4454%`, `50.5551%`, and `51.0888%`. These are not convergence
results: every subcell repeats the full-block `fractureSpacingLz=3.05 m` closure.

Removing the non-initial fracture composition constraints leaves the old scalar
case unchanged, so fixed pure-gas composition is redundant here. The exchange
transmissibility is matrix-permeability controlled; its Jacobian contains
`dT/dP_matrix` and zero `dT/dP_fracture`.
