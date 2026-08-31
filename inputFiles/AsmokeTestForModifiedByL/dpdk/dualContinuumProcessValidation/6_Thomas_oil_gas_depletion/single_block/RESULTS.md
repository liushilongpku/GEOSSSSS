# Canonical Validation Results

The checked-in result artifacts in this directory were generated from the
strict `1x1x1` matrix and `1x1x1` fracture canonical cases using the
pressure-dependent vertical-equilibrium pseudocapillary table.

## Key Results

| Case | 0.5-year recovery | 2.5-year recovery | Final `So/Sg/Sw` |
| --- | ---: | ---: | --- |
| GEOS, Thomas GDP on | 27.3364% | 45.5423% | 0.40507/0.39479/0.20014 |
| GEOS, GDP off | 13.2576% | 18.7370% | 0.60449/0.19537/0.20014 |
| Independent oracle, fracture `Pc=0` | 27.3141% | 45.4514% | `Sg=0.39554` at 2.5 years |
| Thomas Fig. 4 | 27.0% | 46.0% | - |

The GDP-on GEOS run used 3,963 time steps, with zero time-step cuts and zero
discarded nonlinear iterations. Its 2.5-year recovery differs from Thomas by
`-0.4577 pp`. The maximum GEOS-oracle recovery difference across all 11 stored
times is `0.0909 pp`.

## Result Files

- `analysis/recovery_results.csv`: canonical GDP-on GEOS history and Thomas comparison.
- `analysis/gdp_off/recovery_results.csv`: canonical GDP-off comparison.
- `analysis/thomas_single_cell_oracle_ve.csv`: independent oracle histories for
  zero and endpoint fracture capillary pressure.
- `analysis/geos_oracle_recovery_comparison.csv`: pointwise GDP-on GEOS versus
  zero-fracture-`Pc` oracle comparison.
- `figures/recovery_comparison.png`: canonical GDP-on recovery plot.
- `figures/recovery_comparison.pdf`: vector version of the canonical plot.
- `figures/gdp_off/recovery_comparison.png`: GDP-off diagnostic plot.
- `figures/gdp_off/recovery_comparison.pdf`: vector version of the GDP-off plot.

## Key Inputs And Implementations

- `thomas_singleblock_gas_oil_gravity_drainage_gdp_on.xml`: canonical GDP-on deck.
- `thomas_singleblock_gas_oil_gravity_drainage.xml`: canonical GDP-off deck.
- `generate_vertical_equilibrium_pseudocapillary.py`: rebuilds `Pcgo(Sg,P)` from
  Thomas Tables 2 and 3 and checks the initial-pressure Table 4 reconstruction.
- `tables/ve_pseudo_pc_sg_axis.txt`: gas-saturation axis.
- `tables/ve_pseudo_pc_pressure_axis.txt`: absolute-pressure axis.
- `tables/ve_pseudo_pc_values.txt`: two-dimensional pseudocapillary values.
- Matrix three-phase relative permeability uses the existing GEOS
  `TableRelativePermeability` model with the modified Stone-II interpolator.
- `thomas_single_cell_oracle.py`: independent single-cell black-oil oracle.
- `analyze_results.py`: computes oil-component mass recovery and plots Fig. 4 comparison.
- `reproduce.py`: checks all case and repository dependencies; with `--run`,
  executes both canonical cases and verifies their key recoveries.
- `REPRODUCIBILITY.md`: detailed derivation and reproduction procedure.
- `report_zh.md`: Chinese technical report and limitations.

The complete GEOS simulation directories remain external under `/tmp`; only the
compact, reproducible result artifacts are stored here.
