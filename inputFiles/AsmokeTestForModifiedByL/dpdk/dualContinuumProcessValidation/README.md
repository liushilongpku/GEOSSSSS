<!-- Purpose: define the isolated validation package for dual-continuum pressure, gravity, and imbibition processes. -->

# Dual-continuum process validation

This directory contains the inputs, reference calculations, run commands, and
case-by-case reports for the three processes affected by the dual-continuum
cross-flow potential:

| Case | Process | What is isolated | Reference type |
| --- | --- | --- | --- |
| A0 | pressure-driven exchange | pressure contrast only; no gravity, capillary pressure, or mechanics | exact two-compartment solution |
| G0 | gravity-driven exchange | gravity-density contrast only; equal initial pressure | exact GDP sign and hydrostatic equilibrium |
| I0 | capillary imbibition | saturation-dependent capillary pressure; zero gravity | GEOS history, capillary-potential, and conservation checks |
| J0 | potential/Jacobian contract | pressure, saturation, gravity, and capillary terms | finite-difference algebraic check |
| P1 | coupled poromechanics pressure response | existing DPDP Mandel problem | published analytical solution and existing workflow |

The cases are deliberately separated. A passing result means only the
mechanism named in that row has passed the stated checks; it does not prove all
multiphase, gravity, and poromechanics combinations.

## Evidence standard

Each case report must contain:

1. validation purpose;
2. validation target;
3. why the case can validate that target;
4. governing equations and parameters;
5. exact or semi-analytic reference and pass criterion;
6. the executable command and the actual output location;
7. numerical result, error, and limitations.

Each case keeps its own reference/check/plot scripts. A0/G0/I0 histories and
the four P1 Mandel outputs are archived under their respective case-local
`runs/20260821/` directories.

The three I0 `.txt` files are legacy CO2-brine parameter inputs. Their parser
currently accepts model records only and does not accept comment lines, so
their purpose is documented by the I0 XML and report rather than by an
in-file comment that would make the input invalid.

The P1 package is self-contained: its relative Mandel loading tables,
analytical reference data, problem description, four input decks, run
artifacts, and comparison scripts are all under `P1_dpdp_mandel/`. The
`archive/` subdirectory preserves a duplicate run directory that was found
outside this package during the final audit; it is historical provenance and
is not searched by the active plotting workflow.

## Case status

| Case | Input present | Reference check | GEOS run | Report |
| --- | --- | --- | --- | --- |
| A0 | yes | equal- and unequal-volume pressure exchange passed | passed: GEOS exit code 0 | [综合中文报告](A0_pressure_exchange_combined/report_zh.md) |
| G0 | yes | passed | passed: GEOS exit code 0 | [report](G0_gravity_exchange/report.md) |
| I0 | yes | passed: dynamic history and conservation | passed: GEOS exit code 0, 10449 accepted steps, zero cuts | [report](I0_capillary_exchange/report.md) |
| J0 | yes | passed | not applicable | [report](J0_potgrad_derivative/report.md) |
| P1 | four inputs and analytical workflow | qualified quantitative comparison | passed runtime; qualified analytical agreement | [report](P1_dpdp_mandel/report.md) |

The status table is intentionally conservative and is updated only after a
command has been run and its output inspected.

The complete case-by-case explanation is in
[`summary.md`](summary.md). It states the purpose, target, justification,
reference, criterion, result, and limitations for every case.

Each case owns its input, scripts, run archive, and figures. Regenerate figures
independently with:

```bash
python3 A0_pressure_exchange/scripts/plot_A0.py
python3 G0_gravity_exchange/scripts/plot_G0.py
python3 I0_capillary_exchange/scripts/plot_I0.py
python3 J0_potgrad_derivative/scripts/plot_J0.py
bash P1_dpdp_mandel/scripts/plot_P1.sh
```

The independent checks are `A0_pressure_exchange/scripts/check_A0.py`,
`G0_gravity_exchange/scripts/check_G0.py`,
`I0_capillary_exchange/scripts/check_I0.py`, and
`bash J0_potgrad_derivative/scripts/run_J0.sh`. Their outputs stay in the
corresponding case directories.

Each command writes only to its own `figures/` directory. The A0, G0, and I0
histories used by the plots are archived under their respective
`runs/20260821/` directories; P1 keeps its four run trees and analytical
workflow under `P1_dpdp_mandel/`.
