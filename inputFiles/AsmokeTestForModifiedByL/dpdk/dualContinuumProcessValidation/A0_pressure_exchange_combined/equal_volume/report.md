<!-- Purpose: record the completed A0 pressure-exchange validation and its reproducible evidence. -->

# A0: pressure-driven dual-continuum exchange

## Validation purpose

连续时间参考的完整推导见
[`continuous_reference_derivation.md`](continuous_reference_derivation.md)。

Verify the most elementary pressure-driven exchange process in the dual-
continuum solver before adding gravity, capillary pressure, or mechanics.
This isolates the cross-flow residual sign, the equal-and-opposite matrix and
fracture contributions, and the time integration of the exchange term.

## Validation target

For two closed, co-located, single-phase cells, verify that:

- the pressure contrast decays monotonically toward zero;
- exchange removes mass from the high-pressure continuum and adds it to the
  low-pressure continuum;
- the stored pressure/volume average is conserved;
- the measured backward-Euler decay factor agrees with the coefficient implied
  by the input and the dual-continuum shape-factor calculation.

## Physical picture

The high-pressure matrix drains into the low-pressure fracture through the
cross-flow connection. Matrix pressure falls and fracture pressure rises;
because the system is closed, both approach the common `1.5 MPa` equilibrium.

This case does **not** validate capillary pressure, gravity drainage, or the
changed compositional `PotGrad` kernel. Those mechanisms are covered by I0,
G0, and J0/P1 respectively.

## Why this can validate the target

There is one matrix cell and one fracture cell, equal volume and equal fluid
properties, no external boundary, no gravity, and no capillary model. The only
nonzero inter-continuum term is therefore the pressure exchange term. The
continuum equations reduce to

\[
C_m \frac{dp_m}{dt}=-A(p_m-p_f), \qquad
C_f \frac{dp_f}{dt}=+A(p_m-p_f).
\]

For equal storage, the exact continuous solution has

\[
p_m-p_f=\Delta p_0\exp[-2A t/C],
\qquad p_m+p_f=\text{constant}.
\]

GEOS advances this linearized system with backward Euler, so for `dt=0.02 s`
the expected ratio is

\[
q=\frac{\Delta p^{n+1}}{\Delta p^n}
 =\frac{1}{1+2A\,dt/C}.
\]

For this deck:

- `V=1 m^3`, `k=1e-12 m^2`, and the three shape-factor terms give
  `W=4V/Lx^2+4V/Ly^2+4V/Lz^2=1.2e-3`;
- `T=kW=1.2e-15 m^3`;
- single-phase mobility is `rho/mu=1.0e6`, hence `A=T*rho/mu=1.2e-9`;
- `C=phi*rho*(c_fluid+c_porosity)=0.2*1000*(1e-12+1e-12)=4e-10`;
- therefore `q=1/(1+2*(1.2e-9)*0.02/(4e-10))=0.8928571429`.

## Input and command

Input: [`A0_pressure_exchange.xml`](A0_pressure_exchange.xml)

Command used:

```text
timeout 180s build-ubuntu-lsl-release/bin/geosx \
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/A0_pressure_exchange/A0_pressure_exchange.xml \
  -o inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/A0_pressure_exchange/runs/20260821
```

GEOS completed with exit code `0`, 100 time steps, zero time-step cuts, and
100 successful nonlinear/linear iterations reported by `dualFlow`. The run
reached `t=2.0 s`; the archived periodic history contains samples through
`t=1.98 s` because of the final event ordering.

The archived time-history files are:

- [`runs/20260821/A0_matrix_pressure_history.hdf5`](runs/20260821/A0_matrix_pressure_history.hdf5)
- [`runs/20260821/A0_fracture_pressure_history.hdf5`](runs/20260821/A0_fracture_pressure_history.hdf5)

The first recorded row is the post-initial-condition state labelled `t=0` by
the output event; it is not interpreted as the raw field-specification state.

## Numerical result

The resolved pressure curves now cover the two-second run with `dt=0.02 s`
and archived samples through `t=1.98 s`.
The first-step decay factor is `0.8928571053`, compared with the backward-Euler
reference `0.8928571429`, for a relative error of `4.21e-8`. The maximum
relative error of the pressure-contrast curve while the reference remains
above `100` times the observed GEOS subtraction floor is `4.39e-7` (through
`t=1.16 s`). The pressure contrast decreases from about `8.93e5 Pa` to
`11.97 Pa` in the archived history, while the mean-pressure drift reaches `0.132 Pa`, or
`8.81e-8` relative to `1.5 MPa`.

## Validation figure

- [PNG: pressure exchange and backward-Euler reference](figures/A0_pressure_exchange.png)
- [PDF: pressure exchange and backward-Euler reference](figures/A0_pressure_exchange.pdf)

Regenerate with `python3 scripts/plot_A0.py`; check with `python3 scripts/check_A0.py`.

The figure contains one pressure-versus-time panel. It overlays the GEOS
matrix/fracture histories with both the backward-Euler pressure references and
the continuous exponential references (`tau=0.167 s`); the simple reference
formulas are shown inside the panel. Correct assembly is indicated by matrix pressure
falling while fracture pressure rises, both histories following their
references over the resolved interval, and near-conservation of mean pressure.

| quantity | result | criterion | status |
| --- | ---: | ---: | --- |
| `q` from first two stored rows | `0.8928571053` | `0.8928571429` | pass |
| relative `q` error | `4.21e-8` | `< 1e-5` | pass |
| maximum mean-pressure drift | `0.132 Pa` | relative `< 1e-6` | pass |
| relative mean-pressure drift | `8.81e-8` | `< 1e-6` | pass |
| resolved-curve maximum relative error through `t=1.16 s` | `4.39e-7` | `< 1e-4` | pass |
| last archived pressure contrast at `t=1.98 s` | `11.97 Pa` | follows reference | diagnostic |
| first five contrast magnitudes | `8.9286e5, 7.9719e5, 7.1178e5, 6.3552e5, 5.6743e5 Pa` | decreasing | pass |

The small `q` discrepancy and the slow mean-pressure drift are consistent
with the weak nonlinear compressible fluid/porosity update and the fact that
the reference uses a linearized storage coefficient. The shorter time step
resolves the physical transient; the final pressure-difference floor is still
a subtraction-resolution diagnostic rather than a physical equilibrium error.

## Conclusion and limitation

**A0 passes in the numerically resolved interval.** The isolated
dual-continuum pressure exchange has the expected direction, equal-and-
opposite transfer, backward-Euler decay, and conservation until the pressure
difference reaches the floating-point subtraction floor. It is a prerequisite
test for the later process cases, but by itself it says nothing about the sign
of the capillary or gravity additions in the compositional `PotGrad`
implementation.
