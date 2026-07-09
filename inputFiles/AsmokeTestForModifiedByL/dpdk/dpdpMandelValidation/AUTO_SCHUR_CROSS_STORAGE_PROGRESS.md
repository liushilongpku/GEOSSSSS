# Auto Schur cross-storage calibration progress

Date: 2026-06-11

## Goal

Remove the empirical XML setting

```xml
crossStorageOffDiagScale="0.911"
```

from `DPDP_N2_dispdriven_fim_eff.xml` by computing the missing discrete mechanics Schur contribution from the assembled FIM matrix:

```text
K_pp,eff = K_pp - K_pu K_uu^{-1} K_up
```

The intended result was to replace a user supplied off-diagonal storage scale with an automatic correction based on the actual discrete mechanics block.

## Attempted implementation

The experimental code was placed in:

```text
src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp
```

The attempted path did the following:

1. Override `applyBoundaryConditions()` for the single-phase FIM dual-continuum poromechanics solver.
2. Apply mechanics boundary conditions first, so the displacement block `K_uu` includes Dirichlet constraints.
3. Before flow boundary conditions, compute a single-rank local Schur estimate:

   ```text
   mechanicsOff = K_pm,u K_uu^{-1} K_u,pf
   correctedOffdiag = targetOff + mechanicsOff
   ```

4. Add `deltaOffdiag = correctedOffdiag - assembledOffdiag` to the matrix-fracture pressure off-diagonal rows.

Two calibration modes were tested or partially tested:

- Single-cell fracture pressure impulse.
- Mandel-like drainage pressure mode using `cos(pi*x/(2L))`.

A third mode was started but not validated:

- Volume-weighted modal projection of both matrix and fracture pressure rows.

## Verified numerical result before rollback

The compiled Mandel-like mode ran and printed:

```text
poroSolver: auto cross-storage offdiag correction applied:
  assembled = -5.0877e-10
  corrected = -5.08671e-10
  delta     =  9.87441e-14
```

The desired value corresponding to the previously validated XML setting is approximately:

```text
-5.0877e-10 * 0.911 = -4.635e-10
```

Therefore the automatic Schur estimate did not reproduce the observed discrete correction. It effectively left the off-diagonal storage at the unscaled value.

## Current conclusion

The attempted implementation removes the explicit XML scale, but it still reduces the nonlocal Schur effect to one scalar `correctedOffdiag`. That is not a full removal of scaling; it is an automatic equivalent scalar calibration.

The tested scalar calibration is not good enough:

- It computes a near-continuum Schur contribution.
- It does not recover the observed Q1/Mandel discrete under-cancellation represented by `0.911`.
- It is single-rank only and uses a dense local `K_uu`, so it is not a production-quality distributed implementation.

The robust version should not be a single global off-diagonal scale. It should either:

1. Assemble a low-rank/modal approximation to `K_pu K_uu^{-1} K_up`, or
2. Build a proper distributed Schur correction/operator for the pressure block.

Until that is implemented, the validated `crossStorageOffDiagScale="0.911"` path should be kept for this verification deck.

## Restored baseline

After recording this progress, the code and main verification input were restored to the validated baseline:

- `SinglePhaseDualContinuum.cpp` still uses `DualContinuumCrossFlow::crossStorageOffDiagScale`.
- `DPDP_N2_dispdriven_fim_eff.xml` uses `crossStorageOffDiagScale="0.911"`.
- No automatic Schur correction remains active in `DualContinuumPoromechanicsSolverBase.hpp`.

