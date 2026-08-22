<!-- Purpose: record the algebraic finite-difference check of the capillary potential derivative. -->

# J0: `PotGrad` derivative contract

## Validation purpose

Check the local derivative signs used when the capillary potential is
differentiated with respect to matrix and fracture phase volume fractions.
This is a small algebraic contract test for the changed `PotGrad` logic.

## Validation target

For the reduced potential

```text
Phi(Sm, Sf) = p_m - (a + b Sm) - p_f + (a + b Sf),
```

the expected derivatives are

```text
dPhi/dSm = -b,
dPhi/dSf = +b.
```

The sign is the relevant target: increasing the matrix capillary pressure
reduces the matrix-side potential, while increasing the fracture capillary
pressure increases the potential.

## Why this can validate the target

The finite-difference test removes flash, mobility, gravity, mesh, and Newton
effects and evaluates the exact local algebra represented by the kernel. It is
therefore sensitive to a swapped sign or to assigning the matrix derivative
to the fracture support point. It is deliberately a unit-level check; it
cannot detect an incorrect field mapping or a runtime path that never calls
the derivative.

## Physical picture

J0 is a local potential-direction check rather than a runtime transport case.
Increasing matrix phase fraction raises matrix capillary pressure and lowers
the potential; increasing fracture phase fraction raises the potential. The
two curves therefore have negative and positive slopes, respectively.

## Reference and command

The independent reference calculation is in
[`scripts/check_J0.py`](scripts/check_J0.py).

```text
bash inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/J0_potgrad_derivative/scripts/run_J0.sh
```

The central difference uses `Sm=0.30`, `Sf=0.70`, `a=2000 Pa`, `b=8000 Pa`,
and `h=1e-4`.

## Numerical result

## Validation figure

- [PNG: local potential and finite-difference points](figures/J0_potgrad_derivative.png)
- [PDF: local potential and finite-difference points](figures/J0_potgrad_derivative.pdf)

Regenerate with `python3 scripts/plot_J0.py`.

图中直线是约化势函数的代数基准，橙色点是中心有限差分采样点。正确的
局部组装应使基质曲线斜率为 `-8000 Pa`、裂缝曲线斜率为 `+8000 Pa`；
本次有限差分结果正是 `-8000` 和 `+8000 Pa`，相对误差小于 `1e-10`。
这只能证明局部导数符号契约，不能证明运行时字段映射或 kernel 调用路径。

```text
dPhi/dSm = -8.000000000e+03 Pa
dPhi/dSf =  8.000000000e+03 Pa
```

Both relative errors are below `1e-10`.

The archived check output is [`runs/20260821/check_J0.txt`](runs/20260821/check_J0.txt).

## Conclusion and limitation

**J0 passes as an algebraic derivative check.** It supports the local signs
used by `PotGrad`, but it is not a direct GEOS runtime validation and cannot
override the unresolved dynamic result of I0.
