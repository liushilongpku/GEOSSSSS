# DPDP_N2 Mandel — oscillation diagnosis & fixes

Reference: Mehrabian & Abousleiman (2014), *Generalized Biot's theory and Mandel's
problem of multiple-porosity and multiple-permeability poroelasticity*, JGR Solid Earth.
(`problem_description/mandel_problem2014.md`)

## What was wrong (root causes of the oscillation)

`DPDP_N2_stressload.xml` ran to completion but the pressures oscillated in time
(fracture pressure swinging ±20 kPa, matrix pressure on a growing ~3 s wave).
Three compounding causes:

1. **Sequential coupling never iterated.**
   The deck used `sequentialConvergenceCriterion="SolutionIncrements"`. That
   criterion is computed by `FlowSolverBase::saveSequentialIterationState`, which
   the outer loop only calls when the *flow subsolver's own* `couplingType()` is
   `Sequential` (`CoupledSolver.hpp:534`). Here the flow solver is the **nested
   `dualFlow`** (a `FullyImplicit` `CoupledSolver`), so that call never fires:
   `m_sequentialPresChange` stays `0`, the log prints
   `Max pressure change during outer iteration: 0.000 Pa`, and every step
   "converges in 1 iteration". The poromechanics coupling was therefore a
   single-pass staggered scheme → temporal instability.

2. **`RigidBoundary` could not converge.**
   It applied the platen load as a plain Neumann force, then enforced the rigid
   constraint by a **post-solve projection** (`enforceConstraint` overwrote the
   top-face z-displacements with their average, inside `applySystemSolution`
   every Newton iteration). The projection knocks the iterate off equilibrium, so
   the solid residual plateaus (`8.55e+01`) and never converges — under
   FullyImplicit *or* properly-iterated Sequential. The original only "ran"
   because cause #1 accepted a single unconverged pass.

3. **Wrong time scale in the deck comment.**
   The header claimed `τ=3 → t=308412 s` (t0≈1.03e5 s). The correct
   characteristic time is `t0 = a²·ΣS/Σκ ≈ 4.0 s` (dominated by the fast
   fracture conductivity). With t0=4 s the simulation to 308412 s covers
   τ≈0…77000, i.e. the whole analytical curve.

There is also a **deeper instability in the custom dual-continuum poromechanical
coupling itself** (see "Remaining issue").

## Fix 1 — displacement-controlled verification (canonical Mandel method)

`DPDP_N2_dispdriven.xml` (new). Replaces the `RigidBoundary` load with a
**prescribed uniform top z-displacement** equal to the analytical N=2 platen
displacement `loadFunction0000000(t)` (already present in `DPDP_N2.xml` but left
unused). This is exactly how GEOS's own Mandel verification works
(`inputFiles/poromechanics/PoroElastic_Mandel_base.xml`): a Dirichlet BC realizes
the rigid platen exactly and converges trivially. Runs to t=308412 s with no
convergence failures and no dt cuts.

Comparison vs digitized Fig.5c (`analitical_result/GEOS_vs_Fig5c_dispdriven.png`),
using t0=4.0 s and normalizing each pressure by its Skempton t=0⁺ value
(p_m0=4.55e5, p_f0=4.88e5):

* **Fracture pressure: good qualitative match** — drains on the right time scale
  (τ≈0.1–1), ≈0 by τ≈3.
* **Matrix pressure: drainage timing matches (τ≈3000–30000) but the level is
  ~50% too high** (numerical Mandel-Cryer overshoot ≈1.38 vs analytical ≈1.06,
  plateau 1.38 vs 0.91). The matrix absorbs too much of the load shed by the
  draining fracture → the effective-Biot / stress-partitioning between continua
  is off (the area touched by the recent `K_eff` / effective-Biot commits), plus
  the early-time coupling oscillation.

## Fix 2 — RigidBoundary force-control made convergent (code)

`src/coreComponents/fieldSpecification/RigidBoundary.{hpp,cpp}` and
`SolidMechanicsLagrangianFEM.{hpp,cpp}`.

The post-solve projection is replaced by an in-Jacobian **face-local
graph-Laplacian penalty** (`RigidBoundary::applyRigidConstraint`): for each
boundary face, node component-DOFs are driven toward the face mean
(`r_a += -β(u_a - mean)`, Jacobian `-β(δ_ab − 1/n)`, β = 1e4·|diagonal|). The
uniform (collective) mode is in the penalty null space, so the platen is still
free to translate under the applied load, while within-face differences are
suppressed and shared nodes chain faces into global platen rigidity. All coupled
DOFs share an element, so the entries already exist in the sparsity pattern. The
convergence-breaking projection in `applySystemSolution` was removed.

Validation (pure elastic, isolates the BC from the poromechanics coupling):

| test | converges | top u_z uniformity |
|---|---|---|
| `RigidBoundary_elastic_test.xml` (homogeneous) | yes, R→1e-8 | exact (spread ~1e-20) |
| `RigidBoundary_hetero_test.xml` (100× stiffness contrast) | yes, R→1e-8 | flat platen, residual non-uniformity ≪1% |

(Old code: homogeneous converged trivially because the free solution is already
uniform; the heterogeneous/poro cases — where the projection actually moves the
solution — could not converge.)

## Remaining issue (not fixed here)

The dual-continuum poromechanical **coupling is itself unstable**, independent of
the loading method:
* Sequential + displacement BC converges per step (ResidualNorm) but shows a
  growing ~3 s temporal oscillation in t≈1–10 s that is highly dt-sensitive
  (at t≈2 s: dt=0.1 → p_f≈+1.2e5; dt=0.01 → p_f≈−2.6e4) — a genuine instability,
  not under-resolution.
* `couplingType="FullyImplicit"` **diverges at the first loaded step** (the solid
  residual drops once then plateaus) → the hand-written monolithic Jacobian in
  `DualContinuumPoromechanicsSolverBase::assembleSystem` is incomplete/incorrect.

These live in the custom fixed-stress porosity-update / monolithic-assembly code
(`DualContinuumPoromechanicsSolverBase.hpp`, the `K_eff` / effective-Biot
stress-increment correction). They are the next thing to fix to make the
force-driven `stressload` case match the analytical solution quantitatively.

## Files

* `DPDP_N2_dispdriven.xml` — working displacement-controlled verification (full curve).
* `RigidBoundary_elastic_test.xml`, `RigidBoundary_hetero_test.xml` — BC unit tests.
* `analitical_result/GEOS_vs_Fig5c_dispdriven.png` — numerical vs analytical Fig.5c.
* `DPDP_N2_stressload_fixed.xml` — stressload deck with `ResidualNorm` +
  more outer iterations (forces real sequential iteration; still limited by the
  coupling instability above).

---

## Fix 3 — dual-continuum coupling rewritten (shared volumetric strain)

After confirming the **single-porosity** path is correct (standard
`SinglePhasePoromechanics` matches the analytical Mandel decay to 3 digits via
`dpdkHMValidation/SP_Mandel_calibrated.xml`), the bug was isolated to the
dual-continuum coupling with an **identical-material test**: with matrix=fracture
the two continua diverged completely (p_m ballooned to 1.39 MPa and never
drained; p_f frozen) — physically impossible.

Root cause (`DualContinuumPoromechanicsSolverBase::mapSolutionBetweenSolvers`):
the matrix and fracture porosity (fixed-stress) updates were fed by a convoluted,
**asymmetric** stress-increment reconstruction (`ratio*S + ratio*α*Δp_eq − α*Δp_m`,
then `*=K_eff·v/K … /=`), reconstructing strain from the matrix stress by dividing
by `K_m` while the mechanics used `K_eff`. The two continua were therefore coupled
to the mechanics inconsistently.

Rewrite: compute the **shared** volumetric strain increment once,
`δε_v = (δσ_v + α_m·δp_eq)/K_eff`, then give each continuum its effective-Biot
strain coupling by setting `avgStress_i = v_i·K_eff·δε_v` (so GEOS' BiotPorosity
fixed-stress term `α_i·avgStress_i/K_i = ᾱ_i·δε_v`), using each continuum's own
restored `K_i`, `α_i`. Matrix and fracture are now treated symmetrically; the
K-swap round-trips are gone.

Validation:
* **Symmetric degenerate limit** (`DPDP_N2_symmetric_degradation.xml`, identical
  materials, v_m=v_f=0.5): `p_m == p_f` to machine precision (rel.diff ~1e-16).
  (Before: off by 3×.)
* **Full N=2** (`DPDP_N2_dispdriven.xml`): the early-time temporal **oscillation
  is eliminated** — fracture pressure now drains smoothly/monotonically
  (t=1→10 s: 424→3.6 kPa, no ±20 kPa swings). Matrix Mandel-Cryer overshoot
  dropped 1.38→1.28. See `analitical_result/GEOS_vs_Fig5c_fixedcoupling.png`.

Remaining quantitative gap (matrix plateau ≈1.28 vs analytical 0.91; fracture
drains slightly slow): most likely the **Reuss/series effective drained modulus**
`K_eff=(v_m/K_m+v_f/K_f)⁻¹` (softest estimate — probably not the paper's
Appendix-A modulus) and the within-/cross-continuum storage moduli `M̄_ij`. These
are modulus-calibration refinements, not stability bugs.

---

## Fix 4 — validated analytical reference; GEOS matrix error localized to the storage matrix

Using the full paper PDF (`problem_description/*.pdf`), I built and **validated** an
independent dual-porosity Mandel solver, `dpdkHMValidation/dpdp_mandel_analytical.py`
(mpmath, de Hoog inversion). The key bug was a **sign error** converting the
Appendix-A compliance `a_ij` to the eq.8 Biot modulus: the constant-strain modulus is
`1/Mbar_ij = a_ij - abar_i abar_j / Kbar` (MINUS, a Legendre transform), not `+`.
Single-porosity check: `alpha/(B K) - alpha^2/K = (alpha-phi)/Ks + phi/Kf` (the
standard 1/M). With the `+` sign the storage was ~16x too large, killing the
Mandel-Cryer effect.

After the fix the analytical **matches the digitized Fig5c** (matrix bumps to ~1.07
then settles ~0.89; fracture drains on the right timescale) — see
`analitical_result/analytical_corrected_vs_digitized.png`. The correct characteristic
time is **t0 = a^2*trace(S)/sum(kappa) = 10.5 s** (the broken script gave 4 s; the XML
comment said 1.03e5 s). The correct top displacement settles to u_z ~ 2.70e-5 m
(matching `mandel_input_tables/u_z_analitical.csv`; the XML's `loadFunction0000000`,
which only reached 2.21e-5, was the wrong one).

Conclusions for GEOS:
* Re-running GEOS disp-driven with the **correct** loadFunction
  (`DPDP_N2_dispdriven_correctLF.xml`) STILL gives matrix ~1.39 vs analytical/paper
  ~0.89 — so the matrix over-prediction is **not** the loadFunction; it is GEOS's
  dual-continuum **storage model**. See `analitical_result/GEOS_vs_analytical_vs_digitized.png`.
* GEOS uses per-continuum `BiotPorosity` with the intrinsic constant-strain modulus
  `1/M_i = (alpha_i-phi_i)/Ks + phi_i/Kf`. The correct multi-porosity storage matrix
  (Mehrabian eq.22) is `S_ij = 1/Mbar_ij + abar_i cm_j`, with
  `1/Mbar_ij = delta_ij v_i alpha_i/(B_i K_i) - abar_i abar_j/Kbar`. This has
  (a) volume-fraction weighting, (b) a `-abar_i^2/Kbar` diagonal reduction, and
  (c) a non-zero off-diagonal `-abar_1 abar_2/Kbar` coupling the matrix fluid content
  to the fracture pressure — none of which the per-continuum `BiotPorosity` reproduces.
  This is why GEOS keeps the matrix over-pressurized instead of letting it settle.

**Next GEOS fix (well-defined now):** give the dual-continuum flow the full effective
storage matrix `S_ij` (eq.22) — i.e. each continuum's accumulation must include the
v-weighting, the `-abar_i^2/Kbar` diagonal correction, and the off-diagonal
`-abar_i abar_j/Kbar * dp_j` cross term — rather than per-continuum intrinsic
`BiotPorosity`. Validate against `dpdp_mandel_analytical.py` (matrix -> 0.89 plateau).

---

## Fix 5 — analytical fully validated (pressure + stress); GEOS storage-matrix fix specified

`dpdkHMValidation/dpdp_mandel_analytical.py` now reproduces BOTH Fig5c (pressure)
and Fig5d (stress) of the paper. Two bugs were fixed:
1. `1/Mbar_ij = a_ij - abar_i abar_j/Kbar` (constant-strain Legendre transform; was `+`).
2. Center stress coupling coefficient is `2G*kappa_c/(1-2nu)` (eq.34), not `2G*kappa_c`
   (the extra `(1-2nu)` is correct only for the integrated stress BC that closes f);
   and Fig5d normalizes by `Pc`, not by sigma0. With both fixes the stress shows the
   correct two-hump-return-to-1.0 shape. See `analitical_result/analytical_vs_digitized_full.png`.
   The script `solve()`/`plot_vs_digitized()` sweep tau to 1e6.

GEOS storage-matrix fix (CONFIRMED by emulation): replacing the analytical's `1/Mbar`
with GEOS's per-continuum intrinsic `1/M_i = (alpha_i-phi_i)/Ks + phi_i/Kf` (diagonal,
no off-diagonal) reproduces GEOS's matrix plateau ~1.35-1.40 exactly; the correct
`1/Mbar` gives 0.89. BOTH the diagonal correction and the off-diagonal are required
(diagonal-only gives 1.83 — worse). The fix for the dual-continuum flow accumulation:

  zeta_i = abar_i * eps + sum_j (1/Mbar_ij) p_j ,   with
  1/Mbar_ij = delta_ij * v_i alpha_i/(B_i K_i)  -  abar_i abar_j / Kbar
  (abar_i = Kbar v_i alpha_i/K_i,  Kbar = (v1/K1+v2/K2)^-1,  B_i = Skempton)

i.e. each continuum's accumulation needs (a) diagonal storage `1/Mbar_ii` instead of
intrinsic `1/M_i`, and (b) an OFF-DIAGONAL cross-storage `-abar_i abar_j/Kbar * dp_j/dt`
coupling its fluid content to the OTHER continuum's pressure (an inter-continuum
storage term, analogous to but distinct from the Gamma cross-flow). The strain
coupling `abar_i*eps` is already handled (Fix 3). Acceptance test: matrix plateau ->
0.89 vs `dpdp_mandel_analytical.py`. Implementation point: the dual-continuum flow
accumulation/coupling assembly (alongside `DualContinuumCrossFlow`, which already
assembles the 2-dof inter-continuum Jacobian for Gamma).

---

## Fix 6 — multi-porosity storage matrix implemented in GEOS (matrix 1.39 -> 0.93)

Implemented the effective storage matrix M_bar in the dual-continuum flow accumulation:
- `DualContinuumCrossFlow`: new optional `fractureVolumeFraction` (v_f) input + accessor;
  `DualContinuumFlowSolverBase::getFractureVolumeFraction()` delegates to it.
- `SinglePhaseDualContinuum::assembleCouplingTerms`: when v_f>0, after the Gamma
  cross-flow it adds the cross-storage correction per matrix/fracture element pair:
    row_i += rho_i V [ (1/Mbar_ii - 1/M_i)(p_i-p_i_n) + (1/Mbar_ij)(p_j-p_j_n) ]
  with 1/Mbar_ii = v_i(1/M_i + alpha_i^2/K_i) - abar_i^2/Kbar, 1/Mbar_ij = -abar_i abar_j/Kbar,
  abar_i = Kbar v_i alpha_i/K_i, Kbar = (v_m/K_m + v_f/K_f)^-1, c_f from the fluid dDensity/density.
  Both the diagonal correction and the off-diagonal (matrix<->fracture pressure) Jacobian/
  residual contributions are assembled (entries already in the sparsity from the Gamma block).

Input: add `fractureVolumeFraction="0.03"` to the `<DualContinuumCrossFlow>` block
(done in DPDP_N2_dispdriven*.xml).

Validation (DPDP_N2_dispdriven_correctLF.xml, t0=10.5 s) vs analytical + digitized Fig5c:
  tau   matrix: GEOS / analytical / paper      fracture: GEOS / paper
  0.3   0.967 / 0.967 / 0.970                  0.470 / 0.302
  1     0.930 / 0.888 / 0.910                  0.008 / 0.012
  10-100 0.937-0.940 / ~0.89 / 0.91-0.92       ~0    / ~0
The matrix plateau is now ~0.93 (was 1.39, ~55% high) -> matches the paper's ~0.91.
See analitical_result/GEOS_crossStorage_vs_analytical_vs_digitized.png.

Residual minor gaps (secondary): the early Mandel-Cryer bump is slightly damped in GEOS
(1.01 vs 1.07) and the late matrix drainage (tau>300) lags a little -- attributable to the
flux volume-weighting (kappa_i = v_i k_i) not yet applied (cancels in the timescale but
shifts the tail). The dominant storage error is resolved.
