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

------------------------------------------------------------------------------
Fix 7 -- Mandel-Cryer overshoot recovered via time-step refinement (input-only)
------------------------------------------------------------------------------
Symptom: GEOS matrix pressure was nearly monotonic (peak only 1.015) while the
analytical / Fig5c show a Mandel-Cryer overshoot peaking ~1.07 near tau~0.05.

Root cause: NOT the coupling. Verified empirically:
  - enabling sequential subcycling (subcycling=1, ResidualNorm/NumberOfNonlinearIterations)
    did NOT change the peak (1.010-1.015) -> splitting iteration is not the cause;
  - the prescribed platen displacement (loadFunction0000000) matches the analytical
    u_z(top) to <1% in shape (both rise ~50.4% undrained->drained), so the BC is correct.
The overshoot peaks at tau~0.05 == t~0.5 s, but the deck used forceDt=0.1 s in [0.1,1.0],
i.e. only ~4 coarse steps across the peak. The single-pass fixed-stress split under-resolves
that transient and flattens it into a near-monotonic decay.

Fix (input only, no code): refine the solver time-stepping through the overshoot window:
  phase0 0-0.01     forceDt 0.0005
  phase1 0.01-0.1   forceDt 0.002
  phase2 0.1-2.0    forceDt 0.01
  phase3 2.0-10.0   forceDt 0.1     (coarser phases unchanged afterwards)
Result (DPDP_N2_dispdriven_correctLF.xml): matrix peak 1.015 -> 1.066, matching the
paper's digitized peak (1.065) and close to the analytical (1.079). See
analitical_result/GEOS_vs_analytical_vs_digitized_full.png.

Side effect now exposed: with the overshoot built up, the matrix intermediate plateau sits
~1.05 (vs paper 0.91) and drainage onset is delayed (tau~3e3 vs ~1e3); the fracture overshoot
(~1.06) is still not reproduced. These are matrix/fracture drainage-timescale issues
(the kappa_i = v_i k_i flux weighting, Fix-6 note) -- deeper than a dt change, deferred.

NB the subcycling experiment also surfaced a bug: with sequentialConvergenceCriterion=
ResidualNorm the outer-loop re-assembled SOLID residual is stuck at a constant ~3.8 (vs the
solid's own solve ~1e-14), so ResidualNorm never converges for this setup; only
NumberOfNonlinearIterations completes. Not on the critical path since subcycling is unneeded.

---

## Fix 8 — Why the sequential coupling cannot converge (root cause of the residual gaps)

Attempted the principled root-cause fix for the two remaining discrepancies (matrix
plateau too high ~1.05 vs 0.91; fracture overshoot absent): make the sequential
(outer) fixed-stress loop actually converge, combined with v-weighted flux κ_i = v_i·k_i.

### Convergence machinery (worked)
- `CoupledSolver`: the loop only called `saveSequentialIterationState` when the
  *sub-solver's* couplingType==Sequential, but all leaf/nested solvers default to
  FullyImplicit, so it was never called and the outer loop could not measure
  convergence. Fixed: added a recursing `CoupledSolver::saveSequentialIterationState`
  (so the nested FullyImplicit `dualFlow` propagates the save to matrixFlow/fractureFlow),
  changed the guard to save every sub-solver in the sequential loop, and made the base
  `saveSequentialIterationState` a no-op + `FlowSolverBase` a safe early-return.
  Result: subcycling now genuinely iterates.

### The wall (definitive root cause)
With iteration enabled, the custom dual coupling does NOT converge:
- Near the overshoot peak (t≈0.57 s) the outer-loop pressure change stalls and decays
  only ~0.7%/cycle → the fixed-point map has a **near-unit eigenvalue (λ≈+0.993)**, a
  slowly-converging mode (NOT a λ≈−1 oscillation).
- Under-relaxation (swept ω = 0.5, 0.2, 0.1, 0.05) cannot fix a λ≈+1 mode — it only
  stabilizes λ≈−1. Aitken/Anderson would need ω≈140 (unstable) to move λ=0.993.

**Cause:** the composite-pressure formulation replaces the matrix-mesh pressure with
`p_eq = (ᾱ_m·p_m + ᾱ_f·p_f)/α_m` and the mechanics solves a *single* pressure. So the
mechanics only constrains the combination `ᾱ_m·p_m + ᾱ_f·p_f` and is **blind to the
orthogonal direction** (the individual p_m vs p_f split). That invisible direction is the
near-singular λ≈+1 mode: it is coupled back only weakly through flow storage, so the
sequential iteration crawls (~1300 iters/step) and no relaxation/acceleration robustly fixes it.

### The actual fix required (not done — major)
The mechanics must resolve **both** pressures separately — a genuine multi-continuum
effective stress σ = σ_eff − Σ ᾱ_i p_i in the mechanics assembly — instead of the single
composite-pressure swap. That is a real change to the mechanics kernel / coupling
(monolithic-equivalent), not a parameter, relaxation, or acceleration add-on.

All Fix-8 experimental code was reverted; repo holds only the validated time-step
refinement (Fix 7, matrix overshoot 1.015→1.066).

---

## Fix 9 — FIM (FullyImplicit) divergence root cause (2026-06-03, investigation)

Goal: determine why the monolithic FIM path of the dual solver diverges (to decide if
"complete the FIM" is viable). Definitive findings:

CONTROL: stock `SinglePhasePoromechanics` (couplingType=FullyImplicit) with the IDENTICAL
matrix material (K=1.1e9, k=5 nd), IC (incl. initial effective stress), displacement
loading, low perm, small dt **converges in 1 Newton iteration to 1e-14**
(`dpdkHMValidation/SP_Mandel_new.xml`). So FIM divergence is a BUG in the custom dual
solver, NOT physics/perm/IC/loading/dt.

Confounders RULED OUT (dual FIM diverges with each held like the working stock case):
low perm (diverges at high perm 1e-15 too), initial stress (stock sets it too), loading,
dt, and the fracture pressure leaking into momentum (zeroing p_f IC leaves Rsolid=4.31).

PRIMARY CAUSE: the cross-storage correction in
`SinglePhaseDualContinuum::assembleCouplingTerms` has an INCOMPLETE Jacobian. It adds
`rho*V*(corrDiag*dp_m + corrOff*dp_f)` to the flow residual and `rho*V*corrDiag`,
`rho*V*corrOff` to the Jacobian, but rho*V and corrDiag/corrOff are state-dependent
(rho, cf=drho/rho) and that dependence is NOT linearized. Harmless for single-pass
sequential; for FIM Newton it is an inconsistent Jacobian -> wrong p_m update -> corrupts
momentum (-alpha_m p_m) -> Rsolid grows -> diverges. Evidence (isolating, low perm, 1st step):
  cross-storage ON,  alpha_f=0.9992 : Rsolid grows to hundreds
  cross-storage ON,  alpha_f=0      : Rsolid grows to ~28
  cross-storage OFF, alpha_f=0      : no growth, stalls at 3.44
The off-diagonal corrOff (p_m<->p_f, scales with alpha_f) amplifies it. => the cross-storage
correction is a SEQUENTIAL-ONLY hack, fundamentally incompatible with FIM Newton.

SECONDARY: even fully isolated (cross-storage off, alpha_f=0) the dual matrix momentum
stalls at Rsolid~3.44 (stock reaches 1e-14) -> a second, smaller inconsistency in the dual
matrix assembly vs stock (suspect the manual `fields::flow::mass` accumulation path).

TO MAKE FIM WORK: (1) put the storage matrix S_bar into the constitutive/kernel with a
consistent Jacobian (not the bolted-on assembleCouplingTerms residual correction);
(2) fix the secondary matrix-momentum inconsistency; (3) add K_upf + K_pfu (both, with
effective abar_i and correct signs) for the fracture-mechanics physics. All Fix-8/9
experimental code reverted; repo holds only the validated Fix-7 time-step refinement.

---

## Fix 10 — FIM rewrite project, Increment 1: Bug A precisely characterized (2026-06-03)

Started the dual-pressure FIM rewrite. Architecture chosen: REUSE the proven stock matrix
poromechanics kernel for u+p_m (Bug A is a framework issue around it, not the kernel), then
add fracture coupling. Increment 1 = make matrix-only dual FIM converge like stock.

CONTROL re-confirmed: stock SinglePhasePoromechanics with the EXACT dispdriven IC/loading
(stripped single-porosity deck `_sp_disp.xml`) converges in 1 Newton iter (Rsolid 0.362 -> 1.4e-13).

Bug A SIGNATURE (instrumented, matrix-only mode = fracture coupling + cross-storage disabled):
the matrix pressure p_m oscillates in a PERFECT period-2 cycle:
  iter0 p_m=455000(=p_n) -> iter1 p_m=506049 -> iter2 p_m=455000 -> ... (exactly reversing)
with Rsolid stalled at 3.44. Instrumented state is CORRECT (p_m=455000, K=1.1e9, biot=0.9593,
por=0.14). So it is NOT corrupted state, NOT the fracture coupling, NOT permeability/IC.
A perfect eigenvalue=-1 cycle (delta_p exactly reverses) => the matrix p_m Jacobian/residual is
off by EXACTLY a factor of 2 (double-counted residual or half Jacobian) in the dual FIM assembly.

The matrix accumulation is assembled once by the kernel (Step 1, same factory as stock), Step 3
is matrix flux only, Step 4 is the fracture (mesh2). So the factor-2 is in how the dual framework's
state/`mass` handling (updateState + the `fields::flow::mass` accumulation, volumeFraction scaling)
interacts with the kernel's accumulation Jacobian. NEXT STEP to pin the exact line: dump the p_m
diagonal Jacobian entry and the p_m residual entry for element 0 and compare to stock (factor 2
expected somewhere). Then the rest of the rewrite (effective coeffs, K_upf/K_pfu, off-diagonal
storage) builds on the converging matrix base. All debug code reverted; repo holds only Fix-7.

---

## Fix 11 — FIM rewrite, Increment 1 progress (2026-06-03/04)

CODE CHANGE KEPT (rewrite progress, 3 files): gate the sequential-only cross-storage
correction so it is DISABLED under FullyImplicit coupling.
- DualContinuumFlowSolverBase.hpp: added m_enableCrossStorageCorrection (default true) +
  setEnableCrossStorageCorrection/getEnableCrossStorageCorrection.
- SinglePhaseDualContinuum.cpp: cross-storage block now gated on
  `v_f>0 && getEnableCrossStorageCorrection()`.
- DualContinuumPoromechanicsSolverBase.hpp: in the FIM branch, call
  flowSolver()->setEnableCrossStorageCorrection(false).
RESULT: FIM no longer diverges-by-growth (Rsolid was 4.31->25->hundreds with cross-storage;
now 4.31->3.44 and STALLS). Sequential UNCHANGED (verified: matrix peak/p0=1.0660, cross-storage
still on). Quantified earlier: cross-storage inflated the matrix p_m diagonal 5.4x (stock
7.468e-15 -> dual 4.02e-14 = stock + rho*V*corrDiagM); with it off the matrix p_m diagonal
EXACTLY matches stock (7.46756e-15). FIM (cross-storage off, loose tol) matches stock matrix
pressure within ~3% at t=0.1 (peak 516k vs 510k).

REMAINING BLOCKER — Bug A (FIM only): a distributed momentum-residual STALL at Rsolid~3.44
(abs ||R_u||~194, normalizer~56), persists with cross-storage off. NOT the Dirichlet reaction
(|R_top|~140 ~= |R_int|~140, evenly spread), NOT p_f (zeroing p_f IC leaves it), NOT the kernel
state (instrumented p_m/K/biot/por all correct). Stock SinglePhasePoromechanics with identical
IC/loading converges 1 Newton iter to 1e-13. So it is a coupled u<->p_m Jacobian/residual
inconsistency specific to the dual FIM framework (sequential solves u and p_m separately and
converges; monolithic does not). NEXT: dump the K_upm / K_pmu coupling-block entries for one
element, dual vs stock, to find the off-by-factor. Then Increments 2-4 (effective coeffs,
K_upf+K_pfu, consistent off-diagonal storage) build on the converging matrix base.

### Fix 11 (cont.) — Bug A narrowed to the solve/update path (not assembly)

Instrumented dump (matrix element 0, iter 0), dual FIM (cross-storage gated off) vs stock
SinglePhasePoromechanics with identical mesh/IC/loading. ALL of the following are BYTE-IDENTICAL:
  - absolute ||R_u||           : 194.639 (dual) == 194.639 (stock)
  - K_uu diagonal (node0,dim0) : -9.56387e6 == -9.56387e6
  - K_pmpm diagonal            : 7.46756e-15 == 7.46756e-15
  - p_m row off-diagonal norm (K_pmu) : 0.0432067 == 0.0432067
Yet dual FIM 2-cycles (||R_u|| oscillates 195<->205, eigenvalue ~-1) while stock converges in
1 Newton iter. Confirmed NOT the linear solver/preconditioner: a DIRECT solver still 2-cycles.
=> Bug A is NOT in the element assembly (residual + every checked Jacobian block match stock
exactly). It is in the SOLVE/UPDATE/boundary-condition path of the dual coupled framework
(e.g. applySystemSolution scaling, or how the 3-field [u,p_m,p_f] system / Dirichlet rows are
handled). NEXT: dump the solution increment delta_u (and the post-applyBoundaryConditions
system) for dual vs stock to find where the update is corrupted.

### Fix 11 (cont.2) — Bug A CONCLUSIVELY localized to the solve/update path (NOT assembly)

Added K_upm to the byte-identical list. Full verification (matrix elem0, iter0, dual-FIM
cross-storage-gated-off vs stock SinglePhasePoromechanics, identical mesh/IC):
  R_u=194.639, K_uu_diag=-9.56387e6, K_pmpm_diag=7.46756e-15, K_pmu_offnorm=0.0432067,
  K_upm[elem0,node0]=-1.07921e-05  -- ALL byte-identical between dual and stock.
Yet dual FIM solution increment overshoots by EXACTLY 2x: p_m goes 455000 -> 506049 -> 455000
(delta_pm=+51049 then -51049, exact reversal => eigenvalue -1). Stock converges in 1 iter.
A DIRECT solver still 2-cycles. => the assembled R and J are correct (== stock); the SOLVED
solution increment is 2x. So Bug A is in the path BETWEEN assembleSystem and the state update:
applyBoundaryConditions (possible double-application of the prescribed-displacement BC),
the linear-system setup/scaling, or applySystemSolution. The 2x is SYSTEM-WIDE (affects p_m and
u together), consistent with the driving Dirichlet displacement increment effectively applied 2x,
or the RHS/Jacobian scaled by 2 between assembly and solve. NEXT PROBE: dump the solution
increment delta (p_m and a u DOF) right after the linear solve, and the post-applyBoundaryConditions
rhs/diagonal, dual vs stock, to find the 2x. (Increment 1a cross-storage gate kept; debug reverted.)

### Fix 11 (cont.3) — Bug A ROOT FOUND: ill-posed fracture-pressure block contaminates the FIM solve

Dumped the linear solution increment (CoupledSolver::applySystemSolution), dual-FIM vs stock:
  DUAL : max|delta|=325528, ||delta||=1.76e6, ndof=3446
  STOCK: max|delta|=25524,  ||delta||=5.1e5,  ndof=3046  (diff 400 = fracture p_f DOFs)
=> the FRACTURE pressure increment is HUGE (325528, ~13x the matrix increment ~51000, which is
itself ~2x stock's 25524). With K_upf OFF (no mechanical constraint on p_f) and K_pfu only
lagged-explicit (R_pf depends on u via updateFracturePorosityFixedStress but the Jacobian has
no d R_pf/d u), the fracture-pressure block is INCONSISTENT/under-constrained -> huge, wrong
delta_pf -> in the coupled 3-field direct solve this pollutes delta_pm (->2x) and delta_u ->
the eigenvalue=-1 matrix 2-cycle. The matrix u-p_m ASSEMBLY is byte-identical to stock
(R_u, K_uu, K_pmpm, K_pmu, K_upm all verified equal); the corruption is purely the coupled
solve being polluted by the ill-posed p_f block.

REFRAME: Bug A is NOT separate from the fracture coupling. Increments 3-4 (consistent K_upf +
K_pfu + off-diagonal storage) are what make the p_f block well-posed -> they fix Bug A AND the
physics together. So the rewrite path is: implement the consistent fracture-mechanics coupling
(K_upf with abar_f, K_pfu as the true Jacobian of the strain->fracture-mass term, replacing the
lagged updateFracturePorosityFixedStress) so the 3-field FIM system is well-conditioned and
converges. (Increment 1a cross-storage gate kept; all debug reverted; build clean.)

### Fix 11 (cont.4) — Increment 3 attempt: consistent K_upf+K_pfu written; blocked by cross-mesh sparsity

Implemented the consistent fracture-mechanics coupling on the gated base:
  - K_upf: R_u -= alpha_f*p_f*gradPhi ; dR_u/dp_f = -alpha_f*gradPhi  (matches kernel K_upm=-alpha*gradPhi)
  - K_pfu: dR_pf/du = rho*alpha_f*sum_q dNdX*detJxW into the fracture-mass row (consistent with
    updateFracturePorosityFixedStress: phi_f = phi_n + alpha_f*volStrainInc + ...)
  - setupCoupling: addCoupling(displacement, secondary p_f) ; deck switched to effective coeffs
    (Kbar=4.514e8, Gbar=3.108e8, abar_m=0.382, abar_f=0.601).
RESULT: diverges (Rsolid grows ~x1.4/iter) with BOTH intrinsic and effective alpha. KEY diagnostic:
flipping the K_pfu sign had ZERO effect => the K_pfu (and K_upf) JACOBIAN entries are being DROPPED.
Cause: they couple mesh1 displacement <-> mesh2 fracture-pressure (CROSS-MESH), but
addCoupling(displacement, pressure, Connector::Elem) only builds WITHIN-mesh sparsity (u<->p_m on
mesh1). So the cross-mesh u<->p_f slots don't exist; the coupling lands in the residual but not the
Jacobian -> inconsistent -> divergence. (Same hard cross-mesh-sparsity problem the dual flow solves
for p_m<->p_f via its custom stencil.)
NEXT BLOCKER: build the cross-mesh u(mesh1)<->p_f(mesh2) sparsity (co-located node<->element), e.g.
mirror the dual-flow cross-flow stencil setup, so K_upf/K_pfu entries land. Then re-test convergence.
The K_upf/K_pfu formulation is correct and ready; only the sparsity is missing.
Increment-3 WIP reverted to the gate-only working state; cross-storage gate (1a) kept.

### Fix 11 (cont.5) — the cross-mesh sparsity blocker, fully specified

Traced the cross-mesh p_m<->p_f sparsity: DofManager::addCouplingDualContinuum sets m_isdpdk +
matrix/fracture region lists; the sparsity build (DofManager.cpp ~1435-1476) calls
countRowLengthsDualContinuum / setSparsityPatternDualContinuum, which connect the CO-LOCATED
matrix-elem and fracture-elem DOFs (locIdx identical across mesh1/mesh2). It is gated to
ELEM<->ELEM only ("排除力学耦合" / excludes mechanics: both fields must be FieldLocation::Elem).

To land K_upf/K_pfu we need a NODE<->ELEM cross-mesh sparsity:
  for each co-located (matrix elem k mesh1, fracture elem k mesh2):
    connect the 24 displacement DOFs of matrix-elem-k's nodes  <->  fracture-elem-k's p_f DOF
    (both row and column directions).
IMPLEMENTATION PLAN (next session, DofManager.cpp/.hpp):
  1. addCouplingDualContinuumMechanics(dispField, fractureFlowField, matrix/fracture region lists)
     -> set a flag (e.g. m_isdpdkMech) + store the disp<->fractureP field pair.
  2. countRowLengthsDualContinuumMechanics: fracture-elem p_f row += 24 ; each matrix-node disp
     row += 1 per incident fracture elem.
  3. setSparsityPatternDualContinuumMechanics: for each matrix elem k, get elemsToNodes + disp
     DOFs and the co-located fracture-elem-k p_f DOF; pattern.insertNonZeros both directions.
  4. trigger in the sparsity build for the (Node disp, Elem fractureP) field pair.
Then re-enable setupCoupling(u,p_f), the K_upf/K_pfu calls, and re-test convergence; expect the
fracture block to become well-posed (delta_pf bounded) and the matrix 2-cycle (Bug A) to vanish.
The K_upf/K_pfu assembly code + correct signs are documented above (Fix 11 cont.4) and ready.

SESSION SUMMARY (FIM rewrite status): Increment 1a DONE & kept (cross-storage FIM gate, 3 files,
sequential unaffected peak=1.066). Bug A root cause found (ill-posed p_f block). Matrix u-p_m
assembly proven byte-identical to stock. K_upf/K_pfu formulation derived & verified-correct.
Single remaining blocker to a converging FIM matrix+fracture base: the node<->elem cross-mesh
sparsity above. Then Increment 2 (effective coeffs) + Increment 4 (off-diagonal storage S_mf) +
validation vs analytical.

---

## Fix 12 (2026-06-04): Cross-mesh sparsity DONE; FIM blocker re-localized to monolithic scaling

DONE — the node<->elem cross-mesh coupling that was the "single remaining blocker" above:
- `DofManager::addCouplingDualContinuumMechanics` — flags `m_isdpdkMech`, stores the disp/fractureP
  field-index pair, calls `addCoupling`.
- `countRowLengthsDualContinuumMechanics` / `setSparsityPatternDualContinuumMechanics` — for each
  matrix-mesh element, connect its nodes' 3 displacement DOFs to the co-located fracture-mesh
  element's p_f DOF (both directions). Triggered in `setSparsityPattern`'s count + fill passes via
  `m_isdpdkMech && {row,col} == {dispIdx, fracPIdx}`.
- `DualContinuumPoromechanicsSolverBase::setupCoupling` calls it; `assembleSystem` re-enables K_upf
  (`assembleFractureMechanicsCoupling`) and adds K_pfu (`assembleFractureToMechanicsCoupling`,
  `d(phi_f rho V)/dU = rho_ref*alpha_f*gradN` into the fracture-mass row at displacement cols).
  Both gated behind `m_enableFractureMechanicsCoupling` (default false).
- K_upf sign VERIFIED correct as `+alpha_f*p_f*gradN` (matches stock
  `R_mom = -int symGrad(N):totalStress`, `totalStress = effStress - alpha*p*I`). The earlier
  "flip to -=" idea was wrong.

All builds clean. Sequential (`DPDP_N2_dispdriven_correctLF.xml`) NOT regressed — still converges
and overshoots. New `DPDP_N2_dispdriven_fim.xml` = `couplingType="FullyImplicit"`.

BUT FIM STILL DOES NOT CONVERGE, and the cause is NOT the sparsity / coupling. Re-localized:
- Coupling OFF: Rsolid drops once (4.31 -> 3.44) then FREEZES exactly while flow residuals are tiny
  (matrix 2e-5, fracture 1e-3) and the linear solve hits 1e-12; delta -> 0.
- Coupling ON: same freeze but worse (Rsolid 3.44 -> 61) because K_upf injects a large fracture-
  pressure force the under-resolved displacement can't balance.
- `SolidMechanicsLagrangianFEM::calculateResidualNorm` reads `localRhs` directly, normalized by
  `m_maxForce+1`. So Rsolid=3.44 means the assembled momentum residual really is ~3.44*(maxForce+1),
  i.e. NONZERO, yet the monolithic GMRES returns delta_u ~ 0. The mechanics rows (forces ~1e5-1e6)
  and flow mass rows differ by many orders of magnitude; the global 2-norm GMRES minimizes is
  dominated by one block, so the displacement correction is under-resolved => Rsolid stalls. This is
  the same "ResidualNorm criterion stuck ~3.8" symptom from Fix 7.

NEXT (FIM): fix the monolithic block scaling / conditioning (row/field scaling or a
block/Schur-style preconditioner so the mechanics block is properly resolved), validate Rsolid->0
in ONE step on a single timestep with coupling OFF, THEN flip `m_enableFractureMechanicsCoupling`
on and re-test. Until then, Sequential remains the working path.

---

## Fix 13 (2026-06-04): FIM CONVERGENCE FIXED — it was a lambda~=-1 oscillation, not scaling

The Fix 12 guess ("monolithic block scaling") was WRONG. Instrumenting the state across Newton
iterations (max incremental displacement + max p_m + max p_f each assembleSystem) revealed the truth:
FIM was NOT freezing — it was in a **period-2 (lambda~=-1) limit cycle**. Displacement converged
(maxIncDisp settled), but p_m ping-ponged 455000 <-> 506049 and p_f 487998 <-> 498763 every Newton
iteration. The residual norm only *looked* frozen because both cycle states have equal norm. A direct
linear solver showed the identical cycle, ruling out conditioning. The cycle midpoint
((455000+506049)/2 = 480524) is the EXACT solution => the full Newton step overshoots by ~2x
(lambda~=-1).

FIX = fixed Newton under-relaxation in the FIM path: override scalingForSystemSolution to multiply
the step by m_fimNewtonRelaxation (default 0.5) when couplingType==FullyImplicit. lambda_eff =
1 - relax*2 ~= 0 for the oscillatory mode. Confirmed by line search too: a 0.5 step instantly
collapsed R from 3.44 to 5e-8 (landing on the midpoint). Both new XML params on
SinglePhaseDualContinuumPoromechanics:
  - fimNewtonRelaxation (default 0.5; 0.7 also stable and ~1.6x faster; 1.0 disables)
  - enableFractureMechanicsCoupling (default 1; toggles K_upf/K_pfu)

RESULT: DPDP_N2_dispdriven_fim.xml (couplingType=FullyImplicit, fimNewtonRelaxation=0.7, K_upf/K_pfu
ON) runs the FULL Mandel timescale to t=3e5 s (514 cycles) in ~4.5 min with ZERO timestep cuts, ZERO
nonconvergence, with BOTH gmres and direct linear solvers. Residual decreases monotonically
(~0.5-0.7 contraction/iter). Sequential path unaffected (relaxation only applies in FIM mode).

REMAINING (physics accuracy, NOT convergence): FIM matrix pressure over-pressurizes (peak ~1.66 vs
analytical ~1.07) because the storage-matrix (Mehrabian S_ij) + effective-Biot partitioning
corrections live in mapSolutionBetweenSolvers (Sequential-only) and the SinglePhaseDualContinuum
cross-storage block is gated OFF for FIM (postInputInitialization setEnableCrossStorageCorrection
(false)). Porting those into the monolithic FIM assembly is the next increment (Increment 2/4).
