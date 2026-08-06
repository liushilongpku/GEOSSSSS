/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumPoromechanicsSolverBase.hpp
 *
 * @brief A coupled solver that binds dual continuum flow solvers with solid mechanics for poromechanics in dual continuum models.
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMPOROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumFlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsKernelsDispatchTypeList.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/FieldIdentifiers.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/logger/Logger.hpp"
#include <cmath>
#include <vector>

namespace geos
{

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER, typename MECHANICS_SOLVER = SolidMechanicsLagrangianFEM >
class DualContinuumPoromechanicsSolverBase : public PoromechanicsSolver< DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >, MECHANICS_SOLVER >
{
public:

  using Base = PoromechanicsSolver< DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >, MECHANICS_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /// Compile-time flag: true when the flow continua are compositional multiphase.
  /// Drives if-constexpr dispatch in the FIM assembly path so the single-phase code
  /// is preserved byte-for-byte when PRIMARY_FLOW_SOLVER = SinglePhaseBase.
  static constexpr bool isMultiphaseFlow =
    std::is_base_of_v< CompositionalMultiphaseBase, PRIMARY_FLOW_SOLVER >;

  enum class SolverType : integer
  {
    DualContinuumFlow = 0,
    SolidMechanics = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "dualcontinuumporomechanics"; }

  /**
   * @brief main constructor for DualContinuumPoromechanicsSolverBase Objects
   * @param name the name of this instantiation of DualContinuumPoromechanicsSolverBase in the repository
   * @param parent the parent group of this instantiation of DualContinuumPoromechanicsSolverBase
   */
  DualContinuumPoromechanicsSolverBase( const string & name,
                                        dataRepository::Group * const parent )
    : Base( name, parent )
  {
    this->registerWrapper( viewKeyStruct::fractureVolumeFractionString(), &m_fractureVolumeFraction ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( -1.0 ).
      setDescription( "Fracture volume fraction v_f = V_f / (V_m+V_f). "
                      "Default (<0) computes v_f from mesh element volumes. "
                      "Set explicitly when the mesh does not reflect the "
                      "physical volume fractions (e.g. dual-continuum models "
                      "with co-located meshes)." );

    this->registerWrapper( viewKeyStruct::fimNewtonRelaxationString(), &m_fimNewtonRelaxation ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 0.5 ).
      setDescription( "FullyImplicit Newton step under-relaxation factor (0,1]. The monolithic "
                      "dual-porosity pressure blocks produce a period-2 (lambda~=-1) Newton "
                      "oscillation; 0.5 damps it (midpoint = exact solution). Set 1.0 to disable." );

    this->registerWrapper( viewKeyStruct::sequentialPressureRelaxationString(), &m_sequentialPressureRelaxation ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 1.0 ).
      setDescription( "Sequential outer-loop pressure under-relaxation factor (0,1]. Values below 1 "
                      "mix the dual-flow pressure solution with the previous outer iterate before "
                      "mapping pressure to mechanics. Set 1.0 to disable." );

    this->registerWrapper( viewKeyStruct::enableFractureMechanicsCouplingString(), &m_enableFractureMechanicsCoupling ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 1 ).
      setDescription( "Assemble the explicit fracture<->mechanics Jacobian coupling (K_upf / K_pfu) "
                      "in the FullyImplicit path. Requires the cross-mesh node<->elem sparsity." );

    this->registerWrapper( viewKeyStruct::enableFracturePorosityStrainCouplingString(), &m_enableFracturePorosityStrainCoupling ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 1 ).
      setDescription( "Diagnostic FullyImplicit switch for the mechanics->fracture-mass half of the "
                      "fracture-mechanics coupling. When 0, K_upf is kept but the fracture porosity "
                      "strain injection and K_pfu are disabled." );

    this->registerWrapper( viewKeyStruct::enableFimCrossStorageString(), &m_enableFimCrossStorage ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 1 ).
      setDescription( "Enable the multi-porosity (Mehrabian S_ij) cross-storage correction in the "
                      "FullyImplicit path (adds the matrix<->fracture off-diagonal storage that "
                      "drives the Mandel-Cryer matrix decay)." );

    this->registerWrapper( viewKeyStruct::logFimCouplingDiagnosticsString(), &m_logFimCouplingDiagnostics ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 0 ).
      setDescription( "Diagnostic switch for FullyImplicit dual-continuum poromechanics. When nonzero, "
                      "logs norms of hand-assembled FIM coupling blocks and pressure updates "
                      "(K_upf, cross-storage, dp_m/dp_f). Does not change the assembled system." );

    this->registerWrapper( viewKeyStruct::logSequentialMassDiagnosticsString(), &m_logSequentialMassDiagnostics ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 0 ).
      setDescription( "Diagnostic switch for sequential dual-continuum poromechanics. When nonzero, "
                      "logs matrix and fracture CO2 inventories before and after the mechanics-to-flow "
                      "porosity/state updates. Does not change the solution." );

    this->registerWrapper( viewKeyStruct::useIntrinsicInputString(), &m_useIntrinsicInput ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 1 ).
      setDescription( "Mechanics/storage input mode for the dual-continuum effective medium. "
                      "1 (default): the deck sets the INTRINSIC matrix/fracture moduli/Biot on the constitutive "
                      "models; the solver uses its default homogenization algorithm "
                      "(Reuss Kbar/Gbar, abar_i=Kbar*v_i*alpha_i/K_i). FullyImplicit and Sequential "
                      "both homogenize once at initialization. "
                      "Storage is reconstructed from intrinsic parameters. "
                      "0: the deck sets all homogenized EFFECTIVE moduli/Biot and the direct effective "
                      "storage coefficients on DualContinuumCrossFlow; the solver does not homogenize "
                      "these values again. Porosity and permeability are always intrinsic inputs and "
                      "are volume-fraction scaled by the dual-continuum flow solver in both modes." );

    this->registerWrapper( viewKeyStruct::autoInitializeStressString(), &m_autoInitializeStress ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setApplyDefaultValue( 0 ).
      setDescription( "Automatically initialize the matrix effective stress to the dual-continuum total "
                      "Biot stress sigma' = (alpha_m*p_m + alpha_f*p_f)*I so the initial total stress is "
                      "in equilibrium (Rsolid~0 at t=0), the same way single-porosity poromechanics starts "
                      "balanced. Removes the need to set matrixSolid_stress by hand. Assumes isotropic "
                      "initial stress and zero external load; for anisotropic/loaded initial states set "
                      "this 0 and prescribe matrixSolid_stress manually. Default 0 (backward compatible); "
                      "the compositional dual solver defaults it to 1." );
  }

  // Override postInputInitialization to include checks for dual continuum
  virtual void postInputInitialization() override
  {
    Base::postInputInitialization();

    setMGRStrategy();

    GEOS_THROW_IF( m_useIntrinsicInput != 0 && m_useIntrinsicInput != 1,
                   GEOS_FMT( "{}: useIntrinsicInput must be 0 (effective input, no homogenization) "
                             "or 1 (intrinsic input, apply default homogenization), got {}.",
                             this->getName(), m_useIntrinsicInput ),
                   InputError );
    GEOS_THROW_IF( m_sequentialPressureRelaxation <= 0.0 || m_sequentialPressureRelaxation > 1.0,
                   GEOS_FMT( "{}: {} must be in (0, 1], got {}.",
                             this->getName(), viewKeyStruct::sequentialPressureRelaxationString(),
                             m_sequentialPressureRelaxation ),
                   InputError );

    // Additional checks for dual continuum, e.g., ensure flow solvers are compatible
    GEOS_THROW_IF( this->flowSolver()->primarySolver()->isThermal() != this->flowSolver()->secondarySolver()->isThermal(),
                   GEOS_FMT( "{} {}: Primary and secondary flow solvers must have the same thermal setting",
                             this->getCatalogName(), this->getName() ),
                   InputError );

    validateDualContinuumVolumeFractions();
    validateMaterialInputMode( this->template getGroupByPath< DomainPartition >( "/Problem/domain" ) );

    // Enable fixed-stress poromechanics update BEFORE registerDataOnMesh is called.
    // registerDataOnMesh runs on sub-solvers before the coupled solver, so the flags
    // must be set during postInputInitialization for pressure_k to be registered.
    if( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
    {
      GEOS_THROW_IF( this->getNonlinearSolverParameters().m_subcyclingOption == 0,
                     GEOS_FMT( "{} {}: dual-continuum poromechanics Sequential coupling requires "
                               "subcycling=\"1\". A single-pass split leaves the fixed-stress "
                               "stabilization term in the accepted flow solution and changes the "
                               "effective storage/drainage time scale.",
                               this->getCatalogName(), this->getName() ),
                     InputError );
      GEOS_THROW_IF( this->getNonlinearSolverParameters().sequentialConvergenceCriterion() !=
                     NonlinearSolverParameters::SequentialConvergenceCriterion::SolutionIncrements,
                     GEOS_FMT( "{} {}: dual-continuum poromechanics Sequential coupling requires "
                               "sequentialConvergenceCriterion=\"SolutionIncrements\". ResidualNorm is "
                               "ambiguous after the dual-continuum pressure/state mapping and can report "
                               "false convergence while the matrix/fracture pressures are still changing.",
                               this->getCatalogName(), this->getName() ),
                     InputError );
      this->solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
      this->flowSolver()->enableFixedStressPoromechanicsUpdate();
    }
    else
    {
      // FullyImplicit: the multi-porosity (Mehrabian S_ij) cross-storage correction adds the
      // v-weighted diagonal + the matrix<->fracture off-diagonal storage that lets the matrix
      // shed pressure as the fracture drains (the Mandel-Cryer decay). Its Jacobian is assembled
      // consistently and the lambda~=-1 oscillation it used to trigger is now damped by the FIM
      // Newton under-relaxation (scalingForSystemSolution). Keep it on for FIM, gated by
      // m_enableFimCrossStorage so it can be toggled for debugging.
      this->flowSolver()->setEnableCrossStorageCorrection( m_enableFimCrossStorage != 0 );
    }
  }

  // Override setupDofs to handle dual continuum DOFs
  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    // Setup DOFs for dual continuum flow and solid mechanics
    this->flowSolver()->setupDofs( domain, dofManager );//双重介质的耦合关系已经在DualContinuumFlowSolverBase的setupDofs中设置了，这里直接调用即可
    this->solidMechanicsSolver()->setupDofs( domain, dofManager );
    this->setupCoupling( domain, dofManager );//这里仅需要设置流体和固体的耦合关系
  }

  // Setup force coupling between mechanics and flow fields
  virtual void setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                              DofManager & dofManager ) const override
  {
    // Setup coupling between solid mechanics displacement and both dual-continuum flow fields.
    // 只需要设置一个基质的流固耦合关系
    dofManager.addCoupling( fields::solidMechanics::totalDisplacement::key(),
                            PRIMARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString(),
                            DofManager::Connector::Elem );

    // Fracture (secondary mesh) mechanics<->flow coupling: the displacement DOFs (Node, mesh1)
    // and the fracture-pressure DOFs (Elem, mesh2) live on different, co-located meshes, so the
    // ordinary within-mesh addCoupling cannot create the u<->p_f sparsity. Use the dedicated
    // cross-mesh node<->elem coupling so the K_upf / K_pfu Jacobian entries are not dropped.
    string_array const & matrixRegionList =
      this->flowSolver()->template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList =
      this->flowSolver()->template getReference< string_array >( "fractureRegionList" );

    dofManager.addCouplingDualContinuumMechanics( matrixRegionList,
                                                  fractureRegionList,
                                                  fields::solidMechanics::totalDisplacement::key(),
                                                  SECONDARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString(),
                                                  DofManager::Connector::Elem );
  }

  virtual void setMGRStrategy()
  {
    LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();

    if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
      return;

    GEOS_THROW_IF( this->m_isThermal,
                   GEOS_FMT( "{}: MGR strategy is not implemented for thermal {}",
                             this->getName(), this->getCatalogName() ),
                   InputError );

    linearSolverParameters.mgr.separateComponents = true;
    linearSolverParameters.dofsPerNode = 3;

    if constexpr ( isMultiphaseFlow )
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphaseDualContinuumPoromechanics;
    }
    else
    {
      GEOS_THROW( GEOS_FMT( "{}: MGR strategy is not implemented for {}",
                            this->getName(), this->getCatalogName() ),
                  InputError );
    }

    GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolver,
                           GEOS_FMT( "{}: MGR strategy set to {}", this->getName(),
                                     EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString(
                                       linearSolverParameters.mgr.strategy ) ) );
  }

  // TODO@LSL: Monolithic poromechanics kernel for the MATRIX continuum (u + p_m).
  //
  // WHY MONOLITHIC (for matrix only):
  //   The block-by-block approach (flow-solver → mech-solver → coupling-kernel) runs
  //   three separate element traversals.  Each traversal evaluates the constitutive
  //   model at a slightly different state point during the Newton iteration, producing
  //   Jacobian blocks that are individually correct but mutually inconsistent.
  //
  //   When the matrix permeability is very low (k_matrix ≪ 1e-13 m²), the pressure
  //   diagonal block K_pp shrinks to ~10⁻¹⁴ while the off-diagonal coupling
  //   K_pu·K_uu⁻¹·K_up stays at ~10⁻¹⁰.  The Schur complement S = K_pp − K_pu·K_uu⁻¹·K_up
  //   becomes dominated by the inconsistent off-diagonal terms, can flip sign,
  //   and causes the Newton solver to diverge.
  //
  //   This monolithic kernel fuses the matrix momentum balance + matrix mass balance
  //   in a SINGLE quadrature-point loop.  One constitutive evaluation per integration
  //   point produces all blocks (K_uu, K_upm, K_pmu, K_pmpm) from the SAME state,
  //   guaranteeing a self-consistent Schur complement irrespective of permeability.
  //
  // WHAT REMAINS BLOCK-BY-BLOCK:
  //   - Fracture flow (K_pfpf): still assembled by SinglePhaseFVM.  Fracture permeability
  //     is high (> 1e-15 m²) so its K_pp is well-conditioned.
  //   - Cross-flow coupling (K_pmpf, K_pfpm): assembled by DualContinuumFVM.
  //   - Fracture-mechanics coupling (K_upf, K_pfu): omitted (α_f ≈ 0 for most shales).
  //
  // The monolithic kernel reuses the existing SinglePhasePoromechanics kernel,
  // which is already well-tested for single-porosity poromechanics, via the
  // PoromechanicsSolver::assemblyLaunch infrastructure.

  virtual void assembleSystem( real64 const time_n,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override
  {
    // ---- Step 0: Map fracture data (p_f, α_f, DOF#) from mesh2 to mesh1 ----
    mapFractureDataToMatrix( domain, dofManager );

    // ---- Step 1: Monolithic matrix kernel (u + p_m [+ components]) ----
    // Assembles K_uu, K_upm, K_pmu, K_pmpm in one quadrature-point loop
    // on mesh1 (matrix) regions only.  For multiphase flow the same monolithic
    // matrix block is assembled with the compositional poromechanics kernel
    // (couples displacement with pressure + all component DOFs).
    this->template forDiscretizationOnMeshTargets<>(
      domain.getMeshBodies(),
      [&]( string const & meshBodyName,
           MeshLevel & mesh,
           string_array const & regionNames )
    {
      if( meshBodyName != "mesh1" ) return;   // matrix mesh only

      if constexpr ( isMultiphaseFlow )
      {
        string const flowDofKey =
          dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

        this->template assemblyLaunch<
          PoromechanicsKernelsDispatchTypeList,
          poromechanicsKernels::MultiphasePoromechanicsKernelFactory >(
            mesh, dofManager, regionNames,
            viewKeyStruct::porousMaterialNamesString(),
            localMatrix, localRhs, dt,
            flowDofKey,
            this->flowSolver()->primarySolver()->numFluidComponents(),
            this->flowSolver()->primarySolver()->numFluidPhases(),
            this->flowSolver()->primarySolver()->useSimpleAccumulation(),
            this->flowSolver()->primarySolver()->useTotalMassEquation(),
            this->m_performStressInitialization,
            FlowSolverBase::viewKeyStruct::fluidNamesString() );
      }
      else
      {
        string const flowDofKey =
          dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

        this->template assemblyLaunch<
          PoromechanicsKernelsDispatchTypeList,
          poromechanicsKernels::SinglePhasePoromechanicsKernelFactory >(
            mesh, dofManager, regionNames,
            viewKeyStruct::porousMaterialNamesString(),
            localMatrix, localRhs, dt,
            flowDofKey,
            this->m_performStressInitialization,
            FlowSolverBase::viewKeyStruct::fluidNamesString() );
      }
    } );

    // ---- Step 1.5: Update fracture porosity with matrix strain (fixed-stress) ----
    updateFracturePorosityFixedStress( domain, time_n + dt );

    // ---- Step 2: Fracture-mechanics coupling (K_upf) ----
    // Fracture pressure -> displacement residual (+alpha_f * p_f * gradN), consistent with the
    // GEOS total-stress convention. The cross-mesh u<->p_f sparsity is now provided by
    // DofManager::addCouplingDualContinuumMechanics, so these entries are no longer dropped.
    if( m_enableFractureMechanicsCoupling )
      assembleFractureMechanicsCoupling( domain, dofManager, localMatrix, localRhs );

    // ---- Step 2b: Mechanics->fracture-mass coupling (K_pfu) ----
    // Displacement -> fracture mass Jacobian. Makes the monolithic Newton consistent with the
    // strain term injected into the fracture porosity (Step 1.5), instead of relying on a
    // Picard-style porosity refresh alone. Single-phase: d(phi_f rho V)/dU = rho alpha_f gradN.
    // Multiphase: per-component d(phi_f * compDens_c * V)/dU = compDens_c * alpha_f * gradN
    // (with the useTotalMassEquation row transform applied to match the fracture accumulation).
    if( m_enableFractureMechanicsCoupling && m_enableFracturePorosityStrainCoupling )
    {
      if constexpr ( isMultiphaseFlow )
        assembleFractureToMechanicsCouplingMultiphase( domain, dofManager, localMatrix, localRhs );
      else
        assembleFractureToMechanicsCoupling( domain, dofManager, localMatrix, localRhs );
    }

    // ---- Step 3: Matrix face-based flux terms ----
    this->flowSolver()->primarySolver()->assembleFluxTerms(
      dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4: Fracture flow assembly (K_pfpf accumulation + face fluxes) ----
    this->flowSolver()->secondarySolver()->assembleSystem(
      time_n, dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4: Cross-flow coupling (K_pmpf, K_pfpm) ----
    this->flowSolver()->assembleCouplingTerms(
      time_n, dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4b: Multi-porosity (Mehrabian) cross-storage (multiphase) ----
    // Inter-continuum storage coupling d(M_i,c)/d(p_j): the matrix fluid content responds to the
    // fracture pressure (and vice-versa) through the shared solid skeleton. The single-phase path
    // assembles this inside SinglePhaseDualContinuum::assembleCouplingTerms; the compositional flow
    // coupling only does the transfer flux, so for multiphase we add it here.
    if constexpr ( isMultiphaseFlow )
    {
      if( m_enableFimCrossStorage )
        assembleFimCrossStorageMultiphase( domain, dofManager, localMatrix, localRhs );
    }
  }

  // FIM Newton under-relaxation. The monolithic dual-porosity pressure blocks produce a
  // period-2 (lambda~=-1) Newton oscillation: the full step overshoots the solution by ~2x and
  // the iteration ping-pongs between two states whose midpoint is the exact solution. A fixed
  // step scaling of m_fimNewtonRelaxation (~0.5) collapses that mode in one step (lambda_eff =
  // 1 - relax*2 ~= 0) and keeps all other modes contracting. Applied only in FullyImplicit mode.
  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override
  {
    real64 const baseScaling = Base::scalingForSystemSolution( domain, dofManager, localSolution );
    // Apply oscillation detection from PhysicsSolverBase (not called by sub-solver overrides)
    real64 const oscScaling = PhysicsSolverBase::scalingForSystemSolution( domain, dofManager, localSolution );
    real64 scalingFactor = baseScaling * oscScaling;
    if( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      scalingFactor *= m_fimNewtonRelaxation;
    }
    if( m_logFimCouplingDiagnostics != 0 &&
        this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      logFimPressureUpdateDiagnostics( domain, dofManager, localSolution, scalingFactor );
    }
    return scalingFactor;
  }

  virtual void updateState( DomainPartition & domain ) override
  {
    Base::updateState( domain );

    if constexpr ( isMultiphaseFlow )
    {
      if( this->getNonlinearSolverParameters().couplingType() ==
          NonlinearSolverParameters::CouplingType::FullyImplicit )
      {
        enforceFractureCompositionalVolumeClosure( domain, false );
      }
    }
  }

  // Override mapSolutionBetweenSolvers to handle dual-mesh sequential coupling.
  // 1. After flow: add -α_f·p_f·I to stored effective stress (so mechanics kernel
  //    computes correct total stress without needing K_upf).
  // 2. After mechanics: restore stress, average mean stress on mesh1 only,
  //    update porosity/permeability on mesh1 only.
  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override
  {
    bool const useCompositePressure = useCompositePressureForSequential();
    if( solverType == static_cast< integer >( SolverType::DualContinuumFlow ) )
    {
      // The dual-flow subsolver must remain FullyImplicit internally so matrix and
      // fracture pressures are solved together. Therefore the generic outer
      // CoupledSolver does not call saveSequentialIterationState() for it. Do it
      // here before mapping pressures to mechanics so pressure_k records the previous
      // accepted outer iterate and SolutionIncrements sees the applied flow change.
      relaxSequentialFlowPressure( domain );
      this->flowSolver()->saveSequentialIterationState( domain );
      copyFracturePressureToMesh1( domain );
      if( useCompositePressure )
      {
        swapToCompositePressure( domain );
      }
      else
      {
        m_tempCompositePressure.clear();
        m_tempCompositePressure_n.clear();
      }
      Base::updateBulkDensity( domain );
      if( m_logSequentialMassDiagnostics != 0 )
      {
        logSequentialMassInventory( domain, "afterFlow" );
      }
    }
    else if( solverType == static_cast< integer >( SolverType::SolidMechanics )
             && !this->m_performStressInitialization )
    {
      if( useCompositePressure )
      {
        restoreCompositePressure( domain );
      }

      // Average mean total stress on mesh1 only (mesh2 has no FE)
      MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
      MeshLevel & meshLevel = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
      using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
      DualFlow & dualFlow = *this->flowSolver();
      string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
      m_tempVolStrainIncr.clear();
      meshLevel.getElemManager().forElementSubRegions< CellElementSubRegion >( matrixRegions,
        [&]( localIndex const, auto & subRegion )
      {
        string const & solidName = subRegion.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & solid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );

        arrayView1d< real64 > const averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();

        finiteElement::FiniteElementBase & subRegionFE =
          subRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&]( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );
          computeAverageVolumetricStrainIncrement< FE_TYPE >(
            meshLevel, subRegion, finiteElement, averageMeanTotalStressIncrement_k );
        } );

        if( subRegion.hasWrapper( viewKeyStruct::fracturePressureString() ) )
        {
          arrayView1d< real64 const > const K_view = solid.getBulkModulus();
          arrayView1d< real64 const > const alpha_view = solid.getBiotCoefficient();
          arrayView1d< real64 const > const pressure =
            subRegion.template getField< fields::flow::pressure >();
          arrayView1d< real64 const > const pressure_n =
            subRegion.template getField< fields::flow::pressure_n >();
          for( localIndex k = 0; k < subRegion.size(); ++k )
          {
            real64 const delta_eps_v = averageMeanTotalStressIncrement_k[k];
            real64 const totalStressPressureIncrement = alpha_view[k] * ( pressure[k] - pressure_n[k] );
            if( useCompositePressure )
            {
              // Intrinsic-input Sequential path: mechanics used p_eq and temporary K_eff.
              // Feed each intrinsic continuum with v_i*K_eff*dEps - alpha_i*dp_i so the fixed-stress
              // porosity update receives the mean total stress increment.
              arrayView1d< real64 const > const K_eff_view =
                subRegion.template getReference< array1d< real64 > >( viewKeyStruct::effectiveBulkModulusString() );
              averageMeanTotalStressIncrement_k[k] =
                ( 1.0 - m_fractureVolumeFraction ) * K_eff_view[k] * delta_eps_v
                - totalStressPressureIncrement;
            }
            else
            {
              // Effective-input Sequential path: mechanics used physical p_m/p_f with
              // effective Biot coefficients directly.
              averageMeanTotalStressIncrement_k[k] = K_view[k] * delta_eps_v - totalStressPressureIncrement;
            }
            m_tempVolStrainIncr.push_back( delta_eps_v );
          }
        }
        if( m_logSequentialMassDiagnostics != 0 )
        {
          logSequentialMassInventory( domain, "beforeMatrixPorosity" );
        }
        this->flowSolver()->updatePorosityAndPermeability( subRegion );
        if( m_logSequentialMassDiagnostics != 0 )
        {
          logSequentialMassInventory( domain, "afterMatrixPorosity" );
        }
        this->flowSolver()->primarySolver()->updateFluidState( subRegion );
        if( m_logSequentialMassDiagnostics != 0 )
        {
          logSequentialMassInventory( domain, "afterMatrixFluidState" );
        }
        this->updateBulkDensity( subRegion );
      } );
      if( m_logSequentialMassDiagnostics != 0 )
      {
        logSequentialMassInventory( domain, "afterMatrixState" );
      }
      m_tempCompositePressure.clear();
      m_tempCompositePressure_n.clear();

      // Copy matrix avgStressIncr to fracture so the fracture flow
      // solver's fixed-stress porosity update includes the mechanics strain.
      {
        MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
        MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
        DualFlow & ff = *this->flowSolver();
        string_array const & matRegions = ff.template getReference< string_array >( "matrixRegionList" );
        string_array const & fracRegions = ff.template getReference< string_array >( "fractureRegionList" );
        localIndex strainIdx = 0;

        for( size_t iPair = 0; iPair < matRegions.size(); ++iPair )
        {
          ElementRegionBase & matRegion = meshLevel.getElemManager().getRegion( matRegions[ iPair ] );
          ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fracRegions[ iPair ] );

          matRegion.forElementSubRegionsIndex< CellElementSubRegion >(
            [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
          {
            if( subRegIdx >= fracRegion.numSubRegions() ) return;
            CellElementSubRegion & fracSubReg =
              dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

            string const & fracSolidName =
              fracSubReg.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
            constitutive::CoupledSolidBase & fracSolid =
              this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
            arrayView1d< real64 > const fracAvgStressIncr = fracSolid.getAverageMeanTotalStressIncrement_k();
            arrayView1d< real64 const > const K_f = fracSolid.getBulkModulus();
            arrayView1d< real64 const > const alpha_f = fracSolid.getBiotCoefficient();
            arrayView1d< real64 const > const p_f = fracSubReg.template getField< fields::flow::pressure >();
            arrayView1d< real64 const > const p_f_n = fracSubReg.template getField< fields::flow::pressure_n >();

            // Fracture receives the SAME shared volumetric strain increment dEps_v as
            // the matrix. Both effective-input and intrinsic-input Sequential use the
            // homogenized fracture constitutive Biot as abar_f, so avgStress_f=K_f*dEps_v
            // gives abar_f*dEps_v. In both cases subtract
            // alpha_f*dp_f because BiotPorosity expects a mean total stress increment.
            if( matSubReg.hasWrapper( viewKeyStruct::effectiveBulkModulusString() ) )
            {
              arrayView1d< real64 const > const K_eff_frac =
                matSubReg.getReference< array1d< real64 > >( viewKeyStruct::effectiveBulkModulusString() );
              arrayView1d< localIndex const > const matrixToFracture =
                matSubReg.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );
              real64 const v_f = m_fractureVolumeFraction;
              for( localIndex k = 0; k < matSubReg.size(); ++k )
              {
                localIndex const kf = matrixToFracture[k];
                GEOS_ERROR_IF( kf < 0 || kf >= fracSubReg.size(),
                               "Invalid dual-continuum fracture stress copy connectivity for matrix subregion "
                               << matSubReg.getName() << ", local element " << k
                               << ": mapped fracture element " << kf
                               << " is outside fracture subregion " << fracSubReg.getName()
                               << " size " << fracSubReg.size() );
                GEOS_ERROR_IF( strainIdx >= static_cast< localIndex >( m_tempVolStrainIncr.size() ),
                               "Missing volumetric strain increment for dual-continuum fracture mapping." );
                real64 const delta_eps_v = m_tempVolStrainIncr[strainIdx++];
                real64 const pressureStressIncrement = alpha_f[kf] * ( p_f[kf] - p_f_n[kf] );
                fracAvgStressIncr[kf] =
                  useCompositePressure
                  ? v_f * K_eff_frac[k] * delta_eps_v - pressureStressIncrement
                  : K_f[kf] * delta_eps_v - pressureStressIncrement;
              }
            }
            if( m_logSequentialMassDiagnostics != 0 )
            {
              logSequentialMassInventory( domain, "beforeFracturePorosity" );
            }
            this->flowSolver()->secondarySolver()->updatePorosityAndPermeability( fracSubReg );
            if( m_logSequentialMassDiagnostics != 0 )
            {
              logSequentialMassInventory( domain, "afterFracturePorosity" );
            }
            this->flowSolver()->secondarySolver()->updateFluidState( fracSubReg );
            if( m_logSequentialMassDiagnostics != 0 )
            {
              logSequentialMassInventory( domain, "afterFractureFluidState" );
            }
          } );
        }
      }
      if( m_logSequentialMassDiagnostics != 0 )
      {
        logSequentialMassInventory( domain, "afterFractureState" );
      }
    }
  }

  virtual bool checkSequentialConvergence( integer const cycleNumber,
                                           integer const iter,
                                           real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition & domain ) override
  {
    NonlinearSolverParameters const & params = this->getNonlinearSolverParameters();
    if( params.m_subcyclingOption == 0 ||
        params.sequentialConvergenceCriterion() !=
        NonlinearSolverParameters::SequentialConvergenceCriterion::ResidualNorm )
    {
      return Base::checkSequentialConvergence( cycleNumber, iter, time_n, dt, domain );
    }

    GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence, GEOS_FMT( "  Iteration {:2}: outer-loop convergence check", iter + 1 ) );

    auto computeResidualNorm = [&]( auto * solver )
    {
      solver->getLocalMatrix().toViewConstSizes().zero();
      solver->getSystemRhs().zero();
      arrayView1d< real64 > const localRhs = solver->getSystemRhs().open();

      solver->assembleSystem( time_n,
                              dt,
                              domain,
                              solver->getDofManager(),
                              solver->getLocalMatrix().toViewConstSizes(),
                              localRhs );
      solver->applyBoundaryConditions( time_n,
                                       dt,
                                       domain,
                                       solver->getDofManager(),
                                       solver->getLocalMatrix().toViewConstSizes(),
                                       localRhs );
      solver->getSystemRhs().close();

      return solver->calculateResidualNorm( time_n,
                                            dt,
                                            domain,
                                            solver->getDofManager(),
                                            solver->getSystemRhs().values() );
    };

    real64 residualNorm = 0.0;

    // Flow residual is defined on the physical matrix/fracture pressures after
    // restoreCompositePressure.
    real64 const flowResidualNorm = computeResidualNorm( this->flowSolver() );
    residualNorm += flowResidualNorm * flowResidualNorm;

    // Mechanics residual must use the same pressure representation as the actual
    // sequential mechanics step.
    copyFracturePressureToMesh1( domain );
    bool const useCompositePressure = useCompositePressureForSequential();
    if( useCompositePressure )
    {
      swapToCompositePressure( domain );
    }
    Base::updateBulkDensity( domain );

    real64 const solidResidualNorm = computeResidualNorm( this->solidMechanicsSolver() );
    residualNorm += solidResidualNorm * solidResidualNorm;

    if( useCompositePressure )
    {
      restoreCompositePressure( domain );
    }
    Base::updateBulkDensity( domain );

    residualNorm = sqrt( residualNorm );
    GEOS_LOG_LEVEL_RANK_0( logInfo::ResidualNorm,
                           GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
    this->getConvergenceStats().setResidualValue( "R", residualNorm );
    this->updateAndWriteConvergenceStep( time_n, dt, cycleNumber, iter );

    bool const isConverged = ( residualNorm < params.m_newtonTol );
    if( isConverged )
    {
      GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                             GEOS_FMT( "***** The iterative coupling has converged in {} iteration(s) *****", iter + 1 ) );
    }
    return isConverged;
  }

private:

  bool useCompositePressureForSequential() const
  {
    // Intrinsic Sequential now uses the same initialization-time homogenized mechanics state as
    // FullyImplicit. The old runtime composite-pressure path is kept in the source for reference
    // but is disabled because it is not equivalent to the effective-input validation deck.
    return false;
  }

  void relaxSequentialFlowPressure( DomainPartition & domain ) const
  {
    if( m_sequentialPressureRelaxation >= 1.0 ||
        this->getNonlinearSolverParameters().couplingType() !=
        NonlinearSolverParameters::CouplingType::Sequential )
    {
      return;
    }

    GEOS_ERROR_IF( m_sequentialPressureRelaxation <= 0.0,
                   this->getName() << ": " << viewKeyStruct::sequentialPressureRelaxationString()
                                   << " must be in (0, 1]." );

    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    MeshLevel & meshLevel1 = domain.getMeshBody( "mesh1" ).getMeshLevels().template getGroup< MeshLevel >( 0 );
    meshLevel1.getElemManager().forElementSubRegions< CellElementSubRegion >( matrixRegions,
      [&]( localIndex const, CellElementSubRegion & subRegion )
    {
      relaxSubRegionPressureToPreviousOuterIteration( subRegion );
      this->flowSolver()->primarySolver()->updatePorosityAndPermeability( subRegion );
      this->flowSolver()->primarySolver()->updateFluidState( subRegion );
    } );

    MeshLevel & meshLevel2 = domain.getMeshBody( "mesh2" ).getMeshLevels().template getGroup< MeshLevel >( 0 );
    meshLevel2.getElemManager().forElementSubRegions< CellElementSubRegion >( fractureRegions,
      [&]( localIndex const, CellElementSubRegion & subRegion )
    {
      relaxSubRegionPressureToPreviousOuterIteration( subRegion );
      this->flowSolver()->secondarySolver()->updatePorosityAndPermeability( subRegion );
      this->flowSolver()->secondarySolver()->updateFluidState( subRegion );
    } );
  }

  void relaxSubRegionPressureToPreviousOuterIteration( ElementSubRegionBase & subRegion ) const
  {
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
    arrayView1d< real64 > const pressure = subRegion.template getField< fields::flow::pressure >();
    arrayView1d< real64 const > const pressure_k = subRegion.template getField< fields::flow::pressure_k >();
    real64 const relaxation = m_sequentialPressureRelaxation;

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        pressure[k] = pressure_k[k] + relaxation * ( pressure[k] - pressure_k[k] );
      }
    } );
  }

  template< typename FE_TYPE >
  void computeAverageVolumetricStrainIncrement(
    MeshLevel & mesh,
    CellElementSubRegion const & subRegion,
    FE_TYPE const & finiteElement,
    arrayView1d< real64 > const averageVolumetricStrainIncrement ) const
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
    fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const uhat =
      nodeManager.getField< fields::solidMechanics::incrementalDisplacement >();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList().toViewConst();

    typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
    finiteElement::FiniteElementBase::initialize< FE_TYPE >(
      nodeManager, mesh.getEdgeManager(), mesh.getFaceManager(), subRegion, meshData );

    constexpr localIndex numNodesPerElem = FE_TYPE::maxSupportPoints;
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename FE_TYPE::StackVariables feStack;
      finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
      localIndex const numSupportPoints = finiteElement.template numSupportPoints< FE_TYPE >( feStack );

      real64 xLocal[ numNodesPerElem ][ 3 ] = {};
      real64 uhatLocal[ numNodesPerElem ][ 3 ] = {};
      for( localIndex a = 0; a < numSupportPoints; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( integer dim = 0; dim < 3; ++dim )
        {
          xLocal[ a ][ dim ] = X[ nodeIndex ][ dim ];
          uhatLocal[ a ][ dim ] = uhat[ nodeIndex ][ dim ];
        }
      }

      real64 volume = 0.0;
      real64 integratedVolumetricStrainIncrement = 0.0;
      for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
      {
        real64 dNdX[ numNodesPerElem ][ 3 ];
        real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, feStack, dNdX );

        real64 strainIncrement[ 6 ] = {};
        FE_TYPE::symmetricGradient( dNdX, uhatLocal, strainIncrement );
        integratedVolumetricStrainIncrement +=
          detJxW * ( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );
        volume += detJxW;
      }
      averageVolumetricStrainIncrement[k] =
        LvArray::math::abs( volume ) > 0.0 ? integratedVolumetricStrainIncrement / volume : 0.0;
    } );
  }

  void logFimPressureUpdateDiagnostics( DomainPartition & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localSolution,
                                        real64 const scalingFactor )
  {
    if constexpr ( !isMultiphaseFlow )
    {
      GEOS_UNUSED_VAR( domain, dofManager, localSolution, scalingFactor );
      return;
    }
    else
    {
      if( domain.getMeshBodies().numSubGroups() < 2 ) return;

      MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
      MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
      MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
      MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
      using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
      DualFlow & dualFlow = *this->flowSolver();
      string_array const & matrixRegions   = dualFlow.template getReference< string_array >( "matrixRegionList" );
      string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );
      string const dofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );
      globalIndex const rankOffset = dofManager.rankOffset();

      real64 maxAbsRawDpM = 0.0;
      real64 maxAbsRawDpF = 0.0;
      real64 maxAbsScaledDpM = 0.0;
      real64 maxAbsScaledDpF = 0.0;
      real64 minPredPM = LvArray::NumericLimits< real64 >::max;
      real64 minPredPF = LvArray::NumericLimits< real64 >::max;
      real64 minPM = LvArray::NumericLimits< real64 >::max;
      real64 minPF = LvArray::NumericLimits< real64 >::max;

      for( size_t ir = 0; ir < matrixRegions.size(); ++ir )
      {
        ElementRegionBase const & matReg  = meshLevel1.getElemManager().getRegion( matrixRegions[ir] );
        ElementRegionBase const & fracReg = meshLevel2.getElemManager().getRegion( fractureRegions[ir] );

        matReg.template forElementSubRegionsIndex< CellElementSubRegion >(
          [&]( localIndex const isr, CellElementSubRegion const & matSR )
        {
          if( isr >= fracReg.numSubRegions() ) return;
          CellElementSubRegion const & fracSR =
            dynamic_cast< CellElementSubRegion const & >( fracReg.getSubRegion( isr ) );

          arrayView1d< globalIndex const > const dofM =
            matSR.template getReference< array1d< globalIndex > >( dofKey );
          arrayView1d< globalIndex const > const dofF =
            fracSR.template getReference< array1d< globalIndex > >( dofKey );
          arrayView1d< integer const > const ghostM = matSR.ghostRank();
          arrayView1d< integer const > const ghostF = fracSR.ghostRank();
          arrayView1d< real64 const > const pM = matSR.template getField< fields::flow::pressure >();
          arrayView1d< real64 const > const pF = fracSR.template getField< fields::flow::pressure >();
          arrayView1d< real64 const > pressureScalingM;
          arrayView1d< real64 const > pressureScalingF;
          bool const hasLocalPressureScalingM = matSR.hasWrapper( fields::flow::pressureScalingFactor::key() );
          bool const hasLocalPressureScalingF = fracSR.hasWrapper( fields::flow::pressureScalingFactor::key() );
          if( hasLocalPressureScalingM )
          {
            pressureScalingM = matSR.template getField< fields::flow::pressureScalingFactor >();
          }
          if( hasLocalPressureScalingF )
          {
            pressureScalingF = fracSR.template getField< fields::flow::pressureScalingFactor >();
          }

          for( localIndex k = 0; k < matSR.size(); ++k )
          {
            if( ghostM[k] < 0 )
            {
              localIndex const lid = LvArray::integerConversion< localIndex >( dofM[k] - rankOffset );
              if( lid >= 0 && lid < localSolution.size() )
              {
                real64 const raw = localSolution[lid];
                real64 const localScale = hasLocalPressureScalingM ? pressureScalingM[k] : 1.0;
                real64 const scaled = scalingFactor * localScale * raw;
                maxAbsRawDpM = LvArray::math::max( maxAbsRawDpM, LvArray::math::abs( raw ) );
                maxAbsScaledDpM = LvArray::math::max( maxAbsScaledDpM, LvArray::math::abs( scaled ) );
                minPM = LvArray::math::min( minPM, pM[k] );
                minPredPM = LvArray::math::min( minPredPM, pM[k] + scaled );
              }
            }
          }

          for( localIndex kf = 0; kf < fracSR.size(); ++kf )
          {
            if( ghostF[kf] < 0 )
            {
              localIndex const lid = LvArray::integerConversion< localIndex >( dofF[kf] - rankOffset );
              if( lid >= 0 && lid < localSolution.size() )
              {
                real64 const raw = localSolution[lid];
                real64 const localScale = hasLocalPressureScalingF ? pressureScalingF[kf] : 1.0;
                real64 const scaled = scalingFactor * localScale * raw;
                maxAbsRawDpF = LvArray::math::max( maxAbsRawDpF, LvArray::math::abs( raw ) );
                maxAbsScaledDpF = LvArray::math::max( maxAbsScaledDpF, LvArray::math::abs( scaled ) );
                minPF = LvArray::math::min( minPF, pF[kf] );
                minPredPF = LvArray::math::min( minPredPF, pF[kf] + scaled );
              }
            }
          }
        } );
      }

      maxAbsRawDpM = MpiWrapper::max( maxAbsRawDpM );
      maxAbsRawDpF = MpiWrapper::max( maxAbsRawDpF );
      maxAbsScaledDpM = MpiWrapper::max( maxAbsScaledDpM );
      maxAbsScaledDpF = MpiWrapper::max( maxAbsScaledDpF );
      minPM = MpiWrapper::min( minPM );
      minPF = MpiWrapper::min( minPF );
      minPredPM = MpiWrapper::min( minPredPM );
      minPredPF = MpiWrapper::min( minPredPF );

      GEOS_LOG_RANK_0( GEOS_FMT( "{}: FIM pressure-update diagnostics: "
                                 "scale={:.6e}, max|raw dpM|={:.6e}, max|scaled dpM|={:.6e}, "
                                 "min pM {:.6e}->{:.6e}; "
                                 "max|raw dpF|={:.6e}, max|scaled dpF|={:.6e}, min pF {:.6e}->{:.6e}",
                                 this->getName(), scalingFactor,
                                 maxAbsRawDpM, maxAbsScaledDpM, minPM, minPredPM,
                                 maxAbsRawDpF, maxAbsScaledDpF, minPF, minPredPF ) );
    }
  }

  void validateDualContinuumVolumeFractions() const
  {
    real64 const flowFractureVolumeFraction = this->flowSolver()->getFractureVolumeFraction();
    bool const mechanicsHasVolumeFraction = m_fractureVolumeFraction >= 0.0;
    bool const flowHasVolumeFraction = flowFractureVolumeFraction >= 0.0;

    if( mechanicsHasVolumeFraction )
    {
      GEOS_THROW_IF( m_fractureVolumeFraction <= 0.0 || m_fractureVolumeFraction >= 1.0,
                     GEOS_FMT( "{}: fractureVolumeFraction on the dual-continuum poromechanics solver "
                               "must be in (0,1), got {}",
                               this->getName(), m_fractureVolumeFraction ),
                     InputError );
    }

    if( flowHasVolumeFraction )
    {
      GEOS_THROW_IF( flowFractureVolumeFraction <= 0.0 || flowFractureVolumeFraction >= 1.0,
                     GEOS_FMT( "{}: DualContinuumCrossFlow fractureVolumeFraction must be in (0,1), got {}",
                               this->getName(), flowFractureVolumeFraction ),
                     InputError );
    }

    GEOS_THROW_IF( mechanicsHasVolumeFraction != flowHasVolumeFraction,
                   GEOS_FMT( "{}: fractureVolumeFraction must be set consistently in both the "
                             "dual-continuum poromechanics solver and its DualContinuumCrossFlow block. "
                             "The poromechanics value controls effective mechanics/Biot terms, while "
                             "the cross-flow value controls pore-volume scaling and multi-porosity "
                             "storage. Got poromechanics={}, DualContinuumCrossFlow={}.",
                             this->getName(), m_fractureVolumeFraction, flowFractureVolumeFraction ),
                   InputError );

    if( mechanicsHasVolumeFraction )
    {
      GEOS_THROW_IF( std::fabs( m_fractureVolumeFraction - flowFractureVolumeFraction ) > 1e-12,
                     GEOS_FMT( "{}: inconsistent fractureVolumeFraction values. The dual-continuum "
                               "poromechanics solver has {}, but DualContinuumCrossFlow has {}. "
                               "Both must represent the same REV fracture volume fraction v_f.",
                               this->getName(), m_fractureVolumeFraction, flowFractureVolumeFraction ),
                     InputError );
    }

    if( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      bool const fimRequiresVolumeFraction = ( m_useIntrinsicInput != 0 ) || ( m_enableFimCrossStorage != 0 );
      GEOS_THROW_IF( fimRequiresVolumeFraction && !mechanicsHasVolumeFraction,
                     GEOS_FMT( "{}: FullyImplicit dual-continuum poromechanics with useIntrinsicInput={} "
                               "and enableFimCrossStorage={} requires an explicit fractureVolumeFraction "
                               "in both the poromechanics solver and DualContinuumCrossFlow. This is the "
                               "REV volume fraction v_f, not the intrinsic fracture porosity.",
                               this->getName(), m_useIntrinsicInput, m_enableFimCrossStorage ),
                     InputError );
    }
  }

  void validateMaterialInputMode( DomainPartition & domain ) const
  {
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow const & dualFlow = *this->flowSolver();
    real64 const intrKm = dualFlow.getIntrinsicMatrixBulkModulus();
    real64 const intrAm = dualFlow.getIntrinsicMatrixBiot();
    real64 const intrKf = dualFlow.getIntrinsicFractureBulkModulus();
    real64 const intrAf = dualFlow.getIntrinsicFractureBiot();
    bool const hasAnyIntrinsicStorageInput = dualFlow.hasExplicitIntrinsicStorageInput();
    bool const hasEffectiveStorageInput = dualFlow.hasEffectiveStorageInput();
    bool const hasExplicitIntrinsics = intrKm > 0.0 && intrAm >= 0.0 && intrKf > 0.0 && intrAf >= 0.0;

    GEOS_THROW_IF( m_useIntrinsicInput == 0 && hasAnyIntrinsicStorageInput,
                   GEOS_FMT( "{}: useIntrinsicInput=0 is the effective mechanics/storage mode. Do not provide "
                             "internal intrinsic storage reconstruction state in this mode; provide "
                             "effectiveMatrixStorage/effectiveFractureStorage/effectiveCrossStorage instead.",
                             this->getName() ),
                   InputError );
    GEOS_THROW_IF( m_useIntrinsicInput == 0 && !hasEffectiveStorageInput,
                   GEOS_FMT( "{}: useIntrinsicInput=0 is the effective mechanics/storage mode and requires "
                             "direct effective storage inputs on DualContinuumCrossFlow: "
                             "effectiveMatrixStorage, effectiveFractureStorage, and effectiveCrossStorage.",
                             this->getName() ),
                   InputError );
    GEOS_THROW_IF( m_useIntrinsicInput != 0 && hasEffectiveStorageInput,
                   GEOS_FMT( "{}: useIntrinsicInput=1 is the intrinsic mechanics/storage mode. Do not provide "
                             "direct effective storage inputs in this mode; storage is reconstructed "
                             "from the intrinsic constitutive parameters.",
                             this->getName() ),
                   InputError );
    GEOS_THROW_IF( m_useIntrinsicInput != 0 && hasAnyIntrinsicStorageInput && !hasExplicitIntrinsics,
                   GEOS_FMT( "{}: useIntrinsicInput=1 reconstructs storage from intrinsic material models. "
                             "Partial internal intrinsic storage state is not supported.",
                             this->getName() ),
                   InputError );
    if( !hasExplicitIntrinsics )
    {
      return;
    }

    auto relClose = []( real64 const a, real64 const b )
    {
      real64 const scale = LvArray::math::max( real64( 1.0 ), LvArray::math::max( LvArray::math::abs( a ),
                                                                                  LvArray::math::abs( b ) ) );
      return LvArray::math::abs( a - b ) <= 1e-6 * scale;
    };

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    for( size_t i = 0; i < matrixRegions.size(); ++i )
    {
      ElementRegionBase & matRegion = meshLevel1.getElemManager().getRegion( matrixRegions[i] );
      ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fractureRegions[i] );

      matRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
      {
        if( subRegIdx >= fracRegion.numSubRegions() ) return;
        CellElementSubRegion & fracSubReg =
          dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

        string const & matSolidName = matSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & matSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
        string const & fracSolidName = fracSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & fracSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );

        arrayView1d< real64 const > const K_m = matSolid.getBulkModulus();
        arrayView1d< real64 const > const A_m = matSolid.getBiotCoefficient();
        arrayView1d< real64 const > const K_f = fracSolid.getBulkModulus();
        arrayView1d< real64 const > const A_f = fracSolid.getBiotCoefficient();
        if( matSubReg.size() == 0 || fracSubReg.size() == 0 ) return;

        bool const looksIntrinsic =
          relClose( K_m[0], intrKm ) && relClose( A_m[0], intrAm ) &&
          relClose( K_f[0], intrKf ) && relClose( A_f[0], intrAf );

        GEOS_THROW_IF( m_useIntrinsicInput == 0 && looksIntrinsic,
                       GEOS_FMT( "{}: useIntrinsicInput=0 means the constitutive K/Biot values are "
                                 "already EFFECTIVE and will not be homogenized. However, subregion "
                                 "'{}'/'{}' has constitutive values matching the intrinsic* values on "
                                 "DualContinuumCrossFlow. If the XML is giving intrinsic material "
                                 "parameters, set useIntrinsicInput=\"1\"; otherwise provide effective "
                                 "Kbar/Gbar/abar in the constitutive block.",
                                 this->getName(), matSubReg.getName(), fracSubReg.getName() ),
                       InputError );
        GEOS_THROW_IF( m_useIntrinsicInput != 0 && !looksIntrinsic,
                       GEOS_FMT( "{}: useIntrinsicInput=1 means the constitutive K/Biot values are "
                                 "INTRINSIC material parameters. They must match the intrinsic* values "
                                 "on DualContinuumCrossFlow when those intrinsic* values are specified. "
                                 "Subregion '{}'/'{}' does not match; either correct the intrinsic "
                                 "material parameters or set useIntrinsicInput=\"0\" for effective input.",
                                 this->getName(), matSubReg.getName(), fracSubReg.getName() ),
                       InputError );
      } );
    }
  }

  // Copy p_f from mesh2 to mesh1 fracture wrappers directly (bypasses
  // DofManager, needed in sequential mode where mapFractureDataToMatrix
  // is not called).
  void copyFracturePressureToMesh1( DomainPartition & domain )
  {
    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    for( size_t i = 0; i < matrixRegions.size(); ++i )
    {
      ElementRegionBase & matrixRegion = meshLevel1.getElemManager().getRegion( matrixRegions[i] );
      ElementRegionBase & fractureRegion = meshLevel2.getElemManager().getRegion( fractureRegions[i] );

      matrixRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matrixSubRegion )
      {
        if( subRegIdx >= fractureRegion.numSubRegions() ) return;
        CellElementSubRegion & fractureSubRegion =
          dynamic_cast< CellElementSubRegion & >( fractureRegion.getSubRegion( subRegIdx ) );

        arrayView1d< real64 > const p_f_wrapper =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fracturePressureString() );
        arrayView1d< real64 const > const p_f_mesh2 = fractureSubRegion.getField< fields::flow::pressure >();

        arrayView1d< real64 > const alpha_f_wrapper =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 > const fractureMechanicsScale =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );
        string const & fracturePorousName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & fractureSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fractureSubRegion, fracturePorousName );
        arrayView1d< real64 const > const alpha_f_mesh2 = fractureSolid.getBiotCoefficient();
        arrayView1d< localIndex const > const matrixToFracture =
          matrixSubRegion.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fractureSubRegion.size(),
                         "Invalid dual-continuum copyFracturePressureToMesh1 connectivity for matrix subregion "
                         << matrixSubRegion.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fractureSubRegion.getName()
                         << " size " << fractureSubRegion.size() );
          p_f_wrapper[ k ] = p_f_mesh2[ kf ];
          alpha_f_wrapper[ k ] = alpha_f_mesh2[ kf ];
          fractureMechanicsScale[ k ] = 1.0;
        }
      } );
    }

  }

  // useIntrinsicInput=1 path: homogenize the INTRINSIC matrix/fracture constitutive parameters
  // into the effective medium ONCE at initialization for both FullyImplicit and Sequential.
  // This writes the effective values permanently so the mechanics kernel and K_upf use the homogenized
  // response, while the intrinsics are pushed to DualContinuumCrossFlow so the multi-porosity
  // storage (Step 4b) keeps using the intrinsic 1/M_bar. No-op when the flag is off.
  //   Kbar = (v_m/K_m + v_f/K_f)^-1 ,  Gbar = (v_m/G_m + v_f/G_f)^-1   (Reuss)
  //   abar_i = Kbar * v_i * alpha_i / K_i
  // Safe wrt other physics: only runs for SinglePhaseDualContinuumPoromechanics and only writes
  // this solver's own matrix/fracture constitutive instances.
  // Automatically set the matrix initial EFFECTIVE stress to the dual-continuum total Biot stress
  //   sigma'_0 = (alpha_m * p_m,0 + alpha_f * p_f,0) * I   (isotropic, GEOS sign convention)
  // so the initial TOTAL stress sigma = sigma' - alpha_m*p_m - alpha_f*p_f = 0 is in equilibrium and
  // the first Newton iteration starts balanced (Rsolid ~ 0), the same way single-porosity poromechanics
  // does (there the matrix Biot term is incremental from a zero initial effective stress). Without this
  // the co-located dual mesh must set matrixSolid_stress by hand. Assumes isotropic initial stress and
  // no external load; for anisotropic/loaded states set autoInitializeStress=0 and prescribe manually.
  void autoInitializeEffectiveStress( DomainPartition & domain )
  {
    if( m_autoInitializeStress == 0 ) return;
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions   = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    bool anyInitialized = false;
    for( size_t i = 0; i < matrixRegions.size(); ++i )
    {
      ElementRegionBase & matRegion  = meshLevel1.getElemManager().getRegion( matrixRegions[i] );
      ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fractureRegions[i] );

      matRegion.template forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
      {
        if( subRegIdx >= fracRegion.numSubRegions() ) return;
        CellElementSubRegion & fracSubReg =
          dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

        arrayView1d< real64 const > const p_m = matSubReg.template getField< fields::flow::pressure >();
        arrayView1d< real64 const > const p_f = fracSubReg.template getField< fields::flow::pressure >();

        string const & matSolidName = matSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
        arrayView1d< real64 const > const alpha_m = matSolid.getBiotCoefficient();

        string const & fracSolidName = fracSubReg.getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
        arrayView1d< real64 const > const alpha_f = fracSolid.getBiotCoefficient();
        arrayView1d< localIndex const > const matrixToFracture =
          matSubReg.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        // matrix elastic solid (carries the stress tensor used by the mechanics kernel)
        string const & matElasticName = matSolid.getReference< string >(
          constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
        constitutive::SolidBase & matElastic =
          matSolid.getParent().template getGroup< constitutive::SolidBase >( matElasticName );
        arrayView3d< real64, solid::STRESS_USD > const stress = matElastic.getStress();

        // Respect a user-prescribed initial stress: if matrixSolid_stress was set by a
        // FieldSpecification (already applied at this point), do NOT clobber it. This keeps
        // anisotropic / geostatic / gravity-loaded initial states under user control.
        real64 maxAbsStress = 0.0;
        for( localIndex k = 0; k < matSubReg.size(); ++k )
          for( localIndex q = 0; q < stress.size( 1 ); ++q )
            for( integer c = 0; c < 6; ++c )
              maxAbsStress = LvArray::math::max( maxAbsStress, LvArray::math::abs( stress[k][q][c] ) );
        if( maxAbsStress > 1.0 )  // ~0 Pa tolerance; nonzero => user prescribed it
        {
          GEOS_LOG_RANK_0( this->getName() << ": autoInitializeStress=1 but a nonzero initial "
                           "matrix stress is already prescribed; leaving it unchanged." );
          return;
        }

        // The fracture pressure enters the matrix momentum balance only through K_upf
        // (assembleFractureMechanicsCoupling). The fracture Biot coefficient stored on the
        // constitutive model is the effective dual-continuum coefficient, so no additional
        // volume-fraction scaling is applied here.
        bool const includeFracture = ( m_enableFractureMechanicsCoupling != 0 );
        for( localIndex k = 0; k < matSubReg.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fracSubReg.size(),
                         "Invalid dual-continuum autoInitializeEffectiveStress connectivity for matrix subregion "
                         << matSubReg.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fracSubReg.getName()
                         << " size " << fracSubReg.size() );
          real64 const sigma0 =
            alpha_m[k] * p_m[k] + ( includeFracture ? alpha_f[kf] * p_f[kf] : 0.0 );
          for( localIndex q = 0; q < stress.size( 1 ); ++q )
          {
            stress[k][q][0] = sigma0; stress[k][q][1] = sigma0; stress[k][q][2] = sigma0;
            stress[k][q][3] = 0.0;    stress[k][q][4] = 0.0;    stress[k][q][5] = 0.0;
          }
        }
        matElastic.saveConvergedState();  // copy to old stress so the increment starts from sigma0
        anyInitialized = true;
      } );
    }
    if( anyInitialized )
      GEOS_LOG_RANK_0( this->getName() << ": auto-initialized matrix effective stress to "
                       "sigma' = alpha_m*p_m + alpha_f*p_f for dual-continuum equilibrium." );
  }

  // Multi-porosity (Mehrabian 2014) cross-storage for the compositional path.
  // Each continuum's fluid content is coupled to BOTH pressures through the shared solid skeleton:
  //   d(phi_i)/d(p_i) gets the effective-medium value Sbar_ii, and a NEW off-diagonal Sbar_ij
  //   couples continuum i to the other continuum's pressure p_j. The compositional accumulation only
  //   put the intrinsic skeleton diagonal (alpha-phi)/Ks on each continuum, and NO off-diagonal, so
  //   here we add the correction. Per component c the pore-volume change V*(Sbar_ii*dp_i + Sbar_ij*dp_j)
  //   carries component mass compDens_c, so:
  //     row(i,c) += compDens_i,c * V * ( corrDiag_i * dp_i + corrOff * dp_j ).
  //   Jacobian: pressure columns p_i (corrDiag_i) and p_j (corrOff). The compDens self-derivative
  //   (~ proportional to dp) is omitted in this first version (small; the nonlinear safeguards absorb
  //   it). The useTotalMassEquation row transform is replicated. Skeleton-only storage (no fluid
  //   compressibility term) because the compositional flash already carries fluid compressibility.
  //   Disabled when v_f<=0 or m_enableFimCrossStorage==0; single-porosity limit (v_f->0) => zero.
  void assembleFimCrossStorageMultiphase( DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();

    real64 const v_f = dualFlow.getFractureVolumeFraction();
    if( !( v_f > 0.0 ) ) return;   // single-porosity limit / unset => no cross-storage
    real64 const v_m = 1.0 - v_f;
    real64 const offScale = dualFlow.getCrossStorageOffDiagScale();
    bool const useEffectiveStorageInput = dualFlow.hasEffectiveStorageInput();
    real64 const effectiveStorageM = dualFlow.getEffectiveMatrixStorage();
    real64 const effectiveStorageF = dualFlow.getEffectiveFractureStorage();
    real64 const effectiveStorageMF = dualFlow.getEffectiveCrossStorage();
    GEOS_ERROR_IF( useEffectiveStorageInput && !( effectiveStorageM > 0.0 && effectiveStorageF > 0.0 ),
                   "DualContinuumCrossFlow effective storage input requires positive "
                   "effectiveMatrixStorage and effectiveFractureStorage." );

    // Intrinsic (true physical) Biot + drained bulk modulus for the M_bar storage; fall back to the
    // material value when unset (<0).
    real64 const intrMatA  = dualFlow.getIntrinsicMatrixBiot();
    real64 const intrMatK  = dualFlow.getIntrinsicMatrixBulkModulus();
    real64 const intrFracA = dualFlow.getIntrinsicFractureBiot();
    real64 const intrFracK = dualFlow.getIntrinsicFractureBulkModulus();

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    string const dofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );
    globalIndex const rankOffset = dofManager.rankOffset();
    integer const numComp      = this->flowSolver()->primarySolver()->numFluidComponents();
    integer const useTotalMass = this->flowSolver()->primarySolver()->useTotalMassEquation();
    constexpr integer maxNumComp = 16;

    string_array const & matrixRegions   = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    real64 maxAbsStorageM = 0.0;
    real64 maxAbsStorageF = 0.0;
    real64 maxAbsCorrDiagM = 0.0;
    real64 maxAbsCorrDiagF = 0.0;
    real64 maxAbsCorrOff = 0.0;
    real64 maxAbsDpM = 0.0;
    real64 maxAbsDpF = 0.0;
    real64 sumAbsStorageM = 0.0;
    real64 sumAbsStorageF = 0.0;
    real64 minPM = LvArray::NumericLimits< real64 >::max;
    real64 minPF = LvArray::NumericLimits< real64 >::max;
    real64 maxPM = -LvArray::NumericLimits< real64 >::max;
    real64 maxPF = -LvArray::NumericLimits< real64 >::max;

    for( size_t ir = 0; ir < matrixRegions.size(); ++ir )
    {
      ElementRegionBase const & matReg  = meshLevel1.getElemManager().getRegion( matrixRegions[ir] );
      ElementRegionBase const & fracReg = meshLevel2.getElemManager().getRegion( fractureRegions[ir] );

      matReg.template forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const isr, CellElementSubRegion const & matSR )
      {
        if( isr >= fracReg.numSubRegions() ) return;
        CellElementSubRegion const & fracSR =
          dynamic_cast< CellElementSubRegion const & >( fracReg.getSubRegion( isr ) );

        arrayView1d< globalIndex const > const dofM = matSR.template getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< globalIndex const > const dofF = fracSR.template getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< integer const > const ghostM = matSR.ghostRank();
        arrayView1d< integer const > const ghostF = fracSR.ghostRank();
        arrayView1d< real64 const > const pM  = matSR.template getField< fields::flow::pressure >();
        arrayView1d< real64 const > const pMn = matSR.template getField< fields::flow::pressure_n >();
        arrayView1d< real64 const > const pF  = fracSR.template getField< fields::flow::pressure >();
        arrayView1d< real64 const > const pFn = fracSR.template getField< fields::flow::pressure_n >();
        arrayView1d< real64 const > const volM = matSR.getElementVolume();
        arrayView1d< real64 const > const volF = fracSR.getElementVolume();
        auto const compDensM = matSR.template getField< fields::flow::globalCompDensity >();
        auto const compDensF = fracSR.template getField< fields::flow::globalCompDensity >();

        string const & solidMn = matSR.template getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        string const & solidFn = fracSR.template getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & csM = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSR, solidMn );
        constitutive::CoupledSolidBase const & csF = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSR, solidFn );
        arrayView1d< real64 const > const Km = csM.getBulkModulus();
        arrayView1d< real64 const > const Kf = csF.getBulkModulus();
        arrayView1d< real64 const > const aM = csM.getBiotCoefficient();
        arrayView1d< real64 const > const aF = csF.getBiotCoefficient();
        arrayView1d< real64 const > const phiM = csM.getReferencePorosity();
        arrayView1d< real64 const > const phiF = csF.getReferencePorosity();
        arrayView1d< real64 const > const KsM =
          dynamic_cast< constitutive::BiotPorosity const & >( csM.getBasePorosityModel() ).getGrainBulkModulus();
        arrayView1d< real64 const > const KsF =
          dynamic_cast< constitutive::BiotPorosity const & >( csF.getBasePorosityModel() ).getGrainBulkModulus();
        arrayView1d< localIndex const > const matrixToFracture =
          matSR.template getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        for( localIndex k = 0; k < matSR.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fracSR.size(),
                         "Invalid dual-continuum cross-storage connectivity for matrix subregion "
                         << matSR.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fracSR.getName()
                         << " size " << fracSR.size() );
          real64 const aMI = ( intrMatA  > 0.0 ) ? intrMatA  : aM[k];
          real64 const aFI = ( intrFracA > 0.0 ) ? intrFracA : aF[kf];
          real64 const KmI = ( intrMatK  > 0.0 ) ? intrMatK  : Km[k];
          real64 const KfI = ( intrFracK > 0.0 ) ? intrFracK : Kf[kf];

          real64 SbarMM = effectiveStorageM;
          real64 SbarFF = effectiveStorageF;
          real64 corrOff = effectiveStorageMF;
          if( !useEffectiveStorageInput )
          {
            real64 const Kbar = 1.0 / ( v_m/KmI + v_f/KfI );
            real64 const abm  = Kbar*v_m*aMI/KmI;
            real64 const abf  = Kbar*v_f*aFI/KfI;
            // phiM/phiF have already been scaled by dual-continuum netToGross, so they are REV pore
            // fractions (v_i * phi_i). The intrinsic skeleton storage needs the continuum-local
            // porosity phi_i; recover it here before applying Mehrabian's M_bar formula.
            real64 const phiMI = phiM[k] / v_m;
            real64 const phiFI = phiF[kf] / v_f;
            // skeleton-only intrinsic storage (no fluid compressibility; multiphase flash carries it)
            real64 const invMmI = (aMI-phiMI)/KsM[k];
            real64 const invMfI = (aFI-phiFI)/KsF[kf];
            SbarMM = v_m*(invMmI + aMI*aMI/KmI) - abm*abm/Kbar;
            SbarFF = v_f*(invMfI + aFI*aFI/KfI) - abf*abf/Kbar;
            corrOff = -abm*abf/Kbar;
          }
          // what the compositional accumulation already put on the diagonal: skeleton (alpha-phi)/Ks
          real64 const invMmMat = (aM[k]-phiM[k])/KsM[k];
          real64 const invMfMat = (aF[kf]-phiF[kf])/KsF[kf];
          real64 const corrDiagM = SbarMM - invMmMat;
          real64 const corrDiagF = SbarFF - invMfMat;
          corrOff *= offScale;  // symmetric off-diagonal

          real64 const dpM = pM[k]-pMn[k];
          real64 const dpF = pF[kf]-pFn[kf];

          maxAbsCorrDiagM = LvArray::math::max( maxAbsCorrDiagM, LvArray::math::abs( corrDiagM ) );
          maxAbsCorrDiagF = LvArray::math::max( maxAbsCorrDiagF, LvArray::math::abs( corrDiagF ) );
          maxAbsCorrOff = LvArray::math::max( maxAbsCorrOff, LvArray::math::abs( corrOff ) );
          maxAbsDpM = LvArray::math::max( maxAbsDpM, LvArray::math::abs( dpM ) );
          maxAbsDpF = LvArray::math::max( maxAbsDpF, LvArray::math::abs( dpF ) );
          minPM = LvArray::math::min( minPM, pM[k] );
          minPF = LvArray::math::min( minPF, pF[kf] );
          maxPM = LvArray::math::max( maxPM, pM[k] );
          maxPF = LvArray::math::max( maxPF, pF[kf] );

          // per-component coefficients (compDens_c), with optional total-mass row transform
          real64 coeffM[ maxNumComp ] = {};
          real64 coeffF[ maxNumComp ] = {};
          if( useTotalMass )
          {
            real64 sM = 0.0, sF = 0.0;
            for( integer c = 0; c < numComp; ++c ) { sM += compDensM[k][c]; sF += compDensF[kf][c]; }
            coeffM[0] = sM; coeffF[0] = sF;
            for( integer i = 1; i < numComp; ++i ) { coeffM[i] = compDensM[k][i-1]; coeffF[i] = compDensF[kf][i-1]; }
          }
          else
          {
            for( integer c = 0; c < numComp; ++c ) { coeffM[c] = compDensM[k][c]; coeffF[c] = compDensF[kf][c]; }
          }

          // matrix continuum mass rows: vs p_m (diag) and p_f (off)
          if( ghostM[k] < 0 )
          {
            real64 const V = volM[k];
            real64 const storageIncrement = V * ( corrDiagM * dpM + corrOff * dpF );
            maxAbsStorageM = LvArray::math::max( maxAbsStorageM, LvArray::math::abs( storageIncrement ) );
            sumAbsStorageM += LvArray::math::abs( storageIncrement );
            globalIndex const colP_self  = dofM[k];       // matrix pressure DOF (offset 0)
            globalIndex const colP_other = dofF[kf];      // fracture pressure DOF (offset 0)
            for( integer i = 0; i < numComp; ++i )
            {
              localIndex const row = LvArray::integerConversion< localIndex >( dofM[k] + i - rankOffset );
              if( row < 0 || row >= localMatrix.numRows() ) continue;
              localRhs[row] += coeffM[i] * storageIncrement;
              globalIndex cols[2] = { colP_self, colP_other };
              real64 vals[2] = { coeffM[i]*V*corrDiagM, coeffM[i]*V*corrOff };
              localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );

              if( useTotalMass && i == 0 )
              {
                for( integer jc = 0; jc < numComp; ++jc )
                {
                  globalIndex const colC = dofM[k] + jc + 1;
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &colC, &storageIncrement, 1 );
                }
              }
              else
              {
                integer const compIndex = useTotalMass ? i - 1 : i;
                globalIndex const colC = dofM[k] + compIndex + 1;
                localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &colC, &storageIncrement, 1 );
              }
            }
          }
          // fracture continuum mass rows: vs p_f (diag) and p_m (off)
          if( ghostF[kf] < 0 )
          {
            real64 const V = volF[kf];
            real64 const storageIncrement = V * ( corrDiagF * dpF + corrOff * dpM );
            maxAbsStorageF = LvArray::math::max( maxAbsStorageF, LvArray::math::abs( storageIncrement ) );
            sumAbsStorageF += LvArray::math::abs( storageIncrement );
            globalIndex const colP_self  = dofF[kf];      // fracture pressure DOF (offset 0)
            globalIndex const colP_other = dofM[k];       // matrix pressure DOF (offset 0)
            for( integer i = 0; i < numComp; ++i )
            {
              localIndex const row = LvArray::integerConversion< localIndex >( dofF[kf] + i - rankOffset );
              if( row < 0 || row >= localMatrix.numRows() ) continue;
              localRhs[row] += coeffF[i] * storageIncrement;
              globalIndex cols[2] = { colP_self, colP_other };
              real64 vals[2] = { coeffF[i]*V*corrDiagF, coeffF[i]*V*corrOff };
              localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );

              if( useTotalMass && i == 0 )
              {
                for( integer jc = 0; jc < numComp; ++jc )
                {
                  globalIndex const colC = dofF[kf] + jc + 1;
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &colC, &storageIncrement, 1 );
                }
              }
              else
              {
                integer const compIndex = useTotalMass ? i - 1 : i;
                globalIndex const colC = dofF[kf] + compIndex + 1;
                localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &colC, &storageIncrement, 1 );
              }
            }
          }
        }
      } );
    }

    if( m_logFimCouplingDiagnostics != 0 )
    {
      maxAbsStorageM = MpiWrapper::max( maxAbsStorageM );
      maxAbsStorageF = MpiWrapper::max( maxAbsStorageF );
      maxAbsCorrDiagM = MpiWrapper::max( maxAbsCorrDiagM );
      maxAbsCorrDiagF = MpiWrapper::max( maxAbsCorrDiagF );
      maxAbsCorrOff = MpiWrapper::max( maxAbsCorrOff );
      maxAbsDpM = MpiWrapper::max( maxAbsDpM );
      maxAbsDpF = MpiWrapper::max( maxAbsDpF );
      sumAbsStorageM = MpiWrapper::sum( sumAbsStorageM );
      sumAbsStorageF = MpiWrapper::sum( sumAbsStorageF );
      minPM = MpiWrapper::min( minPM );
      minPF = MpiWrapper::min( minPF );
      maxPM = MpiWrapper::max( maxPM );
      maxPF = MpiWrapper::max( maxPF );

      GEOS_LOG_RANK_0( GEOS_FMT( "{}: FIM cross-storage diagnostics: "
                                 "max|storage_m|={:.6e}, sum|storage_m|={:.6e}, "
                                 "max|storage_f|={:.6e}, sum|storage_f|={:.6e}, "
                                 "max|corrDiagM|={:.6e}, max|corrDiagF|={:.6e}, max|corrOff|={:.6e}, "
                                 "max|dpM|={:.6e}, max|dpF|={:.6e}, "
                                 "pM=[{:.6e},{:.6e}], pF=[{:.6e},{:.6e}]",
                                 this->getName(),
                                 maxAbsStorageM, sumAbsStorageM,
                                 maxAbsStorageF, sumAbsStorageF,
                                 maxAbsCorrDiagM, maxAbsCorrDiagF, maxAbsCorrOff,
                                 maxAbsDpM, maxAbsDpF,
                                 minPM, maxPM, minPF, maxPF ) );
    }
  }

  void enforceFractureCompositionalVolumeClosure( DomainPartition & domain,
                                                  bool const saveConvergedState )
  {
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    integer const numComp = this->flowSolver()->secondarySolver()->numFluidComponents();
    integer const numPhase = this->flowSolver()->secondarySolver()->numFluidPhases();

    for( string const & fractureRegionName : fractureRegions )
    {
      ElementRegionBase & fracReg = meshLevel2.getElemManager().getRegion( fractureRegionName );
      fracReg.template forElementSubRegions< CellElementSubRegion >(
        [&]( CellElementSubRegion & fracSR )
      {
        if( !fracSR.hasWrapper( fields::flow::globalCompDensity::key() ) )
        {
          return;
        }

        arrayView2d< real64, compflow::USD_COMP > const compDens =
          fracSR.template getField< fields::flow::globalCompDensity >();
        arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
          fracSR.template getField< fields::flow::phaseVolumeFraction >();
        arrayView1d< integer const > const ghostRank = fracSR.ghostRank();

        constexpr integer maxClosureIter = 6;
        bool updatedFluidState = false;
        for( integer closureIter = 0; closureIter < maxClosureIter; ++closureIter )
        {
          RAJA::ReduceMax< parallelDeviceReduce, real64 > maxVolumeClosureError( 0.0 );
          forAll< parallelDevicePolicy<> >( fracSR.size(),
            [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            if( ghostRank[k] >= 0 )
            {
              return;
            }

            real64 phaseVolFracSum = 0.0;
            for( integer ip = 0; ip < numPhase; ++ip )
            {
              phaseVolFracSum += phaseVolFrac[k][ip];
            }

            real64 const volumeClosureError = LvArray::math::abs( 1.0 - phaseVolFracSum );
            maxVolumeClosureError.max( volumeClosureError );

            if( phaseVolFracSum > 1.0e-12 && volumeClosureError > 1.0e-10 )
            {
              real64 const densityScale = 1.0 / phaseVolFracSum;
              for( integer ic = 0; ic < numComp; ++ic )
              {
                compDens[k][ic] *= densityScale;
              }
            }
          } );

          if( maxVolumeClosureError.get() <= 1.0e-10 )
          {
            break;
          }

          this->flowSolver()->secondarySolver()->updateFluidState( fracSR );
          updatedFluidState = true;
        }
        if( saveConvergedState )
        {
          if( !updatedFluidState )
          {
            this->flowSolver()->secondarySolver()->updateFluidState( fracSR );
          }
          this->flowSolver()->secondarySolver()->saveConvergedState( fracSR );
        }
      } );
    }

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::globalCompDensity::key() },
                                     fractureRegions );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, meshLevel2, domain.getNeighbors(), true );

    for( string const & fractureRegionName : fractureRegions )
    {
      ElementRegionBase & fracReg = meshLevel2.getElemManager().getRegion( fractureRegionName );
      fracReg.template forElementSubRegions< CellElementSubRegion >(
        [&]( CellElementSubRegion & fracSR )
      {
        this->flowSolver()->secondarySolver()->updateFluidState( fracSR );
        if( saveConvergedState )
        {
          this->flowSolver()->secondarySolver()->saveConvergedState( fracSR );
        }
      } );
    }
  }

  void assembleFractureVolumeClosureJacobianMultiphase( DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix )
  {
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    string const dofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );
    globalIndex const rankOffset = dofManager.rankOffset();
    integer const numComp = this->flowSolver()->secondarySolver()->numFluidComponents();
    integer const numPhase = this->flowSolver()->secondarySolver()->numFluidPhases();
    constexpr integer maxNumComp = 16;
    GEOS_THROW_IF( numComp > maxNumComp,
                   GEOS_FMT( "{}: fracture volume-closure Jacobian supports at most {} components, got {}",
                             this->getName(), maxNumComp, numComp ),
                   std::runtime_error );

    for( string const & fractureRegionName : fractureRegions )
    {
      ElementRegionBase & fracReg = meshLevel2.getElemManager().getRegion( fractureRegionName );
      fracReg.template forElementSubRegions< CellElementSubRegion >(
        [&]( CellElementSubRegion & fracSR )
      {
        arrayView1d< globalIndex const > const dofNumber =
          fracSR.template getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< integer const > const ghostRank = fracSR.ghostRank();
        arrayView1d< real64 const > const volume = fracSR.getElementVolume();
        arrayView2d< real64 const, compflow::USD_COMP > const compDens =
          fracSR.template getField< fields::flow::globalCompDensity >();
        arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
          fracSR.template getField< fields::flow::phaseVolumeFraction >();

        string const & solidName =
          fracSR.template getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & solid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSR, solidName );
        arrayView2d< real64 const > const porosity_n = solid.getPorosity_n();

        forAll< parallelDevicePolicy<> >( fracSR.size(),
          [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          if( ghostRank[k] >= 0 )
          {
            return;
          }

          real64 totalCompDensity = 0.0;
          for( integer ic = 0; ic < numComp; ++ic )
          {
            totalCompDensity += compDens[k][ic];
          }

          if( totalCompDensity <= 1.0e-12 )
          {
            return;
          }

          real64 phaseVolFracSum = 0.0;
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            phaseVolFracSum += phaseVolFrac[k][ip];
          }

          real64 const oldPoreVolume = volume[k] * porosity_n[k][0];
          real64 const dVolumeConstraint_dCompDensity =
            -oldPoreVolume * phaseVolFracSum / totalCompDensity;

          globalIndex cols[maxNumComp]{};
          real64 vals[maxNumComp]{};
          for( integer ic = 0; ic < numComp; ++ic )
          {
            cols[ic] = dofNumber[k] + ic + 1;
            vals[ic] = dVolumeConstraint_dCompDensity;
          }

          localIndex const row =
            LvArray::integerConversion< localIndex >( dofNumber[k] + numComp - rankOffset );
          if( row >= 0 && row < localMatrix.numRows() )
          {
            localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( row, cols, vals, numComp );
          }
        } );
      } );
    }
  }

  void computeEffectiveFromIntrinsic( DomainPartition & domain )
  {
    if( m_useIntrinsicInput == 0 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    real64 intrKm = -1.0, intrAm = -1.0, intrKf = -1.0, intrAf = -1.0;

    for( size_t i = 0; i < matrixRegions.size(); ++i )
    {
      ElementRegionBase & matRegion = meshLevel1.getElemManager().getRegion( matrixRegions[i] );
      ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fractureRegions[i] );

      matRegion.template forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
      {
        if( subRegIdx >= fracRegion.numSubRegions() ) return;
        CellElementSubRegion & fracSubReg =
          dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

        // matrix constitutive (intrinsic in → effective out)
        string const & matSolidName = matSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
        string const & matElasticName = matSolid.getReference< string >(
          constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
        constitutive::ElasticIsotropic & matElastic =
          matSolid.getParent().template getGroup< constitutive::ElasticIsotropic >( matElasticName );
        arrayView1d< real64 > const K_m = matElastic.bulkModulus();
        arrayView1d< real64 > const G_m = matElastic.shearModulus();
        constitutive::BiotPorosity & matPor =
          dynamic_cast< constitutive::BiotPorosity & >( matSolid.getBasePorosityModel() );
        arrayView1d< real64 > const alpha_m = matPor.getBiotCoefficientWritable();
        typename constitutive::BiotPorosity::KernelWrapper matPorUpdates = matPor.createKernelUpdates();

        // fracture constitutive
        string const & fracSolidName = fracSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
        string const & fracElasticName = fracSolid.getReference< string >(
          constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
        constitutive::ElasticIsotropic & fracElastic =
          fracSolid.getParent().template getGroup< constitutive::ElasticIsotropic >( fracElasticName );
        arrayView1d< real64 > const K_f = fracElastic.bulkModulus();
        arrayView1d< real64 > const G_f = fracElastic.shearModulus();
        constitutive::BiotPorosity & fracPor =
          dynamic_cast< constitutive::BiotPorosity & >( fracSolid.getBasePorosityModel() );
        arrayView1d< real64 > const alpha_f = fracPor.getBiotCoefficientWritable();
        typename constitutive::BiotPorosity::KernelWrapper fracPorUpdates = fracPor.createKernelUpdates();
        arrayView1d< real64 > const effectiveFractureBiot =
          matSubReg.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );

        bool const useMeshVolumes = ( m_fractureVolumeFraction < 0.0 );
        arrayView1d< real64 const > V_m, V_f;
        if( useMeshVolumes ) { V_m = matSubReg.getElementVolume(); V_f = fracSubReg.getElementVolume(); }
        arrayView1d< localIndex const > const matrixToFracture =
          matSubReg.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        for( localIndex k = 0; k < matSubReg.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fracSubReg.size(),
                         "Invalid dual-continuum computeEffectiveFromIntrinsic connectivity for matrix subregion "
                         << matSubReg.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fracSubReg.getName()
                         << " size " << fracSubReg.size() );
          real64 const KmI = K_m[k];
          real64 const GmI = G_m[k];
          real64 const aMI = alpha_m[k];
          real64 const KfI = K_f[kf];
          real64 const GfI = G_f[kf];
          real64 const aFI = alpha_f[kf];

          // capture a representative intrinsic set (homogeneous material) for the storage params
          if( intrKm < 0.0 ) { intrKm = KmI; intrAm = aMI; intrKf = KfI; intrAf = aFI; }

          real64 const v_f = useMeshVolumes ? V_f[kf] / ( V_m[k] + V_f[kf] ) : m_fractureVolumeFraction;
          real64 const v_m = 1.0 - v_f;

          real64 const Kbar = 1.0 / ( v_m / KmI + v_f / KfI );
          real64 const Gbar = 1.0 / ( v_m / GmI + v_f / GfI );
          real64 const abar_m = Kbar * v_m * aMI / KmI;
          real64 const abar_f = Kbar * v_f * aFI / KfI;

          K_m[k] = Kbar;  G_m[k] = Gbar;   // matrix mechanics → effective drained moduli
          matPorUpdates.updateBiotCoefficientAndAssignModuli( k, Kbar, Gbar );
          alpha_m[k] = abar_m;             // matrix Biot → effective (kernel K_upm)
          K_f[kf] = Kbar;  G_f[kf] = Gbar;  // fracture flow/porosity base state → effective
          fracPorUpdates.updateBiotCoefficientAndAssignModuli( kf, Kbar, Gbar );
          alpha_f[kf] = abar_f;
          effectiveFractureBiot[k] = abar_f; // fracture mechanics Biot → effective (K_upf / K_pfu)
        }
      } );
    }

    // Push intrinsics to the multi-porosity storage (Step 4b uses intrinsic 1/M_bar).
    if( intrKm > 0.0 )
    {
      dualFlow.setIntrinsicMatrixBulkModulus( intrKm );
      dualFlow.setIntrinsicMatrixBiot( intrAm );
      dualFlow.setIntrinsicFractureBulkModulus( intrKf );
      dualFlow.setIntrinsicFractureBiot( intrAf );
      GEOS_LOG_RANK_0( this->getName() << ": useIntrinsicInput → homogenized effective medium from intrinsics ("
                       << "K_m=" << intrKm << ", alpha_m=" << intrAm
                       << ", K_f=" << intrKf << ", alpha_f=" << intrAf << ")" );
    }
  }

  // Replace matrix pressure (and pressure_n) on mesh1 with equivalent single-porosity
  // values so the mechanics kernel — which uses the matrix constitutive model (K_m, α_m)
  // — produces the correct dual-porosity composite response:
  //
  //   K_eff   = (v_m/K_m + v_f/K_f)^{-1}
  //   α_eff_i = K_eff · v_i · α_i / K_i
  //   p_eq    = (α_eff_m·p_m + α_eff_f·p_f) / α_m
  //
  // Stress:          σ_biot = -α_m·p_eq  = -(α_eff_m·p_m + α_eff_f·p_f)           ✓
  // Stress increment: α_m·(p_eq-p_eq_n)  = α_eff_m·(p_m-p_m_n)+α_eff_f·(p_f-p_f_n) ✓
  //
  // K_eff is stored on the subRegion for post-correction of meanTotalStressIncrement
  // (whose effective-stress part uses K_m instead of K_eff).
  //
  // The fracture Biot wrapper is zeroed to avoid double-counting.
  void swapToCompositePressure( DomainPartition & domain )
  {
    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    for( size_t i = 0; i < matrixRegions.size(); ++i )
    {
      ElementRegionBase & matRegion = meshLevel1.getElemManager().getRegion( matrixRegions[i] );
      ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fractureRegions[i] );

      matRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
      {
        if( subRegIdx >= fracRegion.numSubRegions() ) return;
        if( !matSubReg.hasWrapper( viewKeyStruct::fracturePressureString() ) ) return;

        CellElementSubRegion & fracSubReg =
          dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

        arrayView1d< real64 > const p_m = matSubReg.getField< fields::flow::pressure >();
        arrayView1d< real64 > const p_m_n = matSubReg.getField< fields::flow::pressure_n >();
        arrayView1d< real64 const > const p_f_wrapper = matSubReg.getReference< array1d< real64 > >(
          viewKeyStruct::fracturePressureString() );
        arrayView1d< real64 const > const p_f_n = fracSubReg.getField< fields::flow::pressure_n >();
        arrayView1d< real64 > const alpha_f_wrapper = matSubReg.getReference< array1d< real64 > >(
          viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 > const K_eff_wrapper = matSubReg.getReference< array1d< real64 > >(
          viewKeyStruct::effectiveBulkModulusString() );

        // Matrix drained bulk modulus and Biot coefficient
        string const & matSolidName = matSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
        arrayView1d< real64 const > const K_m = matSolid.getBulkModulus();
        arrayView1d< real64 const > const alpha_m = matSolid.getBiotCoefficient();

        // Get writable views on matrix ElasticIsotropic to swap in K_eff, G_eff
        string const & matElasticName = matSolid.getReference< string >(
          constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
        constitutive::ElasticIsotropic & matElastic =
          matSolid.getParent().template getGroup< constitutive::ElasticIsotropic >( matElasticName );
        arrayView1d< real64 > const K_m_writable = matElastic.bulkModulus();
        arrayView1d< real64 > const G_m_writable = matElastic.shearModulus();

        // Fracture drained bulk modulus and Biot coefficient
        string const & fracSolidName = fracSubReg.getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
        arrayView1d< real64 const > const K_f = fracSolid.getBulkModulus();
        arrayView1d< real64 const > const alpha_f = fracSolid.getBiotCoefficient();

        // Get fracture shear modulus for G_eff
        string const & fracElasticName = fracSolid.getReference< string >(
          constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
        constitutive::ElasticIsotropic const & fracElastic =
          fracSolid.getParent().template getGroup< constitutive::ElasticIsotropic >( fracElasticName );
        arrayView1d< real64 const > const G_f = fracElastic.shearModulus();

        // Use user-specified volume fraction or compute from mesh element volumes
        bool const useMeshVolumes = (m_fractureVolumeFraction < 0.0);
        arrayView1d< real64 const > V_m, V_f;
        if( useMeshVolumes )
        {
          V_m = matSubReg.getElementVolume();
          V_f = fracSubReg.getElementVolume();
        }
        arrayView1d< localIndex const > const matrixToFracture =
          matSubReg.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        for( localIndex k = 0; k < matSubReg.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fracSubReg.size(),
                         "Invalid dual-continuum swapToCompositePressure connectivity for matrix subregion "
                         << matSubReg.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fracSubReg.getName()
                         << " size " << fracSubReg.size() );
          m_tempPm.push_back( p_m[k] );
          m_tempPm_n.push_back( p_m_n[k] );
          m_tempAlphaF.push_back( alpha_f_wrapper[k] );

          real64 const v_f = useMeshVolumes
            ? V_f[kf] / (V_m[k] + V_f[kf])
            : m_fractureVolumeFraction;
          real64 const v_m = 1.0 - v_f;

          // K_eff = (v_m/K_m + v_f/K_f)^{-1}  (Reuss average, Eq. A20)
          real64 const K_eff_inv = v_m / K_m[k] + v_f / K_f[kf];
          real64 const K_eff = 1.0 / K_eff_inv;

          // G_eff = (v_m/G_m + v_f/G_f)^{-1}  (Reuss average for shear)
          real64 const G_eff_inv = v_m / G_m_writable[k] + v_f / G_f[kf];
          real64 const G_eff = 1.0 / G_eff_inv;

          // α_eff_i = K_eff · v_i · α_i / K_i  (Eq. A21)
          real64 const alpha_eff_m = K_eff * v_m * alpha_m[k] / K_m[k];
          real64 const alpha_eff_f = K_eff * v_f * alpha_f[kf] / K_f[kf];

          // p_eq = (α_eff_m · p_m + α_eff_f · p_f) / α_m
          // p_eq_n = (α_eff_m · p_m_n + α_eff_f · p_f_n) / α_m
          real64 const p_eq = (alpha_eff_m * p_m[k] + alpha_eff_f * p_f_wrapper[k]) / alpha_m[k];
          real64 const p_eq_n = (alpha_eff_m * p_m_n[k] + alpha_eff_f * p_f_n[kf]) / alpha_m[k];

          m_tempCompositePressure.push_back( p_eq );
          m_tempCompositePressure_n.push_back( p_eq_n );
          K_eff_wrapper[k] = K_eff;

          // Save original matrix moduli and swap to effective
          m_tempKm.push_back( K_m_writable[k] );
          m_tempGm.push_back( G_m_writable[k] );
          K_m_writable[k] = K_eff;
          G_m_writable[k] = G_eff;

          p_m[k] = p_eq;
          p_m_n[k] = p_eq_n;
          // Zero out fracture Biot wrapper to avoid double-counting in the mechanics kernel.
          // The kernel has a legacy fracture Biot branch (-α_f·p_f·I) that was added before
          // the composite-pressure approach.  In sequential coupling the fracture contribution
          // is already accounted for through p_eq → α_eff_f·p_f, so the wrapper must be
          // disabled here.  It MAY still be needed by the monolithic path which does not call
          // swapToCompositePressure — verify before removing.
          alpha_f_wrapper[k] = 0.0;
        }
      } );
    }
  }

  void restoreCompositePressure( DomainPartition & domain )
  {
    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );

    size_t idx = 0;
    meshLevel1.getElemManager().forElementSubRegions< CellElementSubRegion >( matrixRegions,
      [&]( localIndex const, CellElementSubRegion & matrixSubRegion )
    {
      if( !matrixSubRegion.hasWrapper( viewKeyStruct::fracturePressureString() ) ) return;

      arrayView1d< real64 > const p_m = matrixSubRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 > const p_m_n = matrixSubRegion.getField< fields::flow::pressure_n >();
      arrayView1d< real64 > const alpha_f = matrixSubRegion.getReference< array1d< real64 > >(
        viewKeyStruct::fractureBiotCoefficientString() );

      // Restore matrix solid moduli
      string const & matSolidName = matrixSubRegion.template getReference< string >(
        Base::viewKeyStruct::porousMaterialNamesString() );
      constitutive::CoupledSolidBase & matSolid =
        this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matrixSubRegion, matSolidName );
      string const & matElasticName = matSolid.getReference< string >(
        constitutive::CoupledSolidBase::viewKeyStruct::solidModelNameString() );
      constitutive::ElasticIsotropic & matElastic =
        matSolid.getParent().template getGroup< constitutive::ElasticIsotropic >( matElasticName );
      arrayView1d< real64 > const K_view = matElastic.bulkModulus();
      arrayView1d< real64 > const G_view = matElastic.shearModulus();

      for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
      {
        p_m[k]   = m_tempPm[idx];
        p_m_n[k] = m_tempPm_n[idx];
        alpha_f[k] = m_tempAlphaF[idx];
        K_view[k] = m_tempKm[idx];
        G_view[k] = m_tempGm[idx];
        ++idx;
      }
    } );
    m_tempPm.clear();
    m_tempPm_n.clear();
    m_tempAlphaF.clear();
    m_tempKm.clear();
    m_tempGm.clear();
  }

  // Override assembleCouplingTerms — poromechanics coupling is handled by derived class.
  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    // DualContinuumFlowSolverBase::assembleCouplingTerms handles cross-flow within the flow
    // step. Here we only handle the flow-mechanics coupling.
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( domain );
    GEOS_UNUSED_VAR( dofManager );
    GEOS_UNUSED_VAR( localMatrix );
    GEOS_UNUSED_VAR( localRhs );
  }

private:

  // Map fracture pressure, Biot coeff, DOF# from mesh2 to mesh1
  void mapFractureDataToMatrix( DomainPartition & domain, DofManager const & dofManager )
  {
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    ElementRegionManager & elemManager1 = meshLevel1.getElemManager();
    ElementRegionManager & elemManager2 = meshLevel2.getElemManager();

    string const flowDofKey = isMultiphaseFlow
      ? dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() )
      : dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    using DualFlowSolver = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlowSolver & dualFlow = *this->flowSolver();
    string_array const & matrixRegionList = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList = dualFlow.template getReference< string_array >( "fractureRegionList" );
    if( matrixRegionList.empty() ) return;

    FieldIdentifiers fieldsToBeSync;
    if constexpr ( isMultiphaseFlow )
    {
      fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                          fields::flow::globalCompDensity::key() },
                                        fractureRegionList );
    }
    else
    {
      fieldsToBeSync.addElementFields( { fields::flow::pressure::key() },
                                        fractureRegionList );
    }
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, meshLevel2, domain.getNeighbors(), true );

    for( size_t iPair = 0; iPair < matrixRegionList.size(); ++iPair )
    {
      string const & matrixRegionName   = matrixRegionList[ iPair ];
      string const & fractureRegionName = fractureRegionList[ iPair ];

      ElementRegionBase & matrixRegion   = elemManager1.getRegion( matrixRegionName );
      ElementRegionBase & fractureRegion = elemManager2.getRegion( fractureRegionName );

      matrixRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matrixSubRegion )
      {
        if( subRegIdx >= fractureRegion.numSubRegions() ) return;
        CellElementSubRegion & fractureSubRegion =
          dynamic_cast< CellElementSubRegion & >( fractureRegion.getSubRegion( subRegIdx ) );

        arrayView1d< real64 > const fracturePressure =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fracturePressureString() );
        arrayView1d< real64 > const fractureBiotCoeff =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< globalIndex > const fractureDofNumber =
          matrixSubRegion.getReference< array1d< globalIndex > >( viewKeyStruct::fractureDofNumberString() );
        arrayView1d< real64 > const fractureMechanicsScale =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );

        arrayView1d< real64 const > const p_f = fractureSubRegion.getField< fields::flow::pressure >();
        arrayView1d< globalIndex const > const p_f_dof =
          fractureSubRegion.getReference< array1d< globalIndex > >( flowDofKey );

        string const & fracturePorousName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & fractureSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fractureSubRegion, fracturePorousName );
        arrayView1d< real64 const > const alpha_f = fractureSolid.getBiotCoefficient();
        arrayView1d< localIndex const > const matrixToFracture =
          matrixSubRegion.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );

        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fractureSubRegion.size(),
                         "Invalid dual-continuum mapFractureDataToMatrix connectivity for matrix subregion "
                         << matrixSubRegion.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fractureSubRegion.getName()
                         << " size " << fractureSubRegion.size() );
          fracturePressure[ k ]  = p_f[ kf ];
          if( !( m_useIntrinsicInput != 0 &&
                 this->getNonlinearSolverParameters().couplingType() ==
                 NonlinearSolverParameters::CouplingType::FullyImplicit ) )
          {
            fractureBiotCoeff[ k ] = alpha_f[ kf ];
          }
          fractureDofNumber[ k ] = p_f_dof[ kf ];
          fractureMechanicsScale[ k ] = 1.0;
        }
      } );
    }
  }

  // Assemble K_upf: fracture pressure contribution to displacement residual
  void assembleFractureMechanicsCoupling( DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    this->template forDiscretizationOnMeshTargets<>(
      domain.getMeshBodies(),
      [&]( string const & meshBodyName, MeshLevel & mesh, string_array const & regionNames )
    {
      if( meshBodyName != "mesh1" ) return;

      string const dispDofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
      NodeManager const & nodeManager = mesh.getNodeManager();
      if( !nodeManager.hasWrapper( dispDofKey ) ) return;

      globalIndex const rankOffset = dofManager.rankOffset();
      arrayView1d< globalIndex const > const dispDofNumber =
        nodeManager.getReference< globalIndex_array >( dispDofKey );

      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >(
        regionNames,
        [&]( localIndex const, CellElementSubRegion & subRegion )
      {
        if( !subRegion.hasWrapper( viewKeyStruct::fracturePressureString() ) ) return;

        arrayView1d< real64 const > const fracturePressure =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::fracturePressureString() );
        arrayView1d< real64 const > const fractureBiotCoeff =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 const > const fractureMechanicsScale =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );
        arrayView1d< globalIndex const > const fractureDofNumber =
          subRegion.getReference< array1d< globalIndex > >( viewKeyStruct::fractureDofNumberString() );

        finiteElement::FiniteElementBase & subRegionFE =
          subRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&]( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );
          constexpr localIndex numNodesPerElem = FE_TYPE::maxSupportPoints;
          constexpr integer numDofPerNode = 3;

          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList().toViewConst();
          arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
          typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
          finiteElement::FiniteElementBase::initialize< FE_TYPE >(
            nodeManager, mesh.getEdgeManager(), mesh.getFaceManager(), subRegion, meshData );

          RAJA::ReduceMax< parallelDeviceReduce, real64 > fractureMaxForce( 0.0 );
          RAJA::ReduceMax< parallelDeviceReduce, real64 > fractureMaxJac( 0.0 );
          RAJA::ReduceSum< parallelDeviceReduce, real64 > fractureAbsForce( 0.0 );
          RAJA::ReduceSum< parallelDeviceReduce, real64 > fractureAbsJac( 0.0 );

          forAll< parallelDevicePolicy<> >( subRegion.size(),
            [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            typename FE_TYPE::StackVariables feStack;
            finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
            localIndex const numSupportPoints = finiteElement.template numSupportPoints< FE_TYPE >( feStack );

            globalIndex localRowDofIndex[ numNodesPerElem * numDofPerNode ];
            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              localIndex const nodeIndex = elemsToNodes( k, a );
              for( integer dim = 0; dim < numDofPerNode; ++dim )
                localRowDofIndex[ numDofPerNode * a + dim ] = dofNumber[ nodeIndex ] + dim;
            }

            real64 const alpha_f = fractureMechanicsScale[ k ] * fractureBiotCoeff[ k ];
            real64 const p_f     = fracturePressure[ k ];
            if( fractureMechanicsScale[ k ] <= 0.0 || fractureDofNumber[ k ] < 0 )
            {
              return;
            }

            real64 localResidualMomentum[ numNodesPerElem * numDofPerNode ] = {};
            real64 dLocalResMomentum_dFracPressure[ numNodesPerElem * numDofPerNode ] = {};

            for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
            {
              real64 dNdX[ numNodesPerElem ][ 3 ], xLocal[ numNodesPerElem ][ 3 ];
              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                localIndex const nodeIndex = elemsToNodes( k, a );
                for( integer dim = 0; dim < numDofPerNode; ++dim )
                  xLocal[ a ][ dim ] = X[ nodeIndex ][ dim ];
              }
              real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, feStack, dNdX );

              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                for( integer dim = 0; dim < numDofPerNode; ++dim )
                {
                  real64 const gradPhi = dNdX[ a ][ dim ];
                  localResidualMomentum[ numDofPerNode * a + dim ] += alpha_f * p_f * gradPhi * detJxW;
                  dLocalResMomentum_dFracPressure[ numDofPerNode * a + dim ] += alpha_f * gradPhi * detJxW;
                }
              }
            }

            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              for( integer dim = 0; dim < numDofPerNode; ++dim )
              {
                localIndex const localRow = numDofPerNode * a + dim;
                localIndex const dof = LvArray::integerConversion< localIndex >( localRowDofIndex[ localRow ] - rankOffset );
                if( dof < 0 || dof >= localMatrix.numRows() ) continue;

                // Single-phase: the modified monolithic SinglePhasePoromechanics kernel already adds
                // the fracture-pressure term (sigma_tot includes -alpha_f*p_f) to the residual, so
                // K_upf contributes the Jacobian only. The STANDARD multiphase MultiphasePoromechanics
                // kernel knows nothing about the fracture, so for the compositional path K_upf must add
                // the residual (+alpha_f*p_f*gradN) here too, otherwise the Jacobian has no matching
                // residual (inconsistent Newton) and the fracture pressure never loads the mechanics.
                if constexpr ( isMultiphaseFlow )
                {
                  RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], localResidualMomentum[ localRow ] );
                  fractureMaxForce.max( fabs( localResidualMomentum[ localRow ] ) );
                  fractureAbsForce += fabs( localResidualMomentum[ localRow ] );
                }
                fractureMaxJac.max( fabs( dLocalResMomentum_dFracPressure[ localRow ] ) );
                fractureAbsJac += fabs( dLocalResMomentum_dFracPressure[ localRow ] );
                localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                  dof, &fractureDofNumber[ k ], &dLocalResMomentum_dFracPressure[ localRow ], 1 );
              }
            }
          } ); // forAll

          if constexpr ( isMultiphaseFlow )
          {
            this->solidMechanicsSolver()->getMaxForce() =
              LvArray::math::max( this->solidMechanicsSolver()->getMaxForce(), fractureMaxForce.get() );
          }
          if( m_logFimCouplingDiagnostics != 0 )
          {
            GEOS_LOG_RANK_0( GEOS_FMT( "{}: K_upf diagnostics on subRegion '{}': "
                                       "max|R_upf|={:.6e}, sum|R_upf|={:.6e}, "
                                       "max|dR_upf/dp_f|={:.6e}, sum|dR_upf/dp_f|={:.6e}",
                                       this->getName(), subRegion.getName(),
                                       fractureMaxForce.get(), fractureAbsForce.get(),
                                       fractureMaxJac.get(), fractureAbsJac.get() ) );
          }
        } ); // dispatch3D
      } ); // forElementSubRegions
    } ); // forDiscretizationOnMeshTargets
  }

  // Assemble K_pfu: displacement contribution to the fracture mass-balance Jacobian.
  // The fracture porosity follows the fixed-stress rule phi_f = phi_f,n + alpha_f*dEps_v + ...,
  // so the fracture fluid mass M_f = phi_f*rho_f*V depends on the matrix volumetric strain:
  //   dM_f/dU_{a,dim} = rho_f * alpha_f * (dN_a/dx_dim) * detJxW   (sum over quadrature points).
  // This is the off-diagonal that makes the monolithic Newton consistent; the corresponding
  // residual is already produced by updateFracturePorosityFixedStress + the fracture mass assembly.
  void assembleFractureToMechanicsCoupling( DomainPartition & domain,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
  {
    GEOS_UNUSED_VAR( localRhs );
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    this->template forDiscretizationOnMeshTargets<>(
      domain.getMeshBodies(),
      [&]( string const & meshBodyName, MeshLevel & mesh, string_array const & regionNames )
    {
      if( meshBodyName != "mesh1" ) return;

      string const dispDofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
      NodeManager const & nodeManager = mesh.getNodeManager();
      if( !nodeManager.hasWrapper( dispDofKey ) ) return;

      globalIndex const rankOffset = dofManager.rankOffset();

      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >(
        regionNames,
        [&]( localIndex const, CellElementSubRegion & subRegion )
      {
        if( !subRegion.hasWrapper( viewKeyStruct::fracturePressureString() ) ) return;

        arrayView1d< real64 const > const fractureBiotCoeff =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 const > const fractureMechanicsScale =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );
        arrayView1d< globalIndex const > const fractureDofNumber =
          subRegion.getReference< array1d< globalIndex > >( viewKeyStruct::fractureDofNumberString() );

        finiteElement::FiniteElementBase & subRegionFE =
          subRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&]( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );
          constexpr localIndex numNodesPerElem = FE_TYPE::maxSupportPoints;
          constexpr integer numDofPerNode = 3;

          // Reference fluid density for the fixed-stress mass linearization (matches the
          // rho_ref used in updateFracturePorosityFixedStress / fracture mass assembly).
          real64 constexpr rho_ref = 1000.0;

          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList().toViewConst();
          arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
          typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
          finiteElement::FiniteElementBase::initialize< FE_TYPE >(
            nodeManager, mesh.getEdgeManager(), mesh.getFaceManager(), subRegion, meshData );

          forAll< parallelDevicePolicy<> >( subRegion.size(),
            [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            typename FE_TYPE::StackVariables feStack;
            finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
            localIndex const numSupportPoints = finiteElement.template numSupportPoints< FE_TYPE >( feStack );

            globalIndex localColDofIndex[ numNodesPerElem * numDofPerNode ];
            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              localIndex const nodeIndex = elemsToNodes( k, a );
              for( integer dim = 0; dim < numDofPerNode; ++dim )
                localColDofIndex[ numDofPerNode * a + dim ] = dofNumber[ nodeIndex ] + dim;
            }

            real64 const alpha_f = fractureMechanicsScale[ k ] * fractureBiotCoeff[ k ];
            if( fractureMechanicsScale[ k ] <= 0.0 || fractureDofNumber[ k ] < 0 )
            {
              return;
            }

            real64 dFracMass_dU[ numNodesPerElem * numDofPerNode ] = {};

            for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
            {
              real64 dNdX[ numNodesPerElem ][ 3 ], xLocal[ numNodesPerElem ][ 3 ];
              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                localIndex const nodeIndex = elemsToNodes( k, a );
                for( integer dim = 0; dim < numDofPerNode; ++dim )
                  xLocal[ a ][ dim ] = X[ nodeIndex ][ dim ];
              }
              real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, feStack, dNdX );

              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                for( integer dim = 0; dim < numDofPerNode; ++dim )
                {
                  dFracMass_dU[ numDofPerNode * a + dim ] += rho_ref * alpha_f * dNdX[ a ][ dim ] * detJxW;
                }
              }
            }

            // Row = fracture mass DOF; columns = element displacement DOFs.
            localIndex const fracRow =
              LvArray::integerConversion< localIndex >( fractureDofNumber[ k ] - rankOffset );
            if( fracRow < 0 || fracRow >= localMatrix.numRows() ) return;

            localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
              fracRow, localColDofIndex, dFracMass_dU, numSupportPoints * numDofPerNode );
          } ); // forAll
        } ); // dispatch3D
      } ); // forElementSubRegions
    } ); // forDiscretizationOnMeshTargets
  }

  // Multiphase K_pfu: displacement -> fracture component-mass Jacobian.
  // Each fracture component mass is M_c = phi_f * compDens_c * V_elem, and the fixed-stress
  // fracture porosity carries the matrix volumetric strain
  //   phi_f = phi_f,n + alpha_f * dEps_v + (1/N) * dp_f,
  // so the missing off-diagonal (u not in the fracture flux stencil) is
  //   dM_c/dU_{a,dim} = compDens_c * alpha_f * (dN_a/dx_dim) * detJxW   (summed over q).
  // The useTotalMassEquation row transform (row 0 <- sum of component rows, others shifted ahead
  // by one; see CompositionalMultiphaseUtilities) is replicated so the entries land on exactly the
  // rows the fracture accumulation assembled. The volume-balance row (numComp) is porosity-free.
  void assembleFractureToMechanicsCouplingMultiphase( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs )
  {
    GEOS_UNUSED_VAR( localRhs );  // Jacobian-only; residual comes from the fracture accumulation
    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    string const dispDofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
    NodeManager const & nodeManager = meshLevel1.getNodeManager();
    if( !nodeManager.hasWrapper( dispDofKey ) ) return;
    string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );
    globalIndex const rankOffset = dofManager.rankOffset();
    arrayView1d< globalIndex const > const dispDofNumber =
      nodeManager.getReference< globalIndex_array >( dispDofKey );

    using DualFlowSolver = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlowSolver & dualFlow = *this->flowSolver();
    string_array const & matrixRegionList   = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList = dualFlow.template getReference< string_array >( "fractureRegionList" );
    if( matrixRegionList.empty() ) return;

    integer const numComp     = this->flowSolver()->secondarySolver()->numFluidComponents();
    integer const useTotalMass = this->flowSolver()->secondarySolver()->useTotalMassEquation();

    for( size_t iPair = 0; iPair < matrixRegionList.size(); ++iPair )
    {
      ElementRegionBase & matrixRegion   = meshLevel1.getElemManager().getRegion( matrixRegionList[ iPair ] );
      ElementRegionBase & fractureRegion = meshLevel2.getElemManager().getRegion( fractureRegionList[ iPair ] );

      matrixRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matrixSubRegion )
      {
        if( subRegIdx >= fractureRegion.numSubRegions() ) return;
        CellElementSubRegion & fractureSubRegion =
          dynamic_cast< CellElementSubRegion & >( fractureRegion.getSubRegion( subRegIdx ) );

        if( !matrixSubRegion.hasWrapper( viewKeyStruct::fractureBiotCoefficientString() ) ) return;
        arrayView1d< real64 const > const alpha_f =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 const > const fractureMechanicsScale =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );
        arrayView1d< globalIndex const > const fractureDofNumber =
          fractureSubRegion.getReference< array1d< globalIndex > >( flowDofKey );

        auto const compDens = fractureSubRegion.getField< fields::flow::globalCompDensity >();
        arrayView1d< localIndex const > const matrixToFracture =
          matrixSubRegion.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );
        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fractureSubRegion.size(),
                         "Invalid dual-continuum K_pfu connectivity for matrix subregion "
                         << matrixSubRegion.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fractureSubRegion.getName()
                         << " size " << fractureSubRegion.size() );
        }

        finiteElement::FiniteElementBase & subRegionFE =
          matrixSubRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&]( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );
          constexpr localIndex numNodesPerElem = FE_TYPE::maxSupportPoints;
          constexpr integer numDofPerNode = 3;
          constexpr integer maxNumComp = 16;

          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = matrixSubRegion.nodeList().toViewConst();
          typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
          finiteElement::FiniteElementBase::initialize< FE_TYPE >(
            nodeManager, meshLevel1.getEdgeManager(), meshLevel1.getFaceManager(), matrixSubRegion, meshData );

          integer const nc = numComp;
          integer const totalMass = useTotalMass;

          forAll< parallelDevicePolicy<> >( matrixSubRegion.size(),
            [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            typename FE_TYPE::StackVariables feStack;
            finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
            localIndex const kf = matrixToFracture[k];
            localIndex const numSupportPoints = finiteElement.template numSupportPoints< FE_TYPE >( feStack );
            integer const numDispDof = numSupportPoints * numDofPerNode;

            globalIndex localColDofIndex[ numNodesPerElem * numDofPerNode ];
            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              localIndex const ni = elemsToNodes( k, a );
              for( integer dim = 0; dim < numDofPerNode; ++dim )
                localColDofIndex[ numDofPerNode * a + dim ] = dispDofNumber[ ni ] + dim;
            }

            // geometric factor g[a*ND+dim] = alpha_f,strain * sum_q (dN_a/dx_dim) * detJxW
            real64 g[ numNodesPerElem * numDofPerNode ] = {};
            real64 const alpha = fractureMechanicsScale[ k ] * alpha_f[ k ];
            for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
            {
              real64 dNdX[ numNodesPerElem ][ 3 ], xLocal[ numNodesPerElem ][ 3 ];
              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                localIndex const ni = elemsToNodes( k, a );
                for( integer dim = 0; dim < 3; ++dim ) xLocal[ a ][ dim ] = X[ ni ][ dim ];
              }
              real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, feStack, dNdX );
              for( localIndex a = 0; a < numSupportPoints; ++a )
                for( integer dim = 0; dim < numDofPerNode; ++dim )
                  g[ numDofPerNode * a + dim ] += alpha * dNdX[ a ][ dim ] * detJxW;
            }

            // per-row coefficient zeta_c = compDens_c, with optional total-mass transform
            real64 coeff[ maxNumComp ] = {};
            if( totalMass )
            {
              real64 sum = 0.0;
              for( integer c = 0; c < nc; ++c ) sum += compDens[ kf ][ c ];
              coeff[ 0 ] = sum;
              for( integer i = 1; i < nc; ++i ) coeff[ i ] = compDens[ kf ][ i - 1 ];
            }
            else
            {
              for( integer c = 0; c < nc; ++c ) coeff[ c ] = compDens[ kf ][ c ];
            }

            // write each component-mass equation row (volume-balance row numComp is porosity-free)
            for( integer i = 0; i < nc; ++i )
            {
              globalIndex const rowDof = fractureDofNumber[ kf ] + i;
              localIndex const fracRow =
                LvArray::integerConversion< localIndex >( rowDof - rankOffset );
              if( fracRow < 0 || fracRow >= localMatrix.numRows() ) continue;

              real64 dRow[ numNodesPerElem * numDofPerNode ];
              for( integer j = 0; j < numDispDof; ++j ) dRow[ j ] = coeff[ i ] * g[ j ];

              arraySlice1d< globalIndex const > const cols = localMatrix.getColumns( fracRow );
              arraySlice1d< real64 > const vals = localMatrix.getEntries( fracRow );
              localIndex const nnz = localMatrix.numNonZeros( fracRow );
              for( integer j = 0; j < numDispDof; ++j )
              {
                for( localIndex a = 0; a < nnz; ++a )
                {
                  if( cols[a] == localColDofIndex[j] )
                  {
                    RAJA::atomicAdd< parallelDeviceAtomic >( &vals[a], dRow[j] );
                    break;
                  }
                }
              }
            }
          } ); // forAll
        } ); // dispatch3D
      } ); // forElementSubRegionsIndex
    }
  }

  // Update fracture porosity with matrix strain increment using fixed-stress formula.
  // Displacement→fracture-pressure coupling is implicit through porosity→mass in flow assembly.
  // Update fracture porosity using matrix volumetric strain increment.
  // The matrix displacement → fracture porosity coupling is implicit:
  // porosity update → flow solver mass → consistent K_pfpf and residual.
  void updateFracturePorosityFixedStress( DomainPartition & domain,
                                          real64 const & currentTime )
  {
    // Multiphase flow: delegate to the compositional variant, which updates the
    // fracture porosity with the same fixed-stress formula but lets the multiphase
    // accumulation kernel recompute the per-component masses (no single-phase
    // fields::flow::mass write). The single-phase body below is never instantiated
    // for the multiphase template (discarded if-constexpr branch is not odr-used).
    if constexpr ( isMultiphaseFlow )
    {
      updateFracturePorosityFixedStressMultiphase( domain, currentTime );
      return;
    }

    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    using DualFlowSolver = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlowSolver & dualFlow = *this->flowSolver();
    string_array const & matrixRegionList   = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList = dualFlow.template getReference< string_array >( "fractureRegionList" );
    if( matrixRegionList.empty() ) return;

    NodeManager const & nodeManager = meshLevel1.getNodeManager();
    if( !nodeManager.hasWrapper( fields::solidMechanics::incrementalDisplacement::key() ) ) return;
    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const incDisp =
      nodeManager.getField< fields::solidMechanics::incrementalDisplacement >();

    for( size_t iPair = 0; iPair < matrixRegionList.size(); ++iPair )
    {
      ElementRegionBase & matrixRegion   = meshLevel1.getElemManager().getRegion( matrixRegionList[ iPair ] );
      ElementRegionBase & fractureRegion = meshLevel2.getElemManager().getRegion( fractureRegionList[ iPair ] );

      matrixRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matrixSubRegion )
      {
        if( subRegIdx >= fractureRegion.numSubRegions() ) return;
        CellElementSubRegion & fractureSubRegion =
          dynamic_cast< CellElementSubRegion & >( fractureRegion.getSubRegion( subRegIdx ) );

        if( !matrixSubRegion.hasWrapper( viewKeyStruct::fractureBiotCoefficientString() ) ) return;
        arrayView1d< real64 const > const alpha_f =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 const > const fractureMechanicsScale =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );

        // Get fracture porosity on mesh2 for direct update
        string const & fracSolidName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          fractureSubRegion.getConstitutiveModel< constitutive::CoupledSolidBase >( fracSolidName );
        arrayView2d< real64 const > const phi_n_frac = fracSolid.getPorosity_n();
        arrayView2d< real64 > const phi_frac = fracSolid.getPorosity();
        arrayView2d< real64 > const dPhiFrac_dPres =
          fracSolid.getBasePorosityModel().getField< fields::porosity::dPorosity_dPressure >().reference();
        arrayView2d< real64 > const dPhiFrac_dTemp =
          fracSolid.getBasePorosityModel().getField< fields::porosity::dPorosity_dTemperature >().reference();

        // Get fracture mass and dMass for direct update (needed because
        // updateState runs before assembly, so mass must be updated here)
        arrayView1d< real64 > const mass2 = fractureSubRegion.getField< fields::flow::mass >();
        arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const dMass2 =
          fractureSubRegion.getField< fields::flow::dMass >();
        arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const rho_f =
          fractureSubRegion.getConstitutiveModel< constitutive::SingleFluidBase >(
            fractureSubRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() ) ).density();
        arrayView1d< real64 const > const volume = fractureSubRegion.getElementVolume();

        // Get fracture pressure for Biot porosity update
        arrayView1d< real64 const > const p_f = fractureSubRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const p_f_n = fractureSubRegion.getField< fields::flow::pressure_n >();

        // Get grain bulk modulus and reference porosity for the Biot formula
        constitutive::BiotPorosity const & fracBiotPorosity =
            dynamic_cast< constitutive::BiotPorosity const & >( fracSolid.getBasePorosityModel() );
        arrayView1d< real64 const > const grainBulkModulus = fracBiotPorosity.getGrainBulkModulus();
        arrayView1d< real64 const > const refPorosity = fracSolid.getReferencePorosity();

        // Get fluid density derivative for dMass/dP computation
        arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dFluidDensity =
            fractureSubRegion.getConstitutiveModel< constitutive::SingleFluidBase >(
                fractureSubRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() ) ).dDensity();
        arrayView1d< localIndex const > const matrixToFracture =
          matrixSubRegion.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );
        arrayView1d< integer const > const matrixGhostRank = matrixSubRegion.ghostRank();
        arrayView1d< integer const > const fractureGhostRank = fractureSubRegion.ghostRank();
        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fractureSubRegion.size(),
                         "Invalid dual-continuum single-phase fracture porosity connectivity for matrix subregion "
                         << matrixSubRegion.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fractureSubRegion.getName()
                         << " size " << fractureSubRegion.size() );
        }

        string const & matrixSolidName =
          matrixSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matrixSolid =
          matrixSubRegion.getConstitutiveModel< constitutive::CoupledSolidBase >( matrixSolidName );
        arrayView3d< real64 const, solid::STRESS_USD > const matrixEffectiveStress = matrixSolid.getEffectiveStress();
        localIndex const numMatrixStressQ = matrixEffectiveStress.size( 1 );

        // Get matrix FE space for strain computation
        finiteElement::FiniteElementBase & subRegionFE =
          matrixSubRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        constitutive::ConstitutivePassThru< constitutive::CoupledSolidBase >::execute(
          fracSolid, [&]( auto & castedFracSolid )
        {
          typename TYPEOFREF( castedFracSolid ) ::KernelWrapper fracSolidWrapper =
            castedFracSolid.createKernelUpdates();

          finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
          dispatch3D( subRegionFE, [&]( auto const finiteElement )
          {
            using FE_TYPE = decltype( finiteElement );
            constexpr localIndex NN = FE_TYPE::maxSupportPoints;
            constexpr integer ND = 3;

            arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
            arrayView2d< localIndex const, cells::NODE_MAP_USD > const e2n = matrixSubRegion.nodeList().toViewConst();
            typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
            finiteElement::FiniteElementBase::initialize< FE_TYPE >(
              nodeManager, meshLevel1.getEdgeManager(), meshLevel1.getFaceManager(), matrixSubRegion, meshData );

            forAll< parallelDevicePolicy<> >( matrixSubRegion.size(),
              [=] GEOS_HOST_DEVICE ( localIndex const k )
            {
              typename FE_TYPE::StackVariables feStack;
              finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
              localIndex const kf = matrixToFracture[k];
              if( matrixGhostRank[k] >= 0 || fractureGhostRank[kf] >= 0 )
              {
                return;
              }

              real64 uhat[ NN ][ ND ];
              for( localIndex a = 0; a < NN; ++a )
              { localIndex const ni = e2n( k, a );
                for( integer d = 0; d < ND; ++d ) uhat[ a ][ d ] = incDisp[ ni ][ d ]; }

              real64 volStrainIncSum = 0.0, elemVol = 0.0;
              for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
              {
                real64 dNdX[ NN ][ ND ], xLoc[ NN ][ ND ];
                for( localIndex a = 0; a < NN; ++a )
                { localIndex const ni = e2n( k, a );
                  for( integer d = 0; d < ND; ++d ) xLoc[ a ][ d ] = X[ ni ][ d ]; }
                real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLoc, feStack, dNdX );
                real64 strainInc[ 6 ] = {};
                FE_TYPE::symmetricGradient( dNdX, uhat, strainInc );
                volStrainIncSum += ( strainInc[0] + strainInc[1] + strainInc[2] ) * detJxW;
                elemVol += detJxW;
              }
              real64 const volStrainInc = ( elemVol > 0.0 ) ? volStrainIncSum / elemVol : 0.0;

              if( kf < phi_frac.size( 0 ) )
              {
                real64 const delta_p = p_f[ kf ] - p_f_n[ kf ];
                real64 const scaledAlpha = fractureMechanicsScale[ k ] * alpha_f[ k ];
                real64 const biotSkeletonModulusInverse =
                  ( scaledAlpha - refPorosity[ kf ] ) / grainBulkModulus[ kf ];
                real64 const phi_new =
                  phi_n_frac[ kf ][ 0 ] + scaledAlpha * volStrainInc + biotSkeletonModulusInverse * delta_p;
                phi_frac[ kf ][ 0 ] = phi_new;
                dPhiFrac_dPres[ kf ][ 0 ] = biotSkeletonModulusInverse;
                dPhiFrac_dTemp[ kf ][ 0 ] = 0.0;

                real64 const V = volume[ kf ];
                real64 const rho = rho_f[ kf ][ 0 ];
                mass2[ kf ] = phi_new * rho * V;

                real64 const drho_dp_val = dFluidDensity[ kf ][ 0 ][ 0 ];
                dMass2[ kf ][ 0 ] = (biotSkeletonModulusInverse * rho + phi_new * drho_dp_val) * V;

                real64 effectiveStress[6] = {};
                for( localIndex q = 0; q < numMatrixStressQ; ++q )
                {
                  for( integer component = 0; component < 6; ++component )
                  {
                    effectiveStress[component] += matrixEffectiveStress[k][q][component];
                  }
                }
                real64 const invNumMatrixStressQ = 1.0 / ( numMatrixStressQ > 0 ? numMatrixStressQ : 1 );
                for( integer component = 0; component < 6; ++component )
                {
                  effectiveStress[component] *= invNumMatrixStressQ;
                }
                fracSolidWrapper.updatePermeabilityFromEffectiveStress( kf, 0, effectiveStress, currentTime );
              }
            } ); // forAll
          } ); // dispatch3D
        } ); // ConstitutivePassThru
      } ); // forElementSubRegionsIndex
    }
  }

  // Multiphase fixed-stress fracture-porosity update.
  // Same fixed-stress rule as the single-phase variant
  //   phi_f = phi_f,n + alpha_f * dEps_v + (alpha_f - phi_ref)/K_s * dp_f
  // and its pressure derivative. The compositional
  // fracture flow solver's accumulation (Step 4) then recomputes the per-component
  // masses and Jacobian from this porosity — so there is no single-phase
  // fields::flow::mass / dMass manipulation here.
  void updateFracturePorosityFixedStressMultiphase( DomainPartition & domain,
                                                    real64 const & currentTime )
  {
    // The strain->fracture-porosity injection and its Jacobian counterpart (K_pfu) are
    // toggled together by m_enableFractureMechanicsCoupling so the monolithic Newton stays
    // consistent. When the explicit fracture<->mechanics coupling is off (default for the
    // first multiphase version), we also skip the Jacobian-less strain injection: the
    // fracture continuum is then coupled to the matrix only through cross-flow, and its
    // porosity follows pressure alone (with a consistent Jacobian from the fracture flow
    // accumulation). When K_pfu is implemented and enabled, the strain term is applied here.
    if( !m_enableFractureMechanicsCoupling || !m_enableFracturePorosityStrainCoupling ) return;

    if( domain.getMeshBodies().numSubGroups() < 2 ) return;

    MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
    MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
    MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );

    using DualFlowSolver = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlowSolver & dualFlow = *this->flowSolver();
    string_array const & matrixRegionList   = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList = dualFlow.template getReference< string_array >( "fractureRegionList" );
    if( matrixRegionList.empty() ) return;

    NodeManager const & nodeManager = meshLevel1.getNodeManager();
    if( !nodeManager.hasWrapper( fields::solidMechanics::incrementalDisplacement::key() ) ) return;
    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const incDisp =
      nodeManager.getField< fields::solidMechanics::incrementalDisplacement >();

    for( size_t iPair = 0; iPair < matrixRegionList.size(); ++iPair )
    {
      ElementRegionBase & matrixRegion   = meshLevel1.getElemManager().getRegion( matrixRegionList[ iPair ] );
      ElementRegionBase & fractureRegion = meshLevel2.getElemManager().getRegion( fractureRegionList[ iPair ] );

      matrixRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const subRegIdx, CellElementSubRegion & matrixSubRegion )
      {
        if( subRegIdx >= fractureRegion.numSubRegions() ) return;
        CellElementSubRegion & fractureSubRegion =
          dynamic_cast< CellElementSubRegion & >( fractureRegion.getSubRegion( subRegIdx ) );

        if( !matrixSubRegion.hasWrapper( viewKeyStruct::fractureBiotCoefficientString() ) ) return;
        arrayView1d< real64 const > const alpha_f =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() );
        arrayView1d< real64 const > const fractureMechanicsScale =
          matrixSubRegion.getReference< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() );

        string const & fracSolidName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          fractureSubRegion.getConstitutiveModel< constitutive::CoupledSolidBase >( fracSolidName );
        arrayView2d< real64 const > const phi_n_frac = fracSolid.getPorosity_n();
        arrayView2d< real64 > const phi_frac = fracSolid.getPorosity();
        arrayView2d< real64 > const dPhiFrac_dPres =
          fracSolid.getBasePorosityModel().getField< fields::porosity::dPorosity_dPressure >().reference();
        arrayView2d< real64 > const dPhiFrac_dTemp =
          fracSolid.getBasePorosityModel().getField< fields::porosity::dPorosity_dTemperature >().reference();

        arrayView1d< real64 const > const p_f = fractureSubRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const p_f_n = fractureSubRegion.getField< fields::flow::pressure_n >();

        constitutive::BiotPorosity const & fracBiotPorosity =
          dynamic_cast< constitutive::BiotPorosity const & >( fracSolid.getBasePorosityModel() );
        arrayView1d< real64 const > const grainBulkModulus = fracBiotPorosity.getGrainBulkModulus();
        arrayView1d< real64 const > const refPorosity = fracSolid.getReferencePorosity();
        arrayView1d< localIndex const > const matrixToFracture =
          matrixSubRegion.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );
        arrayView1d< integer const > const matrixGhostRank = matrixSubRegion.ghostRank();
        arrayView1d< integer const > const fractureGhostRank = fractureSubRegion.ghostRank();
        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          localIndex const kf = matrixToFracture[k];
          GEOS_ERROR_IF( kf < 0 || kf >= fractureSubRegion.size(),
                         "Invalid dual-continuum fracture porosity connectivity for matrix subregion "
                         << matrixSubRegion.getName() << ", local element " << k
                         << ": mapped fracture element " << kf
                         << " is outside fracture subregion " << fractureSubRegion.getName()
                         << " size " << fractureSubRegion.size() );
        }

        string const & matrixSolidName =
          matrixSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matrixSolid =
          matrixSubRegion.getConstitutiveModel< constitutive::CoupledSolidBase >( matrixSolidName );
        arrayView3d< real64 const, solid::STRESS_USD > const matrixEffectiveStress = matrixSolid.getEffectiveStress();
        localIndex const numMatrixStressQ = matrixEffectiveStress.size( 1 );

        finiteElement::FiniteElementBase & subRegionFE =
          matrixSubRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        constitutive::ConstitutivePassThru< constitutive::CoupledSolidBase >::execute(
          fracSolid, [&]( auto & castedFracSolid )
        {
          typename TYPEOFREF( castedFracSolid ) ::KernelWrapper fracSolidWrapper =
            castedFracSolid.createKernelUpdates();

          finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
          dispatch3D( subRegionFE, [&]( auto const finiteElement )
          {
            using FE_TYPE = decltype( finiteElement );
            constexpr localIndex NN = FE_TYPE::maxSupportPoints;
            constexpr integer ND = 3;

            arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
            arrayView2d< localIndex const, cells::NODE_MAP_USD > const e2n = matrixSubRegion.nodeList().toViewConst();
            typename FE_TYPE::template MeshData< CellElementSubRegion > meshData;
            finiteElement::FiniteElementBase::initialize< FE_TYPE >(
              nodeManager, meshLevel1.getEdgeManager(), meshLevel1.getFaceManager(), matrixSubRegion, meshData );

            forAll< parallelDevicePolicy<> >( matrixSubRegion.size(),
              [=] GEOS_HOST_DEVICE ( localIndex const k )
            {
              typename FE_TYPE::StackVariables feStack;
              finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
              localIndex const kf = matrixToFracture[k];
              if( matrixGhostRank[k] >= 0 || fractureGhostRank[kf] >= 0 )
              {
                return;
              }

              real64 uhat[ NN ][ ND ];
              for( localIndex a = 0; a < NN; ++a )
              { localIndex const ni = e2n( k, a );
                for( integer d = 0; d < ND; ++d ) uhat[ a ][ d ] = incDisp[ ni ][ d ]; }

              real64 volStrainIncSum = 0.0, elemVol = 0.0;
              for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
              {
                real64 dNdX[ NN ][ ND ], xLoc[ NN ][ ND ];
                for( localIndex a = 0; a < NN; ++a )
                { localIndex const ni = e2n( k, a );
                  for( integer d = 0; d < ND; ++d ) xLoc[ a ][ d ] = X[ ni ][ d ]; }
                real64 const detJxW = finiteElement.template getGradN< FE_TYPE >( k, q, xLoc, feStack, dNdX );
                real64 strainInc[ 6 ] = {};
                FE_TYPE::symmetricGradient( dNdX, uhat, strainInc );
                volStrainIncSum += ( strainInc[0] + strainInc[1] + strainInc[2] ) * detJxW;
                elemVol += detJxW;
              }
              real64 const volStrainInc = ( elemVol > 0.0 ) ? volStrainIncSum / elemVol : 0.0;

              if( kf < phi_frac.size( 0 ) )
              {
                real64 const delta_p = p_f[ kf ] - p_f_n[ kf ];
                real64 const scaledAlpha = fractureMechanicsScale[ k ] * alpha_f[ k ];
                real64 const biotSkeletonModulusInverse =
                  ( scaledAlpha - refPorosity[ kf ] ) / grainBulkModulus[ kf ];
                phi_frac[ kf ][ 0 ] = phi_n_frac[ kf ][ 0 ] + scaledAlpha * volStrainInc + biotSkeletonModulusInverse * delta_p;
                dPhiFrac_dPres[ kf ][ 0 ] = biotSkeletonModulusInverse;
                dPhiFrac_dTemp[ kf ][ 0 ] = 0.0;

                real64 effectiveStress[6] = {};
                for( localIndex q = 0; q < numMatrixStressQ; ++q )
                {
                  for( integer component = 0; component < 6; ++component )
                  {
                    effectiveStress[component] += matrixEffectiveStress[k][q][component];
                  }
                }
                real64 const invNumMatrixStressQ = 1.0 / ( numMatrixStressQ > 0 ? numMatrixStressQ : 1 );
                for( integer component = 0; component < 6; ++component )
                {
                  effectiveStress[component] *= invNumMatrixStressQ;
                }
                fracSolidWrapper.updatePermeabilityFromEffectiveStress( kf, 0, effectiveStress, currentTime );
              }
            } ); // forAll
          } ); // dispatch3D
        } ); // ConstitutivePassThru
      } ); // forElementSubRegionsIndex
    }
  }

  // Override registerDataOnMesh to ensure enableFixedStressPoromechanicsUpdate
  // is called BEFORE sub-solvers' registerDataOnMesh (which checks the flag
  // to determine whether to register pressure_k for Sequential coupling).
  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override
  {
    // Enable fixed-stress update BEFORE Base::registerDataOnMesh so that
    // sub-solvers register pressure_k during their own registerDataOnMesh.
    if( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
    {
      this->solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
      this->flowSolver()->enableFixedStressPoromechanicsUpdate();
    }
    Base::registerDataOnMesh( meshBodies );

    // Register mapped fracture fields on mesh1 element subregions
    if( meshBodies.hasGroup( "mesh1" ) )
    {
      MeshBody & mesh1 = meshBodies.getGroup< MeshBody >( "mesh1" );
      MeshLevel & meshLevel = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
      ElementRegionManager & elemManager = meshLevel.getElemManager();

      for( localIndex iRegion = 0; iRegion < elemManager.numRegions(); ++iRegion )
      {
        ElementRegionBase & region = elemManager.getRegion( iRegion );
        region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
        {
          subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::fracturePressureString() )
            .setApplyDefaultValue( 0.0 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() )
            .setApplyDefaultValue( 0.0 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< globalIndex > >( viewKeyStruct::fractureDofNumberString() )
            .setApplyDefaultValue( -1 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::fractureMechanicsScaleString() )
            .setApplyDefaultValue( 0.0 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::effectiveBulkModulusString() )
            .setApplyDefaultValue( -1 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
        } );
      }
    }
  }

public:

  virtual void implicitStepComplete( real64 const & time, real64 const & dt,
                                     DomainPartition & domain ) override
  {
    // For sequential coupling, fracture porosity was already updated in
    // mapSolution(SolidMechanics). Calling updateFracturePorosityFixedStress
    // here would manually overwrite mass_f (absorbing strain into mass_n).
    // For monolithic coupling, the manual update in assembleSystem is needed.
    bool const isSequential =
      this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential;
    if( !isSequential )
    {
      updateFracturePorosityFixedStress( domain, time );
    }
    Base::implicitStepComplete( time, dt, domain );
  }

  // Override initializePostInitialConditionsPreSubGroups to setup dual continuum
  virtual void initializePostInitialConditionsPreSubGroups() override
  {
    Base::initializePostInitialConditionsPreSubGroups();

    // Setup dual continuum cross flow
    DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
    MeshBody & matrix = domain.getMeshBody("mesh1");
    MeshBody & fracture = domain.getMeshBody("mesh2");
    MeshLevel & primaryMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & secondaryMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );
    this->flowSolver()->initializePostInitialConditionsPreSubGroups();

    // Set volumeFraction on matrix (v_m) and fracture (v_f) subregions.
    // The flow solver's mass accumulation (updateMass) does φ * ρ * V_elem,
    // but the physical pore volume per dual-porosity REV is φ * v * V_elem.
    // Scaling the element volume by volumeFraction corrects this.

    GEOS_UNUSED_VAR( primaryMesh, secondaryMesh );
  }

  // Override initializePostInitialConditionsPostSubGroups for gravity setup
  virtual void initializePostInitialConditionsPostSubGroups() override
  {
    Base::initializePostInitialConditionsPostSubGroups();

    DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
    MeshBody & matrix = domain.getMeshBody("mesh1");
    MeshBody & fracture = domain.getMeshBody("mesh2");
    MeshLevel & primaryMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & secondaryMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );
    this->flowSolver()->initializePostInitialConditionsPostSubGroups();
    GEOS_UNUSED_VAR( primaryMesh, secondaryMesh );

    // If requested, homogenize intrinsic constitutive parameters into the effective medium now
    // (FullyImplicit only; no-op by default). Done after constitutive post-init so the moduli /
    // Biot arrays are populated from their defaults before we overwrite them.
    computeEffectiveFromIntrinsic( domain );

    // If requested, set the matrix initial effective stress to the dual-continuum total Biot stress
    // so the model starts in mechanical equilibrium (no manual matrixSolid_stress needed). Done after
    // computeEffectiveFromIntrinsic so the (possibly homogenized) Biot coefficients are final.
    autoInitializeEffectiveStress( domain );

    if constexpr ( isMultiphaseFlow )
    {
      if( this->getNonlinearSolverParameters().couplingType() ==
          NonlinearSolverParameters::CouplingType::FullyImplicit )
      {
        enforceFractureCompositionalVolumeClosure( domain, true );
      }
    }
  }

private:

  // Report the inventory at the sequential split points.  The stages are recorded
  // because mechanics changes porosity between the flow and mechanics sub-solves,
  // and updateFluidState() can then rewrite compAmount from phi*V*rho.  Comparing
  // afterFlow, afterPorosity, and afterFluidState therefore identifies the exact
  // stage that changes the inventory, while matrix/fracture and total sums distinguish
  // continuum redistribution from a non-conservative change in global mass.  The values
  // are reconstructed from the same compAmount field used by the accumulation kernel,
  // so this diagnostic observes the split coupling without perturbing the solve.
  void logSequentialMassInventory( DomainPartition & domain, char const * stage ) const
  {
    if constexpr ( !isMultiphaseFlow )
    {
      GEOS_UNUSED_VAR( domain, stage );
      return;
    }

    MeshBody & matrix = domain.getMeshBody( "mesh1" );
    MeshBody & fracture = domain.getMeshBody( "mesh2" );
    MeshLevel & matrixMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & fractureMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );

    using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlow & dualFlow = *this->flowSolver();
    string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

    real64 co2M = 0.0;
    real64 co2MN = 0.0;
    real64 co2DerivedM = 0.0;
    real64 poreM = 0.0;
    real64 co2F = 0.0;
    real64 co2FN = 0.0;
    real64 co2DerivedF = 0.0;
    real64 poreF = 0.0;

    CompositionalMultiphaseFormulationType const formulationType =
      dualFlow.primarySolver()->template getReference< CompositionalMultiphaseFormulationType >(
        CompositionalMultiphaseBase::viewKeyStruct::formulationTypeString() );

    auto accumulate = [&]( MeshLevel & mesh,
                           string_array const & regions,
                           real64 & co2,
                           real64 & co2N,
                           real64 & co2Derived,
                           real64 & poreVolume )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >(
        regions,
        [&]( localIndex const, CellElementSubRegion & subRegion )
      {
        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
        arrayView1d< real64 const > const volume = subRegion.getElementVolume();
        string const & solidName = subRegion.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & porousMaterial =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();
        arrayView2d< real64 const, compflow::USD_COMP > const compAmount =
          subRegion.template getField< fields::flow::compAmount >();
        arrayView2d< real64 const, compflow::USD_COMP > const compAmountN =
          subRegion.template getField< fields::flow::compAmount_n >();

        arrayView2d< real64 const, compflow::USD_COMP > compDens;
        arrayView2d< real64 const, compflow::USD_COMP > compFrac;
        arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > totalDens;
        if( formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
        {
          compFrac = subRegion.template getField< fields::flow::globalCompFraction >();
          string const & fluidName = subRegion.template getReference< string >(
            CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
          constitutive::MultiFluidBase const & fluid =
            this->template getConstitutiveModel< constitutive::MultiFluidBase >( subRegion, fluidName );
          totalDens = fluid.totalDensity();
        }
        else
        {
          compDens = subRegion.template getField< fields::flow::globalCompDensity >();
        }

        for( localIndex k = 0; k < subRegion.size(); ++k )
        {
          if( ghostRank[k] < 0 )
          {
            co2 += compAmount[k][0];
            co2N += compAmountN[k][0];
            real64 const currentComponentDensity =
              formulationType == CompositionalMultiphaseFormulationType::OverallComposition
              ? compFrac[k][0] * totalDens[k][0]
              : compDens[k][0];
            co2Derived += porosity[k][0] * volume[k] * currentComponentDensity;
            poreVolume += porosity[k][0] * volume[k];
          }
        }
      } );
    };

    accumulate( matrixMesh, matrixRegions, co2M, co2MN, co2DerivedM, poreM );
    accumulate( fractureMesh, fractureRegions, co2F, co2FN, co2DerivedF, poreF );

    co2M = MpiWrapper::sum( co2M );
    co2MN = MpiWrapper::sum( co2MN );
    co2DerivedM = MpiWrapper::sum( co2DerivedM );
    poreM = MpiWrapper::sum( poreM );
    co2F = MpiWrapper::sum( co2F );
    co2FN = MpiWrapper::sum( co2FN );
    co2DerivedF = MpiWrapper::sum( co2DerivedF );
    poreF = MpiWrapper::sum( poreF );

    real64 const total = co2M + co2F;
    real64 const totalN = co2MN + co2FN;
    real64 const derivedTotal = co2DerivedM + co2DerivedF;
    real64 const avgDensM = poreM > 0.0 ? co2M / poreM : 0.0;
    real64 const avgDensF = poreF > 0.0 ? co2F / poreF : 0.0;

    GEOS_LOG_RANK_0( GEOS_FMT( "{}: SEQ-MASS stage={} co2M={:.12e} co2MN={:.12e} dM={:.12e} "
                               "derivedM={:.12e} storedMinusDerivedM={:.12e} poreM={:.12e} cbarM={:.12e} "
                               "co2F={:.12e} co2FN={:.12e} dF={:.12e} derivedF={:.12e} "
                               "storedMinusDerivedF={:.12e} poreF={:.12e} cbarF={:.12e} "
                               "total={:.12e} totalN={:.12e} dTotal={:.12e} derivedTotal={:.12e} "
                               "storedMinusDerivedTotal={:.12e}",
                               this->getName(), stage,
                               co2M, co2MN, co2M - co2MN, co2DerivedM, co2M - co2DerivedM, poreM, avgDensM,
                               co2F, co2FN, co2F - co2FN, co2DerivedF, co2F - co2DerivedF, poreF, avgDensF,
                               total, totalN, total - totalN, derivedTotal, total - derivedTotal ) );
  }

  // Accessor for dual continuum flow solver
  DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER > * dualContinuumFlowSolver() const
  {
    return std::get< toUnderlying( SolverType::DualContinuumFlow ) >( m_solvers );
  }

  // Override updateBulkDensity if needed for dual continuum
  virtual void updateBulkDensity( ElementSubRegionBase & subRegion ) override
  {
    // Default implementation does nothing. Concrete derived solvers should
    // provide an appropriate implementation when needed.
    GEOS_UNUSED_VAR( subRegion );
  }

  // Additional methods for dual continuum specific functionality
  real64 getFracSpacingLz() const
  {
    return this->flowSolver()->getFracSpacingLz();
  }

  void updateGravityPressure( MeshLevel & meshMatrix, MeshLevel & meshFracture, real64 const & gravityCoefficient )
  {
    this->flowSolver()->updateGravityPressure( meshMatrix, meshFracture, gravityCoefficient );
  }

public:

  struct viewKeyStruct : Base::viewKeyStruct
  {
    static constexpr char const * fracturePressureString() { return "fracturePressure"; }
    static constexpr char const * fractureBiotCoefficientString() { return "fractureBiotCoefficient"; }
    static constexpr char const * fractureDofNumberString() { return "fractureDofNumber"; }
    static constexpr char const * fractureMechanicsScaleString() { return "fractureMechanicsScale"; }
    static constexpr char const * fractureVolumeFractionString() { return "fractureVolumeFraction"; }
    static constexpr char const * fimNewtonRelaxationString() { return "fimNewtonRelaxation"; }
    static constexpr char const * sequentialPressureRelaxationString() { return "sequentialPressureRelaxation"; }
    static constexpr char const * enableFractureMechanicsCouplingString() { return "enableFractureMechanicsCoupling"; }
    static constexpr char const * enableFracturePorosityStrainCouplingString()
    { return "enableFracturePorosityStrainCoupling"; }
    static constexpr char const * enableFimCrossStorageString() { return "enableFimCrossStorage"; }
    static constexpr char const * logFimCouplingDiagnosticsString() { return "logFimCouplingDiagnostics"; }
    static constexpr char const * logSequentialMassDiagnosticsString() { return "logSequentialMassDiagnostics"; }
    static constexpr char const * useIntrinsicInputString() { return "useIntrinsicInput"; }
    static constexpr char const * autoInitializeStressString() { return "autoInitializeStress"; }
    static constexpr char const * effectiveBulkModulusString() { return "effectiveBulkModulus"; }
    static constexpr char const * volumeFractionString() { return "volumeFraction"; }
  };

  /// User-specified fracture volume fraction; (<0) → computed from mesh volumes
  real64 m_fractureVolumeFraction;

  // TODO: replace with registered fields for GPU compatibility
  std::vector< real64 > m_tempPm;
  std::vector< real64 > m_tempPm_n;
  std::vector< real64 > m_tempAlphaF;
  std::vector< real64 > m_tempKm;                   // original K_m before swap
  std::vector< real64 > m_tempGm;                   // original G_m before swap
  std::vector< real64 > m_tempCompositePressure;    // p_eq during mechanics step
  std::vector< real64 > m_tempVolStrainIncr;        // shared dEps_v fed to both continua
  std::vector< real64 > m_tempCompositePressure_n;  // p_eq_n during mechanics step

  // Toggle for the explicit fracture<->mechanics Jacobian coupling (K_upf / K_pfu) in the FIM
  // path. Used to isolate the effect of those terms during convergence debugging.
  integer m_enableFractureMechanicsCoupling = 1;

  integer m_enableFracturePorosityStrainCoupling = 1;

  // Fixed Newton step relaxation applied in FullyImplicit (FIM) mode to damp the period-2
  // (lambda~=-1) pressure-block oscillation. 0.5 is near-optimal (midpoint = exact solution);
  // 1.0 disables relaxation.
  real64 m_fimNewtonRelaxation = 0.5;

  real64 m_sequentialPressureRelaxation = 1.0;

  // Toggle the multi-porosity (Mehrabian S_ij) cross-storage correction in the FIM path.
  integer m_enableFimCrossStorage = 1;

  integer m_logFimCouplingDiagnostics = 0;

  integer m_logSequentialMassDiagnostics = 0;

  // Input mode for dual-continuum mechanics/storage parameters (default 1 = intrinsic input).
  // Porosity and permeability are always intrinsic input and are scaled by the continuum volume
  // fractions in the dual-continuum flow solver.
  //   1 (intrinsic): the deck sets the INTRINSIC matrix/fracture moduli/Biot (K_m,G_m,alpha_m /
  //                  K_f,G_f,alpha_f) on the constitutive models. At initialization the solver
  //                  homogenizes them for FullyImplicit and Sequential. Removes the need to
  //                  pre-compute effective parameters by hand.
  //   0 (effective): the deck sets the HOMOGENIZED effective moduli/Biot and direct effective
  //                  storage coefficients. Nothing is auto-computed and Sequential uses the
  //                  physical matrix/fracture pressures directly.
  integer m_useIntrinsicInput = 1;

  // Auto-initialize matrix effective stress to the dual-continuum total Biot stress so the initial
  // total stress is in equilibrium (Rsolid~0 at t=0), mirroring single-porosity poromechanics.
  // Default 0 (backward compatible); the compositional dual solver defaults it to 1.
  integer m_autoInitializeStress = 0;

};

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMPOROMECHANICSSOLVER_HPP_
