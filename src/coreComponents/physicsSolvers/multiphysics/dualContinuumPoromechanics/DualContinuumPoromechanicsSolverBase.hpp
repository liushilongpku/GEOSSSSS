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
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsKernelsDispatchTypeList.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/logger/Logger.hpp"

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
    // Additional initialization if needed, e.g., for dual continuum specific parameters
  }

  // Override postInputInitialization to include checks for dual continuum
  virtual void postInputInitialization() override
  {
    Base::postInputInitialization();

    // Additional checks for dual continuum, e.g., ensure flow solvers are compatible
    GEOS_THROW_IF( this->flowSolver()->primarySolver()->isThermal() != this->flowSolver()->secondarySolver()->isThermal(),
                   GEOS_FMT( "{} {}: Primary and secondary flow solvers must have the same thermal setting",
                             this->getCatalogName(), this->getName() ),
                   InputError );

    // Enable fixed-stress poromechanics update BEFORE registerDataOnMesh is called.
    // registerDataOnMesh runs on sub-solvers before the coupled solver, so the flags
    // must be set during postInputInitialization for pressure_k to be registered.
    if( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
    {
      this->solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
      this->flowSolver()->enableFixedStressPoromechanicsUpdate();
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
    
    // 后续有需要再设置裂缝的流固耦合
    // dofManager.addCoupling( fields::solidMechanics::totalDisplacement::key(),
    //                         SECONDARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString(),
    //                         DofManager::Connector::Elem );
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
    GEOS_UNUSED_VAR( time_n );

    // ---- Step 1: Monolithic matrix kernel (u + p_m) ----
    // Assembles K_uu, K_upm, K_pmu, K_pmpm in one quadrature-point loop
    // on mesh1 (matrix) regions only.
    string const flowDofKey =
      dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    this->template forDiscretizationOnMeshTargets<>(
      domain.getMeshBodies(),
      [&]( string const & meshBodyName,
           MeshLevel & mesh,
           string_array const & regionNames )
    {
      if( meshBodyName != "mesh1" ) return;   // matrix mesh only

      this->template assemblyLaunch<
        PoromechanicsKernelsDispatchTypeList,
        poromechanicsKernels::SinglePhasePoromechanicsKernelFactory >(
          mesh, dofManager, regionNames,
          viewKeyStruct::porousMaterialNamesString(),
          localMatrix, localRhs, dt,
          flowDofKey,
          this->m_performStressInitialization,
          FlowSolverBase::viewKeyStruct::fluidNamesString() );
    } );

    // ---- Step 2: Matrix face-based flux terms ----
    this->flowSolver()->primarySolver()->assembleFluxTerms(
      dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 3: Fracture flow assembly (K_pfpf accumulation + face fluxes) ----
    this->flowSolver()->secondarySolver()->assembleSystem(
      time_n, dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4: Cross-flow coupling (K_pmpf, K_pfpm) ----
    this->flowSolver()->assembleCouplingTerms(
      time_n, dt, domain, dofManager, localMatrix, localRhs );
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

  struct viewKeyStruct : Base::viewKeyStruct
  {
    // Add any additional view keys if needed
  };

};

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMPOROMECHANICSSOLVER_HPP_