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

  // Override assembleCouplingTerms to include dual continuum coupling
  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    // DualContinuumFlowSolverBase implements its own coupling in assembleSystem (through its own assembleCouplingTerms),
    // so we do not directly call its protected assembleCouplingTerms here.

    // Base class (Poromechanics) coupling is still applied.
    Base::assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
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