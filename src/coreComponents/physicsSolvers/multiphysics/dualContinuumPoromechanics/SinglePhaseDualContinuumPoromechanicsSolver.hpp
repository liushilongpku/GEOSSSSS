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
 * @file SinglePhaseDualContinuumPoromechanicsSolver.hpp
 *
 * @brief Single-phase dual continuum poromechanics solver.
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"

namespace geos
{

class SinglePhaseDualContinuumPoromechanicsSolver : public DualContinuumPoromechanicsSolverBase< SinglePhaseBase, SinglePhaseBase, SolidMechanicsLagrangianFEM >
{
public:

  using Base = DualContinuumPoromechanicsSolverBase< SinglePhaseBase, SinglePhaseBase, SolidMechanicsLagrangianFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "singlephasedualcontinuumporomechanics"; }

  /**
   * @brief main constructor for SinglePhaseDualContinuumPoromechanicsSolver Objects
   * @param name the name of this instantiation of SinglePhaseDualContinuumPoromechanicsSolver in the repository
   * @param parent the parent group of this instantiation of SinglePhaseDualContinuumPoromechanicsSolver
   */
  SinglePhaseDualContinuumPoromechanicsSolver( const string & name,
                                               dataRepository::Group * const parent )
    : Base( name, parent )
  {}

  /**
   * @brief default destructor
   */
  ~SinglePhaseDualContinuumPoromechanicsSolver() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name
   */
  static string catalogName() { return "SinglePhaseDualContinuumPoromechanics"; }

  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override
  {
    Base::initializePreSubGroups();

    this->getWrapper< integer >( Base::viewKeyStruct::isThermalString() ).setApplyDefaultValue( 0 );
  }

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    Base::setupDofs( domain, dofManager );//Base dpdk-mech
    // 在基类中已经设置了耦合关系，在此无需重复。
    // Set up coupling between flow and mechanics
    // this->setupCoupling( domain, dofManager );
  }

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override
  {
    // Call base setup coupling for dual continuum flow
    Base::setupCoupling( domain, dofManager );//which is couple solver

    // TODO@LSL Additional coupling setup if needed between flow and mechanics
    // For now, rely on base class
  }

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    // Assemble coupling terms for dual continuum flow
    Base::assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );

    // Assemble coupling terms between flow and solid mechanics
    assemblePoromechanicsCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  virtual void updateBulkDensity( ElementSubRegionBase & subRegion ) override
  {
    // Use the same bulk density update as SinglePhasePoromechanics: launch kernel
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    constitutive::SingleFluidBase const & fluid = this->template getConstitutiveModel< constitutive::SingleFluidBase >( subRegion, fluidName );

    string const & solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
    constitutive::CoupledSolidBase const & solid = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );

    poromechanicsKernels::
      SinglePhaseBulkDensityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( fluid,
                                                 solid,
                                                 subRegion );
  }

protected:

  void assemblePoromechanicsCouplingTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
  {
    // Placeholder: no additional coupling is enforced here.
    // This function exists to provide an entry point for future coupling implementations.
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( domain );
    GEOS_UNUSED_VAR( dofManager );
    GEOS_UNUSED_VAR( localMatrix );
    GEOS_UNUSED_VAR( localRhs );
  }

  struct viewKeyStruct : Base::viewKeyStruct
  {
    // Add any specific view keys if needed
  };

};

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_