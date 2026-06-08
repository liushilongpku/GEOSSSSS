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
 * @file CompositionalMultiphaseDualContinuumPoromechanicsSolver.hpp
 *
 * @brief Compositional multiphase dual continuum poromechanics solver (FIM).
 *
 * Extends the dual-porosity (DPDP) hydro-mechanical coupling to compositional
 * multiphase flow. The matrix continuum (mesh1) is assembled monolithically with
 * the standard GEOS multiphase poromechanics kernel (u + p + components), the
 * fracture continuum (mesh2) flow is assembled by the compositional dual flow
 * solver, and the matrix<->fracture cross-flow / cross-storage is handled by
 * DualContinuumCrossFlow. The fracture<->mechanics fixed-stress coupling is
 * applied through the fracture porosity update (see
 * DualContinuumPoromechanicsSolverBase::updateFracturePorosityFixedStressMultiphase).
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"

namespace geos
{

class CompositionalMultiphaseDualContinuumPoromechanicsSolver
  : public DualContinuumPoromechanicsSolverBase< CompositionalMultiphaseFVM, CompositionalMultiphaseFVM, SolidMechanicsLagrangianFEM >
{
public:

  using Base = DualContinuumPoromechanicsSolverBase< CompositionalMultiphaseFVM, CompositionalMultiphaseFVM, SolidMechanicsLagrangianFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "compositionalmultiphasedualcontinuumporomechanics"; }

  /**
   * @brief main constructor
   * @param name the name of this instantiation in the repository
   * @param parent the parent group of this instantiation
   */
  CompositionalMultiphaseDualContinuumPoromechanicsSolver( const string & name,
                                                           dataRepository::Group * const parent )
    : Base( name, parent )
  {}

  ~CompositionalMultiphaseDualContinuumPoromechanicsSolver() override = default;

  /**
   * @brief name of the solver in the object catalog
   */
  static string catalogName() { return "CompositionalMultiphaseDualContinuumPoromechanics"; }

  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override
  {
    Base::initializePreSubGroups();

    // Isothermal by default.
    this->getWrapper< integer >( Base::viewKeyStruct::isThermalString() ).setApplyDefaultValue( 0 );

    // Fracture<->mechanics coupling is fully consistent for the compositional path:
    // K_upf (alpha_f*p_f in the displacement residual + its p_f Jacobian) and the
    // per-component K_pfu (d(phi_f*compDens_c*V)/dU, with the useTotalMassEquation row
    // transform) are both assembled, together with the strain term in the fracture
    // porosity (updateFracturePorosityFixedStressMultiphase). This makes the monolithic
    // dual-porosity Newton consistent (true FIM) so it converges without step damping.
    this->getWrapper< integer >( Base::viewKeyStruct::enableFractureMechanicsCouplingString() ).setApplyDefaultValue( 1 );

    // Automatically initialize the matrix effective stress to the dual-continuum total Biot stress
    // sigma' = alpha_m*p_m + alpha_f*p_f so the model starts in mechanical equilibrium (Rsolid~0 at
    // t=0), the same way single-porosity poromechanics does -- no manual matrixSolid_stress needed.
    this->getWrapper< integer >( Base::viewKeyStruct::autoInitializeStressString() ).setApplyDefaultValue( 1 );
  }

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    Base::setupDofs( domain, dofManager );
  }

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override
  {
    Base::setupCoupling( domain, dofManager );
  }

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    // Poromechanics coupling (K_upm, K_pmu) is assembled by the monolithic
    // MultiphasePoromechanics kernel in DualContinuumPoromechanicsSolverBase::assembleSystem.
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( domain );
    GEOS_UNUSED_VAR( dofManager );
    GEOS_UNUSED_VAR( localMatrix );
    GEOS_UNUSED_VAR( localRhs );
  }

  virtual void updateBulkDensity( ElementSubRegionBase & subRegion ) override
  {
    // Same bulk-density update as MultiphasePoromechanics: launch the multiphase
    // bulk-density kernel using the fluid + coupled-solid models on the subregion.
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    constitutive::MultiFluidBase const & fluid = this->template getConstitutiveModel< constitutive::MultiFluidBase >( subRegion, fluidName );

    string const & solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
    constitutive::CoupledSolidBase const & solid = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );

    poromechanicsKernels::
      MultiphaseBulkDensityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( this->flowSolver()->primarySolver()->numFluidComponents(),
                                                 fluid,
                                                 solid,
                                                 subRegion );
  }

protected:

  struct viewKeyStruct : Base::viewKeyStruct
  {};

};

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_
