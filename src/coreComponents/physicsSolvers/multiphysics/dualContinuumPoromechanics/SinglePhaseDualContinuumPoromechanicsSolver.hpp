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
#include "finiteElement/FiniteElementDispatch.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"
#include "mesh/CellElementSubRegion.hpp"

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
    // Poromechanics coupling (K_upm, K_pmu) is now handled by the monolithic
    // SinglePhasePoromechanics kernel launched in
    // DualContinuumPoromechanicsSolverBase::assembleSystem, Step 1.
    // The block-by-block coupling kernel (PoromechanicsCouplingKernel) is
    // retained below under #if 0 for reference / regression comparison.
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( domain );
    GEOS_UNUSED_VAR( dofManager );
    GEOS_UNUSED_VAR( localMatrix );
    GEOS_UNUSED_VAR( localRhs );
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

#if 0
  // Block-by-block coupling kernel — retained for reference / regression.
  // Replaced by monolithic SinglePhasePoromechanics kernel in
  // DualContinuumPoromechanicsSolverBase::assembleSystem.
  void assemblePoromechanicsCouplingTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs );
#endif

  struct viewKeyStruct : Base::viewKeyStruct
  {
    // Add any specific view keys if needed
  };

};

#if 0
// ============================================================================
// BLOCK-BY-BLOCK COUPLING KERNEL — retained for reference / regression.
// This kernel was used in the original block-by-block assembly approach.
// It has been replaced by the monolithic SinglePhasePoromechanics kernel
// launched in DualContinuumPoromechanicsSolverBase::assembleSystem.
//
// To re-activate:
//   1. Change #if 0 → #if 1 below.
//   2. In DualContinuumPoromechanicsSolverBase::assembleSystem, revert to
//      the block-by-block version (see git history).
// ============================================================================

namespace
{

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class PoromechanicsCouplingKernel
{
public:

  static constexpr integer numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  static constexpr integer numDofPerNode = 3;

  PoromechanicsCouplingKernel( NodeManager const & nodeManager,
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               SUBREGION_TYPE const & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               arrayView1d< globalIndex const > const dispDofNumber,
                               arrayView1d< globalIndex const > const flowDofNumber,
                               arrayView1d< real64 const > const biotCoefficient,
                               arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const fluidDensity,
                               real64 const dt,
                               globalIndex const dofRankOffset,
                               CRSMatrixView< real64, globalIndex const > const localMatrix )
    : m_finiteElementSpace( finiteElementSpace ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_X( nodeManager.referencePosition() ),
    m_dispDofNumber( dispDofNumber ),
    m_flowDofNumber( flowDofNumber ),
    m_biotCoefficient( biotCoefficient ),
    m_fluidDensity( fluidDensity ),
    m_dt( dt ),
    m_dofRankOffset( dofRankOffset ),
    m_matrix( localMatrix )
  {
    finiteElement::FiniteElementBase::
      initialize< FE_TYPE >( nodeManager,
                             edgeManager,
                             faceManager,
                             elementSubRegion,
                             m_meshData );
  }

  struct StackVariables
  {
    real64 xLocal[ FE_TYPE::maxSupportPoints ][ 3 ];
    typename FE_TYPE::StackVariables feStack;

    globalIndex localRowDofIndex[ FE_TYPE::maxSupportPoints * numDofPerNode ];
    globalIndex localPressureDofIndex;

    real64 dLocalResidualMomentum_dPressure[ FE_TYPE::maxSupportPoints * numDofPerNode ]{};
    real64 dLocalResidualMass_dDisplacement[ FE_TYPE::maxSupportPoints * numDofPerNode ]{};
  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );

    for( localIndex a = 0; a < numSupportPoints; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );

      for( localIndex dim = 0; dim < numDofPerNode; ++dim )
      {
        stack.xLocal[ a ][ dim ] = m_X[ nodeIndex ][ dim ];
        stack.localRowDofIndex[ numDofPerNode * a + dim ] =
          m_dispDofNumber[ nodeIndex ] + dim;
      }
    }

    stack.localPressureDofIndex = m_flowDofNumber[ k ];

    for( localIndex i = 0; i < FE_TYPE::maxSupportPoints * numDofPerNode; ++i )
    {
      stack.dLocalResidualMomentum_dPressure[ i ] = 0.0;
      stack.dLocalResidualMass_dDisplacement[ i ] = 0.0;
    }
  }

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 dNdX[ FE_TYPE::maxSupportPoints ][ 3 ]{};
    real64 const detJxW =
      m_finiteElementSpace.template getGradN< FE_TYPE >( k, q,
                                                         stack.xLocal,
                                                         stack.feStack,
                                                         dNdX );

    real64 const alpha = m_biotCoefficient[ k ];
    real64 const rho   = m_fluidDensity[ k ][ 0 ];
    real64 const invDt = 1.0 / m_dt;

    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );

    for( localIndex a = 0; a < numSupportPoints; ++a )
    {
      for( localIndex dim = 0; dim < numDofPerNode; ++dim )
      {
        real64 const gradPhi = dNdX[ a ][ dim ];

        // K_upm: dR_momentum / d(pressure) = -alpha * gradPhi * detJw
        stack.dLocalResidualMomentum_dPressure[ numDofPerNode * a + dim ]
          -= alpha * gradPhi * detJxW;

        // K_pmu: dR_mass / d(displacement) = rho * alpha / dt * gradPhi * detJw
        stack.dLocalResidualMass_dDisplacement[ numDofPerNode * a + dim ]
          += rho * alpha * invDt * gradPhi * detJxW;
      }
    }
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const k,
                 StackVariables & stack ) const
  {
    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    integer const numDispDof = numSupportPoints * numDofPerNode;

    // --- Assemble K_upm: displacement-equation rows ← pressure column ---
    for( localIndex localNode = 0; localNode < numSupportPoints; ++localNode )
    {
      for( localIndex dim = 0; dim < numDofPerNode; ++dim )
      {
        localIndex const dof =
          LvArray::integerConversion< localIndex >(
            stack.localRowDofIndex[ numDofPerNode * localNode + dim ]
            - m_dofRankOffset );

        if( dof < 0 || dof >= m_matrix.numRows() )
        {
          continue;
        }

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
          dof,
          &stack.localPressureDofIndex,
          &stack.dLocalResidualMomentum_dPressure[ numDofPerNode * localNode + dim ],
          1 );
      }
    }

    // --- Assemble K_pmu: pressure-equation row ← displacement columns ---
    localIndex const dof =
      LvArray::integerConversion< localIndex >(
        stack.localPressureDofIndex - m_dofRankOffset );

    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >(
        dof,
        stack.localRowDofIndex,
        stack.dLocalResidualMass_dDisplacement,
        numDispDof );
    }
  }

  template< typename POLICY,
            typename KERNEL_TYPE >
  static void
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q = 0; q < KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      kernelComponent.complete( k, stack );
    } );
  }

protected:

  FE_TYPE const & m_finiteElementSpace;

  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type >
  const m_elemsToNodes;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  arrayView1d< globalIndex const > const m_dispDofNumber;

  arrayView1d< globalIndex const > const m_flowDofNumber;

  arrayView1d< real64 const > const m_biotCoefficient;

  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID >
  const m_fluidDensity;

  real64 const m_dt;

  globalIndex const m_dofRankOffset;

  CRSMatrixView< real64, globalIndex const > const m_matrix;

  typename FE_TYPE::template MeshData< SUBREGION_TYPE > m_meshData;
};

} // anonymous namespace

// -----------------------------------------------------------------------------
// assemblePoromechanicsCouplingTerms (block-by-block)
// -----------------------------------------------------------------------------

void
SinglePhaseDualContinuumPoromechanicsSolver::
assemblePoromechanicsCouplingTerms( real64 const time_n,
                                    real64 const dt,
                                    DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( localRhs );  // coupling terms only — no RHS assembly here

  // Use const_cast because the GEOS data repository requires non-const access
  // even for read-only operations (getConstitutiveModel, getReference, etc.).
  // This matches the pattern used by all other assembly functions in the codebase.
  DomainPartition & nonConstDomain = const_cast< DomainPartition & >( domain );

  // Iterate over poromechanics target regions.
  // NOTE: displacement DOFs live only on mesh1 (the matrix mesh). The fracture
  // mesh (mesh2) has pressure DOFs but NO displacement DOFs. Therefore we
  // specifically avoid assembling coupling terms on mesh2 regions.
  // The K_upm block couples displacement (mesh1) ↔ matrix pressure (mesh1).
  // The K_upf block (displacement ↔ fracture pressure) is deferred (alpha_f ≈ 0).
  // The K_pmpf / K_pfpm blocks (cross-flow) are handled by DualContinuumFVM.

  this->template forDiscretizationOnMeshTargets<>(
    nonConstDomain.getMeshBodies(),
    [&]( string const & meshBodyName,
         MeshLevel & mesh,
         string_array const & regionNames )
  {
    // Skip mesh bodies that have no displacement DOFs
    string const dispDofKey =
      dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
    NodeManager const & nodeManager = mesh.getNodeManager();
    if( !nodeManager.hasWrapper( dispDofKey ) )
    {
      return;  // mesh2 (fracture) has no displacement DOFs — skip
    }
    GEOS_UNUSED_VAR( meshBodyName );

    arrayView1d< globalIndex const > const dispDofNumber =
      nodeManager.getReference< globalIndex_array >( dispDofKey );

    string const flowDofKey =
      dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >(
      regionNames,
      [&]( localIndex const,
           CellElementSubRegion & subRegion )
    {
      // ---- Pressure DOF numbers ----
      arrayView1d< globalIndex const > const flowDofNumber =
        subRegion.template getReference< array1d< globalIndex > >( flowDofKey );

      // ---- Biot coefficient (matrix rock) ----
      string const & solidName =
        subRegion.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
      constitutive::CoupledSolidBase const & solidModel =
        this->template getConstitutiveModel< constitutive::CoupledSolidBase >(
          subRegion, solidName );
      arrayView1d< real64 const > const biotCoefficient =
        solidModel.getBiotCoefficient();

      // ---- Fluid density (matrix fluid) ----
      string const & fluidName =
        subRegion.template getReference< string >(
          FlowSolverBase::viewKeyStruct::fluidNamesString() );
      constitutive::SingleFluidBase const & fluid =
        this->template getConstitutiveModel< constitutive::SingleFluidBase >(
          subRegion, fluidName );
      arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID >
      const fluidDensity = fluid.density();

      // ---- FE space (same discretization as solid mechanics) ----
      finiteElement::FiniteElementBase & subRegionFE =
        subRegion.template getReference< finiteElement::FiniteElementBase >(
          this->solidMechanicsSolver()->getDiscretizationName() );

      // Dispatch on the element type and launch the coupling kernel
      finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
      dispatch3D( subRegionFE,
                  [&]( auto const finiteElement )
      {
        using FE_TYPE = decltype( finiteElement );

        PoromechanicsCouplingKernel< CellElementSubRegion, FE_TYPE >
        kernel( nodeManager,
                mesh.getEdgeManager(),
                mesh.getFaceManager(),
                subRegion,
                finiteElement,
                dispDofNumber,
                flowDofNumber,
                biotCoefficient,
                fluidDensity,
                dt,
                dofManager.rankOffset(),
                localMatrix );

        PoromechanicsCouplingKernel< CellElementSubRegion, FE_TYPE >::
        template kernelLaunch< parallelDevicePolicy<> >(
          subRegion.size(), kernel );
      } );
    } ); // forElementSubRegions
  } ); // forDiscretizationOnMeshTargets
}

#endif // #if 0 — block-by-block coupling kernel

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEDUALCONTINUUMPOROMECHANICSSOLVER_HPP_