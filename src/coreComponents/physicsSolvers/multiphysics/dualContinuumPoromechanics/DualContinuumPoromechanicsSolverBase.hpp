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
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/logger/Logger.hpp"
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

    // 裂缝的流固耦合关系 — K_upf disabled for numerical stability
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

    // ---- Step 0: Map fracture data (p_f, α_f, DOF#) from mesh2 to mesh1 ----
    mapFractureDataToMatrix( domain, dofManager );

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

    // ---- Step 1.5: Update fracture porosity with matrix strain (fixed-stress) ----
    updateFracturePorosityFixedStress( domain );

    // ---- Step 2: Fracture-mechanics coupling (K_upf) ----
    // K_upf disabled: fracture pressure→displacement coupling is numerically unstable.
    // Displacement→fracture coupling (K_pfu) is implicit through porosity/mass update.
    // assembleFractureMechanicsCoupling( domain, dofManager, localMatrix, localRhs );

    // ---- Step 3: Matrix face-based flux terms ----
    this->flowSolver()->primarySolver()->assembleFluxTerms(
      dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4: Fracture flow assembly (K_pfpf accumulation + face fluxes) ----
    this->flowSolver()->secondarySolver()->assembleSystem(
      time_n, dt, domain, dofManager, localMatrix, localRhs );

    // ---- Step 4: Cross-flow coupling (K_pmpf, K_pfpm) ----
    this->flowSolver()->assembleCouplingTerms(
      time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  // Override mapSolutionBetweenSolvers to handle dual-mesh sequential coupling.
  // 1. After flow: add -α_f·p_f·I to stored effective stress (so mechanics kernel
  //    computes correct total stress without needing K_upf).
  // 2. After mechanics: restore stress, average mean stress on mesh1 only,
  //    update porosity/permeability on mesh1 only.
  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override
  {
    if( solverType == static_cast< integer >( SolverType::DualContinuumFlow ) )
    {
      copyFracturePressureToMesh1( domain );
      {
        // DEBUG: show p_m and p_f before composite swap
        MeshBody & m1 = domain.getMeshBody( "mesh1" );
        MeshBody & m2 = domain.getMeshBody( "mesh2" );
        MeshLevel & ml1 = m1.getMeshLevels().getGroup< MeshLevel >( 0 );
        MeshLevel & ml2 = m2.getMeshLevels().getGroup< MeshLevel >( 0 );
        auto & df = *this->flowSolver();
        auto const & matRegs = df.template getReference< string_array >( "matrixRegionList" );
        auto const & fracRegs = df.template getReference< string_array >( "fractureRegionList" );
        if( matRegs.size() > 0 && ml1.getElemManager().getRegion(matRegs[0]).numSubRegions() > 0 )
        {
          auto & msr = dynamic_cast<CellElementSubRegion&>(ml1.getElemManager().getRegion(matRegs[0]).getSubRegion(0));
          auto & fsr = dynamic_cast<CellElementSubRegion&>(ml2.getElemManager().getRegion(fracRegs[0]).getSubRegion(0));
          printf("[DBG] Flow→Mech: elem0 p_m=%.2f p_f=%.2f\n",
                 msr.getField<fields::flow::pressure>()[0],
                 fsr.getField<fields::flow::pressure>()[0]);
        }
      }
      swapToCompositePressure( domain );
      Base::updateBulkDensity( domain );
    }
    else if( solverType == static_cast< integer >( SolverType::SolidMechanics )
             && !this->m_performStressInitialization )
    {
      // Restore p_m BEFORE the matrix pass so updatePorosityAndPermeability
      // uses the original matrix pressure, not the composite p_eq.
      // Otherwise the matrix uses p_eq (which includes p_f) while the
      // fracture uses p_f, creating a ~0.5·biotSkelInv·(p_m-p_f) difference
      // in porosity that amplifies into divergence over subsequent steps.
      restoreCompositePressure( domain );

      // Average mean total stress on mesh1 only (mesh2 has no FE)
      MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
      MeshLevel & meshLevel = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
      using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
      DualFlow & dualFlow = *this->flowSolver();
      string_array const & matrixRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
      meshLevel.getElemManager().forElementSubRegions< CellElementSubRegion >( matrixRegions,
        [&]( localIndex const, auto & subRegion )
      {
        string const & solidName = subRegion.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & solid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );

        arrayView2d< real64 const > const meanTotalStressIncrement_k = solid.getMeanTotalStressIncrement_k();
        arrayView1d< real64 > const averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();

        finiteElement::FiniteElementBase & subRegionFE =
          subRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&]( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );
          AverageOverQuadraturePoints1DKernelFactory::
            createAndLaunch< CellElementSubRegion, FE_TYPE, parallelDevicePolicy<> >(
              meshLevel.getNodeManager(), meshLevel.getEdgeManager(), meshLevel.getFaceManager(),
              subRegion, finiteElement, meanTotalStressIncrement_k, averageMeanTotalStressIncrement_k );
        } );

        this->flowSolver()->updatePorosityAndPermeability( subRegion );
        this->updateBulkDensity( subRegion );
      } );

      // Copy matrix avgStressIncr to fracture so the fracture flow
      // solver's fixed-stress porosity update includes the mechanics strain.
      {
        MeshBody & mesh1 = domain.getMeshBody( "mesh1" );
        MeshBody & mesh2 = domain.getMeshBody( "mesh2" );
        MeshLevel & meshLevel1 = mesh1.getMeshLevels().getGroup< MeshLevel >( 0 );
        MeshLevel & meshLevel2 = mesh2.getMeshLevels().getGroup< MeshLevel >( 0 );
        using DualFlow = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
        DualFlow & dualFlow = *this->flowSolver();
        string_array const & matRegions = dualFlow.template getReference< string_array >( "matrixRegionList" );
        string_array const & fracRegions = dualFlow.template getReference< string_array >( "fractureRegionList" );

        for( size_t iPair = 0; iPair < matRegions.size(); ++iPair )
        {
          ElementRegionBase & matRegion = meshLevel1.getElemManager().getRegion( matRegions[ iPair ] );
          ElementRegionBase & fracRegion = meshLevel2.getElemManager().getRegion( fracRegions[ iPair ] );

          matRegion.forElementSubRegionsIndex< CellElementSubRegion >(
            [&]( localIndex const subRegIdx, CellElementSubRegion & matSubReg )
          {
            if( subRegIdx >= fracRegion.numSubRegions() ) return;
            CellElementSubRegion & fracSubReg =
              dynamic_cast< CellElementSubRegion & >( fracRegion.getSubRegion( subRegIdx ) );

            string const & matSolidName = matSubReg.template getReference< string >(
              Base::viewKeyStruct::porousMaterialNamesString() );
            constitutive::CoupledSolidBase & matSolid =
              this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
            arrayView1d< real64 const > const matAvgStressIncr = matSolid.getAverageMeanTotalStressIncrement_k();

            string const & fracSolidName =
              fracSubReg.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
            constitutive::CoupledSolidBase & fracSolid =
              this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
            arrayView1d< real64 > const fracAvgStressIncr = fracSolid.getAverageMeanTotalStressIncrement_k();

            for( localIndex k = 0; k < matSubReg.size(); ++k )
            {
              fracAvgStressIncr[ k ] = matAvgStressIncr[ k ];
            }

            // DEBUG before
            {
              arrayView2d<real64 const> phi_f = fracSolid.getPorosity();
              printf("[DBG] frac elem0 BEFORE updatePorosity: phi=%.6e avgStress=%.6e p=%.2f\n",
                     phi_f[0][0], fracAvgStressIncr[0],
                     fracSubReg.getField<fields::flow::pressure>()[0]);
            }
            this->flowSolver()->secondarySolver()->updatePorosityAndPermeability( fracSubReg );
            {
              arrayView2d<real64 const> phi_f = fracSolid.getPorosity();
              printf("[DBG] frac elem0 AFTER  updatePorosity: phi=%.6e\n", phi_f[0][0]);
            }
          } );
        }
      }
    }
  }

private:

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
        string const & fracturePorousName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & fractureSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fractureSubRegion, fracturePorousName );
        arrayView1d< real64 const > const alpha_f_mesh2 = fractureSolid.getBiotCoefficient();

        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          p_f_wrapper[ k ] = p_f_mesh2[ k ];
          alpha_f_wrapper[ k ] = alpha_f_mesh2[ k ];
        }
      } );
    }
  }

  // Replace matrix pressure on mesh1 with equivalent single-porosity pressure
  // derived from Mehrabian & Abousleiman (2014), Appendix A, Eqs. A20-A27:
  //
  //   K_eff   = (v_m/K_m + v_f/K_f)^{-1}                     (Eq. A20, Reuss average)
  //   α_eff_i = K_eff · v_i · α_i / K_i                      (Eq. A21)
  //   p_eq    = α_eff_m · p_m + α_eff_f · p_f
  //
  // p_eq replaces p_m so the mechanics kernel (which multiplies by α_m) computes:
  //   σ = σ'_drained - α_m · p_eq = σ'_drained - (α_m·α_eff_m·p_m + α_m·α_eff_f·p_f)
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
        arrayView1d< real64 const > const p_f_wrapper = matSubReg.getReference< array1d< real64 > >(
          viewKeyStruct::fracturePressureString() );
        arrayView1d< real64 > const alpha_f_wrapper = matSubReg.getReference< array1d< real64 > >(
          viewKeyStruct::fractureBiotCoefficientString() );

        // Element volumes for volume fractions
        arrayView1d< real64 const > const V_m = matSubReg.getElementVolume();
        arrayView1d< real64 const > const V_f = fracSubReg.getElementVolume();

        // Matrix drained bulk modulus and Biot coefficient
        string const & matSolidName = matSubReg.template getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & matSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSubReg, matSolidName );
        arrayView1d< real64 const > const K_m = matSolid.getBulkModulus();
        arrayView1d< real64 const > const alpha_m = matSolid.getBiotCoefficient();

        // Fracture drained bulk modulus and Biot coefficient
        string const & fracSolidName = fracSubReg.getReference< string >(
          Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSubReg, fracSolidName );
        arrayView1d< real64 const > const K_f = fracSolid.getBulkModulus();
        arrayView1d< real64 const > const alpha_f = fracSolid.getBiotCoefficient();

        for( localIndex k = 0; k < matSubReg.size(); ++k )
        {
          m_tempPm.push_back( p_m[k] );
          m_tempAlphaF.push_back( alpha_f_wrapper[k] );

          real64 const V_total = V_m[k] + V_f[k];
          real64 const v_m = V_m[k] / V_total;
          real64 const v_f = V_f[k] / V_total;

          // K_eff = (v_m/K_m + v_f/K_f)^{-1}  (Reuss average, Eq. A20)
          real64 const K_eff_inv = v_m / K_m[k] + v_f / K_f[k];
          real64 const K_eff = 1.0 / K_eff_inv;

          // α_eff_i = K_eff · v_i · α_i / K_i  (Eq. A21)
          real64 const alpha_eff_m = K_eff * v_m * alpha_m[k] / K_m[k];
          real64 const alpha_eff_f = K_eff * v_f * alpha_f[k] / K_f[k];

          // p_eq = α_eff_m · p_m + α_eff_f · p_f  (no division by α_m)
          p_m[k] = alpha_eff_m * p_m[k] + alpha_eff_f * p_f_wrapper[k];
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
      arrayView1d< real64 > const alpha_f = matrixSubRegion.getReference< array1d< real64 > >(
        viewKeyStruct::fractureBiotCoefficientString() );

      for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
      {
        p_m[k]   = m_tempPm[idx];
        alpha_f[k] = m_tempAlphaF[idx];
        ++idx;
      }
    } );
    m_tempPm.clear();
    m_tempAlphaF.clear();
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

    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    using DualFlowSolver = DualContinuumFlowSolverBase< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
    DualFlowSolver & dualFlow = *this->flowSolver();
    string_array const & matrixRegionList = dualFlow.template getReference< string_array >( "matrixRegionList" );
    string_array const & fractureRegionList = dualFlow.template getReference< string_array >( "fractureRegionList" );
    if( matrixRegionList.empty() ) return;

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

        arrayView1d< real64 const > const p_f = fractureSubRegion.getField< fields::flow::pressure >();
        arrayView1d< globalIndex const > const p_f_dof =
          fractureSubRegion.getReference< array1d< globalIndex > >( flowDofKey );

        string const & fracturePorousName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase const & fractureSolid =
          this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fractureSubRegion, fracturePorousName );
        arrayView1d< real64 const > const alpha_f = fractureSolid.getBiotCoefficient();

        for( localIndex k = 0; k < matrixSubRegion.size(); ++k )
        {
          fracturePressure[ k ]  = p_f[ k ];
          fractureBiotCoeff[ k ] = alpha_f[ k ];
          fractureDofNumber[ k ] = p_f_dof[ k ];
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

            real64 localResidualMomentum[ numNodesPerElem * numDofPerNode ] = {};
            real64 dLocalResMomentum_dFracPressure[ numNodesPerElem * numDofPerNode ] = {};
            real64 const alpha_f = fractureBiotCoeff[ k ];
            real64 const p_f     = fracturePressure[ k ];

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

                RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], localResidualMomentum[ localRow ] );
                localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                  dof, &fractureDofNumber[ k ], &dLocalResMomentum_dFracPressure[ localRow ], 1 );
              }
            }
          } ); // forAll
        } ); // dispatch3D
      } ); // forElementSubRegions
    } ); // forDiscretizationOnMeshTargets
  }

  // Update fracture porosity with matrix strain increment using fixed-stress formula.
  // Displacement→fracture-pressure coupling is implicit through porosity→mass in flow assembly.
  // Update fracture porosity using matrix volumetric strain increment.
  // The matrix displacement → fracture porosity coupling is implicit:
  // porosity update → flow solver mass → consistent K_pfpf and residual.
  void updateFracturePorosityFixedStress( DomainPartition & domain )
  {
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

        // Get fracture porosity on mesh2 for direct update
        string const & fracSolidName =
          fractureSubRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & fracSolid =
          fractureSubRegion.getConstitutiveModel< constitutive::CoupledSolidBase >( fracSolidName );
        arrayView2d< real64 const > const phi_n_frac = fracSolid.getPorosity_n();
        arrayView2d< real64 > const phi_frac = fracSolid.getPorosity();

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

        // Get matrix FE space for strain computation
        finiteElement::FiniteElementBase & subRegionFE =
          matrixSubRegion.template getReference< finiteElement::FiniteElementBase >(
            this->solidMechanicsSolver()->getDiscretizationName() );

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

            if( k < phi_frac.size( 0 ) )
            {
              real64 const delta_p = p_f[ k ] - p_f_n[ k ];
              real64 const biotSkeletonModulusInverse = (alpha_f[ k ] - refPorosity[ k ]) / grainBulkModulus[ k ];
              real64 const phi_new = phi_n_frac[ k ][ 0 ] + alpha_f[ k ] * volStrainInc + biotSkeletonModulusInverse * delta_p;
              phi_frac[ k ][ 0 ] = phi_new;

              real64 const V = volume[ k ];
              real64 const rho = rho_f[ k ][ 0 ];
              mass2[ k ] = phi_new * rho * V;

              real64 const drho_dp_val = dFluidDensity[ k ][ 0 ][ 0 ];
              dMass2[ k ][ 0 ] = (biotSkeletonModulusInverse * rho + phi_new * drho_dp_val) * V;
            }
          } ); // forAll
        } ); // dispatch3D
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
            .setApplyDefaultValue( -1 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::fractureBiotCoefficientString() )
            .setApplyDefaultValue( -1 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );

          subRegion.registerWrapper< array1d< globalIndex > >( viewKeyStruct::fractureDofNumberString() )
            .setApplyDefaultValue( -1 )
            .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
            .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
        } );
      }
    }
  }

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
      updateFracturePorosityFixedStress( domain );
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

public:

  struct viewKeyStruct : Base::viewKeyStruct
  {
    static constexpr char const * fracturePressureString() { return "fracturePressure"; }
    static constexpr char const * fractureBiotCoefficientString() { return "fractureBiotCoefficient"; }
    static constexpr char const * fractureDofNumberString() { return "fractureDofNumber"; }
  };

  // TODO: replace with registered fields for GPU compatibility
  std::vector< real64 > m_tempPm;
  std::vector< real64 > m_tempAlphaF;

};

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMPOROMECHANICSSOLVER_HPP_