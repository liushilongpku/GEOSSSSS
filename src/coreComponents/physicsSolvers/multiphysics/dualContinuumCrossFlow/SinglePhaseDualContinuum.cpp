/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumFVM.cpp
 */

#include "SinglePhaseDualContinuum.hpp"

// #include "common/TimingMacros.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/singlePhase/CrossFlowComputeKernel.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/singlePhase/ThermalCrossFlowComputeKernel.hpp"
#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
namespace geos
{
  using namespace dataRepository;
  using namespace fields;

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::
    assembleCouplingTerms(real64 const time_n,
                          real64 const dt,
                          DomainPartition const &domain,
                          DofManager const &dofManager,
                          CRSMatrixView<real64, globalIndex const> const &localMatrix,
                          arrayView1d<real64> const &localRhs)
  {
    GEOS_LOG("befor the assemble");
    //this->printCRSMatrix(this->m_localMatrix);

    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR(time_n);
    GEOS_UNUSED_VAR(dt);
    // 注册传导率场
    // 根据 aperture 计算传导率
    // 当前使用常数
    NumericalMethodsManager const &numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const &fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const &fluxApprox = fvManager.getFluxApproximation(Base::primarySolver()->getDiscretizationName());
    //FluxApproximationBase const &secondaryFluxApprox = fvManager.getFluxApproximation(Base::secondarySolver()->getDiscretizationName());
    //FluxApproximationBase const *const fluxApproxArray[] = {&primaryFluxApprox,&secondaryFluxApprox};

    string const &dofKey = dofManager.getKey(SinglePhaseBase::viewKeyStruct::elemDofFieldString());

    // single element-based DOF field (two components per node: 0=primary,1=secondary)
    globalIndex const rankOffset = dofManager.rankOffset();

    std::vector<MeshLevel const*> meshLevelPtrs;


    forDiscretizationOnMeshTargets(domain.getMeshBodies(), [&](string const &,
                                                               MeshLevel const &mesh,
                                                               string_array const &)
    {
      meshLevelPtrs.push_back(&mesh);
    });

    if(meshLevelPtrs.size() < 2)// 先判断vector至少有2个元素，避免越界
    {
      GEOS_ERROR("The dual continuum flow solver requires at least two meshes");
    }
    else
    {
      MeshLevel const * matrixMeshPtr = meshLevelPtrs[0];
      MeshLevel const * fractureMeshPtr = meshLevelPtrs[1];

      typename TYPEOFREF( this->getStencil() ) ::KernelWrapper stencilWrapper = this->getStencil().createKernelWrapper();
      GEOS_UNUSED_VAR(stencilWrapper);
      if (Base::primarySolver()->isThermal() && Base::secondarySolver()->isThermal())
      {
        singlePhaseThermalDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
                                                  dofKey,
                                                  fluxApprox.getName(),
                                                  matrixMeshPtr->getElemManager(),
                                                  fractureMeshPtr->getElemManager(),
                                                  stencilWrapper,
                                                  this->getGravityDrainageFlag(),
                                                  this->gravityVector()[2],
                                                  this->getFracSpacingLz(),
                                                  this->getInterporosityExchangeCoefficient(),
                                                  dt,
                                                  localMatrix.toViewConstSizes(),
                                                  localRhs.toView());
      }
      else if (!(Base::primarySolver()->isThermal()) && !(Base::secondarySolver()->isThermal()))
      {
        singlePhaseDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
                                                  dofKey,
                                                  fluxApprox.getName(),
                                                  matrixMeshPtr->getElemManager(),
                                                  fractureMeshPtr->getElemManager(),
                                                  stencilWrapper,
                                                  this->getGravityDrainageFlag(),
                                                  this->gravityVector()[2],
                                                  this->getFracSpacingLz(),
                                                  this->getInterporosityExchangeCoefficient(),
                                                  dt,
                                                  localMatrix.toViewConstSizes(),
                                                  localRhs.toView());
      } else
      {
        GEOS_ERROR(
            "The dual continuum flow solver requires both primary and secondary solvers to be either thermal or non-thermal,in this problem,\n"
            "the primary solver is " << (Base::primarySolver()->isThermal() ? "thermal" : "non-thermal") << ",\n"
                                                                                                            "the secondary solver is "<< (Base::secondarySolver()->isThermal() ? "thermal" : "non-thermal") << ".\n");
      }

      // ---- Multi-porosity effective storage matrix (M_bar) cross-storage ----
      // GEOS's per-continuum BiotPorosity gives each continuum its intrinsic
      // 1/M_i; the Mehrabian (2014) multi-porosity storage (eq.22 / Appendix A)
      // requires  1/Mbar_ij = delta_ij v_i alpha_i/(B_i K_i) - abar_i abar_j/Kbar,
      // with abar_i = Kbar v_i alpha_i/K_i, Kbar = (v_m/K_m + v_f/K_f)^-1.
      // We add the correction (1/Mbar - 1/M) to the dual-flow accumulation so
      // continuum i's fluid content is coupled to BOTH pressures:
      //   row_i += rho_i V [ (1/Mbar_ii - 1/M_i)(p_i-p_i_n) + 1/Mbar_ij (p_j-p_j_n) ].
      // The off-diagonal (matrix<->fracture pressure) is essential; it lets the
      // matrix shed pressure as the fracture drains (Mandel-Cryer).  Enabled when
      // the fracture volume fraction v_f is provided on DualContinuumCrossFlow.
      real64 const v_f = this->getFractureVolumeFraction();
      if( v_f > 0.0 && this->getEnableCrossStorageCorrection() )
      {
        real64 const v_m = 1.0 - v_f;
        string_array const & matRegions  = this->template getReference< string_array >( "matrixRegionList" );
        string_array const & fracRegions = this->template getReference< string_array >( "fractureRegionList" );
        ElementRegionManager const & matEM  = matrixMeshPtr->getElemManager();
        ElementRegionManager const & fracEM = fractureMeshPtr->getElemManager();
        for( size_t ir = 0; ir < matRegions.size(); ++ir )
        {
          ElementRegionBase const & matReg  = matEM.getRegion( matRegions[ir] );
          ElementRegionBase const & fracReg = fracEM.getRegion( fracRegions[ir] );
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
            arrayView1d< real64 const > const vol = matSR.getElementVolume();

            string const & solidMn = matSR.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
            string const & solidFn = fracSR.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
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

            string const & fluidMn = matSR.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
            string const & fluidFn = fracSR.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
            constitutive::SingleFluidBase const & flM = this->template getConstitutiveModel< constitutive::SingleFluidBase >( matSR, fluidMn );
            constitutive::SingleFluidBase const & flF = this->template getConstitutiveModel< constitutive::SingleFluidBase >( fracSR, fluidFn );
            arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const densM = flM.density();
            arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const densF = flF.density();
            arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensM = flM.dDensity();
            arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensF = flF.dDensity();

            for( localIndex k = 0; k < matSR.size(); ++k )
            {
              real64 const Kbar = 1.0/( v_m/Km[k] + v_f/Kf[k] );
              real64 const abm  = Kbar*v_m*aM[k]/Km[k];
              real64 const abf  = Kbar*v_f*aF[k]/Kf[k];
              real64 const cfM  = dDensM[k][0][0]/densM[k][0];
              real64 const cfF  = dDensF[k][0][0]/densF[k][0];
              real64 const invMm = (aM[k]-phiM[k])/KsM[k] + phiM[k]*cfM;
              real64 const invMf = (aF[k]-phiF[k])/KsF[k] + phiF[k]*cfF;
              // 1/Mbar_ii = a_ii - abar_i^2/Kbar, a_ii = v_i(1/M_i + alpha_i^2/K_i)
              real64 const corrDiagM = v_m*(invMm + aM[k]*aM[k]/Km[k]) - abm*abm/Kbar - invMm;
              real64 const corrDiagF = v_f*(invMf + aF[k]*aF[k]/Kf[k]) - abf*abf/Kbar - invMf;
              real64 const corrOff   = -abm*abf/Kbar;
              real64 const dpM = pM[k]-pMn[k];
              real64 const dpF = pF[k]-pFn[k];

              if( ghostM[k] < 0 )
              {
                localIndex const row = LvArray::integerConversion< localIndex >( dofM[k]-rankOffset );
                real64 const rv = densM[k][0]*vol[k];
                localRhs[row] += rv*( corrDiagM*dpM + corrOff*dpF );
                globalIndex cols[2] = { dofM[k], dofF[k] };
                real64 vals[2] = { rv*corrDiagM, rv*corrOff };
                localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );
              }
              if( ghostF[k] < 0 )
              {
                localIndex const row = LvArray::integerConversion< localIndex >( dofF[k]-rankOffset );
                real64 const rv = densF[k][0]*vol[k];
                localRhs[row] += rv*( corrDiagF*dpF + corrOff*dpM );
                globalIndex cols[2] = { dofF[k], dofM[k] };
                real64 vals[2] = { rv*corrDiagF, rv*corrOff };
                localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );
              }
            }
          } );
        }
      }
    }
    GEOS_LOG("after the assemble");
    //this->printCRSMatrix(this->m_localMatrix);
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::DualContinuumFVM(const string &name, dataRepository::Group *parent)
      : Base(name, parent)
  {
    LinearSolverParameters &linParams = m_linearSolverParameters.get();
    linParams.multiscale.fieldName = SinglePhaseBase::viewKeyStruct::elemDofFieldString();
    linParams.multiscale.label = "dualflow";
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::postInputInitialization()
  {
    Base::postInputInitialization();
  }
/*
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::setupCoupling(DomainPartition const &domain,
                                                                                   DofManager &dofManager) const
  {
    // ensure element-based coupling (two components per element) has sparsity
    GEOS_LOG(SinglePhaseBase::viewKeyStruct::elemDofFieldString());

    // Get supports from both solvers
    stdVector< DofManager::FieldSupport > supports;
    auto const & primaryTargets = Base::primarySolver()->getMeshTargets();
    auto const & secondaryTargets = Base::secondarySolver()->getMeshTargets();
    for( auto const & p : primaryTargets )
    {
      MeshBody const & meshBody = domain.getMeshBody( p.first.first );
      MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
      std::set< string > regionNames( p.second.begin(), p.second.end() );
      supports.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
    }
    for( auto const & p : secondaryTargets )
    {
      MeshBody const & meshBody = domain.getMeshBody( p.first.first );
      MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
      std::set< string > regionNames( p.second.begin(), p.second.end() );
      supports.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
    }

    dofManager.addCouplingDualContinuum( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                                         SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                                         DofManager::Connector::Elem,
                                         supports );
  }
*/
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::initializePreSubGroups()
  {
    Base::initializePreSubGroups();
    // Ensure discretization is valid for each underlying flow solver if needed
  }
  /*
  template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::setupDofs( DomainPartition const & domain, DofManager & dofManager ) const
  {
    // Let the base DualContinuumFlowSolverBase call each sub-solver's setupDofs
    Base::setupDofs( domain, dofManager );
  }
  */
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::setupSystem(DomainPartition &domain,
                                                                                 DofManager &dofManager,
                                                                                 CRSMatrix<real64, globalIndex> &localMatrix,
                                                                                 ParallelVector &rhs,
                                                                                 ParallelVector &solution,
                                                                                 bool const setSparsity)
  {
    static bool connectivityRegistered = false;
    if( !connectivityRegistered )
    {
      this->registerMeshConnectivity( domain );
      connectivityRegistered = true;
    }

    // TracAI: Print registered connectivity values to terminal
    if( connectivityRegistered )
    {
      this->printRegisteredConnectivityValues( domain );
    }

    GEOS_MARK_FUNCTION;
    PhysicsSolverBase::setupSystem(domain, dofManager, localMatrix, rhs, solution, setSparsity);

    if (!m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct)
    {
      m_precond = createPreconditioner(domain);
    }
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  std::unique_ptr<PreconditionerBase<LAInterface>> DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::createPreconditioner(DomainPartition &domain) const
  {
    LinearSolverParameters const &linParams = m_linearSolverParameters.get();
    switch (linParams.preconditionerType)
    {
    case LinearSolverParameters::PreconditionerType::multiscale:
    {
      return std::make_unique<MultiscalePreconditioner<LAInterface>>(linParams, domain);
    }
    default:
    {
      return PhysicsSolverBase::createPreconditioner(domain);
    }
    }
  }

  /*
  template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                  real64 const dt,
                                                                                  DomainPartition & domain,
                                                                                  DofManager const & dofManager,
                                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                  arrayView1d< real64 > const & localRhs )
  {
    GEOS_MARK_FUNCTION;
    this->primarySolver()->assembleAccumulationTerms(domain,
                                                      dofManager,
                                                      localMatrix,
                                                      localRhs );
    this->secondarySolver()->assembleAccumulationTerms(domain,
                                                     dofManager,
                                                     localMatrix,
                                                     localRhs );
    if( 0 )
    {
      this->primarySolver()->assembleStabilizedFluxTerms( dt,
                                                       domain,
                                                       dofManager,
                                                       localMatrix,
                                                       localRhs );
    }
    else
    {
      this->primarySolver()->assembleFluxTerms( dt,
                                             domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );
    }
    // Step 3: compute the fluxes (face-based contributions)
    if( 0 )
    {
      this->secondarySolver()->assembleStabilizedFluxTerms( dt,
                                                       domain,
                                                       dofManager,
                                                       localMatrix,
                                                       localRhs );
    }
    else
    {
      this->secondarySolver()->assembleFluxTerms( dt,
                                             domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );
    }
  }
  */
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::updateState( geos::DomainPartition & domain )
  {
    Base::updateState(domain);
    if( this->getGravityDrainageFlag() )
    {
      // Update gravity drainage pressure for dual continuum flow
      GEOS_LOG("Updating gravity drainage pressure for dual continuum flow");
      
      // Get primary and secondary solvers
      PRIMARY_FLOW_SOLVER const * primarySolver = this->primarySolver();
      
      // Get meshes for matrix and fracture
      std::vector<MeshLevel const*> meshLevelPtrs;
      this->forDiscretizationOnMeshTargets(
        domain.getMeshBodies(),
        [&]( string const &,
             MeshLevel const & mesh,
             string_array const & )
        {
          meshLevelPtrs.push_back( &mesh );
        }
      );
      
      if( meshLevelPtrs.size() >= 2 )
      {
        MeshLevel  &  matrixMesh = const_cast<MeshLevel&>(* meshLevelPtrs[0]);
        MeshLevel  &  fractureMesh = const_cast<MeshLevel&>(* meshLevelPtrs[1]);

        real64 const gravityCoefficient = primarySolver->gravityVector()[2];
        Base::updateGravityPressure(matrixMesh,fractureMesh, gravityCoefficient);
      }
    }
  }

  // Explicit instantiation for default template
  template class DualContinuumFVM<SinglePhaseBase, SinglePhaseBase>;
  namespace
  { // Register the solver so it can be used from XML

    typedef DualContinuumFVM<SinglePhaseBase, SinglePhaseBase> DualContinuumFVM;
    REGISTER_CATALOG_ENTRY(PhysicsSolverBase, DualContinuumFVM, string const &, Group *const)
  }

} // namespace geos
