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

#include "physicsSolvers/multiphysics/dualContinuumCrossFlowComputeKernels/CrossFlowComputeKernel.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlowComputeKernels/ThermalCrossFlowComputeKernel.hpp"

namespace geos
{

  using namespace dataRepository;
  using namespace fields;

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::assembleCouplingTerms(real64 const time_n,
                                                                                           real64 const dt,
                                                                                           DomainPartition const &domain,
                                                                                           DofManager const &dofManager,
                                                                                           CRSMatrixView<real64, globalIndex const> const &localMatrix,
                                                                                           arrayView1d<real64> const &localRhs)
  {
    GEOS_MARK_FUNCTION;

    // NOTE: This is a scaffold implementation. A real dual-continuum model should:
    // - read transfer coefficients (e.g. inter-continuum transmissibility) from fields/constitutive models
    // - loop over mesh elements/faces and assemble off-diagonal blocks coupling primary<->secondary pressures
    // - add contributions to `localMatrix` and `localRhs` following the chosen discretization

    // Minimal implementation: per-element linear transfer
    //  T * (p_secondary - p_primary) added to equations as
    //  primary eqn: +T*p_primary - T*p_secondary
    //  secondary eqn: -T*p_primary + T*p_secondary
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

    if (meshLevelPtrs.size() >= 2)
    { // 关键：先判断vector至少有2个元素，避免越界
      MeshLevel const * matrixMeshPtr = meshLevelPtrs[0];
      MeshLevel const * fractureMeshPtr = meshLevelPtrs[1];

      fluxApprox.forAllStencils( * matrixMeshPtr,[&](auto &stencil)
      {
        typename TYPEOFREF(stencil)::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

        if (Base::primarySolver()->isThermal() && Base::secondarySolver()->isThermal())
        {
          GEOS_LOG("thermal conduct");
          /*
          thermalSinglePhaseDualContinuumKernels::
          CrossFlowComputeKernelFactory::
          createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                    dofKey,
                                                    fluxApprox.getName(),
                                                    matrixMeshPtr->getElemManager(),
                                                    fractureMeshPtr->getElemManager(),
                                                    stencilWrapper,
                                                    dt,
                                                    localMatrix.toViewConstSizes(),
                                                    localRhs.toView());
          */
        }//flux = sigma * transmisibility * ( p - p )
        //
        else if (!(Base::primarySolver()->isThermal()) && !(Base::secondarySolver()->isThermal()))
        {
          singlePhaseDualContinuumKernels::
          CrossFlowComputeKernelFactory::
          createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                    dofKey,
                                                    fluxApprox.getName(),
                                                    matrixMeshPtr->getElemManager(),
                                                    fractureMeshPtr->getElemManager(),
                                                    stencilWrapper,
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
      });
    }
    else
    {
      GEOS_ERROR("The dual continuum flow solver requires at least two meshes");
    }
/*

    */
    // step1:add transfer items
    //   primary eqn: +T*p_primary - T*p_secondary
    //   secondary = -primary
    //   按照ptpt来，对于p，
    //   非对角项的编号应该是：上三角(elemTotalNum * numDof +elemIdx *numDof ,elemIdx * numDof)
    //   and 下三角：             (elemIdx * numDof, elemTotalNum * numDof +elemIdx *numDof)

    // step2:fix the var of matrix and the fracture
    //  + transfer

    // 从两个solver中获取需要的pressure
    // 从基质中获得其渗透率
    // 计算面积，传导率，形状因子

    // primarySolver()->getMeshTargets();
    // secondarySolver()->getMeshTargets();

    /*
    // local copy of transfer coefficient
    real64 const T = m_transferCoefficient;
    if( std::abs( T ) < 1e-18 )
    {
      return; // nothing to assemble
    }

    this->template forDiscretizationOnMeshTargets<>( domain.getMeshBodies(), [&]( string const &,
                                                                                  MeshLevel const & mesh,
                                                                                  string_array const & regionNames )
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();

      string const localDofKey = dofKey;

      // element dof numbers accessed per subregion below
      (void) elemManager; // silence unused variable warning in some builds

      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                 CellElementSubRegion const & subRegion )
      {
        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
        arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( localDofKey );

        // serial host loop for stability of matrix updates
        forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
        {
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          globalIndex const base = dofNumber[ei];
          globalIndex const col0 = base + 0;
          globalIndex const col1 = base + 1;
          globalIndex const row0 = base - rankOffset + 0;
          globalIndex const row1 = base - rankOffset + 1;

          globalIndex cols[2] = { col0, col1 };
          real64 vals0[2] = { T, -T };
          real64 vals1[2] = { -T, T };

          localMatrix.addToRow< RAJA::seq_atomic >( row0, cols, vals0, 2 );
          localMatrix.addToRow< RAJA::seq_atomic >( row1, cols, vals1, 2 );
        } );
      } );
    } );*/
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
    // Let the base DualContinuumFlowSolver call each sub-solver's setupDofs
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

  // Explicit instantiation for default template
  template class DualContinuumFVM<SinglePhaseBase, SinglePhaseBase>;
  namespace
  { // Register the solver so it can be used from XML

    typedef DualContinuumFVM<SinglePhaseBase, SinglePhaseBase> DualContinuumFVM;
    REGISTER_CATALOG_ENTRY(PhysicsSolverBase, DualContinuumFVM, string const &, Group *const)
  }

} // namespace geos
