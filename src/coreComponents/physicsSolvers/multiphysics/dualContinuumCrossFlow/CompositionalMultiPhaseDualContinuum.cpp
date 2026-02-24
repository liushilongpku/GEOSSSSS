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

#include "CompositionalMultiPhaseDualContinuum.hpp"

// #include "common/TimingMacros.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/CrossFlowComputeKernel.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/ThermalCrossFlowComputeKernel.hpp"
namespace geos
{
  using namespace dataRepository;
  using namespace fields;

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  void CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::
    assembleCouplingTerms(real64 const time_n,
                          real64 const dt,
                          DomainPartition const &domain,
                          DofManager const &dofManager,
                          CRSMatrixView<real64, globalIndex const> const &localMatrix,
                          arrayView1d<real64> const &localRhs)
  {
    GEOS_LOG("befor the assemble");
    this->printCRSMatrix(this->m_localMatrix);

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

    string const &dofKey = dofManager.getKey(CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString());

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
        //TODO@LSL 需要针对多相流动的情况更新一个窜流项的kernel
        singlePhaseThermalDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
                                                  dofKey,
                                                  fluxApprox.getName(),
                                                  matrixMeshPtr->getElemManager(),
                                                  fractureMeshPtr->getElemManager(),
                                                  stencilWrapper,
                                                  dt,
                                                  localMatrix.toViewConstSizes(),
                                                  localRhs.toView());
      }
      else if (!(Base::primarySolver()->isThermal()) && !(Base::secondarySolver()->isThermal()))
      {
        //TODO@LSL 需要针对多相流动的情况更新一个窜流项的kernel
        singlePhaseDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
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
    }
    GEOS_LOG("after the assemble");
    this->printCRSMatrix(this->m_localMatrix);
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::CompositionalMultiPhaseDualContinuumFVM(const string &name, dataRepository::Group *parent)
      : Base(name, parent)
  {
    LinearSolverParameters &linParams = m_linearSolverParameters.get();
    linParams.multiscale.fieldName = CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString();
    linParams.multiscale.label = "dualflow";
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::postInputInitialization()
  {
    Base::postInputInitialization();
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::initializePreSubGroups()
  {
    Base::initializePreSubGroups();
    // Ensure discretization is valid for each underlying flow solver if needed
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::setupSystem(DomainPartition &domain,
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
  std::unique_ptr<PreconditionerBase<LAInterface>> CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::createPreconditioner(DomainPartition &domain) const
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


  // Explicit instantiation for default template
  template class CompositionalMultiPhaseDualContinuumFVM<CompositionalMultiphaseBase, CompositionalMultiphaseBase>;
  namespace
  { // Register the solver so it can be used from XML

    typedef CompositionalMultiPhaseDualContinuumFVM<CompositionalMultiphaseBase, CompositionalMultiphaseBase> CompositionalMultiPhaseDualContinuumFVM;
    REGISTER_CATALOG_ENTRY(PhysicsSolverBase, CompositionalMultiPhaseDualContinuumFVM, string const &, Group *const)
  }

} // namespace geos
