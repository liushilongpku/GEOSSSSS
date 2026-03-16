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


#include "kernels/compositionalMultiPhase/FluxComputeKernel.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/FluxComputeKernelBase.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/C1PPUPhaseFlux.hpp"
#include "kernels/compositionalMultiPhase/ThermalFluxComputeKernel.hpp"
#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"

//#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/FluxComputeKernel.hpp"
//#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/Thermal"
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


  //TODO@LSL 扩散与弥散作用强不强？？

  GEOS_MARK_FUNCTION;

  using namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels;

  BitFlags< KernelFlags > kernelFlags;
  bool hasPrimaryCapPressure = this->primarySolver()->ishasCapPressure();
  bool hasSecondaryCapPressure = this->secondarySolver()->ishasCapPressure();
  if( hasPrimaryCapPressure || hasSecondaryCapPressure )
  {
    if( !hasPrimaryCapPressure || !hasSecondaryCapPressure )
    {
      GEOS_ERROR("Both primary and secondary solvers must have capillary pressure models if either one has it. "
                 "If you only want capillary pressure in one solver, set the capillary pressure values to zero in the other solver.");
    }
    kernelFlags.set( KernelFlags::CapPressure );
  }
  if( this->primarySolver()->ishasDiffusion()  )
    kernelFlags.set( KernelFlags::Diffusion );
  if( this->primarySolver()->ishasDispersion()  )
    kernelFlags.set( KernelFlags::Dispersion );
  if( this->primarySolver()->isuseTotalMassEquation()  )
    kernelFlags.set( KernelFlags::TotalMassEquation );
  if( this->primarySolver()->getgravityDensityScheme() == GravityDensityScheme::PhasePresence && this->primarySolver()->getgravityDensityScheme() == GravityDensityScheme::PhasePresence)
    kernelFlags.set( KernelFlags::CheckPhasePresenceInGravity );
  if( Base::getGravityDrainageFlag() )
    kernelFlags.set( KernelFlags::GravityDrainage);


//    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
//                                                                 MeshLevel const & mesh,
//                                                                 string_array const & )
//    {
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->primarySolver()->getDiscretizationName() );

  //寻找不同网格下的MeshLevel指针，后续需要传入kernel
  std::vector<MeshLevel const*> meshLevelPtrs;
  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
                                 [&](string const &,
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

    typename TYPEOFREF( this->getStencil() )::KernelWrapper stencilWrapper = this->getStencil().createKernelWrapper();
    GEOS_UNUSED_VAR( stencilWrapper );


    auto const & upwindingParams = fluxApprox.upwindingParams();
    if( upwindingParams.upwindingScheme == UpwindingScheme::C1PPU &&
        isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities::epsC1PPU > 0 )
      kernelFlags.set( KernelFlags::C1PPU );
    else if( upwindingParams.upwindingScheme == UpwindingScheme::IHU )
      kernelFlags.set( KernelFlags::IHU );

    string const & elemDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );
    /*
     *  ├── 公式类型 (m_formulationType)
        │   ├── Ov (Overall Composition)
        │   └── 其他类型
        │       ├── 热效应 (m_isThermal)
        │       │   ├── 是：计算非等温对流通量
        │       │   └── 否：计算等温对流通量
        │       │       ├── DBC (m_dbcParams.useDBC)
        │       │       │   ├── 是：带DBC的等温对流通量
        │       │       │   └── 否：标准等温对流通量
        │       └── 扩散/弥散 (m_hasDiffusion || m_hasDispersion)
        │           ├── 是：计算扩散/弥散通量
        │           │   ├── 热效应 (m_isThermal)
        │           │   │   ├── 是：非等温扩散/弥散通量
        │           │   │   └── 否：等温扩散/弥散通量
     */
//        fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
//        {
//          typename TYPEOFREF( stencil )::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

    // Convective flux

    if( false )
    {
      // isothermal only for now
      // Overall Composition
      // TODO@LSL 这也是为什么做不了多相带热的问题的原因
    }
    else
    {
      if( this->primarySolver()->isThermal() && this->secondarySolver()->isThermal() )//thermal
      {
        thermalDualContinuumCompositionalMultiPhaseCrossFlowKernels::
        FluxComputeKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( this->primarySolver()->numFluidComponents(),
                                                   this->primarySolver()->numFluidPhases(),
                                                   dofManager.rankOffset(),
                                                   dofManager.rankOffset(),
                                                   elemDofKey,
                                                   kernelFlags,
                                                   this->primarySolver()->getName(),
                                                   this->secondarySolver()->getName(),
                                                   matrixMeshPtr->getElemManager(),
                                                   fractureMeshPtr->getElemManager(),
                                                   stencilWrapper,
                                                   dt,
                                                   localMatrix.toViewConstSizes(),
                                                   localRhs.toView() );
      }
      else
      {
        if( false )
        {
          //DBC
        }
        else if( ! this->primarySolver()->isThermal() &&  ! this->secondarySolver()->isThermal())
        {
          isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels::
          FluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( this->primarySolver()->numFluidComponents(),
                                                     this->primarySolver()->numFluidPhases(),
                                                     dofManager.rankOffset(),
                                                     dofManager.rankOffset(),//TODO@LSL 裂缝与基质的总偏移量是相同的
                                                     elemDofKey,
                                                     kernelFlags,
                                                     this->primarySolver()->getName(),
                                                     this->secondarySolver()->getName(),
                                                     matrixMeshPtr->getElemManager(),
                                                     fractureMeshPtr->getElemManager(),
                                                     stencilWrapper,
                                                     dt,
                                                     localMatrix.toViewConstSizes(),
                                                     localRhs.toView() );
        }
        else
        {
          GEOS_ERROR("primiarySolver and the secondary Solver should have the same thermal");
        }
      }
    }
//        } );
  }
//    } );
  //this->printCRSMatrix(this->m_localMatrix);

}

/*
    GEOS_LOG("befor the assemble");
    this->printCRSMatrix(this->m_localMatrix);

    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR(time_n);
    GEOS_UNUSED_VAR(dt);

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
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),//TODO@LSL 添加组分数量与相的数量
                                                  dofManager.rankOffset(),
                                                  dofKey,//TODO@LSL 检查与elemdofkey差别
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
*/


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

template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
void CompositionalMultiPhaseDualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::updateState( geos::DomainPartition & domain )
{
  Base::updateState(domain);
  if( this->getGravityDrainageFlag() )
  {
    // Update gravity drainage pressure for compositional dual continuum flow
    GEOS_LOG("Updating gravity drainage pressure for compositional dual continuum flow");
    
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
      MeshLevel const * matrixMeshPtr = meshLevelPtrs[0];
      MeshLevel const * fractureMeshPtr = meshLevelPtrs[1];
      
      ElementRegionManager const & matrixElemManager = matrixMeshPtr->getElemManager();
      ElementRegionManager const & fractureElemManager = fractureMeshPtr->getElemManager();
      
      // Get gravity coefficient from primary solver
      real64 gravityCoefficient = primarySolver->gravityVector()[2];  // z-component
      MeshLevel  &  matrixMesh = const_cast<MeshLevel&>(* meshLevelPtrs[0]);
      MeshLevel  &  fractureMesh = const_cast<MeshLevel&>(* meshLevelPtrs[1]);

      Base::updateGravityPressure(matrixMesh,fractureMesh, gravityCoefficient);


      /*
      // Get fracture spacing Lz
      real64 Lz = this->getFracSpacingLz();
      
      GEOS_LOG("Gravity coefficient: " << gravityCoefficient << ", Lz: " << Lz);
      
      // TODO: Implement gravity drainage pressure update for compositional dual continuum
      // This requires accessing gravityDrainagePressure constitutive models
      // and calling updateState() with density data from both matrix and fracture regions
       */
    }
  }
}
// Explicit instantiation for default template
template class CompositionalMultiPhaseDualContinuumFVM<CompositionalMultiphaseFVM, CompositionalMultiphaseFVM>;
namespace
{ // Register the solver so it can be used from XML

typedef CompositionalMultiPhaseDualContinuumFVM<CompositionalMultiphaseFVM, CompositionalMultiphaseFVM> CompositionalMultiPhaseDualContinuumFVM;
REGISTER_CATALOG_ENTRY(PhysicsSolverBase, CompositionalMultiPhaseDualContinuumFVM, string const &, Group *const)
}

} // namespace geos
