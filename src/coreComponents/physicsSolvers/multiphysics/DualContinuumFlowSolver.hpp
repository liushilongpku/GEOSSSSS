/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumFlowSolver.hpp
 *
 * @brief A coupled solver that binds two flow solvers for dual-continuum/dual-porosity-style models.
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseDpdk.hpp"
#include "mesh/DomainPartition.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{
    namespace stabilization
    {
        enum class StabilizationType : integer
        {
            None,
            Global,
            Local,
        };

        ENUM_STRINGS( StabilizationType,
                      "None",
                      "Global",
                      "Local" );
    }

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER >
class DualContinuumFlowSolver : public CoupledSolver< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >
{
public:

  using Base = CoupledSolver< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Primary = 0,
    Secondary = 1
  };

  static string coupledSolverAttributePrefix() { return "dualcontinuum"; }

  DualContinuumFlowSolver( string const & name, dataRepository::Group * parent )
    : Base( name, parent )
  {
    // dual-flow has two scalar dofs per node (primary + secondary), set a reasonable default
    LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();
    linearSolverParameters.dofsPerNode = 2;
    linearSolverParameters.multiscale.label = "dualflow";
    // default transfer coefficient (per-element scalar). Users may override in input

    // m_transferCoefficient = 0.0;

    this->registerWrapper( viewKeyStruct::transferCoefficientString(), &m_transferCoefficient ).
      setApplyDefaultValue( 0.0 ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "Inter-continuum transfer coefficient (per-element scalar, default 0)" );

  }

  virtual void postInputInitialization() override {
    Base::postInputInitialization();
    // add validation of solver types if needed
    GEOS_LOG("some thermal check");
  }

  /**
   * @brief accessor for the pointer to the primary flow solver
   * @return a pointer to the primary flow solver
   */
  PRIMARY_FLOW_SOLVER * primarySolver() const
  {
    return std::get< toUnderlying( SolverType::Primary ) >( m_solvers );
  }
  /**
   * @brief accessor for the pointer to the secondary flow solver
   * @return a pointer to the secondary flow solver
   */
  SECONDARY_FLOW_SOLVER * secondarySolver() const
  {
    return std::get< toUnderlying( SolverType::Secondary ) >( m_solvers );
  }

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    // ensure both flow solvers register their dofs
    primarySolver()->setupDofs( domain, dofManager );
    secondarySolver()->setupDofs( domain, dofManager );

    // allow base to set any additional coupling DOFs/attributes
    this->setupCoupling( domain, dofManager );
  }

  virtual void setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                              DofManager & dofManager ) const override
  {
    // ensure element-based coupling (two components per element) has sparsity
    dofManager.addCoupling( SinglePhaseBaseDpdk::viewKeyStruct::elemDofFieldString(),
                            SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                            DofManager::Connector::Elem );
  }

  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    GEOS_MARK_FUNCTION;

    GEOS_LOG_LEVEL_RANK_0( logInfo::Coupling, GEOS_FMT( "{}: assembling dual-continuum coupling terms", this->getName() ) );

    // NOTE: This is a scaffold implementation. A real dual-continuum model should:
    // - read transfer coefficients (e.g. inter-continuum transmissibility) from fields/constitutive models
    // - loop over mesh elements/faces and assemble off-diagonal blocks coupling primary<->secondary pressures
    // - add contributions to `localMatrix` and `localRhs` following the chosen discretization

    // Minimal implementation: per-element linear transfer
    //  T * (p_secondary - p_primary) added to equations as
    //  primary eqn: +T*p_primary - T*p_secondary
    //  secondary eqn: -T*p_primary + T*p_secondary

    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );

    // single element-based DOF field (two components per node: 0=primary,1=secondary)
    string const dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
    globalIndex const rankOffset = dofManager.rankOffset();

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
    } );
  }

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override
  {
    // Apply solution for both primary and secondary solvers to their respective meshes
    primarySolver()->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
    secondarySolver()->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
  }

  virtual void applyBoundaryConditions( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override
  {
    // Apply boundary conditions through both solvers since they target different meshes
    // This prevents duplicate boundary condition applications in multi-mesh coupled solvers.
    // Previously, only the primary solver was called, leading to missing BCs for secondary mesh.
    // Now, both primary and secondary solvers apply their respective BCs to avoid NaN solutions.
    primarySolver()->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    secondarySolver()->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

protected:

    template< typename TYPE_LIST,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
    real64 assemblyLaunch( MeshLevel & mesh,
                           DofManager const & dofManager,
                           string_array const & regionNames,
                           string const & materialNamesString,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           real64 const dt,
                           PARAMS && ... params )
    {
      GEOS_MARK_FUNCTION;
      //additional mesh manger is needed
      NodeManager const & nodeManager = mesh.getNodeManager();

      string const dofKey = dofManager.getKey( fields::flow::pressure::key() );
      arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

      real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( this->gravityVector() );

      KERNEL_WRAPPER kernelWrapper( dofNumber,
                                    dofManager.rankOffset(),
                                    localMatrix,
                                    localRhs,
                                    dt,
                                    gravityVectorData,
                                    std::forward< PARAMS >( params )... );
      return 0;
      /*
      return finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< >,
              TYPE_LIST >( mesh,
                           regionNames,
                           this->solidMechanicsSolver()->getDiscretizationName(),
                           materialNamesString,
                           kernelWrapper );*/
    }

    stabilization::StabilizationType m_stabilizationType;

private:
  real64 m_transferCoefficient;

  struct viewKeyStruct : Base::viewKeyStruct
  {
    static constexpr char const * transferCoefficientString() { return "transferCoefficient"; }
    // Add any DualContinuum-specific keys here if needed in the future
  };

}; // class DualContinuumFlowSolver

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_
