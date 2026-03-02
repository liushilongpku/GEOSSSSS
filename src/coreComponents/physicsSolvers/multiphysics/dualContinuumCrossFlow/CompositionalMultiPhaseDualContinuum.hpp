/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumFVM.hpp
 *
 * @brief Dual-continuum finite-volume solver that delegates to two `CompositionalMultiPhase` flow solvers.
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHASE_FLUIDFLOW_DUALCONTINUUMFVM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHASE_FLUIDFLOW_DUALCONTINUUMFVM_HPP_

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumFlowSolver.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
namespace geos
{

template< typename PRIMARY_FLOW_SOLVER = CompositionalMultiphaseFVM, typename SECONDARY_FLOW_SOLVER = CompositionalMultiphaseFVM >
class CompositionalMultiPhaseDualContinuumFVM : public DualContinuumFlowSolver<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>
{
public:

  using Base = DualContinuumFlowSolver<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>;
  using Base::forDiscretizationOnMeshTargets;
  using Base::m_discretizationName;
  using Base::m_linearSolverParameters;
  using Base::m_nonlinearSolverParameters;
  using Base::m_precond;

  static string catalogName() { return "CompositionalMultiPhaseDualContinuumFVM"; }
  string getCatalogName() const override { return catalogName(); }

  CompositionalMultiPhaseDualContinuumFVM( const string & name, dataRepository::Group * parent );
  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override;

  virtual void assembleCouplingTerms( real64 const time_n,
                                    real64 const dt,
                                    DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs ) override;

protected:

private:

};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHASE_FLUIDFLOW_DUALCONTINUUMFVM_HPP_
