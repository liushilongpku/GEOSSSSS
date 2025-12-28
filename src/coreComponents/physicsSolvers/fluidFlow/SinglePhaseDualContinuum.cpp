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

//#include "common/TimingMacros.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::DualContinuumFVM( const string & name, dataRepository::Group * parent )
  : Base( name, parent )
{
  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.multiscale.fieldName = SinglePhaseBase::viewKeyStruct::elemDofFieldString();
  linParams.multiscale.label = "dualflow";
}

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER >
void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::postInputInitialization()
{
  Base::postInputInitialization();
}

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER >
void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                               DofManager & dofManager ) const
{
  dofManager.addCoupling( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::initializePreSubGroups()
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
template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::setupSystem( DomainPartition & domain,
                                            DofManager & dofManager,
                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                            ParallelVector & rhs,
                                            ParallelVector & solution,
                                            bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  PhysicsSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    m_precond = createPreconditioner( domain );
  }
}

template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
std::unique_ptr< PreconditionerBase< LAInterface > > DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::createPreconditioner( DomainPartition & domain ) const
{
  LinearSolverParameters const & linParams = m_linearSolverParameters.get();
  switch( linParams.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::multiscale:
    {
      return std::make_unique< MultiscalePreconditioner< LAInterface > >( linParams, domain );
    }
    default:
    {
      return PhysicsSolverBase::createPreconditioner( domain );
    }
  }
}

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


// Explicit instantiation for default template
template class DualContinuumFVM< SinglePhaseBase, SinglePhaseBase  >;
namespace
{// Register the solver so it can be used from XML

    typedef DualContinuumFVM< SinglePhaseBase, SinglePhaseBase  > DualContinuumFVM;
    REGISTER_CATALOG_ENTRY( PhysicsSolverBase, DualContinuumFVM, string const &, Group * const )
}

} // namespace geos


