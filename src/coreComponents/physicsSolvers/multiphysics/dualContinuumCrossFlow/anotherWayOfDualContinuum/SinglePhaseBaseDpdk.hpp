//
// Created by hello on 2025/12/22.
// This File is not in unused.

//
/**
 * @file SinglePhaseBaseDpdk.hpp
 */
#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEDPDK_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEDPDK_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{
namespace constitutive
{
  class ConstitutiveBase;
}


class SinglePhaseBaseDpdk : public SinglePhaseBase
{
public:
  SinglePhaseBaseDpdk(const string & name,
                      Group * const parent);

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
  static constexpr char const * elemDofFieldString() { return "singlePhaseVariablesDpdk"; }
  };
virtual void
applyBoundaryConditions( real64 const time_n,
                         real64 const dt,
                         DomainPartition & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) override;

virtual void
applyDirichletBC( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override;

virtual void assembleAccumulationTerms( DomainPartition & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs ) override;

template< typename SUBREGION_TYPE >
void accumulationAssemblyLaunchDPDK( DofManager const & dofManager,
                                 SUBREGION_TYPE const & subRegion,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

private:
};



template< typename SUBREGION_TYPE >
void SinglePhaseBaseDpdk::accumulationAssemblyLaunchDPDK( DofManager const & dofManager,
                                                  SUBREGION_TYPE const & subRegion,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  if( m_isThermal )
  {
    thermalSinglePhaseBaseKernels::
    AccumulationKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                               dofKey,
                                               subRegion,
                                               localMatrix,
                                               localRhs );
  }
  else
  {
    singlePhaseBaseKernels::
    AccumulationKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                               dofKey,
                                               subRegion,
                                               localMatrix,
                                               localRhs );
  }
}
}/* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEDPDK_HPP_
