

/**
 * @file FluxComputeKernelBase.cpp
 */

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/FluxComputeKernelBase.hpp"

#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{
using namespace fields;

namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels
{

/******************************** FluxComputeKernelBase ********************************/

FluxComputeKernelBase::FluxComputeKernelBase( integer const numPhases,
                                              globalIndex const rankOffset_m,
                                              globalIndex const rankOffset_f,
                                              DofNumberAccessor const & dofNumberAccessor_m,
                                              CompFlowAccessors const & compFlowAccessors_m,
                                              MultiFluidAccessors const & multiFluidAccessors_m,
                                              DofNumberAccessor const & dofNumberAccessor_f,
                                              CompFlowAccessors const & compFlowAccessors_f,
                                              MultiFluidAccessors const & multiFluidAccessors_f,
                                              GravityDrainagePressureAccessors const & gravityDrainagePressureAccessors_m,
                                              real64 const dt,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs,
                                              BitFlags< KernelFlags > kernelFlags )
  : m_numPhases( numPhases ),
  m_rankOffset_m( rankOffset_m ),
  m_rankOffset_f( rankOffset_f ),
  m_dt( dt ),
  m_dofNumber_m( dofNumberAccessor_m.toNestedViewConst() ),
  m_dofNumber_f( dofNumberAccessor_f.toNestedViewConst() ),
  m_ghostRank_m( compFlowAccessors_m.get( ghostRank {} ) ),
  m_ghostRank_f( compFlowAccessors_f.get( ghostRank {} ) ),
  m_gravCoef_m( compFlowAccessors_m.get( flow::gravityCoefficient {} ) ),
  m_gravCoef_f( compFlowAccessors_f.get( flow::gravityCoefficient {} ) ),
  m_pres_m( compFlowAccessors_m.get( flow::pressure {} ) ),
  m_pres_f( compFlowAccessors_f.get( flow::pressure {} ) ),
  m_phaseVolFrac_m( compFlowAccessors_m.get( flow::phaseVolumeFraction {} ) ),
  m_phaseVolFrac_f( compFlowAccessors_f.get( flow::phaseVolumeFraction {} ) ),
  m_dPhaseVolFrac_m( compFlowAccessors_m.get( flow::dPhaseVolumeFraction {} ) ),
  m_dPhaseVolFrac_f( compFlowAccessors_f.get( flow::dPhaseVolumeFraction {} ) ),
  m_dCompFrac_dCompDens_m( compFlowAccessors_m.get( flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
  m_dCompFrac_dCompDens_f( compFlowAccessors_f.get( flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
  m_phaseCompFrac_m( multiFluidAccessors_m.get( multifluid::phaseCompFraction {} ) ),
  m_phaseCompFrac_f( multiFluidAccessors_f.get( multifluid::phaseCompFraction {} ) ),
  m_dPhaseCompFrac_m( multiFluidAccessors_m.get( multifluid::dPhaseCompFraction {} ) ),
  m_dPhaseCompFrac_f( multiFluidAccessors_f.get( multifluid::dPhaseCompFraction {} ) ),
  m_gravityDrainagePressure(gravityDrainagePressureAccessors_m.get( gravdrainage::gravityDrainagePressure {} )),
  m_localMatrix( localMatrix ),
  m_localRhs( localRhs ),
  m_kernelFlags( kernelFlags )
{}

} // namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels

} // namespace geos
