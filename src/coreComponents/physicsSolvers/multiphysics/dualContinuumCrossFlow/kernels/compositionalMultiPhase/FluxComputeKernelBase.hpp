

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_ISOTHERMALDUALCONTINUUMCOMPOSITIONALMULTIPHASECROSSFLOW_KERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_ISOTHERMALDUALCONTINUUMCOMPOSITIONALMULTIPHASECROSSFLOW_KERNELBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"


namespace geos
{

namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels
{

enum class KernelFlags
{
  /// Flag to specify whether capillary pressure is used or not
  CapPressure = 1 << 0, // 1
  /// Flag to specify whether diffusion is used or not
  Diffusion = 1 << 1, // 2
  /// Flag to specify whether dispersion is used or not
  Dispersion = 1 << 2, // 4
  /// Flag indicating whether total mass equation is formed or not
  TotalMassEquation = 1 << 3, // 8
  /// Flag indicating whether gravity treatment is checking phase presence or not
  CheckPhasePresenceInGravity = 1 << 4, // 16
  /// Flag indicating whether C1-PPU is used or not
  C1PPU = 1 << 5, // 32
  /// Flag indicating whether IHU is used or not
  IHU = 1 << 6, // 64
  /// Flag indicating whether HU 2-phase simplified version is used or not
  HU2PH = 1 << 7 // 128
};

/******************************** FluxComputeKernelBase ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 */
class FluxComputeKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::gravityCoefficient,
                      fields::flow::pressure,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::phaseMobility,
                      fields::flow::dPhaseMobility >;
  using MultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::dPhaseMassDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using CapPressureAccessors =
    StencilMaterialAccessors< constitutive::CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] m_rankOffset_m the offset of my MPI rank for matrix
   * @param[in] f_rankOffset_f the offset of my MPI rank for fracture
   * @param[in] dofNumberAccessor_m accessor for the dof numbers of matrix
   * @param[in] compFlowAccessors_m accessor for wrappers registered by the solver for matrix
   * @param[in] multiFluidAccessors_m accessor for wrappers registered by the multifluid model for matrix
   * @param[in] dofNumberAccessor_f accessor for the dof numbers of fracture
   * @param[in] compFlowAccessors_f accessor for wrappers registered by the solver for fracture
   * @param[in] multiFluidAccessors_f accessor for wrappers registered by the multifluid model for fracture
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed all together
   */
  FluxComputeKernelBase( integer const numPhases,
                         globalIndex const m_rankOffset_m,
                         globalIndex const m_rankOffset_f,
                         DofNumberAccessor const & dofNumberAccessor_m,
                         CompFlowAccessors const & compFlowAccessors_m,
                         MultiFluidAccessors const & multiFluidAccessors_m,
                         DofNumberAccessor const & dofNumberAccessor_f,
                         CompFlowAccessors const & compFlowAccessors_f,
                         MultiFluidAccessors const & multiFluidAccessors_f, 
                         real64 const dt,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         BitFlags< KernelFlags > kernelFlags );

protected:

  /// Number of fluid phases
  //裂缝与基质的相的数量是一样的
  //TODO@LSL 后续可以做不一样的，是一个可行的创新
  integer const m_numPhases;

  /// Offset for my MPI rank
  //_m 表示矩阵的偏移
  //_f 表示裂缝的偏移
  globalIndex const m_rankOffset_m;
  globalIndex const m_rankOffset_f;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber_m;
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber_f;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank_m;
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank_f;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef_m;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef_f;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres_m;
  ElementViewConst< arrayView1d< real64 const > > const m_pres_f;

  /// Views on phase volume fractions
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseVolFrac_m;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseVolFrac_f;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac_m;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac_f;

  /// Views on derivatives of comp fractions
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens_m;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens_f;

  /// Views on phase component fractions
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const m_phaseCompFrac_m;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const m_phaseCompFrac_f;
  ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac_m;
  ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac_f;

  // Residual and jacobian
  // 矩阵与向量不需要区分，因为是同一个矩阵

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< KernelFlags > const m_kernelFlags;
};

//TODO@LSL 这里需要对双重介质网格作适配
namespace helpers
{
template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void calculateMeanDensity( localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  integer const checkPhasePresenceInGravity,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac_m,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac_f,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens_m,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens_f,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens_m,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens_f,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens_m,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens_f,
                                  real64 & densMean, real64 (& dDensMean_dPres)[numFluxSupportPoints], real64 (& dDensMean_dComp)[numFluxSupportPoints][numComp] )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  densMean = 0;
  integer denom = 0;
  real64 dDens_dC[numComp]{};
  for( localIndex i = 0; i < numFluxSupportPoints; ++i )
  {
    localIndex const er = seri[i];
    localIndex const esr = sesri[i];
    localIndex const ei = sei[i];

    bool const phaseExists = (i == 0 ? phaseVolFrac_m[er][esr][ei][ip] > 0 : phaseVolFrac_f[er][esr][ei][ip] > 0);
    if( checkPhasePresenceInGravity && !phaseExists )
    {
      dDensMean_dPres[i] = 0.0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dComp[i][jc] = 0.0;
      }
      continue;
    }

    // density
    real64 const density = (i == 0 ? phaseMassDens_m[er][esr][ei][0][ip] : phaseMassDens_f[er][esr][ei][0][ip]);
    real64 const dDens_dPres = (i == 0 ? dPhaseMassDens_m[er][esr][ei][0][ip][Deriv::dP] : dPhaseMassDens_f[er][esr][ei][0][ip][Deriv::dP]);

    applyChainRule( numComp,
                    (i == 0 ? dCompFrac_dCompDens_m[er][esr][ei] : dCompFrac_dCompDens_f[er][esr][ei]),
                    (i == 0 ? dPhaseMassDens_m[er][esr][ei][0][ip] : dPhaseMassDens_f[er][esr][ei][0][ip]),
                    dDens_dC,
                    Deriv::dC );

    // average density and derivatives
    densMean += density;
    dDensMean_dPres[i] = dDens_dPres;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dDensMean_dComp[i][jc] = dDens_dC[jc];
    }
    denom++;
  }
  if( denom > 1 )
  {
    densMean /= denom;
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      dDensMean_dPres[i] /= denom;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dComp[i][jc] /= denom;
      }
    }
  }
}

}

} // namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMCOMPOSITIONALMULTIPHASECROSSFLOW_KERNELBASE_HPP
