
/**
 * @file PPUPhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/Layouts.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "PotGrad.hpp"
#include "PhaseComponentFlux.hpp"

namespace geos
{

namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct PPUPhaseFlux
{
  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param numPhase number of phases
   * @param ip phase index
   * @param hasCapPressure flag indicating if there is capillary pressure
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param trans transmissibility at the connection
   * @param dTrans_dPres derivative of transmissibility wrt pressure
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param phaseMob phase mobility
   * @param dPhaseMob derivative of phase mobility wrt pressure, temperature, comp density
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param phaseCapPressure phase capillary pressure
   * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   * @param dPhaseFlux_dTrans derivative of phase flux wrt transmissibility
   */
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           integer const checkPhasePresenceInGravity,
           integer const hasGravityDraingae,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           real64 const ( &trans )[2],
           real64 const ( &dTrans_dPres )[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres_m,
           ElementViewConst< arrayView1d< real64 const > > const & pres_f,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef_m,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef_f,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob_m,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob_f,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_m,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_f,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac_m,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac_f,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_m,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_f,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens_m,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens_f,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens_m,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens_f,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens_m,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens_f,
           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure_m,
           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure_f,
           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac_m,
           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac_f,
           ElementViewConst< arrayView2d< real64 const >> const & gravityDrainagePressure_m,
           real64 & potGrad,
           real64 & phaseFlux,
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp],
           real64 & dPhaseFlux_dTrans )
  {
    // assign to zero
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] = 0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] = 0;
      }
    }

    real64 dPotGrad_dTrans = 0;
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, numFluxSupportPoints >( numPhase, ip, hasCapPressure, checkPhasePresenceInGravity,hasGravityDraingae,
                                                       seri, sesri, sei, trans, dTrans_dPres, pres_m, pres_f,
                                                       gravCoef_m, gravCoef_f, phaseVolFrac_m, phaseVolFrac_f, dPhaseVolFrac_m, dPhaseVolFrac_f, dCompFrac_dCompDens_m, dCompFrac_dCompDens_f,
                                                       phaseMassDens_m, phaseMassDens_f, dPhaseMassDens_m, dPhaseMassDens_f,
                                                       phaseCapPressure_m, phaseCapPressure_f, dPhaseCapPressure_dPhaseVolFrac_m, dPhaseCapPressure_dPhaseVolFrac_f,
                                                       gravityDrainagePressure_m,
                                                       potGrad, dPotGrad_dTrans, dPresGrad_dP,
                                                       dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    // *** upwinding ***

    // choose upstream cell
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;
    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    // Use matrix properties for first support point, fracture for second
    auto const & phaseMob = (k_up == 0) ? phaseMob_m : phaseMob_f;
    auto const & dPhaseMob = (k_up == 0) ? dPhaseMob_m : dPhaseMob_f;

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

    // compute phase flux using upwind mobility
    phaseFlux = mobility * potGrad;

    dPhaseFlux_dTrans = mobility * dPotGrad_dTrans;

    // pressure gradient depends on all points in the stencil
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += mobility * (dPresGrad_dP[ke] - dGravHead_dP[ke]);
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += mobility * (dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc]);
      }
    }

    real64 const dMob_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dMob_dC = dPhaseMob[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dMob_dC[Deriv::dC+jc] * potGrad;
    }
  }
};

} // namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP
