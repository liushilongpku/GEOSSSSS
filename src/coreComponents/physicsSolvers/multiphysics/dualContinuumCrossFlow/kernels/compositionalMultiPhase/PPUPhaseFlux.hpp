
/**
 * @file PPUPhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/Layouts.hpp"
#include "constitutive/relativePermeability/Layouts.hpp"
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
  GEOS_HOST_DEVICE
  static localIndex selectMobilitySupport( real64 const potGrad,
                                           integer const matrixControlledExchangeUpwinding )
  {
    return matrixControlledExchangeUpwinding || potGrad >= 0.0 ? 0 : 1;
  }

  GEOS_HOST_DEVICE
  static localIndex selectCompositionUpstream( real64 const potGrad )
  {
    return potGrad >= 0.0 ? 0 : 1;
  }

  GEOS_HOST_DEVICE
  static real64 fractureCoverage( real64 const potGrad,
                                  integer const matrixControlledExchangeUpwinding,
                                  real64 const fracturePhaseVolumeFraction )
  {
    return matrixControlledExchangeUpwinding && potGrad < 0.0 ? fracturePhaseVolumeFraction : 1.0;
  }

  GEOS_HOST_DEVICE
  static real64 reverseExchangeRelPerm( real64 const configuredRelPerm,
                                        real64 const currentMatrixRelPerm )
  {
    return configuredRelPerm >= 0.0 ? configuredRelPerm : currentMatrixRelPerm;
  }

  GEOS_HOST_DEVICE
  static real64 reverseExchangeMobility( real64 const relPerm,
                                         real64 const phaseDensity,
                                         real64 const phaseViscosity,
                                         real64 const fractureCoverage )
  {
    return phaseViscosity > 0.0
           ? relPerm * phaseDensity / phaseViscosity * fractureCoverage
           : 0.0;
  }

  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           integer const checkPhasePresenceInGravity,
           integer const hasGravityDraingae,
           integer const matrixControlledExchangeUpwinding,
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
           ElementViewConst< arrayView3d< real64 const > > const & gravityDrainagePressure_m,
           real64 & potGrad,
           real64 & phaseFlux,
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp],
           real64 & dPhaseFlux_dTrans )
  {
    GEOS_UNUSED_VAR( matrixControlledExchangeUpwinding );
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] = 0.0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] = 0.0;
      }
    }

    real64 dPotGrad_dTrans = 0.0;
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, numFluxSupportPoints >(
      numPhase, ip, hasCapPressure, checkPhasePresenceInGravity, hasGravityDraingae,
      seri, sesri, sei, trans, dTrans_dPres, pres_m, pres_f, gravCoef_m, gravCoef_f,
      phaseVolFrac_m, phaseVolFrac_f, dPhaseVolFrac_m, dPhaseVolFrac_f,
      dCompFrac_dCompDens_m, dCompFrac_dCompDens_f,
      phaseMassDens_m, phaseMassDens_f, dPhaseMassDens_m, dPhaseMassDens_f,
      phaseCapPressure_m, phaseCapPressure_f,
      dPhaseCapPressure_dPhaseVolFrac_m, dPhaseCapPressure_dPhaseVolFrac_f,
      gravityDrainagePressure_m, potGrad, dPotGrad_dTrans, dPresGrad_dP,
      dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    localIndex const k_up = potGrad >= 0.0 ? 0 : 1;
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];
    auto const & phaseMob = k_up == 0 ? phaseMob_m : phaseMob_f;
    auto const & dPhaseMob = k_up == 0 ? dPhaseMob_m : dPhaseMob_f;
    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

    phaseFlux = mobility * potGrad;
    dPhaseFlux_dTrans = mobility * dPotGrad_dTrans;
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += mobility * (dPresGrad_dP[ke] - dGravHead_dP[ke]);
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += mobility * (dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc]);
      }
    }

    dPhaseFlux_dP[k_up] += dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP] * potGrad;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC+jc] * potGrad;
    }
  }

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
  computeMatrixControlled( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           integer const checkPhasePresenceInGravity,
           integer const hasGravityDraingae,
           integer const matrixControlledExchangeUpwinding,
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
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens_f,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens_f,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseVisc_f,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseVisc_f,
           ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const & phaseRelPerm_m,
           ElementViewConst< arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > > const & dPhaseRelPerm_dPhaseVolFrac_m,
           arrayView1d< real64 const > const & matrixControlledReverseExchangeRelPerm,
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
           ElementViewConst< arrayView3d< real64 const >> const & gravityDrainagePressure_m,
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

    // choose mobility support for standard PPU and matrix-to-fracture flow
    localIndex const k_up = selectMobilitySupport( potGrad, matrixControlledExchangeUpwinding );
    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    // Use matrix properties for first support point, fracture for second
    auto const & phaseMob = (k_up == 0) ? phaseMob_m : phaseMob_f;
    auto const & dPhaseMob = (k_up == 0) ? dPhaseMob_m : dPhaseMob_f;

    localIndex const er_f = seri[1];
    localIndex const esr_f = sesri[1];
    localIndex const ei_f = sei[1];
    bool const matrixControlledReverse = matrixControlledExchangeUpwinding && potGrad < 0.0;
    real64 const coverage = fractureCoverage( potGrad, matrixControlledExchangeUpwinding,
                                              phaseVolFrac_f[er_f][esr_f][ei_f][ip] );
    real64 baseMobility = phaseMob[er_up][esr_up][ei_up][ip];
    if( matrixControlledReverse )
    {
      localIndex const er_m = seri[0];
      localIndex const esr_m = sesri[0];
      localIndex const ei_m = sei[0];
      real64 const relPerm = reverseExchangeRelPerm(
        matrixControlledReverseExchangeRelPerm[ip], phaseRelPerm_m[er_m][esr_m][ei_m][0][ip] );
      baseMobility = reverseExchangeMobility(
        relPerm, phaseDens_f[er_f][esr_f][ei_f][0][ip],
        phaseVisc_f[er_f][esr_f][ei_f][0][ip], 1.0 );
    }
    real64 const mobility = baseMobility * coverage;

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

    if( !matrixControlledReverse )
    {
      real64 const dMob_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dMob_dC =
        dPhaseMob[er_up][esr_up][ei_up][ip];
      dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[k_up][jc] += dMob_dC[Deriv::dC+jc] * potGrad;
      }
    }
    else
    {
      localIndex const er_m = seri[0];
      localIndex const esr_m = sesri[0];
      localIndex const ei_m = sei[0];
      real64 const configuredRelPerm = matrixControlledReverseExchangeRelPerm[ip];
      real64 const relPerm = reverseExchangeRelPerm(
        configuredRelPerm, phaseRelPerm_m[er_m][esr_m][ei_m][0][ip] );
      real64 const density = phaseDens_f[er_f][esr_f][ei_f][0][ip];
      real64 const viscosity = phaseVisc_f[er_f][esr_f][ei_f][0][ip];
      real64 const inverseViscosity = viscosity > 0.0 ? 1.0 / viscosity : 0.0;
      real64 const pvtMobility = density * inverseViscosity;

      real64 dDens_dC[numComp]{};
      real64 dVisc_dC[numComp]{};
      applyChainRule( numComp, dCompFrac_dCompDens_f[er_f][esr_f][ei_f],
                      dPhaseDens_f[er_f][esr_f][ei_f][0][ip], dDens_dC, Deriv::dC );
      applyChainRule( numComp, dCompFrac_dCompDens_f[er_f][esr_f][ei_f],
                      dPhaseVisc_f[er_f][esr_f][ei_f][0][ip], dVisc_dC, Deriv::dC );

      real64 dRelPerm_dP = 0.0;
      real64 dRelPerm_dC[numComp]{};
      if( configuredRelPerm < 0.0 )
      {
        for( integer jp = 0; jp < numPhase; ++jp )
        {
          real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac_m[er_m][esr_m][ei_m][0][ip][jp];
          dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_m[er_m][esr_m][ei_m][jp][Deriv::dP];
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_m[er_m][esr_m][ei_m][jp][Deriv::dC+jc];
          }
        }
      }

      dPhaseFlux_dP[0] += coverage * pvtMobility * dRelPerm_dP * potGrad;
      real64 const dPvtMobility_dP =
        inverseViscosity * (dPhaseDens_f[er_f][esr_f][ei_f][0][ip][Deriv::dP] -
                            pvtMobility * dPhaseVisc_f[er_f][esr_f][ei_f][0][ip][Deriv::dP]);
      dPhaseFlux_dP[1] +=
        (coverage * relPerm * dPvtMobility_dP +
         baseMobility * dPhaseVolFrac_f[er_f][esr_f][ei_f][ip][Deriv::dP]) * potGrad;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[0][jc] += coverage * pvtMobility * dRelPerm_dC[jc] * potGrad;
        real64 const dPvtMobility_dC =
          inverseViscosity * (dDens_dC[jc] - pvtMobility * dVisc_dC[jc]);
        dPhaseFlux_dC[1][jc] +=
          (coverage * relPerm * dPvtMobility_dC +
           baseMobility * dPhaseVolFrac_f[er_f][esr_f][ei_f][ip][Deriv::dC+jc]) * potGrad;
      }
    }
  }
};

} // namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUX_HPP
