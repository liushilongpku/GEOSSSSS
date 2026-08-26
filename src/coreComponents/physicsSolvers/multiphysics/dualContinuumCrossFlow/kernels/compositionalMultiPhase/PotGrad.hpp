
/**
 * @file PotGrad.hpp
 */

#ifndef GEOS_CROSSFLOW_COMPOSITIONAL_POTGRAD_HPP
#define GEOS_CROSSFLOW_COMPOSITIONAL_POTGRAD_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/Layouts.hpp"
#include "mesh/ElementRegionManager.hpp"

#include <algorithm>
#include <cmath>


namespace geos
{

namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct PotGrad
{
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute ( integer const numPhase,
            integer const ip,
            integer const hasCapPressure,
            integer const checkPhasePresenceInGravity,
            integer const hasGravityDraingae,
            localIndex const ( &seri )[numFluxSupportPoints],
            localIndex const ( &sesri )[numFluxSupportPoints],
            localIndex const ( &sei )[numFluxSupportPoints],
            real64 const ( &trans )[numFluxSupportPoints],
            real64 const ( &dTrans_dPres )[numFluxSupportPoints],
            ElementViewConst< arrayView1d< real64 const > > const & pres_m,
            ElementViewConst< arrayView1d< real64 const > > const & pres_f,
            ElementViewConst< arrayView1d< real64 const > > const & gravCoef_m,
            ElementViewConst< arrayView1d< real64 const > > const & gravCoef_f,
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
            ElementViewConst< arrayView2d< real64 >> const & gravityDrainagePressure_m,
            real64 & potGrad,
            real64 & dPotGrad_dTrans,
            real64 ( & dPresGrad_dP )[numFluxSupportPoints],
            real64 ( & dPresGrad_dC )[numFluxSupportPoints][numComp],
            real64 ( & dGravHead_dP )[numFluxSupportPoints],
            real64 ( & dGravHead_dC )[numFluxSupportPoints][numComp] )
  {
    // assign derivatives arrays to zero
    for( integer i = 0; i < numFluxSupportPoints; ++i )
    {
      dPresGrad_dP[i] = 0.0;
      dGravHead_dP[i] = 0.0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] = 0.0;
        dGravHead_dC[i][jc] = 0.0;
      }
    }

    real64 presGrad = 0.0;
    real64 dPresGrad_dTrans = 0.0;
    real64 gravHead = 0.0;
    real64 dGravHead_dTrans = 0.0;
    real64 dCapPressure_dC[numComp]{};

    // create local work arrays
    real64 densMean = 0.0;
    real64 dDensMean_dP[numFluxSupportPoints]{};
    real64 dDensMean_dC[numFluxSupportPoints][numComp]{};
    isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels::helpers::
      calculateMeanDensity( ip, seri, sesri, sei,
                            checkPhasePresenceInGravity,
                            phaseVolFrac_m, phaseVolFrac_f,
                            dCompFrac_dCompDens_m, dCompFrac_dCompDens_f,
                            phaseMassDens_m, phaseMassDens_f,
                            dPhaseMassDens_m, dPhaseMassDens_f,
                            densMean, dDensMean_dP, dDensMean_dC );

    /// compute the TPFA potential difference
    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // capillary pressure
      real64 capPressure     = 0.0;
      real64 dCapPressure_dP = 0.0;

      for( integer ic = 0; ic < numComp; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( hasCapPressure )
      {
        if( i == 0 )
        {
          capPressure = phaseCapPressure_m[er][esr][ei][0][ip];//TODO@LSL 裂缝的 capillary pressure与基质不同

          for( integer jp = 0; jp < numPhase; ++jp )
          {
            real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac_m[er][esr][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_m[er][esr][ei][jp][Deriv::dP];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_m[er][esr][ei][jp][Deriv::dC+jc];
            }
          }
        }
        else
        {
          capPressure = phaseCapPressure_f[er][esr][ei][0][ip];

          for( integer jp = 0; jp < numPhase; ++jp )
          {
            real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac_f[er][esr][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_f[er][esr][ei][jp][Deriv::dP];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_f[er][esr][ei][jp][Deriv::dC+jc];
            }
          }
        }
      }

      if( i==0 ) //判断是应该使用那种物质中的场
      {
        real64 gravityDrainagePressure = 0.0;
        real64 dGDP_dP[numFluxSupportPoints]{};
        real64 dGDP_dC[numFluxSupportPoints][numComp]{};
        if( hasGravityDraingae )
        {
          gravityDrainagePressure = gravityDrainagePressure_m[er][esr][ei][0];
          // GDP = g_z * Lz/2 * (rho_f_mix - rho_mix), where
          //   rho_side_mix = sum_ip S_side,ip * rho_side,ip.
          // Recover the local density-difference scale from the stored GDP value to
          // build the implicit dGDP/dP and dGDP/dC Jacobian contributions. The scale
          // is protected against a vanishing density difference (GDP close to zero),
          // where the derivative is negligible anyway.
          real64 rho_m_mix = 0.0;
          real64 rho_f_mix = 0.0;
          for( integer jp = 0; jp < numPhase; ++jp )
          {
            rho_m_mix += phaseVolFrac_m[er][esr][ei][jp] * phaseMassDens_m[er][esr][ei][0][jp];
            rho_f_mix += phaseVolFrac_f[seri[1]][sesri[1]][sei[1]][jp] * phaseMassDens_f[seri[1]][sesri[1]][sei[1]][0][jp];
          }
          real64 const densityDiff = rho_f_mix - rho_m_mix;
          real64 gdpScale = 0.0;
          real64 const tolScale = std::max( std::fabs( rho_f_mix ), std::fabs( rho_m_mix ) );
          if( std::fabs( densityDiff ) > 1.0e-12 * std::max( 1.0, tolScale ) )
          {
            gdpScale = gravityDrainagePressure / densityDiff;
          }
          // dGDP/dx = gdpScale * d(rho_f_mix - rho_mix)/dx
          // matrix side (support point 0): d rho_m_mix/dx
          for( integer jp = 0; jp < numPhase; ++jp )
          {
            real64 const s_m = phaseVolFrac_m[er][esr][ei][jp];
            real64 const rho_m = phaseMassDens_m[er][esr][ei][0][jp];
            dGDP_dP[0] += -gdpScale * ( dPhaseVolFrac_m[er][esr][ei][jp][Deriv::dP] * rho_m
                                        + s_m * dPhaseMassDens_m[er][esr][ei][0][jp][Deriv::dP] );
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dGDP_dC[0][jc] += -gdpScale * ( dPhaseVolFrac_m[er][esr][ei][jp][Deriv::dC+jc] * rho_m
                                              + s_m * dPhaseMassDens_m[er][esr][ei][0][jp][Deriv::dC+jc] );
            }
          }
          // fracture side (support point 1): d rho_f_mix/dx
          for( integer jp = 0; jp < numPhase; ++jp )
          {
            real64 const s_f = phaseVolFrac_f[seri[1]][sesri[1]][sei[1]][jp];
            real64 const rho_f = phaseMassDens_f[seri[1]][sesri[1]][sei[1]][0][jp];
            dGDP_dP[1] += gdpScale * ( dPhaseVolFrac_f[seri[1]][sesri[1]][sei[1]][jp][Deriv::dP] * rho_f
                                       + s_f * dPhaseMassDens_f[seri[1]][sesri[1]][sei[1]][0][jp][Deriv::dP] );
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dGDP_dC[1][jc] += gdpScale * ( dPhaseVolFrac_f[seri[1]][sesri[1]][sei[1]][jp][Deriv::dC+jc] * rho_f
                                             + s_f * dPhaseMassDens_f[seri[1]][sesri[1]][sei[1]][0][jp][Deriv::dC+jc] );
            }
          }
        }
        // GDP is added to matrix-side pressure (Kazemi model: P_grav = |g * (rho_f - rho_m) * Lz / 2|).
        // GDP is now included implicitly: its dGDP/dP and dGDP/dC are added to the
        // matrix/fracture pressure and component Jacobian entries below.
        real64 const dP = pres_m[er][esr][ei] - capPressure + gravityDrainagePressure;
        presGrad += trans[i] * dP;
        dPresGrad_dTrans += dP;
        dPresGrad_dP[i] += trans[i] * ( 1 - dCapPressure_dP + dGDP_dP[i] ) + dTrans_dPres[i] * dP;
        // GDP also depends on the fracture-side (support point 1) state.
        dPresGrad_dP[1] += trans[i] * dGDP_dP[1];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPresGrad_dC[i][jc] += -trans[i] * dCapPressure_dC[jc] + trans[i] * dGDP_dC[i][jc];
        }
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPresGrad_dC[1][jc] += trans[i] * dGDP_dC[1][jc];
        }

        real64 const gC = gravCoef_m[er][esr][ei];
        real64 const gravD = trans[i] * gC;
        real64 const dGravD_dTrans = gC;
        real64 const dGravD_dP = dTrans_dPres[i] * gC;

        // the density used in the potential difference is always a mass density
        // unlike the density used in the phase mobility, which is a mass density
        // if useMass == 1 and a molar density otherwise
        gravHead += densMean * gravD;
        dGravHead_dTrans += densMean * dGravD_dTrans;

        // need to add contributions from both cells the mean density depends on
        for( integer j = 0; j < numFluxSupportPoints; ++j )
        {
          dGravHead_dP[j] += dDensMean_dP[j] * gravD;
          if( j == i )
          {
            dGravHead_dP[j] += dGravD_dP * densMean;
          }
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
          }
        }
      }
      else if(i == 1)//TODO@LSL 当裂缝网格与基质网格的网格中心不同时，会产生由重力带来的差异，需要检查这种差异
      {
        //std::cout << "phase " << ip << " in block "<< i << " has pressure of " << pres_m[er][esr][ei]<< " and the cappres is :"<< capPressure<< std::endl;
        real64 const dP = pres_f[er][esr][ei] - capPressure;
        presGrad += -trans[i] * dP;
        dPresGrad_dTrans -= dP;
        dPresGrad_dP[i] += -trans[i] * ( 1 - dCapPressure_dP ) - dTrans_dPres[i] * dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPresGrad_dC[i][jc] += trans[i] * dCapPressure_dC[jc];
        }

        real64 const gC = gravCoef_f[er][esr][ei];
        real64 const gravD = -trans[i] * gC;
        real64 const dGravD_dTrans = -gC;
        real64 const dGravD_dP = -dTrans_dPres[i] * gC;

        // the density used in the potential difference is always a mass density
        // unlike the density used in the phase mobility, which is a mass density
        // if useMass == 1 and a molar density otherwise
        gravHead += densMean * gravD;
        dGravHead_dTrans += densMean * dGravD_dTrans;

        // need to add contributions from both cells the mean density depends on
        for( integer j = 0; j < numFluxSupportPoints; ++j )
        {
          dGravHead_dP[j] += dDensMean_dP[j] * gravD;
          if( j == i )
          {
            dGravHead_dP[j] += dGravD_dP * densMean;
          }
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
          }
        }
      }
      else
      {
        GEOS_ERROR("there should be 2 element in crossflow stencils");
      }
    }

    // compute phase potential gradient
    potGrad = presGrad - gravHead;
    dPotGrad_dTrans = dPresGrad_dTrans - dGravHead_dTrans;
  }

};

} // namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities

} // namespace geos

#endif // GEOS_CROSSFLOW_COMPOSITIONAL_POTGRAD_HPP
