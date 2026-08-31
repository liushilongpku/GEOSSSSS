/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_THOMASGASOILGRAVITYDRAINAGEPRESSURE_HPP_
#define GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_THOMASGASOILGRAVITYDRAINAGEPRESSURE_HPP_

#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Thomas gas/oil gravity balance for a dual-continuum matrix block.
 *
 * The unresolved gas/oil gravity pressure is assigned to the oil phase:
 * GDP_o = (rho_o-rho_g) |g_z| Lz/2, while gas and water GDP are zero.
 */
class ThomasGasOilGravityDrainagePressure : public GravityDrainagePressureBase
{
public:
  ThomasGasOilGravityDrainagePressure( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "ThomasGasOilGravityDrainagePressure"; }
  string getCatalogName() const override { return catalogName(); }

  void setupGravityDrainagePressureFromPhaseMassDensities(
    arrayView3d< real64 const > const matrixPhaseMassDensity,
    arrayView2d< real64 const > const matrixPhaseVolumeFraction,
    arrayView3d< real64 const > const fracturePhaseMassDensity,
    arrayView2d< real64 const > const fracturePhaseVolumeFraction,
    string_array const & phaseNames,
    real64 gravityCoefficient,
    real64 Lz ) const override;

  void updateState( real64 gravityCoefficient,
                    arrayView1d< real64 const > const & densityMatrix,
                    arrayView1d< real64 const > const & densityFracture,
                    localIndex numElements,
                    real64 Lz ) override;
};

} // namespace constitutive
} // namespace geos

#endif
