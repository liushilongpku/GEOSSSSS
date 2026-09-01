/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalPerPhaseGravityDrainagePressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_COMPOSITIONALPERPHASEGRAVITYDRAINAGEPRESSURE_HPP_
#define GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_COMPOSITIONALPERPHASEGRAVITYDRAINAGEPRESSURE_HPP_

#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Compositional gravity-drainage pressure evaluated independently for each phase.
 *
 * The sign is determined by the mixture-density contrast and the magnitude uses the
 * matrix phase mass density:
 * GDP_alpha = sign(rho_m_mix-rho_f_mix) |g_z| rho_m,alpha Lz/2.
 */
class CompositionalPerPhaseGravityDrainagePressure : public GravityDrainagePressureBase
{
public:
  CompositionalPerPhaseGravityDrainagePressure( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "CompositionalPerPhaseGravityDrainagePressure"; }
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
