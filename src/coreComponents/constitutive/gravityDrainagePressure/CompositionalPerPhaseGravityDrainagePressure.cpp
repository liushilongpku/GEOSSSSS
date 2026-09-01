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
 * @file CompositionalPerPhaseGravityDrainagePressure.cpp
 */

#include "CompositionalPerPhaseGravityDrainagePressure.hpp"

#include "common/GeosxMacros.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

CompositionalPerPhaseGravityDrainagePressure::CompositionalPerPhaseGravityDrainagePressure(
  string const & name, Group * const parent )
  : GravityDrainagePressureBase( name, parent )
{}

void CompositionalPerPhaseGravityDrainagePressure::setupGravityDrainagePressureFromPhaseMassDensities(
  arrayView3d< real64 const > const matrixPhaseMassDensity,
  arrayView2d< real64 const > const matrixPhaseVolumeFraction,
  arrayView3d< real64 const > const fracturePhaseMassDensity,
  arrayView2d< real64 const > const fracturePhaseVolumeFraction,
  string_array const & phaseNames,
  real64 const gravityCoefficient,
  real64 const Lz ) const
{
  GEOS_UNUSED_VAR( phaseNames );
  GEOS_THROW_IF_NE_MSG( matrixPhaseMassDensity.size( 2 ), fracturePhaseMassDensity.size( 2 ),
                        "Matrix and fracture phase-density extents must match.", InputError );
  GEOS_THROW_IF_NE_MSG( matrixPhaseMassDensity.size( 2 ), matrixPhaseVolumeFraction.size( 1 ),
                        "Matrix phase-density and saturation extents must match.", InputError );
  GEOS_THROW_IF_NE_MSG( fracturePhaseMassDensity.size( 2 ), fracturePhaseVolumeFraction.size( 1 ),
                        "Fracture phase-density and saturation extents must match.", InputError );
  GEOS_THROW_IF( m_gravityDrainagePressure.size( 2 ) < matrixPhaseMassDensity.size( 2 ),
                 GEOS_FMT( "{} has insufficient phase storage.", getFullName() ), InputError );

  localIndex const numElements = m_gravityDrainagePressure.size( 0 );
  localIndex const numPhases = matrixPhaseMassDensity.size( 2 );
  real64 const gravityMagnitude = std::abs( gravityCoefficient );
  real64 const phaseHeight = 0.5 * Lz;
  auto gravityDrainagePressure = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 matrixMixtureDensity = 0.0;
    real64 fractureMixtureDensity = 0.0;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      matrixMixtureDensity += matrixPhaseVolumeFraction[ei][ip] * matrixPhaseMassDensity[ei][0][ip];
      fractureMixtureDensity += fracturePhaseVolumeFraction[ei][ip] * fracturePhaseMassDensity[ei][0][ip];
    }

    real64 const densityContrast = matrixMixtureDensity - fractureMixtureDensity;
    real64 const direction = densityContrast > 0.0 ? 1.0 : densityContrast < 0.0 ? -1.0 : 0.0;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      gravityDrainagePressure[ei][0][ip] = direction * gravityMagnitude *
                                           matrixPhaseMassDensity[ei][0][ip] * phaseHeight;
    }
  } );
}

void CompositionalPerPhaseGravityDrainagePressure::updateState(
  real64 const gravityCoefficient,
  arrayView1d< real64 const > const & densityMatrix,
  arrayView1d< real64 const > const & densityFracture,
  localIndex const numElements,
  real64 const Lz )
{
  GEOS_UNUSED_VAR( gravityCoefficient, densityMatrix, densityFracture, numElements, Lz );
  GEOS_ERROR( GEOS_FMT( "{} requires compositional phase mass densities for initialization.", getFullName() ) );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalPerPhaseGravityDrainagePressure, string const &, Group * const )

} // namespace constitutive
} // namespace geos
