/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ThomasGasOilGravityDrainagePressure.hpp"

#include "common/GeosxMacros.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ThomasGasOilGravityDrainagePressure::ThomasGasOilGravityDrainagePressure( string const & name,
                                                                          Group * const parent )
  : GravityDrainagePressureBase( name, parent )
{}

void ThomasGasOilGravityDrainagePressure::setupGravityDrainagePressureFromPhaseMassDensities(
  arrayView3d< real64 const > const matrixPhaseMassDensity,
  arrayView2d< real64 const > const matrixPhaseVolumeFraction,
  arrayView3d< real64 const > const fracturePhaseMassDensity,
  arrayView2d< real64 const > const fracturePhaseVolumeFraction,
  string_array const & phaseNames,
  real64 const gravityCoefficient,
  real64 const Lz ) const
{
  GEOS_UNUSED_VAR( matrixPhaseVolumeFraction, fracturePhaseMassDensity, fracturePhaseVolumeFraction );

  integer oilIndex = -1;
  integer gasIndex = -1;
  for( integer ip = 0; ip < LvArray::integerConversion< integer >( phaseNames.size() ); ++ip )
  {
    if( phaseNames[ip] == "oil" )
    {
      oilIndex = ip;
    }
    else if( phaseNames[ip] == "gas" )
    {
      gasIndex = ip;
    }
  }
  GEOS_THROW_IF( oilIndex < 0 || gasIndex < 0,
                 GEOS_FMT( "{} requires configured oil and gas phases.", getFullName() ), InputError );
  GEOS_THROW_IF_NE_MSG( matrixPhaseMassDensity.size( 2 ),
                        LvArray::integerConversion< localIndex >( phaseNames.size() ),
                        "Phase-density and phase-name extents must match.", InputError );
  GEOS_THROW_IF( m_gravityDrainagePressure.size( 2 ) < matrixPhaseMassDensity.size( 2 ),
                 GEOS_FMT( "{} has insufficient phase storage.", getFullName() ), InputError );

  localIndex const numElements = m_gravityDrainagePressure.size( 0 );
  localIndex const numPhases = matrixPhaseMassDensity.size( 2 );
  real64 const gravityMagnitude = std::abs( gravityCoefficient );
  real64 const drainageHeight = 0.5 * Lz;
  auto gravityDrainagePressure = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      gravityDrainagePressure[ei][0][ip] = 0.0;
    }
    gravityDrainagePressure[ei][0][oilIndex] =
      ( matrixPhaseMassDensity[ei][0][oilIndex] - matrixPhaseMassDensity[ei][0][gasIndex] )
      * gravityMagnitude * drainageHeight;
  } );
}

void ThomasGasOilGravityDrainagePressure::updateState( real64 const gravityCoefficient,
                                                        arrayView1d< real64 const > const & densityMatrix,
                                                        arrayView1d< real64 const > const & densityFracture,
                                                        localIndex const numElements,
                                                        real64 const Lz )
{
  GEOS_UNUSED_VAR( gravityCoefficient, densityMatrix, densityFracture, numElements, Lz );
  GEOS_ERROR( GEOS_FMT( "{} is only valid for compositional gas/oil flow.", getFullName() ) );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThomasGasOilGravityDrainagePressure, string const &, Group * const )

} // namespace constitutive
} // namespace geos
