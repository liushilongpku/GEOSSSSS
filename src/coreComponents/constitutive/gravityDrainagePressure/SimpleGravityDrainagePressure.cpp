/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SimpleGravityDrainagePressure.cpp
 */

#include "SimpleGravityDrainagePressure.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

SimpleGravityDrainagePressure::SimpleGravityDrainagePressure( string const & name,
                                                              dataRepository::Group * const parent )
  : GravityDrainagePressureBase( name, parent )
{
  // No additional parameters for simple model
}

void SimpleGravityDrainagePressure::setupGravityDrainagePressure(arrayView2d< real64 const> const matrixFluidDensity, arrayView2d< real64 const > const fractureFluidDensity, real64 gravityCoefficient, real64 Lz ) const
{

  localIndex const numE = m_gravityDrainagePressure.size( 0 );
  auto gravDrainPressView = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const densityDifference = fractureFluidDensity[ei][0] - matrixFluidDensity[ei][0];
    gravDrainPressView[ei][0] = std::abs( gravityCoefficient * densityDifference * Lz / 2 );
  } );
}

void SimpleGravityDrainagePressure::setupGravityDrainagePressure(arrayView3d< real64 const> const matrixFluidDensity,
                                                                  arrayView2d< real64 const> const matrixPhaseVolumeFraction,
                                                                  arrayView2d< real64 const > const fractureTotalDensity,
                                                                  real64 gravityCoefficient,
                                                                  real64 Lz ) const
{
  // GDP is stored as a single scalar per element (phase index 0).
  // Per-phase GDP values are not currently used by the flow kernels.
  GEOS_UNUSED_VAR( matrixPhaseVolumeFraction );
  localIndex const numE = m_gravityDrainagePressure.size(0);
  auto gravDrainPressView = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    // Use the first phase's density difference as representative GDP
    gravDrainPressView[ei][0] = std::abs( gravityCoefficient * ( matrixFluidDensity[ei][0][0] - fractureTotalDensity[ei][0] ) * Lz / 2 );
  } );
}

void SimpleGravityDrainagePressure::updateState( real64 const gravityCoefficient,
                                                  arrayView1d< real64 const > const & densityMatrix,
                                                  arrayView1d< real64 const > const & densityFracture,
                                                  localIndex const numElements,
                                                  real64 const Lz )
{
  auto gravDrainPressView = m_gravityDrainagePressure.toView();
  real64 const halfLz = 0.5 * Lz;

  forAll< parallelDevicePolicy<> >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    real64 const densityDifference = densityFracture[k] - densityMatrix[k];
    gravDrainPressView[k][0] = std::abs( gravityCoefficient * densityDifference * halfLz );
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SimpleGravityDrainagePressure, string const &, Group * const )

}

}
