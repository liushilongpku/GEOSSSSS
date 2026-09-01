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
    // SIGNED gravity drainage pressure (Kazemi). Do NOT take the absolute value:
    // the sign carries the drive direction. With gravityCoefficient = g_z (< 0) and
    // (rho_f - rho_m), this evaluates to |g|*(rho_m - rho_f)*Lz/2, i.e. > 0 when the
    // matrix fluid is denser than the fracture fluid. A positive value is added to the
    // matrix-side potential (see PotGrad), so the denser matrix fluid drains to the
    // fracture (and the lighter fracture fluid imbibes into the matrix) - the
    // conventional gravity-drainage direction. When the contrast reverses, so does GDP.
    real64 const densityDifference = fractureFluidDensity[ei][0] - matrixFluidDensity[ei][0];
    gravDrainPressView[ei][0][0] = gravityCoefficient * densityDifference * Lz / 2;
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
    // SIGNED GDP, consistent with the (rho_f - rho_m) convention used elsewhere.
    // NOTE: this 3-arg overload only has the matrix per-phase density available, so it
    // falls back to the first phase as a representative density. The compositional
    // setup path instead passes matrix/fracture TOTAL densities to the 2-arg overload,
    // which is dimensionally consistent and should be preferred.
    gravDrainPressView[ei][0][0] =
      gravityCoefficient * ( fractureTotalDensity[ei][0] - matrixFluidDensity[ei][0][0] ) * Lz / 2;
  } );
}

void SimpleGravityDrainagePressure::setupGravityDrainagePressureFromPhaseMassDensities(
  arrayView3d< real64 const > const matrixPhaseMassDensity,
  arrayView2d< real64 const > const matrixPhaseVolumeFraction,
  arrayView3d< real64 const > const fracturePhaseMassDensity,
  arrayView2d< real64 const > const fracturePhaseVolumeFraction,
  string_array const & phaseNames,
  real64 gravityCoefficient,
  real64 Lz ) const
{
  GEOS_UNUSED_VAR( phaseNames );
  localIndex const numE = m_gravityDrainagePressure.size( 0 );
  localIndex const numPhase = matrixPhaseMassDensity.size( 2 );
  auto gravDrainPressView = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 matrixMassDensity = 0.0;
    real64 fractureMassDensity = 0.0;
    for( localIndex ip = 0; ip < numPhase; ++ip )
    {
      matrixMassDensity += matrixPhaseVolumeFraction[ei][ip] * matrixPhaseMassDensity[ei][0][ip];
      fractureMassDensity += fracturePhaseVolumeFraction[ei][ip] * fracturePhaseMassDensity[ei][0][ip];
    }

    // Keep the original scalar closure available to every phase in the phase-aware field.
    real64 const gravityDrainagePressure =
      gravityCoefficient * ( fractureMassDensity - matrixMassDensity ) * Lz / 2;
    for( localIndex ip = 0; ip < gravDrainPressView.size( 2 ); ++ip )
    {
      gravDrainPressView[ei][0][ip] = gravityDrainagePressure;
    }
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
    // SIGNED GDP (see setupGravityDrainagePressure): the sign encodes the drive direction.
    real64 const densityDifference = densityFracture[k] - densityMatrix[k];
    gravDrainPressView[k][0][0] = gravityCoefficient * densityDifference * halfLz;
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SimpleGravityDrainagePressure, string const &, Group * const )

}

}
