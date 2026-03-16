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

void SimpleGravityDrainagePressure::setupGravityDrainagePressure(arrayView2d< real64 const> const matrixFluidDensity, arrayView2d< real64 const > const fractureFluidDensity, real64 Lz ) const
{

  // 这里实际上不用做什么，因为在allocateConstitutiveData分配数据的时候已经有了初始值,
  // 在此仅是对初始化方式的一种测试。
  localIndex const numE = m_gravityDrainagePressure.size();
  auto gravDrainPressView = m_gravityDrainagePressure.toView();

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for(localIndex phaseIndex = 0; phaseIndex < matrixFluidDensity[ei].size(); phaseIndex++ )
    {
      gravDrainPressView[ei][phaseIndex] = std::abs( matrixFluidDensity[ei][0] - fractureFluidDensity[ei][0] ) * Lz / 2;
    }
  } );
}

void SimpleGravityDrainagePressure::setupGravityDrainagePressure(arrayView3d< real64 const> const matrixFluidDensity,
                                                                 arrayView2d< real64 const> const matrixPhaseVolumeFraction,
                                                                 arrayView2d< real64 const > const fractureTotalDensity,
                                                                 real64 gravityCoefficient,
                                                                 real64 Lz ) const
{

  // 这里实际上不用做什么，因为在allocateConstitutiveData分配数据的时候已经有了初始值,
  // 在此仅是对初始化方式的一种测试。
  localIndex const numE = m_gravityDrainagePressure.size(0);
  auto gravDrainPressView = m_gravityDrainagePressure.toView();
  localIndex const phaseNumber = matrixFluidDensity.size(2);
  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for(localIndex phaseIndex = 0; phaseIndex < phaseNumber; phaseIndex++ )
    {
      //fractureTotalDensity 尽管是个二维的数据但是第二维度长度是1
      //gravDrainPressView[ei][phaseIndex] = gravityCoefficient * matrixPhaseVolumeFraction[ei][phaseIndex] * (std::abs( matrixFluidDensity[ei][0][phaseIndex] - fractureTotalDensity[ei][0] ) * Lz / 2);
      gravDrainPressView[ei][phaseIndex] =std::abs( gravityCoefficient * ( matrixFluidDensity[ei][0][phaseIndex] - fractureTotalDensity[ei][0] ) * Lz / 2);
    }
  } );
}

void SimpleGravityDrainagePressure::updateState( real64 const gravityCoefficient,
                                                  arrayView1d< real64 const > const & densityMatrix,
                                                  arrayView1d< real64 const > const & densityFracture,
                                                  localIndex const numElements,
                                                  real64 const Lz )
{
  GEOS_ASSERT( m_gravityDrainagePressure.size() >= numElements );
  GEOS_ASSERT( densityMatrix.size() >= numElements );
  GEOS_ASSERT( densityFracture.size() >= numElements );

  // Compute gravity drainage pressure for each element
  // P_grav = (rho_fracture - rho_matrix) * g * Lz / 2
  /*
  real64 const halfLz = 0.5 * Lz;

  for( localIndex k = 0; k < numElements; ++k )
  {
    real64 const densityDifference = densityFracture[k] - densityMatrix[k];
    m_gravityDrainagePressure[k][] = densityDifference * gravityCoefficient * halfLz;
  }
  */
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SimpleGravityDrainagePressure, string const &, Group * const )

}

}
