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
 * @file GravityDrainagePressureBase.cpp
 */

#include "GravityDrainagePressureBase.hpp"
#include "constitutive/gravityDrainagePressure/GravityDrainagePressureFields.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

namespace constitutive
{

GravityDrainagePressureBase::GravityDrainagePressureBase( string const & name,
                                                          dataRepository::Group * const parent )
  : ConstitutiveBase( name, parent )
{
  // Register the gravity drainage pressure field
  registerField< fields::gravdrainage::gravityDrainagePressure >( &m_gravityDrainagePressure );
}

void GravityDrainagePressureBase::setupGravityDrainagePressure(arrayView2d< real64 const > const matrixFluidDensity, arrayView2d< real64 const > const fractureFluidDensity, real64 m_fracSpacingLz ) const
{
  // Default implementation: do nothing, let derived classes override if needed
}

void GravityDrainagePressureBase::setupGravityDrainagePressure(arrayView3d< real64 const > const matrixFluidDensity, arrayView3d< real64 const > const fractureFluidDensity, real64 m_fracSpacingLz ) const
{
  // Default implementation: do nothing, let derived classes override if needed
}

void GravityDrainagePressureBase::setupGravityDrainagePressure(arrayView3d< real64 const> const matrixFluidDensity,
                                   arrayView2d< real64 const> const matrixPhaseVolumeFraction,
                                   arrayView2d< real64 const > const fractureFluidDensity,
                                   real64 gravityCoefficient,
                                   real64 Lz ) const
{
  // Default implementation: do nothing, let derived classes override if needed
}

void GravityDrainagePressureBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                            localIndex const numPts )
{
  int phaseNumber = 3 ;//TODO@LSL 直接分配三个相的数量，不管实际问题是什么
  m_gravityDrainagePressure.resize( numPts, phaseNumber);
  m_gravityDrainagePressure.zero();

  //需要放在设置值的后边，才能分配到constitutive中的subgroup中。
  ConstitutiveBase::allocateConstitutiveData( parent, numPts );

}

}

}
