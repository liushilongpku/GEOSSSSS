/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "constitutive/gravityDrainagePressure/ThomasGasOilGravityDrainagePressure.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::constitutive;
using namespace geos::dataRepository;

namespace
{

class ThomasGasOilGravityDrainagePressureTest : public ::testing::Test
{
protected:
  ThomasGasOilGravityDrainagePressureTest()
    : m_node(),
    m_parent( "parent", m_node ),
    m_model( "gravityDrainagePressure", &m_parent )
  {
    m_parent.resize( 1 );
    m_model.allocateConstitutiveData( m_parent, 1 );
  }

  void check( string_array const & phaseNames, real64 const gravityCoefficient )
  {
    localIndex const numPhases = LvArray::integerConversion< localIndex >( phaseNames.size() );
    array3d< real64 > matrixDensity( 1, 1, numPhases );
    array3d< real64 > fractureDensity( 1, 1, numPhases );
    array2d< real64 > matrixSaturation( 1, numPhases );
    array2d< real64 > fractureSaturation( 1, numPhases );
    integer oilIndex = -1;
    integer gasIndex = -1;

    for( integer ip = 0; ip < numPhases; ++ip )
    {
      matrixDensity[0][0][ip] = phaseNames[ip] == "oil" ? 800.0 : phaseNames[ip] == "gas" ? 100.0 : 1000.0;
      fractureDensity[0][0][ip] = 50.0;
      matrixSaturation[0][ip] = 1.0 / numPhases;
      fractureSaturation[0][ip] = 1.0 / numPhases;
      oilIndex = phaseNames[ip] == "oil" ? ip : oilIndex;
      gasIndex = phaseNames[ip] == "gas" ? ip : gasIndex;
    }

    m_model.setupGravityDrainagePressureFromPhaseMassDensities(
      matrixDensity.toViewConst(), matrixSaturation.toViewConst(),
      fractureDensity.toViewConst(), fractureSaturation.toViewConst(),
      phaseNames, gravityCoefficient, 4.0 );

    arrayView3d< real64 const > const gravityDrainagePressure = m_model.gravityDrainagePressure();
    EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][oilIndex], 14000.0 );
    EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][gasIndex], 0.0 );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      if( ip != oilIndex )
      {
        EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][ip], 0.0 );
      }
    }
  }

  conduit::Node m_node;
  Group m_parent;
  ThomasGasOilGravityDrainagePressure m_model;
};

TEST_F( ThomasGasOilGravityDrainagePressureTest, computesOilGasBalance )
{
  check( { "oil", "gas", "water" }, -10.0 );
}

TEST_F( ThomasGasOilGravityDrainagePressureTest, isIndependentOfGravityAndPhaseOrdering )
{
  check( { "water", "gas", "oil" }, 10.0 );
}

TEST_F( ThomasGasOilGravityDrainagePressureTest, rejectsMissingGasPhase )
{
  string_array const phaseNames = { "oil", "water" };
  array3d< real64 > density( 1, 1, 2 );
  array2d< real64 > saturation( 1, 2 );
  EXPECT_THROW(
    m_model.setupGravityDrainagePressureFromPhaseMassDensities(
      density.toViewConst(), saturation.toViewConst(), density.toViewConst(), saturation.toViewConst(),
      phaseNames, -9.81, 3.05 ),
    InputError );
}

} // namespace
