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

#include "constitutive/gravityDrainagePressure/CompositionalPerPhaseGravityDrainagePressure.hpp"
#include "constitutive/gravityDrainagePressure/SimpleGravityDrainagePressure.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::constitutive;
using namespace geos::dataRepository;

namespace
{

class GravityDrainagePressureTest : public ::testing::Test
{
protected:
  GravityDrainagePressureTest()
    : m_node(),
    m_parent( "parent", m_node ),
    m_phaseModel( "phaseModel", &m_parent ),
    m_simpleModel( "simpleModel", &m_parent )
  {
    m_parent.resize( 1 );
    m_phaseModel.allocateConstitutiveData( m_parent, 1 );
    m_simpleModel.allocateConstitutiveData( m_parent, 1 );
  }

  conduit::Node m_node;
  Group m_parent;
  CompositionalPerPhaseGravityDrainagePressure m_phaseModel;
  SimpleGravityDrainagePressure m_simpleModel;
};

TEST_F( GravityDrainagePressureTest, computesOneHeadPerPhaseFromMixtureDirection )
{
  string_array const phaseNames = { "oil", "gas", "water" };
  array3d< real64 > matrixDensity( 1, 1, 3 );
  array3d< real64 > fractureDensity( 1, 1, 3 );
  array2d< real64 > matrixSaturation( 1, 3 );
  array2d< real64 > fractureSaturation( 1, 3 );
  for( localIndex ip = 0; ip < 3; ++ip )
  {
    matrixDensity[0][0][ip] = ip == 0 ? 800.0 : ip == 1 ? 100.0 : 1000.0;
    fractureDensity[0][0][ip] = 50.0;
    matrixSaturation[0][ip] = 1.0 / 3.0;
    fractureSaturation[0][ip] = 1.0 / 3.0;
  }

  m_phaseModel.setupGravityDrainagePressureFromPhaseMassDensities(
    matrixDensity.toViewConst(), matrixSaturation.toViewConst(),
    fractureDensity.toViewConst(), fractureSaturation.toViewConst(),
    phaseNames, -10.0, 4.0 );

  arrayView3d< real64 const > const gravityDrainagePressure = m_phaseModel.gravityDrainagePressure();
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][0], 16000.0 );
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][1], 2000.0 );
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][2], 20000.0 );
}

TEST_F( GravityDrainagePressureTest, simpleModelUsesOneScalarForEveryPhase )
{
  array3d< real64 > matrixDensity( 1, 1, 3 );
  array3d< real64 > fractureDensity( 1, 1, 3 );
  array2d< real64 > matrixSaturation( 1, 3 );
  array2d< real64 > fractureSaturation( 1, 3 );
  for( localIndex ip = 0; ip < 3; ++ip )
  {
    matrixDensity[0][0][ip] = ip == 0 ? 800.0 : ip == 1 ? 100.0 : 1000.0;
    fractureDensity[0][0][ip] = 50.0;
    matrixSaturation[0][ip] = 1.0 / 3.0;
    fractureSaturation[0][ip] = 1.0 / 3.0;
  }

  m_simpleModel.setupGravityDrainagePressureFromPhaseMassDensities(
    matrixDensity.toViewConst(), matrixSaturation.toViewConst(),
    fractureDensity.toViewConst(), fractureSaturation.toViewConst(),
    { "oil", "gas", "water" }, -10.0, 4.0 );

  arrayView3d< real64 const > const gravityDrainagePressure = m_simpleModel.gravityDrainagePressure();
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][0], gravityDrainagePressure[0][0][1] );
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][1], gravityDrainagePressure[0][0][2] );
}

TEST_F( GravityDrainagePressureTest, reversesAllPhaseHeadsForNegativeMixtureContrast )
{
  array3d< real64 > matrixDensity( 1, 1, 3 );
  array3d< real64 > fractureDensity( 1, 1, 3 );
  array2d< real64 > matrixSaturation( 1, 3 );
  array2d< real64 > fractureSaturation( 1, 3 );
  for( localIndex ip = 0; ip < 3; ++ip )
  {
    matrixDensity[0][0][ip] = 50.0;
    fractureDensity[0][0][ip] = ip == 0 ? 800.0 : ip == 1 ? 100.0 : 1000.0;
    matrixSaturation[0][ip] = 1.0 / 3.0;
    fractureSaturation[0][ip] = 1.0 / 3.0;
  }

  m_phaseModel.setupGravityDrainagePressureFromPhaseMassDensities(
    matrixDensity.toViewConst(), matrixSaturation.toViewConst(),
    fractureDensity.toViewConst(), fractureSaturation.toViewConst(),
    { "oil", "gas", "water" }, -10.0, 4.0 );

  arrayView3d< real64 const > const gravityDrainagePressure = m_phaseModel.gravityDrainagePressure();
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][0], -1000.0 );
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][1], -1000.0 );
  EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][2], -1000.0 );
}

TEST_F( GravityDrainagePressureTest, returnsZeroHeadsForEqualMixtureDensities )
{
  array3d< real64 > matrixDensity( 1, 1, 3 );
  array3d< real64 > fractureDensity( 1, 1, 3 );
  array2d< real64 > matrixSaturation( 1, 3 );
  array2d< real64 > fractureSaturation( 1, 3 );
  for( localIndex ip = 0; ip < 3; ++ip )
  {
    matrixDensity[0][0][ip] = fractureDensity[0][0][ip] = 50.0 + 25.0 * ip;
    matrixSaturation[0][ip] = fractureSaturation[0][ip] = 1.0 / 3.0;
  }

  m_phaseModel.setupGravityDrainagePressureFromPhaseMassDensities(
    matrixDensity.toViewConst(), matrixSaturation.toViewConst(),
    fractureDensity.toViewConst(), fractureSaturation.toViewConst(),
    { "oil", "gas", "water" }, -10.0, 4.0 );

  arrayView3d< real64 const > const gravityDrainagePressure = m_phaseModel.gravityDrainagePressure();
  for( localIndex ip = 0; ip < 3; ++ip )
  {
    EXPECT_DOUBLE_EQ( gravityDrainagePressure[0][0][ip], 0.0 );
  }
}

} // namespace
