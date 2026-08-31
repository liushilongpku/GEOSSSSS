/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/FluxComputeKernelBase.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/PPUPhaseFlux.hpp"

#include <gtest/gtest.h>

using geos::isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities::PPUPhaseFlux;

TEST( DualContinuumExchangeUpwinding, preservesStandardPPUByDefault )
{
  EXPECT_EQ( PPUPhaseFlux::selectMobilitySupport( 1.0, 0 ), 0 );
  EXPECT_EQ( PPUPhaseFlux::selectMobilitySupport( 0.0, 0 ), 0 );
  EXPECT_EQ( PPUPhaseFlux::selectMobilitySupport( -1.0, 0 ), 1 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::fractureCoverage( -1.0, 0, 0.25 ), 1.0 );
}

TEST( DualContinuumExchangeUpwinding, separatesMobilityFromTransportUpstream )
{
  EXPECT_EQ( PPUPhaseFlux::selectMobilitySupport( 1.0, 1 ), 0 );
  EXPECT_EQ( PPUPhaseFlux::selectMobilitySupport( -1.0, 1 ), 0 );
  EXPECT_EQ( PPUPhaseFlux::selectCompositionUpstream( 1.0 ), 0 );
  EXPECT_EQ( PPUPhaseFlux::selectCompositionUpstream( -1.0 ), 1 );
}

TEST( DualContinuumExchangeUpwinding, appliesFracturePhaseCoverageOnlyToReverseExchange )
{
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::fractureCoverage( 1.0, 1, 0.25 ), 1.0 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::fractureCoverage( -1.0, 1, 0.25 ), 0.25 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::fractureCoverage( -1.0, 1, 0.0 ), 0.0 );
}

TEST( DualContinuumExchangeUpwinding, selectsConfiguredOrCurrentMatrixRelativePermeability )
{
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::reverseExchangeRelPerm( 0.42, 0.01 ), 0.42 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::reverseExchangeRelPerm( -1.0, 0.01 ), 0.01 );
}

TEST( DualContinuumExchangeUpwinding, formsReverseMobilityFromMatrixRelPermAndFracturePVT )
{
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::reverseExchangeMobility( 0.42, 100.0, 0.02, 0.25 ), 525.0 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::reverseExchangeMobility( 0.42, 100.0, 0.02, 0.0 ), 0.0 );
  EXPECT_DOUBLE_EQ( PPUPhaseFlux::reverseExchangeMobility( 0.42, 100.0, 0.0, 0.25 ), 0.0 );
}
