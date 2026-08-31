/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */

#include "RelpermDriverRunTest.hpp"
#include "constitutive/relativePermeability/ThreePhaseTableRelativePermeability.hpp"

namespace geos
{
template void RelpermDriver::runTest< geos::constitutive::ThreePhaseTableRelativePermeability >(
  geos::constitutive::ThreePhaseTableRelativePermeability &, arrayView2d< real64 > const & );
}
