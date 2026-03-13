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
 * @file GravityDrainagePressureFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREFIELDS_HPP_
#define GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace gravdrainage
{

DECLARE_FIELD( gravityDrainagePressure,
               "gravityDrainagePressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Gravity drainage pressure for dual continuum model" );

}

}

}

#endif // GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREFIELDS_HPP_
