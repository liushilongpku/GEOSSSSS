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
 * @file FlowSolverBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

DECLARE_FIELD( pressure,
               "pressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure" );

DECLARE_FIELD( pressure_n,
               "pressure_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Pressure at the previous converged time step" );

DECLARE_FIELD( pressure_k,
               "pressure_k",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Pressure at the previous sequential iteration" );

DECLARE_FIELD( initialPressure,
               "initialPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial pressure" );

DECLARE_FIELD( deltaPressure,
               "deltaPressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Delta pressure: current pressure - initial pressure" );

DECLARE_FIELD( facePressure,
               "facePressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face pressure" );

DECLARE_FIELD( facePressure_n,
               "facePressure_n",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Face pressure at the previous converged time step" );

DECLARE_FIELD( pressureGradient,
               "pressureGradient",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure gradient" );

DECLARE_FIELD( temperature,
               "temperature",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Temperature" );

DECLARE_FIELD( temperature_n,
               "temperature_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Temperature at the previous converged time step" );

DECLARE_FIELD( temperature_k,
               "temperature_k",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Temperature at the previous sequential iteration" );

DECLARE_FIELD( initialTemperature,
               "initialTemperature",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial temperature" );

DECLARE_FIELD( faceTemperature,
               "faceTemperature",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face temperature" );

DECLARE_FIELD( netToGross,
               "netToGross",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Net to gross" );

DECLARE_FIELD( deltaVolume,
               "deltaVolume",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Accumulated volume updates" );

DECLARE_FIELD( aperture0,
               "aperture_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial aperture" );

DECLARE_FIELD( hydraulicAperture,
               "hydraulicAperture",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Hydraulic aperture" );

DECLARE_FIELD( gravityCoefficient,
               "gravityCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Gravity coefficient (dot product of gravity acceleration by gravity vector)" );

DECLARE_FIELD( mimGravityCoefficient,
               "mimGravityCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Mimetic gravity coefficient" );

DECLARE_FIELD( macroElementIndex,
               "macroElementIndex",
               array1d< integer >,
               -1,
               LEVEL_1,
               WRITE_AND_READ,
               "Index of the macroelement for a given element" );

DECLARE_FIELD( bcPressure,
               "bcPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition pressure" );

DECLARE_FIELD( bcTemperature,
               "bcTemperature",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition temperature" );

DECLARE_FIELD( elementStabConstant,
               "elementStabConstant",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "Stabilization constant for pressure jump stabilization" );

DECLARE_FIELD( transMultiplier,
               "permeabilityTransMultiplier",
               array1d< real64 >,
               1,
               LEVEL_0,
               WRITE_AND_READ,
               "Permeability transmissibility multipliers" );

DECLARE_FIELD( pressureScalingFactor,
               "pressureScalingFactor",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for pressure" );

DECLARE_FIELD( temperatureScalingFactor,
               "temperatureScalingFactor",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for temperature" );

DECLARE_FIELD( energy,
               "energy",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Energy" );

DECLARE_FIELD( energy_n,
               "energy_n",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Energy at the previous converged time step" );

DECLARE_FIELD( fractureCreationTime,
               "fractureCreationTime",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "The creation time for the fracture cell." );

DECLARE_FIELD( pressure_DPDK,
               "pressure_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure (DPDK)" );

DECLARE_FIELD( pressure_n_DPDK,
               "pressure_n_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Pressure at the previous converged time step (DPDK)" );

DECLARE_FIELD( pressure_k_DPDK,
               "pressure_k_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Pressure at the previous sequential iteration (DPDK)" );

DECLARE_FIELD( initialPressure_DPDK,
               "initialPressure_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial pressure (DPDK)" );

DECLARE_FIELD( deltaPressure_DPDK,
               "deltaPressure_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Delta pressure: current pressure - initial pressure (DPDK)" );

DECLARE_FIELD( facePressure_DPDK,
               "facePressure_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face pressure (DPDK)" );

DECLARE_FIELD( facePressure_n_DPDK,
               "facePressure_n_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Face pressure at the previous converged time step (DPDK)" );

DECLARE_FIELD( pressureGradient_DPDK,
               "pressureGradient_DPDK",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure gradient (DPDK)" );

DECLARE_FIELD( temperature_DPDK,
               "temperature_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Temperature (DPDK)" );

DECLARE_FIELD( temperature_n_DPDK,
               "temperature_n_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Temperature at the previous converged time step (DPDK)" );

DECLARE_FIELD( temperature_k_DPDK,
               "temperature_k_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Temperature at the previous sequential iteration (DPDK)" );

DECLARE_FIELD( initialTemperature_DPDK,
               "initialTemperature_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial temperature (DPDK)" );

DECLARE_FIELD( faceTemperature_DPDK,
               "faceTemperature_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face temperature (DPDK)" );

DECLARE_FIELD( netToGross_DPDK,
               "netToGross_DPDK",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Net to gross (DPDK)" );

DECLARE_FIELD( deltaVolume_DPDK,
               "deltaVolume_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Accumulated volume updates (DPDK)" );

DECLARE_FIELD( aperture0_DPDK,
               "aperture_n_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial aperture (DPDK)" );

DECLARE_FIELD( hydraulicAperture_DPDK,
               "hydraulicAperture_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Hydraulic aperture (DPDK)" );

DECLARE_FIELD( gravityCoefficient_DPDK,
               "gravityCoefficient_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Gravity coefficient (dot product of gravity acceleration by gravity vector) (DPDK)" );

DECLARE_FIELD( mimGravityCoefficient_DPDK,
               "mimGravityCoefficient_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Mimetic gravity coefficient (DPDK)" );

DECLARE_FIELD( macroElementIndex_DPDK,
               "macroElementIndex_DPDK",
               array1d< integer >,
               -1,
               LEVEL_1,
               WRITE_AND_READ,
               "Index of the macroelement for a given element (DPDK)" );

DECLARE_FIELD( bcPressure_DPDK,
               "bcPressure_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition pressure (DPDK)" );

DECLARE_FIELD( bcTemperature_DPDK,
               "bcTemperature_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition temperature (DPDK)" );

DECLARE_FIELD( elementStabConstant_DPDK,
               "elementStabConstant_DPDK",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "Stabilization constant for pressure jump stabilization (DPDK)" );

DECLARE_FIELD( transMultiplier_DPDK,
               "permeabilityTransMultiplier_DPDK",
               array1d< real64 >,
               1,
               LEVEL_0,
               WRITE_AND_READ,
               "Permeability transmissibility multipliers (DPDK)" );

DECLARE_FIELD( pressureScalingFactor_DPDK,
               "pressureScalingFactor_DPDK",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for pressure (DPDK)" );

DECLARE_FIELD( temperatureScalingFactor_DPDK,
               "temperatureScalingFactor_DPDK",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for temperature (DPDK)" );

DECLARE_FIELD( energy_DPDK,
               "energy_DPDK",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Energy (DPDK)" );

DECLARE_FIELD( energy_n_DPDK,
               "energy_n_DPDK",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Energy at the previous converged time step (DPDK)" );

DECLARE_FIELD( fractureCreationTime_DPDK,
               "fractureCreationTime_DPDK",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "The creation time for the fracture cell (DPDK)." );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_
