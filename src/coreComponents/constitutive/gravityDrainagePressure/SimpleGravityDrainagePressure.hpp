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
 * @file SimpleGravityDrainagePressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_SIMPLEGRAVITYDRAINAGEPRESSURE_HPP_
#define GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_SIMPLEGRAVITYDRAINAGEPRESSURE_HPP_

#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Simple gravity drainage pressure model
 *
 * Computes gravity drainage pressure as:
 * P_grav = rho * g * Lz / 2
 * 
 * where rho is the density difference between fracture and matrix,
 * g is gravitational acceleration, and Lz is the fracture spacing in z-direction
 */
class SimpleGravityDrainagePressure : public GravityDrainagePressureBase
{
public:

  SimpleGravityDrainagePressure( string const & name,
                                 dataRepository::Group * const parent );

  static string catalogName() { return "SimpleGravityDrainagePressure"; }

  virtual string getCatalogName() const override { return catalogName(); }


  void setupGravityDrainagePressure(arrayView2d< real64 const> const matrixFluidDensity, arrayView2d< real64 const> const fractureFluidDensity, real64 gravityCoefficient, real64 m_fracSpacingLz ) const override;

  //void setupGravityDrainagePressure(arrayView3d< real64 const> const matrixFluidDensity, arrayView3d< real64 const> const fractureFluidDensity, real64 m_fracSpacingLz ) const override;
  void setupGravityDrainagePressure(arrayView3d< real64 const> const matrixFluidDensity,
                                                                   arrayView2d< real64 const> const matrixPhaseVolumeFraction,
                                                                   arrayView2d< real64 const > const fractureFluidDensity,
                                                                   real64 gravityCoefficient,
                                                                   real64 Lz ) const override;

  void setupGravityDrainagePressureFromPhaseMassDensities(
    arrayView3d< real64 const > const matrixPhaseMassDensity,
    arrayView2d< real64 const > const matrixPhaseVolumeFraction,
    arrayView3d< real64 const > const fracturePhaseMassDensity,
    arrayView2d< real64 const > const fracturePhaseVolumeFraction,
    real64 gravityCoefficient,
    real64 Lz ) const override;

  /**
   * @brief Update gravity drainage pressure
   * @param[in] gravityCoefficient gravity coefficient (g in m/s^2)
   * @param[in] densityMatrix density in matrix elements
   * @param[in] densityFracture density in fracture elements
   * @param[in] numElements number of elements to process
   * @param[in] Lz fracture spacing in z-direction
   */
  virtual void updateState( real64 const gravityCoefficient,
                            arrayView1d< real64 const > const & densityMatrix,
                            arrayView1d< real64 const > const & densityFracture,
                            localIndex const numElements,
                            real64 const Lz ) override;

};

}

}

#endif // GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_SIMPLEGRAVITYDRAINAGEPRESSURE_HPP_
