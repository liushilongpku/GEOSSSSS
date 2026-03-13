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
 * @file GravityDrainagePressureBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREBASE_HPP_
#define GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Gravity drainage pressure wrapper for updating state
 */
class GravityDrainagePressureUpdate
{
public:

  /**
   * @brief Get gravity drainage pressure field
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 > gravityDrainagePressure() { return m_gravityDrainagePressure; }

  /**
   * @brief Get gravity drainage pressure field (const)
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > gravityDrainagePressure() const { return m_gravityDrainagePressure; }

  /**
   * @brief Get number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_gravityDrainagePressure.size(); }

protected:

  GravityDrainagePressureUpdate( arrayView1d< real64 > const & gravityDrainagePressure )
    : m_gravityDrainagePressure( gravityDrainagePressure )
  {}

  arrayView1d< real64 > m_gravityDrainagePressure;

private:

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       real64 const gravityDrainagePres ) const = 0;
};

/**
 * @brief Base class for gravity drainage pressure models
 */
class GravityDrainagePressureBase : public ConstitutiveBase
{
public:


  GravityDrainagePressureBase( string const & name,
                               dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  /**
   * @brief Get gravity drainage pressure field
   */
  arrayView1d< real64 > gravityDrainagePressure() { return m_gravityDrainagePressure; }

  /**
   * @brief Get gravity drainage pressure field (const)
   */
  arrayView1d< real64 const > gravityDrainagePressure() const { return m_gravityDrainagePressure; }

  /**
   * @brief Update gravity drainage pressure
   */

  virtual void updateState( real64 const gravityCoefficient,
                            arrayView1d< real64 const > const & densityMatrix,
                            arrayView1d< real64 const > const & densityFracture,
                            localIndex const numElements,
                            real64 const Lz ) = 0;

  virtual void setupGravityDrainagePressure() const;

protected:

  array1d< real64 > m_gravityDrainagePressure;
};

}

}

#endif // GEOS_CONSTITUTIVE_GRAVITYDRAINAGEPRESSURE_GRAVITYDRAINAGEPRESSUREBASE_HPP_
