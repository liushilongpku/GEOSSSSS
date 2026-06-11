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
 * @file MeanEffectiveStressPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_MEANEFFECTIVESTRESSPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_MEANEFFECTIVESTRESSPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geos
{
namespace constitutive
{

class MeanEffectiveStressPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  MeanEffectiveStressPermeabilityUpdate( R1Tensor const initialPermeability,
                                         real64 const stressSensitivity,
                                         real64 const referenceMeanEffectiveStress,
                                         arrayView3d< real64 > const & permeability,
                                         arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_initialPermeability( initialPermeability ),
    m_stressSensitivity( stressSensitivity ),
    m_referenceMeanEffectiveStress( referenceMeanEffectiveStress )
  {}

  GEOS_HOST_DEVICE
  virtual void updateFromMeanEffectiveStress( localIndex const k,
                                              localIndex const q,
                                              real64 const & meanEffectiveStress ) const override final
  {
    real64 const multiplier =
      exp( -m_stressSensitivity * ( meanEffectiveStress - m_referenceMeanEffectiveStress ) );

    for( localIndex dim = 0; dim < m_permeability.size( 2 ); ++dim )
    {
      m_permeability[k][q][dim] = m_initialPermeability[dim] * multiplier;
    }
  }

private:

  /// Permeability components at the reference mean effective stress.
  R1Tensor m_initialPermeability;

  /// Stress sensitivity coefficient beta_sigma [1/Pa].
  real64 m_stressSensitivity;

  /// Reference compression-positive mean effective stress [Pa].
  real64 m_referenceMeanEffectiveStress;
};

class MeanEffectiveStressPermeability : public PermeabilityBase
{
public:

  MeanEffectiveStressPermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "MeanEffectiveStressPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  virtual void initializeState() const override final;

  using KernelWrapper = MeanEffectiveStressPermeabilityUpdate;

  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_initialPermeabilityComponents,
                          m_stressSensitivity,
                          m_referenceMeanEffectiveStress,
                          m_permeability,
                          m_dPerm_dPressure );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * initialPermeabilityComponentsString()
    { return "initialPermeabilityComponents"; }
    static constexpr char const * stressSensitivityString()
    { return "stressSensitivity"; }
    static constexpr char const * referenceMeanEffectiveStressString()
    { return "referenceMeanEffectiveStress"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the reference mean effective stress.
  R1Tensor m_initialPermeabilityComponents;

  /// Stress sensitivity coefficient beta_sigma [1/Pa].
  real64 m_stressSensitivity;

  /// Reference compression-positive mean effective stress [Pa].
  real64 m_referenceMeanEffectiveStress;
};

} /* namespace constitutive */
} /* namespace geos */

#endif // GEOS_CONSTITUTIVE_PERMEABILITY_MEANEFFECTIVESTRESSPERMEABILITY_HPP_
