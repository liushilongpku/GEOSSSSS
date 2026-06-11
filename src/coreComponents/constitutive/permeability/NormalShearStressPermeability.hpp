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
 * @file NormalShearStressPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_NORMALSHEARSTRESSPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_NORMALSHEARSTRESSPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geos
{
namespace constitutive
{

class NormalShearStressPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  NormalShearStressPermeabilityUpdate( R1Tensor const initialPermeability,
                                       R1Tensor const fractureNormal,
                                       real64 const normalStressSensitivity,
                                       real64 const shearStressSensitivity,
                                       real64 const referenceNormalStress,
                                       real64 const referenceShearStress,
                                       real64 const frictionCoefficient,
                                       real64 const cohesion,
                                       real64 const slipReferenceStress,
                                       real64 const slipEnhancementCoefficient,
                                       real64 const slipEnhancementExponent,
                                       real64 const decayStartTime,
                                       real64 const timeDecayRate,
                                       arrayView3d< real64 > const & permeability,
                                       arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_initialPermeability( initialPermeability ),
    m_fractureNormal( fractureNormal ),
    m_normalStressSensitivity( normalStressSensitivity ),
    m_shearStressSensitivity( shearStressSensitivity ),
    m_referenceNormalStress( referenceNormalStress ),
    m_referenceShearStress( referenceShearStress ),
    m_frictionCoefficient( frictionCoefficient ),
    m_cohesion( cohesion ),
    m_slipReferenceStress( slipReferenceStress ),
    m_slipEnhancementCoefficient( slipEnhancementCoefficient ),
    m_slipEnhancementExponent( slipEnhancementExponent ),
    m_decayStartTime( decayStartTime ),
    m_timeDecayRate( timeDecayRate )
  {}

  GEOS_HOST_DEVICE
  virtual void updateFromEffectiveStress( localIndex const k,
                                          localIndex const q,
                                          real64 const ( & effectiveStress )[6] ) const override final
  {
    updateFromEffectiveStress( k, q, effectiveStress, 0.0 );
  }

  GEOS_HOST_DEVICE
  virtual void updateFromEffectiveStress( localIndex const k,
                                          localIndex const q,
                                          real64 const ( & effectiveStress )[6],
                                          real64 const & currentTime ) const override final
  {
    real64 const traction[3] =
    {
      effectiveStress[0] * m_fractureNormal[0]
      + effectiveStress[5] * m_fractureNormal[1]
      + effectiveStress[4] * m_fractureNormal[2],
      effectiveStress[5] * m_fractureNormal[0]
      + effectiveStress[1] * m_fractureNormal[1]
      + effectiveStress[3] * m_fractureNormal[2],
      effectiveStress[4] * m_fractureNormal[0]
      + effectiveStress[3] * m_fractureNormal[1]
      + effectiveStress[2] * m_fractureNormal[2]
    };
    real64 const normalTraction = traction[0] * m_fractureNormal[0]
                                  + traction[1] * m_fractureNormal[1]
                                  + traction[2] * m_fractureNormal[2];
    real64 const tractionNormSquared = traction[0] * traction[0]
                                       + traction[1] * traction[1]
                                       + traction[2] * traction[2];
    real64 const shearStressSquared = tractionNormSquared - normalTraction * normalTraction;
    real64 const shearStress = sqrt( shearStressSquared > 0.0 ? shearStressSquared : 0.0 );
    real64 const normalStress = -normalTraction;

    updateFromNormalAndShearStressAndTime( k, q, normalStress, shearStress, currentTime );
  }

  GEOS_HOST_DEVICE
  virtual void updateFromNormalAndShearStress( localIndex const k,
                                               localIndex const q,
                                               real64 const & normalStress,
                                               real64 const & shearStress ) const override final
  {
    updateFromNormalAndShearStressAndTime( k, q, normalStress, shearStress, 0.0 );
  }

  GEOS_HOST_DEVICE
  void updateFromNormalAndShearStressAndTime( localIndex const k,
                                              localIndex const q,
                                              real64 const & normalStress,
                                              real64 const & shearStress,
                                              real64 const & currentTime ) const
  {
    real64 const normalMultiplier =
      exp( -m_normalStressSensitivity * ( normalStress - m_referenceNormalStress ) );
    real64 const shearMultiplier =
      exp( m_shearStressSensitivity * ( shearStress - m_referenceShearStress ) );

    real64 slipMultiplier = 1.0;
    real64 const slipFunction = shearStress - m_frictionCoefficient * normalStress - m_cohesion;
    if( slipFunction > 0.0 && m_slipEnhancementCoefficient > 0.0 && m_slipReferenceStress > 0.0 )
    {
      slipMultiplier +=
        m_slipEnhancementCoefficient * pow( slipFunction / m_slipReferenceStress, m_slipEnhancementExponent );
    }

    real64 const timeMultiplier =
      currentTime <= m_decayStartTime ? 1.0 : exp( -m_timeDecayRate * ( currentTime - m_decayStartTime ) );

    real64 const multiplier = normalMultiplier * shearMultiplier * slipMultiplier * timeMultiplier;

    for( localIndex dim = 0; dim < m_permeability.size( 2 ); ++dim )
    {
      m_permeability[k][q][dim] = m_initialPermeability[dim] * multiplier;
    }
  }

private:

  /// Permeability components at the reference stress state.
  R1Tensor m_initialPermeability;

  /// Unit normal of the fracture plane in global coordinates.
  R1Tensor m_fractureNormal;

  /// Normal stress sensitivity coefficient beta_n [1/Pa].
  real64 m_normalStressSensitivity;

  /// Shear stress sensitivity coefficient beta_tau [1/Pa].
  real64 m_shearStressSensitivity;

  /// Reference compression-positive normal effective stress [Pa].
  real64 m_referenceNormalStress;

  /// Reference effective shear stress magnitude [Pa].
  real64 m_referenceShearStress;

  /// Friction coefficient for the slip activation function.
  real64 m_frictionCoefficient;

  /// Cohesion for the slip activation function [Pa].
  real64 m_cohesion;

  /// Reference stress used to normalize positive slip tendency [Pa].
  real64 m_slipReferenceStress;

  /// Slip enhancement coefficient A_s [-].
  real64 m_slipEnhancementCoefficient;

  /// Slip enhancement exponent m [-].
  real64 m_slipEnhancementExponent;

  /// Time at which long-term permeability decay starts [s].
  real64 m_decayStartTime;

  /// Long-term permeability decay rate lambda_t [1/s].
  real64 m_timeDecayRate;
};

class NormalShearStressPermeability : public PermeabilityBase
{
public:

  NormalShearStressPermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "NormalShearStressPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  virtual void initializeState() const override final;

  using KernelWrapper = NormalShearStressPermeabilityUpdate;

  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_initialPermeabilityComponents,
                          m_fractureNormal,
                          m_normalStressSensitivity,
                          m_shearStressSensitivity,
                          m_referenceNormalStress,
                          m_referenceShearStress,
                          m_frictionCoefficient,
                          m_cohesion,
                          m_slipReferenceStress,
                          m_slipEnhancementCoefficient,
                          m_slipEnhancementExponent,
                          m_decayStartTime,
                          m_timeDecayRate,
                          m_permeability,
                          m_dPerm_dPressure );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * initialPermeabilityComponentsString()
    { return "initialPermeabilityComponents"; }
    static constexpr char const * fractureNormalString()
    { return "fractureNormal"; }
    static constexpr char const * normalStressSensitivityString()
    { return "normalStressSensitivity"; }
    static constexpr char const * shearStressSensitivityString()
    { return "shearStressSensitivity"; }
    static constexpr char const * referenceNormalStressString()
    { return "referenceNormalStress"; }
    static constexpr char const * referenceShearStressString()
    { return "referenceShearStress"; }
    static constexpr char const * frictionCoefficientString()
    { return "frictionCoefficient"; }
    static constexpr char const * cohesionString()
    { return "cohesion"; }
    static constexpr char const * slipReferenceStressString()
    { return "slipReferenceStress"; }
    static constexpr char const * slipEnhancementCoefficientString()
    { return "slipEnhancementCoefficient"; }
    static constexpr char const * slipEnhancementExponentString()
    { return "slipEnhancementExponent"; }
    static constexpr char const * decayStartTimeString()
    { return "decayStartTime"; }
    static constexpr char const * timeDecayRateString()
    { return "timeDecayRate"; }
  };

  R1Tensor const & fractureNormal() const
  { return m_fractureNormal; }

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the reference stress state.
  R1Tensor m_initialPermeabilityComponents;

  /// Unit normal of the fracture plane in global coordinates.
  R1Tensor m_fractureNormal;

  /// Normal stress sensitivity coefficient beta_n [1/Pa].
  real64 m_normalStressSensitivity;

  /// Shear stress sensitivity coefficient beta_tau [1/Pa].
  real64 m_shearStressSensitivity;

  /// Reference compression-positive normal effective stress [Pa].
  real64 m_referenceNormalStress;

  /// Reference effective shear stress magnitude [Pa].
  real64 m_referenceShearStress;

  /// Friction coefficient for the slip activation function.
  real64 m_frictionCoefficient;

  /// Cohesion for the slip activation function [Pa].
  real64 m_cohesion;

  /// Reference stress used to normalize positive slip tendency [Pa].
  real64 m_slipReferenceStress;

  /// Slip enhancement coefficient A_s [-].
  real64 m_slipEnhancementCoefficient;

  /// Slip enhancement exponent m [-].
  real64 m_slipEnhancementExponent;

  /// Time at which long-term permeability decay starts [s].
  real64 m_decayStartTime;

  /// Long-term permeability decay rate lambda_t [1/s].
  real64 m_timeDecayRate;
};

} /* namespace constitutive */
} /* namespace geos */

#endif // GEOS_CONSTITUTIVE_PERMEABILITY_NORMALSHEARSTRESSPERMEABILITY_HPP_
