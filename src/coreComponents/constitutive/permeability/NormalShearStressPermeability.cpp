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
 * @file NormalShearStressPermeability.cpp
 */

#include "NormalShearStressPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

NormalShearStressPermeability::NormalShearStressPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::initialPermeabilityComponentsString(), &m_initialPermeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Initial xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::fractureNormalString(), &m_fractureNormal ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fracture-plane normal vector in global coordinates." );

  registerWrapper( viewKeyStruct::normalStressSensitivityString(), &m_normalStressSensitivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal effective stress sensitivity coefficient beta_n [1/Pa]." );

  registerWrapper( viewKeyStruct::shearStressSensitivityString(), &m_shearStressSensitivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear effective stress sensitivity coefficient beta_tau [1/Pa]." );

  registerWrapper( viewKeyStruct::referenceNormalStressString(), &m_referenceNormalStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference compression-positive normal effective stress [Pa]." );

  registerWrapper( viewKeyStruct::referenceShearStressString(), &m_referenceShearStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference effective shear stress magnitude [Pa]." );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.6 ).
    setDescription( "Friction coefficient used in the slip activation function." );

  registerWrapper( viewKeyStruct::cohesionString(), &m_cohesion ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Cohesion used in the slip activation function [Pa]." );

  registerWrapper( viewKeyStruct::slipReferenceStressString(), &m_slipReferenceStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e6 ).
    setDescription( "Reference stress used to normalize positive slip tendency [Pa]." );

  registerWrapper( viewKeyStruct::slipEnhancementCoefficientString(), &m_slipEnhancementCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Slip enhancement coefficient A_s. A value of zero disables slip enhancement." );

  registerWrapper( viewKeyStruct::slipEnhancementExponentString(), &m_slipEnhancementExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Slip enhancement exponent." );

  registerWrapper( viewKeyStruct::decayStartTimeString(), &m_decayStartTime ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Time at which long-term permeability decay starts [s]." );

  registerWrapper( viewKeyStruct::timeDecayRateString(), &m_timeDecayRate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Long-term permeability decay rate lambda_t [1/s]. A value of zero disables time decay." );
}

void NormalShearStressPermeability::postInputInitialization()
{
  GEOS_ERROR_IF( m_normalStressSensitivity < 0.0,
                 getDataContext() << ": normalStressSensitivity must be nonnegative." );
  GEOS_ERROR_IF( m_shearStressSensitivity < 0.0,
                 getDataContext() << ": shearStressSensitivity must be nonnegative." );
  GEOS_ERROR_IF( m_frictionCoefficient < 0.0,
                 getDataContext() << ": frictionCoefficient must be nonnegative." );
  GEOS_ERROR_IF( m_cohesion < 0.0,
                 getDataContext() << ": cohesion must be nonnegative." );
  GEOS_ERROR_IF( m_slipReferenceStress <= 0.0,
                 getDataContext() << ": slipReferenceStress must be positive." );
  GEOS_ERROR_IF( m_slipEnhancementCoefficient < 0.0,
                 getDataContext() << ": slipEnhancementCoefficient must be nonnegative." );
  GEOS_ERROR_IF( m_slipEnhancementExponent < 0.0,
                 getDataContext() << ": slipEnhancementExponent must be nonnegative." );
  GEOS_ERROR_IF( m_timeDecayRate < 0.0,
                 getDataContext() << ": timeDecayRate must be nonnegative." );

  real64 const normSquared = m_fractureNormal[0] * m_fractureNormal[0]
                             + m_fractureNormal[1] * m_fractureNormal[1]
                             + m_fractureNormal[2] * m_fractureNormal[2];
  GEOS_ERROR_IF( normSquared <= 0.0,
                 getDataContext() << ": fractureNormal must have nonzero magnitude." );

  real64 const invNorm = 1.0 / sqrt( normSquared );
  for( integer dim = 0; dim < 3; ++dim )
  {
    m_fractureNormal[dim] *= invNorm;
  }
}

void NormalShearStressPermeability::allocateConstitutiveData( Group & parent,
                                                              localIndex const numPts )
{
  PermeabilityBase::allocateConstitutiveData( parent, numPts );

  GEOS_UNUSED_VAR( numPts );
  integer constexpr numQuad = 1;

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim = 0; dim < 3; ++dim )
      {
        m_permeability[ei][q][dim] = m_initialPermeabilityComponents[dim];
      }
    }
  }
}

void NormalShearStressPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1;

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_initialPermeabilityComponents[0],
                                     m_initialPermeabilityComponents[1],
                                     m_initialPermeabilityComponents[2] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim = 0; dim < 3; ++dim )
      {
        if( permView[ei][q][dim] < 0 )
        {
          permView[ei][q][dim] = permComponents[dim];
        }
      }
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, NormalShearStressPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geos */
