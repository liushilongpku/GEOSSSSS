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
 * @file MeanEffectiveStressPermeability.cpp
 */

#include "MeanEffectiveStressPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

MeanEffectiveStressPermeability::MeanEffectiveStressPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::initialPermeabilityComponentsString(), &m_initialPermeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Initial xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::stressSensitivityString(), &m_stressSensitivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Stress sensitivity coefficient beta_sigma [1/Pa]." );

  registerWrapper( viewKeyStruct::referenceMeanEffectiveStressString(), &m_referenceMeanEffectiveStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference compression-positive mean effective stress [Pa]." );
}

void MeanEffectiveStressPermeability::postInputInitialization()
{
  GEOS_ERROR_IF( m_stressSensitivity < 0.0,
                 getDataContext() << ": stressSensitivity must be nonnegative." );
}

void MeanEffectiveStressPermeability::allocateConstitutiveData( Group & parent,
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

void MeanEffectiveStressPermeability::initializeState() const
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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MeanEffectiveStressPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geos */
