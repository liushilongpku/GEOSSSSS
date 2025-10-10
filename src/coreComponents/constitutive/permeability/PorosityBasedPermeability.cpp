//
// Created by lsl on 2025/10/9.
//


/**
 * @file PorosityBasedPermeability.cpp
 */
//new tap a tap tap tap



//new tap2
#include "PorosityBasedPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PorosityBasedPermeability::PorosityBasedPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::permeabilityComponentsString(), &m_permeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz components of a diagonal permeability tensor." );
}

std::unique_ptr< ConstitutiveBase >
PorosityBasedPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void PorosityBasedPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  integer const numQuad = 1; // NOTE: enforcing 1 quadrature point

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      m_permeability[ei][q][0] =  m_permeabilityComponents[0];
      m_permeability[ei][q][1] =  m_permeabilityComponents[1];
      m_permeability[ei][q][2] =  m_permeabilityComponents[2];
    }
  }
}

void PorosityBasedPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_permeabilityComponents[0],
                                     m_permeabilityComponents[1],
                                     m_permeabilityComponents[2] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 3; ++dim )
      {
        // The default value is -1 so if it still -1 it needs to be set to something physical
        if( permView[ei][q][dim] < 0 )
        {
          permView[ei][q][dim] =  permComponents[dim];
        }
      }
    }
  } );
}

void PorosityBasedPermeability::postInputInitialization()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorosityBasedPermeability, string const &, Group * const )

}
} /* namespace geos */