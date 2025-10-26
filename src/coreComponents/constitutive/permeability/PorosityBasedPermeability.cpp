//
// Created by lsl on 2025/10/9.
//


/**
 * @file PorosityBasedPermeability.cpp
 */

#include "PorosityBasedPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

//NOTE@LSL: add attribute like this
PorosityBasedPermeability::PorosityBasedPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::particleDiameterString(), &m_particleDiameter ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Diameter of the spherical particles." );

  registerWrapper( viewKeyStruct::sphericityString(), &m_sphericity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Sphericity of the particles." );

  registerWrapper( viewKeyStruct::anisotropyString(), &m_anisotropy ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 1.0, 1.0, 1.0 } ).
    setDescription( "Anisotropy factors for three permeability components." );

  registerWrapper( viewKeyStruct::dPerm_dPorosityString(), &m_dPerm_dPorosity );
}

std::unique_ptr< ConstitutiveBase >
PorosityBasedPermeability::deliverClone( string const & name,
                                        Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void PorosityBasedPermeability::allocateConstitutiveData( Group & parent,
                                                         localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dPorosity.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );
}

void PorosityBasedPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_particleDiameter,
                                     m_particleDiameter,
                                     m_particleDiameter };

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
