//
// Created by lsl on 2025/10/9.
//

/**
 * @file PorosityBasedPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_POROSITYBASEDPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_POROSITYBASEDPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geos
{
namespace constitutive
{

class PorosityBasedPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

    PorosityBasedPermeabilityUpdate(arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure)
      : PermeabilityBaseUpdate( permeability, dPerm_dPressure )
    {}

private:

};

class PorosityBasedPermeability : public PermeabilityBase
{
public:

  PorosityBasedPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PorosityBasedPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PorosityBasedPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * permeabilityComponentsString() { return "permeabilityComponents"; }
  } viewKeys;

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  R1Tensor m_permeabilityComponents;

};
}/*namespace constitutive*/

}/*namespace geos*/


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_POROSITYBASEDPERMEABILITY_HPP_