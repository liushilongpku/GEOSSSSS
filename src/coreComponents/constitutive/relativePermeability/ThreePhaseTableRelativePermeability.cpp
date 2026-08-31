/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ThreePhaseTableRelativePermeability.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ThreePhaseTableRelativePermeability::ThreePhaseTableRelativePermeability( string const & name, Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::waterRelPermTableNameString(), &m_waterRelPermTableName ).setInputFlag( InputFlags::REQUIRED );
  registerWrapper( viewKeyStruct::gasRelPermTableNameString(), &m_gasRelPermTableName ).setInputFlag( InputFlags::REQUIRED );
  registerWrapper( viewKeyStruct::oilRelPermTableNameString(), &m_oilRelPermTableName ).setInputFlag( InputFlags::REQUIRED );
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction )
    .setInputFlag( InputFlags::FALSE ).setSizedFromParent( 0 );
}

void ThreePhaseTableRelativePermeability::postInputInitialization()
{
  RelativePermeabilityBase::postInputInitialization();
  GEOS_THROW_IF( numFluidPhases() != 3, GEOS_FMT( "{}: exactly three phases are required", getFullName() ), InputError );
  FunctionManager const & functions = FunctionManager::getInstance();
  for( string const & name : {m_waterRelPermTableName, m_gasRelPermTableName, m_oilRelPermTableName} )
  {
    GEOS_THROW_IF( !functions.hasGroup( name ), GEOS_FMT( "{}: table {} was not found", getFullName(), name ), InputError );
  }
}

void ThreePhaseTableRelativePermeability::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();
  m_phaseMinVolumeFraction.resize( MAX_NUM_PHASES );
  m_phaseMinVolumeFraction.zero();
}

ThreePhaseTableRelativePermeability::KernelWrapper ThreePhaseTableRelativePermeability::createKernelWrapper()
{
  FunctionManager const & functions = FunctionManager::getInstance();
  m_water = functions.getGroup< TableFunction >( m_waterRelPermTableName ).createKernelWrapper();
  m_gas = functions.getGroup< TableFunction >( m_gasRelPermTableName ).createKernelWrapper();
  m_oil = functions.getGroup< TableFunction >( m_oilRelPermTableName ).createKernelWrapper();
  return KernelWrapper( m_water, m_gas, m_oil, m_phaseTypes, m_phaseOrder,
                        m_phaseRelPerm, m_dPhaseRelPerm_dPhaseVolFrac, m_phaseTrappedVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThreePhaseTableRelativePermeability, string const &, Group * const )

} // namespace constitutive
} // namespace geos
