/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_CONSTITUTIVE_THREEPHASETABLERELATIVEPERMEABILITY_HPP
#define GEOS_CONSTITUTIVE_THREEPHASETABLERELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

class ThreePhaseTableRelativePermeability : public RelativePermeabilityBase
{
public:
  ThreePhaseTableRelativePermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "ThreePhaseTableRelativePermeability"; }
  string getCatalogName() const override { return catalogName(); }

  class KernelWrapper final : public RelativePermeabilityBaseUpdate
  {
  public:
    KernelWrapper( TableFunction::KernelWrapper const & water,
                   TableFunction::KernelWrapper const & gas,
                   TableFunction::KernelWrapper const & oil,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                   arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrappedVolFrac )
      : RelativePermeabilityBaseUpdate( phaseTypes, phaseOrder, phaseRelPerm,
                                        dPhaseRelPerm_dPhaseVolFrac, phaseTrappedVolFrac ),
      m_water( water ), m_gas( gas ), m_oil( oil ) {}

    GEOS_HOST_DEVICE
    void update( localIndex const k, localIndex const q,
                 arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
    {
      arraySlice1d< real64, relperm::USD_RELPERM - 2 > const relPerm = m_phaseRelPerm[k][q];
      arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const derivatives = m_dPhaseRelPerm_dPhaseVolFrac[k][q];
      LvArray::forValuesInSlice( derivatives, []( real64 & value ) { value = 0.0; } );

      integer const ipWater = m_phaseOrder[PhaseType::WATER];
      integer const ipGas = m_phaseOrder[PhaseType::GAS];
      integer const ipOil = m_phaseOrder[PhaseType::OIL];
      real64 derivative[1] = {0.0};
      relPerm[ipWater] = m_water.compute( &(phaseVolFraction)[ipWater], derivative );
      derivatives[ipWater][ipWater] = derivative[0];
      relPerm[ipGas] = m_gas.compute( &(phaseVolFraction)[ipGas], derivative );
      derivatives[ipGas][ipGas] = derivative[0];
      relPerm[ipOil] = m_oil.compute( &(phaseVolFraction)[ipGas], derivative );
      derivatives[ipOil][ipGas] = derivative[0];

      for( integer ip = 0; ip < numPhases(); ++ip )
      {
        m_phaseTrappedVolFrac[k][q][ip] = 0.0;
      }
    }

  private:
    TableFunction::KernelWrapper m_water;
    TableFunction::KernelWrapper m_gas;
    TableFunction::KernelWrapper m_oil;
  };

  using KernelWrapperType = KernelWrapper;
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * waterRelPermTableNameString() { return "waterRelPermTableName"; }
    static constexpr char const * gasRelPermTableNameString() { return "gasRelPermTableName"; }
    static constexpr char const * oilRelPermTableNameString() { return "oilRelPermTableName"; }
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
  };

  arrayView1d< real64 const > getPhaseMinVolumeFraction() const override { return m_phaseMinVolumeFraction; }
  real64 getWettingPhaseMinVolumeFraction() const override { return m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::WATER]]; }
  real64 getNonWettingMinVolumeFraction() const override { return m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::GAS]]; }

private:
  void postInputInitialization() override;
  void initializePreSubGroups() override;

  string m_waterRelPermTableName;
  string m_gasRelPermTableName;
  string m_oilRelPermTableName;
  array1d< real64 > m_phaseMinVolumeFraction;
  TableFunction::KernelWrapper m_water;
  TableFunction::KernelWrapper m_gas;
  TableFunction::KernelWrapper m_oil;
};

} // namespace constitutive
} // namespace geos

#endif
