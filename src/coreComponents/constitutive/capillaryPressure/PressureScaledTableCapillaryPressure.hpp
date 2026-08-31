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
 * @file PressureScaledTableCapillaryPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_PRESSURESCALEDTABLECAPILLARYPRESSURE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_PRESSURESCALEDTABLECAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Table-based capillary pressure with optional pressure dependence.
 *
 * This model behaves like TableCapillaryPressure, but gas/oil capillary
 * pressure can either be scaled by a pressure-dependent factor:
 *
 *   Pc_go = sigma(P) / sigma_ref * Pc_go_table   (Thomas 1983 Eq. 28)
 *
 * where sigma(P) is provided as a table, or supplied directly as a
 * two-dimensional Pc_go(Sg, P) table for pressure-dependent pseudofunctions.
 *
 * It is provided as a separate model (instead of extending TableCapillaryPressure)
 * so that the existing TableCapillaryPressure behavior is left completely untouched.
 */
class PressureScaledTableCapillaryPressure : public CapillaryPressureBase
{
public:

  /// order of the phase properties for three-phase flow
  struct ThreePhasePairPhaseType
  {
    enum : integer
    {
      INTERMEDIATE_WETTING = 0,   ///< index for intermediate-wetting
      INTERMEDIATE_NONWETTING = 1 ///< index for intermediate-non-wetting
    };
  };

  PressureScaledTableCapillaryPressure( std::string const & name, dataRepository::Group * const parent );

  static std::string catalogName() { return "PressureScaledTableCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
public:

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                    arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac,
                    arrayView1d< real64 const > const & pressure,
                    bool const hasPressureScaling,
                    TableFunction::KernelWrapper const & pressureScalingWrapper,
                    bool const hasPressureDependentTable,
                    TableFunction::KernelWrapper const & pressureDependentTableWrapper );

    GEOS_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  real64 const pressure,
                  arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                  arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Array of kernel wrappers for the capillary pressures
    /// Is of size 1 for two-phase flow, and of size 2 for three-phase flow
    arrayView1d< TableFunction::KernelWrapper const > const m_capPresKernelWrappers;

    /// Cell-wise reference phase pressure (read from the parent sub-region)
    arrayView1d< real64 const > const m_pressure;

    /// Whether the gas/oil capillary pressure is scaled by a pressure-dependent factor
    bool m_hasPressureScaling;

    /// Kernel wrapper for the pressure-dependent scaling factor sigma(P)/sigma_ref
    TableFunction::KernelWrapper m_pressureScalingWrapper;

    /// Whether gas/oil capillary pressure is supplied directly as Pc(Sg, pressure)
    bool m_hasPressureDependentTable;

    /// Kernel wrapper for the pressure-dependent gas/oil capillary-pressure table
    TableFunction::KernelWrapper m_pressureDependentTableWrapper;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * wettingNonWettingCapPresTableNameString() { return "wettingNonWettingCapPressureTableName"; }
    static constexpr char const * wettingIntermediateCapPresTableNameString() { return "wettingIntermediateCapPressureTableName"; }
    static constexpr char const * nonWettingIntermediateCapPresTableNameString() { return "nonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * capPresWrappersString() { return "capPresWrappers"; }
    /// Optional pressure-dependent scaling factor table sigma(P)/sigma_ref for the gas/oil capillary pressure.
    static constexpr char const * pressureScalingTableNameString() { return "pressureScalingTableName"; }
    /// Optional two-dimensional gas/oil capillary-pressure table Pc(Sg, pressure).
    static constexpr char const * pressureDependentTableNameString() { return "pressureDependentTableName"; }

  };

private:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /// Capillary pressure table names (one for each phase in the wetting-non-wetting pair)
  string m_wettingNonWettingCapPresTableName;

  /// Capillary pressure table names (one for each phase in the wetting-intermediate pair)
  string m_wettingIntermediateCapPresTableName;

  /// Capillary pressure table names (one for each phase in the non-wetting-intermediate pair)
  string m_nonWettingIntermediateCapPresTableName;

  /// Capillary pressure kernel wrapper for the first pair (wetting-intermediate if NP=3, wetting-non-wetting otherwise)
  array1d< TableFunction::KernelWrapper > m_capPresKernelWrappers;

  /// Optional name of the pressure-dependent scaling factor table sigma(P)/sigma_ref (empty = no scaling)
  string m_pressureScalingTableName;

  /// Optional name of a two-dimensional gas/oil capillary-pressure table Pc(Sg, pressure)
  string m_pressureDependentTableName;

};

GEOS_HOST_DEVICE
inline void
PressureScaledTableCapillaryPressure::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           real64 const pressure,
           arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
           arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  if( m_hasPressureDependentTable && ipGas >= 0 )
  {
    if( ipWater >= 0 && ipOil >= 0 )
    {
      using TPT = PressureScaledTableCapillaryPressure::ThreePhasePairPhaseType;
      phaseCapPres[ipWater] =
        m_capPresKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipWater],
                                                                    &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
    }
    real64 const coordinates[2] = { phaseVolFraction[ipGas], pressure };
    real64 derivatives[2] = { 0.0, 0.0 };
    phaseCapPres[ipGas] = -m_pressureDependentTableWrapper.compute( coordinates, derivatives );
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] = -derivatives[0];
  }

  if( !m_hasPressureDependentTable && ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    using TPT = PressureScaledTableCapillaryPressure::ThreePhasePairPhaseType;

    // water-oil capillary pressure
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipWater],
                                                                  &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );

    // gas-oil capillary pressure
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipGas],
                                                                     &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( !m_hasPressureDependentTable && ipWater < 0 )
  {
    // put capillary pressure on the non-wetting phase
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipGas],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipOil < 0 || ipGas < 0 )
  {
    // put capillary pressure on the wetting phase
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
  }

  // Optional pressure-dependent scaling of the gas/oil capillary pressure:
  // Pc_go = sigma(P)/sigma_ref * Pc_go_table (Thomas 1983 Eq. 28).
  if( m_hasPressureScaling && !m_hasPressureDependentTable && ipGas >= 0 )
  {
    real64 const scaling = m_pressureScalingWrapper.compute( &pressure );
    phaseCapPres[ipGas] *= scaling;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= scaling;
  }
}

GEOS_HOST_DEVICE
inline void
PressureScaledTableCapillaryPressure::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geos::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_pressure[k],
           m_phaseCapPressure[k][q],
           m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
}

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_PRESSURESCALEDTABLECAPILLARYPRESSURE_HPP
