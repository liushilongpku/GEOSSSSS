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
 * @file MultiphaseDualContinuumPoromechanics.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEDUALCONTINUUMPOROMECHANICS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEDUALCONTINUUMPOROMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MGR strategy for compositional multiphase dual-continuum poromechanics.
 *
 * The dual-continuum poromechanics solver registers fields in this order:
 *   - compositional variables on both matrix and fracture mesh supports
 *   - nodal displacements
 *
 * This strategy first eliminates displacement labels using a mechanics AMG F-solver, then reduces both
 * compositional continua to their pressure label. The final coarse system contains pressure DOFs on both the matrix
 * and fracture mesh supports and is solved with BoomerAMG.
 */
class MultiphaseDualContinuumPoromechanics : public MGRStrategyBase< 3 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit MultiphaseDualContinuumPoromechanics( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] +
                                                                 numComponentsPerField[1] ) )
  {
    GEOS_ERROR_IF( numComponentsPerField.size() != 2,
                   "MultiphaseDualContinuumPoromechanics MGR requires exactly two fields: "
                   "dual-continuum flow variables and displacement." );

    HYPRE_Int const numFlowLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    HYPRE_Int const displacementLabelStart = numFlowLabels;

    GEOS_ERROR_IF( numFlowLabels < 2 || numComponentsPerField[1] != 3,
                   "MultiphaseDualContinuumPoromechanics MGR expects compositional flow variables "
                   "with at least pressure plus one component, followed by 3 displacement components." );

    // Level 0: eliminate displacement degrees of freedom.
    m_labels[0].resize( displacementLabelStart );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );

    // Level 1: eliminate the local volume-constraint component.
    m_labels[1].resize( numFlowLabels - 1 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 0 );

    // Level 2: eliminate the remaining component variables and keep pressure.
    m_labels[2].push_back( 0 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    // Level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::none;
    m_levelInterpType[2]          = MGRInterpolationType::injection;
    m_levelRestrictType[2]        = MGRRestrictionType::blockColLumped;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[2] = 1;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams MGR configuration parameters
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, mgrData );
    setMechanicsFSolver( precond, mgrData, mgrParams.separateComponents );
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEDUALCONTINUUMPOROMECHANICS_HPP_*/
