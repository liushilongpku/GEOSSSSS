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

#ifndef GEOS_FIELDSPECIFICATION_RIGIDBOUNDARY_HPP
#define GEOS_FIELDSPECIFICATION_RIGIDBOUNDARY_HPP

#include "FieldSpecificationBase.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/NodeManager.hpp"

namespace geos
{

class RigidBoundary : public FieldSpecificationBase
{
public:
  RigidBoundary( string const & name, Group * parent );
  RigidBoundary() = delete;
  virtual ~RigidBoundary() = default;
  RigidBoundary( RigidBoundary const & ) = delete;
  RigidBoundary( RigidBoundary && ) = default;
  RigidBoundary & operator=( RigidBoundary const & ) = delete;
  RigidBoundary & operator=( RigidBoundary && ) = delete;

  static string catalogName() { return "RigidBoundary"; }

  /// Apply distributed pressure load to RHS.
  void applyLoad( real64 const time,
                  arrayView1d< globalIndex const > const nodeDofNumber,
                  globalIndex const dofRankOffset,
                  FaceManager const & faceManager,
                  NodeManager const & nodeManager,
                  SortedArrayView< localIndex const > const & targetSet,
                  arrayView1d< real64 > const & localRhs ) const;

  /// Apply the rigid-platen load AND enforce the equal-displacement (rigid
  /// frictionless platen) constraint directly in the linear system, via a
  /// penalty that ties every boundary-node component-DOF to the (lagged)
  /// weighted-mean boundary displacement.  Unlike enforceConstraint(), this is
  /// assembled into the Jacobian/residual so Newton converges (the post-solve
  /// projection used previously prevented convergence for non-uniform fields).
  void applyRigidConstraint( real64 const time,
                             arrayView1d< globalIndex const > const nodeDofNumber,
                             globalIndex const dofRankOffset,
                             FaceManager const & faceManager,
                             NodeManager const & nodeManager,
                             SortedArrayView< localIndex const > const & targetSet,
                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                             arrayView1d< real64 > const & localRhs ) const;

  /// Force all boundary nodes to have the same displacement (modifies displacement field directly).
  /// Legacy post-solve projection; retained for reference but no longer used in the implicit solve.
  void enforceConstraint( FaceManager const & faceManager,
                          NodeManager & nodeManager,
                          SortedArrayView< localIndex const > const & targetSet ) const;

  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {
  };

protected:
  virtual void postInputInitialization() override;
};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_RIGIDBOUNDARY_HPP */
