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

#include "RigidBoundary.hpp"

#include "common/MpiWrapper.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionBase.hpp"
#include "functions/FunctionManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

RigidBoundary::RigidBoundary( string const & name, Group * parent ):
  FieldSpecificationBase( name, parent )
{
  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Constrained coordinate component: 0=x, 1=y, 2=z." );

  getWrapper< real64 >( FieldSpecificationBase::viewKeyStruct::scaleString() ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Applied pressure (Pa)." );
}


void RigidBoundary::postInputInitialization()
{
  int const component = getComponent();
  GEOS_ERROR_IF( component < 0 || component > 2,
                 getDataContext() << ": " << viewKeyStruct::componentString()
                                 << " must be 0, 1, or 2." );
}


void RigidBoundary::applyLoad( real64 const time,
                               arrayView1d< globalIndex const > const nodeDofNumber,
                               globalIndex const dofRankOffset,
                               FaceManager const & faceManager,
                               NodeManager const & nodeManager,
                               SortedArrayView< localIndex const > const & targetSet,
                               arrayView1d< real64 > const & localRhs ) const
{
  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();

  int const component = getComponent();

  real64 pressureMagnitude = 0.0;
  if( getFunctionName().empty() )
  {
    pressureMagnitude = getScale();
  }
  else
  {
    FunctionManager const & functionManager = FunctionManager::getInstance();
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( getFunctionName() );
    if( function.isFunctionOfTime() == 2 )
    {
      pressureMagnitude = getScale() * function.evaluate( &time );
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": Spatial functions are not supported." );
    }
  }

  // Per-node weights and total area
  std::map< localIndex, real64 > nodeWeight;
  real64 localTotalArea = 0.0;

  for( localIndex i = 0; i < targetSet.size(); ++i )
  {
    localIndex const kf = targetSet[ i ];
    localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
    real64 const area = faceArea[ kf ];
    localTotalArea += area;
    real64 const nodeContrib = area / numNodes;
    for( localIndex a = 0; a < numNodes; ++a )
    {
      nodeWeight[ faceToNodeMap( kf, a ) ] += nodeContrib;
    }
  }

  real64 const globalTotalArea = MpiWrapper::allReduce( localTotalArea, MpiWrapper::Reduction::Sum );
  GEOS_ERROR_IF( globalTotalArea <= 0.0,
                 getDataContext() << ": RigidBoundary has zero total area." );

  real64 localTotalWeight = 0.0;
  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    if( ghostRank[ nodeIdx ] < 0 ) localTotalWeight += weight;
  }
  real64 const globalTotalWeight = MpiWrapper::allReduce( localTotalWeight, MpiWrapper::Reduction::Sum );
  GEOS_ERROR_IF( globalTotalWeight <= 0.0,
                 getDataContext() << ": RigidBoundary has zero total weight." );

  real64 const totalForce = pressureMagnitude * globalTotalArea;

  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    globalIndex const dof = nodeDofNumber[ nodeIdx ] + component;
    localIndex const localRow = LvArray::integerConversion< localIndex >( dof - dofRankOffset );
    if( localRow >= 0 && localRow < localRhs.size() )
    {
      localRhs[ localRow ] += totalForce * weight / globalTotalWeight;
    }
  }
}


void RigidBoundary::applyRigidConstraint( real64 const time,
                                          arrayView1d< globalIndex const > const nodeDofNumber,
                                          globalIndex const dofRankOffset,
                                          FaceManager const & faceManager,
                                          NodeManager const & nodeManager,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) const
{
  // Penalty factor relative to the assembled diagonal stiffness.  Large enough
  // to tie the boundary DOFs together, small enough to keep the system well
  // conditioned for the iterative linear solver.
  constexpr real64 penaltyFactor = 1.0e4;

  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();
  solidMechanics::arrayViewConst2dLayoutTotalDisplacement const disp =
    nodeManager.getField< solidMechanics::totalDisplacement >();

  int const component = getComponent();

  // ---- applied pressure magnitude (constant or function of time) ----
  real64 pressureMagnitude = 0.0;
  if( getFunctionName().empty() )
  {
    pressureMagnitude = getScale();
  }
  else
  {
    FunctionManager const & functionManager = FunctionManager::getInstance();
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( getFunctionName() );
    GEOS_ERROR_IF( function.isFunctionOfTime() != 2,
                   getDataContext() << ": Only functions of time are supported." );
    pressureMagnitude = getScale() * function.evaluate( &time );
  }

  // ---- tributary node weights and total area ----
  std::map< localIndex, real64 > nodeWeight;
  real64 localTotalArea = 0.0;
  for( localIndex i = 0; i < targetSet.size(); ++i )
  {
    localIndex const kf = targetSet[ i ];
    localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
    localTotalArea += faceArea[ kf ];
    real64 const nodeContrib = faceArea[ kf ] / numNodes;
    for( localIndex a = 0; a < numNodes; ++a )
    {
      nodeWeight[ faceToNodeMap( kf, a ) ] += nodeContrib;
    }
  }
  real64 const globalTotalArea = MpiWrapper::allReduce( localTotalArea, MpiWrapper::Reduction::Sum );
  GEOS_ERROR_IF( globalTotalArea <= 0.0, getDataContext() << ": RigidBoundary has zero total area." );

  // total owned weight (for load distribution)
  real64 localTotalWeight = 0.0;
  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    if( ghostRank[ nodeIdx ] < 0 )
      localTotalWeight += weight;
  }
  real64 const globalTotalWeight = MpiWrapper::allReduce( localTotalWeight, MpiWrapper::Reduction::Sum );
  GEOS_ERROR_IF( globalTotalWeight <= 0.0, getDataContext() << ": RigidBoundary has zero total weight." );

  real64 const totalForce = pressureMagnitude * globalTotalArea;

  // (1) distributed external load on each owned boundary node
  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    if( ghostRank[ nodeIdx ] >= 0 )
      continue;
    globalIndex const dof = nodeDofNumber[ nodeIdx ] + component;
    localIndex const localRow = LvArray::integerConversion< localIndex >( dof - dofRankOffset );
    if( localRow >= 0 && localRow < localMatrix.numRows() )
      localRhs[ localRow ] += totalForce * weight / globalTotalWeight;
  }

  // (2) rigid (equal-displacement) constraint as a face-local graph-Laplacian
  //     penalty.  For each boundary face, drive its node component-DOFs toward
  //     the face mean: r_a += -beta*(u_a - mean_face), with Jacobian
  //     d r_a / d u_b = -beta*(delta_ab - 1/n) for nodes a,b of the face.
  //     The uniform (collective) mode is in the null space of this penalty, so
  //     the platen is still free to translate under the applied load, while
  //     within-face differences are suppressed; shared nodes chain faces
  //     together to enforce global platen uniformity.  All coupled DOFs (a,b)
  //     belong to the same element, so the entries already exist in the
  //     sparsity pattern.  localMatrix == -K, so penalty terms are added with a
  //     leading minus sign to match the residual Jacobian convention.
  for( localIndex i = 0; i < targetSet.size(); ++i )
  {
    localIndex const kf = targetSet[ i ];
    localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
    if( numNodes < 2 )
      continue;

    // face-mean displacement (current iterate) and a representative penalty scale
    real64 faceMeanDisp = 0.0;
    for( localIndex a = 0; a < numNodes; ++a )
      faceMeanDisp += disp( faceToNodeMap( kf, a ), component );
    faceMeanDisp /= numNodes;

    for( localIndex a = 0; a < numNodes; ++a )
    {
      localIndex const nodeA = faceToNodeMap( kf, a );
      if( ghostRank[ nodeA ] >= 0 )
        continue;
      globalIndex const dofA = nodeDofNumber[ nodeA ] + component;
      localIndex const rowA = LvArray::integerConversion< localIndex >( dofA - dofRankOffset );
      if( rowA < 0 || rowA >= localMatrix.numRows() )
        continue;

      arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( rowA );
      arraySlice1d< real64 > const entries = localMatrix.getEntries( rowA );
      localIndex const numEntries = localMatrix.numNonZeros( rowA );

      // penalty magnitude scaled by the assembled diagonal stiffness of this row
      real64 diagonal = 0.0;
      for( localIndex j = 0; j < numEntries; ++j )
      {
        if( columns[ j ] == dofA ) { diagonal = entries[ j ]; break; }
      }
      real64 const beta = penaltyFactor * LvArray::math::abs( diagonal );
      if( beta <= 0.0 )
        continue;

      // residual: r_a += -beta*(u_a - faceMean)
      localRhs[ rowA ] += -beta * ( disp( nodeA, component ) - faceMeanDisp );

      // Jacobian: add -beta*(delta_ab - 1/n) to (rowA, dof_b) for each face node b
      for( localIndex b = 0; b < numNodes; ++b )
      {
        globalIndex const dofB = nodeDofNumber[ faceToNodeMap( kf, b ) ] + component;
        real64 const coeff = -beta * ( ( a == b ? 1.0 : 0.0 ) - 1.0 / numNodes );
        for( localIndex j = 0; j < numEntries; ++j )
        {
          if( columns[ j ] == dofB ) { entries[ j ] += coeff; break; }
        }
      }
    }
  }
}


void RigidBoundary::enforceConstraint( FaceManager const & faceManager,
                                       NodeManager & nodeManager,
                                       SortedArrayView< localIndex const > const & targetSet ) const
{
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();
  solidMechanics::arrayView2dLayoutTotalDisplacement const disp =
    nodeManager.getField< solidMechanics::totalDisplacement >();

  int const component = getComponent();

  // Collect boundary nodes with tributary weights
  std::map< localIndex, real64 > nodeWeight;
  for( localIndex i = 0; i < targetSet.size(); ++i )
  {
    localIndex const kf = targetSet[ i ];
    localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
    real64 const nodeContrib = faceManager.faceArea()[ kf ] / numNodes;
    for( localIndex a = 0; a < numNodes; ++a )
    {
      nodeWeight[ faceToNodeMap( kf, a ) ] += nodeContrib;
    }
  }

  // Compute global weighted average displacement
  real64 localWeightedDisp = 0.0;
  real64 localTotalWeight = 0.0;

  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    if( ghostRank[ nodeIdx ] < 0 )
    {
      localWeightedDisp += weight * disp( nodeIdx, component );
      localTotalWeight += weight;
    }
  }

  real64 const globalWeightedDisp = MpiWrapper::allReduce( localWeightedDisp, MpiWrapper::Reduction::Sum );
  real64 const globalTotalWeight  = MpiWrapper::allReduce( localTotalWeight, MpiWrapper::Reduction::Sum );

  if( globalTotalWeight <= 0.0 ) return;

  real64 const uAvg = globalWeightedDisp / globalTotalWeight;

  // Set all boundary node displacements to uAvg
  for( auto const & [nodeIdx, weight] : nodeWeight )
  {
    GEOS_UNUSED_VAR( weight );
    disp( nodeIdx, component ) = uAvg;
  }
}


REGISTER_CATALOG_ENTRY( FieldSpecificationBase, RigidBoundary, string const &, Group * const )

} /* namespace geos */
