//
// Created by hello on 2026/1/19.
//


#include "physicsSolvers/multiphysics/dualContinuumCrossFlowComputeKernels/DualContinuumStencil.hpp"


namespace geos
{
DualContinuumStencil::DualContinuumStencil()
  : StencilBase()
{
  m_reciprocalFractureSpacingArray.resize(0,3);
}
  //end DualContinuumWrapper constructor nothing to construct
void DualContinuumStencil::reserve( localIndex const size )
{
  StencilBase::reserve( size );

  m_reciprocalFractureSpacingArray.reserve( 3 * size );
  m_transMultiplier.reserve( size );
}

void DualContinuumStencil::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOS_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in TPFA stencil should be 2" );

  localIndex const oldSize = m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;

  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, 1 );
  //当基质的渗透率远小于裂缝的渗透率时，仅需要对基质本身的传导系数进行设置。


  for( localIndex a=0; a<numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
  }
  m_weights( oldSize, 0 ) = weights[0];
  m_connectorIndices[connectorIndex] = oldSize;
}

void DualContinuumStencil::addVectors( real64 const & transMultiplier )
{
  localIndex const oldSize = m_transMultiplier.size( 0 );
  localIndex const newSize = oldSize + 1;

  m_transMultiplier.resize( newSize );
  m_transMultiplier[oldSize] = transMultiplier;
}


DualContinuumStencil::KernelWrapper
DualContinuumStencil::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights,
           m_reciprocalFractureSpacingArray,
           m_transMultiplier};
}

DualContinuumStencilWrapper::
DualContinuumStencilWrapper( IndexContainerType const & elementRegionIndices,
                               IndexContainerType const & elementSubRegionIndices,
                               IndexContainerType const & elementIndices,
                               WeightContainerType const & weights,
                               arrayView2d< real64 > const & reciprocalFractureSpacingArray,
                               arrayView1d< real64 > const & transMultiplier )
  : StencilWrapperBase( elementRegionIndices,
                        elementSubRegionIndices,
                        elementIndices,
                        weights ),
    m_reciprocalFractureSpacingArray(reciprocalFractureSpacingArray),
    m_transMultiplier( transMultiplier )
{}

}//end namespace geos

