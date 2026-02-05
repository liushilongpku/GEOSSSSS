//
// Created by hello on 2026/1/19.
//


#include "DualContinuumStencil.hpp"


namespace geos
{
DualContinuumStencil::DualContinuumStencil()
  : StencilBase()
{
}
  //end DualContinuumWrapper constructor nothing to construct
void DualContinuumStencil::reserve( localIndex const size )
{
  StencilBase::reserve( size );
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
  m_weights.resize( newSize, 3 );
  //当基质的渗透率远小于裂缝的渗透率时，仅需要对基质本身的传导系数进行设置。


  for( localIndex a=0; a<numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
  }
  for(localIndex dim=0; dim < 3; ++dim )
  {
    m_weights( oldSize, dim ) = weights[dim];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}
/*目前好像不需要这些内容
void DualContinuumStencil::addVectors( arrayView2d< real64 > const & transMultiplier )
{
  localIndex const oldSize = m_transMultiplier.size( 0 );
  localIndex const newSize = oldSize + 1;

  m_transMultiplier.resize( newSize );
  m_transMultiplier[oldSize] = transMultiplier;
}
*/

DualContinuumStencil::KernelWrapper
DualContinuumStencil::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights};
}

DualContinuumStencilWrapper::
DualContinuumStencilWrapper( IndexContainerType const & elementRegionIndices,
                               IndexContainerType const & elementSubRegionIndices,
                               IndexContainerType const & elementIndices,
                               WeightContainerType const & weights )
  : StencilWrapperBase( elementRegionIndices,
                        elementSubRegionIndices,
                        elementIndices,
                        weights )
{}

}//end namespace geos

