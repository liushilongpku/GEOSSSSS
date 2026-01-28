//
// Created by hello on 2026/1/19.
//

#ifndef GEOSX_DUALCONTINUUMSTENCIL_HPP
#define GEOSX_DUALCONTINUUMSTENCIL_HPP

#include "StencilBase.hpp"

namespace geos
{
class DualContinuumStencilWrapper : public StencilWrapperBase< TwoPointStencilTraits >
{
public:
  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::ElementViewConst< VIEWTYPE >;


  DualContinuumStencilWrapper(IndexContainerType const & elementRegionIndices,
                              IndexContainerType const & elementSubRegionIndices,
                              IndexContainerType const & elementIndices,
                              WeightContainerType const & weights,
                              arrayView2d< real64 > const & m_reciprocalFractureSpacingArray,
                              arrayView1d< real64 > const & transMultiplier );


  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                       real64  & weight,
                       real64  & dWeight_dVar) const;

  localIndex size() const
  {
    return m_elementRegionIndices.size( 0 );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex stencilSize( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numPointsInFlux( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

  typename TwoPointStencilTraits::IndexContainerViewConstType m_MeshIndices;

private:

  arrayView2d< real64 > m_reciprocalFractureSpacingArray;
  arrayView1d< real64 > m_transMultiplier;

};//end class DualContinuumStencilWrapper

class DualContinuumStencil final : public StencilBase< TwoPointStencilTraits, DualContinuumStencil >
{
public:
  DualContinuumStencil();

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  void addVectors( real64 const & transMultiplier );

  virtual localIndex size() const override
  { return m_elementRegionIndices.size( 0 ); }

  virtual void reserve( localIndex const size ) override;

  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }


  using KernelWrapper = DualContinuumStencilWrapper;
  KernelWrapper createKernelWrapper() const;

  typename TwoPointStencilTraits::IndexContainerType m_MeshIndices;

private:

  array2d< real64 > m_reciprocalFractureSpacingArray;
  array1d< real64 > m_transMultiplier;
};
//end class DualContinuumWrapper

GEOS_HOST_DEVICE
inline void
DualContinuumStencilWrapper::
computeWeights( localIndex const iconn,
                CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                real64 & weight,
                real64 & dWeight_dVar  ) const//导数项也只有一个分量，因为传导系数考虑为只和基质的渗透率有关，而基质的渗透率仅与其自身压力
{
  localIndex const mi = m_MeshIndices[iconn][0];
  localIndex const er = m_elementRegionIndices[iconn][0];
  localIndex const esr = m_elementSubRegionIndices[iconn][0];
  localIndex const ei = m_elementIndices[iconn][0];
  //在这里，仅仅使用基质的渗透率，并且假设第一个是基质渗透率
  // Coeff = Weight * Perm
  weight = LvArray::tensorOps::AiBi< 3 >( m_reciprocalFractureSpacingArray[iconn], coefficient[er][esr][ei][0] );

  // B. 计算导数 (链式法则)
  // d(Coeff)/dVar = d(Weight * Perm)/dVar
  // 因为 Weight 是几何常数 (预计算好的)，所以：
  // dCoeff/dVar = Weight * dPerm/dVar

  // 注意：这里是一个标量乘向量的操作（如果 dPerm_dVar 是个梯度向量）
  // 或者简单的标量乘法
  dWeight_dVar = LvArray::tensorOps::AiBi< 3 >( m_reciprocalFractureSpacingArray[iconn], dCoeff_dVar[er][esr][ei][0] );
}

}//end namespace geos



#endif //GEOSX_DUALCONTINUUMSTENCIL_HPP
