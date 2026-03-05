#include "DualContinuumCrossFlow.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementRegionBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumFlowSolver.hpp"
namespace geos
{

using namespace dataRepository;

DualContinuumCrossFlow::DualContinuumCrossFlow( string const & name,
                                                Group * const parent )
  : Group( name, parent )
{

  setInputFlags( InputFlags::OPTIONAL );
  
  // Register parameters to be read from XML
  registerWrapper( viewKeyStruct::fractureSpacingLxString(), &m_fracSpacingLx ).
                                                                                 setInputFlag( InputFlags::REQUIRED ).
                                                                                 setDescription( "Fracture spacing in X direction" );

  registerWrapper( viewKeyStruct::fractureSpacingLyString(), &m_fracSpacingLy ).
                                                                                 setInputFlag( InputFlags::REQUIRED ).
                                                                                 setDescription( "Fracture spacing in Y direction" );

  registerWrapper( viewKeyStruct::fractureSpacingLzString(), &m_fracSpacingLz ).
                                                                                 setInputFlag( InputFlags::REQUIRED ).
                                                                                 setDescription( "Fracture spacing in Z direction" );

  registerWrapper( viewKeyStruct::matrixRegionList(), &m_matrixRegionList ).
                                                                             setInputFlag( InputFlags::REQUIRED ).
                                                                             setDescription( "List of matrix regions" );

  registerWrapper( viewKeyStruct::fractureRegionList(), &m_fractureRegionList ).
                                                                                 setInputFlag( InputFlags::REQUIRED ).
                                                                                 setDescription( "List of fracture regions" );

  // Register the stencil so it shows up in output/checkpoint

  registerWrapper( viewKeyStruct::DualContinuumStencilString(), &m_stencil ).
                          setRestartFlags( RestartFlags::NO_WRITE );

}

localIndex DualContinuumCrossFlow::findRegionIndexInList( string const & regionName )
{
  auto found = find( m_matrixRegionList.begin(),
                     m_matrixRegionList.end(),
                     regionName );    // std::find 接收迭代器：begin()（首元素）、end()（尾后）

  return found != m_matrixRegionList.end() ? distance( m_matrixRegionList.begin(), found ) : -1;    // 迭代器转换为索引：distance(起始迭代器, 目标迭代器)
}


localIndex DualContinuumCrossFlow::findRegionIndexInRegionManager( ElementRegionManager const & elemManager,
                                                                   string const & regionName )
{
  localIndex idx = 0;
  elemManager.forElementRegions( [&]( ElementRegionBase const & elemRegion )
                                 {
                                   if( elemRegion.getName() == regionName && idx < elemRegion.getNumberOfElements() )
                                   {
                                     return;
                                   }
                                   idx++;
                                 } );

  return idx;
}

void DualContinuumCrossFlow::initialize( MeshLevel & meshMatrix,
                                         MeshLevel & meshFracture )
{
  if( m_matrixRegionList.size() == 0 || m_fractureRegionList.size() == 0 )
  {
    GEOS_ERROR( "Matrix region list or fracture region list is empty." );
    return;
  }
  // 1. 分别获取两个网格的 Region 管理器
  ElementRegionManager const & elemManagerMatrix = meshMatrix.getElemManager();
  ElementRegionManager const & elemManagerFracture = meshFracture.getElemManager();

  string meshMatrixBodyName = meshMatrix.getParent().getParent().getName();
  string meshFractureBodyName = meshFracture.getParent().getParent().getName();

  // 2. 预估 Stencil 大小 (Calculate total size)
  localIndex totalConnections = 0;
  for( auto const & regionName: m_matrixRegionList )
  {
    // 需要先找到 Region 才能获取 size
    // 为了安全，这里也使用辅助查找逻辑，或者假设名字存在
    // 这里简化处理，直接通过 Manager 获取
    // 注意：getAttribute/size 等方法建议在 CellElementRegion 上调用
    ElementRegionBase const & regionBase = elemManagerMatrix.getRegion( regionName );
    if( auto const * cellRegion = dynamic_cast< CellElementRegion const * >( &regionBase ) )
    {
      totalConnections += cellRegion->getNumberOfElements();
    }
    else
    {
      GEOS_ERROR( "Matrix region " << regionName << " is not a CellElementRegion." );
    }
  }
  m_stencil.reserve( totalConnections );


  // --- D. 获取体积属性 ---
  // 注意：getReference 返回的是数组的引用，不是标量
  // 假设属性名为 "elementVolume" 或 "volume"，请根据你的 GEOS 版本确认

  real64 const invLx2 = ( m_fracSpacingLx > 0 ) ? 1.0 / ( m_fracSpacingLx * m_fracSpacingLx ) : 0.0;
  real64 const invLy2 = ( m_fracSpacingLy > 0 ) ? 1.0 / ( m_fracSpacingLy * m_fracSpacingLy ) : 0.0;
  real64 const invLz2 = ( m_fracSpacingLz > 0 ) ? 1.0 / ( m_fracSpacingLz * m_fracSpacingLz ) : 0.0;

  localIndex ConnIdx = 0;

  // 3.循环所有 Region，构建 Stencil
  localIndex regionMatrixIdx = 0;
  elemManagerMatrix.forElementRegions( [&]( ElementRegionBase const & elemRegionMatrix ){

    //根据matrix的名字找到对应的fracture名字与region
    string regionName = elemRegionMatrix.getName();
    localIndex couplingRegionIndexInList = findRegionIndexInList( regionName );
    string fractureRegionName = m_fractureRegionList[couplingRegionIndexInList];
    localIndex regionFractureIdx = findRegionIndexInRegionManager( elemManagerFracture, fractureRegionName );

    ElementRegionBase const & elemRegionFracture = elemManagerFracture.getRegion( fractureRegionName );

    // 校验大小是否匹配
    if( elemRegionMatrix.getNumberOfElements() != elemRegionFracture.getNumberOfElements() ){
      GEOS_ERROR( "Region size mismatch between matrix region " << regionName << " and fracture region " << fractureRegionName );
      return;
    }


    localIndex subRegionIdx = 0;//默认matrix 与 fracture 的subregion相同
    elemRegionMatrix.forElementSubRegions( [&]( ElementSubRegionBase const & elementSubRegionMatrix ){
      //目前只需要获得matrix的体积属性,因此无需获取fracture的subregion
      //ElementSubRegionBase const & elementSubRegionFracture = elemRegionFracture.getSubRegions()[subRegionIdx];
      auto const & cellVolumeArrayViewMatrix = elementSubRegionMatrix.getReference< array1d< real64 > >( "elementVolume" );

      for( localIndex i = 0; i < elementSubRegionMatrix.size(); i++ )
      {
        localIndex regionIndices[2] = { regionMatrixIdx, regionFractureIdx };
        localIndex subRegionIndices[2] = { subRegionIdx, subRegionIdx }; // Assuming default subregion 0
        localIndex elementIndices[2] = { i, i }; // 1-to-1 mapping

        // Compute Geometric Weights [Wx, Wy, Wz]
        // W = 4 * V / L^2
        // 修正：从数组中获取第 i 个单元的体积
        real64 const Volume = cellVolumeArrayViewMatrix[i];
        real64 weights[3];
        weights[0] = 4.0 * Volume * invLx2;
        weights[1] = 4.0 * Volume * invLy2;
        weights[2] = 4.0 * Volume * invLz2;

        // Add to Stencil
        m_stencil.add( 2, regionIndices, subRegionIndices, elementIndices, weights, ConnIdx );
        ConnIdx++;
      }
      subRegionIdx++;
    } );
    regionMatrixIdx++;
  } );
}

// Example implementation of assembly (you need to adapt types to match SinglePhaseBase)
template< typename SOLVER_TRAITS >
void DualContinuumCrossFlow::assembleCouplingTerms( MeshLevel & mesh,
                                                    DofManager const & dofManager,
                                                    typename SOLVER_TRAITS::PhysicsData & physicsData,
                                                    typename SOLVER_TRAITS::ResidualType & residual ) const
{
  // Implementation depends on Solver Physics Data structure
  // Logic:
  // 1. auto wrapper = m_stencil.createKernelWrapper();
  // 2. RAJA::forall( ... ) {
  //      wrapper.computeTransferCoeff(i, perm, dPerm, T, dT);
  //      Flux = T * (P_frac - P_mat);
  //      atomicAdd( res_mat, -Flux );
  //      atomicAdd( res_frac, +Flux );
  // }
}
//using DualContinuumSinglePhaseSolver = DualContinuumFlowSolver< SinglePhaseBase, SinglePhaseBase >;

//REGISTER_CATALOG_ENTRY( DualContinuumSinglePhaseSolver, DualContinuumCrossFlow, string const &, Group * const )

} // namespace geos

