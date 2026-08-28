#include "DualContinuumCrossFlow.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementRegionBase.hpp"
#include "mesh/CellElementRegion.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumFlowSolverBase.hpp"
#include "constitutive/gravityDrainagePressure/SimpleGravityDrainagePressure.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

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
  this->registerWrapper( viewKeyStruct::gravityDrainageFlag(), &m_gravityDrainageFlag ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDefaultValue( 0 ).
        setDescription( "flag of gravity drainage" );

  registerWrapper( viewKeyStruct::shapeFactorTypeString(), &m_shapeFactorType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Matrix-fracture shape factor formulation. Valid options:\n* " +
                    EnumStrings< ShapeFactorType >::concat( "\n* " ) ).
    setApplyDefaultValue( ShapeFactorType::Kazemi );

  registerWrapper( viewKeyStruct::shapeFactorValueString(), &m_shapeFactorValue ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Directly specified matrix-fracture shape factor sigma [m^-2], "
                    "used when shapeFactorType=direct." );

  registerWrapper( viewKeyStruct::interporosityExchangeCoefficientString(),
                   &m_interporosityExchangeCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direct interporosity exchange coefficient Gamma [Pa^{-1} s^{-1}]. "
                    "When > 0, bypasses the Kazemi shape-factor formula and uses "
                    "transmissibility = Gamma * mu * V_element." );

  registerWrapper( viewKeyStruct::fractureVolumeFractionString(),
                   &m_fractureVolumeFraction ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fracture (secondary continuum) volume fraction v_f used by the "
                    "multi-porosity effective storage matrix. It is required for "
                    "co-located dual-continuum meshes; the matrix fraction is "
                    "v_m=1-v_f and v_m+v_f=1." );

  registerWrapper( viewKeyStruct::effectiveMatrixStorageString(), &m_effectiveMatrixStorage ).
    setApplyDefaultValue( 0.0 ).setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direct effective matrix storage coefficient Sbar_mm [1/Pa] for the multi-porosity "
                    "storage matrix. When positive, bypasses the intrinsic-property reconstruction of "
                    "Sbar_mm. Single-phase inputs include the fluid compressibility contribution; "
                    "compositional inputs should provide the skeleton part used by the compositional "
                    "cross-storage correction." );
  registerWrapper( viewKeyStruct::effectiveFractureStorageString(), &m_effectiveFractureStorage ).
    setApplyDefaultValue( 0.0 ).setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direct effective fracture storage coefficient Sbar_ff [1/Pa] for the multi-porosity "
                    "storage matrix. When positive, bypasses the intrinsic-property reconstruction of "
                    "Sbar_ff. Single-phase inputs include the fluid compressibility contribution; "
                    "compositional inputs should provide the skeleton part used by the compositional "
                    "cross-storage correction." );
  registerWrapper( viewKeyStruct::effectiveCrossStorageString(), &m_effectiveCrossStorage ).
    setApplyDefaultValue( 0.0 ).setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direct effective off-diagonal storage coefficient Sbar_mf=Sbar_fm [1/Pa] for the "
                    "multi-porosity storage matrix. This value is usually negative. It is still multiplied "
                    "by crossStorageOffDiagScale, matching the intrinsic reconstruction path." );

  registerWrapper( viewKeyStruct::crossStorageOffDiagScaleString(), &m_crossStorageOffDiagScale ).
    setApplyDefaultValue( 1.0 ).setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Scale on the multi-porosity off-diagonal storage (1 = paper bulk-Kbar value). "
                    "Values other than 1 are retained for compatibility with historical diagnostic "
                    "scans and should not be used as a physical correction without separate evidence." );
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
  GEOS_ERROR_IF( !elemManager.hasRegion( regionName ),
                 GEOS_FMT( "Region '{}' was not found in '{}'", regionName, elemManager.getName() ) );
  return elemManager.getRegion( regionName ).getIndexInParent();
}

void DualContinuumCrossFlow::setupCrossFlow( DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                             MeshLevel & meshMatrix,
                                             MeshLevel & meshFracture )
{


  if( m_matrixRegionList.size() == 0 || m_fractureRegionList.size() == 0 )
  {
    GEOS_ERROR( "Matrix region list or fracture region list is empty." );
    return;
  }

  // 0.1 分别获取两个网格的 Region 管理器
  ElementRegionManager const & elemManagerMatrix = meshMatrix.getElemManager();
  ElementRegionManager const & elemManagerFracture = meshFracture.getElemManager();

  string meshMatrixBodyName = meshMatrix.getParent().getParent().getName();
  string meshFractureBodyName = meshFracture.getParent().getParent().getName();



  // 2 构建窜流项stencil
  // 2.1 预估 Stencil 大小 (Calculate total size)
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
      cellRegion->forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
      {
        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
        for( localIndex i = 0; i < subRegion.size(); ++i )
        {
          if( ghostRank[i] < 0 )
          {
            ++totalConnections;
          }
        }
      } );
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

  // 2.2 循环所有 Region，构建 Stencil
  elemManagerMatrix.forElementRegions( [&]( ElementRegionBase const & elemRegionMatrix )
  {

    //根据matrix的名字找到对应的fracture名字与region
    string regionName = elemRegionMatrix.getName();
    localIndex couplingRegionIndexInList = findRegionIndexInList( regionName );
    if( couplingRegionIndexInList < 0 )
    {
      return;
    }
    string fractureRegionName = m_fractureRegionList[couplingRegionIndexInList];
    localIndex regionFractureIdx = findRegionIndexInRegionManager( elemManagerFracture, fractureRegionName );
    localIndex const regionMatrixIdx = elemRegionMatrix.getIndexInParent();
    ElementRegionBase const & elemRegionFracture = elemManagerFracture.getRegion( fractureRegionName );
    // 校验大小是否匹配
    if( elemRegionMatrix.getNumberOfElements() != elemRegionFracture.getNumberOfElements() )
    {
      GEOS_ERROR( "Region size mismatch between matrix region " << regionName << " and fracture region " << fractureRegionName );
      return;
    }
    elemRegionMatrix.forElementSubRegions( [&]( ElementSubRegionBase const & elementSubRegionMatrix )
    {
      localIndex const subRegionIdx = elementSubRegionMatrix.getIndexInParent();//默认matrix 与 fracture 的subregion相同
      ElementSubRegionBase const & elementSubRegionFracture = elemRegionFracture.getSubRegion( subRegionIdx );
      auto const & cellVolumeArrayViewMatrix = elementSubRegionMatrix.getReference< array1d< real64 > >( "elementVolume" );
      arrayView1d< integer const > const matrixGhostRank = elementSubRegionMatrix.ghostRank();
      arrayView1d< localIndex const > const matrixToFractureConnectivity =
        elementSubRegionMatrix.getReference< array1d< localIndex > >( "mesh1ToMesh2Connectivity" );
      for( localIndex i = 0; i < elementSubRegionMatrix.size(); i++ )
      {
        if( matrixGhostRank[i] >= 0 )
        {
          continue;
        }
        localIndex const fractureElementIndex = matrixToFractureConnectivity[i];
        GEOS_ERROR_IF( fractureElementIndex < 0 || fractureElementIndex >= elementSubRegionFracture.size(),
                       "Invalid dual-continuum matrix-to-fracture connectivity for matrix region "
                       << regionName << ", subregion " << elementSubRegionMatrix.getName()
                       << ", local element " << i << ": mapped fracture element "
                       << fractureElementIndex << " is outside fracture subregion "
                       << elementSubRegionFracture.getName() << " size "
                       << elementSubRegionFracture.size() );
        localIndex regionIndices[2] = { regionMatrixIdx, regionFractureIdx };
        localIndex subRegionIndices[2] = { subRegionIdx, subRegionIdx }; // Assuming default subregion 0
        localIndex elementIndices[2] = { i, fractureElementIndex };
        // Compute Geometric Weights [Wx, Wy, Wz]
        // 薄板形状的形状因子，乘体积是因为算这个网格的交换而不是每单位体积的交换
        // The shape factor formulation is selectable:
        //   - direct: sigma supplied directly [m^-2] (isotropic), W = (sigma/3) * V
        //   - kazemi: sigma = 4(1/Lx^2 + 1/Ly^2 + 1/Lz^2), W_i = 4 * V / L_i^2
        // 修正：从数组中获取第 i 个单元的体积
        real64 const Volume = cellVolumeArrayViewMatrix[i];
        real64 shapeFactory[3];
        if( m_shapeFactorType == ShapeFactorType::Direct )
        {
          real64 const sigmaPerDir = m_shapeFactorValue / 3.0;
          shapeFactory[0] = sigmaPerDir * Volume;
          shapeFactory[1] = sigmaPerDir * Volume;
          shapeFactory[2] = sigmaPerDir * Volume;
        }
        else // Kazemi
        {
          shapeFactory[0] = 4.0 * Volume * invLx2;
          shapeFactory[1] = 4.0 * Volume * invLy2;
          shapeFactory[2] = 4.0 * Volume * invLz2;
        }
        // Add to Stencil
        m_stencil.add( 2, regionIndices, subRegionIndices, elementIndices, shapeFactory, ConnIdx );
        ConnIdx++;
      }
    } );
  } );


  // 重力排驱压力初始化
  // 1.1 初始化重力排驱压力,在flowsolver中后初始化
  // setupGravityDrainagePressure(meshMatrix,meshFracture,9.81);
}


void DualContinuumCrossFlow::setupGravityDrainagePressure( MeshLevel & meshMatrix,
                                                           MeshLevel & fractureMatrix,
                                                           real64 const & gravityCoefficient)
{
  ElementRegionManager & elemManagerMatrix = meshMatrix.getElemManager();
  ElementRegionManager & elemManagerFracture = fractureMatrix.getElemManager();

/*
  typename TYPEOFREF(m_stencil)::KernelWrapper stencilWrapper = m_stencil.createKernelWrapper();


  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity >;

  //在多相计算中可能需要饱和度

  string const nameMatrix = "matrix";
  string const nameFracture = "fracture";

  SinglePhaseFluidAccessors matrixFluidAccessors(elemRegionMatrix,nameMatrix);
  SinglePhaseFluidAccessors fractureFluidAccessors(elemRegionFracture,nameFracture);

  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  forAll< parallelDevicePolicy<> >( stencilWrapper.size() ,[=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {


  });

*/


  if( m_gravityDrainageFlag )
  {
    elemManagerMatrix.forElementRegions( [&]( ElementRegionBase const & elemRegionMatrix ){
      //根据matrix的名字找到对应的fracture名字与region
      string regionName = elemRegionMatrix.getName();
      localIndex couplingRegionIndexInList = findRegionIndexInList( regionName );
      if( couplingRegionIndexInList < 0 ) return;
      string fractureRegionName = m_fractureRegionList[couplingRegionIndexInList];

      ElementRegionBase const & elemRegionFracture = elemManagerFracture.getRegion( fractureRegionName );

      elemRegionMatrix.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const subRegionIndex, CellElementSubRegion const & matrixSubRegion )
      {
        string const & fluidNameMatrix = matrixSubRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

        dataRepository::Group const & matrixConstitutiveModels = matrixSubRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

        ElementSubRegionBase const & fractureSubRegion = elemRegionFracture.getSubRegion( subRegionIndex );
        string const & fluidNameFracture = fractureSubRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

        dataRepository::Group const & fractureConstitutiveModels = fractureSubRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

        dataRepository::Group const & constitutiveModels = matrixSubRegion.getConstitutiveModels();
        constitutiveModels.forSubGroups< GravityDrainagePressureBase >( [&]( GravityDrainagePressureBase const & gdModel ){

          // Check if fluid is compositional multiphase (MultiFluidBase) or single-phase (SingleFluidBase)
          constitutive::MultiFluidBase const * multiFluidMatrix = matrixConstitutiveModels.getGroupPointer< constitutive::MultiFluidBase >( fluidNameMatrix );
          if( multiFluidMatrix )
          {
            arrayView3d< real64 const > const matrixPhaseMassDensity = multiFluidMatrix->phaseMassDensity();
            arrayView2d< real64 const > const matrixPhaseVolumeFraction =
              matrixSubRegion.getField< fields::flow::phaseVolumeFraction >();

            constitutive::MultiFluidBase const & fluidFracture = fractureConstitutiveModels.getGroup< constitutive::MultiFluidBase >( fluidNameFracture );
            arrayView3d< real64 const > const fracturePhaseMassDensity = fluidFracture.phaseMassDensity();
            arrayView2d< real64 const > const fracturePhaseVolumeFraction =
              fractureSubRegion.getField< fields::flow::phaseVolumeFraction >();

            gdModel.setupGravityDrainagePressureFromPhaseMassDensities(
              matrixPhaseMassDensity,
              matrixPhaseVolumeFraction,
              fracturePhaseMassDensity,
              fracturePhaseVolumeFraction,
              gravityCoefficient,
              m_fracSpacingLz );
          }
          else
          {
            constitutive::SingleFluidBase const * singleFluidMatrix = matrixConstitutiveModels.getGroupPointer< constitutive::SingleFluidBase >( fluidNameMatrix );
            constitutive::SingleFluidBase const * singleFluidFracture = fractureConstitutiveModels.getGroupPointer< constitutive::SingleFluidBase >( fluidNameFracture );
            if( singleFluidMatrix && singleFluidFracture )
            {
              gdModel.setupGravityDrainagePressure( singleFluidMatrix->density(), singleFluidFracture->density(), gravityCoefficient, m_fracSpacingLz );
            }
          }
        });
      });
    } );
  }
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
//using DualContinuumSinglePhaseSolver = DualContinuumFlowSolverBase< SinglePhaseBase, SinglePhaseBase >;

//REGISTER_CATALOG_ENTRY( DualContinuumSinglePhaseSolver, DualContinuumCrossFlow, string const &, Group * const )

} // namespace geos
