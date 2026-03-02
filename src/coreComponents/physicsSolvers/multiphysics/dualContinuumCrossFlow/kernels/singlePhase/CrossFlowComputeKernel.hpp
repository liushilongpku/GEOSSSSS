

/**
 * @file CrossFlowComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASE_CROSSFLOWCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASE_CROSSFLOWCOMPUTEKERNEL_HPP

#include "CrossFlowComputeKernelBase.hpp"

namespace geos
{

namespace singlePhaseDualContinuumKernels
{

/**
 * @class FluxComputeKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class CrossFlowComputeKernel : public CrossFlowComputeKernelBase
{
public:

/// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

/// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

/// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

/// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

/// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

/**
 * @brief Constructor for the kernel interface
 * @param[in] rankOffset the offset of my MPI rank
 * @param[in] stencilWrapper reference to the stencil wrapper
 * @param[in] dofNumberAccessor
 * @param[in] singlePhaseFlowAccessors
 * @param[in] singlePhaseFluidAccessors
 * @param[in] permeabilityAccessors
 * @param[in] dt time step size
 * @param[inout] localMatrix the local CRS matrix
 * @param[inout] localRhs the local right-hand side vector
 */
  CrossFlowComputeKernel( globalIndex const rankOffsetM,
                          globalIndex const rankOffsetF,
                          STENCILWRAPPER const & stencilWrapper,
                          DofNumberAccessor const & dofMatrixNumberAccessor,
                          DofNumberAccessor const & dofFractureNumberAccessor,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                          PermeabilityAccessors const & permeabilityAccessors,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors_fracture,
                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors_fracture,
                          PermeabilityAccessors const & permeabilityAccessors_fracture,
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
    : CrossFlowComputeKernelBase( rankOffsetM,
                                  rankOffsetF,
                                  dofMatrixNumberAccessor,
                                  dofFractureNumberAccessor,
                                  singlePhaseFlowAccessors,
                                  singlePhaseFluidAccessors,
                                  permeabilityAccessors,
                                  singlePhaseFlowAccessors_fracture,
                                  singlePhaseFluidAccessors_fracture,
                                  permeabilityAccessors_fracture,
                                  dt,
                                  localMatrix,
                                  localRhs ),
      m_stencilWrapper( stencilWrapper ),
      m_seri( stencilWrapper.getElementRegionIndices() ),
      m_sesri( stencilWrapper.getElementSubRegionIndices() ),
      m_sei( stencilWrapper.getElementIndices() )
  {}
/**
 * @struct StackVariables
 * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
 */
  struct StackVariables
  {
  public:

    /**
   * @brief Constructor for the stack variables
   * @param[in] size size of the stencil for this connection
   * @param[in] numElems number of elements for this connection
   */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
        numFluxElems( numElems ),
        dofColIndices( size * numDof ),
        localFlux( numElems * numEqn ),
        localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;

    /// Number of elements for a given connection
    localIndex const numFluxElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres{} ;
    //real64 dTrans_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDof > dofColIndices;

    /// Storage for the face local residual vector (all equations except volume balance)
    stackArray1d< real64, maxNumElems * numEqn > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > localFluxJacobian;

  };

  /**
 * @brief Getter for the stencil size at this connection
 * @param[in] iconn the connection index
 * @return the size of the stencil at this connection
 */
  GEOS_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

  /**
 * @brief Getter for the number of elements at this connection
 * @param[in] iconn the connection index
 * @return the number of elements at this connection
 */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
 * @brief Performs the setup phase for the kernel.
 * @param[in] iconn the connection index
 * @param[in] stack the stack variables
 */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face

      globalIndex const offsetM = m_dofMatrixNumber[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )];
      globalIndex const offsetF = m_dofFractureNumber[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )];

      //这里处理
      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[jdof] = offsetM + jdof;
        stack.dofColIndices[numDof + jdof] = offsetF + jdof;
      }

  }

  /**
 * @brief Compute the local flux contributions to the residual and Jacobian
 * @tparam FUNC the type of the function that can be used to customize the computation of the flux
 * @param[in] iconn the connection index
 * @param[inout] stack the stack variables
 * @param[in] NoOpFunc the function used to customize the computation of the flux
 */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const
  {
    //TODO@LSL：这里删除了用于多个面的connectionIndex，因为在Flux的计算中需要考虑裂缝面单元和基质单元的关系，双重介质不需要
    // first, compute the transmissibilities at this face
    // 事实上这里应该也需要修改 传进去的不应该只有m_permeability和m_dPerm_dPres
    // m_permeability_fracture和m_dPerm_dPres_fracture也应该参与计算
    // 现在看应该不用修改，因为只需要计算基质的传导率而裂缝的传导率在窜流项的计算中几乎不起作用。
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );
    localIndex k[2];
    k[0]=0;//同一个连接关系中元素的索引，在此stencil中仅是0与1，所以不循环
    k[1]=1;
    real64 fluxVal = 0.0;
    real64 dFlux_dTrans = 0.0;
    real64 alpha = 0.0;
    real64 mobility = 0.0;
    real64 potGrad = 0.0;
    real64 trans =  stack.transmissibility ;//快速修复，实际上应该改为trans也是real64
    real64 dTrans =  stack.dTrans_dPres;
    real64 dFlux_dP[2] = {0.0, 0.0};
    localIndex const regionIndex[2]    = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
    localIndex const subRegionIndex[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
    localIndex const elementIndex[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

    geos::CrossFlowComputeKernelBase::
    computeSinglePhaseCrossFlow( regionIndex,
                                 subRegionIndex,
                                 elementIndex,
                                 trans,
                                 dTrans,
                                 m_pres,
                                 m_gravCoef,
                                 m_dens,
                                 m_dDens,
                                 m_mob,
                                 m_dMob,
                                 m_pres_fracture,
                                 m_gravCoef_fracture,
                                 m_dens_fracture,
                                 m_dDens_fracture,
                                 m_mob_fracture,
                                 m_dMob_fracture,
                                 alpha,
                                 mobility,
                                 potGrad,
                                 fluxVal,
                                 dFlux_dP,
                                 dFlux_dTrans );

    // populate local flux vector and derivatives
    stack.localFlux[k[0]*numEqn] += m_dt * fluxVal;//matrix 是正值
    stack.localFlux[k[1]*numEqn] -= m_dt * fluxVal;//fracture 是负值

    for( integer ke = 0; ke < 2; ++ke )
    {
      localIndex const localDofIndexPres = k[ke] * numDof;

      //TODO@LSL check 这里与flux中的符号相反，但是我的传导结果也是对的，不知道哪里反了
      //现在我改为与flux相同的模式试一下，为什么没有影响 QAQ ？
      stack.localFluxJacobian[k[0]*numEqn][localDofIndexPres] += m_dt * dFlux_dP[ke];
      stack.localFluxJacobian[k[1]*numEqn][localDofIndexPres] -= m_dt * dFlux_dP[ke];
    }
    kernelOp( k, regionIndex, subRegionIndex, elementIndex, alpha, mobility, potGrad, fluxVal, dFlux_dP );

  }

  /**
 * @brief Performs the complete phase for the kernel.
 * @param[in] iconn the connection index
 * @param[inout] stack the stack variables
 */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = NoOpFunc{} ) const
  {
    // add contribution to residual and jacobian into:
    // - the mass balance equation
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels

    //判断是否为边界的虚拟单元
      if( m_ghostRank[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )] < 0 )
      {
        globalIndex const globalRow = m_dofMatrixNumber[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffsetM );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow );
        GEOS_LOG("matrix ");
        GEOS_LOG("the connection id is " << iconn);
        GEOS_LOG("the global row number is " << globalRow);
        GEOS_LOG("the local row number is " << localRow);
        GEOS_LOG("the dof column indices are " << stack.dofColIndices);
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow], stack.localFlux[0 * numEqn] );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            stack.dofColIndices.data(),
                                                                            stack.localFluxJacobian[1 * numEqn].dataIfContiguous(),
                                                                            stack.stencilSize * numDof );

        // call the lambda to assemble additional terms, such as thermal terms
        kernelOp( 0, localRow );
      }
    if( m_ghostRank[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )] < 0 )
    {
      globalIndex const globalRow = m_dofFractureNumber[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffsetF );
      GEOS_ASSERT_GE( localRow, 0 );
      GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow );
      GEOS_LOG("frac ");
      GEOS_LOG("the connection id is " << iconn);
      GEOS_LOG("the global row number is " << globalRow);
      GEOS_LOG("the local row number is " << localRow);
      GEOS_LOG("the dof column indices are " << stack.dofColIndices);
      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow], stack.localFlux[1 * numEqn] );
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                          stack.dofColIndices.data(),
                                                                          stack.localFluxJacobian[0 * numEqn].dataIfContiguous(),
                                                                          stack.stencilSize * numDof );

      // call the lambda to assemble additional terms, such as thermal terms
      kernelOp( 1, localRow );
    }

  }

  /**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 * @tparam KERNEL_TYPE the kernel type
 * @param[in] numConnections the number of connections
 * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
 */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }

protected:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};


/**
 * @class FluxComputeKernelFactory
 */
class CrossFlowComputeKernelFactory
{
public:

/**
 * @brief Create a new kernel and launch
 * @tparam POLICY the policy used in the RAJA kernel
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @param[in] rankOffset the offset of my MPI rank
 * @param[in] dofKey string to get the element degrees of freedom numbers
 * @param[in] solverName name of the solver (to name accessors)
 * @param[in] elemManager reference to the element region manager
 * @param[in] stencilWrapper reference to the stencil wrapper
 * @param[in] dt time step size
 * @param[inout] localMatrix the local CRS matrix
 * @param[inout] localRhs the local right-hand side vector
 */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffsetM,
                   globalIndex const rankOffsetF,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & matrixElemManager,
                   ElementRegionManager const & fractureElemManager,
                   DualContinuumStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {

    integer constexpr NUM_EQN = 1;
    integer constexpr NUM_DOF = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofMatrixNumberAccessor =
      matrixElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofFractureNumberAccessor =
      fractureElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofMatrixNumberAccessor.setName( solverName + "/accessors/" + dofKey );
    dofFractureNumberAccessor.setName( solverName + "/accessors/" + dofKey );


    using kernelType = CrossFlowComputeKernel< NUM_EQN, NUM_DOF, DualContinuumStencilWrapper >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( matrixElemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( matrixElemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( matrixElemManager, solverName );
    typename kernelType::SinglePhaseFlowAccessors flowAccessors_fracture( fractureElemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors_fracture( fractureElemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors_fracture( fractureElemManager, solverName );

    kernelType kernel( rankOffsetM,rankOffsetF, stencilWrapper, dofMatrixNumberAccessor,
                       dofFractureNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors,
                       flowAccessors_fracture,fluidAccessors_fracture,permAccessors_fracture,
                       dt, localMatrix, localRhs );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace singlePhaseDualContinuumKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASE_CROSSFLOWCOMPUTEKERNEL_HPP
