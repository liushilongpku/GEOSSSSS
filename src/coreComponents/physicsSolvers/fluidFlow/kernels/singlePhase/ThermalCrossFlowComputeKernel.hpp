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
 * @file ThermalFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/multiphysics/dualContinuumCrossFlowComputeKernels/CrossFlowComputeKernel.hpp"

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace singlePhaseThermalDualContinuumKernels
{
/******************************** FluxComputeKernel ********************************/

/**
 * @class FluxComputeKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class CrossFlowComputeKernel : public singlePhaseDualContinuumKernels::CrossFlowComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = singlePhaseDualContinuumKernels::CrossFlowComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;

  using AbstractBase::m_rankOffsetM;
  using AbstractBase::m_dofMatrixNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens;

  using AbstractBase::m_rankOffsetF;
  using AbstractBase::m_dofFractureNumber;
  using AbstractBase::m_gravCoef_fracture;
  using AbstractBase::m_mob_fracture;
  using AbstractBase::m_dMob_fracture;
  using AbstractBase::m_dens_fracture;
  using AbstractBase::m_dDens_fracture;



  using Base = singlePhaseDualContinuumKernels::CrossFlowComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity,
                              fields::thermalconductivity::dEffectiveConductivity_dT >;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] singlePhaseFlowAccessors accessor for wrappers registered by the solver
   * @param[in] thermalSinglePhaseFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the single fluid model
   * @param[in] thermalSinglePhaseFluidAccessors accessor for *thermal* wrappers registered by the single fluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] thermalConductivityAccessors accessor for wrappers registered by the thermal conductivity model
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
                          ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                          ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                          PermeabilityAccessors const & permeabilityAccessors,
                          ThermalConductivityAccessors const & thermalConductivityAccessors,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors_fracture,
                          ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors_fracture,
                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors_fracture,
                          ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors_fracture,
                          PermeabilityAccessors const & permeabilityAccessors_fracture,
                          ThermalConductivityAccessors const & thermalConductivityAccessors_fracture,
                          real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffsetM,
            rankOffsetF,
            stencilWrapper,
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
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_temp_fracture( thermalSinglePhaseFlowAccessors_fracture.get( fields::flow::temperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy {} ) ),
    m_enthalpy_fracture( thermalSinglePhaseFluidAccessors_fracture.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy_fracture( thermalSinglePhaseFluidAccessors_fracture.get( fields::singlefluid::dEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_thermalConductivity_fracture( thermalConductivityAccessors_fracture.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dThermalCond_dT( thermalConductivityAccessors.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) ),
    m_dThermalCond_dT_fracture( thermalConductivityAccessors_fracture.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      energyFlux( 0.0 ),
      dEnergyFlux_dP( size ),
      dEnergyFlux_dT( size )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    // Thermal transmissibility
    //real64 thermalTransmissibility[maxNumConns][2]{};
    real64 thermalTransmissibility{};//初始化为0.0，裂缝和基质共用一个传导率，后续可以再区分
    //裂缝的传导是否单独享有一个传导系数？

    /// Derivatives of thermal transmissibility with respect to temperature
    //real64 dThermalTrans_dT[maxNumConns][2]{};//同样的
    real64 dThermalTrans_dT{};
    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux;
    /// Derivatives of energy fluxes wrt pressure
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dP;
    /// Derivatives of energy fluxes wrt temperature
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dT;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of energyFlux  and its derivatives (including derivatives wrt temperature)
    //
    // Computing dFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, fluxVal, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                                 localIndex const (&seri)[2],
                                                 localIndex const (&sesri)[2],
                                                 localIndex const (&sei)[2],
                                                 real64 const alpha,
                                                 real64 const mobility,
                                                 real64 const & potGrad,
                                           real64 const & fluxVal,
                                           real64 const (&dFlux_dP)[2] )
    {
      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      real64 dDensMean_dT[2]{0.0, 0.0};

      real64 const trans[2] = { stack.transmissibility, stack.transmissibility};//TODO@LSL传导率应该不一样，在此优先简化保证正确


      dDensMean_dT[0] = 0.5 * m_dDens[seri[0]][sesri[0]][sei[0]][0][DerivOffset::dT];
      dDensMean_dT[1] = 0.5 * m_dDens_fracture[seri[1]][sesri[1]][sei[1]][0][DerivOffset::dT];

      // Step 2: compute the derivatives of the potential difference wrt temperature
      //***** calculation of flux *****

      real64 dGravHead_dT[2]{0.0, 0.0};

      // compute potential difference

      // compute derivative of gravity potential difference wrt temperature
      real64 const gravD_matrix = trans[0] * m_gravCoef[seri[0]][sesri[0]][sei[0]];
      real64 const gravD_fracture = trans[1] * m_gravCoef_fracture[seri[1]][sesri[1]][sei[1]];

      for( integer i = 0; i < 2; ++i )
      {
        dGravHead_dT[i] += dDensMean_dT[i] *(gravD_matrix + gravD_fracture);
      }


      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      real64 dFlux_dT[2]{0.0, 0.0};

      // Step 3.1: compute the derivative of flux wrt temperature
      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] -= dGravHead_dT[ke];
      }

      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] *= mobility;
      }

      real64 dMob_dT[2]{};

      if( alpha <= 0.0  )
      {
        //localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );//根据alpha确定上游
        dMob_dT[0] = m_dMob[seri[0]][sesri[0]][sei[0]][DerivOffset::dT];
      }
      else if(alpha >= 1.0)
      {
        dMob_dT[1] = m_dMob_fracture[seri[1]][sesri[1]][sei[1]][DerivOffset::dT];
      }
      else
      {
        dMob_dT[0] = alpha * m_dMob[seri[0]][sesri[0]][sei[0]][DerivOffset::dT];
        dMob_dT[1] = (1.0 - alpha )* m_dMob_fracture[seri[1]][sesri[1]][sei[1]][DerivOffset::dT];
      }

      // add contribution from upstream cell mobility derivatives
      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] += dMob_dT[ke] * potGrad;//potGrad 未经验证
      }

      // add dFlux_dTemp to localFluxJacobian 这里貌似只添加了流量对温度的导数，所以没有k[0]*numEqn坐标没有改变，而只是改变了第二个索引
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - 1;
        stack.localFluxJacobian[k[0]*numEqn][localDofIndexTemp] += m_dt * dFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn][localDofIndexTemp] -= m_dt * dFlux_dT[ke];
      }

      // Step 4: compute the enthalpy flux
      real64 enthalpy = 0.0;
      real64 dEnthalpy_dP[2]{0.0, 0.0};
      real64 dEnthalpy_dT[2]{0.0, 0.0};

      //TODO@LSL 后续需要检查，目前理解这里的alpha是0（基质）贡献的占比，因此当alpha<0时，全部由裂缝
      if(  alpha <= 0.0 )
      {
        localIndex const k_up = 1;

        enthalpy = m_enthalpy_fracture[seri[k_up]][sesri[k_up]][sei[k_up]][0];//裂缝介质
        dEnthalpy_dP[k_up] = m_dEnthalpy_fracture[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dP];
        dEnthalpy_dT[k_up] = m_dEnthalpy_fracture[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dT];
      }
      else if( alpha >= 1.0)
      {
        localIndex const k_up = 0;

        enthalpy = m_enthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0];//基质介质
        dEnthalpy_dP[k_up] = m_dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dP];
        dEnthalpy_dT[k_up] = m_dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dT];
      }
      else
      {
        real64 const mobWeights[2] = { alpha, 1.0 - alpha };

        enthalpy = mobWeights[0] * m_enthalpy[seri[0]][sesri[0]][sei[0]][0] + mobWeights[1] * m_enthalpy_fracture[seri[1]][sesri[1]][sei[1]][0];
        dEnthalpy_dP[0] = mobWeights[0] * m_dEnthalpy[seri[0]][sesri[0]][sei[0]][0][DerivOffset::dP];
        dEnthalpy_dP[1] = mobWeights[1] * m_dEnthalpy_fracture[seri[1]][sesri[1]][sei[1]][0][DerivOffset::dP];
        dEnthalpy_dT[0] = mobWeights[0] * m_dEnthalpy[seri[0]][sesri[0]][sei[0]][0][DerivOffset::dT];
        dEnthalpy_dT[1] = mobWeights[1] * m_dEnthalpy_fracture[seri[1]][sesri[1]][sei[1]][0][DerivOffset::dT];
        
      }

      stack.energyFlux += fluxVal * enthalpy;

      for( integer ke = 0; ke < 2; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += dFlux_dP[ke] * enthalpy;
        stack.dEnergyFlux_dT[ke] += dFlux_dT[ke] * enthalpy;
      }

      for( integer ke = 0; ke < 2; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += fluxVal * dEnthalpy_dP[ke];
        stack.dEnergyFlux_dT[ke] += fluxVal * dEnthalpy_dT[ke];
      }

    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    // We follow how the thermal compositional multi-phase solver does to update the thermal transmissibility
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_dThermalCond_dT,
                                     stack.thermalTransmissibility,
                                     stack.dThermalTrans_dT );

    localIndex k[2];
    k[0]=0;
    k[1]=1;
    //TODO@LSL 第一个是基质传导系数，第二个是裂缝传导系数，裂缝实际上不拥有固态的传导率，应该使用液相的等效传导系数去替代
    //在这里还没有做，暂时使用相同的传导系数，后续需要完善一下
    //在基质中，应该是固相和液相两种的综合传导率，check是否是这样
    //此外还需要检查一下传导系数是否已经包含了传导方向，这对后边计算energyflux有重要作用。
    real64 const thermalTrans[2] = { stack.thermalTransmissibility, stack.thermalTransmissibility };
    real64 const dThermalTrans_dT[2] = { stack.dThermalTrans_dT, stack.dThermalTrans_dT };

    localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
    localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
    localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

    // Step 2: compute temperature difference at the interface
    for( integer ke = 0; ke < 2; ++ke )
    {
      localIndex const er  = seri[ke];
      localIndex const esr = sesri[ke];
      localIndex const ei  = sei[ke];
    }
    //这里采用+ 的原因可能是因为传导系数thermalTrans已经包含了流动方向的影响
    //这里仍然需要保留+=，以包含前面的热对流
    //测试一下第二个裂缝部分的值减小，懵一个
    //TODO@LSL对了，那应该是原来的thermalTrans中含有方向，但是现在我直接抓来的数据没有，后续需要测试一下
    stack.energyFlux += (thermalTrans[0] * m_temp[seri[0]][sesri[0]][sei[0]] - thermalTrans[1] * m_temp_fracture[seri[1]][sesri[1]][sei[1]]);
    stack.dEnergyFlux_dT[0] += thermalTrans[0] + dThermalTrans_dT[0] * m_temp[seri[0]][sesri[0]][sei[0]];
    stack.dEnergyFlux_dT[1] += thermalTrans[1] + dThermalTrans_dT[1] * m_temp_fracture[seri[1]][sesri[1]][sei[1]];
    

    // add energyFlux and its derivatives to localFlux and localFluxJacobian
    stack.localFlux[k[0]*numEqn + numEqn - 1] += m_dt * stack.energyFlux;
    stack.localFlux[k[1]*numEqn + numEqn - 1] -= m_dt * stack.energyFlux;

    for( integer ke = 0; ke < 2; ++ke )
    {
      integer const localDofIndexPres = k[ke] * numDof;
      stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexPres] =  m_dt * stack.dEnergyFlux_dP[ke];
      stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[ke];
      integer const localDofIndexTemp = localDofIndexPres + numDof - 1;
      stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexTemp] =  m_dt * stack.dEnergyFlux_dT[ke];
      stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexTemp] = -m_dt * stack.dEnergyFlux_dT[ke];
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Call Case::complete to assemble the mass balance equations
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      // Different from the one in compositional multi-phase flow, which has a volume balance eqn.
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn-1], stack.localFlux[i * numEqn + numEqn-1] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn-1,
                                                                                        stack.dofColIndices.data(),
                                                                                        stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
                                                                                        stack.stencilSize * numDof );

    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;
  ElementViewConst< arrayView1d< real64 const > > const m_temp_fracture;


  /// Views on enthalpies
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_enthalpy;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dEnthalpy;
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_enthalpy_fracture;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dEnthalpy_fracture;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity_fracture;

  /// View on derivatives of thermal conductivity w.r.t. temperature
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT;
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT_fracture;

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
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( globalIndex const rankOffsetM,
                   globalIndex const rankOffsetF,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & matrixElemManager,
                   ElementRegionManager const & fractureElemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;
    integer constexpr NUM_EQN = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      matrixElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor_fracture =
      fractureElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using KernelType = CrossFlowComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename KernelType::SinglePhaseFlowAccessors flowAccessors( matrixElemManager, solverName );
    typename KernelType::SinglePhaseFlowAccessors flowAccessors_fracture( fractureElemManager, solverName );

    typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( matrixElemManager, solverName );
    typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors_fracture( fractureElemManager, solverName );

    typename KernelType::SinglePhaseFluidAccessors fluidAccessors( matrixElemManager, solverName );
    typename KernelType::SinglePhaseFluidAccessors fluidAccessors_fracture( fractureElemManager, solverName );

    typename KernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( matrixElemManager, solverName );
    typename KernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors_fracture( fractureElemManager, solverName );

    typename KernelType::PermeabilityAccessors permAccessors( matrixElemManager, solverName );
    typename KernelType::PermeabilityAccessors permAccessors_fracture( fractureElemManager, solverName );

    typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( matrixElemManager, solverName );
    typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors_fracture( fractureElemManager, solverName );

    KernelType kernel( rankOffsetM,rankOffsetF, stencilWrapper, dofNumberAccessor,dofNumberAccessor_fracture,
                       flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                       permAccessors, thermalConductivityAccessors,
                       flowAccessors_fracture, thermalFlowAccessors_fracture, fluidAccessors_fracture, thermalFluidAccessors_fracture,
                       permAccessors_fracture,thermalConductivityAccessors_fracture,
                       dt, localMatrix, localRhs );
    KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace thermalSinglePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALFLUXCOMPUTEKERNEL_HPP
