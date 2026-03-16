
/**
 * @file FluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHASE_COMPOSITIONAL_CROSSFLOWCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHASE_COMPOSITIONAL_CROSSFLOWCOMPUTEKERNEL_HPP

#include "FluxComputeKernelBase.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/KernelLaunchSelectors.hpp"
#include "PPUPhaseFlux.hpp"
#include "C1PPUPhaseFlux.hpp"
#include "IHUPhaseFlux.hpp"
#include "HU2PhaseFlux.hpp"
#include "PhaseComponentFlux.hpp"

namespace geos
{

namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernels
{

/**
 * @class FluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public FluxComputeKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations (all of them, except the volume balance equation)
  static constexpr integer numEqn = NUM_DOF-1;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /// Number of flux support points (hard-coded for TFPA)
  static constexpr integer numFluxSupportPoints = 2;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] m_rankOffset_m the offset of my MPI rank for matrix
   * @param[in] m_rankOffset_f the offset of my MPI rank for fracture
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor_m accessor for the dof numbers of matrix
   * @param[in] compFlowAccessors_m accessor for wrappers registered by the solver for matrix
   * @param[in] multiFluidAccessors_m accessor for wrappers registered by the multifluid model for matrix
   * @param[in] dofNumberAccessor_f accessor for the dof numbers of fracture
   * @param[in] compFlowAccessors_f accessor for wrappers registered by the solver for fracture
   * @param[in] multiFluidAccessors_f accessor for wrappers registered by the multifluid model for fracture
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  FluxComputeKernel( integer const numPhases,
                     globalIndex const rankOffset_m,
                     globalIndex const rankOffset_f,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor_m,
                     CompFlowAccessors const & compFlowAccessors_m,
                     MultiFluidAccessors const & multiFluidAccessors_m,
                     DofNumberAccessor const & dofNumberAccessor_f,
                     CompFlowAccessors const & compFlowAccessors_f,
                     MultiFluidAccessors const & multiFluidAccessors_f,
                     CapPressureAccessors const & capPressureAccessors_m,
                     CapPressureAccessors const & capPressureAccessors_f,
                     PermeabilityAccessors const & permeabilityAccessors_m,
                     PermeabilityAccessors const & permeabilityAccessors_f,
                     GravityDrainagePressureAccessors const & gravityDrainagePressureAccessors_m,
                     real64 const dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs,
                     BitFlags< KernelFlags > kernelFlags )
    : FluxComputeKernelBase( numPhases,
                             rankOffset_m,
                             rankOffset_f,
                             dofNumberAccessor_m,
                             compFlowAccessors_m,
                             multiFluidAccessors_m,
                             dofNumberAccessor_f,
                             compFlowAccessors_f,
                             multiFluidAccessors_f,
                             gravityDrainagePressureAccessors_m,
                             dt,
                             localMatrix,
                             localRhs,
                             kernelFlags ),
    m_permeability_m( permeabilityAccessors_m.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres_m( permeabilityAccessors_m.get( fields::permeability::dPerm_dPressure {} ) ),
    m_permeability_f( permeabilityAccessors_f.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres_f( permeabilityAccessors_f.get( fields::permeability::dPerm_dPressure {} ) ),
    m_phaseMob_m( compFlowAccessors_m.get( fields::flow::phaseMobility {} ) ),
    m_dPhaseMob_m( compFlowAccessors_m.get( fields::flow::dPhaseMobility {} ) ),
    m_phaseMob_f( compFlowAccessors_f.get( fields::flow::phaseMobility {} ) ),
    m_dPhaseMob_f( compFlowAccessors_f.get( fields::flow::dPhaseMobility {} ) ),
    m_phaseMassDens_m( multiFluidAccessors_m.get( fields::multifluid::phaseMassDensity {} ) ),
    m_dPhaseMassDens_m( multiFluidAccessors_m.get( fields::multifluid::dPhaseMassDensity {} ) ),
    m_phaseMassDens_f( multiFluidAccessors_f.get( fields::multifluid::phaseMassDensity {} ) ),
    m_dPhaseMassDens_f( multiFluidAccessors_f.get( fields::multifluid::dPhaseMassDensity {} ) ),
    m_phaseCapPressure_m( capPressureAccessors_m.get( fields::cappres::phaseCapPressure {} ) ),
    m_dPhaseCapPressure_dPhaseVolFrac_m( capPressureAccessors_m.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
    m_phaseCapPressure_f( capPressureAccessors_f.get( fields::cappres::phaseCapPressure {} ) ),
    m_dPhaseCapPressure_dPhaseVolFrac_f( capPressureAccessors_f.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
    m_gravityDrainagePressure_m( gravityDrainagePressureAccessors_m.get( fields::gravdrainage::gravityDrainagePressure {} ) ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  { }

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
      numConnectedElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;
    /// Number of elements connected at a given connection
    localIndex const numConnectedElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres{};

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
  inline
  localIndex stencilSize( localIndex const iconn ) const { return m_sei[iconn].size(); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  inline
  localIndex numPointsInFlux( localIndex const iconn ) const { return m_stencilWrapper.numPointsInFlux( iconn ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face

      globalIndex const offset_m = m_dofNumber_m[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )];
      globalIndex const offset_f = m_dofNumber_f[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )];
      //分别对matrix和fracture的dof进行赋值

    for( integer jdof = 0; jdof < numDof; ++jdof )
    {
      stack.dofColIndices[ jdof ] = offset_m + jdof;
      stack.dofColIndices[numDof + jdof] = offset_f + jdof;

    }

  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    using namespace isothermalDualContinuumCompositionalMultiPhaseCrossFlowKernelUtilities;

    // first, compute the transmissibilities at this face
    // 在这里仅使用matrix的permeability来计算transmissibility
    // 因为在计算窜流相的过程中fracture几乎不对传导率有贡献
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability_m,
                                     m_dPerm_dPres_m,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );


    localIndex k[numFluxSupportPoints];
    localIndex connectionIndex = 0;
    for( k[0] = 0; k[0] < stack.numConnectedElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numConnectedElems; ++k[1] )
      {
        /// cell indices
        localIndex const seri[numFluxSupportPoints]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[numFluxSupportPoints] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[numFluxSupportPoints]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // clear working arrays
        real64 compFlux[numComp]{};
        real64 dCompFlux_dP[numFluxSupportPoints][numComp]{};
        real64 dCompFlux_dC[numFluxSupportPoints][numComp][numComp]{};
        real64 dCompFlux_dTrans[numComp]{};

        real64 const trans[numFluxSupportPoints] = { stack.transmissibility,
                                                     stack.transmissibility };

        real64 const dTrans_dPres[numFluxSupportPoints] = { stack.dTrans_dPres,
                                                            stack.dTrans_dPres };

        //***** calculation of flux *****
        // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // create local work arrays
          real64 potGrad = 0.0;
          real64 phaseFlux = 0.0;
          real64 dPhaseFlux_dP[numFluxSupportPoints]{};
          real64 dPhaseFlux_dC[numFluxSupportPoints][numComp]{};
          real64 dPhaseFlux_dTrans = 0.0; // not really used

          if( m_kernelFlags.isSet( KernelFlags::C1PPU ) )
          {
            C1PPUPhaseFlux::compute< numComp, numFluxSupportPoints >
              ( m_numPhases,
              ip,
              m_kernelFlags.isSet( KernelFlags::CapPressure ),
              m_kernelFlags.isSet( KernelFlags::CheckPhasePresenceInGravity ),
              m_kernelFlags.isSet( KernelFlags::GravityDrainage),
              seri, sesri, sei,
              trans,
              dTrans_dPres,
              m_pres_m, m_pres_f, // Using matrix pressure for first support point, fracture for second
              m_gravCoef_m, m_gravCoef_f, // Using matrix gravity coefficient for first support point, fracture for second
              m_phaseMob_m, m_phaseMob_f, m_dPhaseMob_m, m_dPhaseMob_f, // Using matrix phase mobilities for first support point, fracture for second
              m_phaseVolFrac_m, m_phaseVolFrac_f, m_dPhaseVolFrac_m, m_dPhaseVolFrac_f, // Using matrix phase volume fractions for first support point, fracture for second
              m_dCompFrac_dCompDens_m, m_dCompFrac_dCompDens_f, // Using matrix comp frac derivatives for first support point, fracture for second
              m_phaseMassDens_m, m_phaseMassDens_f, m_dPhaseMassDens_m, m_dPhaseMassDens_f, // Using matrix phase mass densities for first support point, fracture for second
              m_phaseCapPressure_m, m_phaseCapPressure_f, m_dPhaseCapPressure_dPhaseVolFrac_m, m_dPhaseCapPressure_dPhaseVolFrac_f, // Using matrix cap pressure for first support point, fracture for second
              m_gravityDrainagePressure,
              potGrad,
              phaseFlux,
              dPhaseFlux_dP,
              dPhaseFlux_dC );
          }
          else if( m_kernelFlags.isSet( KernelFlags::IHU ) )
          {
            IHUPhaseFlux::compute< numComp, numFluxSupportPoints >
              ( m_numPhases,
              ip,
              m_kernelFlags.isSet( KernelFlags::CapPressure ),
              m_kernelFlags.isSet( KernelFlags::CheckPhasePresenceInGravity ),
              m_kernelFlags.isSet( KernelFlags::GravityDrainage),
              seri, sesri, sei,
              trans,
              dTrans_dPres,
              m_pres_m, m_pres_f, // Using matrix pressure for first support point, fracture for second
              m_gravCoef_m, m_gravCoef_f, // Using matrix gravity coefficient for first support point, fracture for second
              m_phaseMob_m, m_phaseMob_f, m_dPhaseMob_m, m_dPhaseMob_f, // Using matrix phase mobilities for first support point, fracture for second
              m_phaseVolFrac_m, m_phaseVolFrac_f, m_dPhaseVolFrac_m, m_dPhaseVolFrac_f, // Using matrix phase volume fractions for first support point, fracture for second
              m_dCompFrac_dCompDens_m, m_dCompFrac_dCompDens_f, // Using matrix comp frac derivatives for first support point, fracture for second
              m_phaseMassDens_m, m_phaseMassDens_f, m_dPhaseMassDens_m, m_dPhaseMassDens_f, // Using matrix phase mass densities for first support point, fracture for second
              m_phaseCapPressure_m, m_phaseCapPressure_f, m_dPhaseCapPressure_dPhaseVolFrac_m, m_dPhaseCapPressure_dPhaseVolFrac_f, // Using matrix cap pressure for first support point, fracture for second
              m_gravityDrainagePressure_m,
              potGrad,
              phaseFlux,
              dPhaseFlux_dP,
              dPhaseFlux_dC );
          }
          else if( m_kernelFlags.isSet( KernelFlags::HU2PH ) )
          {
            HU2PhaseFlux::compute< numComp, numFluxSupportPoints >
              ( m_numPhases,
              ip,
              m_kernelFlags.isSet( KernelFlags::CapPressure ),
              m_kernelFlags.isSet( KernelFlags::CheckPhasePresenceInGravity ),
              m_kernelFlags.isSet( KernelFlags::GravityDrainage),
              seri, sesri, sei,
              trans,
              dTrans_dPres,
              m_pres_m, m_pres_f, // Using matrix pressure for first support point, fracture for second
              m_gravCoef_m, m_gravCoef_f, // Using matrix gravity coefficient for first support point, fracture for second
              m_phaseMob_m, m_phaseMob_f, m_dPhaseMob_m, m_dPhaseMob_f, // Using matrix phase mobilities for first support point, fracture for second
              m_phaseVolFrac_m, m_phaseVolFrac_f, m_dPhaseVolFrac_m, m_dPhaseVolFrac_f, // Using matrix phase volume fractions for first support point, fracture for second
              m_dCompFrac_dCompDens_m, m_dCompFrac_dCompDens_f, // Using matrix comp frac derivatives for first support point, fracture for second
              m_phaseMassDens_m, m_phaseMassDens_f, m_dPhaseMassDens_m, m_dPhaseMassDens_f, // Using matrix phase mass densities for first support point, fracture for second
              m_phaseCapPressure_m, m_phaseCapPressure_f, m_dPhaseCapPressure_dPhaseVolFrac_m, m_dPhaseCapPressure_dPhaseVolFrac_f, // Using matrix cap pressure for first support point, fracture for second
              m_gravityDrainagePressure,
              potGrad,
              phaseFlux,
              dPhaseFlux_dP,
              dPhaseFlux_dC );
          }
          else
          {
            PPUPhaseFlux::compute< numComp, numFluxSupportPoints >
              ( m_numPhases,
              ip,
              m_kernelFlags.isSet( KernelFlags::CapPressure ),
              m_kernelFlags.isSet( KernelFlags::CheckPhasePresenceInGravity ),
              m_kernelFlags.isSet( KernelFlags::GravityDrainage),
              seri, sesri, sei,
              trans,
              dTrans_dPres,
              m_pres_m, m_pres_f, // Using matrix pressure for first support point, fracture for second
              m_gravCoef_m, m_gravCoef_f, // Using matrix gravity coefficient for first support point, fracture for second
              m_phaseMob_m, m_phaseMob_f, m_dPhaseMob_m, m_dPhaseMob_f, // Using matrix phase mobilities for first support point, fracture for second
              m_phaseVolFrac_m, m_phaseVolFrac_f, m_dPhaseVolFrac_m, m_dPhaseVolFrac_f, // Using matrix phase volume fractions for first support point, fracture for second
              m_dCompFrac_dCompDens_m, m_dCompFrac_dCompDens_f, // Using matrix comp frac derivatives for first support point, fracture for second
              m_phaseMassDens_m, m_phaseMassDens_f, m_dPhaseMassDens_m, m_dPhaseMassDens_f, // Using matrix phase mass densities for first support point, fracture for second
              m_phaseCapPressure_m, m_phaseCapPressure_f, m_dPhaseCapPressure_dPhaseVolFrac_m, m_dPhaseCapPressure_dPhaseVolFrac_f, // Using matrix cap pressure for first support point, fracture for second
              m_gravityDrainagePressure,
              potGrad,
              phaseFlux,
              dPhaseFlux_dP,
              dPhaseFlux_dC,
              dPhaseFlux_dTrans );
          }

          // choose upstream cell for composition upwinding
          localIndex k_up;

          // distribute on phaseComponentFlux here
          //TODO@LSL: 计算相中各组分的通量，并计算其对压力、组分浓度和传导率的导数。这里需要根据phaseFlux的符号来确定上游单元是k[0]还是k[1]，并使用对应的compFrac和其导数进行计算。

          if(phaseFlux >= 0)
          {
            k_up = 0;
            PhaseComponentFlux::compute( ip, k_up, seri, sesri, sei,
                                         m_phaseCompFrac_m, m_dPhaseCompFrac_m, m_dCompFrac_dCompDens_m, // Using matrix comp fractions
                                         phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, dPhaseFlux_dTrans,
                                         compFlux, dCompFlux_dP, dCompFlux_dC, dCompFlux_dTrans );
          }
          else
          {
            k_up = 1;
            PhaseComponentFlux::compute( ip, k_up, seri, sesri, sei,
                                         m_phaseCompFrac_f, m_dPhaseCompFrac_f, m_dCompFrac_dCompDens_f, // Using fracture comp fractions
                                         phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, dPhaseFlux_dTrans,
                                         compFlux, dCompFlux_dP, dCompFlux_dC, dCompFlux_dTrans );
          }


          // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
          // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
          // 钩子函数，用于在计算完phaseComponentFlux后，进行一些额外的计算
          compFluxKernelOp( ip, m_kernelFlags.isSet( KernelFlags::CheckPhasePresenceInGravity ),
                            k, seri, sesri, sei, connectionIndex,
                            k_up, seri[k_up], sesri[k_up], sei[k_up], potGrad,
                            phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

        }   // loop over phases

        /// populate local flux vector and derivatives

//                节点0压力  节点0组分0  节点0组分1  节点1压力  节点1组分0  节点1组分1
// 节点0压力 已组装   ∂R0/∂P0   ∂R0/∂C0,0  ∂R0/∂C0,1  ∂R0/∂P1   ∂R0/∂C1,0  ∂R0/∂C1,1
// 节点0组分0        ∂R0/∂P0   ∂R0/∂C0,0  ∂R0/∂C0,1  ∂R0/∂P1   ∂R0/∂C1,0  ∂R0/∂C1,1
// 节点0组分1        ∂R1/∂P0   ∂R1/∂C0,0  ∂R1/∂C0,1  ∂R1/∂P1   ∂R1/∂C1,0  ∂R1/∂C1,1
// 节点1压力 已组装   ∂R2/∂P0   ∂R2/∂C0,0  ∂R2/∂C0,1  ∂R2/∂P1   ∂R2/∂C1,0  ∂R2/∂C1,1
// 节点1组分0        ∂R2/∂P0   ∂R2/∂C0,0  ∂R2/∂C0,1  ∂R2/∂P1   ∂R2/∂C1,0  ∂R2/∂C1,1
// 节点1组分1        ∂R3/∂P0   ∂R3/∂C0,0  ∂R3/∂C0,1  ∂R3/∂P1   ∂R3/∂C1,0  ∂R3/∂C1,1

        for( integer ic = 0; ic < numComp; ++ic )
        {
          integer const eqIndex0 = k[0] * numEqn + ic;
          integer const eqIndex1 = k[1] * numEqn + ic;

          stack.localFlux[eqIndex0]  +=  m_dt * compFlux[ic];
          stack.localFlux[eqIndex1]  -=  m_dt * compFlux[ic];

          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            localIndex const localDofIndexPres = k[ke] * numDof;
            stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dCompFlux_dP[ke][ic];
            stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dCompFlux_dP[ke][ic];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
              stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dCompFlux_dC[ke][ic][jc];
              stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dCompFlux_dC[ke][ic][jc];
            }
          }
        }
        connectionIndex++;
      }   // loop over k[1]
    }   // loop over k[0]

  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;
      // TODO@LSL: Need to determine if this is matrix or fracture and use the appropriate accessor

    if( m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // Apply equation/variable change transformation(s)
      stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numEqn, numDof * stack.stencilSize, stack.numConnectedElems,
                                                               stack.localFluxJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, numEqn, stack.numConnectedElems,
                                                                 stack.localFlux );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels

      // TODO@LSL: Need to determine if this is matrix or fracture and use the appropriate accessor
      // For now, using matrix accessor as default
      if( m_ghostRank_m[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber_m[m_seri( iconn, 0 )][m_sesri( iconn, 0 )][m_sei( iconn, 0 )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset_m );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic],
                           stack.localFlux[ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        assemblyKernelOp( 0, localRow );
      }
      if( m_ghostRank_f[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber_f[m_seri( iconn, 1 )][m_sesri( iconn, 1 )][m_sei( iconn, 1 )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset_f );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic],
                           stack.localFlux[numEqn + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
                         ( localRow + ic,
                           stack.dofColIndices.data(),
                           stack.localFluxJacobian[ numEqn + ic].dataIfContiguous(),
                           stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        assemblyKernelOp( 1, localRow );
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

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > const m_permeability_m;
  ElementViewConst< arrayView3d< real64 const > > const m_dPerm_dPres_m;
  ElementViewConst< arrayView3d< real64 const > > const m_permeability_f;
  ElementViewConst< arrayView3d< real64 const > > const m_dPerm_dPres_f;

  /// Views on phase mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob_m;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob_m;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob_f;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob_f;

  /// Views on phase mass densities
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseMassDens_m;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseMassDens_m;
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseMassDens_f;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseMassDens_f;

  /// Views on phase capillary pressure
  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const m_phaseCapPressure_m;
  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac_m;
  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const m_phaseCapPressure_f;
  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac_f;

  /// 重力项
  ElementViewConst< arrayView2d< real64 const > > const m_gravityDrainagePressure_m;

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
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] m_rankOffset_m the offset of my MPI rank for matrix
   * @param[in] m_rankOffset_f the offset of my MPI rank for fracture
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] kernelFlags flags packed together
   * @param[in] primarySolverName name of the primary solver (to name accessors)
   * @param[in] secondarySolverName name of the secondary solver (to name accessors)
   * @param[in] matrixElemManager reference to the matrix element region manager
   * @param[in] fractureElemManager reference to the fracture element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const m_rankOffset_m,
                   globalIndex const m_rankOffset_f,
                   string const & dofKey,
                   BitFlags< KernelFlags > kernelFlags,
                   string const & primarySolverName,
                   string const & secondarySolverName,
                   ElementRegionManager const & matrixElemManager,
                   ElementRegionManager const & fractureElemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor_m =
        matrixElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor_m.setName( primarySolverName + "/accessors/" + dofKey + "_m" );

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor_f =
        fractureElemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor_f.setName( secondarySolverName + "/accessors/" + dofKey + "_f" );

      using kernelType = FluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename kernelType::CompFlowAccessors compFlowAccessors_m( matrixElemManager, primarySolverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors_m( matrixElemManager, primarySolverName );
      typename kernelType::CompFlowAccessors compFlowAccessors_f( fractureElemManager, secondarySolverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors_f( fractureElemManager, secondarySolverName );
      typename kernelType::CapPressureAccessors capPressureAccessors_m( matrixElemManager, primarySolverName );
      typename kernelType::CapPressureAccessors capPressureAccessors_f( fractureElemManager, secondarySolverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors_m( matrixElemManager, primarySolverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors_f( fractureElemManager, secondarySolverName );
      typename kernelType::GravityDrainagePressureAccessors gravityDrainagePressureAccessors_m(matrixElemManager,primarySolverName);

      kernelType kernel( numPhases, m_rankOffset_m, m_rankOffset_f, stencilWrapper, dofNumberAccessor_m,
                         compFlowAccessors_m, multiFluidAccessors_m, dofNumberAccessor_f, compFlowAccessors_f, multiFluidAccessors_f,
                         capPressureAccessors_m, capPressureAccessors_f, permeabilityAccessors_m, permeabilityAccessors_f,
                         gravityDrainagePressureAccessors_m,
                         dt, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUXCOMPUTEKERNEL_HPP
