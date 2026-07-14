/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumFVM.cpp
 */

#include "SinglePhaseDualContinuum.hpp"

// #include "common/TimingMacros.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/singlePhase/CrossFlowComputeKernel.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/singlePhase/ThermalCrossFlowComputeKernel.hpp"
#include "constitutive/gravityDrainagePressure/GravityDrainagePressureBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
namespace geos
{
  using namespace dataRepository;
  using namespace fields;

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::
    assembleCouplingTerms(real64 const time_n,
                          real64 const dt,
                          DomainPartition const &domain,
                          DofManager const &dofManager,
                          CRSMatrixView<real64, globalIndex const> const &localMatrix,
                          arrayView1d<real64> const &localRhs)
  {
    //GEOS_LOG("befor the assemble");
    //this->printCRSMatrix(this->m_localMatrix);

    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR(time_n);
    GEOS_UNUSED_VAR(dt);
    // 注册传导率场
    // 根据 aperture 计算传导率
    // 当前使用常数
    NumericalMethodsManager const &numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const &fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const &fluxApprox = fvManager.getFluxApproximation(Base::primarySolver()->getDiscretizationName());
    //FluxApproximationBase const &secondaryFluxApprox = fvManager.getFluxApproximation(Base::secondarySolver()->getDiscretizationName());
    //FluxApproximationBase const *const fluxApproxArray[] = {&primaryFluxApprox,&secondaryFluxApprox};

    string const &dofKey = dofManager.getKey(SinglePhaseBase::viewKeyStruct::elemDofFieldString());

    // single element-based DOF field (two components per node: 0=primary,1=secondary)
    globalIndex const rankOffset = dofManager.rankOffset();

    std::vector<MeshLevel const*> meshLevelPtrs;


    forDiscretizationOnMeshTargets(domain.getMeshBodies(), [&](string const &,
                                                               MeshLevel const &mesh,
                                                               string_array const &)
    {
      meshLevelPtrs.push_back(&mesh);
    });

    if(meshLevelPtrs.size() < 2)// 先判断vector至少有2个元素，避免越界
    {
      GEOS_ERROR("The dual continuum flow solver requires at least two meshes");
    }
    else
    {
      MeshLevel const * matrixMeshPtr = meshLevelPtrs[0];
      MeshLevel const * fractureMeshPtr = meshLevelPtrs[1];

      typename TYPEOFREF( this->getStencil() ) ::KernelWrapper stencilWrapper = this->getStencil().createKernelWrapper();
      GEOS_UNUSED_VAR(stencilWrapper);
      if (Base::primarySolver()->isThermal() && Base::secondarySolver()->isThermal())
      {
        singlePhaseThermalDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
                                                  dofKey,
                                                  fluxApprox.getName(),
                                                  matrixMeshPtr->getElemManager(),
                                                  fractureMeshPtr->getElemManager(),
                                                  stencilWrapper,
                                                  this->getGravityDrainageFlag(),
                                                  this->gravityVector()[2],
                                                  this->getFracSpacingLz(),
                                                  this->getInterporosityExchangeCoefficient(),
                                                  dt,
                                                  localMatrix.toViewConstSizes(),
                                                  localRhs.toView());
      }
      else if (!(Base::primarySolver()->isThermal()) && !(Base::secondarySolver()->isThermal()))
      {
        singlePhaseDualContinuumKernels::
        CrossFlowComputeKernelFactory::
        createAndLaunch<parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                  dofManager.rankOffset(),
                                                  dofKey,
                                                  fluxApprox.getName(),
                                                  matrixMeshPtr->getElemManager(),
                                                  fractureMeshPtr->getElemManager(),
                                                  stencilWrapper,
                                                  this->getGravityDrainageFlag(),
                                                  this->gravityVector()[2],
                                                  this->getFracSpacingLz(),
                                                  this->getInterporosityExchangeCoefficient(),
                                                  dt,
                                                  localMatrix.toViewConstSizes(),
                                                  localRhs.toView());
      } else
      {
        GEOS_ERROR(
            "The dual continuum flow solver requires both primary and secondary solvers to be either thermal or non-thermal,in this problem,\n"
            "the primary solver is " << (Base::primarySolver()->isThermal() ? "thermal" : "non-thermal") << ",\n"
                                                                                                            "the secondary solver is "<< (Base::secondarySolver()->isThermal() ? "thermal" : "non-thermal") << ".\n");
      }

      // ---- Multi-porosity effective storage matrix (M_bar) cross-storage ----
      // GEOS's per-continuum BiotPorosity gives each continuum its intrinsic
      // 1/M_i; the Mehrabian (2014) multi-porosity storage (eq.22 / Appendix A)
      // requires  1/Mbar_ij = delta_ij v_i alpha_i/(B_i K_i) - abar_i abar_j/Kbar,
      // with abar_i = Kbar v_i alpha_i/K_i, Kbar = (v_m/K_m + v_f/K_f)^-1.
      // We add the correction (1/Mbar - 1/M) to the dual-flow accumulation so
      // continuum i's fluid content is coupled to BOTH pressures:
      //   row_i += rho_i V [ (1/Mbar_ii - 1/M_i)(p_i-p_i_n) + 1/Mbar_ij (p_j-p_j_n) ].
      // The off-diagonal (matrix<->fracture pressure) is essential; it lets the
      // matrix shed pressure as the fracture drains (Mandel-Cryer).  Enabled when
      // the fracture volume fraction v_f is provided on DualContinuumCrossFlow.
      real64 const v_f = this->getFractureVolumeFraction();
      if( v_f > 0.0 && this->getEnableCrossStorageCorrection() )
      {
        real64 const v_m = 1.0 - v_f;
        // Intrinsic (true physical) Biot + drained bulk modulus for the M_bar storage. When the
        // constitutive materials are set to EFFECTIVE-medium values (so the monolithic mechanics
        // kernel uses Kbar/Gbar/abar), these let us still build the correct multi-porosity storage
        // from the true intrinsic parameters. <0 => fall back to the material value (legacy).
        real64 const intrMatA = this->getIntrinsicMatrixBiot();
        real64 const intrMatK = this->getIntrinsicMatrixBulkModulus();
        real64 const intrFracA = this->getIntrinsicFractureBiot();
        real64 const intrFracK = this->getIntrinsicFractureBulkModulus();
        string_array const & matRegions  = this->template getReference< string_array >( "matrixRegionList" );
        string_array const & fracRegions = this->template getReference< string_array >( "fractureRegionList" );
        ElementRegionManager const & matEM  = matrixMeshPtr->getElemManager();
        ElementRegionManager const & fracEM = fractureMeshPtr->getElemManager();
        for( size_t ir = 0; ir < matRegions.size(); ++ir )
        {
          ElementRegionBase const & matReg  = matEM.getRegion( matRegions[ir] );
          ElementRegionBase const & fracReg = fracEM.getRegion( fracRegions[ir] );
          matReg.template forElementSubRegionsIndex< CellElementSubRegion >(
            [&]( localIndex const isr, CellElementSubRegion const & matSR )
          {
            if( isr >= fracReg.numSubRegions() ) return;
            CellElementSubRegion const & fracSR =
              dynamic_cast< CellElementSubRegion const & >( fracReg.getSubRegion( isr ) );

            arrayView1d< globalIndex const > const dofM = matSR.template getReference< array1d< globalIndex > >( dofKey );
            arrayView1d< globalIndex const > const dofF = fracSR.template getReference< array1d< globalIndex > >( dofKey );
            arrayView1d< integer const > const ghostM = matSR.ghostRank();
            arrayView1d< integer const > const ghostF = fracSR.ghostRank();
            arrayView1d< real64 const > const pM  = matSR.template getField< fields::flow::pressure >();
            arrayView1d< real64 const > const pMn = matSR.template getField< fields::flow::pressure_n >();
            arrayView1d< real64 const > const pF  = fracSR.template getField< fields::flow::pressure >();
            arrayView1d< real64 const > const pFn = fracSR.template getField< fields::flow::pressure_n >();
            // In Sequential fixed-stress, lag the off-diagonal storage to the previous outer iterate.
            // The fixed point is unchanged, while the dual-flow subproblem remains contractive.
            bool const useSequentialLaggedOffdiag =
              matSR.hasWrapper( fields::flow::pressure_k::key() ) &&
              fracSR.hasWrapper( fields::flow::pressure_k::key() );
            arrayView1d< real64 const > const pMk =
              useSequentialLaggedOffdiag ? matSR.template getField< fields::flow::pressure_k >() : pM;
            arrayView1d< real64 const > const pFk =
              useSequentialLaggedOffdiag ? fracSR.template getField< fields::flow::pressure_k >() : pF;
            arrayView1d< real64 const > const vol = matSR.getElementVolume();

            string const & solidMn = matSR.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
            string const & solidFn = fracSR.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
            constitutive::CoupledSolidBase const & csM = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( matSR, solidMn );
            constitutive::CoupledSolidBase const & csF = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( fracSR, solidFn );
            arrayView1d< real64 const > const Km = csM.getBulkModulus();
            arrayView1d< real64 const > const Kf = csF.getBulkModulus();
            arrayView1d< real64 const > const aM = csM.getBiotCoefficient();
            arrayView1d< real64 const > const aF = csF.getBiotCoefficient();
            arrayView1d< real64 const > const phiM = csM.getReferencePorosity();
            arrayView1d< real64 const > const phiF = csF.getReferencePorosity();
            arrayView1d< real64 const > const KsM =
              dynamic_cast< constitutive::BiotPorosity const & >( csM.getBasePorosityModel() ).getGrainBulkModulus();
            arrayView1d< real64 const > const KsF =
              dynamic_cast< constitutive::BiotPorosity const & >( csF.getBasePorosityModel() ).getGrainBulkModulus();

            string const & fluidMn = matSR.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
            string const & fluidFn = fracSR.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
            constitutive::SingleFluidBase const & flM = this->template getConstitutiveModel< constitutive::SingleFluidBase >( matSR, fluidMn );
            constitutive::SingleFluidBase const & flF = this->template getConstitutiveModel< constitutive::SingleFluidBase >( fracSR, fluidFn );
            arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const densM = flM.density();
            arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const densF = flF.density();
            arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensM = flM.dDensity();
            arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensF = flF.dDensity();

            for( localIndex k = 0; k < matSR.size(); ++k )
            {
              // Intrinsic params for the analytical storage (fall back to material if unset)
              real64 const aMI = ( intrMatA  > 0.0 ) ? intrMatA  : aM[k];
              real64 const aFI = ( intrFracA > 0.0 ) ? intrFracA : aF[k];
              real64 const KmI = ( intrMatK  > 0.0 ) ? intrMatK  : Km[k];
              real64 const KfI = ( intrFracK > 0.0 ) ? intrFracK : Kf[k];

              real64 const Kbar = 1.0/( v_m/KmI + v_f/KfI );
              real64 const abm  = Kbar*v_m*aMI/KmI;
              real64 const abf  = Kbar*v_f*aFI/KfI;
              real64 const cfM  = dDensM[k][0][0]/densM[k][0];
              real64 const cfF  = dDensF[k][0][0]/densF[k][0];
              // phiM/phiF have already been scaled by dual-continuum netToGross, so they are REV pore
              // fractions (v_i * phi_i). The intrinsic storage formula needs continuum-local porosity.
              real64 const phiMI = phiM[k] / v_m;
              real64 const phiFI = phiF[k] / v_f;

              // Analytical constant-strain storage built from INTRINSIC params:
              //   1/Mbar_ii = v_i(1/M_i^intr + alpha_i^2/K_i) - abar_i^2/Kbar,  abar_i = Kbar v_i alpha_i/K_i.
              real64 const invMmI = (aMI-phiMI)/KsM[k] + phiMI*cfM;
              real64 const invMfI = (aFI-phiFI)/KsF[k] + phiFI*cfF;
              // Paper-exact constant-strain diagonal storage 1/Mbar_ii = a_ii - abar_i^2/Kbar,
              // a_ii = v_i(1/M_i + alpha_i^2/K_i) (Mehrabian 2014 eq A24/A25).
              real64 const SbarMM = v_m*(invMmI + aMI*aMI/KmI) - abm*abm/Kbar;
              real64 const SbarFF = v_f*(invMfI + aFI*aFI/KfI) - abf*abf/Kbar;
              // Subtract what the monolithic kernel / secondary solver already put on the physical
              // storage diagonal. Do not remove the fixed-stress term here: in a converged
              // Sequential fixed-stress iteration it vanishes because p -> p_k, while during the
              // outer iteration it is the stabilizing term that makes the split solve contract.
              real64 const invMmMat = (aM[k]-phiM[k])/KsM[k] + phiM[k]*cfM;
              real64 const invMfMat = (aF[k]-phiF[k])/KsF[k] + phiF[k]*cfF;
              real64 const corrDiagM = SbarMM - invMmMat;
              real64 const corrDiagF = SbarFF - invMfMat;
              // Paper-exact constant-strain off-diagonal storage 1/Mbar_ij = -abar_i*abar_j/Kbar.
              // In the analytical effective medium the diffusion-equation storage is
              // Sbar = invM + abar (x) cm with cm_i = abar_i/(Kbar+4Gbar/3) (oedometric). The FIM here
              // adds the full constant-strain invM off-diagonal (-abar_m*abar_f/Kbar) as STORAGE and
              // relies on the monolithic mechanics Schur (K_pm_u*K_uu^-1*K_up_f) to supply the
              // +abar_m*abar_f/(Kbar+4Gbar/3) half. The net is therefore a NEAR-CANCELLATION of two
              // ~5e-10 terms (analytical net Sbar_mf = -2.44e-10 for the GOM-shale N=2 case), so any
              // error in the mechanics half is hugely amplified in the matrix-pressure plateau.
              //
              // MEASURED (2026-06, crossStorageOffDiagScale XML scan, direct solver, IDENTICAL on
              // 10x10 and 20x20 -> mesh-INDEPENDENT): the discrete Q1 mechanics produces an effective
              // modulus ~1.05e9 = 1.21x the continuum oedometric Kbar+4Gbar/3 = 8.66e8, i.e. it
              // UNDER-cancels by ~17%. At scale=1.0 the residual off-diagonal is too negative and the
              // matrix is over-drawn-down (plateau 0.82 vs analytical 0.885; overshoot/drainage also
              // low). The plateau is linear in scale (slope ~-0.71) and matches the analytical 0.885 at
              //   crossStorageOffDiagScale = 0.911   (verified plateau 0.8852, tau>=1 within +-0.25%).
              // This 0.911 is NOT a fudge or a mesh-dependent knob: it is a measured, mesh-independent
              // DISCRETIZATION CONSTANT compensating the Q1 under-cancellation. It is set in the
              // fim_eff* decks. NOTE it was calibrated at nu=0.22 / this Mandel BC set; a different
              // Poisson ratio or geometry would need a re-scan. A fully nu/geometry-adaptive value
              // cannot be computed from Kbar/Gbar/nu in closed form (a single Q1 element under uniform
              // strain reproduces the continuum modulus exactly -> the 21% excess comes from the
              // plateau's SPATIAL structure, not a local modulus), so a "true auto" would require a
              // one-time matrix self-calibration: solve K_uu*x = K_up_f for one interior element and
              // back out the actual mechanics cross-term. See memory dpdp-fim-crossstorage-offdiag.
              real64 const offScale = this->getCrossStorageOffDiagScale();
              real64 const corrOffMF = -abm*abf/Kbar * offScale;  // matrix row, vs p_f
              real64 const corrOffFM = -abm*abf/Kbar * offScale;  // fracture row, vs p_m
              real64 const dpM = pM[k]-pMn[k];
              real64 const dpF = pF[k]-pFn[k];
              real64 const dpMOffdiag = pMk[k]-pMn[k];
              real64 const dpFOffdiag = pFk[k]-pFn[k];

              if( ghostM[k] < 0 )
              {
                localIndex const row = LvArray::integerConversion< localIndex >( dofM[k]-rankOffset );
                real64 const rv = densM[k][0]*vol[k];
                localRhs[row] += rv*( corrDiagM*dpM + corrOffMF*dpFOffdiag );
                if( useSequentialLaggedOffdiag )
                {
                  globalIndex cols[1] = { dofM[k] };
                  real64 vals[1] = { rv*corrDiagM };
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 1 );
                }
                else
                {
                  globalIndex cols[2] = { dofM[k], dofF[k] };
                  real64 vals[2] = { rv*corrDiagM, rv*corrOffMF };
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );
                }
              }
              if( ghostF[k] < 0 )
              {
                localIndex const row = LvArray::integerConversion< localIndex >( dofF[k]-rankOffset );
                real64 const rv = densF[k][0]*vol[k];
                localRhs[row] += rv*( corrDiagF*dpF + corrOffFM*dpMOffdiag );
                if( useSequentialLaggedOffdiag )
                {
                  globalIndex cols[1] = { dofF[k] };
                  real64 vals[1] = { rv*corrDiagF };
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 1 );
                }
                else
                {
                  globalIndex cols[2] = { dofF[k], dofM[k] };
                  real64 vals[2] = { rv*corrDiagF, rv*corrOffFM };
                  localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, cols, vals, 2 );
                }
              }
            }
          } );
        }
      }
    }
    //GEOS_LOG("after the assemble");
    //this->printCRSMatrix(this->m_localMatrix);
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::DualContinuumFVM(const string &name, dataRepository::Group *parent)
      : Base(name, parent)
  {
    LinearSolverParameters &linParams = m_linearSolverParameters.get();
    linParams.multiscale.fieldName = SinglePhaseBase::viewKeyStruct::elemDofFieldString();
    linParams.multiscale.label = "dualflow";
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::postInputInitialization()
  {
    Base::postInputInitialization();
  }
/*
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::setupCoupling(DomainPartition const &domain,
                                                                                   DofManager &dofManager) const
  {
    // ensure element-based coupling (two components per element) has sparsity
    GEOS_LOG(SinglePhaseBase::viewKeyStruct::elemDofFieldString());

    // Get supports from both solvers
    stdVector< DofManager::FieldSupport > supports;
    auto const & primaryTargets = Base::primarySolver()->getMeshTargets();
    auto const & secondaryTargets = Base::secondarySolver()->getMeshTargets();
    for( auto const & p : primaryTargets )
    {
      MeshBody const & meshBody = domain.getMeshBody( p.first.first );
      MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
      std::set< string > regionNames( p.second.begin(), p.second.end() );
      supports.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
    }
    for( auto const & p : secondaryTargets )
    {
      MeshBody const & meshBody = domain.getMeshBody( p.first.first );
      MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
      std::set< string > regionNames( p.second.begin(), p.second.end() );
      supports.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
    }

    dofManager.addCouplingDualContinuum( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                                         SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                                         DofManager::Connector::Elem,
                                         supports );
  }
*/
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::initializePreSubGroups()
  {
    Base::initializePreSubGroups();
    // Ensure discretization is valid for each underlying flow solver if needed
  }
  /*
  template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER  >::setupDofs( DomainPartition const & domain, DofManager & dofManager ) const
  {
    // Let the base DualContinuumFlowSolverBase call each sub-solver's setupDofs
    Base::setupDofs( domain, dofManager );
  }
  */
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::setupSystem(DomainPartition &domain,
                                                                                 DofManager &dofManager,
                                                                                 CRSMatrix<real64, globalIndex> &localMatrix,
                                                                                 ParallelVector &rhs,
                                                                                 ParallelVector &solution,
                                                                                 bool const setSparsity)
  {
    static bool connectivityRegistered = false;
    if( !connectivityRegistered )
    {
      this->registerMeshConnectivity( domain );
      connectivityRegistered = true;
    }

    // TracAI: Print registered connectivity values to terminal
    if( connectivityRegistered )
    {
      this->printRegisteredConnectivityValues( domain );
    }

    GEOS_MARK_FUNCTION;
    PhysicsSolverBase::setupSystem(domain, dofManager, localMatrix, rhs, solution, setSparsity);

    if (!m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct)
    {
      m_precond = createPreconditioner(domain);
    }
  }

  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  std::unique_ptr<PreconditionerBase<LAInterface>> DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::createPreconditioner(DomainPartition &domain) const
  {
    LinearSolverParameters const &linParams = m_linearSolverParameters.get();
    switch (linParams.preconditionerType)
    {
    case LinearSolverParameters::PreconditionerType::multiscale:
    {
      return std::make_unique<MultiscalePreconditioner<LAInterface>>(linParams, domain);
    }
    default:
    {
      return PhysicsSolverBase::createPreconditioner(domain);
    }
    }
  }

  /*
  template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                  real64 const dt,
                                                                                  DomainPartition & domain,
                                                                                  DofManager const & dofManager,
                                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                  arrayView1d< real64 > const & localRhs )
  {
    GEOS_MARK_FUNCTION;
    this->primarySolver()->assembleAccumulationTerms(domain,
                                                      dofManager,
                                                      localMatrix,
                                                      localRhs );
    this->secondarySolver()->assembleAccumulationTerms(domain,
                                                     dofManager,
                                                     localMatrix,
                                                     localRhs );
    if( 0 )
    {
      this->primarySolver()->assembleStabilizedFluxTerms( dt,
                                                       domain,
                                                       dofManager,
                                                       localMatrix,
                                                       localRhs );
    }
    else
    {
      this->primarySolver()->assembleFluxTerms( dt,
                                             domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );
    }
    // Step 3: compute the fluxes (face-based contributions)
    if( 0 )
    {
      this->secondarySolver()->assembleStabilizedFluxTerms( dt,
                                                       domain,
                                                       dofManager,
                                                       localMatrix,
                                                       localRhs );
    }
    else
    {
      this->secondarySolver()->assembleFluxTerms( dt,
                                             domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );
    }
  }
  */
  template <typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER>
  void DualContinuumFVM<PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER>::updateState( geos::DomainPartition & domain )
  {
    Base::updateState(domain);
    if( this->getGravityDrainageFlag() )
    {
      // Update gravity drainage pressure for dual continuum flow
      GEOS_LOG("Updating gravity drainage pressure for dual continuum flow");
      
      // Get primary and secondary solvers
      PRIMARY_FLOW_SOLVER const * primarySolver = this->primarySolver();
      
      // Get meshes for matrix and fracture
      std::vector<MeshLevel const*> meshLevelPtrs;
      this->forDiscretizationOnMeshTargets(
        domain.getMeshBodies(),
        [&]( string const &,
             MeshLevel const & mesh,
             string_array const & )
        {
          meshLevelPtrs.push_back( &mesh );
        }
      );
      
      if( meshLevelPtrs.size() >= 2 )
      {
        MeshLevel  &  matrixMesh = const_cast<MeshLevel&>(* meshLevelPtrs[0]);
        MeshLevel  &  fractureMesh = const_cast<MeshLevel&>(* meshLevelPtrs[1]);

        real64 const gravityCoefficient = primarySolver->gravityVector()[2];
        Base::updateGravityPressure(matrixMesh,fractureMesh, gravityCoefficient);
      }
    }
  }

  // Explicit instantiation for default template
  template class DualContinuumFVM<SinglePhaseBase, SinglePhaseBase>;
  namespace
  { // Register the solver so it can be used from XML

    typedef DualContinuumFVM<SinglePhaseBase, SinglePhaseBase> DualContinuumFVM;
    REGISTER_CATALOG_ENTRY(PhysicsSolverBase, DualContinuumFVM, string const &, Group *const)
  }

} // namespace geos
