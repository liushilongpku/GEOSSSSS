/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DualContinuumFlowSolverBase.hpp
 *
 * @brief A coupled solver that binds two flow solvers for dual-continuum/dual-porosity-style models.
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/SinglePhaseBaseDpdk.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/FieldIdentifiers.hpp"
#include "mesh/InterObjectRelation.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "codingUtilities/Utilities.hpp"
#include <map>
#include <tuple>
#include "DualContinuumCrossFlow.hpp"
namespace geos
{
    namespace stabilization
    {
  enum class StabilizationType : integer;
}
template< typename PRIMARY_FLOW_SOLVER, typename SECONDARY_FLOW_SOLVER >
class DualContinuumFlowSolverBase : public CoupledSolver< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >
{
public:

  using Base = CoupledSolver< PRIMARY_FLOW_SOLVER, SECONDARY_FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Primary = 0,
    Secondary = 1
  };

  static string coupledSolverAttributePrefix() { return "dualcontinuum"; }

  DualContinuumFlowSolverBase( string const & name, dataRepository::Group * parent )
    : Base( name, parent ),
    m_crossFlow( viewKeyStruct::DualContinuumCrossFlow(), this )

  {
    // dual-flow has two scalar dofs per node (primary + secondary), set a reasonable default
    LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();
    linearSolverParameters.dofsPerNode = 2;
    linearSolverParameters.multiscale.label = "dualflow";
    // default transfer coefficient (per-element scalar). Users may override in input

    // m_transferCoefficient = 0.0;

    this->registerWrapper( viewKeyStruct::transferCoefficientString(), &m_transferCoefficient ).
      setApplyDefaultValue( 0.0 ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "Inter-continuum transfer coefficient (per-element scalar, default 0)" );

    // TracAI: Added to store dual continuum region pairs from XML
    this->registerWrapper( viewKeyStruct::matrixRegionList(), &m_matrixRegionList ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "List of matrix regions" );

    this->registerGroup( viewKeyStruct::DualContinuumCrossFlow(), &m_crossFlow );

    this->registerWrapper( viewKeyStruct::fractureRegionList(), &m_fractureRegionList ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "List of fracture regions" );

  }

  virtual void postInputInitialization() override {
    Base::postInputInitialization();
    /*
     *     this->getMeshTargets();

    std::set< string > targetMeshBodiesSet;
    for( auto const & pair : primarySolver()->getMeshTargets() )
    {
      targetMeshBodiesSet.insert( pair.first.first );
    }

    std::set< string > targetMeshBodiesSetFracture;
    for( auto const & pair : secondarySolver()->getMeshTargets() )
    {
      targetMeshBodiesSetFracture.insert( pair.first.first );
    }
    */
    //DomainPartition& domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
    //MeshBody & matrix = domain.getMeshBody("mesh1");
    //MeshBody & fracture = domain.getMeshBody("mesh2");
    //MeshLevel & primaryMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    //MeshLevel & secondaryMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );

    //m_crossFlow.initialize( primaryMesh, secondaryMesh );
    // add validation of solver types if needed
  }

  virtual void initializePostInitialConditionsPreSubGroups() override
  {
    Base::initializePostInitialConditionsPreSubGroups();
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      solver->initializePostInitialConditionsPreSubGroupsPublic();
    } );

    DomainPartition& domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
    MeshBody & matrix = domain.getMeshBody("mesh1");
    MeshBody & fracture = domain.getMeshBody("mesh2");
    MeshLevel & primaryMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & secondaryMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );
    this->registerMeshConnectivity( domain );
    m_crossFlow.setupCrossFlow( domain, primaryMesh, secondaryMesh );//TODO@LSL 这块需要对多region进行支持，可能会报错
  };

  virtual void initializePostInitialConditionsPostSubGroups() override
  {
    Base::initializePostInitialConditionsPostSubGroups();
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      solver->initializePostInitialConditionsPostSubGroupsPublic();
    } );

    DomainPartition& domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
    MeshBody & matrix = domain.getMeshBody("mesh1");
    MeshBody & fracture = domain.getMeshBody("mesh2");
    MeshLevel & primaryMesh = matrix.getMeshLevels().getGroup< MeshLevel >( 0 );
    MeshLevel & secondaryMesh = fracture.getMeshLevels().getGroup< MeshLevel >( 0 );
    m_crossFlow.setupGravityDrainagePressure(primaryMesh, secondaryMesh, this->gravityVector()[2]);
  }
  /**
   * @brief accessor for the pointer to the primary flow solver
   * @return a pointer to the primary flow solver
   */
  PRIMARY_FLOW_SOLVER * primarySolver() const
  {
    return std::get< toUnderlying( SolverType::Primary ) >( m_solvers );
  }
  /**
   * @brief accessor for the pointer to the secondary flow solver
   * @return a pointer to the secondary flow solver
   */
  SECONDARY_FLOW_SOLVER * secondarySolver() const
  {
    return std::get< toUnderlying( SolverType::Secondary ) >( m_solvers );
  }

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    // ensure both flow solvers register their dofs
    primarySolver()->setupDofs( domain, dofManager );
    secondarySolver()->setupDofs( domain, dofManager );

    // allow base to set any additional coupling DOFs/attributes
    this->setupCoupling( domain, dofManager );
  }

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override
  {
    // ensure element-based coupling (two components per element) has sparsity

    // Get supports from both solvers
    /*
    stdVector< DofManager::FieldSupport > supports;
    auto const & primaryTargets = primarySolver()->getMeshTargets();
    auto const & secondaryTargets = secondarySolver()->getMeshTargets();
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
    */
    dofManager.addCouplingDualContinuum( m_matrixRegionList,
                                         m_fractureRegionList,
                                         PRIMARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString(),
                                         SECONDARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString(),
                                         DofManager::Connector::Elem);
  }

  void initializeState( DomainPartition & domain )
  {
    primarySolver()->initializeState( domain );
    secondarySolver()->initializeState( domain );
  }

  void assembleFluxTerms( real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const
  {
    primarySolver()->assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );
    secondarySolver()->assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );
  }

  void assembleStabilizedFluxTerms( real64 const dt,
                                    DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs ) const
  {
    primarySolver()->assembleStabilizedFluxTerms( dt, domain, dofManager, localMatrix, localRhs );
    secondarySolver()->assembleStabilizedFluxTerms( dt, domain, dofManager, localMatrix, localRhs );
  }

  real64 updateFluidState( ElementSubRegionBase & subRegion ) const
  {
    return primarySolver()->updateFluidState( subRegion );
  }

  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
  {
    primarySolver()->updateSolidInternalEnergyModel( dataGroup );
    secondarySolver()->updateSolidInternalEnergyModel( dataGroup );
  }

  /**
   * @brief Get the fracture spacing Lz
   * @return The fracture spacing in z-direction
   */
  real64 getFracSpacingLz() const
  {
    return m_crossFlow.getFracSpacingLz();
  }

  /**
   * @brief Get the direct interporosity exchange coefficient Gamma
   * @return Gamma in [Pa^{-1} s^{-1}]; 0 = use Kazemi shape factor
   */
  real64 getInterporosityExchangeCoefficient() const
  {
    return m_crossFlow.getInterporosityExchangeCoefficient();
  }

  /// Fracture volume fraction v_f (from DualContinuumCrossFlow); <0 = unset
  real64 getFractureVolumeFraction() const
  {
    return m_crossFlow.getFractureVolumeFraction();
  }

  /// Intrinsic-parameter accessors for the FIM multi-porosity storage (<0 = use material value)
  real64 getIntrinsicMatrixBiot() const { return m_crossFlow.getIntrinsicMatrixBiot(); }
  real64 getIntrinsicMatrixBulkModulus() const { return m_crossFlow.getIntrinsicMatrixBulkModulus(); }
  real64 getIntrinsicFractureBiot() const { return m_crossFlow.getIntrinsicFractureBiot(); }
  real64 getIntrinsicFractureBulkModulus() const { return m_crossFlow.getIntrinsicFractureBulkModulus(); }
  real64 getCrossStorageOffDiagScale() const { return m_crossFlow.getCrossStorageOffDiagScale(); }

  /// Setters used by the useIntrinsicInput path (auto-derived intrinsics for FIM storage)
  void setIntrinsicMatrixBiot( real64 const v ) { m_crossFlow.setIntrinsicMatrixBiot( v ); }
  void setIntrinsicMatrixBulkModulus( real64 const v ) { m_crossFlow.setIntrinsicMatrixBulkModulus( v ); }
  void setIntrinsicFractureBiot( real64 const v ) { m_crossFlow.setIntrinsicFractureBiot( v ); }
  void setIntrinsicFractureBulkModulus( real64 const v ) { m_crossFlow.setIntrinsicFractureBulkModulus( v ); }

  /// The multi-porosity cross-storage correction in assembleCouplingTerms is a
  /// SEQUENTIAL-only construct: it adds a residual term with an incomplete (non-linearized)
  /// Jacobian and inflates the matrix pressure-storage diagonal, which makes a FullyImplicit
  /// Newton diverge. For FIM the matrix storage comes from the monolithic kernel and the
  /// matrix<->fracture cross-storage is added as a consistent term instead. The coupled
  /// poromechanics solver disables this flag when couplingType=FullyImplicit.
  void setEnableCrossStorageCorrection( bool const flag ) { m_enableCrossStorageCorrection = flag; }
  bool getEnableCrossStorageCorrection() const { return m_enableCrossStorageCorrection; }

  // Support for PoromechanicsSolver expectations (delegated to primary/secondary solvers)
  integer isThermal() const
  {
    // Assume both solvers have the same thermal flag
    return primarySolver()->isThermal();
  }

  void enableFixedStressPoromechanicsUpdate()
  {
    primarySolver()->enableFixedStressPoromechanicsUpdate();
    secondarySolver()->enableFixedStressPoromechanicsUpdate();
  }

  void enableJumpStabilization()
  {
    primarySolver()->enableJumpStabilization();
    secondarySolver()->enableJumpStabilization();
  }

  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
  {
    primarySolver()->setKeepVariablesConstantDuringInitStep( keepVariablesConstantDuringInitStep );
    secondarySolver()->setKeepVariablesConstantDuringInitStep( keepVariablesConstantDuringInitStep );
  }

  void updatePorosityAndPermeability( ElementSubRegionBase & subRegion )
  {
    // Forward to the appropriate overload on the underlying solvers depending on subregion type
    if( auto * cellSub = dynamic_cast< CellElementSubRegion * >( &subRegion ) )
    {
      primarySolver()->updatePorosityAndPermeability( *cellSub );
      secondarySolver()->updatePorosityAndPermeability( *cellSub );
    }
    else if( auto * surfSub = dynamic_cast< SurfaceElementSubRegion * >( &subRegion ) )
    {
      primarySolver()->updatePorosityAndPermeability( *surfSub );
      secondarySolver()->updatePorosityAndPermeability( *surfSub );
    }
    else
    {
      GEOS_ERROR( "Unsupported ElementSubRegionBase type in updatePorosityAndPermeability" );
    }
  }

public:
  /**
   * @brief Register mesh connectivity between matrix and fracture regions
   * @param domain The domain partition containing both meshes
   * @note TracAI: Added to support dual continuum flow solver with mesh connectivity
   */
  void registerMeshConnectivity( DomainPartition & domain )
  {
    // Check if matrix and fracture region lists are defined
    if( !m_matrixRegionList.empty() && !m_fractureRegionList.empty() )
    {
      GEOS_LOG( "Registering mesh connectivity based on matrix and fracture region lists" );
      
      // Check if the number of regions in both lists match
      if( m_matrixRegionList.size() != m_fractureRegionList.size() )
      {
        GEOS_ERROR( "Matrix and fracture region lists must have the same number of elements: " << m_matrixRegionList.size() << " matrix regions, " << m_fractureRegionList.size() << " fracture regions provided" );
        return;
      }
      
      // Get mesh1 (primary) and mesh2 (secondary) if available
      if( domain.getMeshBodies().numSubGroups() >= 2 )
      {
        MeshBody & meshBody1 = domain.getMeshBody( 0 );
        MeshBody & meshBody2 = domain.getMeshBody( 1 );
        
        MeshLevel & mesh1 = meshBody1.getMeshLevels().getGroup< MeshLevel >( 0 );
        MeshLevel & mesh2 = meshBody2.getMeshLevels().getGroup< MeshLevel >( 0 );
        
        ElementRegionManager & elemManager1 = mesh1.getGroup< ElementRegionManager >( "ElementRegions" );
        ElementRegionManager & elemManager2 = mesh2.getGroup< ElementRegionManager >( "ElementRegions" );
        
        // Process each pair of regions (matrix <-> fracture)
        for( size_t i = 0; i < m_matrixRegionList.size(); ++i )
        {
          std::string matrixRegion = m_matrixRegionList[i];
          std::string fractureRegion = m_fractureRegionList[i];
          
          GEOS_LOG( "Processing dual continuum pair: " << matrixRegion << " (matrix) <-> " << fractureRegion << " (fracture)" );
          
          // Try to find the matrix region in mesh1
          ElementRegionBase * matrixRegionPtr = nullptr;
          elemManager1.forElementRegions( [&]( ElementRegionBase & region )
          {
            if( region.getName() == matrixRegion )
            {
              matrixRegionPtr = &region;
            }
          } );
          
          // Try to find the fracture region in mesh2
          ElementRegionBase * fractureRegionPtr = nullptr;
          elemManager2.forElementRegions( [&]( ElementRegionBase & region )
          {
            if( region.getName() == fractureRegion )
            {
              fractureRegionPtr = &region;
            }
          } );
          
          // Register connectivity between subregions of the found regions
          if( matrixRegionPtr && fractureRegionPtr )
          {
            if( auto * cellMatrixRegion = dynamic_cast< CellElementRegion * >( matrixRegionPtr ) )
            {
              if( auto * cellFractureRegion = dynamic_cast< CellElementRegion * >( fractureRegionPtr ) )
              {
                // Get subregions from both regions
                stdVector< ElementSubRegionBase * > matrixSubRegions;
                cellMatrixRegion->forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
                {
                  matrixSubRegions.push_back( &subRegion );
                } );
                
                stdVector< ElementSubRegionBase * > fractureSubRegions;
                cellFractureRegion->forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
                {
                  fractureSubRegions.push_back( &subRegion );
                } );
                
                // Register connectivity based on subregion order
                size_t numSubRegions = std::min( matrixSubRegions.size(), fractureSubRegions.size() );
                GEOS_LOG( "Registering connectivity for " << numSubRegions << " subregion pairs between " << matrixRegion << " and " << fractureRegion );
                GEOS_LOG( "Assuming subregions are paired by their registration order" );
                
                for( size_t j = 0; j < numSubRegions; ++j )
                {
                  if( auto * matrixSubRegion = dynamic_cast< CellElementSubRegion * >( matrixSubRegions[j] ) )
                  {
                    if( auto * fractureSubRegion = dynamic_cast< CellElementSubRegion * >( fractureSubRegions[j] ) )
                    {
                      GEOS_LOG( "Pairing subregion " << matrixSubRegion->getName() << " (matrix) with " << fractureSubRegion->getName() << " (fracture)" );
                      
                      // Register connectivity from matrix to fracture
                      localIndex numElements = matrixSubRegion->size();
                      if( !matrixSubRegion->hasWrapper( viewKeyStruct::mesh1ToMesh2ConnectivityString() ) )
                      {
                        matrixSubRegion->registerWrapper< array1d< localIndex > >( viewKeyStruct::mesh1ToMesh2ConnectivityString() )
                          .setApplyDefaultValue( false )
                          .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
                          .setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
                          .setDescription( "Connectivity from matrix elements to fracture elements" )
                          .setRegisteringObjects( this->getName() );
                      }
                      
                      auto & matrixToFractureConnectivity =
                        matrixSubRegion->getReference< array1d< localIndex > >( viewKeyStruct::mesh1ToMesh2ConnectivityString() );
                      matrixToFractureConnectivity.resize( numElements );
                      arrayView1d< globalIndex const > const matrixLocalToGlobal = matrixSubRegion->localToGlobalMap();
                      auto const & fractureGlobalToLocal = fractureSubRegion->globalToLocalMap();
                      for( localIndex k = 0; k < numElements; ++k )
                      {
                        auto const iter = fractureGlobalToLocal.find( matrixLocalToGlobal[k] );
                        GEOS_ERROR_IF( iter == fractureGlobalToLocal.end(),
                                       "Unable to build dual-continuum connectivity from matrix subregion "
                                       << matrixSubRegion->getName() << " to fracture subregion "
                                       << fractureSubRegion->getName() << ": fracture subregion does not contain global element "
                                       << matrixLocalToGlobal[k] << " on this rank." );
                        matrixToFractureConnectivity[k] = iter->second;
                      }
                      
                      // Register connectivity from fracture to matrix
                      numElements = fractureSubRegion->size();
                      if( !fractureSubRegion->hasWrapper( viewKeyStruct::mesh2ToMesh1ConnectivityString() ) )
                      {
                        fractureSubRegion->registerWrapper< array1d< localIndex > >( viewKeyStruct::mesh2ToMesh1ConnectivityString() )
                          .setApplyDefaultValue( false )
                          .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
                          .setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
                          .setDescription( "Connectivity from fracture elements to matrix elements" )
                          .setRegisteringObjects( this->getName() );
                      }
                      
                      auto & fractureToMatrixConnectivity =
                        fractureSubRegion->getReference< array1d< localIndex > >( viewKeyStruct::mesh2ToMesh1ConnectivityString() );
                      fractureToMatrixConnectivity.resize( numElements );
                      arrayView1d< globalIndex const > const fractureLocalToGlobal = fractureSubRegion->localToGlobalMap();
                      auto const & matrixGlobalToLocal = matrixSubRegion->globalToLocalMap();
                      for( localIndex k = 0; k < numElements; ++k )
                      {
                        auto const iter = matrixGlobalToLocal.find( fractureLocalToGlobal[k] );
                        GEOS_ERROR_IF( iter == matrixGlobalToLocal.end(),
                                       "Unable to build dual-continuum connectivity from fracture subregion "
                                       << fractureSubRegion->getName() << " to matrix subregion "
                                       << matrixSubRegion->getName() << ": matrix subregion does not contain global element "
                                       << fractureLocalToGlobal[k] << " on this rank." );
                        fractureToMatrixConnectivity[k] = iter->second;
                      }
                    }
                  }
                }
                
                if( matrixSubRegions.size() != fractureSubRegions.size() )
                {
                  GEOS_LOG( "Warning: Number of subregions differs between " << matrixRegion << " (" << matrixSubRegions.size() << ") and " << fractureRegion << " (" << fractureSubRegions.size() << ")" );
                  GEOS_LOG( "Only the first " << numSubRegions << " subregions will be paired" );
                }
              }
            }
          }
          else
          {
            GEOS_LOG( "Warning: Could not find both regions for pair " << matrixRegion << " <-> " << fractureRegion );
            if( !matrixRegionPtr )
            {
              GEOS_LOG( "  Matrix region " << matrixRegion << " not found in mesh1" );
            }
            if( !fractureRegionPtr )
            {
              GEOS_LOG( "  Fracture region " << fractureRegion << " not found in mesh2" );
            }
          }
        }
        
        GEOS_LOG( "Registered mesh connectivity for " << m_matrixRegionList.size() << " dual continuum region pairs" );
      }
    }
    else
    {
      // Fallback to default behavior if no region pairs are defined
      GEOS_LOG( "No matrix or fracture region lists defined, using default mesh connectivity" );
      
      if( domain.getMeshBodies().numSubGroups() >= 2 )
      {
        // Get mesh1 (primary) and mesh2 (secondary)
        MeshBody & meshBody1 = domain.getMeshBody( 0 );
        MeshBody & meshBody2 = domain.getMeshBody( 1 );
        
        // Get the latest mesh levels using the first available mesh level
        MeshLevel & mesh1 = meshBody1.getMeshLevels().getGroup< MeshLevel >( 0 );
        MeshLevel & mesh2 = meshBody2.getMeshLevels().getGroup< MeshLevel >( 0 );
        
        // Get element managers for both meshes
        ElementRegionManager & elemManager1 = mesh1.getGroup< ElementRegionManager >( "ElementRegions" );
        ElementRegionManager & elemManager2 = mesh2.getGroup< ElementRegionManager >( "ElementRegions" );
        
        // Register connectivity for each CellElementSubRegion in mesh1 to corresponding subregions in mesh2
        // Iterate through all element regions in mesh1
        elemManager1.forElementRegions( [&]( ElementRegionBase & region1 )
        {
          // Check if this is a CellElementRegion
          if( auto * cellRegion1 = dynamic_cast< CellElementRegion * >( &region1 ) )
          {
            // Iterate through all element subregions in this region
            cellRegion1->forElementSubRegions( [&]( ElementSubRegionBase & subRegion1 )
            {
              // Check if this is a CellElementSubRegion
              if( auto * cellSubRegion1 = dynamic_cast< CellElementSubRegion * >( &subRegion1 ) )
              {
                // Register the connectivity field for mesh1 to mesh2
                // We'll use a simple array1d for the connectivity
                localIndex numElements = cellSubRegion1->size();
                
                // Register the connectivity field
                cellSubRegion1->registerWrapper< array1d< localIndex > >( viewKeyStruct::mesh1ToMesh2ConnectivityString() )
                  .setApplyDefaultValue( false )
                  .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
                  .setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
                  .setDescription( "Connectivity from mesh1 elements to mesh2 elements" )
                  .setRegisteringObjects( this->getName() );
                
                // Get the field and fill it with one-to-one mapping
                auto & connectivity = cellSubRegion1->getReference< array1d< localIndex > >( viewKeyStruct::mesh1ToMesh2ConnectivityString() );
                connectivity.resize( numElements );
                for( localIndex i = 0; i < numElements; ++i )
                {
                  connectivity[i] = i;
                }
              }
            } );
          }
        } );
        
        // Register connectivity for each CellElementSubRegion in mesh2 to corresponding subregions in mesh1
        // Iterate through all element regions in mesh2
        elemManager2.forElementRegions( [&]( ElementRegionBase & region2 )
        {
          // Check if this is a CellElementRegion
          if( auto * cellRegion2 = dynamic_cast< CellElementRegion * >( &region2 ) )
          {
            // Iterate through all element subregions in this region
            cellRegion2->forElementSubRegions( [&]( ElementSubRegionBase & subRegion2 )
            {
              // Check if this is a CellElementSubRegion
              if( auto * cellSubRegion2 = dynamic_cast< CellElementSubRegion * >( &subRegion2 ) )
              {
                // Register the connectivity field for mesh2 to mesh1
                // We'll use a simple array1d for the connectivity
                localIndex numElements = cellSubRegion2->size();
                
                // Register the connectivity field
                cellSubRegion2->registerWrapper< array1d< localIndex > >( viewKeyStruct::mesh2ToMesh1ConnectivityString() )
                  .setApplyDefaultValue( false )
                  .setPlotLevel( dataRepository::PlotLevel::NOPLOT )
                  .setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
                  .setDescription( "Connectivity from mesh2 elements to mesh1 elements" )
                  .setRegisteringObjects( this->getName() );
                
                // Get the field and fill it with one-to-one mapping
                auto & connectivity = cellSubRegion2->getReference< array1d< localIndex > >( viewKeyStruct::mesh2ToMesh1ConnectivityString() );
                connectivity.resize( numElements );
                for( localIndex i = 0; i < numElements; ++i )
                {
                  connectivity[i] = i;
                }
              }
            } );
          }
        } );
        
        GEOS_LOG( "Registered mesh connectivity between " << meshBody1.getName() << " and " << meshBody2.getName() );
      }
    }
  }

  /// Assemble coupling blocks between the two flow solvers (exchange/transfer terms)
  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    GEOS_ERROR("should be override");
  }

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override
  {
    // The primary and secondary flow solvers usually register the same element DOF field name
    // (e.g. "compositionalVariables") on different mesh supports. DofManager::addVectorToField()
    // applies an update to every support of that field, so calling both sub-solvers would add the
    // same Newton increment twice. Apply once, then synchronize all dual-flow target fields.
    if( string( PRIMARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString() ) ==
        string( SECONDARY_FLOW_SOLVER::viewKeyStruct::elemDofFieldString() ) )
    {
      primarySolver()->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );

      this->forDiscretizationOnMeshTargets( domain.getMeshBodies(),
        [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
      {
        stdVector< string > syncFieldNames{ fields::flow::pressure::key() };
        if constexpr ( std::is_base_of_v< CompositionalMultiphaseBase, PRIMARY_FLOW_SOLVER > )
        {
          CompositionalMultiphaseFormulationType const formulationType =
            primarySolver()->template getReference< CompositionalMultiphaseFormulationType >(
              CompositionalMultiphaseBase::viewKeyStruct::formulationTypeString() );
          syncFieldNames.emplace_back( formulationType == CompositionalMultiphaseFormulationType::OverallComposition ?
                                       fields::flow::globalCompFraction::key() :
                                       fields::flow::globalCompDensity::key() );
        }

        FieldIdentifiers fieldsToBeSync;
        fieldsToBeSync.addElementFields( syncFieldNames, regionNames );
        CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
      } );
      return;
    }

    primarySolver()->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
    secondarySolver()->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
  }

  virtual void applyBoundaryConditions( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override
  {
    primarySolver()->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    secondarySolver()->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  /**
   * @brief Set up the linear system
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param setSparsity flag to indicate if the sparsity pattern should be set
   * @note TracAI: Overridden to register mesh connectivity
   */
  virtual void setupSystem(  DomainPartition & domain,
                           DofManager & dofManager,
                           CRSMatrix< real64, globalIndex > & localMatrix,
                           ParallelVector & rhs,
                           ParallelVector & solution,
                           bool const setSparsity = true ) override
  {
    // TracAI: Register mesh connectivity on first setup system call
    static bool connectivityRegistered = false;
    if( !connectivityRegistered )
    {
      //this->registerMeshConnectivity( domain );
      //TODO: temporarily disable automatic connectivity registration
      // the fracture and matrix region should be rigorously same
      connectivityRegistered = true;
    }

    // Call base class implementation
    Base::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
  }

  void updateGravityPressure(MeshLevel & meshMatrix, MeshLevel & meshFracture, real64 const & gravityCoefficient)
  {
    m_crossFlow.setupGravityDrainagePressure(meshMatrix, meshFracture, gravityCoefficient);
  }

protected:
  /**
   * @brief Print registered mesh connectivity values to terminal
   * @param domain the domain containing the mesh and fields
   * @note TracAI: Added to print registered connectivity values
   */
  void printRegisteredConnectivityValues( DomainPartition & domain )
  {
    GEOS_LOG( "=== Printing registered mesh connectivity values ===" );
    
    // Print matrix and fracture region pairs if defined
    if( !m_matrixRegionList.empty() && !m_fractureRegionList.empty() )
    {
      GEOS_LOG( "Matrix and fracture region pairs defined in input:" );
      if( m_matrixRegionList.size() == m_fractureRegionList.size() )
      {
        for( size_t i = 0; i < m_matrixRegionList.size(); ++i )
        {
          std::string matrixRegion = m_matrixRegionList[i];
          std::string fractureRegion = m_fractureRegionList[i];
          GEOS_LOG( "  Pair " << (i+1) << ": " << matrixRegion << " (matrix) <-> " << fractureRegion << " (fracture)" );
        }
      }
      else
      {
        GEOS_LOG( "  Invalid region lists (must have the same number of elements)" );
        GEOS_LOG( "  Matrix regions: " << m_matrixRegionList.size() );
        GEOS_LOG( "  Fracture regions: " << m_fractureRegionList.size() );
      }
    }
    
    // Iterate through all mesh bodies to check connectivity
    domain.forMeshBodies( [&]( MeshBody & meshBody )
    {
      MeshLevel & mesh = meshBody.getMeshLevels().getGroup< MeshLevel >( 0 );
      ElementRegionManager & elemManager = mesh.getGroup< ElementRegionManager >( "ElementRegions" );
      
      // Iterate through all element regions in the current mesh body
      elemManager.forElementRegions( [&]( ElementRegionBase & region )
      {
        if( auto * cellRegion = dynamic_cast< CellElementRegion * >( &region ) )
        {
          // Iterate through all element subregions in this region
          cellRegion->forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
          {
            if( auto * cellSubRegion = dynamic_cast< CellElementSubRegion * >( &subRegion ) )
            {
              // Check if the connectivity field exists for mesh1 to mesh2
              if( cellSubRegion->hasWrapper( viewKeyStruct::mesh1ToMesh2ConnectivityString() ) )
              {
                GEOS_LOG( "  " << meshBody.getName() << "/Region: " << cellRegion->getName() << "/SubRegion: " << cellSubRegion->getName() );
                GEOS_LOG( "  Has connectivity: mesh1ToMesh2Connectivity" );
                
                // Get the connectivity field
                auto const & connectivity = cellSubRegion->getReference< array1d< localIndex > >( viewKeyStruct::mesh1ToMesh2ConnectivityString() );
                
                // Print first few values to avoid cluttering the terminal
                localIndex printCount = std::min( localIndex( 5 ), cellSubRegion->size() );
                std::string values;
                for( localIndex i = 0; i < printCount; ++i )
                {
                  values += std::to_string( connectivity[i] ) + " ";
                }
                if( printCount < cellSubRegion->size() )
                {
                  values += "...";
                }
                GEOS_LOG( "  Connectivity values: " << values );
                GEOS_LOG( "  Total elements: " << cellSubRegion->size() );
              }
              
              // Check if the connectivity field exists for mesh2 to mesh1
              if( cellSubRegion->hasWrapper( viewKeyStruct::mesh2ToMesh1ConnectivityString() ) )
              {
                GEOS_LOG( "  " << meshBody.getName() << "/Region: " << cellRegion->getName() << "/SubRegion: " << cellSubRegion->getName() );
                GEOS_LOG( "  Has connectivity: mesh2ToMesh1Connectivity" );
                
                // Get the connectivity field
                auto const & connectivity = cellSubRegion->getReference< array1d< localIndex > >( viewKeyStruct::mesh2ToMesh1ConnectivityString() );
                
                // Print first few values to avoid cluttering the terminal
                localIndex printCount = std::min( localIndex( 5 ), cellSubRegion->size() );
                std::string values;
                for( localIndex i = 0; i < printCount; ++i )
                {
                  values += std::to_string( connectivity[i] ) + " ";
                }
                if( printCount < cellSubRegion->size() )
                {
                  values += "...";
                }
                GEOS_LOG( "  Connectivity values: " << values );
                GEOS_LOG( "  Total elements: " << cellSubRegion->size() );
              }
            }
          } );
        }
      } );
    } );
    
    // Analyze whether numbering is reset per-region or per-subregion
    // We inspect the connectivity arrays registered per subregion. If each subregion's
    // connectivity starts at zero we conclude that numbering is re-counted per-subregion.
    // If the first subregion of a region starts at 0 and subsequent subregions start
    // right after the previous max we conclude numbering is contiguous per-region.

    GEOS_LOG( "--- Analyzing numbering behavior (region vs subregion) ---" );

    // Helper to analyze one connectivity key
    auto analyzeKey = [&]( char const * key )
    {
      // group by mesh body and region
      std::map< std::pair< std::string, std::string >, std::vector< std::tuple< std::string, localIndex, localIndex > > > groups;

      domain.forMeshBodies( [&]( MeshBody & meshBody )
      {
        MeshLevel & mesh = meshBody.getMeshLevels().getGroup< MeshLevel >( 0 );
        ElementRegionManager & elemManager = mesh.getGroup< ElementRegionManager >( "ElementRegions" );

        elemManager.forElementRegions( [&]( ElementRegionBase & region )
        {
          if( auto * cellRegion = dynamic_cast< CellElementRegion * >( &region ) )
          {
            cellRegion->forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
            {
              if( auto * cellSubRegion = dynamic_cast< CellElementSubRegion * >( &subRegion ) )
              {
                if( cellSubRegion->hasWrapper( key ) )
                {
                  auto const & connectivity = cellSubRegion->getReference< array1d< localIndex > >( key );
                  if( connectivity.size() == 0 ) return;

                  localIndex minVal = connectivity[0];
                  localIndex maxVal = connectivity[ connectivity.size() - 1 ];
                  // record subregion name, min and max
                  groups[ { meshBody.getName(), cellRegion->getName() } ].push_back( std::make_tuple( cellSubRegion->getName(), minVal, maxVal ) );
                }
              }
            } );
          }
        } );
      } );

      if( groups.empty() )
      {
        GEOS_LOG( "No connectivity arrays found for key: " << key );
        return;
      }

      // Analyze each (meshBody,region) group
      for( auto const & kv : groups )
      {
        auto const & meshName = kv.first.first;
        auto const & regionName = kv.first.second;
        auto const & infos = kv.second;

        bool allStartAtZero = true;
        bool contiguousAcross = true;

        // check first at zero and contiguous
        localIndex prevMax = -1;
        for( size_t i = 0; i < infos.size(); ++i )
        {
          localIndex minVal = std::get<1>( infos[i] );
          localIndex maxVal = std::get<2>( infos[i] );
          if( minVal != 0 ) allStartAtZero = false;
          if( i == 0 )
          {
            prevMax = maxVal;
            if( minVal != 0 ) contiguousAcross = false; // first doesn't start at zero
          }
          else
          {
            if( minVal != prevMax + 1 ) contiguousAcross = false;
            prevMax = maxVal;
          }
        }

        // Print conclusion
        if( allStartAtZero )
        {
          GEOS_LOG( "Mesh '" << meshName << "' Region '" << regionName << "': 元素编号在每个 subregion 内从 0 重新计数（按 subregion 重置）。" );
        }
        else if( contiguousAcross )
        {
          GEOS_LOG( "Mesh '" << meshName << "' Region '" << regionName << "': 元素编号在 region 内连续（按 region 连续编号）。" );
        }
        else
        {
          GEOS_LOG( "Mesh '" << meshName << "' Region '" << regionName << "': 编号行为混合或不一致（既不是严格按 subregion 重置，也不是按 region 连续）。" );
        }

        // Print a short example for this region
        for( auto const & info : infos )
        {
          GEOS_LOG( "  Subregion '" << std::get<0>( info ) << "' => min: " << std::get<1>( info ) << ", max: " << std::get<2>( info ) );
        }
      }
    };

    // Analyze both directions if present
    analyzeKey( viewKeyStruct::mesh1ToMesh2ConnectivityString() );
    analyzeKey( viewKeyStruct::mesh2ToMesh1ConnectivityString() );

    GEOS_LOG( "=== End of connectivity values ===" );
  }

    template< typename TYPE_LIST,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
    real64 assemblyLaunch( MeshLevel & mesh,
                           DofManager const & dofManager,
                           string_array const & regionNames,
                           string const & materialNamesString,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           real64 const dt,
                           PARAMS && ... params )
    {
      GEOS_MARK_FUNCTION;
      //additional mesh manger is needed
      NodeManager const & nodeManager = mesh.getNodeManager();

      string const dofKey = dofManager.getKey( fields::flow::pressure::key() );
      arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

      real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( this->gravityVector() );

      KERNEL_WRAPPER kernelWrapper( dofNumber,
                                    dofManager.rankOffset(),
                                    localMatrix,
                                    localRhs,
                                    dt,
                                    gravityVectorData,
                                    std::forward< PARAMS >( params )... );
      return 0;
      /*
      return finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< >,
              TYPE_LIST >( mesh,
                           regionNames,
                           this->solidMechanicsSolver()->getDiscretizationName(),
                           materialNamesString,
                           kernelWrapper );*/
    }

    stabilization::StabilizationType m_stabilizationType;
  DualContinuumStencil & getStencil(){return m_crossFlow.getStencil();};

  // Accessor for gravity drainage flag
  int getGravityDrainageFlag() const { return m_crossFlow.m_gravityDrainageFlag; }

  // Added flag to indicate if gravity draiange should be applied
private:
  //std::shared_ptr< DualContinuumCrossFlow > m_crossFlow;
  DualContinuumCrossFlow m_crossFlow;
  real64 m_transferCoefficient;
  bool m_enableCrossStorageCorrection = true;  // sequential-only; disabled for FullyImplicit coupling
  
  // TracAI: Added to store dual continuum region pairs from XML
  string_array m_matrixRegionList;
  string_array m_fractureRegionList;
  
  struct viewKeyStruct : Base::viewKeyStruct
  {
    static constexpr char const * transferCoefficientString() { return "transferCoefficient"; }
    // TracAI: Added connectivity field keys
    static constexpr char const * mesh1ToMesh2ConnectivityString() { return "mesh1ToMesh2Connectivity"; }
    static constexpr char const * mesh2ToMesh1ConnectivityString() { return "mesh2ToMesh1Connectivity"; }
    // TracAI: Added key for dual continuum region pairs
    static constexpr char const * matrixRegionList() { return "matrixRegionList"; }
    static constexpr char const * fractureRegionList() { return "fractureRegionList"; }
    static constexpr char const * DualContinuumCrossFlow() { return "DualContinuumCrossFlow"; }

  };

}; // class DualContinuumFlowSolverBase

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALCONTINUUMFLOWSOLVER_HPP_
