/*
 * ------------------------------------------------------------------------------------------------------------
lsl
2025.11.17
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALPOROSITYDUALPERM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALPOROSITYDUALPERM_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"

namespace geos
{

template< typename FLOW_SOLVER0 = SinglePhaseBase, typename FLOW_SOLVER1 = SinglePhaseBase >
class SinglePhasePoromechanics : public PoromechanicsSolver< FLOW_SOLVER0, FLOW_SOLVER1 >
{
public:

  using Base = PoromechanicsSolver< FLOW_SOLVER0, FLOW_SOLVER1 >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_stabilizationType;
  using Base::m_stabilizationRegionNames;
  using Base::m_stabilizationMultiplier;
  using Base::updateBulkDensity;

  /**
   * @brief main constructor for SinglePhasePoromechanics objects
   * @param name the name of this instantiation of SinglePhasePoromechanics in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanics
   */
  SinglePhasePoromechanics( const string & name,
                            dataRepository::Group * const parent );

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanics object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseBase > ) // special case
    {
      return "SinglePhasePoromechanics";
    }
    else // default
    {
      return FLOW_SOLVER::catalogName() + "Poromechanics";
    }
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void postInputInitialization() override;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  void assembleElementBasedTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  /**@}*/

  struct viewKeyStruct : Base::viewKeyStruct
  {
    /// Flag to indicate if it is the phase-field formulation
    constexpr static char const * damageFlagString() { return "damageFlag"; }
  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  integer m_damageFlag;

  virtual void setMGRStrategy()
  {
    if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}", this->getName(), this->getCatalogName()));
  }

  /**
   * @brief Helper function to recompute the bulk density
   * @param[in] subRegion the element subRegion
   */
  virtual void updateBulkDensity( ElementSubRegionBase & subRegion ) override;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_DUALPOROSITYDUALPERM_HPP_ */
