#ifndef GEOS_DUALCONTINUUMCROSSFLOW_HPP
#define GEOS_DUALCONTINUUMCROSSFLOW_HPP

#include "dataRepository/Group.hpp"
#include "finiteVolume/DualContinuumStencil.hpp"

namespace geos
{
// Forward declarations
class MeshLevel;
class DofManager;

/**
 * @brief Manages the coupling between matrix and fracture regions.
 * Handles parameter reading, stencil building, and kernel launching.
 */
class DualContinuumCrossFlow : public dataRepository::Group
{
public:
  // Constructor
  DualContinuumCrossFlow( string const & name, Group * parent );

  // Find region index in the provided list
  localIndex findRegionIndexInList( string const & regionName );

  // Find region index in the ElementRegionManager
  localIndex findRegionIndexInRegionManager(ElementRegionManager const & elemManager,
                                            string const & regionName );
  // Initialize: Read params and build stencil
  void initialize( MeshLevel & meshMatrix,
                   MeshLevel & meshFracture );

  // Launch the kernel to add transfer terms to residual
  // (Template to support different solver traits if needed)
  template< typename SOLVER_TRAITS >
  void assembleCouplingTerms( MeshLevel & mesh,
                              DofManager const & dofManager,
                              typename SOLVER_TRAITS::PhysicsData & physicsData,
                              typename SOLVER_TRAITS::ResidualType & residual ) const;

  // Accessor for the stencil wrapper
  auto createKernelWrapper() const { return m_stencil.createKernelWrapper(); }

  static string catalogName() { return "DualContinuumCrossFlow"; }

  struct viewKeyStruct : public dataRepository::Group::viewKeyStruct
  {
    // 定义的 XML 参数名
    static constexpr char const * fractureSpacingLxString() { return "fractureSpacingLx"; }
    static constexpr char const * fractureSpacingLyString() { return "fractureSpacingLy"; }
    static constexpr char const * fractureSpacingLzString() { return "fractureSpacingLz"; }
    static constexpr char const * matrixRegionList() { return "matrixRegionList"; }
    static constexpr char const * fractureRegionList() { return "fractureRegionList"; }
    static constexpr char const * DualContinuumStencilString() {return "DualContinuumStencil";}
  };
private:
  // --- Parameters ---
  real64 m_fracSpacingLx;
  real64 m_fracSpacingLy;
  real64 m_fracSpacingLz;

  string_array m_matrixRegionList;
  string_array m_fractureRegionList;

  // --- The Stencil ---
  DualContinuumStencil m_stencil;
};


} // namespace geos

#endif