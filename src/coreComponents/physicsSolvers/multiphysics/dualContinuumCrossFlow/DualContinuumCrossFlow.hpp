#ifndef GEOS_DUALCONTINUUMCROSSFLOW_HPP
#define GEOS_DUALCONTINUUMCROSSFLOW_HPP

#include "dataRepository/Group.hpp"
#include "physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/DualContinuumStencil.hpp"

namespace geos
{
// Forward declarations
class MeshLevel;
class DomainPartition;
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
  // Setup: Read params and build stencil
  void setupCrossFlow( DomainPartition & domain,
                       MeshLevel & meshMatrix,
                       MeshLevel & meshFracture );

  void setupGravityDrainagePressure( MeshLevel & meshMatrix,
                                     MeshLevel & fractureMatrix,
                                     real64 const & gravityCoefficient);

  // Launch the kernel to add transfer terms to residual
  // (Template to support different solver traits if needed)
  template< typename SOLVER_TRAITS >
  void assembleCouplingTerms( MeshLevel & mesh,
                              DofManager const & dofManager,
                              typename SOLVER_TRAITS::PhysicsData & physicsData,
                              typename SOLVER_TRAITS::ResidualType & residual ) const;

  // Accessor for the stencil wrapper
  auto createKernelWrapper() const { return m_stencil.createKernelWrapper(); }

  // Accessor for fracture spacing
  real64 getFracSpacingLz() const { return m_fracSpacingLz; }

  static string catalogName() { return "DualContinuumCrossFlow"; }

  struct viewKeyStruct : public dataRepository::Group::viewKeyStruct
  {
    static constexpr char const * fractureSpacingLxString() { return "fractureSpacingLx"; }
    static constexpr char const * fractureSpacingLyString() { return "fractureSpacingLy"; }
    static constexpr char const * fractureSpacingLzString() { return "fractureSpacingLz"; }
    static constexpr char const * matrixRegionList() { return "matrixRegionList"; }
    static constexpr char const * fractureRegionList() { return "fractureRegionList"; }
    static constexpr char const * DualContinuumStencilString() {return "DualContinuumStencil";}
    static constexpr char const * gravityDrainageFlag() { return "gravityDrainageFlag"; }
    /// Direct interporosity exchange coefficient Gamma [Pa^{-1} s^{-1}].
    /// When > 0, bypasses the Kazemi shape-factor formula and uses
    /// transmissibility = Gamma * mu * V_element.
    static constexpr char const * interporosityExchangeCoefficientString()
    { return "interporosityExchangeCoefficient"; }
    /// Fracture (secondary continuum) volume fraction v_f; matrix v_m = 1 - v_f.
    /// Needed for the multi-porosity effective storage matrix M_bar.
    static constexpr char const * fractureVolumeFractionString()
    { return "fractureVolumeFraction"; }
    /// Intrinsic (true physical) matrix/fracture Biot coefficient and drained bulk modulus.
    /// When > 0 these are used to compute the multi-porosity storage M_bar in the FIM path,
    /// independently of the (possibly effective-medium) constitutive material values that the
    /// monolithic mechanics kernel consumes. <0 means "use the material value".
    static constexpr char const * intrinsicMatrixBiotString() { return "intrinsicMatrixBiot"; }
    static constexpr char const * intrinsicMatrixBulkModulusString() { return "intrinsicMatrixBulkModulus"; }
    static constexpr char const * intrinsicFractureBiotString() { return "intrinsicFractureBiot"; }
    static constexpr char const * intrinsicFractureBulkModulusString() { return "intrinsicFractureBulkModulus"; }
    /// Scale on the multi-porosity off-diagonal storage (1=paper bulk-Kbar value). Accounts for the
    /// incomplete monolithic-Schur cancellation under laterally confined (Mandel) geometry.
    static constexpr char const * crossStorageOffDiagScaleString() { return "crossStorageOffDiagScale"; }
  };

  DualContinuumStencil & getStencil(){return m_stencil;};
  int m_gravityDrainageFlag;

  /// Get the direct interporosity exchange coefficient
  real64 getInterporosityExchangeCoefficient() const
  { return m_interporosityExchangeCoefficient; }

  /// Get the fracture volume fraction (v_f); <0 means "unset"
  real64 getFractureVolumeFraction() const
  { return m_fractureVolumeFraction; }

  /// Intrinsic-parameter accessors for the FIM multi-porosity storage (<0 = use material value)
  real64 getIntrinsicMatrixBiot() const { return m_intrinsicMatrixBiot; }
  real64 getIntrinsicMatrixBulkModulus() const { return m_intrinsicMatrixBulkModulus; }
  real64 getIntrinsicFractureBiot() const { return m_intrinsicFractureBiot; }
  real64 getIntrinsicFractureBulkModulus() const { return m_intrinsicFractureBulkModulus; }
  real64 getCrossStorageOffDiagScale() const { return m_crossStorageOffDiagScale; }

  /// Setters used by the useIntrinsicInput path: the dual-poromechanics solver reads the
  /// intrinsic moduli/Biot off the constitutive models and pushes them here so the FIM
  /// multi-porosity storage (Step 4b) uses the intrinsic values while the mechanics kernel
  /// uses the homogenized effective ones.
  void setIntrinsicMatrixBiot( real64 const v ) { m_intrinsicMatrixBiot = v; }
  void setIntrinsicMatrixBulkModulus( real64 const v ) { m_intrinsicMatrixBulkModulus = v; }
  void setIntrinsicFractureBiot( real64 const v ) { m_intrinsicFractureBiot = v; }
  void setIntrinsicFractureBulkModulus( real64 const v ) { m_intrinsicFractureBulkModulus = v; }

private:
  real64 m_fracSpacingLx;
  real64 m_fracSpacingLy;
  real64 m_fracSpacingLz;

  string_array m_matrixRegionList;
  string_array m_fractureRegionList;

  /// Direct interporosity exchange coefficient [Pa^{-1} s^{-1}]; 0 = use Kazemi
  real64 m_interporosityExchangeCoefficient = 0.0;

  /// Fracture volume fraction v_f (-1 = unset -> cross-storage disabled)
  real64 m_fractureVolumeFraction = -1.0;

  /// Intrinsic matrix/fracture Biot + drained bulk modulus for the FIM M_bar storage
  /// (-1 = unset -> use the constitutive material value)
  real64 m_intrinsicMatrixBiot = -1.0;
  real64 m_intrinsicMatrixBulkModulus = -1.0;
  real64 m_intrinsicFractureBiot = -1.0;
  real64 m_intrinsicFractureBulkModulus = -1.0;

  /// Off-diagonal multi-porosity storage scale (1 = paper bulk-Kbar value)
  real64 m_crossStorageOffDiagScale = 1.0;

  DualContinuumStencil m_stencil;
};


} // namespace geos

#endif