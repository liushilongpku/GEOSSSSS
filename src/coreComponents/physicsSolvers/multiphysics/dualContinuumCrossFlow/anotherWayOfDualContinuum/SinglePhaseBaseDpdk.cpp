//
// Created by hello on 2025/12/22.
// This File is not in unused.
//
#include "SinglePhaseBaseDpdk.hpp"



#include "common/TimingMacros.hpp"

#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"


namespace geos
{

    using namespace dataRepository;
    using namespace constitutive;
    using namespace fields;
    using namespace singlePhaseBaseKernels;

SinglePhaseBaseDpdk::SinglePhaseBaseDpdk( const string & name,
                                  Group * const parent ):
 SinglePhaseBase( name, parent )
{}


void SinglePhaseBaseDpdk::assembleAccumulationTerms( DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
    GEOS_MARK_FUNCTION;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
        mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                SurfaceElementSubRegion >( regionNames,
                                           [&]( localIndex const,
                                                auto & subRegion )
                                           {
                                               accumulationAssemblyLaunchDPDK( dofManager, subRegion, localMatrix, localRhs );
                                           } );
    } );
}

void SinglePhaseBaseDpdk::applyBoundaryConditions( real64 time_n,
                                               real64 dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
    GEOS_MARK_FUNCTION;

    if( m_keepVariablesConstantDuringInitStep )
    {
        // this function is going to force the current flow state to be constant during the time step
        // this is used when the poromechanics solver is performing the stress initialization
        // TODO: in the future, a dedicated poromechanics kernel should eliminate the flow vars to construct a reduced system
        //       which will remove the need for this brittle passing aroung of flag
        keepVariablesConstantDuringInitStep( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
    }
    else
    {
        applySourceFluxBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
        applyDirichletBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
        applyAquiferBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
    }
}
namespace
{

    void applyAndSpecifyFieldValue( real64 const & time_n,
                                    real64 const & dt,
                                    MeshLevel & mesh,
                                    globalIndex const rankOffset,
                                    string const dofKey,
                                    bool const,
                                    integer const idof,
                                    string const fieldKey,
                                    string const boundaryFieldKey,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
    {
        FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

        fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                                 mesh,
                                                 fieldKey,
                                                 [&]( FieldSpecificationBase const & fs,
                                                      string const &,
                                                      SortedArrayView< localIndex const > const & lset,
                                                      ElementSubRegionBase & subRegion,
                                                      string const & )
                                                 {
                                                     // Specify the bc value of the field
                                                     fs.applyFieldValue< FieldSpecificationEqual,
                                                             parallelDevicePolicy<> >( lset,
                                                                                       time_n + dt,
                                                                                       subRegion,
                                                                                       boundaryFieldKey );

                                                     arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
                                                     arrayView1d< globalIndex const > const dofNumber =
                                                             subRegion.getReference< array1d< globalIndex > >( dofKey );
                                                     arrayView1d< real64 const > const bcField =
                                                             subRegion.getReference< array1d< real64 > >( boundaryFieldKey );
                                                     arrayView1d< real64 const > const field =
                                                             subRegion.getReference< array1d< real64 > >( fieldKey );

                                                     forAll< parallelDevicePolicy<> >( lset.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
                                                     {
                                                         localIndex const ei = lset[a];
                                                         if( ghostRank[ei] >= 0 )
                                                         {
                                                             return;
                                                         }

                                                         globalIndex const dofIndex = dofNumber[ei];
                                                         localIndex const localRow = dofIndex - rankOffset;
                                                         real64 rhsValue;

                                                         // Apply field value to the matrix/rhs
                                                         FieldSpecificationEqual::SpecifyFieldValue( dofIndex + idof,
                                                                                                     rankOffset,
                                                                                                     localMatrix,
                                                                                                     rhsValue,
                                                                                                     bcField[ei],
                                                                                                     field[ei] );
                                                         localRhs[localRow + idof] = rhsValue;
                                                     } );
                                                 } );
    }

}

void SinglePhaseBaseDpdk::applyDirichletBC( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )  const {
    GEOS_MARK_FUNCTION;

    string const dofKey = dofManager.getKey(viewKeyStruct::elemDofFieldString());
    globalIndex const rankOffset = dofManager.rankOffset();
    bool const isFirstNonlinearIteration = (m_nonlinearSolverParameters.m_numNewtonIterations == 0);

    forDiscretizationOnMeshTargets(domain.getMeshBodies(), [&](string const &,
                                                               MeshLevel &mesh,
                                                               string_array const &) {
        applyAndSpecifyFieldValue(time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration,
                                                   0, flow::pressure::key(), flow::bcPressure::key(),
                                                   localMatrix, localRhs);
        if (m_isThermal) {
            applyAndSpecifyFieldValue(time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration,
                                                       1, flow::temperature::key(), flow::bcTemperature::key(),
                                                       localMatrix, localRhs);
        }
    });
}



}/* namespace geos */