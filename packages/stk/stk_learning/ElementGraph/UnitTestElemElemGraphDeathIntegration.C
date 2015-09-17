
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/ElemElemGraph.hpp>  // for process_killed_elements
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <vector>                       // for vector
#include "BulkDataElementGraphTester.hpp"
#include "ElementGraphTester.hpp"       // for ElemElemGraphTester
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityVector, PartVector
#include "stk_unit_test_utils/ioUtils.hpp"  // for fill_mesh_using_stk_io
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
namespace stk { namespace mesh { class Part; } }





namespace
{

unsigned count_elements(const stk::mesh::BulkData &bulk, stk::mesh::Part &activePart)
{
    std::vector<unsigned> countsPerEntityType;
    stk::mesh::count_entities(activePart, bulk, countsPerEntityType);
    return countsPerEntityType[stk::topology::ELEMENT_RANK];
}

void deactivateElement(stk::mesh::BulkData &bulk, stk::mesh::Entity elementToDeactivate, stk::mesh::Part *activePart)
{
    bulk.change_entity_parts(elementToDeactivate, {}, {activePart});
}

void deactivateElements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector elementsToDeactivate, stk::mesh::Part *activePart)
{
    for(stk::mesh::Entity element : elementsToDeactivate)
    {
        deactivateElement(bulk, element, activePart);
    }
}

// Test cases: (all cases must involve a "element death" starter)
//  * create neighbor then deactivate => "refinement"
//  * move processor boundary (change_entity_owner) => rebalance

TEST(ElemElemGraph, killAnElementAndUpdateGraphThenKillAnotherElement)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::Part &activePart = meta.declare_part("active");
        stk::mesh::Part &boundaryPart = meta.declare_part("boundary");
        BulkDataElementGraphTester bulk(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulk, comm);

        stk::unit_test_util::put_mesh_into_part(bulk, activePart);

        ElemElemGraphTester elemGraph(bulk);

        stk::mesh::EntityVector elementsDeactivated;
        int procRankToKillElementsOn = 0;
        int procRank = stk::parallel_machine_rank(comm);
        if(procRank == procRankToKillElementsOn)
        {
            int element1Id = 1;
            stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEM_RANK, element1Id);
            elementsDeactivated.push_back(element1);
        }

        bulk.modification_begin();
        deactivateElements(bulk, elementsDeactivated, &activePart);
        bulk.modification_end();

        //doh! so painful... bulk.update_exposed_boundary({element1}, activePart);
        stk::mesh::PartVector newlyCreatedBoundaryFacesGetAddedToTheseParts = { &boundaryPart };
        stk::mesh::PartVector exposedButExistingBoundaryFacesGetAddedToTheseParts = { &boundaryPart };
        process_killed_elements(bulk,
                                elemGraph,
                                elementsDeactivated,
                                activePart,
                                newlyCreatedBoundaryFacesGetAddedToTheseParts,
                                &exposedButExistingBoundaryFacesGetAddedToTheseParts);

        unsigned goldNumElements = 2u;
        if(procRank == procRankToKillElementsOn)
        {
            goldNumElements = 1u;
        }
        EXPECT_EQ(goldNumElements, count_elements(bulk, activePart));

//        add_new_element();
//        update_graph();

        elementsDeactivated.clear();

        if(procRank != procRankToKillElementsOn)
        {
            int element4Id = 4;
            stk::mesh::Entity element4 = bulk.get_entity(stk::topology::ELEM_RANK, element4Id);
            elementsDeactivated.push_back(element4);
        }

        bulk.modification_begin();
        deactivateElements(bulk, elementsDeactivated, &activePart);
        bulk.modification_end();

        process_killed_elements(bulk,
                                elemGraph,
                                elementsDeactivated,
                                activePart,
                                newlyCreatedBoundaryFacesGetAddedToTheseParts,
                                &exposedButExistingBoundaryFacesGetAddedToTheseParts);

        unsigned numElements = stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::ELEM_RANK));
        EXPECT_EQ(1u, numElements);

        unsigned goldNumFaces = 1;
        EXPECT_EQ(goldNumFaces, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::FACE_RANK)));
    }
}

} // namespace
