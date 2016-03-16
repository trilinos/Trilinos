#include "stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>       // for comm_mesh_counts, count_entities
#include <stk_mesh/base/CreateFaces.hpp>       // for comm_mesh_counts, count_entities
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "unit_tests/Setup2Block2HexMesh.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>


namespace
{

TEST(ForEntityFunctionInMeshImplUtils, test_counting_nodes_using_raw_bucket_loops)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        stk::mesh::BulkData bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        {
            unsigned numNodes = 0;
            stk::mesh::impl::for_each_selected_entity_run_no_threads(bulkData, stk::topology::NODE_RANK, metaData.universal_part(),
                [&numNodes](stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
                {
                    stk::mesh::Entity entity = stk::mesh::impl::get_entity(meshIndex);
                    if(mesh.is_valid(entity))
                    {
                        numNodes++;
                    }
                }
            );
            EXPECT_EQ(16u, numNodes);
        }
        {
            unsigned numNodes = 0;
            stk::mesh::impl::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK,
                [&numNodes](const stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
                {
                    stk::mesh::Entity entity = stk::mesh::impl::get_entity(meshIndex);
                    if(mesh.is_valid(entity))
                    {
                        numNodes++;
                    }
                }
            );
            EXPECT_EQ(16u, numNodes);
        }
        {
            const unsigned numNodes = stk::mesh::count_selected_entities(metaData.universal_part(),
                                                                         bulkData.buckets(stk::topology::NODE_RANK));
            ASSERT_EQ(16u, numNodes);
            std::vector<unsigned> localIds(numNodes);
            unsigned index = 0;
            stk::mesh::impl::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK,
                [&index](stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
                {
                    stk::mesh::Entity node = stk::mesh::impl::get_entity(meshIndex);
                    if(mesh.is_valid(node))
                    {
                        mesh.set_local_id(node, index);
                        index++;
                    }
                }
            );
            stk::mesh::impl::for_each_entity_run(bulkData, stk::topology::NODE_RANK,
                [&localIds](const stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
                {
                    stk::mesh::Entity node = stk::mesh::impl::get_entity(meshIndex);
                    if(mesh.is_valid(node))
                    {
                        unsigned localId = mesh.local_id(node);
                        localIds[localId] = localId;
                    }
                }
            );
            for(size_t i=0; i<numNodes; i++)
            {
                EXPECT_EQ(i, localIds[i]);
            }
        }
    }
}

}
