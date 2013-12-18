#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace {

int numProcessors(MPI_Comm comm)
{
  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);
  return numProcs;
}

void checkNodeInSelectedBucket(stk::mesh::Selector selectNode, stk::mesh::EntityId expectedGlobalId, stk::mesh::BulkData &stkMeshBulkData)
{
    const stk::mesh::BucketVector &node1BucketVector = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, selectNode);
    size_t expectedNumBucketsForAnyNode = 1;
    EXPECT_EQ(expectedNumBucketsForAnyNode, node1BucketVector.size());
    size_t expectedNumNodesPerBucket = 1;
    EXPECT_EQ(expectedNumNodesPerBucket, node1BucketVector[0]->size());
    const stk::mesh::Bucket &node1Bucket = *node1BucketVector[0];
    EXPECT_EQ(expectedGlobalId, stkMeshBulkData.identifier(node1Bucket[0]));
}
/*
         7--------------8
         |\             |\
         | \      s2    | \
         |  \           |  \          y
         |   3--------------4      z  |
         |s1 |          |   |       \ |
         |   |          |   |        \|
         5---|----------6   |         o------x
          \  |     s3    \  |
           \ |            \ |
            \|             \|
             1--------------2
*/
TEST(PartToBucket, hexWithThreeSidesets)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (numProcessors(communicator) != 1)
    {
        return;
    }
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);

    //generated-mesh 'sideset:xYz' syntax adds face surfaces on 3 sides of the mesh,
    //(minimum 'x' side, maximum 'y' side, minimum 'z' side)
    //and the IO system will create corresponding stk::mesh::Parts named 'surface_1',
    //'surface_2' and 'surface_3', respectively, which are referenced in code below.
    const std::string generatedMeshSpecification = "generated:1x1x1|sideset:xYz";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Part &surface2Part = *stkMeshMetaData.get_part("surface_2");
    stk::mesh::Part &surface3Part = *stkMeshMetaData.get_part("surface_3");

    stk::mesh::Selector surface1Selector(surface1Part);
    const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1Selector);
    size_t expectedSurface1NodeBuckets = 4;
    EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

    stk::mesh::Selector selectNode1 = surface1Part & !surface2Part & surface3Part;
    stk::mesh::EntityId expectedGlobalId = 1;
    checkNodeInSelectedBucket(selectNode1, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode2 = !surface1Part & !surface2Part & surface3Part;
    expectedGlobalId = 2;
    checkNodeInSelectedBucket(selectNode2, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode3 = surface1Part & surface2Part & surface3Part;
    expectedGlobalId = 3;
    checkNodeInSelectedBucket(selectNode3, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode4 = !surface1Part & surface2Part & surface3Part;
    expectedGlobalId = 4;
    checkNodeInSelectedBucket(selectNode4, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode5 = surface1Part & !surface2Part & !surface3Part;
    expectedGlobalId = 5;
    checkNodeInSelectedBucket(selectNode5, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Part &block1Part = *stkMeshMetaData.get_part("block_1");
    stk::mesh::Selector selectNode6 = block1Part & !surface1Part & !surface2Part & !surface3Part;
    expectedGlobalId = 6;
    checkNodeInSelectedBucket(selectNode6, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode7 = surface1Part & surface2Part & !surface3Part;
    expectedGlobalId = 7;
    checkNodeInSelectedBucket(selectNode7, expectedGlobalId, stkMeshBulkData);

    stk::mesh::Selector selectNode8 = !surface1Part & surface2Part & !surface3Part;
    expectedGlobalId = 8;
    checkNodeInSelectedBucket(selectNode8, expectedGlobalId, stkMeshBulkData);

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 8;
    EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

    const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    size_t expectedElemBuckets = 1;
    EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

    const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
    size_t expectedEdgeBuckets = 0;
    EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

    const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
    size_t expectedFaceBuckets = 3;
    EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());
}

}
