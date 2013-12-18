#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace
{

int numProcessors(MPI_Comm comm)
{
  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);
  return numProcs;
}

TEST(PartToBucket, hex)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (numProcessors(communicator) != 1)
    {
        return;
    }
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x1";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 1;
    EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

    const stk::mesh::Part& locallyOwned = stkMeshMetaData.locally_owned_part();
    const stk::mesh::Part& globallyShared = stkMeshMetaData.globally_shared_part();
    const stk::mesh::Part& universal = stkMeshMetaData.universal_part();

    const stk::mesh::Bucket& nodeBucket = *nodeBuckets[0];

    EXPECT_TRUE(nodeBucket.member(locallyOwned));
    EXPECT_TRUE(nodeBucket.member(universal));
    EXPECT_FALSE(nodeBucket.member(globallyShared));

    const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    size_t expectedElemBuckets = 1;
    EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

    const stk::mesh::Bucket& elemBucket = *elemBuckets[0];

    EXPECT_TRUE(elemBucket.member(locallyOwned));
    EXPECT_TRUE(elemBucket.member(universal));
    EXPECT_FALSE(elemBucket.member(globallyShared));

    const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
    size_t expectedEdgeBuckets = 0;
    EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

    const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
    size_t expectedFaceBuckets = 0;
    EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());
}

TEST(PartToBucket, hexWithSingleSideset)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (numProcessors(communicator) != 1)
    {
        return;
    }
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x1|sideset:X";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector surface1NodesSelector(surface1Part);

    const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1NodesSelector);
    size_t expectedSurface1NodeBuckets = 1;
    EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

    size_t numNodesInBoundary = 4;
    EXPECT_EQ(numNodesInBoundary, surface1NodeBuckets[0]->size());

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 2;
    EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

    const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    size_t expectedElemBuckets = 1;
    EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

    const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
    size_t expectedEdgeBuckets = 0;
    EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

    const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
    size_t expectedFaceBuckets = 1;
    EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());

    size_t numFacesInBucket = 1;
    EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

TEST(PartToBucket, hexWithTwoSidesets)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (numProcessors(communicator) != 1)
    {
        return;
    }
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x1|sideset:XY";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector surface1NodesSelector(surface1Part);

    const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1NodesSelector);
    size_t expectedSurface1NodeBuckets = 2;
    EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 4;
    EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

    const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    size_t expectedElemBuckets = 1;
    EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

    const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
    size_t expectedEdgeBuckets = 0;
    EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

    const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
    size_t expectedFaceBuckets = 2;
    EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());

    size_t numFacesInBucket = 1;
    EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

}
