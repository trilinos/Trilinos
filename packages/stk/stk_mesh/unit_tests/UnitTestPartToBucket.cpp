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

TEST(PartToBucket, simple)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:2x2x2";
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

TEST(PartToBucket, simpleWithSingleSideset)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    const stk::mesh::BucketVector &boundaryNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, boundaryNodesSelector);
    size_t expectedBoundaryNodeBuckets = 1;
    EXPECT_EQ(expectedBoundaryNodeBuckets, boundaryNodeBuckets.size());

    size_t numNodesInBoundary = 9;
    EXPECT_EQ(numNodesInBoundary, boundaryNodeBuckets[0]->size());

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

    size_t numFacesInBucket = 4;
    EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

TEST(PartToBucket, simpleWithTwoSidesets)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:XY";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    const stk::mesh::BucketVector &boundaryNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, boundaryNodesSelector);
    size_t expectedBoundaryNodeBuckets = 2;
    EXPECT_EQ(expectedBoundaryNodeBuckets, boundaryNodeBuckets.size());

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

    size_t numFacesInBucket = 4;
    EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

TEST(PartToBucket, simpleWithSixSidesets)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:xXyYzZ";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector surface1Selector(surface1Part);

    const stk::mesh::BucketVector &boundaryNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1Selector);
    size_t expectedBoundaryNodeBuckets = 9;
    EXPECT_EQ(expectedBoundaryNodeBuckets, boundaryNodeBuckets.size());

    stk::mesh::Part &surface3Part = *stkMeshMetaData.get_part("surface_3");
    stk::mesh::Selector selectTwoAdjacentSurfaces = surface1Part | surface3Part;

    const stk::mesh::BucketVector &twoAdjacentSurfacesNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, selectTwoAdjacentSurfaces);
    size_t expectedBucketsForTwoAdjacentSurfaces = 15;
    EXPECT_EQ(expectedBucketsForTwoAdjacentSurfaces, twoAdjacentSurfacesNodeBuckets.size());

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 27;
    EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

    const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    size_t expectedElemBuckets = 1;
    EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

    const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
    size_t expectedEdgeBuckets = 0;
    EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

    const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
    size_t expectedFaceBuckets = 6;
    EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());
}

}
