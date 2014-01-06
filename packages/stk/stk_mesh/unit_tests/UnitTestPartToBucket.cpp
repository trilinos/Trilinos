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

int procId(MPI_Comm comm)
{
  int procId = 1;
  MPI_Comm_rank(comm, &procId);
  return procId;
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

/*
       C 7--------------8 D
         |\             |\
         | \      s2    | \
         |  \           |  \          y
         | C 3--------------4 D    z  |
         |   |          |   |       \ |
         |   |          | s1|        \|
       A 5---|----------6 B |         o------x
          \  |           \  |
           \ |            \ |
            \|             \|
           A 1--------------2 B
*/

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
    size_t expectedNodeBuckets = 4; // A,B,C,D in above picture
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

/*
           C 11------------12 C
             |\             |\                  Locally Owned   Globally Shared      Ghosted
             | \            | \                     x  x             o  o             *  *
             |  \           |  \                 x        x       o        o       *        *
   proc 1    | B 7--------------8 B             x          x     o          o     *          *
             |   |`         |   |`              x          x     o          o     *          *
             |   | `        |   | `              x        x       o        o       *        *
           C 9---|---------10 C |  `                x  x             o  o             *  *
              \  |   `       \  |   `
               \ |    `       \ |    `
                \|     `       \|     `
               B 5--------------6 B    `
                  `      `       `      `
                   `    B 7--------------8 B
                    `     |\       `     |\                               Locally Owned        Ghosted
                     `    | \       `    | \                                  x  x              *  *
                      `   |  \       `   |  \            y                 x     o  x        *        *
             proc 0    `  | A 3--------------4 A      z  |                x     o GS x      *          *
                        ` |   |        ` |   |         \ |                x     o    x      *          *
                         `|   |         `|   |          \|                 x     o  x        *        *
                        B 5---|----------6 B |           o------x             x  x              *  *
                           \  |           \  |
                            \ |            \ |
                             \|             \|
                            A 1--------------2 A
*/

TEST(PartToBucket, np2TwoHex)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (numProcessors(communicator) != 2)
    {
        return;
    }
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

    const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    size_t expectedNodeBuckets = 3;
    ASSERT_EQ(expectedNodeBuckets, nodeBuckets.size());
    size_t expectedNumNodesInEveryBucket = 4;
    for(size_t i=0; i<nodeBuckets.size(); i++)
    {
        EXPECT_EQ(expectedNumNodesInEveryBucket, nodeBuckets[i]->size());
    }

    const stk::mesh::Part& locallyOwned = stkMeshMetaData.locally_owned_part();
    const stk::mesh::Part& globallyShared = stkMeshMetaData.globally_shared_part();
    stk::mesh::Selector locallyOwnedSelector(locallyOwned);
    stk::mesh::Selector globallySharedSelector(globallyShared);
    stk::mesh::Selector ghostedSelector(!(locallyOwned | globallyShared));
    const stk::mesh::BucketVector &locallyOwnedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, locallyOwnedSelector);
    const stk::mesh::BucketVector &globallySharedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, globallySharedSelector);
    const stk::mesh::BucketVector &ghostedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, ghostedSelector);

    size_t expectedNumLocallyOwnedNodeBuckets = 0;
    size_t expectedNumGloballySharedNodeBuckets = 0;
    size_t expectedNumGhostedNodeBuckets = 0;
    if(procId(communicator) == 0)
    {
        expectedNumLocallyOwnedNodeBuckets = 2; // A, B
        expectedNumGloballySharedNodeBuckets = 1; // B
        expectedNumGhostedNodeBuckets = 1; // C
    }
    else
    {
        expectedNumLocallyOwnedNodeBuckets = 1; // C
        expectedNumGloballySharedNodeBuckets = 1; // B
        expectedNumGhostedNodeBuckets = 1; // A
    }
    EXPECT_EQ(expectedNumLocallyOwnedNodeBuckets, locallyOwnedNodeBuckets.size());
    EXPECT_EQ(expectedNumGloballySharedNodeBuckets, globallySharedNodeBuckets.size());
    EXPECT_EQ(expectedNumGhostedNodeBuckets, ghostedNodeBuckets.size());


    const stk::mesh::BucketVector &locallyOwnedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, locallyOwnedSelector);
    const stk::mesh::BucketVector &globallySharedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, globallySharedSelector);
    const stk::mesh::BucketVector &ghostedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, ghostedSelector);
    size_t expectedNumLocallyOwnedElementBuckets = 1;
    size_t expectedNumGloballySharedElementBuckets = 0;
    size_t expectedNumGhostedElementBuckets = 1;
    EXPECT_EQ(expectedNumLocallyOwnedElementBuckets, locallyOwnedElementBuckets.size());
    EXPECT_EQ(expectedNumGloballySharedElementBuckets, globallySharedElementBuckets.size());
    EXPECT_EQ(expectedNumGhostedElementBuckets, ghostedElementBuckets.size());
}

}
