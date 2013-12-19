#include <gtest/gtest.h>
#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_mesh/base/BulkData.hpp>

TEST(BucketRepositoryTest, createBuckets)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    size_t spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim, stk::mesh::entity_rank_names());

    stk::mesh::PartVector parts;
    parts.push_back(&stkMeshMetaData.universal_part());
    parts.push_back(&stkMeshMetaData.locally_owned_part());
    parts.push_back(&stkMeshMetaData.declare_part("part1"));
    parts.push_back(&stkMeshMetaData.declare_part("part2"));
    stkMeshMetaData.commit();

    stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, comm);
    stk::mesh::impl::EntityRepository entityRepository(stkMeshBulkData);

    stk::mesh::impl::BucketRepository &bucketRepository =
      stk::mesh::impl::Partition::getRepository(stkMeshBulkData);

    stk::mesh::impl::Partition* partition = bucketRepository.get_or_create_partition(stk::topology::NODE_RANK, parts);

    size_t numNodes = 1024;
    for(size_t i=0; i<numNodes; i++)
    {
        stk::mesh::EntityId nodeID = i+1;
        stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK, nodeID);
        std::pair<stk::mesh::Entity,bool> createResult = entityRepository.internal_create_entity(nodeKey);
        stk::mesh::Entity node = createResult.first;
        bool aNewEntityWasCreated = createResult.second;
        EXPECT_TRUE(aNewEntityWasCreated);
        partition->add(node);
    }

    size_t expectedNumBuckets = 2;
    EXPECT_EQ(expectedNumBuckets, partition->num_buckets());

    const stk::mesh::BucketVector & nodeBuckets = bucketRepository.buckets(stk::topology::NODE_RANK);
    EXPECT_EQ(expectedNumBuckets, nodeBuckets.size());

    size_t expectedBucketSize = 512;
    EXPECT_EQ(expectedBucketSize, nodeBuckets[0]->size());
    EXPECT_EQ(expectedBucketSize, nodeBuckets[1]->size());

    stk::mesh::EntityId nodeID = numNodes+1;
    stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK, nodeID);
    std::pair<stk::mesh::Entity,bool> createResult = entityRepository.internal_create_entity(nodeKey);
    stk::mesh::Entity node = createResult.first;
    bool aNewEntityWasCreated = createResult.second;
    EXPECT_TRUE(aNewEntityWasCreated);
    partition->add(node);

    expectedNumBuckets = 3;
    const stk::mesh::BucketVector & newNodeBuckets = bucketRepository.buckets(stk::topology::NODE_RANK);
    EXPECT_EQ(expectedNumBuckets, newNodeBuckets.size());
    EXPECT_EQ(expectedNumBuckets, partition->num_buckets());

    EXPECT_EQ(expectedBucketSize, nodeBuckets[0]->size());
    EXPECT_EQ(expectedBucketSize, nodeBuckets[1]->size());
    expectedBucketSize = 1;
    EXPECT_EQ(expectedBucketSize, nodeBuckets[2]->size());

    size_t expectedBucketCapacity = 512;
    EXPECT_EQ(expectedBucketCapacity, nodeBuckets[0]->capacity());
    EXPECT_EQ(expectedBucketCapacity, nodeBuckets[1]->capacity());
    EXPECT_EQ(expectedBucketCapacity, nodeBuckets[2]->capacity());

}
