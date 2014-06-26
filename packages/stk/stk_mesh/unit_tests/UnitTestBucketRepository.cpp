#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition
#include <utility>                      // for pair
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for BucketVector, OrdinalVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc

TEST(BucketRepositoryTest, createBuckets)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    size_t spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim, stk::mesh::entity_rank_names());

    stk::mesh::OrdinalVector parts;
    parts.push_back(stkMeshMetaData.universal_part().mesh_meta_data_ordinal());
    parts.push_back(stkMeshMetaData.locally_owned_part().mesh_meta_data_ordinal());
    parts.push_back(stkMeshMetaData.declare_part("part1").mesh_meta_data_ordinal());
    parts.push_back(stkMeshMetaData.declare_part("part2").mesh_meta_data_ordinal());
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
