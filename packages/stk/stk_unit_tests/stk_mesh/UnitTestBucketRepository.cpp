// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <utility>                      // for pair
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
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

    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, comm);
    stk::mesh::impl::EntityRepository entityRepository;

    stk::mesh::impl::BucketRepository &bucketRepository = stkMeshBulkData.my_get_bucket_repository();
    stk::mesh::impl::Partition* partition = bucketRepository.get_or_create_partition(stk::topology::NODE_RANK, parts);

    size_t numNodes = 1024;
    for(size_t i=0; i<numNodes; i++)
    {
        stk::mesh::EntityId nodeID = i+1;
        stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK, nodeID);
        std::pair<stk::mesh::entity_iterator,bool> createResult = entityRepository.internal_create_entity(nodeKey);
        bool aNewEntityWasCreated = createResult.second;
        EXPECT_TRUE(aNewEntityWasCreated);
        stk::mesh::Entity node = stkMeshBulkData.my_generate_new_entity();
        stkMeshBulkData.my_set_entity_key(node, nodeKey);
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
    std::pair<stk::mesh::entity_iterator,bool> createResult = entityRepository.internal_create_entity(nodeKey);
    bool aNewEntityWasCreated = createResult.second;
    EXPECT_TRUE(aNewEntityWasCreated);
    stk::mesh::Entity node = stkMeshBulkData.my_generate_new_entity();
    stkMeshBulkData.my_set_entity_key(node, nodeKey);
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
