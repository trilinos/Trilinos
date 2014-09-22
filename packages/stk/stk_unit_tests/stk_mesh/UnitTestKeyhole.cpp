/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/StkMeshIoBroker.hpp"
#include "unit_tests/SetupKeyholeMesh.hpp"

TEST(UnitTestKeyhole, NodeParts_case1)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case1(bulk);

  stk::mesh::Part& shared = meta.globally_shared_part();
  const stk::mesh::BucketVector& shared_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, shared);
  stk::mesh::PartVector blocksA(2);
  blocksA[0] = meta.get_part("block_1");
  blocksA[1] = meta.get_part("block_2");
  unsigned num_shared_nodes = 0;
  for(size_t i=0; i<shared_node_buckets.size(); ++i) {
    num_shared_nodes += shared_node_buckets[i]->size();
    const stk::mesh::Bucket& bucket = *shared_node_buckets[i];
    std::ostringstream oss;
    oss<<"proc "<<bulk.parallel_rank()<<", shared node ids: ";
    for(size_t j=0; j<bucket.size(); ++j) oss <<bulk.identifier(bucket[j])<<" ";
    std::cerr<<oss.str()<<std::endl;
    bool in_both_blocks = bucket.member_all(blocksA);
    EXPECT_TRUE(in_both_blocks);
  }

  const unsigned expected_num_shared_nodes = 2;
  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);

  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(8));
    stk::mesh::Entity node9 = bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(9));

    EXPECT_TRUE(bulk.is_valid(node8));
    EXPECT_TRUE(bulk.is_valid(node9));

    const stk::mesh::Part& aura_part = meta.aura_part();
    EXPECT_TRUE(bulk.bucket(node8).member(aura_part));
    EXPECT_TRUE(bulk.bucket(node9).member(aura_part));

    stk::mesh::Part& block_2 = *meta.get_part("block_2");
    stk::mesh::Part& block_3 = *meta.get_part("block_3");
    stk::mesh::PartVector blocksB(2);
    blocksB[0] = &block_2;
    blocksB[1] = &block_3;
    EXPECT_TRUE(bulk.bucket(node8).member_all(blocksB));
    EXPECT_TRUE(bulk.bucket(node9).member_all(blocksB));
  }
}

TEST(UnitTestKeyhole, NodeParts_case2)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case2(bulk);

  if (bulk.parallel_rank() == 0) {
    stk::mesh::Part& aura = meta.aura_part();
    const stk::mesh::BucketVector& aura_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, aura);
    stk::mesh::PartVector blocks(2);
    blocks[0] = meta.get_part("block_2");
    blocks[1] = meta.get_part("block_3");
    unsigned num_aura_nodes = 0;
    for(size_t i=0; i<aura_node_buckets.size(); ++i) {
      num_aura_nodes += aura_node_buckets[i]->size();
      const stk::mesh::Bucket& bucket = *aura_node_buckets[i];
      std::cerr<<"proc 0, aura node ids: ";
      for(size_t j=0; j<bucket.size(); ++j) std::cerr<<bulk.identifier(bucket[j])<<" ";
      std::cerr<<std::endl;
      bool in_both_blocks = bucket.member_all(blocks);
      EXPECT_TRUE(in_both_blocks);
    }
  
    const unsigned expected_num_aura_nodes = 2;
    EXPECT_EQ(expected_num_aura_nodes, num_aura_nodes);
  }
}

TEST(UnitTestKeyhole, EdgeParts_case1)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case1(bulk);

  stk::mesh::create_edges(bulk);

  //find the edge between nodes 2 and 3.
  stk::mesh::Entity edge = stk::mesh::Entity();
  stk::mesh::EntityId nodeId2 = 2;
  stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, nodeId2);
  unsigned num_edges = bulk.num_edges(node2);
  const stk::mesh::Entity* edges = bulk.begin_edges(node2);
  for(unsigned i=0; i<num_edges; ++i) {
    stk::mesh::Entity this_edge = edges[i];
    const stk::mesh::Entity* edge_nodes = bulk.begin_nodes(this_edge);
    if (bulk.identifier(edge_nodes[0])==2 && bulk.identifier(edge_nodes[1])==3) {
      edge = this_edge;
      break;
    }
  }

  EXPECT_TRUE(bulk.is_valid(edge));
  std::cerr<<"proc "<<bulk.parallel_rank()<<" found edge id="<<bulk.identifier(edge)<<" between nodes 2 and 3"<<std::endl;

  const stk::mesh::Part& block_1 = *meta.get_part("block_1");
  const stk::mesh::Part& block_2 = *meta.get_part("block_2");
  const stk::mesh::Bucket& edge_bucket = bulk.bucket(edge);
  EXPECT_TRUE(edge_bucket.member(block_1));
  EXPECT_FALSE(edge_bucket.member(block_2));

  //find the edge between nodes 8 and 9.
  edge = stk::mesh::Entity();
  stk::mesh::EntityId nodeId8 = 8;
  stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, nodeId8);
  num_edges = bulk.num_edges(node8);
  edges = bulk.begin_edges(node8);
  for(unsigned i=0; i<num_edges; ++i) {
    stk::mesh::Entity this_edge = edges[i];
    const stk::mesh::Entity* edge_nodes = bulk.begin_nodes(this_edge);
    if (bulk.identifier(edge_nodes[0])==8 && bulk.identifier(edge_nodes[1])==9) {
      edge = this_edge;
      break;
    }
  }

  EXPECT_TRUE(bulk.is_valid(edge));
  std::cerr<<"proc "<<bulk.parallel_rank()<<" found edge id="<<bulk.identifier(edge)<<" between nodes 8 and 9"<<std::endl;

  const stk::mesh::Part& block_3 = *meta.get_part("block_3");
  EXPECT_TRUE(bulk.bucket(edge).member(block_2));
  EXPECT_TRUE(bulk.bucket(edge).member(block_3));
}

TEST(UnitTestKeyhole, EdgeParts_case2)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case2(bulk);

  stk::mesh::create_edges(bulk);

  //find the edge between nodes 5 and 6.
  stk::mesh::Entity edge = stk::mesh::Entity();
  stk::mesh::EntityId nodeId5 = 5;
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, nodeId5);
  unsigned num_edges = bulk.num_edges(node5);
  const stk::mesh::Entity* edges = bulk.begin_edges(node5);
  for(unsigned i=0; i<num_edges; ++i) {
    stk::mesh::Entity this_edge = edges[i];
    const stk::mesh::Entity* edge_nodes = bulk.begin_nodes(this_edge);
    if (bulk.identifier(edge_nodes[0])==5 && bulk.identifier(edge_nodes[1])==6) {
      edge = this_edge;
      break;
    }
  }

  EXPECT_TRUE(bulk.is_valid(edge));
  std::cerr<<"proc "<<bulk.parallel_rank()<<" found edge id="<<bulk.identifier(edge)<<" between nodes 5 and 6"<<std::endl;

  const stk::mesh::Part& block_2 = *meta.get_part("block_2");
  const stk::mesh::Part& block_3 = *meta.get_part("block_3");
  const stk::mesh::Bucket& edge_bucket = bulk.bucket(edge);
  EXPECT_TRUE(edge_bucket.member(block_2));
  EXPECT_FALSE(edge_bucket.member(block_3));
}

