// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"

TEST(UnitTestGhosting, WithSharedFiltered)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x6";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();
  
  const stk::mesh::BucketVector& shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, shared_selector);

  size_t num_shared_nodes = 0;
  for(size_t i=0; i<shared_buckets.size(); ++i) {
    num_shared_nodes += shared_buckets[i]->size();
  }

  const size_t expected_num_shared_nodes = 4;
  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);

  stkMeshBulkData.modification_begin();

  //we will attempt to ghost shared-and-owned entities from proc 0 to proc 1.
  //they won't be allowed to be in the ghosting because they are already
  //shared on proc 1.
  int myProc = stkMeshBulkData.parallel_rank();
  int otherProc = 1;
  if (myProc != 0) otherProc = -1;//don't do ghosting if we are proc 1

  stk::mesh::Ghosting& custom_shared_ghosting = stkMeshBulkData.create_ghosting("custom_shared");

  std::vector<stk::mesh::EntityProc> send_shared;
  if (otherProc != -1) {
    for(size_t i=0; i<shared_buckets.size(); ++i) {
      const stk::mesh::Bucket& bucket = *shared_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
        if (bucket.parallel_owner_rank(j) == myProc) {
          send_shared.push_back(std::make_pair(bucket[j], otherProc));
        }
      }
    }
  }

  stkMeshBulkData.change_ghosting(custom_shared_ghosting, send_shared);

  stkMeshBulkData.modification_end();

  std::vector<stk::mesh::EntityProc> custom_shared_send_list;
  custom_shared_ghosting.send_list(custom_shared_send_list);

  size_t expected_send_list_size = 0;
  EXPECT_EQ(expected_send_list_size, custom_shared_send_list.size());
}

TEST(UnitTestGhosting, WithShared)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 3) {
    return;
  }

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x6";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();
  
  const stk::mesh::BucketVector& shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, shared_selector);

  size_t num_shared_nodes = 0;
  for(size_t i=0; i<shared_buckets.size(); ++i) {
    num_shared_nodes += shared_buckets[i]->size();
  }

  int myProc = stkMeshBulkData.parallel_rank();

  size_t expected_num_shared_nodes = 4;
  if (myProc == 1) expected_num_shared_nodes = 8;

  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);

  stkMeshBulkData.modification_begin();

  //we will ghost shared-and-owned entities from proc 0 to proc 2.
  int otherProc = 2;
  if (myProc != 0) otherProc = -1;//don't do ghosting

  stk::mesh::Ghosting& custom_shared_ghosting = stkMeshBulkData.create_ghosting("custom_shared");

  std::vector<stk::mesh::EntityProc> send_shared;
  if (otherProc != -1) {
    for(size_t i=0; i<shared_buckets.size(); ++i) {
      const stk::mesh::Bucket& bucket = *shared_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
        if (bucket.parallel_owner_rank(j) == myProc) {
          send_shared.push_back(std::make_pair(bucket[j], otherProc));
        }
      }
    }
  }

  stkMeshBulkData.change_ghosting(custom_shared_ghosting, send_shared);

  stkMeshBulkData.modification_end();

  std::vector<stk::mesh::EntityProc> custom_shared_send_list;
  custom_shared_ghosting.send_list(custom_shared_send_list);

  EXPECT_EQ(send_shared.size(), custom_shared_send_list.size());
  if (myProc == 0) {
    EXPECT_EQ(num_shared_nodes, custom_shared_send_list.size());
  }
}

