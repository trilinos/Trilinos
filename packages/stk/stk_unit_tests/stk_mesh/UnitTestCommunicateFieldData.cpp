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

//TEST(CommunicateFieldData, pack_ghost)
//{
//  stk::ParallelMachine communicator = MPI_COMM_WORLD;
//
//  int numProcs = stk::parallel_machine_size(communicator);
//  if (numProcs != 2) {
//    return;
//  }
//  int myProc = stk::parallel_machine_rank(communicator);
//
//  const unsigned spatialDim = 2;
//  stk::mesh::MetaData meta(spatialDim);
//  
//  stk::mesh::Field<double>& field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
//  stk::mesh::put_field(field, meta.universal_part());
//
//  stk::mesh::BulkData bulk(meta, communicator);
//
//  setupKeyholeMesh2D_case1(bulk);
//
//  stk::mesh::Part& aura_part = meta.aura_part();
//  const stk::mesh::BucketVector& aura_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, aura_part);
//
//  unsigned num_aura_nodes = 0;
//  for(size_t i=0; i<aura_node_buckets.size(); ++i) {
//    num_aura_nodes += aura_node_buckets[i]->size();
//    stk::mesh::Bucket& bucket = *aura_node_buckets[i];
//    for(size_t j=0; j<bucket.size(); ++j) {
//        stk::mesh::Entity node = bucket[j];
//        double value = myProc*100 + bulk.identifier(node);
//        double* data = stk::mesh::field_data(field, node);
//        *data = value;
//    }
//  }
//
//  unsigned expected_num_aura_nodes = 6;
//  if (myProc == 1) {
//      expected_num_aura_nodes = 2;
//  }
//  EXPECT_EQ(expected_num_aura_nodes, num_aura_nodes);
//
//  const stk::mesh::Ghosting& aura_ghosting = *bulk.ghostings()[stk::mesh::BulkData::AURA];
//
//  stk::mesh::FieldBase* fieldPtr = &field;
//  std::vector<std::vector<unsigned char> > send_data;
//  std::vector<std::vector<unsigned char> > recv_data;
//  stk::mesh::pack_ghost_field_data(bulk, aura_ghosting, 1, &fieldPtr, send_data, recv_data);
//
//  const stk::mesh::VolatileFastGhostCommMapOneRank& ghost_comm_map = bulk.volatile_fast_ghost_comm_map(stk::topology::NODE_RANK);
//  std::vector<unsigned> send_offsets(numProcs, 0);
//  std::vector<unsigned> recv_offsets(numProcs, 0);
//
//  const stk::mesh::BucketVector& all_node_buckets = bulk.buckets(stk::topology::NODE_RANK);
//
//  for(int proc=0; proc<numProcs; ++proc) {
//      for(size_t idata=0; idata<ghost_comm_map[proc].size(); ++idata) {
//          const unsigned ghost_id = ghost_comm_map[proc][idata].ghost_id;
//          const int owner       = ghost_comm_map[proc][idata].owner;
//          const unsigned bucket_id = ghost_comm_map[proc][idata].bucket_id;
//          const unsigned ord    = ghost_comm_map[proc][idata].bucket_ord;
//          EXPECT_EQ(aura_ghosting.ordinal(), ghost_id);
//          const stk::mesh::Bucket& bucket = *all_node_buckets[bucket_id];
//          stk::mesh::Entity node = bucket[ord];
//          double* data = stk::mesh::field_data(field, node);
//          if (owner == myProc) {
//              EXPECT_EQ(*data, *reinterpret_cast<double*>(&send_data[proc][send_offsets[proc]]));
//              send_offsets[proc] += sizeof(double);
//          }
//          else {
//              EXPECT_EQ(*data, *reinterpret_cast<double*>(&recv_data[proc][recv_offsets[proc]]));
//              recv_offsets[proc] += sizeof(double);
//          }
//      }
//  }
//}

TEST(CommunicateFieldData, communicate)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }
  int myProc = stk::parallel_machine_rank(communicator);

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::Field<double>& field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field(field, meta.universal_part());

  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case1(bulk);

  stk::mesh::Part& owned_part = meta.locally_owned_part();
  const stk::mesh::BucketVector& owned_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK,owned_part);

  for(size_t i=0; i<owned_node_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *owned_node_buckets[i];
    for(size_t j=0; j<bucket.size(); ++j) {
        stk::mesh::Entity node = bucket[j];
        double value = myProc*100 + bulk.identifier(node);
        double* data = stk::mesh::field_data(field, node);
        *data = value;
    }
  }

  const stk::mesh::Ghosting& aura_ghosting = *bulk.ghostings()[stk::mesh::BulkData::AURA];

  std::vector<const stk::mesh::FieldBase*> fields(1, &field);
  stk::mesh::communicate_field_data(aura_ghosting, fields);

  stk::mesh::Part& aura_part = meta.aura_part();
  const stk::mesh::BucketVector& aura_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, aura_part);
  for(size_t i=0; i<aura_node_buckets.size(); ++i) {
      stk::mesh::Bucket& bucket = *aura_node_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
          stk::mesh::Entity node = bucket[j];
          int owner = bulk.parallel_owner_rank(node);
          double value = owner*100 + bulk.identifier(node);
          double* data = stk::mesh::field_data(field, node);
          EXPECT_EQ(*data, value);
      }
  }
}



TEST(CommunicateFieldData, communicateMultipleGhostings)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }
  int myProc = stk::parallel_machine_rank(communicator);

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::Field<double>& field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field(field, meta.universal_part());

  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case2(bulk);

  bulk.modification_begin();
  // Create ghosting
  stk::mesh::Ghosting & ghosting_1 = bulk.create_ghosting( "CUSTOM_1" );
  stk::mesh::Ghosting & ghosting_2 = bulk.create_ghosting( "CUSTOM_2" );
  stk::mesh::Ghosting & ghosting_3 = bulk.create_ghosting( "CUSTOM_3" );

  std::vector<stk::mesh::EntityProc> ghosting_1_send;
  std::vector<stk::mesh::EntityProc> ghosting_2_send;
  std::vector<stk::mesh::EntityProc> ghosting_3_send;
  if (myProc == 1)
  {
    const int proc_0 = 0;
    ghosting_1_send.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(7)),proc_0));
    ghosting_2_send.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(8)),proc_0));
    ghosting_3_send.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(9)),proc_0));
    ghosting_3_send.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, stk::mesh::EntityId(7)),proc_0)); // duplicate entity
  }
  bulk.change_ghosting( ghosting_1, ghosting_1_send );
  bulk.change_ghosting( ghosting_2, ghosting_2_send );
  bulk.change_ghosting( ghosting_3, ghosting_3_send );
  bulk.modification_end();

  stk::mesh::Part& owned_part = meta.locally_owned_part();
  const stk::mesh::BucketVector& owned_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK,owned_part);

  for(size_t i=0; i<owned_node_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *owned_node_buckets[i];
    for(size_t j=0; j<bucket.size(); ++j) {
        stk::mesh::Entity node = bucket[j];
        double value = myProc*100 + bulk.identifier(node);
        double* data = stk::mesh::field_data(field, node);
        *data = value;
    }
  }

  stk::mesh::Selector select_ghosted = !meta.locally_owned_part() & !meta.globally_shared_part();
  const stk::mesh::BucketVector& ghosted_buckets = bulk.get_buckets(stk::topology::NODE_RANK,select_ghosted);

  int num_ghosted_nodes = 0;
  for(size_t i=0; i<ghosted_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *ghosted_buckets[i];
    for(size_t j=0; j<bucket.size(); ++j) {
        ++num_ghosted_nodes;
        stk::mesh::Entity node = bucket[j];
        double garbage = -1.2345;
        double* data = stk::mesh::field_data(field, node);
        *data = garbage;
    }
  }
  if (myProc == 0) {
      ASSERT_EQ(5, num_ghosted_nodes);
  }
  else { // myProc == 1
      ASSERT_EQ(2, num_ghosted_nodes);
  }

  std::vector<const stk::mesh::FieldBase*> fields(1, &field);
  stk::mesh::communicate_field_data(bulk, fields);

  for(size_t i=0; i<ghosted_buckets.size(); ++i) {
      stk::mesh::Bucket& bucket = *ghosted_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
          stk::mesh::Entity node = bucket[j];
          int owner = bulk.parallel_owner_rank(node);
          double value = owner*100 + bulk.identifier(node);
          double* data = stk::mesh::field_data(field, node);
          EXPECT_EQ(value, *data);
      }
  }
}
