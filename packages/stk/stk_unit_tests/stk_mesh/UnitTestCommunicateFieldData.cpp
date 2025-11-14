// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "SetupKeyholeMesh.hpp"
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator!, Selector, etc
#include "stk_mesh/base/Types.hpp"      // for EntityProc, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <iostream>                     // for basic_ostream::operator<<
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <vector>                       // for vector
namespace stk { namespace mesh { class Ghosting; } }
using stk::unit_test_util::build_mesh;

TEST(CommunicateFieldData, communicate)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }
  int myProc = stk::parallel_machine_rank(communicator);

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::mesh::Field<double>& field = meta.declare_field<double>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);

  stk::mesh::unit_test::setupKeyholeMesh2D_case1(bulk);

  stk::mesh::Part& owned_part = meta.locally_owned_part();
  const stk::mesh::BucketVector& owned_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK,owned_part);

  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for(size_t i=0; i<owned_node_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *owned_node_buckets[i];
    auto bucketFieldData = fieldData.bucket_values(bucket);
    for(stk::mesh::EntityIdx j : bucket.entities()) {
      stk::mesh::Entity node = bucket[j];
      double value = myProc*100 + bulk.identifier(node);
      bucketFieldData(j,0_comp) = value;
    }
  }

  const stk::mesh::Ghosting& aura_ghosting = *bulk.ghostings()[stk::mesh::BulkData::AURA];

  std::vector<const stk::mesh::FieldBase*> fields(1, &field);
  stk::mesh::communicate_field_data(aura_ghosting, fields);

  stk::mesh::Part& aura_part = meta.aura_part();
  const stk::mesh::BucketVector& aura_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, aura_part);
  fieldData = field.data<stk::mesh::ReadWrite>();
  for(size_t i=0; i<aura_node_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *aura_node_buckets[i];
    auto bucketFieldData = fieldData.bucket_values(bucket);
    for(stk::mesh::EntityIdx j : bucket.entities()) {
      stk::mesh::Entity node = bucket[j];
      int owner = bulk.parallel_owner_rank(node);
      double value = owner*100 + bulk.identifier(node);
      EXPECT_EQ(bucketFieldData(j,0_comp), value);
    }
  }
}



TEST(CommunicateFieldData, communicateMultipleGhostings)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  const int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) { GTEST_SKIP(); }
  const int myProc = stk::parallel_machine_rank(communicator);

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::mesh::Field<double>& field = meta.declare_field<double>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);

  stk::mesh::unit_test::setupKeyholeMesh2D_case2(bulk);

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

  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for(size_t i=0; i<owned_node_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *owned_node_buckets[i];
    auto bucketFieldData = fieldData.bucket_values(bucket);
    for(stk::mesh::EntityIdx j : bucket.entities()) {
      stk::mesh::Entity node = bucket[j];
      double value = myProc*100 + bulk.identifier(node);
      bucketFieldData(j,0_comp) = value;
    }
  }

  stk::mesh::Selector select_ghosted = (!meta.locally_owned_part()) & (!meta.globally_shared_part());
  const stk::mesh::BucketVector& ghosted_buckets = bulk.get_buckets(stk::topology::NODE_RANK,select_ghosted);

  int num_ghosted_nodes = 0;
  for(size_t i=0; i<ghosted_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *ghosted_buckets[i];
    auto bucketFieldData = fieldData.bucket_values(bucket);
    for(stk::mesh::EntityIdx j : bucket.entities()) {
      ++num_ghosted_nodes;
      double garbage = -1.2345;
      bucketFieldData(j,0_comp) = garbage;
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

  fieldData = field.data<stk::mesh::ReadWrite>();
  for(size_t i=0; i<ghosted_buckets.size(); ++i) {
    stk::mesh::Bucket& bucket = *ghosted_buckets[i];
    auto bucketFieldData = fieldData.bucket_values(bucket);
    for(stk::mesh::EntityIdx j : bucket.entities()) {
      stk::mesh::Entity node = bucket[j];
      int owner = bulk.parallel_owner_rank(node);
      double value = owner*100 + bulk.identifier(node);
      EXPECT_EQ(value, bucketFieldData(j,0_comp));
    }
  }
}
