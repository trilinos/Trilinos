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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <vector>                       // for vector, operator==
#include "mpi.h"                        // for MPI_COMM_WORLD, etc

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }

using stk::unit_test_util::build_mesh;

TEST( UnitTestChangeEntityId, change_id_small )
{
  // change_entity_id is a broken concept in parallel:
  //   Only entity owners should have to call change-id. As it stands now,
  //   every process that might know about the entity, even as a ghost,
  //   must make the call. If any sharer/ghoster of the changing entity
  //   forgets to make the call, or changes to wrong id, awful errors ensue
  //   in the next modification_end. The only way to make this safe is to have
  //   change_entity_id be a parallel call similar to change_entity_owner; this
  //   would be very expensive. We need to implement the tracking of exodus-ids
  //   in a different fashion in the long run.

  // Demonstrate how change_entity_id should work

  using namespace stk::mesh;

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned spatial_dim = 2;
  std::shared_ptr<BulkData> bulkPtr = build_mesh(spatial_dim, pm);
  MetaData& meta_data = bulkPtr->mesh_meta_data();
  Part &quad4_part = meta_data.declare_part_with_topology("quad4_Part", stk::topology::QUAD_4_2D);
  meta_data.commit();

  BulkData& mesh = *bulkPtr;
  int p_rank = mesh.parallel_rank();
  int p_numProcs = mesh.parallel_size();
  if (p_numProcs > 1) {
    //change-entity-id test only supported in serial
    return;
  }

  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&quad4_part);

  mesh.modification_begin();

  Entity elem                 = mesh.declare_element(p_rank + 1 /*id*/, elem_parts);

  Entity node0                = mesh.declare_node(p_rank*3 + 5000 /*id*/);
  Entity node1_local_chg_id   = mesh.declare_node(p_rank*3 + 1 /*id*/);
  Entity node2_shared_chg_id  = mesh.declare_node(           2 /*id*/);
  Entity node3                = mesh.declare_node(p_rank*3 + 3 /*id*/);

  mesh.declare_relation(elem, node0              , 0 /*relation ordinal*/);
  mesh.declare_relation(elem, node1_local_chg_id , 1 /*relation ordinal*/);
  mesh.declare_relation(elem, node2_shared_chg_id, 2 /*relation ordinal*/);
  mesh.declare_relation(elem, node3              , 3 /*relation ordinal*/);

  mesh.modification_end();

  mesh.modification_begin();

  const EntityId new_id_for_local_node = 42 + p_rank + 1;
  const EntityId new_id_for_shared_node = 666;

  mesh.change_entity_id(new_id_for_local_node, node1_local_chg_id);
  if (mesh.parallel_owner_rank(node2_shared_chg_id) == p_rank) {
    mesh.change_entity_id(new_id_for_shared_node, node2_shared_chg_id);
  }

  EXPECT_EQ(new_id_for_local_node, mesh.identifier(node1_local_chg_id));

  mesh.modification_end();

  EXPECT_EQ(new_id_for_local_node, mesh.identifier(node1_local_chg_id));
  EXPECT_EQ(new_id_for_shared_node, mesh.identifier(node2_shared_chg_id));
}

TEST( UnitTestChangeEntityId, change_id_large )
{
  using namespace stk::mesh;

  const unsigned NX = 20;
  const unsigned NY = 20;
  const unsigned NZ = 20;
  const unsigned num_elems = NX * NY * NZ;

  fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

  Field<int> & simple_nodal_field = hf.m_meta.declare_field<int>(stk::topology::NODE_RANK, "simple_nodal_field");

  put_field_on_mesh(simple_nodal_field, *hf.m_elem_parts[0], nullptr);

  //create nodal field on hex topo

  hf.m_meta.commit();

  hf.generate_mesh();

  stk::mesh::BulkData & mesh = hf.m_bulk_data;
  if (mesh.parallel_size() > 1) {
    return;
  }

  mesh.modification_begin();

  const BucketVector & node_buckets = mesh.buckets(stk::topology::NODE_RANK);

  for(Bucket * b : node_buckets) {
    int* nodal_field = stk::mesh::field_data( simple_nodal_field, *b );
    for (size_t i =0; i<b->size(); ++i) {
      nodal_field[i] = 1;
    }
  }

  const BucketVector & elem_buckets = mesh.buckets(stk::topology::ELEMENT_RANK);

  std::vector<EntityId> old_ids;
  old_ids.reserve(num_elems);
  for(Bucket * b : elem_buckets) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity e = (*b)[i];
      old_ids.push_back(mesh.identifier(e));
      mesh.change_entity_id( mesh.identifier(e)+num_elems, e);
    }
  }

  mesh.modification_end();

  mesh.modification_begin();
  mesh.modification_end();

  std::vector<EntityId> new_ids_minus_num_elems;
  new_ids_minus_num_elems.reserve(num_elems);
  for(Bucket * b : elem_buckets) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity e = (*b)[i];
      new_ids_minus_num_elems.push_back(mesh.identifier(e)-num_elems);
    }
  }

  EXPECT_TRUE(old_ids == new_ids_minus_num_elems);

  for(Bucket * b : node_buckets) {
    int* nodal_field = stk::mesh::field_data( simple_nodal_field, *b );
    for (size_t i =0; i<b->size(); ++i) {
      EXPECT_TRUE( nodal_field[i] == 1);
    }
  }
}
