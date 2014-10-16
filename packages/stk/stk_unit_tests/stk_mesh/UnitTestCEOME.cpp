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
#include <algorithm>                    // for sort
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>  // for test_shift_ring
#include <unit_tests/Setup8Quad4ProcMesh.hpp>
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
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::BaseEntityRank;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace
{

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }
    virtual ~BulkDataTester()
    {
    }

    void my_update_comm_list_based_on_changes_in_comm_map()
    {
        this->update_comm_list_based_on_changes_in_comm_map();
    }
};

const EntityRank NODE_RANK = stk::topology::NODE_RANK;
const EntityRank EDGE_RANK = stk::topology::EDGE_RANK;
const EntityRank FACE_RANK = stk::topology::FACE_RANK;
const EntityRank ELEM_RANK = stk::topology::ELEM_RANK;

enum EntityStates {
  STATE_VALID,
  STATE_NOT_VALID,
  STATE_OWNED,
  STATE_SHARED,
  STATE_NOT_SHARED,
  STATE_GHOSTED,
  STATE_NOT_GHOSTED
};

bool check_state_CEOME(const stk::mesh::BulkData & mesh, const EntityKey & entityKey, EntityStates state,
                 int p0 = -1, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1, int p5 = -1);

void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem,
                       int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move);
//==============================================================================

void addSharingInfo(stk::mesh::BulkData& bulkData, stk::mesh::Entity entity, stk::mesh::BulkData::GHOSTING_ID ghostingId, int sharingProc )
{
    EXPECT_TRUE(bulkData.entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(ghostingId, sharingProc)));
}

void eraseSharingInfo(stk::mesh::BulkData &bulkData, stk::mesh::Entity entity, stk::mesh::BulkData::GHOSTING_ID ghostingId )
{
    stk::mesh::EntityKey key = bulkData.entity_key(entity);
    bulkData.entity_comm_map_erase(key, *bulkData.ghostings()[ghostingId]);
}

void eraseSharingInfoUsingKey(stk::mesh::BulkData &bulkData, stk::mesh::EntityKey key, stk::mesh::BulkData::GHOSTING_ID ghostingId )
{
    bulkData.entity_comm_map_erase(key, *bulkData.ghostings()[ghostingId]);
}

TEST(CEOME, change_entity_owner_2Elem2ProcMove)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part, element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part, element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  // Check initial state
  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID) );
  }

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(NODE_RANK, 6), 1));
  }

  bulk.change_entity_owner(entity_procs);

  ////////////////////////////////////////////////////////////////////////////

  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  std::vector<bool> nodeSharing;
  if ( p_rank == 0 )
  {
      bool isNodeShared[] = { false, false, true, true, false, false };
      nodeSharing.assign(isNodeShared, isNodeShared+6);
  }
  else
  {
      bool isNodeShared[] = { false, false, true, true, false, false };
      nodeSharing.assign(isNodeShared, isNodeShared+6);
  }

  int otherProc = 1;
  if ( p_rank == 1 )
  {
      otherProc = 0;
  }

//  bulk.modification_begin();

  for (size_t i=0;i<nodeSharing.size();i++)
  {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, i+1);
      eraseSharingInfo(bulk, node, stk::mesh::BulkData::SHARED);
      if ( nodeSharing[i] )
      {
          addSharingInfo(bulk, node, stk::mesh::BulkData::SHARED, otherProc);
      }
  }

  bulk.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  bulk.modification_begin();
  bulk.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  ////////////////////////////////////////////////////////////////////////////

  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
  }

}

TEST(CEOME, change_entity_owner_2Elem2ProcFlip)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData mesh( meta, pm);

  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  mesh.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(mesh, elem_part, element_ids[0], elem_node_ids[0] ) );
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3)), 1);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 1);
  }
  else if (p_rank == 1) {
    elems.push_back(stk::mesh::declare_element(mesh, elem_part, element_ids[1], elem_node_ids[1] ) );
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3)), 0);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 0);
  }

  mesh.modification_end();

  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
  }

  //okay now flip
  //    1/0---4/0---5/1          1/1---4/0---5/0
  //     |     |     |            |     |     |
  //     | 1/0 | 2/1 |       =>   | 1/1 | 2/0 |
  //     |     |     |            |     |     |
  //    2/0---3/0---6/1          2/1---3/0---6/0

  stk::mesh::EntityProcVec entity_procs_flip;
  if (p_rank == 0) {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 2), 1));
  } else {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 2), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 5), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 6), 0));
  }
  mesh.change_entity_owner(entity_procs_flip);

  ////////////////////////////////////////////////////////////////////////////


  std::vector<bool> nodeSharing;
  if ( p_rank == 0 )
  {
      bool isNodeShared[] = { false, false, true, true, false, false };
      nodeSharing.assign(isNodeShared, isNodeShared+6);
  }
  else
  {
      bool isNodeShared[] = { false, false, true, true, false, false };
      nodeSharing.assign(isNodeShared, isNodeShared+6);
  }

  int otherProc = 1;
  if ( p_rank == 1 )
  {
      otherProc = 0;
  }

  for (size_t i=0;i<nodeSharing.size();i++)
  {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
      eraseSharingInfo(mesh, node, stk::mesh::BulkData::SHARED);
      if ( nodeSharing[i] )
      {
          addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
      }
  }

  mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  mesh.modification_begin();
  mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  ////////////////////////////////////////////////////////////////////////////

  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  0) );
  }

}

TEST(CEOME, change_entity_owner_3Elem2ProcMoveRight)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 2)
  {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  EntityVector nodes;
  EntityVector elements;
  if (p_rank == 0) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 1, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 2, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 3, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 4, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 5, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 6, node_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 1, elem_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 2, elem_part));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.declare_relation(elements[1], nodes[2], 0);
    mesh.declare_relation(elements[1], nodes[3], 1);
    mesh.declare_relation(elements[1], nodes[4], 2);
    mesh.declare_relation(elements[1], nodes[5], 3);

    mesh.add_node_sharing(nodes[4], 1);
    mesh.add_node_sharing(nodes[5], 1);
  }
  else if (p_rank == 1) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 5, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 6, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 7, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 8, node_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 3, elem_part));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.add_node_sharing(nodes[0], 0);
    mesh.add_node_sharing(nodes[1], 0);
  }

  mesh.modification_end();


  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED) );
  }


  std::vector<EntityProc> change;
  if(p_rank == 0)
  {
    change.push_back(EntityProc(elements[1], 1));
    change.push_back(EntityProc(nodes[4], 1));
    change.push_back(EntityProc(nodes[5], 1));
  }

  mesh.change_entity_owner(change);

  ////////////////////////////////////////////////////////////////////////////

  //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    size_t numNodes = 8;
    std::vector<bool> nodeSharing;

    if ( p_rank == 0 )
    {
        bool isNodeShared[] = { false, false, true, true, false, false, false, false };
        nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
    }
    else
    {
        bool isNodeShared[] = { false, false, true, true, false, false, false, false };
        nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    int otherProc = 1;
    if ( p_rank == 1 )
    {
        otherProc = 0;
    }

    for (size_t i=0;i<nodeSharing.size();i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
        if ( mesh.is_valid(node) )
        {
            eraseSharingInfo(mesh, node, stk::mesh::BulkData::SHARED);
        }
        if ( nodeSharing[i] )
        {
            addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
        }
    }

    // mesh.update_comm_list_based_on_changes_in_comm_map();

    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED) );
  }
}

TEST(CEOME, change_entity_owner_3Elem2ProcMoveLeft)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 2)
  {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  EntityVector nodes;
  EntityVector elements;
  if (p_rank == 0) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 1, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 2, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 3, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 4, node_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 1, elem_part));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.add_node_sharing(nodes[2], 1);
    mesh.add_node_sharing(nodes[3], 1);
  }
  else if (p_rank == 1) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 3, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 4, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 5, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 6, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 7, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 8, node_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 2, elem_part));
    elements.push_back(mesh.declare_entity(ELEM_RANK, 3, elem_part));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.declare_relation(elements[1], nodes[2], 0);
    mesh.declare_relation(elements[1], nodes[3], 1);
    mesh.declare_relation(elements[1], nodes[4], 2);
    mesh.declare_relation(elements[1], nodes[5], 3);

    mesh.add_node_sharing(nodes[0], 0);
    mesh.add_node_sharing(nodes[1], 0);
  }

  mesh.modification_end();


  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED) );
  }

  std::vector<EntityProc> change;
  if(p_rank == 1)
  {
    change.push_back(EntityProc(elements[0], 0));
    change.push_back(EntityProc(nodes[2], 0));
    change.push_back(EntityProc(nodes[3], 0));
  }

  mesh.change_entity_owner(change);

  ////////////////////////////////////////////////////////////////////////////

  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1


     size_t numNodes = 8;
     std::vector<bool> nodeSharing;

     if ( p_rank == 0 )
     {
         bool isNodeShared[] = { false, false, false, false, true, true, false, false };
         nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
     }
     else
     {
         bool isNodeShared[] = { false, false, false, false, true, true, false, false };
         nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
     }

     EXPECT_EQ(numNodes, nodeSharing.size());

     int otherProc = 1;
     if ( p_rank == 1 )
     {
         otherProc = 0;
     }

     for (size_t i=0;i<nodeSharing.size();i++)
     {
         stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
         if ( mesh.is_valid(node) )
         {
             eraseSharingInfo(mesh, node, stk::mesh::BulkData::SHARED);
         }
         if ( nodeSharing[i] )
         {
             addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
         }
     }

     // mesh.update_comm_list_based_on_changes_in_comm_map();

     mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

     mesh.modification_begin();
     mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

     ////////////////////////////////////////////////////////////////////////////

  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED) );
  }

}

TEST(CEOME, change_entity_owner_4Elem4ProcEdge)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    {|     |          |      |     |    {|     |
  //        | 1/0  | 2/1 | 3/2{| 4/3 |          | 1/0  | 2/1 | 3/0{| 4/3 |
  //        |      |     |    {|     |          |      |     |    {|     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
  //  this edge moves to p0 --^
  //  element 3 moves to proc 0.
  //  nodes 7&8 move to proc 0.
  //  proc 2 forgets everything.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  Part& edge_part = meta_data.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();
  BulkDataTester mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 4)
  {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  EntityKey elem_key_chg_own(ELEM_RANK, 3 /*id*/);
  EntityKey edge_key_chg_own(EDGE_RANK, 1 /*id*/);
  EntityKey node_A_key_chg_own(NODE_RANK, 7 /*id*/);
  EntityKey node_B_key_chg_own(NODE_RANK, 8 /*id*/);
  EntityKey node_C_key(NODE_RANK, 5 /*id*/);
  EntityKey node_D_key(NODE_RANK, 6 /*id*/);

  // Create element
  Entity elem = mesh.declare_entity(ELEM_RANK, p_rank + 1, //elem_id
                                    elem_part);

  // If it is 2nd to last element, it is the one changing
  if(p_rank == 2)
  {
    EXPECT_TRUE(elem_key_chg_own == mesh.entity_key(elem));
  }

  // Create nodes
  EntityVector nodes;
  if(p_rank == 0)
  {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 1, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 2, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 3, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 4, node_part));
  }
  else if(p_rank == 1)
  {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 3, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 4, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 5, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 6, node_part));
  }
  else if(p_rank == 2)
  {
    nodes.push_back(mesh.declare_entity(NODE_RANK, 5, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 6, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 7, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 8, node_part));
  }
  else
  { // p_rank == 3
    nodes.push_back(mesh.declare_entity(NODE_RANK, 7, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 8, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 9, node_part));
    nodes.push_back(mesh.declare_entity(NODE_RANK, 10, node_part));
  }

  // Add element relations to nodes
  unsigned rel_id = 0;
  for(EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
  {
    mesh.declare_relation(elem, *itr, rel_id);
  }

  // Create edge on last two procs

  if(p_rank >= 2)
  {
    Entity edge = mesh.declare_entity(EDGE_RANK, 1, // id
                                      edge_part);
    EXPECT_TRUE(mesh.entity_key(edge) == edge_key_chg_own);

    // Add element relation to edge
    mesh.declare_relation(elem, edge, 1 /*rel-id*/);

    // Add edge relations to nodes
    if(p_rank == 2)
    {
      mesh.declare_relation(edge, nodes[2], 0);
      mesh.declare_relation(edge, nodes[3], 1);
    }
    else
    { // p_rank == 3
      mesh.declare_relation(edge, nodes[0], 0);
      mesh.declare_relation(edge, nodes[1], 1);
    }
  }
  if(p_rank == 0)
  {
    mesh.add_node_sharing(nodes[2], 1);
    mesh.add_node_sharing(nodes[3], 1);
  }
  else if((p_rank == 1) || (p_rank == 2))
  {
    mesh.add_node_sharing(nodes[0], p_rank - 1);
    mesh.add_node_sharing(nodes[1], p_rank - 1);
    mesh.add_node_sharing(nodes[2], p_rank + 1);
    mesh.add_node_sharing(nodes[3], p_rank + 1);
  }
  else
  { // p_rank ==3
    mesh.add_node_sharing(nodes[0], p_rank - 1);
    mesh.add_node_sharing(nodes[1], p_rank - 1);
  }

  mesh.modification_end();

  //test pre-conditions
  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED, 2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED,  2) );
  }
  else if (p_rank == 2) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED,  3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED,  3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED,    2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED,   3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_SHARED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_SHARED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED,  3) );
  }
  else if (p_rank == 3) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED,  3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED,    2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED,   2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_SHARED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_SHARED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED) );
  }

  std::vector<EntityProc> change;
  if(p_rank == 2)
  {
    // Change ownership of changing elem and all entities in it's closure that
    // we own to proc 0.

    Entity changing_elem = mesh.get_entity(elem_key_chg_own);
    ASSERT_TRUE( mesh.is_valid(changing_elem));
    EntityProc eproc(changing_elem, 0 /*new owner*/);
    change.push_back(eproc);

    const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
    {
      stk::mesh::Entity const *to_i = mesh.begin(changing_elem, irank);
      stk::mesh::Entity const *to_e = mesh.end(changing_elem, irank);
      for (; to_i != to_e; ++to_i)
      {
        if (mesh.parallel_owner_rank(*to_i) == p_rank)
        {
          EntityProc eproc_new(*to_i, 0 /*new owner*/);
          change.push_back(eproc_new);
        }
      }
    }
  }

  mesh.change_entity_owner(change);

  ////////////////////////////////////////////////////////////////////////////

  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    {|     |          |      |     |    {|     |
  //        | 1/0  | 2/1 | 3/2{| 4/3 |          | 1/0  | 2/1 | 3/0{| 4/3 |
  //        |      |     |    {|     |          |      |     |    {|     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3


       size_t numNodes = 10;
       std::vector<bool> nodeSharing;
       std::vector<int>  nodeSharingProcs;

       size_t numEdges = 1;
       std::vector<bool> edgeSharing;
       std::vector<int>  edgeSharingProcs;

       std::vector<stk::mesh::Entity> modifiedEntities;

       if ( p_rank == 0 )
       {
           bool isNodeShared[] = { false, false, true, true, true, true, true, true, false, false };
           nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
           int  procIdsNode[]  = {    -1,    -1,    1,    1,    1,    1,    3,    3,    -1,    -1 };
           nodeSharingProcs.assign(procIdsNode, procIdsNode+numNodes);
           bool isEdgeShared[] = { true };
           int  procIdsEdge[] = { 3 };
           edgeSharing.assign(isEdgeShared, isEdgeShared+numEdges);
           edgeSharingProcs.assign(procIdsEdge, procIdsEdge+numEdges);
       }
       else if ( p_rank == 1)
       {
           bool isNodeShared[] = { false, false, true, true, true, true, false, false, false, false };
           nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
           int  procIdsNode[]  = {    -1,    -1,    0,    0,    0,    0,    -1,    -1,    -1,    -1 };
           nodeSharingProcs.assign(procIdsNode, procIdsNode+numNodes);
           bool isEdgeShared[] = { false };
           int  procIdsEdge[] = { -1 };
           edgeSharing.assign(isEdgeShared, isEdgeShared+numEdges);
           edgeSharingProcs.assign(procIdsEdge, procIdsEdge+numEdges);
       }
       else if ( p_rank == 2)
       {
           bool isNodeShared[] = { false, false, false, false, false, false, false, false, false, false };
           nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
           int  procIdsNode[]  = {    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1 };
           nodeSharingProcs.assign(procIdsNode, procIdsNode+numNodes);
           bool isEdgeShared[] = { false };
           int  procIdsEdge[] = { -1 };
           edgeSharing.assign(isEdgeShared, isEdgeShared+numEdges);
           edgeSharingProcs.assign(procIdsEdge, procIdsEdge+numEdges);
       }
       else
       {
           bool isNodeShared[] = { false, false, false, false, false, false, true, true, false, false };
           nodeSharing.assign(isNodeShared, isNodeShared+numNodes);
           int  procIdsNode[]  = {    -1,    -1,    -1,    -1,    -1,    -1,    0,    0,    -1,    -1 };
           nodeSharingProcs.assign(procIdsNode, procIdsNode+numNodes);
           bool isEdgeShared[] = { true };
           int  procIdsEdge[] = { 0 };
           edgeSharing.assign(isEdgeShared, isEdgeShared+numEdges);
           edgeSharingProcs.assign(procIdsEdge, procIdsEdge+numEdges);
       }

       EXPECT_EQ(numNodes, nodeSharing.size());

       for (size_t i=0;i<nodeSharing.size();i++)
       {
           stk::mesh::EntityKey key(stk::topology::NODE_RANK, i+1);
           eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

           if ( nodeSharing[i] )
           {
               stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
               addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeSharingProcs[i]);
               modifiedEntities.push_back(node);
           }
       }

       for (size_t i=0;i<edgeSharing.size();i++)
       {
           stk::mesh::EntityKey key(stk::topology::EDGE_RANK, i+1);
           eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

           if ( edgeSharing[i] )
           {
               stk::mesh::Entity edge = mesh.get_entity(stk::topology::EDGE_RANK, i+1);
               addSharingInfo(mesh, edge, stk::mesh::BulkData::SHARED, edgeSharingProcs[i]);
               modifiedEntities.push_back(edge);
           }
       }

       mesh.my_update_comm_list_based_on_changes_in_comm_map();
       mesh.update_comm_list(modifiedEntities);

       mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

       mesh.modification_begin();
       mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

       ////////////////////////////////////////////////////////////////////////////

  //test post condition
  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED,  3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED,    0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED,   3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_OWNED, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_SHARED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_SHARED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_SHARED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED,  3) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED, 0) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED,  0) );
  }
  else if (p_rank == 2) { //amnesia
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID) );
  }
  else if (p_rank == 3) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED,  3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED,    0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED,   0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_SHARED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_SHARED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5),  STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6),  STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9),  STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED) );
  }

}

TEST(CEOME, change_entity_owner_8Elem4ProcMoveTop)
{
    //
    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
    //
    // This test moves ownership of elements 6 and 7 (as well as their locally-owned
    // nodes) to procs 3 and 0, respectively.
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if(numProcs != 4)
    {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    BulkDataTester mesh(meta, pm);
    const int p_rank = mesh.parallel_rank();

    setup8Quad4ProcMesh2D(mesh);

    // Check the initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED, 1));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED, 1));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED, 2));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED, 3));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED));
    }

    std::vector<stk::mesh::EntityProc> entities_to_move;
    if(mesh.parallel_rank() == 1)
    {
        stk::mesh::Entity elem = mesh.get_entity(ELEM_RANK, 6);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    if(mesh.parallel_rank() == 2)
    {
        stk::mesh::Entity elem = mesh.get_entity(ELEM_RANK, 7);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }

    mesh.change_entity_owner(entities_to_move);

    ////////////////////////////////////////////////////////////////////////////

    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3

    size_t numNodes = 15;
    std::vector<bool> nodeSharing;
    std::vector<std::vector<int> > nodeSharingProcs(15);

    std::vector<stk::mesh::Entity> modifiedEntities;

    if(p_rank == 0)
    {
        bool isNodeShared[] = {
                false, true, false, false, false,
                false, true, true, true, false,
                false, true, true, true, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // node 11: nothing
        // nodes 12, 13, 14
        {
            int procsForNodeSharing[] = {3};
            nodeSharingProcs[12 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[13 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[14 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 15
        // node 6
        // node 7
        {
            int procsForNodeSharing[] = {1, 3};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {1, 2, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {2, 3};
            nodeSharingProcs[9 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 10
        // node 1
        // node 2
        {
            int procsForNodeSharing[] = {1};
            nodeSharingProcs[2 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // nodes 3, 4, 5
    }
    else if(p_rank == 1)
    {
        bool isNodeShared[] = {
                false, true, true, false, false,
                false, true, true, false, false,
                false, false, false, false, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 11-15: nothing
        // nodes 6, 9, 10: nothing
        // node 7
        {
            int procsForNodeSharing[] = {0, 3};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {0, 2, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // nodes 1, 4, 5: nothing
        // node 2
        {
            int procsForNodeSharing[] = {0};
            nodeSharingProcs[2 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 3
        {
            int procsForNodeSharing[] = {2};
            nodeSharingProcs[3 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }
    else if(p_rank == 2)
    {
        bool isNodeShared[] = {
                false, false, true, true, false,
                false, false, true, true, false,
                false, false, false, false, false
        };

        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 11-15: nothing
        // nodes 6, 7, 10: nothing
        // node 8
        {
            int procsForNodeSharing[] = {0, 1, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {0, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // nodes 1, 2, 5: nothing
        // node 3
        {
            int procsForNodeSharing[] = {1};
            nodeSharingProcs[3 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 4
        {
            int procsForNodeSharing[] = {3};
            nodeSharingProcs[4 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }
    else
    {
        bool isNodeShared[] = {
                false, false, false, true, false,
                false, true, true, true, false,
                false, true, true, true, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 1, 2, 3, 5, 6, 10, 11, 15: nothing
        // nodes 12, 13, 14
        {
            int procsForNodeSharing[] = {0};
            nodeSharingProcs[12 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[13 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[14 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 7
        {
            int procsForNodeSharing[] = {0, 1};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {0, 1, 2};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {0, 2};
            nodeSharingProcs[9 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 4
        {
            int procsForNodeSharing[] = {2};
            nodeSharingProcs[4 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);
            for(size_t j = 0; j < nodeSharingProcs[i].size(); j++)
            {
                addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeSharingProcs[i][j]);
            }
            modifiedEntities.push_back(node);
        }
    }

    mesh.my_update_comm_list_based_on_changes_in_comm_map();
    mesh.update_comm_list(modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    // Check the final state
    if(p_rank == 0)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 1, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 2, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED, 3));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED, 0));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED, 0));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 0, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED, 3));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 0, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED));
    }
}

TEST(CEOME, change_entity_owner_4Elem4ProcRotate)
{
    //
    //     id/proc                id/proc
    //      7/3---8/2---9/2        7/2---8/1---9/1
    //       |     |     |          |     |     |
    //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
    //       |     |     |          |     |     |
    //      4/0---5/0---6/1  -->   4/3---5/3---6/0
    //       |     |     |          |     |     |
    //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
    //       |     |     |          |     |     |
    //      1/0---2/0---3/1        1/3---2/3---3/0

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if(numProcs != 4)
    {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    BulkDataTester mesh(meta, pm);
    const int p_rank = mesh.parallel_rank();

    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
    meta.commit();

    // 1 elem-id for each proc
    stk::mesh::EntityId proc_elemIDs[] = {1, 2, 3, 4};

    // list of node-ids for each element
    const int nodesPerElem = 4;
    stk::mesh::EntityId elem_nodeIDs[][nodesPerElem] = {
            {1, 2, 5, 4},
            {2, 3, 6, 5},
            {5, 6, 9, 8},
            {4, 5, 8, 7}
    };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    int shared_nodeIDs_and_procs[][3] = {
            {0, 2, 1}, {0, 5, 1}, {0, 5, 2}, {0, 5, 3}, {0, 4, 3}, // proc 0
            {1, 2, 0}, {1, 6, 2}, {1, 5, 0}, {1, 5, 2}, {1, 5, 3}, // proc 1
            {2, 5, 0}, {2, 5, 1}, {2, 5, 3}, {2, 6, 1}, {2, 8, 3}, // proc 2
            {3, 4, 0}, {3, 5, 0}, {3, 5, 1}, {3, 5, 2}, {3, 8, 2} // proc 3
    };
    int numSharedNodeTriples = 20;

    mesh.modification_begin();

    stk::mesh::EntityId elemID = proc_elemIDs[p_rank];
    stk::mesh::declare_element(mesh, block_1, elemID, elem_nodeIDs[p_rank]);

    for(int proc = 0; proc < numSharedNodeTriples; ++proc)
    {
        if(p_rank == shared_nodeIDs_and_procs[proc][0])
        {
            int nodeID = shared_nodeIDs_and_procs[proc][1];
            int sharingProc = shared_nodeIDs_and_procs[proc][2];
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
            mesh.add_node_sharing(node, sharingProc);
        }
    }

    mesh.modification_end();

    // Check the initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 2));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_GHOSTED));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 3));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_NOT_GHOSTED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 2));
    }

    std::vector<stk::mesh::EntityProc> entities_to_move;
    if(p_rank == 0)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 1)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 2)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
        int dest_proc = 1;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 3)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 4);
        int dest_proc = 2;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }

    mesh.change_entity_owner(entities_to_move);

    std::vector<std::pair<int, int> > entities;

    int numNodes = 9;
    if(p_rank == 0)
    {
        entities.push_back(std::make_pair(2, 3));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 2));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(6, 1));

    }
    else if(p_rank == 1)
    {
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 2));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(6, 0));
        entities.push_back(std::make_pair(8, 2));
    }
    else if(p_rank == 2)
    {
        entities.push_back(std::make_pair(4, 3));
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(8, 1));
    }
    else
    {
        entities.push_back(std::make_pair(2, 0));
        entities.push_back(std::make_pair(4, 2));
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 2));
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    for(int i = 0; i < numNodes; i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
    }

    for(size_t i = 0; i < entities.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, entities[i].first);

        addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, entities[i].second);
        modifiedEntities.push_back(node);
    }

    mesh.my_update_comm_list_based_on_changes_in_comm_map();
    mesh.update_comm_list(modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    // Check the initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 1));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 1));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 1));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_GHOSTED));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_NOT_GHOSTED));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 1));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_SHARED, 3 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_SHARED, 1 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_GHOSTED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 1));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_OWNED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 4), STATE_GHOSTED, 2));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_OWNED, 1));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_SHARED, 0 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_SHARED, 2 ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_NOT_SHARED ));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_NOT_SHARED ));

        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 1), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 2), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 3), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 4), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 5), STATE_NOT_GHOSTED));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 6), STATE_GHOSTED, 0));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 7), STATE_GHOSTED, 2));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 8), STATE_GHOSTED, 1));
        EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 9), STATE_GHOSTED, 1));
    }
}

TEST(CEOME, change_entity_owner_3Elem4Proc1Edge3D)
{
  //  ID.proc
  //                    15.2--------16.2                      15.1--------16.1
  //                     /|          /|                        /|          /|
  //                    / |         / |                       / |         / |
  //                  7.2---------8.2 |                     7.1---------8.1 |
  //                   |  |  3.2   |  |                      |  |  3.1   |  |
  //                   |  |        |  |                      |  |        |  |
  //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
  //         /|        | *|        | /|   -->      /|        | *|        | /|
  //        / |        |* |        |/ |           / |        |* |        |/ |
  //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
  //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
  //       |  |        |  |        |  |          |  |        |  |        |  |
  //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
  //       | /         | /         | /           | /         | /         | /
  //       |/          |/          |/            |/          |/          |/
  //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
  //
  //      (*)edge: 1.0                          (*)edge: 1.1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 4)
  {
    return;
  }

  unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  BulkDataTester mesh(meta, pm);
  const int p_rank = mesh.parallel_rank();

  stk::mesh::Part & elem_part = meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);
  stk::mesh::Part & edge_part = meta.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  meta.commit();

  // 1 elem-id for each proc
  stk::mesh::EntityId proc_elemIDs[] = {1, 2, 3};

  // list of node-ids for each element
  const int nodesPerElem = 8;
  stk::mesh::EntityId elem_nodeIDs[][nodesPerElem] = {
    {1, 2, 5, 4,  9, 10, 13, 12},
    {2, 3, 6, 5, 10, 11, 14, 13},
    {5, 6, 8, 7, 13, 14, 16, 15}
  };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  int shared_nodeIDs_and_procs[][3] = {
      {0, 2, 1}, {0, 5, 1}, {0, 5, 2}, {0, 10, 1}, {0, 13, 1}, {0, 13, 2},                         // proc 0
      {1, 2, 0}, {1, 6, 2}, {1, 5, 0}, {1,  5, 2}, {1, 10, 0}, {1, 14, 2}, {1, 13, 0}, {1, 13, 2}, // proc 1
      {2, 5, 0}, {2, 5, 1}, {2, 6, 1}, {2, 13, 0}, {2, 13, 1}, {2, 14, 1}                          // proc 2
  };
  int numSharedNodeTriples = 20;

  mesh.modification_begin();

  if (p_rank < 3) {
    stk::mesh::EntityId elemID = proc_elemIDs[p_rank];
    Entity elem = stk::mesh::declare_element(mesh, elem_part, elemID, elem_nodeIDs[p_rank]);
    Entity edge = mesh.declare_entity(stk::topology::EDGE_RANK, 1, edge_part);
    std::vector<Entity> nodes;
    nodes.push_back(mesh.get_entity(stk::topology::NODE_RANK, 5));
    nodes.push_back(mesh.get_entity(stk::topology::NODE_RANK, 13));
    stk::mesh::impl::connectEntityToEdge(mesh, elem, edge, nodes);
    mesh.declare_relation(edge, nodes[0], 0);
    mesh.declare_relation(edge, nodes[1], 1);
  }

  for (int proc=0; proc < numSharedNodeTriples; ++proc) {
    if (p_rank == shared_nodeIDs_and_procs[proc][0]) {
      int nodeID = shared_nodeIDs_and_procs[proc][1];
      int sharingProc = shared_nodeIDs_and_procs[proc][2];
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
      mesh.add_node_sharing(node, sharingProc);
    }
  }

  mesh.modification_end();

  // Check the initial state
  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  0   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 1, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  1, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  1, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_GHOSTED,  2) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  0   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 0, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  0, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_SHARED,  2      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  0, 2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_SHARED,  2      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_GHOSTED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_GHOSTED,  2) );
  }
  else if (p_rank == 2) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  0   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  2) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  2) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_GHOSTED) );
  }
  else if (p_rank == 3) {  //knows nothing
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_VALID) );
  }

  std::vector<stk::mesh::EntityProc> entities_to_move;
  if (p_rank == 0) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
    int dest_proc = 3;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);

    elem = mesh.get_entity(stk::topology::EDGE_RANK, 1);
    dest_proc = 1;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
  }
  else if (p_rank == 1) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
    int dest_proc = 0;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  else if (p_rank == 2) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
    int dest_proc = 1;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }

  mesh.change_entity_owner(entities_to_move);

  ////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<int, int> > nodeEntities;
  std::vector<std::pair<int, int> > edgeEntities;

  int numNodes = 16;
  if(p_rank == 0)
  {
      nodeEntities.push_back(std::make_pair(2, 3));
      nodeEntities.push_back(std::make_pair(10, 3));
      nodeEntities.push_back(std::make_pair(5, 1));
      nodeEntities.push_back(std::make_pair(5, 3));
      nodeEntities.push_back(std::make_pair(13, 1));
      nodeEntities.push_back(std::make_pair(13, 3));
      nodeEntities.push_back(std::make_pair(6, 1));
      nodeEntities.push_back(std::make_pair(14, 1));

      edgeEntities.push_back(std::make_pair(1, 1));
      edgeEntities.push_back(std::make_pair(1, 3));
  }
  else if(p_rank == 1)
  {
      nodeEntities.push_back(std::make_pair(5, 0));
      nodeEntities.push_back(std::make_pair(5, 3));
      nodeEntities.push_back(std::make_pair(13, 0));
      nodeEntities.push_back(std::make_pair(13, 3));
      nodeEntities.push_back(std::make_pair(6, 0));
      nodeEntities.push_back(std::make_pair(14, 0));

      edgeEntities.push_back(std::make_pair(1, 0));
      edgeEntities.push_back(std::make_pair(1, 3));
  }
  else if(p_rank == 2)
  {

  }
  else
  {
      nodeEntities.push_back(std::make_pair(2, 0));
      nodeEntities.push_back(std::make_pair(10, 0));
      nodeEntities.push_back(std::make_pair(5, 0));
      nodeEntities.push_back(std::make_pair(5, 1));
      nodeEntities.push_back(std::make_pair(13, 0));
      nodeEntities.push_back(std::make_pair(13, 1));

      edgeEntities.push_back(std::make_pair(1, 0));
      edgeEntities.push_back(std::make_pair(1, 1));
  }

  std::vector<stk::mesh::Entity> modifiedEntities;

  for(int i = 0; i < numNodes; i++)
  {
      stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
      eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
  }

  stk::mesh::EntityKey key(stk::topology::EDGE_RANK, 1);
  eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

  for(size_t i = 0; i < nodeEntities.size(); i++)
  {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeEntities[i].first);

      addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeEntities[i].second);
      modifiedEntities.push_back(node);
  }

  for(size_t i = 0; i < edgeEntities.size(); i++)
  {
      stk::mesh::Entity edge = mesh.get_entity(stk::topology::EDGE_RANK, edgeEntities[i].first);

      addSharingInfo(mesh, edge, stk::mesh::BulkData::SHARED, edgeEntities[i].second);
      modifiedEntities.push_back(edge);
  }

  mesh.my_update_comm_list_based_on_changes_in_comm_map();
  mesh.update_comm_list(modifiedEntities);

  mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  mesh.modification_begin();
  mesh.internal_modification_end_for_change_entity_owner(true, stk::mesh::BulkData::MOD_END_SORT);

  ////////////////////////////////////////////////////////////////////////////

  // Check the final state
  if (p_rank == 0) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  1   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 1, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_SHARED,  3      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  1, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_SHARED,  3      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  1, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_SHARED,  1      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_GHOSTED,  1) );
  }
  else if (p_rank == 1) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_GHOSTED) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  1   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 0, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  0, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  0, 3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_GHOSTED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_GHOSTED) );
  }
  else if (p_rank == 2) {  //knows nothing
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_VALID) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_VALID) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_VALID) );
  }
  else if (p_rank == 3) {
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 2), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::ELEMENT_RANK, 3), STATE_GHOSTED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_OWNED,  1   ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_SHARED, 0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::EDGE_RANK, 1), STATE_NOT_GHOSTED ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_OWNED,  3) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_OWNED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_OWNED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_OWNED,  1) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_SHARED,  0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_SHARED,  0      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_SHARED,  0, 1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_NOT_SHARED      ) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_NOT_SHARED      ) );

    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  1), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  2), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  3), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  4), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  5), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  6), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  7), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  8), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK,  9), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 10), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 11), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 12), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 13), STATE_NOT_GHOSTED) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 14), STATE_GHOSTED,  0) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 15), STATE_GHOSTED,  1) );
    EXPECT_TRUE( check_state_CEOME(mesh, EntityKey(stk::topology::NODE_RANK, 16), STATE_GHOSTED,  1) );
  }
}

bool check_state_CEOME(const stk::mesh::BulkData & mesh, const EntityKey & entityKey, EntityStates state,
                 int p0, int p1, int p2, int p3, int p4, int p5)
{
  // Check to see if the state is as expected for the provided EntityKey.
  //
  // Meaning of the optional processor arguments for various states:
  //     STATE_OWNED: Processor that owns the Entity
  //    STATE_SHARED: List of Processors that we share this entity with
  //   STATE_GHOSTED: Processor that we ghost the Entity from
  //
  std::vector<int> procs;
  if (p0 >= 0) {
    procs.push_back(p0);
  }
  if (p1 >= 0) {
    procs.push_back(p1);
  }
  if (p2 >= 0) {
    procs.push_back(p2);
  }
  if (p3 >= 0) {
    procs.push_back(p3);
  }
  if (p4 >= 0) {
    procs.push_back(p4);
  }
  if (p5 >= 0) {
    procs.push_back(p5);
  }
  std::sort(procs.begin(), procs.end());

  Entity entity = mesh.get_entity(entityKey);
  std::ostringstream oss;

  switch (state) {
    case STATE_VALID:
    {
      if (!procs.empty()) {
        oss << "check_state_CEOME(): Cannot provide processors with validity check." << std::endl;
      }
      if (!mesh.is_valid(entity)) {
        oss << "check_state_CEOME(): Entity " << entityKey << " is not valid when it should have been." << std::endl;
      }
      break;
    }
    case STATE_NOT_VALID:
    {
      if (!procs.empty()) {
        oss << "check_state_CEOME(): Cannot provide processors with STATE_NOT_VALID check." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        oss << "check_state_CEOME(): Entity " << entityKey << " is valid when it shouldn't have been." << std::endl;
      }
      break;
    }
    case STATE_OWNED:
    {
      if (procs.size() != 1u) {
        oss << "check_state_CEOME(): Entities can have only one owner." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        if (procs[0] != mesh.parallel_owner_rank(entity) ) {
          oss << "check_state_CEOME(): Owner of entity " << entityKey << " was proc " << mesh.parallel_owner_rank(entity)
              << " and not proc " << procs[0] << std::endl;
        }
      }
      else {
        oss << "check_state_CEOME(): Can't check ownership of locally-invalid entity." << std::endl;
      }
      break;
    }
    case STATE_SHARED:
    {
      if (procs.empty()) {
        oss << "check_state_CEOME(): Must provide processor(s) with STATE_SHARED check." << std::endl;
      }
      stk::mesh::PairIterEntityComm comm_it = mesh.entity_comm_map_shared(entityKey);
      std::vector<int>::const_iterator procs_it = procs.begin();
      bool lists_match = true;

      if (comm_it.size() != procs.size()) {
        lists_match = false;
      }
      else {
        for ( ; procs_it != procs.end(); ++comm_it, ++procs_it) {
          int comm_proc = comm_it.first->proc;
          int user_proc = *procs_it;
          if (comm_proc != user_proc) {
            lists_match = false;
            break;
          }
        }
      }

      if (!lists_match) {
        oss << "check_state_CEOME(): Entity " << entityKey << " was shared with procs (";
        comm_it = mesh.entity_comm_map_shared(entityKey);
        for ( ; comm_it.first != comm_it.second; ++comm_it) {
          int proc = comm_it.first->proc;
          oss << proc << " ";
        }
        oss << ")" << std::endl
            << "               when it was expected to be shared with procs (";
        procs_it = procs.begin();
        for ( ; procs_it != procs.end(); ++procs_it) {
          oss << *procs_it << " ";
        }
        oss << ")" << std::endl;
      }

      break;
    }
    case STATE_NOT_SHARED:
    {
      if (!procs.empty()) {
        oss << "check_state_CEOME(): Cannot provide processors with STATE_NOT_SHARED check." << std::endl;
      }
      if (!mesh.entity_comm_map_shared(entityKey).empty()) {
        oss << "check_state_CEOME(): Entity " << entityKey << " was shared with procs (";
        stk::mesh::PairIterEntityComm comm_pit = mesh.entity_comm_map_shared(entityKey);
        for ( ; comm_pit.first != comm_pit.second; ++comm_pit) {
          int proc = comm_pit.first->proc;
          oss << proc << " ";
        }
        oss << ") when it shouldn't have been shared." << std::endl;
      }
      break;
    }
    case STATE_GHOSTED:
    {
      if (procs.size() != 1) {
        oss << "check_state_CEOME(): Must provide one processor with STATE_GHOSTED check." << std::endl;
      }
      if (!mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        oss << "check_state_CEOME(): Entity " << entityKey << " was not ghosted from any proc when it should have" << std::endl
            << "               been ghosted from proc " << procs[0] << "." << std::endl;
      }
      else {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        if (owner_rank != procs[0]) {
          oss << "check_state_CEOME(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
              << "               when it should have been ghosted from proc " << procs[0] << "." << std::endl;
        }
      }
      break;
    }
    case STATE_NOT_GHOSTED:
    {
      if (!procs.empty()) {
        oss << "check_state_CEOME(): Cannot provide processors with STATE_NOT_GHOSTED check." << std::endl;
      }
      if (mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        oss << "check_state_CEOME(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
            << "               when it shouldn't have been ghosted." << std::endl;
      }
      break;
    }
  }

  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false;
  }
  else {
    return true;
  }

}

void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem,
                       int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move)
{
    const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
    for(unsigned i = 0; i < bulk.num_nodes(elem); ++i)
    {
        if(bulk.parallel_owner_rank(nodes[i]) == bulk.parallel_rank())
        {
            entities_to_move.push_back(stk::mesh::EntityProc(nodes[i], dest_proc));
        }
    }
}

}

