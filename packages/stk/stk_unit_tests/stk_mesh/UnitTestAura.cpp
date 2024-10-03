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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_mesh/base/Comm.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "UnitTestTextMeshFixture.hpp"

namespace stk { namespace mesh { class Part; } }


namespace stk
{
namespace mesh
{
class FieldBase;
}
}

namespace
{
using stk::unit_test_util::build_mesh;
constexpr stk::mesh::EntityState Unchanged = stk::mesh::EntityState::Unchanged;
constexpr stk::mesh::EntityState Created = stk::mesh::EntityState::Created;
constexpr stk::mesh::EntityState Modified = stk::mesh::EntityState::Modified;
constexpr stk::mesh::EntityState Deleted = stk::mesh::EntityState::Deleted;

void delete_elem2_on_p1(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 1) {
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
    EXPECT_TRUE(bulk.is_valid(elem2));
    bool successfullyDestroyed = bulk.destroy_entity(elem2);
    EXPECT_TRUE(successfullyDestroyed);
  }

  bulk.modification_end();
}

void expect_nodes_1_to_4_in_aura_on_p1(const stk::mesh::BulkData& bulk)
{
  if (bulk.parallel_rank() == 1) {
    for(stk::mesh::EntityId id : {1, 2, 3, 4}) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, id);
      EXPECT_TRUE(bulk.is_valid(node) && bulk.bucket(node).in_aura());
    }
  }
}

void expect_nodes_1_to_8_no_longer_valid_on_p1(const stk::mesh::BulkData& bulk)
{
  if (bulk.parallel_rank() == 1) {
    for(stk::mesh::EntityId id : {1, 2, 3, 4, 5, 6, 7, 8}) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, id);
      EXPECT_FALSE(bulk.is_valid(node));
    }
  }
}

TEST(BulkDataTest, destroyDependentGhostsConnectedToDeletedShared)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x2", bulk);

  expect_nodes_1_to_4_in_aura_on_p1(bulk);
  delete_elem2_on_p1(bulk);
  expect_nodes_1_to_8_no_longer_valid_on_p1(bulk);
}

stk::mesh::Part& setupDavidNobleTestCaseTkt12837(stk::mesh::BulkData& bulk)
{
  //
  //        5____1  1  1____3
  //        |   /  /|\  \   |
  //        |E3/  / | \  \E1|
  //        | /  /E4|E2\  \ |
  //        |/  /___|___\  \|
  //        6   6   4   2   2
  //
  //        P2     P1      P0
  //

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::TRIANGLE_3_2D);
  stk::mesh::Part& nonConformalPart = meta.declare_part("noconform", stk::topology::ELEMENT_RANK);

  meta.commit();

  bulk.modification_begin();

  stk::mesh::EntityIdVector elem1_nodes {1, 2, 3}; // 1
  stk::mesh::EntityIdVector elem2_nodes {1, 4, 2}; // 2
  stk::mesh::EntityIdVector elem3_nodes {6, 1, 5}; // 3
  stk::mesh::EntityIdVector elem4_nodes {6, 4, 1}; // 4

  stk::mesh::EntityId elemId1 = 1; // p0
  stk::mesh::EntityId elemId2 = 2; // p1
  stk::mesh::EntityId elemId3 = 3; // p2
  stk::mesh::EntityId elemId4 = 4; // p1

  if(bulk.parallel_rank() == 0)
  {
    stk::mesh::declare_element(bulk, block_1, elemId1, elem1_nodes);
    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    bulk.add_node_sharing(node1, 1);
    bulk.add_node_sharing(node1, 2);
    bulk.add_node_sharing(node2, 1);
  }
  else if(bulk.parallel_rank() == 1)
  {
    stk::mesh::declare_element(bulk, block_1, elemId2, elem2_nodes);
    stk::mesh::declare_element(bulk, block_1, elemId4, elem4_nodes);

    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
    bulk.add_node_sharing(node1, 2);
    bulk.add_node_sharing(node6, 2);

    bulk.add_node_sharing(node1, 0);
    bulk.add_node_sharing(node2, 0);
  }
  else
  {
    stk::mesh::declare_element(bulk, block_1, elemId3, elem3_nodes);
    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
    bulk.add_node_sharing(node1, 0);
    bulk.add_node_sharing(node1, 1);
    bulk.add_node_sharing(node6, 1);
  }

  bulk.modification_end();

  return nonConformalPart;
}

stk::mesh::Part& setupDavidNobleTestCaseTkt14337(stk::mesh::BulkData& bulk)
{
  //
  //        5____1  1  1____3
  //       /|   /  /|\  \   |
  //      / |E3/  / | \  \E1|
  //     /E5| /  /E4|E2\  \ |
  //    /___|/  /___|___\  \|
  //    7   6   6   4   2   2
  //
  //      P2      P1       P0
  //

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::TRIANGLE_3_2D);
  stk::mesh::Part& nonConformalPart = meta.declare_part("noconform", stk::topology::ELEMENT_RANK);
  stk::mesh::Part& edgePart = meta.declare_part("edge_part", stk::topology::EDGE_RANK);
  stk::mesh::PartVector edgeParts = {&edgePart};

  meta.commit();

  bulk.modification_begin();

  stk::mesh::EntityIdVector elem1_nodes {1, 2, 3}; // 1
  stk::mesh::EntityIdVector elem2_nodes {1, 4, 2}; // 2
  stk::mesh::EntityIdVector elem3_nodes {6, 1, 5}; // 3
  stk::mesh::EntityIdVector elem4_nodes {6, 4, 1}; // 4
  stk::mesh::EntityIdVector elem5_nodes {7, 6, 5}; // 5

  stk::mesh::EntityId elemId1 = 1; // p0
  stk::mesh::EntityId elemId2 = 2; // p1
  stk::mesh::EntityId elemId3 = 3; // p2
  stk::mesh::EntityId elemId4 = 4; // p1
  stk::mesh::EntityId elemId5 = 5; // p2

  if(bulk.parallel_rank() == 0)
  {
    stk::mesh::declare_element(bulk, block_1, elemId1, elem1_nodes);
    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    bulk.add_node_sharing(node1, 1);
    bulk.add_node_sharing(node1, 2);
    bulk.add_node_sharing(node2, 1);
  }
  else if(bulk.parallel_rank() == 1)
  {
    stk::mesh::declare_element(bulk, block_1, elemId2, elem2_nodes);
    stk::mesh::declare_element(bulk, block_1, elemId4, elem4_nodes);

    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
    bulk.add_node_sharing(node1, 2);
    bulk.add_node_sharing(node6, 2);

    bulk.add_node_sharing(node1, 0);
    bulk.add_node_sharing(node2, 0);
  }
  else
  {
    stk::mesh::declare_element(bulk, block_1, elemId3, elem3_nodes);
    stk::mesh::declare_element(bulk, block_1, elemId5, elem5_nodes);
    stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
    bulk.add_node_sharing(node1, 0);
    bulk.add_node_sharing(node1, 1);
    bulk.add_node_sharing(node6, 1);
  }

  bulk.modification_end();

  bulk.initialize_face_adjacent_element_graph();

  bulk.modification_begin();

  if (bulk.parallel_rank() == 1) {
    stk::mesh::Entity elem4 = bulk.get_entity(stk::topology::ELEM_RANK, 4);
    bulk.declare_element_side(elem4, 0, edgeParts);
    bulk.declare_element_side(elem4, 2, edgeParts);
  }
  else if (bulk.parallel_rank() == 2) {
    stk::mesh::Entity elem3 = bulk.get_entity(stk::topology::ELEM_RANK, 3);
    bulk.declare_element_side(elem3, 0, edgeParts);
    bulk.declare_element_side(elem3, 2, edgeParts);
  }

  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Ghosting& ghosting = bulk.create_ghosting("myCustomGhosting");
  stk::mesh::EntityProcVec elemToGhost;
  const int destProc = 0;
  if (bulk.parallel_rank() == 1) {
    elemToGhost.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK, 4), destProc));
  }
  if (bulk.parallel_rank() == 2) {
    elemToGhost.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK, 3), destProc));
  }

  bulk.change_ghosting(ghosting, elemToGhost);
  bulk.modification_end();

  return nonConformalPart;
}

bool isEntityInPart(stk::mesh::BulkData &bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id, const stk::mesh::Part &part)
{
  stk::mesh::Entity entity = bulk.get_entity(rank, id);
  return bulk.bucket(entity).member(part);
}

void put_elements_in_part(stk::mesh::BulkData& bulk,
                          const stk::mesh::EntityIdVector& elemIds,
                          stk::mesh::Part& part)
{
  stk::mesh::PartVector addParts = {&part};
  stk::mesh::PartVector rmParts;
  bulk.modification_begin();

  for(stk::mesh::EntityId elemId : elemIds) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    if (bulk.is_valid(elem) && bulk.parallel_rank() == bulk.parallel_owner_rank(elem)) {
      bulk.change_entity_parts(elem, addParts, rmParts);
    }
  }

  bulk.modification_end();
}

void test_node6_shared_and_aura(const stk::mesh::BulkData& bulk)
{
  stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
  if(bulk.parallel_rank() == 0) {
    EXPECT_TRUE(bulk.bucket(node6).in_aura());
  }
  else if(bulk.parallel_rank() == 1) {
    EXPECT_TRUE(bulk.bucket(node6).owned());
    EXPECT_TRUE(bulk.bucket(node6).shared());
  }
  else {
    EXPECT_TRUE(!bulk.bucket(node6).owned());
    EXPECT_TRUE(bulk.bucket(node6).shared());
  }
}

TEST(BulkDataTest, removeElemPartWithNodeSharedWithOneProcAndAuraToAnotherProc)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(2, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::mesh::Part& nonConformalPart = setupDavidNobleTestCaseTkt12837(bulk);

  {
    put_elements_in_part(bulk, {3}, nonConformalPart);

    EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));
  }

  {
    bulk.modification_begin();
    if(bulk.parallel_rank() == 2)
    {
      stk::mesh::PartVector add_parts;
      stk::mesh::PartVector rm_parts;
      rm_parts.push_back(&nonConformalPart);
      stk::mesh::Entity element_3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
      bulk.change_entity_parts(element_3, add_parts, rm_parts);
    }

    EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));

    bulk.modification_end();

    EXPECT_FALSE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));

    test_node6_shared_and_aura(bulk);
  }
}

TEST(BulkDataTest, removeEdgeFromPartWithNodeSharedWithOneProcAndGhostedToAnotherProc)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(2, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::mesh::Part& nonConformalPart = setupDavidNobleTestCaseTkt14337(bulk);
  stk::mesh::PartVector edgeParts = {bulk.mesh_meta_data().get_part("edge_part")};
  EXPECT_TRUE(!edgeParts.empty());
  EXPECT_TRUE(edgeParts[0] != nullptr);

  {
    put_elements_in_part(bulk, {3,4}, nonConformalPart);

    EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, *edgeParts[0]));
  }

  {
    bulk.modification_begin();
    if(bulk.parallel_rank() == 1)
    {
      stk::mesh::PartVector add_parts;
      stk::mesh::Entity edge31 = bulk.get_entity(stk::topology::EDGE_RANK, 31);
      stk::mesh::Entity edge41 = bulk.get_entity(stk::topology::EDGE_RANK, 41);
      ASSERT_TRUE(bulk.is_valid(edge31));
      ASSERT_TRUE(bulk.is_valid(edge41));
      stk::mesh::Entity elem4 = bulk.get_entity(stk::topology::ELEM_RANK, 4);
      bulk.destroy_relation(elem4, edge31, 2);
      bulk.destroy_entity(edge31);
      bulk.change_entity_parts(edge41, add_parts, edgeParts);
    }
    if(bulk.parallel_rank() == 2)
    {
      stk::mesh::Entity edge31 = bulk.get_entity(stk::topology::EDGE_RANK, 31);
      ASSERT_TRUE(bulk.is_valid(edge31));
      stk::mesh::Entity elem3 = bulk.get_entity(stk::topology::ELEM_RANK, 3);
      bulk.destroy_relation(elem3, edge31, 0);
      bulk.destroy_entity(edge31);
      stk::mesh::PartVector empty_parts;
      bulk.declare_element_side(elem3, 0, empty_parts);
    }

    EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, *edgeParts[0]));

    stk::mesh::EntityProcVec addGhost;
    if (bulk.parallel_rank() == 2) {
      addGhost.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::EDGE_RANK, 31), 0));
    }
    bulk.change_ghosting(*bulk.ghostings()[2], addGhost);

    EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, *edgeParts[0]));
  }
}

void disconnect_elem1_on_proc0(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 0) {
    stk::mesh::EntityId elemId = 1;
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, elemId);

    const unsigned numSharedNodes = 4;
    stk::mesh::EntityId sharedNodeIds[] = {5, 6, 8, 7};
    stk::mesh::ConnectivityOrdinal nodeOrds[] = {4, 5, 6, 7};
    stk::mesh::EntityId newNodeIds[] = {13, 14, 16, 15};

    for(unsigned n=0; n<numSharedNodes; ++n) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, sharedNodeIds[n]);
      bulk.destroy_relation(elem1, node, nodeOrds[n]);
      stk::mesh::Entity newNode = bulk.declare_node(newNodeIds[n]);
      bulk.declare_relation(elem1, newNode, nodeOrds[n]);
    }
  }
  bulk.modification_end();
}

void partially_disconnect_elem1_on_proc0(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 0) {
    stk::mesh::EntityId elemId = 2;
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    bulk.destroy_entity(elem2);

    elemId = 1;
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, elemId);

    const unsigned numSharedNodesToDisconnect = 2;
    stk::mesh::EntityId sharedNodeIds[] = {5, 6};
    stk::mesh::ConnectivityOrdinal nodeOrds[] = {4, 5};
    stk::mesh::EntityId newNodeIds[] = {13, 14};

    for(unsigned n=0; n<numSharedNodesToDisconnect; ++n) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, sharedNodeIds[n]);
      bulk.destroy_relation(elem1, node, nodeOrds[n]);
      bulk.destroy_entity(node);
      stk::mesh::Entity newNode = bulk.declare_node(newNodeIds[n]);
      bulk.declare_relation(elem1, newNode, nodeOrds[n]);
    }
  }
  bulk.modification_end();
}

void disconnect_elem2_on_proc1(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 1) {
    stk::mesh::EntityId elemId = 2;
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, elemId);

    const unsigned numSharedNodes = 4;
    stk::mesh::EntityId sharedNodeIds[] = {5, 6, 8, 7};
    stk::mesh::ConnectivityOrdinal nodeOrds[] = {0, 1, 2, 3};
    stk::mesh::EntityId newNodeIds[] = {13, 14, 16, 15};

    for(unsigned n=0; n<numSharedNodes; ++n) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, sharedNodeIds[n]);
      bulk.destroy_relation(elem2, node, nodeOrds[n]);
      stk::mesh::Entity newNode = bulk.declare_node(newNodeIds[n]);
      bulk.declare_relation(elem2, newNode, nodeOrds[n]);
    }
  }
  bulk.modification_end();
}

void partially_disconnect_elem2_on_proc1(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
    bulk.destroy_entity(elem2);
  }

  if (bulk.parallel_rank() == 1) {
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    bulk.destroy_entity(elem1);

    stk::mesh::EntityId elemId = 2;
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, elemId);

    const unsigned numSharedNodesToDisconnect = 2;
    stk::mesh::EntityId sharedNodeIds[] = {5, 6};
    stk::mesh::ConnectivityOrdinal nodeOrds[] = {0, 1};
    stk::mesh::EntityId newNodeIds[] = {13, 14};

    for(unsigned n=0; n<numSharedNodesToDisconnect; ++n) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, sharedNodeIds[n]);
      bulk.destroy_relation(elem2, node, nodeOrds[n]);
      stk::mesh::Entity newNode = bulk.declare_node(newNodeIds[n]);
      bulk.declare_relation(elem2, newNode, nodeOrds[n]);
    }
  }
  bulk.modification_end();
}

void check_node_states_for_elem(const stk::mesh::BulkData& mesh,
                                stk::mesh::EntityId elemId,
                                stk::mesh::EntityState expectedNodeStates[])
{
  stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elemId);
  ASSERT_TRUE(mesh.is_valid(elem));
  const unsigned numNodes = mesh.num_nodes(elem);
  const stk::mesh::Entity* nodes = mesh.begin_nodes(elem);
  for(unsigned n=0; n<numNodes; ++n) {
    EXPECT_EQ(expectedNodeStates[n], mesh.state(nodes[n]))
        <<"state="<<mesh.state(nodes[n])
       <<" for node "<<mesh.identifier(nodes[n])
      <<" of elem "<<mesh.identifier(elem)
      <<" on proc "<<mesh.parallel_rank()
     <<" expectedNodeState="<<expectedNodeStates[n];
  }
}

void check_elem_state(const stk::mesh::BulkData& mesh,
                      stk::mesh::EntityId elemId,
                      stk::mesh::EntityState expectedState)
{
  stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elemId);
  const bool expectValid = expectedState != Deleted;
  if (expectValid) {
    ASSERT_TRUE(mesh.is_valid(elem))<<"elem "<<elemId;
  }
  EXPECT_EQ(expectedState, mesh.state(elem))
      <<"state="<<mesh.state(elem)
     <<" for elem "<<elemId
    <<" on proc "<<mesh.parallel_rank()
   <<" expectedState="<<expectedState;
}

void confirm_entities_not_valid(const stk::mesh::BulkData& mesh,
                                stk::mesh::EntityRank rank,
                                const stk::mesh::EntityIdVector& ids)
{
  for(stk::mesh::EntityId id : ids) {
    stk::mesh::Entity entity = mesh.get_entity(rank, id);
    EXPECT_FALSE(mesh.is_valid(entity))<<" entity "<<mesh.entity_key(entity)
                                      <<" is valid, expected NOT valid.";
  }
}

void test_aura_disconnect_elem_from_shared_owned_nodes(stk::mesh::BulkData& mesh)
{
  disconnect_elem1_on_proc0(mesh);

  int thisProc = stk::parallel_machine_rank(mesh.parallel());
  stk::mesh::EntityId auraElemId = 2;
  if (thisProc == 1) {
    auraElemId = 1;
  }
  if (thisProc == 0) {
    check_elem_state(mesh, auraElemId, Modified);

    stk::mesh::EntityState expectedAuraElemNodeStates[] = {
      Modified, Modified, Modified, Modified,
      Unchanged, Unchanged, Unchanged, Unchanged
    };
    check_node_states_for_elem(mesh, auraElemId, expectedAuraElemNodeStates);

    stk::mesh::EntityId ownedElemId = 1;
    stk::mesh::Entity ownedElem = mesh.get_entity(stk::topology::ELEM_RANK, ownedElemId);
    EXPECT_EQ(Modified, mesh.state(ownedElem));
    stk::mesh::EntityState expectedNodeStates[] = {
      Unchanged, Unchanged, Unchanged, Unchanged,
      Created, Created, Created, Created
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
  else {
    confirm_entities_not_valid(mesh, stk::topology::ELEM_RANK,
                               stk::mesh::EntityIdVector{auraElemId});
    confirm_entities_not_valid(mesh, stk::topology::NODE_RANK,
                               stk::mesh::EntityIdVector{1, 2, 3, 4});

    stk::mesh::EntityId ownedElemId = 2;
    stk::mesh::EntityState expectedNodeStates[] = {
      Modified, Modified, Modified, Modified,
      Unchanged, Unchanged, Unchanged, Unchanged
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
}

void test_aura_partially_disconnect_elem_from_shared_owned_nodes(stk::mesh::BulkData& mesh)
{
  partially_disconnect_elem1_on_proc0(mesh);

  int thisProc = stk::parallel_machine_rank(mesh.parallel());
  stk::mesh::EntityId auraElemId = 2;
  if (thisProc == 1) {
    auraElemId = 1;
  }
  if (thisProc == 0) {
    check_elem_state(mesh, auraElemId, Created);

    stk::mesh::EntityState expectedAuraElemNodeStates[] = {
      Created, Created, Modified, Modified,
      Modified, Modified, Modified, Modified
    };
    check_node_states_for_elem(mesh, auraElemId, expectedAuraElemNodeStates);

    stk::mesh::EntityId ownedElemId = 1;
    stk::mesh::Entity ownedElem = mesh.get_entity(stk::topology::ELEM_RANK, ownedElemId);
    EXPECT_EQ(Modified, mesh.state(ownedElem));
    stk::mesh::EntityState expectedNodeStates[] = {
      Unchanged, Unchanged, Unchanged, Unchanged,
      Created, Created, Modified, Modified
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
  else {
    check_elem_state(mesh, auraElemId, Created);
    stk::mesh::EntityState expectedAuraElemNodeStates[] = {
      Modified, Modified, Modified, Modified,
      Created, Created, Modified, Modified
    };
    check_node_states_for_elem(mesh, auraElemId, expectedAuraElemNodeStates);

    stk::mesh::EntityId ownedElemId = 2;
    stk::mesh::EntityState expectedNodeStates[] = {
      Modified, Modified, Modified, Modified,
      Modified, Modified, Modified, Modified
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
}

void test_aura_disconnect_elem_from_shared_not_owned_nodes(stk::mesh::BulkData& mesh)
{
  disconnect_elem2_on_proc1(mesh);

  int thisProc = stk::parallel_machine_rank(mesh.parallel());
  stk::mesh::EntityId auraElemId = 2;
  if (thisProc == 1) {
    auraElemId = 1;
  }
  if (thisProc == 0) {
    check_elem_state(mesh, auraElemId, Deleted);
    confirm_entities_not_valid(mesh, stk::topology::NODE_RANK,
                               stk::mesh::EntityIdVector{9, 10, 11, 12});

    stk::mesh::EntityId ownedElemId = 1;
    stk::mesh::Entity ownedElem = mesh.get_entity(stk::topology::ELEM_RANK, ownedElemId);
    EXPECT_EQ(Modified, mesh.state(ownedElem));
    stk::mesh::EntityState expectedNodeStates[] = {
      Unchanged, Unchanged, Unchanged, Unchanged,
      Modified, Modified, Modified, Modified
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
  else {
    confirm_entities_not_valid(mesh, stk::topology::ELEM_RANK,
                               stk::mesh::EntityIdVector{auraElemId});
    confirm_entities_not_valid(mesh, stk::topology::NODE_RANK,
                               stk::mesh::EntityIdVector{1, 2, 3, 4});

    stk::mesh::EntityId ownedElemId = 2;
    stk::mesh::EntityState expectedNodeStates[] = {
      Created, Created, Created, Created,
      Unchanged, Unchanged, Unchanged, Unchanged
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
}

void test_aura_partially_disconnect_elem_from_shared_not_owned_nodes(stk::mesh::BulkData& mesh)
{
  partially_disconnect_elem2_on_proc1(mesh);

  int thisProc = stk::parallel_machine_rank(mesh.parallel());
  stk::mesh::EntityId auraElemId = 2;
  if (thisProc == 1) {
    auraElemId = 1;
  }
  if (thisProc == 0) {
    check_elem_state(mesh, auraElemId, Created);

    stk::mesh::EntityId ownedElemId = 1;
    stk::mesh::Entity ownedElem = mesh.get_entity(stk::topology::ELEM_RANK, ownedElemId);
    EXPECT_EQ(Modified, mesh.state(ownedElem));
    stk::mesh::EntityState expectedNodeStates[] = {
      Modified, Modified, Modified, Modified,
      Modified, Modified, Modified, Modified
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
  else {
    check_elem_state(mesh, auraElemId, Created);

    stk::mesh::EntityId ownedElemId = 2;
    stk::mesh::EntityState expectedNodeStates[] = {
      Created, Created, Modified, Modified,
      Modified, Modified, Modified, Modified
    };
    check_node_states_for_elem(mesh, ownedElemId, expectedNodeStates);
  }
}

class Aura2Hex2Proc : public stk::unit_test_util::MeshFixture
{
public:
  Aura2Hex2Proc()
  {
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
      setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
      //       3----------7----------11
      //      /|         /|         /|
      //     / |        / |        / |
      //    /  |       /  |       /  |
      //   2----------6----------10  |
      //   |   4------|---8------|---12
      //   |  /       |  /       |  /
      //   | /   E1   | /  E2    | /
      //   |/         |/         |/
      //   1----------5----------9
      //       P0         P1
      //  Nodes 5,6,7,8 are shared
      //
    }
  }
};

TEST_F(Aura2Hex2Proc, disconnectElemFromSharedOwnedNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  test_aura_disconnect_elem_from_shared_owned_nodes(get_bulk());
}

TEST_F(Aura2Hex2Proc, partiallyDisconnectElemFromSharedOwnedNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  test_aura_partially_disconnect_elem_from_shared_owned_nodes(get_bulk());
}

TEST_F(Aura2Hex2Proc, disconnectElemFromSharedNotOwnedNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  test_aura_disconnect_elem_from_shared_not_owned_nodes(get_bulk());
}

TEST_F(Aura2Hex2Proc, partiallyDisconnectElemFromSharedNotOwnedNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  test_aura_partially_disconnect_elem_from_shared_not_owned_nodes(get_bulk());
}

void expect_recv_aura(const stk::mesh::BulkData& bulk,
                      stk::mesh::EntityRank rank,
                      const stk::mesh::EntityIdVector& ids)
{
  stk::mesh::Selector select = bulk.mesh_meta_data().aura_part();
  stk::mesh::EntityVector entities;
  stk::mesh::get_selected_entities(select, bulk.buckets(rank), entities);

  EXPECT_EQ(ids.size(), entities.size());

  for(stk::mesh::EntityId id : ids) {
    stk::mesh::Entity entity = bulk.get_entity(rank, id);
    EXPECT_TRUE(bulk.is_valid(entity));
    EXPECT_FALSE(bulk.parallel_owner_rank(entity) == bulk.parallel_rank());
    EXPECT_TRUE(std::binary_search(entities.begin(), entities.end(),
                                   entity, stk::mesh::EntityLess(bulk)))
        <<"P"<<bulk.parallel_rank()<<" expected to find "
       <<bulk.entity_key(entity)<<" in recv-aura but didn't. "
      <<"owned="<<bulk.bucket(entity).owned()
     <<", shared="<<bulk.bucket(entity).shared()<<std::endl;
  }
}

void test_aura_move_elem1_from_proc0_to_proc1(stk::mesh::BulkData& mesh)
{
  const std::string generatedMeshSpec = "generated:1x1x4";
  stk::io::fill_mesh(generatedMeshSpec, mesh);

  //Initial mesh:
  //       3----------7----------11----------15-----------19
  //      /|         /|         /|           /|          /|
  //     / |        / |        / |          / |         / |
  //    /  |       /  |       /  |         /  |        /  |
  //   2----------6----------10-----------14----------18  |
  //   |   4------|---8------|---12-------|--16-------|---20
  //   |  /       |  /       |  /         |  /        |  /
  //   | /   E1   | /  E2    | /    E3    | /    E4   | /
  //   |/         |/         |/           |/          |/
  //   1----------5----------9------------13----------17
  //       P0         P0          P1           P1
  //  Nodes 9,10,11,12 are shared
  //  Elem 3 is aura-ghost on P0 and elem 2 is aura-ghost on P1
  //  Nodes 13-16 are aura-ghosts on P0, nodes 5-8 are aura-ghosts on P1
  //
  int thisProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  {
    stk::mesh::EntityIdVector elemIds[] = {{3}, {2}};
    stk::mesh::EntityIdVector nodeIds[] = {{13,14,15,16}, {5,6,7,8}};
    expect_recv_aura(mesh, stk::topology::ELEM_RANK, elemIds[thisProc]);
    expect_recv_aura(mesh, stk::topology::NODE_RANK, nodeIds[thisProc]);
  }

  //---------------------------------------
  stk::mesh::EntityProcVec elemToMove;
  if (thisProc == 0) {
    stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    elemToMove.push_back(stk::mesh::EntityProc(elem1, 1));
  }

  mesh.change_entity_owner(elemToMove);

  //After change-entity-owner moves elem 1 to P1:
  //       3----------7----------11----------15-----------19
  //      /|         /|         /|           /|          /|
  //     / |        / |        / |          / |         / |
  //    /  |       /  |       /  |         /  |        /  |
  //   2----------6----------10-----------14----------18  |
  //   |   4------|---8------|---12-------|--16-------|---20
  //   |  /       |  /       |  /         |  /        |  /
  //   | /   E1   | /  E2    | /    E3    | /    E4   | /
  //   |/         |/         |/           |/          |/
  //   1----------5----------9------------13----------17
  //       P1         P0          P1           P1
  //  Nodes 1-12 are shared
  //  Elem 2 is aura-ghost on P1 and elems 1 and 3 are aura-ghosts on P0
  //  Nodes 13-16 are aura-ghosts on P0. No nodes are aura-ghosts on P1.
  //
  {
    stk::mesh::EntityIdVector elemIds[] = {{1, 3}, {2}};
    stk::mesh::EntityIdVector nodeIds[] = {{13,14,15,16}, {}};
    expect_recv_aura(mesh, stk::topology::ELEM_RANK, elemIds[thisProc]);
    expect_recv_aura(mesh, stk::topology::NODE_RANK, nodeIds[thisProc]);
  }
}

TEST(BulkData, aura_moveElem1FromProc0ToProc1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD)==2) {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                      .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                                                      .create();
    test_aura_move_elem1_from_proc0_to_proc1(*bulkPtr);
  }
}

TEST(BulkData, aura_moveElem1FromProc0ToProc1_NoUpwardConnectivity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD)==2) {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                      .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                                                      .set_upward_connectivity(false)
                                                      .create();
    test_aura_move_elem1_from_proc0_to_proc1(*bulkPtr);
  }
}

class BulkDataAura : public stk::unit_test_util::MeshFixture
{
public:
  void verify_no_aura()
  {
    verify_aura(false, 0, 0);
  }

  void verify_aura(bool expectAuraOptionIsOn,
                   unsigned expectedNumAuraNodes,
                   unsigned expectedNumAuraElems)
  {
    EXPECT_EQ(expectAuraOptionIsOn, get_bulk().is_automatic_aura_on());

    stk::mesh::Selector selectAura = get_meta().aura_part();
    unsigned numAuraNodes = stk::mesh::count_selected_entities(selectAura, get_bulk().buckets(stk::topology::NODE_RANK));
    unsigned numAuraElements = stk::mesh::count_selected_entities(selectAura, get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(expectedNumAuraNodes, numAuraNodes);
    EXPECT_EQ(expectedNumAuraElems, numAuraElements);
  }
};

TEST_F(BulkDataAura, turnAuraOnAfterConstruction)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {
    return;
  }

  setup_mesh("generated:2x2x2", stk::mesh::BulkData::NO_AUTO_AURA);

  verify_no_aura();

  bool applyImmediately = true;
  get_bulk().set_automatic_aura_option(stk::mesh::BulkData::AUTO_AURA, applyImmediately);

  bool expectAuraOptionIsOn = true;
  unsigned expectedNumAuraNodes = 9;
  unsigned expectedNumAuraElements = 4;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);
}

TEST_F(BulkDataAura, turnAuraOnAfterConstruction_applyAtNextModEnd)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {
    return;
  }

  setup_mesh("generated:2x2x2", stk::mesh::BulkData::NO_AUTO_AURA);

  verify_no_aura();

  bool applyImmediately = false;
  get_bulk().set_automatic_aura_option(stk::mesh::BulkData::AUTO_AURA, applyImmediately);

  bool expectAuraOptionIsOn = true;
  unsigned expectedNumAuraNodes = 0;
  unsigned expectedNumAuraElements = 0;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);

  get_bulk().modification_begin();
  get_bulk().modification_end();

  expectedNumAuraNodes = 9;
  expectedNumAuraElements = 4;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);
}

TEST_F(BulkDataAura, turnAuraOffAfterConstruction)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {
    return;
  }

  setup_mesh("generated:2x2x2", stk::mesh::BulkData::AUTO_AURA);

  bool expectAuraOptionIsOn = true;
  unsigned expectedNumAuraNodes = 9;
  unsigned expectedNumAuraElements = 4;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);

  bool applyImmediately = true;
  get_bulk().set_automatic_aura_option(stk::mesh::BulkData::NO_AUTO_AURA, applyImmediately);

  verify_no_aura();
}

TEST_F(BulkDataAura, turnAuraOffAfterConstruction_applyAtNextModEnd)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {
    return;
  }

  setup_mesh("generated:2x2x2", stk::mesh::BulkData::AUTO_AURA);

  bool expectAuraOptionIsOn = true;
  unsigned expectedNumAuraNodes = 9;
  unsigned expectedNumAuraElements = 4;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);

  bool applyImmediately = false;
  get_bulk().set_automatic_aura_option(stk::mesh::BulkData::NO_AUTO_AURA, applyImmediately);

  expectAuraOptionIsOn = false;
  expectedNumAuraNodes = 9;
  expectedNumAuraElements = 4;
  verify_aura(expectAuraOptionIsOn, expectedNumAuraNodes, expectedNumAuraElements);

  get_bulk().modification_begin();
  get_bulk().modification_end();

  verify_no_aura();
}

class AuraSharedSideMods : public TestTextMeshAura
{
public:
  void verify_num_faces(size_t goldCount)
  {
    EXPECT_EQ(goldCount, stk::mesh::count_entities(get_bulk(), stk::topology::FACE_RANK, get_meta().universal_part()));
  }

  void delete_elem2()
  {
    get_bulk().modification_begin();
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    if (get_bulk().is_valid(elem2) && get_bulk().parallel_owner_rank(elem2) == get_bulk().parallel_rank()) {
      get_bulk().destroy_entity(elem2);
    }
    get_bulk().modification_end();
  }

  void verify_face_owned_on_P0_aura_on_P1(stk::mesh::EntityId faceId)
  {
    stk::mesh::Entity face = get_bulk().get_entity(stk::topology::FACE_RANK, faceId);
    EXPECT_TRUE(get_bulk().is_valid(face));
    if (get_bulk().parallel_rank() == 0) {
      EXPECT_TRUE(get_bulk().parallel_owner_rank(face) == get_bulk().parallel_rank());
    }
    else {
      EXPECT_TRUE(get_bulk().bucket(face).in_aura());
    }
  }

  void verify_nodes_shared(const std::vector<stk::mesh::EntityId>& nodeIds)
  {
    for(stk::mesh::EntityId nodeId : nodeIds) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      EXPECT_TRUE(get_bulk().is_valid(node));
      EXPECT_TRUE(get_bulk().bucket(node).shared());
    }
  }

  void verify_nodes_owned_on_P0_aura_on_P1(const std::vector<stk::mesh::EntityId>& nodeIds)
  {
    const int thisProc = get_bulk().parallel_rank();
    for(stk::mesh::EntityId nodeId : nodeIds) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      EXPECT_TRUE(get_bulk().is_valid(node));
      if (thisProc == 0) {
        EXPECT_TRUE(get_bulk().bucket(node).owned());
        EXPECT_FALSE(get_bulk().bucket(node).in_aura());
      }
      else {
        ASSERT_EQ(1, thisProc);
        EXPECT_FALSE(get_bulk().bucket(node).owned());
        EXPECT_TRUE(get_bulk().bucket(node).in_aura());
      }
    }
  }

  void recreate_elem2_on_P1()
  {
    get_bulk().modification_begin();
    stk::mesh::Entity elem2;
    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Part& hexPart = get_meta().get_topology_root_part(stk::topology::HEX_8);
      stk::mesh::PartVector parts = {&hexPart};
      elem2 = stk::mesh::declare_element(get_bulk(), parts, 2, {5, 6, 7, 8, 9, 10, 11, 12});
    }
    get_bulk().modification_end();
  }

  void declare_side_on_elem2()
  {
    get_bulk().modification_begin();
    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      EXPECT_TRUE(get_bulk().is_valid(elem2));
      const unsigned sideOrdinal = 4;
      get_bulk().declare_element_side<stk::mesh::PartVector>(elem2, sideOrdinal);
    }
    get_bulk().modification_end();
  }

  void declare_relation_on_elem2()
  {
    get_bulk().modification_begin();
    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      EXPECT_TRUE(get_bulk().is_valid(elem2));
      stk::mesh::Entity face16 = get_bulk().get_entity(stk::topology::FACE_RANK, 16);
      const unsigned sideOrdinal = 4;
      get_bulk().declare_relation(elem2, face16, sideOrdinal);
    }
    get_bulk().modification_end();
  }
};

TEST_F(AuraSharedSideMods, sharedFace)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "1, 2, HEX_8,5,6,7,8,9,10,11,12\n"
                         "1, 3, HEX_8,5,13,14,15,16,17,18,19\n"
                         "0, 1, HEX_8,1,2,3,4,5,6,7,8|sideset:data=1,6";
  setup_text_mesh(meshDesc);
  verify_num_elements(3);
  verify_num_faces(1);
}

TEST_F(AuraSharedSideMods, sharedFaceDeleteElemRecreateElem_declareSide)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "1, 2, HEX_8,5,6,7,8,9,10,11,12\n"
                         "1, 3, HEX_8,5,13,14,15,16,17,18,19\n"
                         "0, 1, HEX_8,1,2,3,4,5,6,7,8|sideset:data=1,6";
  setup_text_mesh(meshDesc);
  verify_num_elements(3);
  verify_num_faces(1);

  delete_elem2();
  verify_num_elements(2);
  stk::mesh::EntityId faceId = 16;
  verify_face_owned_on_P0_aura_on_P1(faceId);
  verify_nodes_shared({5});
  verify_nodes_owned_on_P0_aura_on_P1({6, 7, 8});

  recreate_elem2_on_P1();
  declare_side_on_elem2();
  verify_num_faces(1);
  verify_nodes_shared({5, 6, 7, 8});
}

TEST_F(AuraSharedSideMods, sharedFaceDeleteElemRecreateElem_declareRelation)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "1, 2, HEX_8,5,6,7,8,9,10,11,12\n"
                         "1, 3, HEX_8,5,13,14,15,16,17,18,19\n"
                         "0, 1, HEX_8,1,2,3,4,5,6,7,8|sideset:data=1,6";
  setup_text_mesh(meshDesc);
  verify_num_elements(3);
  verify_num_faces(1);

  delete_elem2();
  verify_num_elements(2);
  stk::mesh::EntityId faceId = 16;
  verify_face_owned_on_P0_aura_on_P1(faceId);
  verify_nodes_shared({5});
  verify_nodes_owned_on_P0_aura_on_P1({6, 7, 8});

  recreate_elem2_on_P1();
  declare_relation_on_elem2();
  verify_num_faces(1);
  verify_nodes_shared({5, 6, 7, 8});
}

class AuraTetSide : public TestTextMeshAura
{
public:
  void delete_side_on_p0_check_marking()
  {
    get_bulk().modification_begin();

    stk::mesh::Entity elem2577 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2577);
    stk::mesh::Entity elem2579 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2579);
    stk::mesh::Entity face25772 = get_bulk().get_entity(stk::topology::FACE_RANK, 25772);
    EXPECT_TRUE(get_bulk().is_valid(elem2577));
    EXPECT_TRUE(get_bulk().is_valid(elem2579));
    EXPECT_TRUE(get_bulk().is_valid(face25772));

    if (get_bulk().parallel_rank() == 0) {
      stk::mesh::ConnectivityOrdinal ord = 1;
      EXPECT_TRUE(get_bulk().destroy_relation(elem2577, face25772, ord));
      EXPECT_TRUE(get_bulk().destroy_relation(elem2579, face25772, ord));
      EXPECT_TRUE(get_bulk().destroy_entity(face25772));

      stk::mesh::Entity elem2587 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2587);
      stk::mesh::Entity node484 = get_bulk().get_entity(stk::topology::NODE_RANK, 484);
      stk::mesh::Entity node494 = get_bulk().get_entity(stk::topology::NODE_RANK, 494);
      EXPECT_TRUE(get_bulk().state(elem2587) == stk::mesh::Modified);
      EXPECT_TRUE(get_bulk().state(node484) == stk::mesh::Unchanged);
      EXPECT_TRUE(get_bulk().state(node494) == stk::mesh::Unchanged);
    }

    get_bulk().modification_end();
  }

  void delete_elems_on_p0_check_aura_nodes_deleted()
  {
    stk::mesh::Entity elem2577 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2577);
    stk::mesh::Entity elem2579 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2579);
    EXPECT_TRUE(get_bulk().is_valid(elem2577));
    EXPECT_TRUE(get_bulk().is_valid(elem2579));

    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 0) {
      EXPECT_TRUE(get_bulk().destroy_entity(elem2577));
      EXPECT_TRUE(get_bulk().destroy_entity(elem2579));
    }

    get_bulk().modification_end();

    if (get_bulk().parallel_rank() == 1) {
      elem2577 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2577);
      elem2579 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2579);
      EXPECT_FALSE(get_bulk().is_valid(elem2577));
      EXPECT_FALSE(get_bulk().is_valid(elem2579));

      stk::mesh::Entity node452 = get_bulk().get_entity(stk::topology::NODE_RANK, 452);
      stk::mesh::Entity node464 = get_bulk().get_entity(stk::topology::NODE_RANK, 464);
      stk::mesh::Entity node474 = get_bulk().get_entity(stk::topology::NODE_RANK, 474);
      EXPECT_FALSE(get_bulk().is_valid(node452));
      EXPECT_FALSE(get_bulk().is_valid(node464));
      EXPECT_FALSE(get_bulk().is_valid(node474));
    }
  }

  void delete_elem_on_p1_check_aura_nodes_deleted()
  {
    stk::mesh::Entity elem2601 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2601);
    EXPECT_TRUE(get_bulk().is_valid(elem2601));

    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 1) {
      EXPECT_TRUE(get_bulk().destroy_entity(elem2601));
    }

    get_bulk().modification_end();

    if (get_bulk().parallel_rank() == 1) {
      elem2601 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2601);
      EXPECT_FALSE(get_bulk().is_valid(elem2601));

      stk::mesh::Entity node452 = get_bulk().get_entity(stk::topology::NODE_RANK, 452);
      stk::mesh::Entity node464 = get_bulk().get_entity(stk::topology::NODE_RANK, 464);
      stk::mesh::Entity node474 = get_bulk().get_entity(stk::topology::NODE_RANK, 474);
      EXPECT_FALSE(get_bulk().is_valid(node452));
      EXPECT_FALSE(get_bulk().is_valid(node464));
      EXPECT_FALSE(get_bulk().is_valid(node474));
    }
  }

  void print_face_nodes()
  {
    stk::mesh::Entity face25872 = get_bulk().get_entity(stk::topology::FACE_RANK, 25872);
    EXPECT_TRUE(get_bulk().is_valid(face25872));
    const unsigned numNodes = get_bulk().num_nodes(face25872);
    EXPECT_EQ(3u, numNodes);
    const stk::mesh::Entity* nodes = get_bulk().begin_nodes(face25872);
    std::cerr<<" face 25872 nodes: ";
    for(unsigned i=0; i<numNodes; ++i) {
      std::cerr<<get_bulk().identifier(nodes[i])<<" ";
    }
    std::cerr<<std::endl;
  }

  void delete_aura_elem_on_p1_check_aura_nodes_deleted()
  {
    stk::mesh::Entity elem2587 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2587);
    stk::mesh::Entity elem2601 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2601);
    EXPECT_TRUE(get_bulk().is_valid(elem2587));
    EXPECT_TRUE(get_bulk().is_valid(elem2601));

    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 0) {
      stk::mesh::Entity elem2599 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2599);
      EXPECT_TRUE(get_bulk().destroy_entity(elem2599));
    }
    else {
      EXPECT_TRUE(get_bulk().destroy_entity(elem2601));
    }

    get_bulk().modification_end();

    if (get_bulk().parallel_rank() == 1) {
      elem2587 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2587);
      EXPECT_FALSE(get_bulk().is_valid(elem2587));

      stk::mesh::Entity node484 = get_bulk().get_entity(stk::topology::NODE_RANK, 484);
      stk::mesh::Entity node494 = get_bulk().get_entity(stk::topology::NODE_RANK, 494);
      EXPECT_FALSE(get_bulk().is_valid(node484));
      EXPECT_FALSE(get_bulk().is_valid(node494));
    }
  }

};

TEST_F(AuraTetSide, removeAuraSideOnOwner_checkMarking)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2579, TET_4, 91,452,474,464\n"
                         "0, 2577, TET_4, 50,452,464,474\n"
                         "0, 2587, TET_4, 50,452,484,494\n"
                         "1, 2601, TET_4, 91, 50, 60, 70|sideset:data=2577,2";
  setup_text_mesh(meshDesc);

  delete_side_on_p0_check_marking();
}

TEST_F(AuraTetSide, removeAuraElemsOnNonOwner_checkMarking)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2587, TET_4, 50,452,484,494\n"
                         "0, 2599, TET_4, 50,497,498,499\n"
                         "1, 2601, TET_4, 91, 50, 60, 70|sideset:data=2587,2";
  setup_text_mesh(meshDesc);

  delete_aura_elem_on_p1_check_aura_nodes_deleted();
}

TEST_F(AuraTetSide, printFaceNodes_serial)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2587, TET_4, 50,452,484,494\n"
                         "0, 2599, TET_4, 50,497,498,499\n"
                         "0, 2601, TET_4, 91, 50, 60, 70|sideset:data=2587,2";
  setup_text_mesh(meshDesc);

  print_face_nodes();
}

TEST_F(AuraTetSide, removeElemsOnP0_auraNodesDeletedOnP1)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2579, TET_4, 91,452,474,464\n"
                         "0, 2577, TET_4, 50,452,464,474\n"
                         "1, 2601, TET_4, 91, 50, 60, 70";
  setup_text_mesh(meshDesc);

  delete_elems_on_p0_check_aura_nodes_deleted();
}

TEST_F(AuraTetSide, removeElemOnP1_auraNodesDeletedOnP1)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2579, TET_4, 91,452,474,464\n"
                         "0, 2577, TET_4, 50,452,464,474\n"
                         "1, 2601, TET_4, 91, 50, 60, 70";
  setup_text_mesh(meshDesc);

  delete_elem_on_p1_check_aura_nodes_deleted();
}

TEST_F(AuraTetSide, removeAuraElemOnP2_auraNodesDeletedOnP1)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 2579, TET_4, 91,452,474,464\n"
                         "0, 2577, TET_4, 50,452,464,474\n"
                         "1, 2601, TET_4, 91, 50, 60, 70";
  setup_text_mesh(meshDesc);

  delete_elem_on_p1_check_aura_nodes_deleted();
}

class AuraTetGhosting : public TestTextMeshAura
{
public:

  void make_node110_owned_by_p1()
  {
    stk::mesh::EntityProcVec nodeToMove;
    {
      stk::mesh::Entity node110 = get_bulk().get_entity(stk::topology::NODE_RANK,110);
      if (get_bulk().parallel_rank() == 0) {
        EXPECT_TRUE(get_bulk().bucket(node110).owned());
        int destProc = 1;
        nodeToMove.push_back(stk::mesh::EntityProc(node110, destProc));
      }
    }

    get_bulk().change_entity_owner(nodeToMove);

    {
      stk::mesh::Entity node110 = get_bulk().get_entity(stk::topology::NODE_RANK,110);
      if (get_bulk().parallel_rank() == 1) {
        EXPECT_TRUE(get_bulk().bucket(node110).owned());
      }
    }
  }

  void destroy_elems(const std::vector<std::pair<int,stk::mesh::EntityId>>& procElemIds)
  {
    stk::mesh::EntityVector elemsToDestroy;
    for(const std::pair<int,stk::mesh::EntityId>& procElemId : procElemIds) {
      if (procElemId.first == get_bulk().parallel_rank()) {
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, procElemId.second);
        ASSERT_TRUE(get_bulk().is_valid(elem));
        elemsToDestroy.push_back(elem);
      }
    }

    get_bulk().modification_begin();
    stk::mesh::destroy_elements_no_mod_cycle(get_bulk(), elemsToDestroy, get_meta().universal_part());
    get_bulk().modification_end();
  }

  void custom_ghost(stk::mesh::EntityRank rank, stk::mesh::EntityId entityId, const std::vector<int>& destProcs)
  {
    get_bulk().modification_begin();

    stk::mesh::Ghosting& customGhosting = get_bulk().create_ghosting("myCustomGhosting");

    stk::mesh::Entity entity = get_bulk().get_entity(rank, entityId);
    stk::mesh::EntityProcVec entityToGhost;
    if (get_bulk().is_valid(entity) && get_bulk().parallel_owner_rank(entity)==get_bulk().parallel_rank()) {
      for(int destProc : destProcs) {
        entityToGhost.push_back(stk::mesh::EntityProc(entity, destProc));
      }
    }

    get_bulk().change_ghosting(customGhosting, entityToGhost);

    get_bulk().modification_end();
  }
};

TEST_F(AuraTetGhosting, destroy_elems)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1049, TET_4, 101,140,138,137\n"
                         "0, 1050, TET_4, 140,110,139,136\n"
                         "0, 1051, TET_4, 138,139,108,135\n"
                         "0, 1052, TET_4, 137,136,135,107\n"
                         "0, 1053, TET_4, 138,137,140,136\n"
                         "0, 1054, TET_4, 138,135,137,136\n"
                         "0, 1055, TET_4, 138,139,135,136\n"
                         "0, 1056, TET_4, 138,140,139,136\n"
                         "1, 2000, TET_4, 200,110,162,201\n"
                         "2, 1065, TET_4, 147,108,103,109\n"
                         "2, 1066, TET_4, 105,147,103,109\n"
                         "3, 1108, TET_4, 134,106,133,104\n"
                         "3, 1109, TET_4, 132,133,105,104\n"
                         "3, 1110, TET_4, 132,134,133,104\n"
                         "3, 1111, TET_4, 107,134,132,104\n"
                         "4, 1081, TET_4, 134,136,107,135\n"
                         "4, 1145, TET_4, 135,107,136,134\n"
                         "4, 1146, TET_4, 139,136,110,162\n"
                         "4, 1148, TET_4, 136,135,134,159\n"
                         "4, 1149, TET_4, 136,139,135,159\n"
                         "4, 1150, TET_4, 136,162,139,159\n"
                         "4, 1151, TET_4, 136,134,162,159\n"
                         "5, 1089, TET_4, 108,147,135,159\n"
                         "5, 1090, TET_4, 147,105,132,133\n"
                         "5, 1091, TET_4, 135,132,107,134\n"
                         "5, 1092, TET_4, 159,133,134,106\n"
                         "5, 1093, TET_4, 132,147,133,159\n"
                         "5, 1094, TET_4, 132,135,147,159\n"
                         "5, 1095, TET_4, 132,134,135,159\n"
                         "5, 1096, TET_4, 132,133,134,159\n"
                         "6, 1160, TET_4, 159,133,160,147\n"
                         "6, 1161, TET_4, 159,160,133,106\n"
                         "6, 1162, TET_4, 109,159,160,147\n"
                         "6, 1163, TET_4, 109,159,147,108\n"
                         "6, 1164, TET_4, 105,160,133,147\n"
                         "6, 1165, TET_4, 105,160,147,109\n"
                         "7, 2001, TET_4, 202,162,110,203|sideset:data=1151,4,1146,2,1052,1,1050,1,1081,4,1053,2\n";

  setup_text_mesh(meshDesc);

  destroy_elems({ {0, 1049},
                  {0, 1050},
                  {0, 1051},
                  {0, 1052},
                  {0, 1053},
                  {0, 1054},
                  {0, 1055},
                  {0, 1056},
                  {2, 1065},
                  {2, 1066},
                  {3, 1108},
                  {3, 1109},
                  {3, 1110},
                  {3, 1111},
                  {5, 1089},
                  {5, 1090},
                  {5, 1091},
                  {5, 1092},
                  {5, 1093},
                  {5, 1094},
                  {5, 1095},
                  {5, 1096},
                  {6, 1160},
                  {6, 1161},
                  {6, 1162},
                  {6, 1163},
                  {6, 1164},
                  {6, 1165} });
}

TEST_F(AuraTetGhosting, destroy_elems_custom_ghosting)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1049, TET_4, 101,140,138,137\n"
                         "0, 1050, TET_4, 140,110,139,136\n"
                         "0, 1051, TET_4, 138,139,108,135\n"
                         "0, 1052, TET_4, 137,136,135,107\n"
                         "0, 1053, TET_4, 138,137,140,136\n"
                         "0, 1054, TET_4, 138,135,137,136\n"
                         "0, 1055, TET_4, 138,139,135,136\n"
                         "0, 1056, TET_4, 138,140,139,136\n"
                         "1, 2000, TET_4, 200,110,162,201\n"
                         "2, 1065, TET_4, 147,108,103,109\n"
                         "2, 1066, TET_4, 105,147,103,109\n"
                         "3, 1108, TET_4, 134,106,133,104\n"
                         "3, 1109, TET_4, 132,133,105,104\n"
                         "3, 1110, TET_4, 132,134,133,104\n"
                         "3, 1111, TET_4, 107,134,132,104\n"
                         "4, 1081, TET_4, 134,136,107,135\n"
                         "4, 1145, TET_4, 135,107,136,134\n"
                         "4, 1146, TET_4, 139,136,110,162\n"
                         "4, 1148, TET_4, 136,135,134,159\n"
                         "4, 1149, TET_4, 136,139,135,159\n"
                         "4, 1150, TET_4, 136,162,139,159\n"
                         "4, 1151, TET_4, 136,134,162,159\n"
                         "5, 1089, TET_4, 108,147,135,159\n"
                         "5, 1090, TET_4, 147,105,132,133\n"
                         "5, 1091, TET_4, 135,132,107,134\n"
                         "5, 1092, TET_4, 159,133,134,106\n"
                         "5, 1093, TET_4, 132,147,133,159\n"
                         "5, 1094, TET_4, 132,135,147,159\n"
                         "5, 1095, TET_4, 132,134,135,159\n"
                         "5, 1096, TET_4, 132,133,134,159\n"
                         "6, 1160, TET_4, 159,133,160,147\n"
                         "6, 1161, TET_4, 159,160,133,106\n"
                         "6, 1162, TET_4, 109,159,160,147\n"
                         "6, 1163, TET_4, 109,159,147,108\n"
                         "6, 1164, TET_4, 105,160,133,147\n"
                         "6, 1165, TET_4, 105,160,147,109\n"
                         "7, 2001, TET_4, 202,162,110,203|sideset:data=1151,4,1146,2,1052,1,1050,1,1081,4,1053,2\n";

  setup_text_mesh(meshDesc);

  //This mesh setup provides coverage for some dark corners of BulkData::modification_end.
  //Specifically, node 110 is both aura-ghosted and custom-ghosted from proc 0 to procs 4 and 7,
  //and it is both shared with, and custom-ghosted to proc 1.
  //Similarly, node 135 is shared by procs 0 and 4 and 5, as well as being custom-ghosted
  //from proc 0 to procs 4 and 5.
  //This test is a success if no exception is thrown from BulkData::modification_end.

  custom_ghost(stk::topology::NODE_RANK, 110, {1, 4, 7});
  custom_ghost(stk::topology::NODE_RANK, 135, {4, 5});

  destroy_elems({ {0, 1050},
                  {4, 1081},
                  {4, 1145},
                  {4, 1146},
                  {4, 1148},
                  {4, 1149},
                  {1, 2000} });
}

TEST_F(AuraTetGhosting, destroy_elems_sharing_and_aura_ghosting)
{
  if (get_parallel_size() != 4) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1049, TET_4, 101,140,138,137\n"
                         "0, 1050, TET_4, 140,110,139,136\n"
                         "0, 1051, TET_4, 138,139,108,135\n"
                         "0, 1052, TET_4, 137,136,135,107\n"
                         "1, 2000, TET_4, 200,110,162,201\n"
                         "1, 2001, TET_4, 200,110,165,202\n"
                         "1, 2002, TET_4, 200,203,204,202\n"
                         "2, 1065, TET_4, 200,108,103,104\n"
                         "2, 1066, TET_4, 200,101,105,106\n"
                         "3, 1085, TET_4, 136,189,183,156\n"
                         "3, 1086, TET_4, 136,157,153,186|sideset:data=1050,1,1085,1,1086,2\n";

  //This mesh setup provides coverage for some dark corners of BulkData::modification_end.
  //Specifically, we make node 110 shared by procs 0 and 1, and aura-ghosting to proc 2 and 3.
  //Additionally, we move the ownership of node 110 to proc 1, while the face that it is
  //attached to (face 10501, i.e., the face on side 1 of element 1050) remains owned by proc 0.
  //This test is a success if no exception is thrown from BulkData::modification_end.

  setup_text_mesh(meshDesc);

  make_node110_owned_by_p1();

  destroy_elems({ {2, 1065} });
}

TEST_F(AuraTetGhosting, destroy_elems_sharing_and_aura_and_custom_ghosting)
{
  if (get_parallel_size() != 4) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1049, TET_4, 101,140,138,137\n"
                         "0, 1050, TET_4, 140,110,139,136\n"
                         "0, 1051, TET_4, 138,139,108,135\n"
                         "0, 1052, TET_4, 137,136,135,107\n"
                         "1, 2000, TET_4, 200,110,162,201\n"
                         "1, 2001, TET_4, 200,110,165,202\n"
                         "1, 2002, TET_4, 200,203,204,202\n"
                         "2, 1065, TET_4, 200,108,103,104\n"
                         "2, 1066, TET_4, 200,101,105,106\n"
                         "3, 1085, TET_4, 136,189,183,156\n"
                         "3, 1086, TET_4, 136,157,153,186|sideset:data=1050,1,1085,1,1086,2\n";

  //This mesh setup provides coverage for some dark corners of BulkData::modification_end.
  //Specifically, we make node 110 shared by procs 0 and 1, and aura-ghosting to proc 2 and 3.
  //Additionally, we move the ownership of node 110 to proc 1, while the face that it is
  //attached to (face 10501, i.e., the face on side 1 of element 1050) remains owned by proc 0.
  //Finally, we also custom-ghost elems 1050 and 1051 to procs 2 and 3, as well as node 110 to procs 2 and 3.
  //This test is a success if no exception is thrown from BulkData::modification_end.

  setup_text_mesh(meshDesc);

  make_node110_owned_by_p1();

  custom_ghost(stk::topology::ELEM_RANK, 1050, {2, 3});
  custom_ghost(stk::topology::ELEM_RANK, 1051, {2, 3});
  custom_ghost(stk::topology::NODE_RANK, 110, {2, 3});

  destroy_elems({ {2, 1065} });
}

} // empty namespace

