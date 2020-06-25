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
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
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

constexpr stk::mesh::EntityState Unchanged = stk::mesh::EntityState::Unchanged;
constexpr stk::mesh::EntityState Created   = stk::mesh::EntityState::Created;
constexpr stk::mesh::EntityState Modified  = stk::mesh::EntityState::Modified;

stk::mesh::Part& setupDavidNobleTestCase(stk::mesh::BulkData& bulk)
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

bool isEntityInPart(stk::mesh::BulkData &bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id, const stk::mesh::Part &part)
{
    stk::mesh::Entity entity = bulk.get_entity(rank, id);
    const stk::mesh::PartVector &partsNodes6 = bulk.bucket(entity).supersets();
    bool isInPart = false;
    for(size_t i = 0; i < partsNodes6.size(); ++i)
    {
        if(partsNodes6[i] == &part)
        {
            isInPart = true;
            break;
        }
    }
    return isInPart;
}

TEST(BulkDataTest, testRemovingPartsOnNodeSharedWithOneProcAndAuraToAnotherProc)
{
    //unit test for ticket #12837

    stk::mesh::MetaData meta_data(2);
    stk::mesh::BulkData bulk(meta_data, MPI_COMM_WORLD);

    int num_procs = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if(num_procs == 3)
    {
        stk::mesh::Part& nonConformalPart = setupDavidNobleTestCase(bulk);

        {
            bulk.modification_begin();
            stk::mesh::PartVector add_parts;
            stk::mesh::PartVector rm_parts;
            add_parts.push_back(&nonConformalPart);
            if(bulk.parallel_rank() == 2)
            {
                stk::mesh::Entity element_3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
                bulk.change_entity_parts(element_3, add_parts, rm_parts);
            }
            bulk.modification_end();

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

            EXPECT_TRUE(!isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));

            stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
            if(bulk.parallel_rank() == 0)
            {
                EXPECT_TRUE(bulk.bucket(node6).in_aura());
            }
            else if(bulk.parallel_rank() == 1)
            {
                EXPECT_TRUE(bulk.bucket(node6).owned());
                EXPECT_TRUE(bulk.bucket(node6).shared());
            }
            else
            {
                EXPECT_TRUE(!bulk.bucket(node6).owned());
                EXPECT_TRUE(bulk.bucket(node6).shared());
            }
        }

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
                       <<" on proc "<<mesh.parallel_rank()
                       <<" expectedNodeState="<<expectedNodeStates[n];
  }
}

void check_elem_state(const stk::mesh::BulkData& mesh,
                      stk::mesh::EntityId elemId,
                      stk::mesh::EntityState expectedState)
{
  stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elemId);
  ASSERT_TRUE(mesh.is_valid(elem));
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

TEST(BulkData, aura_disconnectElemOnProcBoundary)
{
  int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  int thisProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  if (numProcs==2)
  {
    stk::mesh::MetaData meta(3);
    //stk::mesh::Part& elemPart = meta.declare_part("myElemPart", stk::topology::ELEM_RANK);
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);

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
    const std::string generatedMeshSpec = "generated:1x1x2";
    stk::io::fill_mesh(generatedMeshSpec, mesh);

    disconnect_elem1_on_proc0(mesh);

    stk::mesh::EntityId auraElemId = 2;
    if (thisProc == 1) {
      auraElemId = 1;
    }
    if (thisProc == 0) {
      check_elem_state(mesh, auraElemId, Modified);

      stk::mesh::EntityState expectedAuraElemNodeStates[] = {
        Modified, Modified, Modified, Modified,
        Modified, Modified, Modified, Modified
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

TEST(BulkData, aura_moveElem1FromProc0ToProc1)
{
  int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  int thisProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  if (numProcs==2)
  {
    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);
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
}

} // empty namespace

