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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Comm.hpp"
#include "stk_mesh/base/DestroyElements.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "UnitTestTextMeshFixture.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

void verify_local_num_aura_entities(const stk::mesh::BulkData& mesh, stk::mesh::EntityRank rank, unsigned expectedNumAuraEntities)
{
  stk::mesh::EntityVector auraEntities;
  stk::mesh::get_entities(mesh, rank, mesh.mesh_meta_data().aura_part(), auraEntities);
  EXPECT_EQ(expectedNumAuraEntities, auraEntities.size());
}

class TestAura2D : public TestTextMeshAura2d
{
public:
  TestAura2D() {}

  void delete_elem5_on_p1()
  {
    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem5 = get_bulk().get_entity(stk::topology::ELEM_RANK,5);
      ASSERT_TRUE(get_bulk().is_valid(elem5));
      EXPECT_TRUE(get_bulk().destroy_entity(elem5));
    }

    get_bulk().modification_end();
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
};

TEST_F(TestAura2D, quadKeyhole)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, QUAD_4_2D, 1,2,6,5\n"
                         "0, 2, QUAD_4_2D, 2,3,7,6\n"
                         "0, 3, QUAD_4_2D, 3,4,8,7\n"
                         "1, 4, QUAD_4_2D, 5,6,10,9\n"
                         "1, 5, QUAD_4_2D, 6,7,11,10\n"
                         "1, 6, QUAD_4_2D, 7,8,12,11\n"
                         "1, 7, QUAD_4_2D, 9,10,14,13\n"
                         "1, 8, QUAD_4_2D, 10,11,15,14\n"
                         "1, 9, QUAD_4_2D, 11,12,16,15\n";

  setup_text_mesh(meshDesc);

  delete_elem5_on_p1();

  if (get_bulk().parallel_rank() == 0) {
    verify_local_num_aura_entities(get_bulk(), stk::topology::ELEM_RANK, 2);
    verify_local_num_aura_entities(get_bulk(), stk::topology::NODE_RANK, 4);
  }
  else {
    verify_local_num_aura_entities(get_bulk(), stk::topology::ELEM_RANK, 3);
    verify_local_num_aura_entities(get_bulk(), stk::topology::NODE_RANK, 4);
  }
}

TEST_F(TestAura2D, triSharingGhosting)
{
  if (get_parallel_size() != 4) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,2,3\n"
                         "0, 2, TRI_3_2D, 3,2,4\n"
                         "0, 8, TRI_3_2D, 3,1,10\n"
                         "1, 3, TRI_3_2D, 2,5,4\n"
                         "1, 4, TRI_3_2D, 2,6,5\n"
                         "1, 9, TRI_3_2D, 4,5,11\n"
                         "2, 5, TRI_3_2D, 2,7,6\n"
                         "2, 6, TRI_3_2D, 1,7,2\n"
                         "3, 7, TRI_3_2D, 6,8,5|sideset:data=3,1,4,1,5,1\n";

  //This test provides coverage for some dark corners of BulkData::modification_end.
  //Specifically, it sets up a combination of sharing and aura-ghosting such that
  //proc 3 has an aura-ghost of node 2, and doesn't know that procs 1 and 2 share
  //that node. Correspondingly, procs 1 and 2 know about each other, but don't know
  //that proc 3 has a recv-ghost of node 2.
  //The handling of this type of situation will get easier if/when we make ghosting
  //info symmetric in BulkData.
 
  setup_text_mesh(meshDesc);

  destroy_elems({ {0,8}, {1,9} });
}

class AuraToSharedToAura : public TestTextMeshAura2d
{
public:
  AuraToSharedToAura() {}

  void create_elem3_p1_and_delete_elem1_p0()
  {
    const stk::mesh::MetaData& meta = get_meta();
    stk::mesh::PartVector triParts = {&meta.get_topology_root_part(stk::topology::TRI_3_2D)};

    const stk::mesh::EntityId elemId = 3;
    stk::mesh::EntityIdVector nodeIds = {1,2,4};

    const int otherProc = 1 - get_bulk().parallel_rank();
    stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

    EXPECT_TRUE(get_bulk().is_valid(node1));
    EXPECT_TRUE(get_bulk().in_ghost(get_bulk().aura_ghosting(), node1));
    EXPECT_TRUE(0 == get_bulk().parallel_owner_rank(node1));

    get_bulk().modification_begin();

    get_bulk().add_node_sharing(node1, otherProc);

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::declare_element(get_bulk(), triParts, elemId, nodeIds);
    }

    bool destroyedOwnedNode = false;
    if (get_bulk().parallel_rank() == 0) {
      stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
      get_bulk().destroy_entity(elem1);
      destroyedOwnedNode = get_bulk().destroy_entity(node1);
    }

    get_bulk().modification_end();

    if (destroyedOwnedNode) {
      //On P0, have to re-get node1 since it was owned, then destroyed,
      //now exists again locally as a recv-aura node.
      node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);
    }
    else {
      //On P1, node1 is still valid: it transitioned from recv-aura to
      //shared, to owned send-aura
    }
    EXPECT_TRUE(get_bulk().is_valid(node1));
    EXPECT_TRUE(get_bulk().in_ghost(get_bulk().aura_ghosting(), node1));
    EXPECT_FALSE(get_bulk().in_shared(node1));
    EXPECT_TRUE(1 == get_bulk().parallel_owner_rank(node1));

    const stk::mesh::ElemElemGraph& eeGraph = get_bulk().get_face_adjacent_element_graph();

    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    ASSERT_TRUE(get_bulk().is_valid(elem3));

    if (get_bulk().parallel_rank() == 0) {
      EXPECT_TRUE(get_bulk().bucket(elem3).in_aura());
      constexpr bool requireValidId = false;
      stk::mesh::impl::LocalId elem3LocalId = eeGraph.get_local_element_id(elem3, requireValidId);
      EXPECT_EQ(stk::mesh::impl::INVALID_LOCAL_ID, elem3LocalId);
    }

    if (get_bulk().parallel_rank() == 1) {
      EXPECT_TRUE(get_bulk().bucket(elem3).owned());
      constexpr bool requireValidId = true;
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      stk::mesh::impl::LocalId elem2LocalId = eeGraph.get_local_element_id(elem2, requireValidId);
      stk::mesh::impl::LocalId elem3LocalId = eeGraph.get_local_element_id(elem3, requireValidId);
      constexpr size_t numConnectedElems = 1;
      EXPECT_EQ(numConnectedElems, eeGraph.get_num_connected_elems(elem3));
      stk::mesh::GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(elem3LocalId);
      ASSERT_EQ(1u, graphEdges.size());
      const stk::mesh::GraphEdge& graphEdge = graphEdges[0];
      EXPECT_EQ(elem3LocalId, graphEdge.elem1());
      EXPECT_EQ(elem2LocalId, graphEdge.elem2());
    }
  }
};

TEST_F(AuraToSharedToAura, makeAuraNodeSharedThenDelete)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,3,2\n"
                         "1, 2, TRI_3_2D, 2,3,4\n";
  setup_text_mesh(meshDesc);
  create_elem3_p1_and_delete_elem1_p0();
}

class Aura2DTri : public TestTextMeshAura2d
{
public:
  Aura2DTri() {}

  void disconnect_node3_on_p1()
  {
    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity node3 = get_bulk().get_entity(stk::topology::NODE_RANK, 3);
      ASSERT_TRUE(get_bulk().is_valid(node3));

      stk::mesh::Entity node6 = get_bulk().declare_node(6);
      stk::mesh::ConnectivityOrdinal ord = 0;

      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      ASSERT_TRUE(get_bulk().is_valid(elem2));

      EXPECT_TRUE(get_bulk().destroy_relation(elem2, node3, ord));
      get_bulk().declare_relation(elem2, node6, ord);
    }
    get_bulk().modification_end();
  }
};

TEST_F(Aura2DTri, destroyRelation_correctAura)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,2,3\n"
                         "1, 2, TRI_3_2D, 3,4,5\n";
  setup_text_mesh(meshDesc);

  verify_local_num_aura_entities(get_bulk(), stk::topology::ELEM_RANK, 1);
  verify_local_num_aura_entities(get_bulk(), stk::topology::NODE_RANK, 2);

  disconnect_node3_on_p1();

  verify_local_num_aura_entities(get_bulk(), stk::topology::ELEM_RANK, 0);
  verify_local_num_aura_entities(get_bulk(), stk::topology::NODE_RANK, 0);
}

class Aura2DTri4Procs : public TestTextMeshAura2d
{
public:
  Aura2DTri4Procs() {}

  void p1_starts_sharing_node2()
  {
    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
      stk::mesh::Entity node8 = get_bulk().get_entity(stk::topology::NODE_RANK, 8);
      stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
      stk::mesh::ConnectivityOrdinal ord = 2;
      get_bulk().destroy_relation(elem4, node8, ord);
      get_bulk().declare_relation(elem4, node2, ord);
    }

    get_bulk().modification_end();
  }

  void p2_starts_sharing_node2()
  {
    get_bulk().modification_begin();

    if (get_bulk().parallel_rank() == 2) {
      stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
      stk::mesh::Entity node10 = get_bulk().get_entity(stk::topology::NODE_RANK, 10);
      stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
      stk::mesh::ConnectivityOrdinal ord = 2;
      get_bulk().destroy_relation(elem4, node10, ord);
      get_bulk().declare_relation(elem4, node2, ord);
    }

    get_bulk().modification_end();
  }
};

TEST_F(Aura2DTri4Procs, addSharerToAlreadySharedAuraNode)
{
  if (get_parallel_size() != 4) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,2,3\n"
                         "1, 3, TRI_3_2D, 1,6,7\n"
                         "1, 4, TRI_3_2D, 1,7,8\n"
                         "2, 5, TRI_3_2D, 9,10,4\n"
                         "3, 2, TRI_3_2D, 2,4,5\n";
  setup_text_mesh(meshDesc);

  stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
  EXPECT_TRUE(get_bulk().is_valid(node2));
  const int proc = get_bulk().parallel_rank();
  EXPECT_EQ(0, get_bulk().parallel_owner_rank(node2));
  if (proc == 1 || proc == 2) {
    EXPECT_TRUE(get_bulk().bucket(node2).in_aura());
    EXPECT_FALSE(get_bulk().bucket(node2).shared());
  }
  else {
    EXPECT_TRUE(get_bulk().bucket(node2).shared());
    EXPECT_FALSE(get_bulk().bucket(node2).in_aura());
  }

  p1_starts_sharing_node2();

  node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
  EXPECT_TRUE(get_bulk().is_valid(node2));
  EXPECT_EQ(0, get_bulk().parallel_owner_rank(node2));

  if (proc == 2) {
    EXPECT_TRUE(get_bulk().bucket(node2).in_aura());
    EXPECT_FALSE(get_bulk().bucket(node2).shared());
  }
  else {
    EXPECT_TRUE(get_bulk().bucket(node2).shared());
    EXPECT_FALSE(get_bulk().bucket(node2).in_aura());
  }

  if (proc == 1 || proc == 3) {
    int ownerProc = 0;
    int otherSharingProc = proc==1 ? 3 : 1;
    EXPECT_TRUE(get_bulk().in_shared(node2, ownerProc));
    EXPECT_TRUE(get_bulk().in_shared(node2, otherSharingProc));
  }
}

TEST_F(Aura2DTri4Procs, addSharerToAlreadySharedAuraNode_addsRemoteGhostingRequestForAlreadyShared)
{
  if (get_parallel_size() != 3) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,2,3\n"
                         "1, 2, TRI_3_2D, 2,4,5\n"
                         "2, 4, TRI_3_2D, 8,9,10\n"
                         "2, 3, TRI_3_2D, 6,7,4\n";
  setup_text_mesh(meshDesc);

  stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
  EXPECT_TRUE(get_bulk().is_valid(node2));
  const int proc = get_bulk().parallel_rank();
  EXPECT_EQ(0, get_bulk().parallel_owner_rank(node2));
  if (proc == 2) {
    EXPECT_TRUE(get_bulk().bucket(node2).in_aura());
    EXPECT_FALSE(get_bulk().bucket(node2).shared());
  }
  else {
    EXPECT_TRUE(get_bulk().bucket(node2).shared());
    EXPECT_FALSE(get_bulk().bucket(node2).in_aura());
  }

  p2_starts_sharing_node2();

  node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
  EXPECT_TRUE(get_bulk().is_valid(node2));
  EXPECT_EQ(0, get_bulk().parallel_owner_rank(node2));

  EXPECT_TRUE(get_bulk().bucket(node2).shared());
  EXPECT_FALSE(get_bulk().bucket(node2).in_aura());

  if (proc == 1 || proc == 2) {
    int ownerProc = 0;
    int otherSharingProc = proc==1 ? 2 : 1;
    EXPECT_TRUE(get_bulk().in_shared(node2, ownerProc));
    EXPECT_TRUE(get_bulk().in_shared(node2, otherSharingProc));
  }
}

} // empty namespace

