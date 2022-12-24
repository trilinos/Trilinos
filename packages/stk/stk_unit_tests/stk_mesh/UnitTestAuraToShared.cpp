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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Comm.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "UnitTestTextMeshFixture.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

class AuraToSharedToAura : public TestTextMeshAura2d
{
public:
  AuraToSharedToAura() {}

  void create_elem3_p1_and_delete_elem1_p0()
  {
    const stk::mesh::MetaData& meta = get_meta();
    stk::mesh::PartVector triParts = {&meta.get_topology_root_part(stk::topology::TRI_3_2D)};

    stk::mesh::EntityId elemId = 3;
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

