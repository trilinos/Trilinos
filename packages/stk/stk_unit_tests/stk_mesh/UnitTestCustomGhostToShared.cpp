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
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/DestroyElements.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "UnitTestTextMeshFixture.hpp"

namespace {

class TestCustomGhostingThenAura2D : public TestTextMesh2d
{
public:
  TestCustomGhostingThenAura2D() {}

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

    stk::mesh::destroy_elements(get_bulk(), elemsToDestroy, get_meta().universal_part());
  }
};

TEST_F(TestCustomGhostingThenAura2D, triSharingGhosting)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, TRI_3_2D, 1,2,3\n"
                         "0, 4, TRI_3_2D, 2,5,3\n"
                         "1, 2, TRI_3_2D, 2,4,5\n"
                         "1, 3, TRI_3_2D, 4,6,7\n"
                         "1, 5, TRI_3_2D, 4,7,5\n";

  setup_text_mesh(meshDesc);

  get_bulk().modification_begin();
  stk::mesh::Ghosting& myGhosting = get_bulk().create_ghosting("myCustomGhosting");
  std::vector<stk::mesh::EntityProc> nodesToGhost;
  if (get_bulk().parallel_rank() == 1) {
    stk::mesh::Entity node4 = get_bulk().get_entity(stk::topology::NODE_RANK,4);
    stk::mesh::Entity node7 = get_bulk().get_entity(stk::topology::NODE_RANK,7);
    nodesToGhost.push_back(stk::mesh::EntityProc(node4,0));
    nodesToGhost.push_back(stk::mesh::EntityProc(node7,0));
  }
  get_bulk().change_ghosting(myGhosting, nodesToGhost);
  get_bulk().modification_end();

  stk::mesh::Entity node4 = get_bulk().get_entity(stk::topology::NODE_RANK,4);
  stk::mesh::Entity node7 = get_bulk().get_entity(stk::topology::NODE_RANK,7);
  EXPECT_TRUE(get_bulk().is_valid(node4));
  EXPECT_TRUE(get_bulk().is_valid(node7));

  constexpr bool applyImmediately = true;
  get_bulk().set_automatic_aura_option(stk::mesh::BulkData::AUTO_AURA, applyImmediately);
  EXPECT_TRUE(get_bulk().is_valid(node4));

  get_bulk().modification_begin();
  if (get_bulk().parallel_rank() == 0) {
    stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK,2);
    EXPECT_TRUE(get_bulk().is_valid(node4));
    stk::mesh::Entity node5 = get_bulk().get_entity(stk::topology::NODE_RANK,5);
    const stk::mesh::MetaData& meta = get_meta();
    stk::mesh::PartVector elemParts =
      {&meta.get_topology_root_part(stk::topology::TRI_3_2D), meta.get_part("block_TRIANGLE_3_2D")};
    EXPECT_TRUE(elemParts[0] != nullptr);
    EXPECT_TRUE(elemParts[1] != nullptr);

    stk::mesh::Entity elem6 = get_bulk().declare_element(6, elemParts);
    get_bulk().declare_relation(elem6, node2, 0);
    get_bulk().declare_relation(elem6, node4, 1);
    get_bulk().declare_relation(elem6, node5, 2);

    stk::mesh::Entity elem7 = get_bulk().declare_element(7, elemParts);
    get_bulk().declare_relation(elem7, node4, 0);
    get_bulk().declare_relation(elem7, node7, 1);
    get_bulk().declare_relation(elem7, node5, 2);
  }

  if (get_bulk().parallel_rank() == 1) {
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK,2);
    EXPECT_TRUE(get_bulk().destroy_entity(elem2));
    stk::mesh::Entity elem5 = get_bulk().get_entity(stk::topology::ELEM_RANK,5);
    EXPECT_TRUE(get_bulk().destroy_entity(elem5));
  }
  get_bulk().modification_end();

  if (get_bulk().parallel_rank() == 0) {
    const bool node4IsSharedAndRecvGhost = get_bulk().bucket(node4).shared() && get_bulk().in_receive_ghost(node4);
    const bool node7IsSharedAndRecvGhost = get_bulk().bucket(node7).shared() && get_bulk().in_receive_ghost(node7);
    EXPECT_TRUE(node4IsSharedAndRecvGhost && node7IsSharedAndRecvGhost);
  }

  EXPECT_NO_THROW(destroy_elems({ {0,6} }));
}

} // anonymous namespace

