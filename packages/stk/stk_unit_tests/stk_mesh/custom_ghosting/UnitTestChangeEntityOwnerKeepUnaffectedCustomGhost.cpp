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

#include <gtest/gtest.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>


class ChangeEntityOwnerTest : public ::testing::Test
{
public:
  ChangeEntityOwnerTest()
    : meta(nullptr), bulk(), customGhosting(nullptr)
  {
  }

  virtual ~ChangeEntityOwnerTest()
  {
  }

  void setup_mesh(stk::ParallelMachine pm)
  {
    stk::mesh::MeshBuilder builder(pm);
    builder.set_spatial_dimension(2);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    bulk = builder.create();
    meta = &(bulk->mesh_meta_data());
    create_1_beam_per_proc();
    custom_ghost_node2_to_proc2();
  }

  void create_1_beam_per_proc()
  {
    const int myProc = stk::parallel_machine_rank(bulk->parallel());
    stk::mesh::Part& beamPart = meta->get_topology_root_part(stk::topology::BEAM_2);

    bulk->modification_begin();

    stk::mesh::EntityId firstNodeId = myProc + 1;
    stk::mesh::Entity node1 = bulk->declare_node(firstNodeId);
    stk::mesh::Entity node2 = bulk->declare_node(firstNodeId+1);
    if (bulk->parallel_rank()==0) {
      int otherProc = 1;
      bulk->add_node_sharing(node2, otherProc);
    }
    if (bulk->parallel_rank()==1) {
      int otherProc = 0;
      bulk->add_node_sharing(node1, otherProc);
      otherProc = 2;
      bulk->add_node_sharing(node2, otherProc);
    }
    if (bulk->parallel_rank()==2) {
      int otherProc = 1;
      bulk->add_node_sharing(node1, otherProc);
    }

    stk::mesh::EntityId beamId = myProc + 1;

    stk::mesh::Entity beam = bulk->declare_element(beamId, stk::mesh::ConstPartVector{&beamPart});
    bulk->declare_relation(beam, node1, 0);
    bulk->declare_relation(beam, node2, 1);

    bulk->modification_end();
  }

  void custom_ghost_node2_to_proc2()
  {
    bulk->modification_begin();
    customGhosting = &bulk->create_ghosting("customGhosting");
    std::vector<stk::mesh::EntityProc> sendGhost;
    if (bulk->parallel_rank() == 0) {
      stk::mesh::EntityId nodeId = 2;
      const int otherProc = 2;
      sendGhost.push_back(stk::mesh::EntityProc(bulk->get_entity(stk::topology::NODE_RANK, nodeId), otherProc));
    }
    bulk->change_ghosting(*customGhosting, sendGhost);
    bulk->modification_end();
  }

  void move_beam2_to_proc0()
  {
    std::vector<stk::mesh::EntityProc> elemToMove;

    if (bulk->parallel_rank()==1) {
      const int otherProc = 0;
      stk::mesh::EntityId beamId = 2;
      elemToMove.push_back(stk::mesh::EntityProc(bulk->get_entity(stk::topology::ELEM_RANK, beamId), otherProc));
    }

    bulk->change_entity_owner(elemToMove);
  }

protected:
  stk::mesh::MetaData* meta = nullptr;
  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::Ghosting * customGhosting;
};

TEST_F(ChangeEntityOwnerTest, keepUnaffectedCustomGhost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

  setup_mesh(MPI_COMM_WORLD);

  stk::mesh::EntityKey node2key(stk::topology::NODE_RANK, 2);
  stk::mesh::Entity node2 = bulk->get_entity(node2key);
  EXPECT_TRUE(bulk->is_valid(node2));

  if (bulk->parallel_rank() == 0) {
    EXPECT_TRUE(bulk->in_send_ghost(*customGhosting, node2key, 2));
  }

  move_beam2_to_proc0();

  if (bulk->parallel_rank() == 0) {
    EXPECT_TRUE(bulk->in_send_ghost(*customGhosting, node2key, 2));
  }
}

