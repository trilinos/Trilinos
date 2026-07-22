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


class ParticleCustomGhostTester : public ::testing::Test
{
public:
  ParticleCustomGhostTester()
    : meta(nullptr), bulk()
  {}

  void initialize(stk::ParallelMachine pm, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    stk::mesh::MeshBuilder builder(pm);
    builder.set_spatial_dimension(3);
    builder.set_aura_option(auraOption);
    bulk = builder.create();
    meta = &(bulk->mesh_meta_data());
    createParts();
  }

  void createParticleOnProc0AndGhostToProc1()
  {
    bulk->modification_begin();
    createParticleOnProc0();
    ghostParticleToProc1();
    bulk->modification_end();
  }

  void createBeamOnProc1AndUpdateSharingOnProc0()
  {
    bulk->modification_begin();
    createBeamOnProc1();
    addNode1Sharing();
    bulk->modification_end();
  }

  void destroyBeamOnProc1AndVerifyParticleIsStillOnProc1()
  {
    destroyBeamOnProc1();
    verifyParticleIsStillOnProc1();
    verifyNode1IsSendGhostAndNotSharedOnProc0();
    verifyNode1IsGhostedAndNotSharedOnProc1();
  }

  void verifyNode1NotShared()
  {
    EXPECT_FALSE(bulk->bucket(node1).shared());
    EXPECT_FALSE(bulk->in_shared(bulk->entity_key(node1)));
  }

  void verifyNode1IsSendGhostAndNotSharedOnProc0()
  {
    if (bulk->parallel_rank() == 0)
    {
      verifyNode1NotShared();
      EXPECT_TRUE(bulk->in_send_ghost(bulk->entity_key(node1)));
    }
  }

  void verifyNode1IsGhostedAndNotSharedOnProc1()
  {
    if (bulk->parallel_rank() == 1)
    {
      verifyNode1NotShared();
      EXPECT_TRUE(bulk->in_receive_custom_ghost(bulk->entity_key(node1)));
    }
  }

  void verifyParticleIsStillOnProc1()
  {
    if (bulk->parallel_rank() == 1)
    {
      particle = bulk->get_entity(stk::topology::ELEMENT_RANK,1);
      EXPECT_TRUE(bulk->is_valid(particle));
      EXPECT_TRUE(bulk->is_valid(node1));
    }
  }

  ~ParticleCustomGhostTester()
  {
  }
private:
  void destroyBeamOnProc1()
  {
    bulk->modification_begin();
    if (bulk->parallel_rank() == 1)
    {
      destroyBeam();
    }
    bulk->modification_end();
  }

  void ghostParticleToProc1()
  {
    stk::mesh::Ghosting & my_ghosting = bulk->create_ghosting("My Custom Ghosting");
    stk::mesh::EntityProcVec ghost_entities;
    if (bulk->parallel_rank() == 0)
      ghost_entities.push_back(stk::mesh::EntityProc(particle,1));
    bulk->change_ghosting(my_ghosting,ghost_entities,{});
  }

  void createParticleOnProc0()
  {
    if (bulk->parallel_rank() == 0)
    {
      particle = bulk->declare_element(1, stk::mesh::ConstPartVector{particle_part});
      node1 = bulk->declare_node(1, stk::mesh::ConstPartVector{node_part});
      bulk->declare_relation(particle,node1,0);
    }
  }

  void createParts()
  {
    node_part = &meta->get_topology_root_part(stk::topology::NODE);
    particle_part = &meta->get_topology_root_part(stk::topology::PARTICLE);
    beam_part = &meta->get_topology_root_part(stk::topology::BEAM_2);
  }

  void createBeamOnProc1()
  {
    if (bulk->parallel_rank() == 1)
      createBeam();
  }

  void declareBeamRelationsToNodes()
  {
    node1 = bulk->get_entity(stk::topology::NODE_RANK,1);
    node2 = bulk->declare_node(2, stk::mesh::ConstPartVector{node_part});
    bulk->declare_relation(beam,node1,0);
    bulk->declare_relation(beam,node2,1);
  }

  void createBeam()
  {
    beam = bulk->declare_element(2, stk::mesh::ConstPartVector{beam_part});
    declareBeamRelationsToNodes();
  }

  void addNode1Sharing()
  {
    addNode1SharingOnProc0();
    addNode1SharingOnProc1();
  }

  void addNode1SharingOnProc1()
  {
    if (bulk->parallel_rank() == 1)
    {
      node1 = bulk->get_entity(stk::topology::NODE_RANK,1);
      bulk->add_node_sharing(node1,0);
    }

  }
  void addNode1SharingOnProc0()
  {
    if (bulk->parallel_rank() == 0)
    {
      bulk->add_node_sharing(node1,1);
    }

  }

  void destroyBeam()
  {
    EXPECT_TRUE(bulk->is_valid(beam));
    EXPECT_TRUE(bulk->destroy_entity(beam));
    EXPECT_FALSE(bulk->is_valid(beam));
  }

  stk::mesh::MetaData* meta = nullptr;
  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::Part * node_part;
  stk::mesh::Part * beam_part;
  stk::mesh::Part * particle_part;
  stk::mesh::Entity beam;
  stk::mesh::Entity particle;
  stk::mesh::Entity node1;
  stk::mesh::Entity node2;

};

TEST_F(ParticleCustomGhostTester, KeepCustomGhostAfterLossOfSharing_ticket_14261)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    initialize(MPI_COMM_WORLD,stk::mesh::BulkData::NO_AUTO_AURA);
    createParticleOnProc0AndGhostToProc1();
    createBeamOnProc1AndUpdateSharingOnProc0();
    destroyBeamOnProc1AndVerifyParticleIsStillOnProc1();
  }
}

