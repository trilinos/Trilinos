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
#include "UnitTestCustomAura.hpp"
#include "stk_tools/mesh_tools/CustomAura.hpp"
#include "stk_util/parallel/CommSparse.hpp"

namespace aura_unit_tests {

std::map<int,std::vector<int>> local_node_comm( stk::mesh::BulkData& meshBulk ) {
  std::map<int, std::vector<int>> localComm;
  stk::mesh::EntityVector nodes;
  meshBulk.get_entities(stk::topology::NODE_RANK,meshBulk.mesh_meta_data().locally_owned_part(),nodes);
  for( auto& node : nodes ) {
    int nodeID = meshBulk.identifier(node);
    std::set<int> pIDs;
    pIDs.emplace(meshBulk.parallel_rank());
    std::vector<int> procList;
    meshBulk.comm_procs(node, procList);
    for( auto p : procList ) {
      pIDs.emplace(p);
    }
    std::vector<int> v(pIDs.begin(),pIDs.end());
    localComm[nodeID] = v;
  }
  return localComm;
}

void FourQuadShellsInSequenceFixture::print_local_node_comm(const int rank) {
  if ( rank == this->get_parallel_rank() ) {
    std::map<int, std::vector<int>> nodeComm = local_node_comm(this->get_bulk());
    std::stringstream pstream;
    pstream << "On Proc " <<  this->get_parallel_rank() << ";" << std::endl;
    for (auto& a : nodeComm) {
      pstream << "goldMap[" << a.first << "] = {";
      for (int p : a.second) {
        pstream << " " << p << ",";
      }
      pstream << "};" << std::endl;
    }
    std::cout << pstream.str();
  }
}

std::map<int,std::vector<int>> pre_aura_node_comm(int numRanks) {
  std::map<int, std::vector<int>> goldMap;
  if ( numRanks == 1 ) {
    goldMap[1] = { 0 };
    goldMap[2] = { 0 };
    goldMap[3] = { 0 };
    goldMap[4] = { 0 };
    goldMap[5] = { 0 };
    goldMap[6] = { 0 };
    goldMap[7] = { 0 };
    goldMap[8] = { 0 };
    goldMap[9] = { 0 };
    goldMap[10] = { 0 };
  } else if ( numRanks == 2 ) {
    goldMap[1] = { 0 };
    goldMap[2] = { 0 };
    goldMap[3] = { 0 };
    goldMap[4] = { 0, 1 };
    goldMap[5] = { 1 };
    goldMap[6] = { 0 };
    goldMap[7] = { 0 };
    goldMap[8] = { 0 };
    goldMap[9] = { 0, 1 };
    goldMap[10] = { 1 };
  } else if ( numRanks == 3 ) {
    goldMap[1] = { 0 };
    goldMap[2] = { 0 };
    goldMap[3] = { 0, 1 };
    goldMap[4] = { 1, 2 };
    goldMap[5] = { 2 };
    goldMap[6] = { 0 };
    goldMap[7] = { 0 };
    goldMap[8] = { 0, 1 };
    goldMap[9] = { 1, 2 };
    goldMap[10] = { 2 };
  } else if ( numRanks == 4) {
    goldMap[1] = { 0 };
    goldMap[2] = { 0, 1 };
    goldMap[3] = { 1, 2 };
    goldMap[4] = { 2, 3 };
    goldMap[5] = { 3 };
    goldMap[6] = { 0 };
    goldMap[7] = { 0, 1 };
    goldMap[8] = { 1, 2 };
    goldMap[9] = { 2, 3 };
    goldMap[10] = { 3 };
  }
  return goldMap;
}

std::map<int,std::vector<int>> post_aura_node_comm(int numRanks, int numRings) {
  std::map<int, std::vector<int>> goldMap;
  if ( numRings == 1 ) {
    if ( numRanks == 1 ) {
      goldMap[1] = { 0 };
      goldMap[2] = { 0 };
      goldMap[3] = { 0 };
      goldMap[4] = { 0 };
      goldMap[5] = { 0 };
      goldMap[6] = { 0 };
      goldMap[7] = { 0 };
      goldMap[8] = { 0 };
      goldMap[9] = { 0 };
      goldMap[10] = { 0 };
    } else if ( numRanks == 2 ) {
      goldMap[1] = { 0 };
      goldMap[2] = { 0, 1 };
      goldMap[3] = { 0, 1 };
      goldMap[4] = { 0, 1 };
      goldMap[5] = { 0, 1 };
      goldMap[6] = { 0 };
      goldMap[7] = { 0, 1 };
      goldMap[8] = { 0, 1 };
      goldMap[9] = { 0, 1 };
      goldMap[10] = { 0, 1 };
    } else if ( numRanks == 3 ) {
      goldMap[1] = { 0, 1 };
      goldMap[2] = { 0, 1 };
      goldMap[3] = { 0, 1, 2 };
      goldMap[4] = { 0, 1, 2 };
      goldMap[5] = { 1, 2 };
      goldMap[6] = { 0, 1 };
      goldMap[7] = { 0, 1 };
      goldMap[8] = { 0, 1, 2 };
      goldMap[9] = { 0, 1, 2 };
      goldMap[10] = { 1, 2 };
    } else if ( numRanks == 4) {
      goldMap[1] = { 0, 1 };
      goldMap[2] = { 0, 1, 2 };
      goldMap[3] = { 0, 1, 2, 3 };
      goldMap[4] = { 1, 2, 3 };
      goldMap[5] = { 2, 3 };
      goldMap[6] = { 0, 1 };
      goldMap[7] = { 0, 1, 2 };
      goldMap[8] = { 0, 1, 2, 3 };
      goldMap[9] = { 1, 2, 3 };
      goldMap[10] = { 2, 3 };
    }
  } else {
    for ( int nID(1); nID <= 10; ++nID ) {
      std::vector<int> pIDs;
      for ( int p(0); p < numRanks; ++p ) {
        pIDs.push_back(p);
      }
      goldMap[nID] = pIDs;
    }
  }
  return goldMap;
}

void test_node_comm_pre_aura( stk::mesh::BulkData& bulk ) {
  std::map<int,std::vector<int>> nodeComm = local_node_comm(bulk);
  std::map<int,std::vector<int>> goldComm = pre_aura_node_comm(bulk.parallel_size());
  for( auto& a : nodeComm ) {
    EXPECT_EQ(a.second,goldComm[a.first]) << " pre aura for node " << a.first;
  }
}

void test_node_comm_post_aura( stk::mesh::BulkData& bulk, int numRings ) {
  std::map<int,std::vector<int>> nodeComm = local_node_comm(bulk);
  std::map<int,std::vector<int>> goldComm = post_aura_node_comm(bulk.parallel_size(),numRings);
  for( auto& a : nodeComm ) {
    EXPECT_EQ(a.second,goldComm[a.first]) << " post aura for node " << a.first;
  }
}

void sym_comm_pre_ghost_check(stk::mesh::BulkData& bulk) {
  stk::mesh::Entity n1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  if ( bulk.parallel_rank() == 0 ) {
    EXPECT_TRUE(n1 != stk::mesh::Entity::InvalidEntity);
  } else {
    EXPECT_TRUE(n1 == stk::mesh::Entity::InvalidEntity);
  }
}

void sym_comm_post_ghost_check(stk::mesh::BulkData& bulk) {
  stk::mesh::Entity n1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(n1 != stk::mesh::Entity::InvalidEntity);

  if ( bulk.parallel_size() == 1 ) {
    return;
  }

  std::vector<int> goldGhostCommProcs;
  for ( int p=0; p < bulk.parallel_size(); ++p ) {
    if (p != bulk.parallel_rank()) {
      goldGhostCommProcs.push_back(p);
    }
  }

  std::vector<int> commProcs;
  bulk.comm_procs(n1,commProcs);
  std::sort(commProcs.begin(),commProcs.end());

  EXPECT_EQ(commProcs,goldGhostCommProcs);

  std::map<stk::mesh::EntityKey, std::set<int>> symmMap;
  stk::tools::populate_symm_comm_map(bulk,symmMap,stk::topology::NODE_RANK);

  // Desired comm_procs
  std::set<int> goldDesiredCommProcs;
  for ( int p(0); p < bulk.parallel_size(); ++p ) {
    goldDesiredCommProcs.insert(p);
  }
  std::set<int>& n1CommProcs = symmMap[bulk.entity_key(n1)];
  EXPECT_EQ(n1CommProcs,goldDesiredCommProcs);

}

TEST_F(FourQuadShellsInSequenceFixture, GenSymmComm) {

  stk::mesh::BulkData& bulk = this->get_bulk();
  sym_comm_pre_ghost_check(bulk);

  // Ghost Node1 to all Procs
  stk::mesh::Ghosting * customGhost = nullptr;
  if (bulk.parallel_size() > 1) {
    bulk.modification_begin();
    customGhost = &bulk.create_ghosting("UnitGenSymmComm");
    stk::mesh::EntityProcVec entitiesToGhost;
    if (bulk.parallel_rank() == 0) {
      stk::mesh::Entity n1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
      for ( int p(1); p < bulk.parallel_size(); ++p ) {
        entitiesToGhost.push_back(stk::mesh::EntityProc(n1,p));
      }
    }
    bulk.change_ghosting(*customGhost, entitiesToGhost);
    bulk.modification_end();
  }

  sym_comm_post_ghost_check(bulk);

  stk::tools::destroy_custom_aura(bulk, customGhost);
}

TEST_F(FourQuadShellsInSequenceFixture, OneRingTest) {

  test_node_comm_pre_aura(this->get_bulk());
  //    std::cout << "Pre Aura:" << std::endl;
  //    for ( int i(0); i < this->get_parallel_size(); ++i ) {
  //      print_local_node_comm(i);
  //    }

  const int num_rings = 1;
  stk::mesh::Ghosting* g = stk::tools::create_ringed_aura(this->get_bulk(), this->get_meta().universal_part(), "n_ring_aura", num_rings);
  if ( this->get_parallel_size() == 1 ) {
    EXPECT_EQ(g,nullptr);
  }

  test_node_comm_post_aura(this->get_bulk(),num_rings);
  //    std::cout << "Post Aura:" << std::endl;
  //    for ( int i(0); i < this->get_parallel_size(); ++i ) {
  //      print_local_node_comm(i);
  //    }

  EXPECT_EQ(1, num_rings);
  stk::tools::destroy_custom_aura(this->get_bulk(), g);
}

TEST_F(FourQuadShellsInSequenceFixture, TwoRingTest) {
  test_node_comm_pre_aura(this->get_bulk());
  const int num_rings = 2;
  stk::mesh::Ghosting* g = stk::tools::create_ringed_aura(this->get_bulk(), this->get_meta().universal_part(), "n_ring_aura", num_rings);
  if ( this->get_parallel_size() == 1 ) {
    EXPECT_EQ(g,nullptr);
  }
  test_node_comm_post_aura(this->get_bulk(),num_rings);
  EXPECT_EQ(2, num_rings);
  stk::tools::destroy_custom_aura(this->get_bulk(), g);
}

TEST_F(FourQuadShellsInSequenceFixture, ThreeRingTest) {
  test_node_comm_pre_aura(this->get_bulk());
  const int num_rings = 3;
  stk::mesh::Ghosting* g = stk::tools::create_ringed_aura(this->get_bulk(), this->get_meta().universal_part(), "n_ring_aura", num_rings);
  if ( this->get_parallel_size() == 1 ) {
    EXPECT_EQ(g,nullptr);
  }
  test_node_comm_post_aura(this->get_bulk(),num_rings);
  EXPECT_EQ(3, num_rings);
  stk::tools::destroy_custom_aura(this->get_bulk(), g);
}

}

