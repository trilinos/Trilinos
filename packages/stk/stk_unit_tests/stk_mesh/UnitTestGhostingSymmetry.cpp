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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include "UnitTestUtils.hpp"
#include <memory>
#include <vector>
#include <map>

using KeyProcMap = std::map<stk::mesh::EntityKey,std::vector<int>>;

template<class DoThis>
void visit_all_comm_entities(const stk::mesh::BulkData& bulk, DoThis doThis)
{
  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
  std::vector<int> commProcs;
  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    const stk::mesh::BucketVector& bkts = bulk.buckets(rank);

    for(const stk::mesh::Bucket* bptr : bkts) {
      for(stk::mesh::Entity ent : *bptr) {
        bulk.comm_procs(ent, commProcs);

        if (!commProcs.empty()) {
          doThis(bulk.entity_key(ent), commProcs);
        }
      }
    }
  }
}

void fill_key_proc_map(const stk::mesh::BulkData& bulk,
                       KeyProcMap& keyProcMap)
{
  visit_all_comm_entities(bulk,
    [&](stk::mesh::EntityKey key, std::vector<int>& commProcs)
    {
      keyProcMap[key] = commProcs;
    });

  //Now communicate to tell each proc which other procs we know
  //about for each entity. This will make keyProcMap globally symmetric
  //because the owner always knows about all comm-procs even if
  //ghost-receivers don't know about each other.

  stk::CommSparse commSparse(bulk.parallel());
  stk::pack_and_communicate(commSparse,
    [&commSparse, &keyProcMap]()
    {
      for(const auto& keyProcs : keyProcMap) {
        const stk::mesh::EntityKey& key = keyProcs.first;
        const std::vector<int>& procs = keyProcs.second;
        for(int p : procs) {
          commSparse.send_buffer(p).pack<stk::mesh::EntityKey>(key);
          commSparse.send_buffer(p).pack(procs);
        }
      }
    });

  for(int p=0; p<commSparse.parallel_size(); ++p) {
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      stk::mesh::EntityKey key;
      buf.unpack<stk::mesh::EntityKey>(key);
      std::vector<int> recvdProcs;
      buf.unpack(recvdProcs);
      KeyProcMap::iterator iter = keyProcMap.find(key);
      STK_ThrowRequireMsg(iter != keyProcMap.end(),"P"<<bulk.parallel_rank()<<" Failed to find "<<key<<" in keyProcMap");
      stk::util::insert_keep_sorted_and_unique(p, recvdProcs);
      std::vector<int>& procs = iter->second;
      stk::util::insert_keep_sorted_and_unique(recvdProcs, procs);
    }
  }
}

bool is_comm_globally_symmetric(const stk::mesh::BulkData& bulk)
{
  std::map<stk::mesh::EntityKey,std::vector<int>> keyProcMap;

  fill_key_proc_map(bulk, keyProcMap);

  stk::CommSparse commSparse(bulk.parallel());
  stk::pack_and_communicate(commSparse,
    [&commSparse, &bulk]()
    {
      visit_all_comm_entities(bulk,
        [&](stk::mesh::EntityKey key, std::vector<int>& commProcs)
        {
          for(int p : commProcs) {
            commSparse.send_buffer(p).pack<stk::mesh::EntityKey>(key);
          }
        });
    });

  bool returnValue = true;

  for(int p=0; p<commSparse.parallel_size(); ++p) {
    if (returnValue == false) {
      break;
    }

    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      stk::mesh::EntityKey key;
      buf.unpack<stk::mesh::EntityKey>(key);
      KeyProcMap::iterator iter = keyProcMap.find(key);
      if (iter == keyProcMap.end()) {
        returnValue = false;
        break;
      }

      std::vector<int>& procs = iter->second;
      std::vector<int>::iterator p_iter = std::find(procs.begin(), procs.end(), p);
      if (p_iter == procs.end()) {
        returnValue = false;
        break;
      }
      procs.erase(p_iter);
    }
  }

  for(const auto& keyProcs : keyProcMap) {
    const std::vector<int>& procs = keyProcs.second;
    const bool onlyContainsSelf = (procs.size() == 1 && procs[0] == bulk.parallel_rank());
    if (!procs.empty() && !onlyContainsSelf) {
      returnValue = false;
      break;
    }
  }

  return stk::is_true_on_all_procs(bulk.parallel(), returnValue);
}

std::shared_ptr<stk::mesh::BulkData>
setup_hex_mesh(MPI_Comm comm, size_t nx, size_t ny, size_t nz,
               bool symmetricGhostInfo = false,
               stk::mesh::BulkData::AutomaticAuraOption /*auraOption*/ =  stk::mesh::BulkData::NO_AUTO_AURA)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                         .set_spatial_dimension(3)
                                         .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                         .set_symmetric_ghost_info(symmetricGhostInfo)
                                         .create();
  stk::mesh::fixtures::HexFixture::fill_mesh(nx,ny,nz, *bulkPtr);
  return bulkPtr;
}

TEST(GhostingSymmetry, ghostNode_fromProc0_toProc1andProc2_notSymmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3);

  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1);

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node1, commProcs);

  const bool isOwned = bulkPtr->bucket(node1).owned();
  size_t numOtherProcs = isOwned ? 2 : 1;
  //because ghost receivers know about owner, not each other

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_FALSE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, ghostNode_fromProc0_toProc1andProc2_symmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo);

  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1);

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node1, commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, ghostNode_fromProc0_toProc1andProc2_symmetric_ghosting)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo);

  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1);

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));

  const stk::mesh::Ghosting* myGhosting = bulkPtr->ghostings().back();
  EXPECT_EQ("myCustomGhosting", myGhosting->name());

  std::vector<int> commProcs;
  bulkPtr->comm_procs(*myGhosting, bulkPtr->entity_key(node1), commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, ghostSharedNode_fromProc0_toProc2_notSymmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3);

  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr);

  stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node5, commProcs);

  const bool isOwned = bulkPtr->bucket(node5).owned();
  size_t numOtherProcs = isOwned ? 2 : 1;
  //because ghost receivers know about owner, not each other
  //also proc1 shares node 5 but doesn't know it's ghosted to proc 2

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_FALSE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, ghostSharedNode_fromProc0_toProc2_symmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo);

  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr);

  stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node5, commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, ghostSharedNode_fromProc0_toProc2_then_change_parts_symmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo);

  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr);

  stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node5, commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));

  stk::mesh::Part& nodePart = bulkPtr->mesh_meta_data().declare_part("myNodePart", stk::topology::NODE_RANK);
  stk::mesh::PartVector addParts = {&nodePart};
  stk::mesh::PartVector rmParts;
  stk::mesh::EntityVector entities;
  if (bulkPtr->parallel_rank() == 0) {
    entities.push_back(bulkPtr->get_entity(stk::topology::NODE_RANK, 5));
  }

  bulkPtr->batch_change_entity_parts(entities, addParts, rmParts);

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, two_separate_mesh_mods_symmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo);

  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr, "node5Ghosting");

  stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node5, commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));

  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1, "node1Ghosting");

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));

  bulkPtr->comm_procs(node1, commProcs);

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

TEST(GhostingSymmetry, two_separate_mesh_mods_with_aura_symmetric)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }

  constexpr bool symmetricGhostInfo = true;
  constexpr stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3, symmetricGhostInfo, auraOption);

  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr, "node5Ghosting");

  stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node5, commProcs);

  const size_t numOtherProcs = 2;

  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));

  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1, "node1Ghosting");
  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 5, "node1Ghosting");

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  node5 = bulkPtr->get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulkPtr->is_valid(node1));
  EXPECT_TRUE(bulkPtr->is_valid(node5));

  bulkPtr->comm_procs(node1, commProcs);
  EXPECT_EQ(numOtherProcs, commProcs.size());
  bulkPtr->comm_procs(node5, commProcs);
  EXPECT_EQ(numOtherProcs, commProcs.size());

  EXPECT_TRUE(is_comm_globally_symmetric(*bulkPtr));
}

