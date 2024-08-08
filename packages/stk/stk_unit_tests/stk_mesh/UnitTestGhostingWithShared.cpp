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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <iostream>                     // for basic_ostream::operator<<
#include <set>                          // for set, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <string>                       // for string
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, operator<<
#include "stk_mesh/base/EntityLess.hpp"  // for EntityLess
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityProc, BucketVector, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/ioUtils.hpp"  // for fill_mesh_using_stk_io
#include "stk_unit_test_utils/BuildMesh.hpp"

using stk::unit_test_util::build_mesh;

TEST(UnitTestGhosting, ThreeElemSendElemWithNonOwnedNodes)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 3) {
    return;
  }

  int procId = stk::parallel_machine_rank(communicator);

  unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::unit_test_util::BulkDataTester bulk(meta, communicator);
  const std::string generatedMeshSpecification = "generated:1x1x3";
  stk::io::fill_mesh(generatedMeshSpecification, bulk);
  bulk.modification_begin();
  stk::mesh::Ghosting& custom_shared_ghosting = bulk.create_ghosting("custom_shared");
  stk::mesh::Ghosting& custom_shared_ghosting2= bulk.create_ghosting("custom_shared");
  EXPECT_EQ(custom_shared_ghosting.ordinal(), custom_shared_ghosting2.ordinal());
  bulk.modification_end();

  stk::mesh::EntityProcVec ownedEntitiesToGhost;

  if (procId == 1)
  {
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
    int destProc = 2;
    ownedEntitiesToGhost.push_back(stk::mesh::EntityProc(elem2, destProc));
  }


  stk::mesh::EntityProcVec entitiesWithClosure;
  bulk.my_add_closure_entities(custom_shared_ghosting, ownedEntitiesToGhost, entitiesWithClosure);
  stk::util::sort_and_unique(entitiesWithClosure, stk::mesh::EntityLess(bulk));

  if (procId == 1)
  {
    std::vector<stk::mesh::EntityKey> gold_keys;
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 7));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 8));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 9));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 10));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 11));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 12));

    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2));

    ASSERT_EQ(gold_keys.size(), entitiesWithClosure.size());

    int otherProc = 2;

    unsigned i=0;
    for(const stk::mesh::EntityProc& entProc : entitiesWithClosure) {
      EXPECT_EQ(gold_keys[i], bulk.entity_key(entProc.first));
      EXPECT_EQ(otherProc, entProc.second);
      ++i;
    }
  }
  else
  {
    ASSERT_TRUE(entitiesWithClosure.empty());
  }

  stk::mesh::impl::move_unowned_entities_for_owner_to_ghost(bulk, entitiesWithClosure);
  stk::util::sort_and_unique(entitiesWithClosure, stk::mesh::EntityLess(bulk));

  if (procId==0)
  {
    std::vector<stk::mesh::EntityKey> gold_keys;
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 7));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 8));

    ASSERT_EQ(gold_keys.size(), entitiesWithClosure.size());

    int otherProc = 2;

    unsigned i=0;
    for(const stk::mesh::EntityProc& entProc : entitiesWithClosure) {
      EXPECT_EQ(gold_keys[i], bulk.entity_key(entProc.first));
      EXPECT_EQ(otherProc, entProc.second);
      ++i;
    }
  }
  else if (procId==1)
  {
    std::vector<stk::mesh::EntityKey> gold_keys;
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 9));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 10));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 11));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 12));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2));

    ASSERT_EQ(gold_keys.size(), entitiesWithClosure.size());

    int otherProc = 2;

    unsigned i=0;
    for(const stk::mesh::EntityProc& entProc : entitiesWithClosure) {
      EXPECT_EQ(gold_keys[i], bulk.entity_key(entProc.first));
      EXPECT_EQ(otherProc, entProc.second);
      ++i;
    }
  }
  else
  {
    ASSERT_TRUE(entitiesWithClosure.empty());
  }

  bulk.modification_begin();

  stk::mesh::Ghosting &ghosting = bulk.create_ghosting("custom ghost unit test");
  bulk.my_ghost_entities_and_fields(ghosting, std::move(entitiesWithClosure));
  bulk.my_internal_modification_end_for_change_ghosting();

  if (procId == 0)
  {
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 5), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 6), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 7), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 8), 2));
    EXPECT_FALSE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2), 2));
  }
  else if (procId == 1)
  {
    EXPECT_FALSE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 5), 2));
    EXPECT_FALSE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 6), 2));
    EXPECT_FALSE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 7), 2));
    EXPECT_FALSE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 8), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 9), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 10), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 11), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 2));
    EXPECT_TRUE(bulk.my_in_send_ghost(ghosting, stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2), 2));
  }
  else if (procId == 2)
  {
    std::vector<stk::mesh::EntityKey> gold_keys;
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 7));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 8));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 9));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 10));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 11));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 12));
    gold_keys.push_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2));

    for (size_t i=0;i<gold_keys.size();++i)
    {
      stk::mesh::Entity entity = bulk.get_entity(gold_keys[i]);
      ASSERT_TRUE(bulk.is_valid(entity));
      ASSERT_TRUE(bulk.in_receive_ghost(ghosting, gold_keys[i]));
    }
  }
}

TEST(UnitTestGhosting, WithSharedFiltered)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x6";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();

  const stk::mesh::BucketVector& shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, shared_selector);

  size_t num_shared_nodes = 0;
  for(size_t i=0; i<shared_buckets.size(); ++i) {
    num_shared_nodes += shared_buckets[i]->size();
  }

  const size_t expected_num_shared_nodes = 4;
  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);

  stkMeshBulkData.modification_begin();

  //we will attempt to ghost shared-and-owned entities from proc 0 to proc 1.
  //they won't be allowed to be in the ghosting because they are already
  //shared on proc 1.
  int myProc = stkMeshBulkData.parallel_rank();
  int otherProc = 1;
  if (myProc != 0) otherProc = -1;//don't do ghosting if we are proc 1

  stk::mesh::Ghosting& custom_shared_ghosting = stkMeshBulkData.create_ghosting("custom_shared");

  std::vector<stk::mesh::EntityProc> send_shared;
  if (otherProc != -1) {
    for(size_t i=0; i<shared_buckets.size(); ++i) {
      const stk::mesh::Bucket& bucket = *shared_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
        if (bucket.parallel_owner_rank(j) == myProc) {
          send_shared.push_back(std::make_pair(bucket[j], otherProc));
        }
      }
    }
  }

  stkMeshBulkData.change_ghosting(custom_shared_ghosting, send_shared);

  stkMeshBulkData.modification_end();

  std::vector<stk::mesh::EntityProc> custom_shared_send_list;
  custom_shared_ghosting.send_list(custom_shared_send_list);

  size_t expected_send_list_size = 0;
  EXPECT_EQ(expected_send_list_size, custom_shared_send_list.size());
}

TEST(UnitTestGhosting, WithShared)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 3) {
    GTEST_SKIP();
  }

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x6";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();

  const stk::mesh::BucketVector& shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, shared_selector);

  size_t num_shared_nodes = 0;
  for(size_t i=0; i<shared_buckets.size(); ++i) {
    num_shared_nodes += shared_buckets[i]->size();
  }

  int myProc = stkMeshBulkData.parallel_rank();

  size_t expected_num_shared_nodes = 4;
  if (myProc == 1) expected_num_shared_nodes = 8;

  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);

  stkMeshBulkData.modification_begin();

  //we will ghost shared-and-owned entities from proc 0 to proc 2.
  int otherProc = 2;
  if (myProc != 0) otherProc = -1;//don't do ghosting

  stk::mesh::Ghosting& custom_shared_ghosting = stkMeshBulkData.create_ghosting("custom_shared");

  std::vector<stk::mesh::EntityProc> send_shared;
  if (otherProc != -1) {
    for(size_t i=0; i<shared_buckets.size(); ++i) {
      const stk::mesh::Bucket& bucket = *shared_buckets[i];
      for(size_t j=0; j<bucket.size(); ++j) {
        if (bucket.parallel_owner_rank(j) == myProc) {
          send_shared.push_back(std::make_pair(bucket[j], otherProc));
        }
      }
    }
  }

  stkMeshBulkData.change_ghosting(custom_shared_ghosting, send_shared);

  if (myProc == otherProc) {
    stk::mesh::Entity node9 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK,9);
    stk::mesh::Entity node10 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK,10);
    stk::mesh::Entity node11 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK,11);
    stk::mesh::Entity node12 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK,12);
    EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(node9));
    EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(node10));
    EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(node11));
    EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(node12));
  }

  stkMeshBulkData.modification_end();

  std::vector<stk::mesh::EntityProc> custom_shared_send_list;
  custom_shared_ghosting.send_list(custom_shared_send_list);

  EXPECT_EQ(send_shared.size(), custom_shared_send_list.size());
  if (myProc == 0) {
    EXPECT_EQ(num_shared_nodes, custom_shared_send_list.size());
  }
}

TEST(UnitTestGhosting, basic1Elem)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::BulkData& bulk = *bulkPtr;
  const std::string generatedMeshSpecification = "generated:1x1x2";
  stk::io::fill_mesh(generatedMeshSpecification, bulk);

  std::vector<stk::mesh::EntityProc> send;
  if (bulk.parallel_rank() == 0) {
    send.emplace_back(bulk.get_entity(stk::topology::ELEM_RANK, 1), 1);
  }
  bulk.modification_begin();
  stk::mesh::Ghosting& ghosting = bulk.create_ghosting("custom1");
  bulk.change_ghosting(ghosting, send, std::vector<stk::mesh::EntityKey>());
  bulk.modification_end();
}

void create_mesh_with_1_tri_per_proc(stk::mesh::BulkData& bulk)
{
  std::vector<stk::mesh::EntityId> nodeIds[] = { {1, 3, 2}, {2, 3, 4} };
  stk::mesh::EntityId elemIds[] = {1, 2};

  stk::mesh::Part& tri2dPart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::TRI_3_2D);
  stk::mesh::PartVector elemParts = {&tri2dPart};
  const int myProc = bulk.parallel_rank();
  const int otherProc = 1 - myProc;

  bulk.modification_begin();
  stk::mesh::declare_element(bulk, elemParts, elemIds[myProc], nodeIds[myProc]);
  bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 2), otherProc);
  bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 3), otherProc);
  bulk.modification_end();
}

void delete_node_2_and_connect_node_5_to_elem_1(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();
{
std::ostringstream os;
os<<"P"<<bulk.parallel_rank()<<"  **** doing test mod ****" << std::endl;
std::cerr<<os.str();
stk::parallel_machine_barrier(bulk.parallel());
}
  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity node5 = bulk.declare_node(5);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);

    bulk.destroy_relation(elem1, node2, 2);
    bulk.destroy_relation(elem2, node2, 0);
    bulk.destroy_entity(elem2);
    bulk.destroy_entity(node2);
    bulk.declare_relation(elem1, node5, 2);
  }
  bulk.modification_end();
}

TEST(UnitTestGhosting, sharedBecomesAuraGhost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;

  /*Create this mesh:

 P0       P1
    2* *2
    /| |\
   / | | \
 1*  | |  *4
   \ | | /
    \| |/
    3* *3

  nodes 2 and 3 are shared, aura ghosts node 1 to P1 and node 4 to P0
*/

  create_mesh_with_1_tri_per_proc(bulk);

  stk::mesh::EntityVector sharedNodes;
  stk::mesh::get_selected_entities(meta.globally_shared_part(),
                                   bulk.buckets(stk::topology::NODE_RANK), sharedNodes);
  EXPECT_EQ(2u, sharedNodes.size());

  stk::mesh::Selector ghosts = !(meta.locally_owned_part() | meta.globally_shared_part());
  stk::mesh::EntityVector ghostNodes;
  stk::mesh::get_selected_entities(ghosts, bulk.buckets(stk::topology::NODE_RANK), ghostNodes);
  EXPECT_EQ(1u, ghostNodes.size());

  /*Now change the mesh to be this:

 P0       P1
    5* *2
    /| |\
   / | | \
 1*  | |  *4
   \ | | /
    \| |/
    3* *3
  i.e., on P0: disconnect node 2 from elem 1 and elem 2, destroy node 2
  node 3 still shared, aura ghosts nodes 1 and 5 to P1 and nodes 2 and 4 to P0
*/

  delete_node_2_and_connect_node_5_to_elem_1(bulk);

  stk::mesh::get_selected_entities(meta.globally_shared_part(),
                                   bulk.buckets(stk::topology::NODE_RANK), sharedNodes);
  EXPECT_EQ(1u, sharedNodes.size());

  stk::mesh::get_selected_entities(ghosts, bulk.buckets(stk::topology::NODE_RANK), ghostNodes);
  EXPECT_EQ(2u, ghostNodes.size());
}

