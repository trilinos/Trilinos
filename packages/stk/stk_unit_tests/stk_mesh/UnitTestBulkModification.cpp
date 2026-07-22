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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include <stddef.h>                                          // for size_t
#include <stdexcept>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>                          // for Entity
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <vector>                                            // for vector, etc
#include "mpi.h"

#include "stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp"
#include "stk_mesh/base/Bucket.hpp"                          // for Bucket, etc
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/Part.hpp"                            // for Part
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_util/util/SortAndUnique.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Bucket;
using stk::mesh::BucketIterator;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::BucketVector;
using stk::mesh::fixtures::RingFixture;

class UnitTestStkMeshBulkModification {
public:
  UnitTestStkMeshBulkModification(stk::ParallelMachine pm)
    : m_comm(pm),
      m_num_procs(stk::parallel_machine_size( m_comm )),
      m_rank(stk::parallel_machine_rank( m_comm )),
      m_ring_mesh(pm)
  { }

  void test_bulkdata_not_synchronized();
  void test_all_local_nodes();
  void test_all_local_elements();
  void test_parallel_consistency();

  BulkData& initialize_ring_fixture()
  {
    m_ring_mesh.m_meta_data.commit();
    BulkData& bulk_data = m_ring_mesh.m_bulk_data;

    bulk_data.modification_begin();
    m_ring_mesh.generate_mesh( );
    STK_ThrowRequire(bulk_data.modification_end());

    m_ring_mesh.fixup_node_ownership( );

    return bulk_data;
  }

  stk::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
  RingFixture m_ring_mesh;
};

namespace {

const EntityRank NODE_RANK = stk::topology::NODE_RANK;

TEST( UnitTestBulkDataNotSynchronized , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_bulkdata_not_synchronized();
}

TEST( UnitTestAllLocalNodes , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_all_local_nodes();
}

TEST( UnitTestAllLocalElements , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_all_local_elements();
}

TEST( UnitTestParallelConsistency , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_parallel_consistency();
}

} //end namespace

void UnitTestStkMeshBulkModification::test_bulkdata_not_synchronized()
{
  BulkData& bulk_data = initialize_ring_fixture();

  bulk_data.modification_begin(); // Intentionally make things unsynced

  std::vector< Entity> entities;
  std::vector< Entity> entities_closure;
  ASSERT_THROW(stk::mesh::find_closure(bulk_data, entities, entities_closure), std::logic_error);
}

void UnitTestStkMeshBulkModification::test_all_local_nodes()
{
  BulkData& bulk_data = initialize_ring_fixture();

  {
    std::vector< Entity> entities;
    std::vector< Entity> entities_closure;
    find_closure(bulk_data, entities, entities_closure);

    // the closure of the an empty set of entities on all procs should be empty
    EXPECT_TRUE(entities_closure.empty());
  }

  {
    // Get a selector for the univeral part (contains local, shared, and ghosted)
    const stk::mesh::Part& universal = m_ring_mesh.m_meta_data.universal_part();
    stk::mesh::Selector universal_selector(universal);

    // Get the buckets that will give us the universal nodes
    BucketVector buckets = bulk_data.get_buckets(NODE_RANK, universal_selector);

    // Get the universal nodes
    std::vector< Entity> universal_entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    stk::util::sort_and_unique(universal_entities, stk::mesh::EntityLess(bulk_data));

    // Get the buckets that will give us the locally used nodes
    stk::mesh::Selector locally_used_selector =
        m_ring_mesh.m_meta_data.locally_owned_part() |
        m_ring_mesh.m_meta_data.globally_shared_part();

    buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, locally_used_selector);

    // Get the locally used nodes
    std::vector< Entity> entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }

    // Get the closure, passing in the locally used nodes on each proc
    std::vector< Entity> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted nodes on this part will be locally used on one of the other
    // procs, so we expect that they will be part of the closure. In other
    // words, the set of nodes returned by find_closure should exactly match
    // the set of universal nodes.
    ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      EXPECT_TRUE(universal_entities[i] == entities_closure[i]);
    }
  }
}

void UnitTestStkMeshBulkModification::test_all_local_elements()
{
  BulkData& bulk_data = initialize_ring_fixture();
  {
    const stk::mesh::Part& universal = m_ring_mesh.m_meta_data.universal_part();
    stk::mesh::Selector universal_selector(universal);

    BucketVector buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, universal_selector);

    // get all the nodes that this process knows about
    std::vector< Entity> universal_entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, universal_selector);

    // get all the elements that this process knows about
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    // universal entities should now have all the universal nodes and elements
    stk::util::sort_and_unique(universal_entities, stk::mesh::EntityLess(bulk_data));

    // get the buckets that we need to traverse to get the locally used elements
    stk::mesh::Selector locally_used_selector =
        m_ring_mesh.m_meta_data.locally_owned_part() |
        m_ring_mesh.m_meta_data.globally_shared_part();

    buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_used_selector);

    // get the locally used elements and store them in entities
    std::vector< Entity> entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }

    // call find_closure, passing in the locally used elements
    std::vector< Entity> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted entities on this proc (element or node) should be contained
    // in the closure of the locally-used element on some other proc, so we
    // expect that they will be part of the closure. In other
    // words, the set of entities returned by find_closure should exactly match
    // the set of universal entities (nodes and elements).
    ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      EXPECT_TRUE(universal_entities[i] == entities_closure[i]);
    }
  }
}

void UnitTestStkMeshBulkModification::test_parallel_consistency()
{
  BulkData& bulk_data = initialize_ring_fixture();

  std::vector< Entity> entities;
  std::vector< Entity> entities_closure;

  // For proc 0 only, add locally used nodes to entities, for all other
  // procs, leave entities empty.
  if (m_rank == 0) {
    stk::mesh::Selector locally_used_selector =
        m_ring_mesh.m_meta_data.locally_owned_part() |
        m_ring_mesh.m_meta_data.globally_shared_part();

    BucketVector const& buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, locally_used_selector);

    for (BucketVector::const_iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }
  }

  // Call find_closure with proc 0 passing in locally-used nodes
  stk::mesh::find_closure(bulk_data, entities, entities_closure);

  // Proc 0 will broadcast the global ids of the nodes it passed to
  // find_closure

  stk::CommBroadcast all(bulk_data.parallel(), 0);

  stk::pack_and_communicate(all, [&](){
    for (Entity entity : entities) {
      all.send_buffer().pack<stk::mesh::EntityKey>(bulk_data.entity_key(entity));
    }
  });

  // clear-out entities and put the nodes that correspond to the keys
  // broadcast by proc 0 into entities.
  entities.clear();
  stk::CommBuffer& buf = all.recv_buffer();
  stk::mesh::EntityKey k ;
  while ( buf.remaining() ) {
    buf.unpack<stk::mesh::EntityKey>(k);
    Entity e = bulk_data.get_entity(k);
    // If a proc is not aware of a key, that means it has no relationship
    // with that entity, so it can ignore it.
    if (bulk_data.is_valid(e)) {
      entities.push_back(e);
    }
  }

  stk::util::sort_and_unique(entities, stk::mesh::EntityLess(bulk_data));

  // If any processor had ghosted nodes that were local to proc 0, those
  // nodes should be in the closure because proc 0 passed them in to
  // find_closure.
  ASSERT_TRUE(entities.size() == entities_closure.size());
  for (size_t i = 0; i < entities_closure.size(); ++i) {
    EXPECT_TRUE(entities[i] == entities_closure[i]);
  }
}
