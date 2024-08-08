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

#include "mpi.h"                        // for MPI_Barrier, MPI_COMM_WORLD, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityCommListInfo.hpp"
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp"  // for RingFixture
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <vector>                       // for vector, etc

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::Entity;
using stk::mesh::EntityKey;
using stk::mesh::EntityProc;
using stk::mesh::fixtures::RingFixture;

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

TEST( UnitTestBoxFixture, verifyRingFixture )
{
  // A unit test to verify the correctness of the RingFixture fixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );
  // Create the ring fixture we'll be testing

  RingFixture fixture(pm);
  MetaData& meta = fixture.m_meta_data;
  BulkData& bulk = fixture.m_bulk_data;

  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

  meta.commit();

  const unsigned p_rank     = bulk.parallel_rank();
  const unsigned p_size     = bulk.parallel_size();
  const unsigned nPerProc   = fixture.m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;
  const unsigned n_extra    = 1 < p_size ? 2 : 0 ;

  bulk.modification_begin();
  fixture.generate_mesh();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used = meta.locally_owned_part() |
      meta.globally_shared_part();
  Selector select_all(  meta.universal_part() );

  std::vector<size_t> local_count;
  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( local_count[NODE_RANK]     , nLocalNode );
  ASSERT_EQ( local_count[element_rank] , nLocalElement );

  std::vector<Entity> all_nodes;
  get_entities( bulk, NODE_RANK, all_nodes);

  unsigned num_selected_nodes =
      count_selected_entities( select_used, bulk.buckets(NODE_RANK) );
  ASSERT_EQ( num_selected_nodes , local_count[NODE_RANK] );

  std::vector<Entity> universal_nodes;
  get_selected_entities(select_all, bulk.buckets(NODE_RANK),
                        universal_nodes );
  ASSERT_EQ( universal_nodes.size() , all_nodes.size() );

  ASSERT_TRUE(bulk.modification_end());

  // Verify declarations and sharing two end nodes:

  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
  ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

  if ( 1 < p_size ) {
    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    Entity const node0 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n0] );
    Entity const node1 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n1] );

    ASSERT_TRUE( bulk.is_valid(node0) );
    ASSERT_TRUE( bulk.is_valid(node1) );

    std::vector<int> node0_shared_procs;
    bulk.comm_shared_procs(bulk.entity_key(node0), node0_shared_procs);
    std::vector<int> node1_shared_procs;
    bulk.comm_shared_procs(bulk.entity_key(node0), node1_shared_procs);
    ASSERT_EQ( node0_shared_procs.size(), 1u );
    ASSERT_EQ( node1_shared_procs.size(), 1u );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;
  bulk.change_entity_owner( change );

  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
  ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

  stk::mesh::count_entities( select_all , bulk , local_count );
  ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode + n_extra );
  ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement + n_extra );

  fixture.fixup_node_ownership();

  // Make sure that element->owner_rank() == element->node[1]->owner_rank()
  if ( 1 < p_size ) {
    stk::mesh::count_entities( select_all , bulk , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode + n_extra );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement + n_extra );

    stk::mesh::count_entities( select_used , bulk , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

    stk::mesh::count_entities( select_owned , bulk , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nPerProc );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nPerProc );
  }
}

}

namespace stk {
namespace unit_test {

void test_shift_ring( RingFixture& ring )
{
  MetaData& meta = ring.m_meta_data;
  stk::unit_test_util::BulkDataTester& bulk = ring.m_bulk_data;

  const int p_rank     = bulk.parallel_rank();
  const int p_size     = bulk.parallel_size();
  const unsigned nPerProc   = ring.m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;

  const unsigned p_send  = ( p_rank + 1 ) % p_size ;
  const unsigned id_send = id_end - 2 ;
  const unsigned id_recv = ( id_begin + id_total - 2 ) % id_total ;

  Selector select_used = meta.locally_owned_part() |
      meta.globally_shared_part();

  std::vector<size_t> local_count ;
  std::vector<EntityProc> change ;

  Entity send_element_1 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_send ] );
  Entity send_element_2 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_send + 1 ] );

  Entity send_node_1 = *(bulk.begin_nodes(send_element_1) + 1);
  Entity send_node_2 = *(bulk.begin_nodes(send_element_2) + 1);

  Entity recv_element_1 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_recv ] );
  Entity recv_element_2 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_recv + 1 ] );

  ASSERT_TRUE( bulk.is_valid(send_element_1) && p_rank == bulk.parallel_owner_rank(send_element_1) );
  ASSERT_TRUE( bulk.is_valid(send_element_2) && p_rank == bulk.parallel_owner_rank(send_element_2) );
  ASSERT_TRUE( !bulk.is_valid(recv_element_1) || p_rank != bulk.parallel_owner_rank(recv_element_1) );
  ASSERT_TRUE( !bulk.is_valid(recv_element_2) || p_rank != bulk.parallel_owner_rank(recv_element_2) );

  if ( p_rank == bulk.parallel_owner_rank(send_node_1) ) {
    EntityProc entry( send_node_1 , p_send );
    change.push_back( entry );
  }
  if ( p_rank == bulk.parallel_owner_rank(send_node_2) ) {
    EntityProc entry( send_node_2 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_element_1 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_element_2 , p_send );
    change.push_back( entry );
  }

  send_element_1 = Entity();
  send_element_2 = Entity();
  send_node_1 = Entity();
  send_node_2 = Entity();
  recv_element_1 = Entity();
  recv_element_2 = Entity();

  bulk.change_entity_owner( change );

  send_element_1 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_send ] );
  send_element_2 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_send + 1 ] );
  recv_element_1 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_recv ] );
  recv_element_2 = bulk.get_entity( stk::topology::ELEMENT_RANK , ring.m_element_ids[ id_recv + 1 ] );

  ASSERT_TRUE( !bulk.is_valid(send_element_1) || p_rank != bulk.parallel_owner_rank(send_element_1) );
  ASSERT_TRUE( !bulk.is_valid(send_element_2) || p_rank != bulk.parallel_owner_rank(send_element_2) );
  ASSERT_TRUE( bulk.is_valid(recv_element_1) && p_rank == bulk.parallel_owner_rank(recv_element_1) );
  ASSERT_TRUE( bulk.is_valid(recv_element_2) && p_rank == bulk.parallel_owner_rank(recv_element_2) );

  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
  ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

  unsigned count_shared = 0 ;
  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i = bulk.my_internal_comm_list().begin() ;
        i != bulk.my_internal_comm_list().end() ; ++i )
  {
    if ( bulk.in_shared( i->key ) ) { ++count_shared ; }
  }
  ASSERT_EQ( count_shared , 2u );

  {
    const EntityKey node_recv = EntityKey(NODE_RANK , ring.m_node_ids[id_recv]);
    const EntityKey node_send = EntityKey(NODE_RANK , ring.m_node_ids[id_send]);

    std::vector<int> node_recv_shared_procs;
    bulk.comm_shared_procs(node_recv,node_recv_shared_procs);
    ASSERT_EQ( node_recv_shared_procs.size(), 1u );
    std::vector<int> node_send_shared_procs;
    bulk.comm_shared_procs(node_send,node_send_shared_procs);
    ASSERT_EQ( node_send_shared_procs.size(), 1u );
  }
}

}
}
