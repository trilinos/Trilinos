/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/fixtures/RingFixture.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <unit_tests/UnitTestModificationEndWrapper.hpp>

using stk::mesh::MetaData;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::Entity;
using stk::mesh::EntityKey;
using stk::mesh::EntityProc;
using stk::mesh::fixtures::RingFixture;

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

STKUNIT_UNIT_TEST( UnitTestBoxFixture, verifyRingFixture )
{
  // A unit test to verify the correctness of the RingFixture fixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create the ring fixture we'll be testing

  RingFixture fixture(pm);
  MetaData& meta = fixture.m_meta_data;
  BulkData& bulk = fixture.m_bulk_data;

  const stk::mesh::EntityRank element_rank = MetaData::ELEMENT_RANK;

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

  std::vector<unsigned> local_count;
  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[NODE_RANK]     , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[element_rank] , nLocalElement );

  std::vector<Entity> all_nodes;
  get_entities( bulk, NODE_RANK, all_nodes);

  unsigned num_selected_nodes =
    count_selected_entities( select_used, bulk.buckets(NODE_RANK) );
  STKUNIT_ASSERT_EQUAL( num_selected_nodes , local_count[NODE_RANK] );

  std::vector<Entity> universal_nodes;
  get_selected_entities(select_all, bulk.buckets(NODE_RANK),
                        universal_nodes );
  STKUNIT_ASSERT_EQUAL( universal_nodes.size() , all_nodes.size() );

  STKUNIT_ASSERT(bulk.modification_end());

  // Verify declarations and sharing two end nodes:

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement );

  if ( 1 < p_size ) {
    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    Entity const node0 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n0] );
    Entity const node1 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n1] );

    STKUNIT_ASSERT( bulk.is_valid(node0) );
    STKUNIT_ASSERT( bulk.is_valid(node1) );

    STKUNIT_ASSERT_EQUAL( bulk.entity_comm_sharing(bulk.entity_key(node0)).size(), 1u );
    STKUNIT_ASSERT_EQUAL( bulk.entity_comm_sharing(bulk.entity_key(node1)).size() , 1u );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( bulk.modification_begin() );
  bulk.change_entity_owner( change );
  STKUNIT_ASSERT( bulk.modification_end());

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement );

  stk::mesh::count_entities( select_all , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode + n_extra );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement + n_extra );

  bulk.modification_begin();
  fixture.fixup_node_ownership();
  STKUNIT_ASSERT(bulk.modification_end());

  // Make sure that element->owner_rank() == element->node[1]->owner_rank()
  if ( 1 < p_size ) {
    stk::mesh::count_entities( select_all , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode + n_extra );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement + n_extra );

    stk::mesh::count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement );

    stk::mesh::count_entities( select_owned , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nPerProc );
    STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nPerProc );
  }
}

}

namespace stk {
namespace unit_test {

void test_shift_ring( RingFixture& ring, bool generate_aura=true )
{
  MetaData& meta = ring.m_meta_data;
  BulkData& bulk = ring.m_bulk_data;

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

  std::vector<unsigned> local_count ;
  std::vector<EntityProc> change ;

  Entity send_element_1 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_send ] );
  Entity send_element_2 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_send + 1 ] );

  Entity send_node_1 = *(bulk.begin_nodes(send_element_1) + 1);
  Entity send_node_2 = *(bulk.begin_nodes(send_element_2) + 1);

  Entity recv_element_1 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_recv ] );
  Entity recv_element_2 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( bulk.is_valid(send_element_1) && p_rank == bulk.parallel_owner_rank(send_element_1) );
  STKUNIT_ASSERT( bulk.is_valid(send_element_2) && p_rank == bulk.parallel_owner_rank(send_element_2) );
  STKUNIT_ASSERT( !bulk.is_valid(recv_element_1) || p_rank != bulk.parallel_owner_rank(recv_element_1) );
  STKUNIT_ASSERT( !bulk.is_valid(recv_element_2) || p_rank != bulk.parallel_owner_rank(recv_element_2) );

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

  STKUNIT_ASSERT( bulk.modification_begin() );
  bulk.change_entity_owner( change );
  STKUNIT_ASSERT( stk::unit_test::modification_end_wrapper( bulk , generate_aura ) );

  send_element_1 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_send ] );
  send_element_2 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_send + 1 ] );
  recv_element_1 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_recv ] );
  recv_element_2 = bulk.get_entity( MetaData::ELEMENT_RANK , ring.m_element_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( !bulk.is_valid(send_element_1) || p_rank != bulk.parallel_owner_rank(send_element_1) );
  STKUNIT_ASSERT( !bulk.is_valid(send_element_2) || p_rank != bulk.parallel_owner_rank(send_element_2) );
  STKUNIT_ASSERT( bulk.is_valid(recv_element_1) && p_rank == bulk.parallel_owner_rank(recv_element_1) );
  STKUNIT_ASSERT( bulk.is_valid(recv_element_2) && p_rank == bulk.parallel_owner_rank(recv_element_2) );

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::NODE_RANK] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[MetaData::ELEMENT_RANK] , nLocalElement );

  unsigned count_shared = 0 ;
  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i = bulk.comm_list().begin() ;
        i != bulk.comm_list().end() ; ++i )
  {
    if ( bulk.in_shared( i->key ) ) { ++count_shared ; }
  }
  STKUNIT_ASSERT_EQUAL( count_shared , 2u );

  {
    const EntityKey node_recv = EntityKey(NODE_RANK , ring.m_node_ids[id_recv]);
    const EntityKey node_send = EntityKey(NODE_RANK , ring.m_node_ids[id_send]);

    STKUNIT_ASSERT_EQUAL( bulk.entity_comm_sharing(node_recv).size(), 1u );
    STKUNIT_ASSERT_EQUAL( bulk.entity_comm_sharing(node_send).size() , 1u );
  }
}

}
}
