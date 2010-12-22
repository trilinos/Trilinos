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
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::Entity;
using stk::mesh::EntityProc;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fem::NODE_RANK;

STKUNIT_UNIT_TEST( UnitTestBoxFixture, verifyRingFixture )
{
  // A unit test to verify the correctness of the RingFixture fixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create the ring fixture we'll be testing

  RingFixture fixture(pm);
  MetaData& meta = fixture.m_meta_data;
  BulkData& bulk = fixture.m_bulk_data;

  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fixture.m_fem);

  meta.commit();

  const unsigned p_rank     = bulk.parallel_rank();
  const unsigned p_size     = bulk.parallel_size();
  const unsigned nPerProc   = fixture.m_num_edge_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;
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
  STKUNIT_ASSERT_EQUAL( local_count[element_rank] , nLocalEdge );

  std::vector<Entity*> all_nodes;
  get_entities( bulk, NODE_RANK, all_nodes);

  unsigned num_selected_nodes =
    count_selected_entities( select_used, bulk.buckets(NODE_RANK) );
  STKUNIT_ASSERT_EQUAL( num_selected_nodes , local_count[NODE_RANK] );

  std::vector<Entity*> universal_nodes;
  get_selected_entities(select_all, bulk.buckets(NODE_RANK),
                        universal_nodes );
  STKUNIT_ASSERT_EQUAL( universal_nodes.size() , all_nodes.size() );

  STKUNIT_ASSERT(bulk.modification_end());

  // Verify declarations and sharing two end nodes:

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge );

  if ( 1 < p_size ) {
    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    Entity * const node0 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n0] );
    Entity * const node1 = bulk.get_entity( NODE_RANK , fixture.m_node_ids[n1] );

    STKUNIT_ASSERT( node0 != NULL );
    STKUNIT_ASSERT( node1 != NULL );

    STKUNIT_ASSERT_EQUAL( node0->sharing().size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( node1->sharing().size() , size_t(1) );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( bulk.modification_begin() );
  bulk.change_entity_owner( change );
  STKUNIT_ASSERT( bulk.modification_end());

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge );

  stk::mesh::count_entities( select_all , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode + n_extra );
  STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge + n_extra );

  bulk.modification_begin();
  fixture.fixup_node_ownership();
  STKUNIT_ASSERT(bulk.modification_end());

  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank()
  if ( 1 < p_size ) {
    stk::mesh::count_entities( select_all , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode + n_extra );
    STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge + n_extra );

    stk::mesh::count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode );
    STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge );

    stk::mesh::count_entities( select_owned , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( local_count[0] , nPerProc );
    STKUNIT_ASSERT_EQUAL( local_count[1] , nPerProc );
  }
}

namespace stk {
namespace unit_test {

void test_shift_ring( RingFixture& ring, bool generate_aura=true )
{
  MetaData& meta = ring.m_meta_data;
  BulkData& bulk = ring.m_bulk_data;

  const unsigned p_rank     = bulk.parallel_rank();
  const unsigned p_size     = bulk.parallel_size();
  const unsigned nPerProc   = ring.m_num_edge_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  const unsigned p_send  = ( p_rank + 1 ) % p_size ;
  const unsigned id_send = id_end - 2 ;
  const unsigned id_recv = ( id_begin + id_total - 2 ) % id_total ;

  Selector select_used = meta.locally_owned_part() |
                         meta.globally_shared_part();

  std::vector<unsigned> local_count ;
  std::vector<EntityProc> change ;

  Entity * send_edge_1 = bulk.get_entity( 1 , ring.m_edge_ids[ id_send ] );
  Entity * send_edge_2 = bulk.get_entity( 1 , ring.m_edge_ids[ id_send + 1 ] );
  Entity * send_node_1 = send_edge_1->relations()[1].entity();
  Entity * send_node_2 = send_edge_2->relations()[1].entity();
  Entity * recv_edge_1 = bulk.get_entity( 1 , ring.m_edge_ids[ id_recv ] );
  Entity * recv_edge_2 = bulk.get_entity( 1 , ring.m_edge_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( NULL != send_edge_1 && p_rank == send_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL != send_edge_2 && p_rank == send_edge_2->owner_rank() );
  STKUNIT_ASSERT( NULL == recv_edge_1 || p_rank != recv_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL == recv_edge_2 || p_rank != recv_edge_2->owner_rank() );

  if ( p_rank == send_node_1->owner_rank() ) {
    EntityProc entry( send_node_1 , p_send );
    change.push_back( entry );
  }
  if ( p_rank == send_node_2->owner_rank() ) {
    EntityProc entry( send_node_2 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_edge_1 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_edge_2 , p_send );
    change.push_back( entry );
  }

  send_edge_1 = NULL ;
  send_edge_2 = NULL ;
  send_node_1 = NULL ;
  send_node_2 = NULL ;
  recv_edge_1 = NULL ;
  recv_edge_2 = NULL ;

  STKUNIT_ASSERT( bulk.modification_begin() );
  bulk.change_entity_owner( change );
  STKUNIT_ASSERT( stk::unit_test::modification_end_wrapper( bulk , generate_aura ) );

  send_edge_1 = bulk.get_entity( 1 , ring.m_edge_ids[ id_send ] );
  send_edge_2 = bulk.get_entity( 1 , ring.m_edge_ids[ id_send + 1 ] );
  recv_edge_1 = bulk.get_entity( 1 , ring.m_edge_ids[ id_recv ] );
  recv_edge_2 = bulk.get_entity( 1 , ring.m_edge_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( NULL == send_edge_1 || p_rank != send_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL == send_edge_2 || p_rank != send_edge_2->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_1 && p_rank == recv_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_2 && p_rank == recv_edge_2->owner_rank() );

  stk::mesh::count_entities( select_used , bulk , local_count );
  STKUNIT_ASSERT_EQUAL( local_count[0] , nLocalNode );
  STKUNIT_ASSERT_EQUAL( local_count[1] , nLocalEdge );

  unsigned count_shared = 0 ;
  for ( std::vector<Entity*>::const_iterator
        i = bulk.entity_comm().begin() ;
        i != bulk.entity_comm().end() ; ++i ) {
    if ( in_shared( **i ) ) { ++count_shared ; }
  }
  STKUNIT_ASSERT_EQUAL( count_shared , 2u );

  {
    Entity * const node_recv = bulk.get_entity( NODE_RANK , ring.m_node_ids[id_recv] );
    Entity * const node_send = bulk.get_entity( NODE_RANK , ring.m_node_ids[id_send] );

    STKUNIT_ASSERT_EQUAL( node_recv->sharing().size() , 1u );
    STKUNIT_ASSERT_EQUAL( node_send->sharing().size() , 1u );
  }
}

}
}
