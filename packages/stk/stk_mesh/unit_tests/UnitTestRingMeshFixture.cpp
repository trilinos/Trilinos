/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <unit_tests/UnitTestRingMeshFixture.hpp>

#include <unit_tests/stk_utest_macros.hpp>
#include <unit_tests/UnitTestBulkData.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

using namespace stk;
using namespace stk::mesh;

RingMeshFixture::RingMeshFixture(
  ParallelMachine pm ,
  unsigned        num_edge_per_proc ,
  bool            use_edge_parts )
  : m_meta_data( fem_entity_rank_names() ),
    m_bulk_data( m_meta_data, pm, 100 ),
    m_edge_parts(),
    m_edge_part_extra( m_meta_data.declare_part( "edge_extra" , Edge ) ),
    m_num_edge_per_proc( num_edge_per_proc ),
    m_node_ids(),
    m_edge_ids()
{
  if ( use_edge_parts ) {
    m_edge_parts.resize( num_edge_per_proc );
    for ( unsigned i = 0 ; i < num_edge_per_proc ; ++i ) {
      std::ostringstream name ;
      name << "EdgePart_" << i ;
      m_edge_parts[i] = & m_meta_data.declare_part( name.str() , Edge );
    }
  }

  m_meta_data.commit();
}

RingMeshFixture::~RingMeshFixture()
{}

void RingMeshFixture::generate_loop( bool generate_aura )
{
  const unsigned p_rank     = m_bulk_data.parallel_rank();
  const unsigned p_size     = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_edge_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;
  const unsigned n_extra    = generate_aura && 1 < p_size ? 2 : 0 ;

  m_node_ids.resize( id_total );
  m_edge_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    m_node_ids[i] = i + 1;
    m_edge_ids[i] = i + 1;
  }

  m_bulk_data.modification_begin();

  // Create a loop of edges:
  {
    const PartVector no_parts ;
    PartVector add_parts ;

    if ( ! m_edge_parts.empty() ) { add_parts.resize(1); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_edge_parts.empty() ) {
        add_parts[0] = m_edge_parts[ i % m_edge_parts.size() ];
      }
      Entity & e_node_0 = m_bulk_data.declare_entity( 0 , m_node_ids[n0] , no_parts );
      Entity & e_node_1 = m_bulk_data.declare_entity( 0 , m_node_ids[n1] , no_parts );
      Entity & e_edge   = m_bulk_data.declare_entity( 1 , m_edge_ids[i] , add_parts );
      m_bulk_data.declare_relation( e_edge , e_node_0 , 0 );
      m_bulk_data.declare_relation( e_edge , e_node_1 , 1 );
    }
  }

  Selector select_owned( m_bulk_data.mesh_meta_data().locally_owned_part() );
  Selector select_used = m_bulk_data.mesh_meta_data().locally_owned_part() |
                         m_bulk_data.mesh_meta_data().globally_shared_part();
  Selector select_all(  m_bulk_data.mesh_meta_data().universal_part() );

  stk::mesh::count_entities( select_used , m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[Node] == nLocalNode );
  STKUNIT_ASSERT( local_count[Edge] == nLocalEdge );

  std::vector<Entity*> all_nodes;
  get_entities( m_bulk_data, Node, all_nodes);

  unsigned num_selected_nodes =
      count_selected_entities( select_used, m_bulk_data.buckets(Node) );
  STKUNIT_ASSERT( num_selected_nodes == local_count[Node] );

  std::vector<Entity*> universal_nodes;
  get_selected_entities( select_all, m_bulk_data.buckets(Node), universal_nodes );
  STKUNIT_ASSERT( universal_nodes.size() == all_nodes.size() );

  STKUNIT_ASSERT( UnitTestBulkData::modification_end( m_bulk_data , generate_aura ) );

  // Verify declarations and sharing two end nodes:

  stk::mesh::count_entities( select_used , m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  if ( 1 < p_size ) {
    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    Entity * const node0 = m_bulk_data.get_entity( Node , m_node_ids[n0] );
    Entity * const node1 = m_bulk_data.get_entity( Node , m_node_ids[n1] );

    STKUNIT_ASSERT( node0 != NULL );
    STKUNIT_ASSERT( node1 != NULL );

    STKUNIT_ASSERT_EQUAL( node0->sharing().size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( node1->sharing().size() , size_t(1) );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( m_bulk_data.modification_begin() );
  m_bulk_data.change_entity_owner( change );
  STKUNIT_ASSERT( UnitTestBulkData::modification_end( m_bulk_data , generate_aura ) );

  stk::mesh::count_entities( select_used , m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  stk::mesh::count_entities( select_all , m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank()
  if ( 1 < p_size ) {
    Entity * const e_node_0 = m_bulk_data.get_entity( 0 , m_node_ids[id_begin] );
    if ( p_rank == e_node_0->owner_rank() ) {
      EntityProc entry ;
      entry.first = e_node_0 ;
      entry.second = ( p_rank + p_size - 1 ) % p_size ;
      change.push_back( entry );
    }
    STKUNIT_ASSERT( m_bulk_data.modification_begin() );
    m_bulk_data.change_entity_owner( change );
    STKUNIT_ASSERT( UnitTestBulkData::modification_end( m_bulk_data , generate_aura ) );

    stk::mesh::count_entities( select_all , m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    stk::mesh::count_entities( select_used , m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    stk::mesh::count_entities( select_owned , m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nPerProc );
    STKUNIT_ASSERT( local_count[1] == nPerProc );
  }
}


void RingMeshFixture::test_shift_loop( bool generate_aura )
{
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned nPerProc = m_num_edge_per_proc ;
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  const unsigned p_send  = ( p_rank + 1 ) % p_size ;
  const unsigned id_send = id_end - 2 ;
  const unsigned id_recv = ( id_begin + id_total - 2 ) % id_total ;

  Selector select_used = m_meta_data.locally_owned_part() |
                         m_meta_data.globally_shared_part();

  std::vector<unsigned> local_count ;
  std::vector<EntityProc> change ;

  Entity * send_edge_1 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_send ] );
  Entity * send_edge_2 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_send + 1 ] );
  Entity * send_node_1 = send_edge_1->relations()[1].entity();
  Entity * send_node_2 = send_edge_2->relations()[1].entity();
  Entity * recv_edge_1 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_recv ] );
  Entity * recv_edge_2 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_recv + 1 ] );

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

  STKUNIT_ASSERT( m_bulk_data.modification_begin() );
  m_bulk_data.change_entity_owner( change );
  STKUNIT_ASSERT( UnitTestBulkData::modification_end( m_bulk_data , generate_aura ) );

  send_edge_1 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_send ] );
  send_edge_2 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_send + 1 ] );
  recv_edge_1 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_recv ] );
  recv_edge_2 = m_bulk_data.get_entity( 1 , m_edge_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( NULL == send_edge_1 || p_rank != send_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL == send_edge_2 || p_rank != send_edge_2->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_1 && p_rank == recv_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_2 && p_rank == recv_edge_2->owner_rank() );

  stk::mesh::count_entities( select_used , m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  unsigned count_shared = 0 ;
  for ( std::vector<Entity*>::const_iterator
        i = m_bulk_data.entity_comm().begin() ;
        i != m_bulk_data.entity_comm().end() ; ++i ) {
    if ( in_shared( **i ) ) { ++count_shared ; }
  }
  STKUNIT_ASSERT( count_shared == 2u );

  {
    Entity * const node_recv = m_bulk_data.get_entity( Node , m_node_ids[id_recv] );
    Entity * const node_send = m_bulk_data.get_entity( Node , m_node_ids[id_send] );

    STKUNIT_ASSERT( node_recv->sharing().size() == 1 );
    STKUNIT_ASSERT( node_send->sharing().size() == 1 );
  }
}

