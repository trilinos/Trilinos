/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

#include <unit_tests/UnitTestBulkData.hpp>
#include <unit_tests/UnitTestRingMeshFixture.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------
// Testing for mesh entities without relations

void UnitTestBulkData::testDestroy_nodes( ParallelMachine pm )
{
  enum { nPerProc = 10 };
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );

  MetaData meta( fem_entity_rank_names() );

  const PartVector no_parts ;

  meta.commit();

  BulkData bulk( meta , pm , 100 );

  // Ids for all entities (all entities have type 0):

  std::vector<EntityId> ids( id_total );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    ids[i] = i + 1;
  }

  // Declare just those entities in my range of ids:

  STKUNIT_ASSERT( bulk.modification_begin() );
  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    bulk.declare_entity( 0 , ids[i] , no_parts );
  }
  STKUNIT_ASSERT( bulk.modification_end() );

  // Verify that I only have entities in my range:

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity * e = bulk.get_entity( 0 , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      STKUNIT_ASSERT( NULL != e );
      STKUNIT_ASSERT( p_rank == e->owner_rank() );
    }
    else {
      STKUNIT_ASSERT( NULL == e );
    }
  }

  // Delete one entity at a time.

  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    Entity * e = bulk.get_entity( 0 , ids[ i ] );

    STKUNIT_ASSERT( NULL != e );

    bulk.modification_begin();
    STKUNIT_ASSERT( bulk.destroy_entity( e ) );
    bulk.modification_end();

    // Due to change logging the previously deleted entity
    // should be gone, but the currently deleted entity
    // should exist in the 'nil' set.

    if ( id_begin < i ) {
      STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[ i - 1 ] ) );
    }

    e = bulk.get_entity( 0 , ids[ i ] );
    STKUNIT_ASSERT( NULL != e );
    STKUNIT_ASSERT( 0 == e->bucket().capacity() );
  }

  std::cout << std::endl
            << "P" << p_rank
            << ": UnitTestBulkData::testDestroy_nodes( NP = "
            << p_size << " ) SUCCESSFULL " << std::endl ;
}

//----------------------------------------------------------------------

void assert_is_destroyed( const Entity * const entity )
{
  STKUNIT_ASSERT( entity == NULL || entity->bucket().capacity() == 0 );
}

void UnitTestBulkData::testDestroy_loop( ParallelMachine pm )
{
  enum { nPerProc = 10 };
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );
  // const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  MetaData meta( fem_entity_rank_names() );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used = meta.locally_owned_part() | meta.globally_shared_part();
  Selector select_all(  meta.universal_part() );

  PartVector no_parts ;

  std::vector<unsigned> local_count ;

  //------------------------------
  { // No ghosting
    const bool aura_flag = false ;

    UnitTestRingMeshFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    mesh.generate_mesh( aura_flag );

    // This process' first element in the loop
    // if a parallel mesh has a shared node
    Entity * edge = mesh.m_bulk_data.get_entity( 1 , mesh.m_edge_ids[ nLocalEdge * p_rank ] );
    Entity * node0 = edge->relations()[0].entity();
    Entity * node1 = edge->relations()[1].entity();

    const size_t node0_edges = node0->relations().size();
    const size_t node1_edges = node1->relations().size();

    STKUNIT_ASSERT( 1 <= node0_edges && node0_edges <= 2 );
    STKUNIT_ASSERT( 1 <= node1_edges && node1_edges <= 2 );

    STKUNIT_ASSERT( node0->relations()[0].entity() == edge ||
                    node0->relations()[1].entity() == edge );

    STKUNIT_ASSERT( node1->relations()[0].entity() == edge ||
                    node1->relations()[1].entity() == edge );

    mesh.m_bulk_data.modification_begin();

    // Destroy the element:
    bool result = mesh.m_bulk_data.destroy_entity( edge );
    STKUNIT_ASSERT( true == result );
    STKUNIT_ASSERT( NULL == edge );

    // Destroy orphanned node:
    if ( node0->relations().size() == 0 ) {
      STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( node0 ) );
      STKUNIT_ASSERT( NULL == node0 );
    }
    if ( node1->relations().size() == 0 ) {
      STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( node1 ) );
      STKUNIT_ASSERT( NULL == node1 );
    }
    STKUNIT_ASSERT( mesh.m_bulk_data.internal_modification_end( aura_flag ) );

    if ( NULL != node0 ) {
      STKUNIT_ASSERT_EQUAL( node0_edges - 1 , node0->relations().size() );
    }
    if ( NULL != node1 ) {
      STKUNIT_ASSERT_EQUAL( node1_edges - 1 , node1->relations().size() );
    }
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    const bool aura_flag = true ;

    UnitTestRingMeshFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    mesh.generate_mesh( aura_flag );

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity * node = mesh.m_bulk_data.get_entity( 0 , mesh.m_node_ids[ nNotOwned ] );

    STKUNIT_ASSERT( node != NULL );
    STKUNIT_ASSERT( p_rank != node->owner_rank() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(2) , node->relations().size() );

    EntityId node_edge_ids[2] ;
    node_edge_ids[0] = node->relations()[0].entity()->identifier();
    node_edge_ids[1] = node->relations()[1].entity()->identifier();

    mesh.m_bulk_data.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached edges.
    // One will be owned and the other ghosted

    while ( node->relations().size() ) {
      Entity * e = node->relations().back().entity();
      STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( e ) );
    }
    STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( node ) );

    STKUNIT_ASSERT( mesh.m_bulk_data.internal_modification_end( aura_flag ) );

    assert_is_destroyed( mesh.m_bulk_data.get_entity(0, mesh.m_node_ids[nNotOwned] ) );
    assert_is_destroyed( mesh.m_bulk_data.get_entity(1, node_edge_ids[0] ) );
    assert_is_destroyed( mesh.m_bulk_data.get_entity(1, node_edge_ids[1] ) );

    // assert that no entities are shared or ghosted
    STKUNIT_ASSERT( mesh.m_bulk_data.entity_comm().empty() );
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    const bool aura_flag = true ;

    UnitTestRingMeshFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    mesh.generate_mesh( aura_flag );

    // The owned shared entity:
    const unsigned nOwned = ( nPerProc * ( p_rank + 1 ) ) % mesh.m_node_ids.size();
    const unsigned nNotOwned = nPerProc * p_rank ;

    Entity * node_owned = mesh.m_bulk_data.get_entity( 0 , mesh.m_node_ids[ nOwned ] );
    Entity * node_not_owned = mesh.m_bulk_data.get_entity( 0 , mesh.m_node_ids[ nNotOwned ] );

    STKUNIT_ASSERT( node_owned != NULL );
    STKUNIT_ASSERT( node_not_owned != NULL );
    STKUNIT_ASSERT( p_rank != node_not_owned->owner_rank() );
    STKUNIT_ASSERT_EQUAL( p_rank , node_owned->owner_rank() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node_owned->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node_not_owned->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(2) , node_owned->relations().size() );

    EntityId node_edge_ids[2] ;
    node_edge_ids[0] = node_owned->relations()[0].entity()->identifier();
    node_edge_ids[1] = node_owned->relations()[1].entity()->identifier();

    mesh.m_bulk_data.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached edges.
    // One will be owned and the other ghosted

    while ( node_owned->relations().size() ) {
      Entity * e = node_owned->relations().back().entity();
      STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( e ) );
    }
    STKUNIT_ASSERT( mesh.m_bulk_data.destroy_entity( node_owned ) );

    STKUNIT_ASSERT( mesh.m_bulk_data.internal_modification_end( aura_flag ) );

    // Ownership of the other process' owned, shared, and destroyed node
    // has been transferred to this process.

    STKUNIT_ASSERT_EQUAL( p_rank , node_not_owned->owner_rank() );
    assert_is_destroyed( mesh.m_bulk_data.get_entity(0, mesh.m_node_ids[ nOwned ] ) );
    assert_is_destroyed( mesh.m_bulk_data.get_entity(1, node_edge_ids[0] ) );
    assert_is_destroyed( mesh.m_bulk_data.get_entity(1, node_edge_ids[1] ) );

    // assert that no entities are shared or ghosted
    STKUNIT_ASSERT( mesh.m_bulk_data.entity_comm().empty() );
  }
}

}
}

