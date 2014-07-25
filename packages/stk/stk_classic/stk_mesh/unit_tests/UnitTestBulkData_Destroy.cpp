/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fixtures/RingFixture.hpp>

#include <unit_tests/UnitTestModificationEndWrapper.hpp>

using stk_classic::mesh::Part;
using stk_classic::mesh::MetaData;
using stk_classic::mesh::BulkData;
using stk_classic::mesh::Entity;
using stk_classic::mesh::Selector;
using stk_classic::mesh::PartVector;
using stk_classic::mesh::EntityId;
using stk_classic::mesh::fixtures::RingFixture;

//----------------------------------------------------------------------
// Testing for mesh entities without relations

STKUNIT_UNIT_TEST(UnitTestingOfBulkData, testDestroy_nodes)
{
  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };
  const unsigned p_rank = stk_classic::parallel_machine_rank( pm );
  const unsigned p_size = stk_classic::parallel_machine_size( pm );
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );

  const int spatial_dimension = 3;
  MetaData meta( stk_classic::mesh::fem::entity_rank_names(spatial_dimension) );

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
}

//----------------------------------------------------------------------

void assert_is_destroyed( const Entity * const entity )
{
  STKUNIT_ASSERT( entity == NULL || entity->bucket().capacity() == 0 );
}

STKUNIT_UNIT_TEST(UnitTestingOfBulkData, testDestory_ring)
{
  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };
  const unsigned p_rank = stk_classic::parallel_machine_rank( pm );
  const unsigned p_size = stk_classic::parallel_machine_size( pm );
  // const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  const int spatial_dimension = 3;
  MetaData meta( stk_classic::mesh::fem::entity_rank_names(spatial_dimension) );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used = meta.locally_owned_part() | meta.globally_shared_part();
  Selector select_all(  meta.universal_part() );

  PartVector no_parts ;

  std::vector<unsigned> local_count ;

  //------------------------------
  { // No ghosting
    const bool aura_flag = false ;

    RingFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    STKUNIT_ASSERT(stk_classic::unit_test::modification_end_wrapper(bulk,
                                                           false /*no aura*/));

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    STKUNIT_ASSERT(stk_classic::unit_test::modification_end_wrapper(bulk,
                                                           false /*no aura*/));

    // This process' first element in the loop
    // if a parallel mesh has a shared node
    Entity * edge = bulk.get_entity( 1 , mesh.m_edge_ids[ nLocalEdge * p_rank ] );
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

    bulk.modification_begin();

    // Destroy the element:
    bool result = bulk.destroy_entity( edge );
    STKUNIT_ASSERT( true == result );
    STKUNIT_ASSERT( NULL == edge );

    // Destroy orphanned node:
    if ( node0->relations().size() == 0 ) {
      STKUNIT_ASSERT( bulk.destroy_entity( node0 ) );
      STKUNIT_ASSERT( NULL == node0 );
    }
    if ( node1->relations().size() == 0 ) {
      STKUNIT_ASSERT( bulk.destroy_entity( node1 ) );
      STKUNIT_ASSERT( NULL == node1 );
    }
    STKUNIT_ASSERT( stk_classic::unit_test::modification_end_wrapper(bulk, aura_flag) );

    if ( NULL != node0 ) {
      STKUNIT_ASSERT_EQUAL( node0_edges - 1 , node0->relations().size() );
    }
    if ( NULL != node1 ) {
      STKUNIT_ASSERT_EQUAL( node1_edges - 1 , node1->relations().size() );
    }
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    STKUNIT_ASSERT( bulk.modification_end() );

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    STKUNIT_ASSERT( bulk.modification_end() );

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity * node = bulk.get_entity( 0 , mesh.m_node_ids[ nNotOwned ] );

    STKUNIT_ASSERT( node != NULL );
    STKUNIT_ASSERT_NE( p_rank , node->owner_rank() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(2) , node->relations().size() );

    EntityId node_edge_ids[2] ;
    node_edge_ids[0] = node->relations()[0].entity()->identifier();
    node_edge_ids[1] = node->relations()[1].entity()->identifier();

    bulk.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached edges.
    // One will be owned and the other ghosted

    while ( node->relations().size() ) {
      Entity * e = node->relations().back().entity();
      STKUNIT_ASSERT( bulk.destroy_entity( e ) );
    }
    STKUNIT_ASSERT( bulk.destroy_entity( node ) );

    STKUNIT_ASSERT( bulk.modification_end() );

    assert_is_destroyed( bulk.get_entity(0, mesh.m_node_ids[nNotOwned] ) );
    assert_is_destroyed( bulk.get_entity(1, node_edge_ids[0] ) );
    assert_is_destroyed( bulk.get_entity(1, node_edge_ids[1] ) );

    // assert that no entities are shared or ghosted
    STKUNIT_ASSERT( bulk.entity_comm().empty() );
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No edge parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    STKUNIT_ASSERT( bulk.modification_end() );

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    STKUNIT_ASSERT( bulk.modification_end() );

    // The owned shared entity:
    const unsigned nOwned = ( nPerProc * ( p_rank + 1 ) ) % mesh.m_node_ids.size();
    const unsigned nNotOwned = nPerProc * p_rank ;

    Entity * node_owned = bulk.get_entity( 0 , mesh.m_node_ids[ nOwned ] );
    Entity * node_not_owned = bulk.get_entity( 0 , mesh.m_node_ids[ nNotOwned ] );

    STKUNIT_ASSERT( node_owned != NULL );
    STKUNIT_ASSERT( node_not_owned != NULL );
    STKUNIT_ASSERT_NE( p_rank , node_not_owned->owner_rank() );
    STKUNIT_ASSERT_EQUAL( p_rank , node_owned->owner_rank() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node_owned->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(1) , node_not_owned->sharing().size() );
    STKUNIT_ASSERT_EQUAL( size_t(2) , node_owned->relations().size() );

    EntityId node_edge_ids[2] ;
    node_edge_ids[0] = node_owned->relations()[0].entity()->identifier();
    node_edge_ids[1] = node_owned->relations()[1].entity()->identifier();

    bulk.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached edges.
    // One will be owned and the other ghosted

    while ( node_owned->relations().size() ) {
      Entity * e = node_owned->relations().back().entity();
      STKUNIT_ASSERT( bulk.destroy_entity( e ) );
    }
    STKUNIT_ASSERT( bulk.destroy_entity( node_owned ) );

    STKUNIT_ASSERT( bulk.modification_end()  );

    // Ownership of the other process' owned, shared, and destroyed node
    // has been transferred to this process.

    STKUNIT_ASSERT_EQUAL( p_rank , node_not_owned->owner_rank() );
    assert_is_destroyed( bulk.get_entity(0, mesh.m_node_ids[ nOwned ] ) );
    assert_is_destroyed( bulk.get_entity(1, node_edge_ids[0] ) );
    assert_is_destroyed( bulk.get_entity(1, node_edge_ids[1] ) );

    // assert that no entities are shared or ghosted
    STKUNIT_ASSERT( bulk.entity_comm().empty() );
  }
}

