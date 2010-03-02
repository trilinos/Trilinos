/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <unit_tests/stk_utest_macros.hpp>
#include <unit_tests/UnitTestMesh.hpp>
#include <unit_tests/UnitTestBoxMeshFixture.hpp>


// UnitTestBulkData_new is the beginnings of a refactoring of the bulk
// data unit test.  It relies on the UnitTestMesh fixture to rapidly
// create a mesh for testing.

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyAssertOwnerDeletedEntity )
{
  stk::unit_test::UnitTestMesh   fixture;

  stk::mesh::BulkData                   &bulk = fixture.nonconst_bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part ();
  stk::mesh::PartVector                  add_part;
  add_part.push_back ( &new_part );

  fixture.generate_boxes ();
  // Find a cell owned by this process
  stk::mesh::Entity  *cell_to_delete = NULL;
  stk::mesh::Entity  *cell_to_delete_copy = NULL;
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture.comm_rank() )
      {
        cell_to_delete = &*cur_entity;
        break;
      }
      cur_entity++;
    }
    cur_bucket++;
  }

  STKUNIT_ASSERT ( cell_to_delete != NULL );
  cell_to_delete_copy = cell_to_delete;
  bulk.modification_begin();
  bulk.destroy_entity ( cell_to_delete );
  STKUNIT_ASSERT_THROW ( bulk.assert_entity_owner_or_not_destroyed ( "method" , *cell_to_delete_copy , fixture.comm_rank()) , std::runtime_error );
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyAssertGoodKey )
{
  stk::unit_test::UnitTestMesh   fixture;

  stk::mesh::BulkData                   &bulk = fixture.nonconst_bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part ();
  stk::mesh::PartVector                  add_part;
  add_part.push_back ( &new_part );

  stk::mesh::EntityKey bad_key1 ( 45 , 1 );  // Bad entity rank
  stk::mesh::EntityKey bad_key2 ( 1 , 0 );   // Bad id

  STKUNIT_ASSERT_THROW ( bulk.assert_good_key ( "method" , bad_key1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW ( bulk.assert_good_key ( "method" , bad_key2 ) , std::runtime_error );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyGetEntityGuards )
{
  stk::unit_test::UnitTestMesh   fixture;

  stk::mesh::BulkData                   &bulk = fixture.nonconst_bulk_data();
  STKUNIT_ASSERT_THROW ( bulk.get_entity ( 1 , 0 ) , std::runtime_error );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyExplicitAddInducedPart )
{
  stk::unit_test::UnitTestMesh   fixture;
  stk::mesh::BulkData     &bulk = fixture.nonconst_bulk_data ();
  stk::mesh::PartVector    empty_vector;
  stk::mesh::PartVector    cell_part_vector;

  stk::mesh::Entity &new_cell = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , empty_vector );
  stk::mesh::Entity &new_node = bulk.declare_entity ( 0 , fixture.comm_rank()+1 , empty_vector );

  bulk.declare_relation ( new_cell , new_node , 1 );

  cell_part_vector.push_back ( &fixture.get_cell_part () );
  bulk.change_entity_parts ( new_cell , cell_part_vector );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_node , cell_part_vector ) , std::runtime_error );
}

/************************
 * This unit test is not possible currently because of the lack of
 * separation between internal part modification routines and public
 * part modification routines.
STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCannotRemoveFromSpecialParts )
{
  stk::unit_test::UnitTestMesh  fixture;
  stk::mesh::BulkData          &bulk = fixture.nonconst_bulk_data();
  stk::mesh::PartVector         test_parts;
  stk::mesh::PartVector         out_parts;
  stk::mesh::PartVector         empty_vector;

  stk::mesh::Entity &new_cell = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , empty_vector );
  test_parts.push_back ( &fixture.meta_data().universal_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.meta_data().locally_owned_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.meta_data().locally_used_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
}
 */


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDefaultPartAddition )
{
  stk::unit_test::UnitTestMesh    fixture;
  stk::mesh::BulkData            &bulk = fixture.nonconst_bulk_data ();

  bulk.modification_begin();
  stk::mesh::Entity &new_cell = fixture.get_new_entity ( 3 , 1 );
  bulk.modification_end();

  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().universal_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().locally_owned_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().locally_used_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangePartsSerial )
{
  stk::unit_test::UnitTestMesh    fixture;
  stk::mesh::BulkData            &bulk = fixture.nonconst_bulk_data ();
  stk::mesh::PartVector           create_parts , remove_parts , add_parts, empty_parts;

  create_parts.push_back ( &fixture.get_test_part() );
  create_parts.push_back ( &fixture.get_part_a_3() );
  remove_parts.push_back ( &fixture.get_part_a_3() );
  add_parts.push_back ( &fixture.get_part_b_superset() );
  add_parts.push_back ( &fixture.get_cell_part() );

  bulk.modification_begin();
  stk::mesh::Entity &new_cell = fixture.get_new_entity ( 3 , 1 );
  bulk.change_entity_parts ( new_cell , create_parts , empty_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , add_parts , remove_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , empty_parts , add_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_cell_part() ) );

  //Verify still a member of default parts
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().universal_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().locally_owned_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.meta_data().locally_used_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyParallelAddParts )
{
  stk::unit_test::UnitTestMesh     fixture;
  stk::mesh::BulkData             &bulk = fixture.nonconst_bulk_data ();
  stk::mesh::PartVector            add_part;

  add_part.push_back ( &fixture.get_part_a_0() );
  fixture.generate_boxes ();

  bulk.modification_begin();
  std::vector<stk::mesh::EntityProc>::const_iterator  cur_entity = bulk.shared_entities().begin();
  while ( cur_entity != bulk.shared_entities().end() )
  {
    if ( cur_entity->first->entity_type() == 0 )
      if ( cur_entity->first->owner_rank() == fixture.comm_rank() )
        bulk.change_entity_parts ( *cur_entity->first , add_part , stk::mesh::PartVector() );
    cur_entity++;
  }
  bulk.modification_end();

  cur_entity = bulk.shared_entities().begin();
  while ( cur_entity != bulk.shared_entities().end() )
  {
    if ( cur_entity->first->entity_type() == 0 )
      STKUNIT_ASSERT ( cur_entity->first->bucket().member ( fixture.get_part_a_0 () ) );
    cur_entity++;
  }
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyInducedMembership )
{
  stk::unit_test::UnitTestMesh     fixture;
  stk::mesh::BulkData             &bulk = fixture.nonconst_bulk_data ();
  stk::mesh::PartVector            create_node_parts , create_cell_parts , empty_parts;

  create_node_parts.push_back ( &fixture.get_part_a_0() );
  create_cell_parts.push_back ( &fixture.get_cell_part() );

  stk::mesh::Entity &node = fixture.get_new_entity ( 0 , 1 );
  stk::mesh::Entity &cell = fixture.get_new_entity ( 3 , 1 );

  bulk.change_entity_parts ( node , create_node_parts , stk::mesh::PartVector () );
  bulk.change_entity_parts ( cell , create_cell_parts , stk::mesh::PartVector () );
  // Add node to cell part
  bulk.declare_relation ( cell , node , 0 );
  bulk.modification_end();

  STKUNIT_ASSERT ( node.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.destroy_relation ( cell , node );
  bulk.modification_end();

  STKUNIT_ASSERT ( !node.bucket().member ( fixture.get_cell_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCanRemoveFromSetWithDifferentRankSubset )
{
  stk::unit_test::UnitTestMesh   fixture;
  stk::mesh::BulkData           &bulk = fixture.nonconst_bulk_data ();
  stk::mesh::PartVector          add_parts , remove_parts, empty_parts;

  add_parts.push_back ( &fixture.get_part_b_3() );
  add_parts.push_back ( &fixture.get_part_a_superset() );

  remove_parts.push_back ( &fixture.get_part_a_superset() );

  stk::mesh::Entity  &e = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , add_parts );
  bulk.modification_end();

  bulk.modification_begin();
  bulk.change_entity_parts ( e , empty_parts , remove_parts );
  bulk.modification_end();

  STKUNIT_ASSERT ( e.bucket().member ( fixture.get_part_b_3() ) );
  STKUNIT_ASSERT ( !e.bucket().member ( fixture.get_part_a_superset() ) );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCommonGhostingName )
{

  stk::unit_test::UnitTestMesh  fixture;
  stk::mesh::BulkData          &bulk = fixture.nonconst_bulk_data ();

  if ( fixture.comm_size() == 1 ) return;

  if ( fixture.comm_rank() == 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 1" ) , std::runtime_error );
  }
  else
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 2" ) , std::runtime_error );
  }
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyTrivialDestroyAllGhostings )
{
  stk::unit_test::UnitTestMesh  fixture;

  if ( fixture.comm_size() == 1 ) return;

  fixture.generate_boxes ();

  stk::mesh::BulkData  &bulk = fixture.nonconst_bulk_data();

  bulk.modification_begin();

  stk::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );

  // Find a cell owned by this process
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  unsigned send_rank = 0;

  std::vector<stk::mesh::EntityProc>  to_send;
  std::vector<stk::mesh::Entity *>    empty_vector;
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        send_rank++;
      }
      cur_entity++;
    }
    cur_bucket++;
  }
  bulk.change_ghosting ( ghosting , to_send , empty_vector );
  bulk.modification_end();

  STKUNIT_ASSERT ( ghosting.send().size() > 0u );
  STKUNIT_ASSERT ( ghosting.receive().size() > 0u );

  bulk.modification_begin();
  bulk.destroy_all_ghosting ();
  bulk.modification_end();

  STKUNIT_ASSERT_EQUAL ( ghosting.send().size() , 0u );
  STKUNIT_ASSERT_EQUAL ( ghosting.receive().size() , 0u );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangeGhostingGuards )
{
  stk::unit_test::UnitTestMesh   fixture1 , fixture2;
  stk::mesh::BulkData   &bulk1 = fixture1.nonconst_bulk_data ();
  stk::mesh::BulkData   &bulk2 = fixture2.nonconst_bulk_data ();

  fixture1.generate_boxes();
  fixture2.generate_boxes();

  bulk1.modification_begin();
  bulk2.modification_begin();

  std::vector<stk::mesh::EntityProc>  to_send;
  std::vector<stk::mesh::Entity *>    empty_vector;
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk1.buckets(3).begin();
  unsigned send_rank = 0;
  while ( cur_bucket != bulk1.buckets(3).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture1.comm_rank() )
      {
        if ( send_rank == fixture1.comm_size() ) send_rank = 0;
        if ( send_rank != fixture1.comm_rank() )
          to_send.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        send_rank++;
      }
      cur_entity++;
    }
    cur_bucket++;
  }

  stk::mesh::Ghosting &ghosting = bulk1.create_ghosting ( "Ghost 1" );
  STKUNIT_ASSERT_THROW ( bulk2.change_ghosting ( ghosting , to_send , empty_vector ) , std::runtime_error );
  STKUNIT_ASSERT_THROW ( bulk1.change_ghosting ( bulk1.shared_aura() , to_send , empty_vector ) , std::runtime_error );
  bulk1.modification_end();
  bulk2.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyOtherGhostingGuards )
{
  stk::unit_test::UnitTestMesh  fixture;
  stk::mesh::BulkData          &bulk = fixture.nonconst_bulk_data ();
  fixture.generate_boxes ();
  bulk.modification_begin();

  std::vector<stk::mesh::EntityProc>  to_send_unowned;
  std::vector<stk::mesh::EntityProc>  empty_send;
  std::vector<stk::mesh::Entity *>    to_remove_not_ghosted;
  std::vector<stk::mesh::Entity *>    empty_remove;
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  unsigned send_rank = 0;
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() != fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send_unowned.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        send_rank++;
      }
      else
      {
        to_remove_not_ghosted.push_back ( &*cur_entity );
      }
      cur_entity++;
    }
    cur_bucket++;
  }

  stk::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );
  if ( to_send_unowned.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove );
  }

  if ( to_remove_not_ghosted.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted );
  }
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyPartsOnCreate )
{
   stk::unit_test::UnitTestMesh    fixture;
   stk::mesh::BulkData           & bulk = fixture.nonconst_bulk_data ();
   stk::mesh::Part               & part_a = fixture.get_part_a_0 ();
   stk::mesh::Part               & part_b = fixture.get_part_b_0 ();

   stk::mesh::PartVector           create_vector;
   create_vector.push_back ( &part_a );

   stk::mesh::Entity &node = bulk.declare_entity ( 0 , fixture.comm_rank()+1 ,create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( node.bucket().member ( part_a ) );

   bulk.modification_begin();
   create_vector.push_back ( &part_b );
   stk::mesh::Entity &node2 = bulk.declare_entity ( 0 , fixture.comm_size() + fixture.comm_rank() + 1 , create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( node2.bucket().member ( part_a ) );
   STKUNIT_ASSERT ( node2.bucket().member ( part_b ) );
}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyBoxGhosting )
{
  const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD );
  if ( 8 < p_size ) { return ; }

  BoxMeshFixture fixture( MPI_COMM_WORLD );

  for ( size_t iz = 0 ; iz < 3 ; ++iz ) {
  for ( size_t iy = 0 ; iy < 3 ; ++iy ) {
  for ( size_t ix = 0 ; ix < 3 ; ++ix ) {
    stk::mesh::Entity * const node = fixture.m_nodes[iz][iy][ix] ;
    STKUNIT_ASSERT( NULL != node );
    if ( NULL != node ) {
      STKUNIT_ASSERT( fixture.m_node_id[iz][iy][ix] == node->identifier() );
      BoxMeshFixture::Scalar * const node_coord =
        stk::mesh::field_data( *fixture.m_coord_field , *node );
      STKUNIT_ASSERT( node_coord != NULL );
      if ( node_coord != NULL ) {
        BoxMeshFixture::Scalar * coord = fixture.m_node_coord[iz][iy][ix] ;
        STKUNIT_ASSERT_EQUAL( node_coord[0] , coord[0] );
        STKUNIT_ASSERT_EQUAL( node_coord[1] , coord[1] );
        STKUNIT_ASSERT_EQUAL( node_coord[2] , coord[2] );
      }
    }
  }
  }
  }

  for ( size_t iz = 0 ; iz < 2 ; ++iz ) {
  for ( size_t iy = 0 ; iy < 2 ; ++iy ) {
  for ( size_t ix = 0 ; ix < 2 ; ++ix ) {
    stk::mesh::Entity * const elem = fixture.m_elems[iz][iy][ix] ;
    STKUNIT_ASSERT( NULL != elem );
    if ( NULL != elem ) {
      stk::mesh::PairIterRelation elem_nodes = elem->relations();
      STKUNIT_ASSERT_EQUAL( 8u , elem_nodes.size() );
      BoxMeshFixture::Scalar ** const elem_node_coord =
        stk::mesh::field_data( *fixture.m_coord_gather_field , *elem );
      for ( size_t j = 0 ; j < elem_nodes.size() ; ++j ) {
        STKUNIT_ASSERT_EQUAL( j , elem_nodes[j].identifier() );
        BoxMeshFixture::Scalar * const node_coord =
          stk::mesh::field_data( *fixture.m_coord_field , *elem_nodes[j].entity() );
        STKUNIT_ASSERT( node_coord == elem_node_coord[ elem_nodes[j].identifier() ] );
      }
      if ( 8u == elem_nodes.size() ) {
        STKUNIT_ASSERT( elem_nodes[0].entity() == fixture.m_nodes[iz][iy][ix] );
        STKUNIT_ASSERT( elem_nodes[1].entity() == fixture.m_nodes[iz][iy][ix+1] );
        STKUNIT_ASSERT( elem_nodes[2].entity() == fixture.m_nodes[iz+1][iy][ix+1] );
        STKUNIT_ASSERT( elem_nodes[3].entity() == fixture.m_nodes[iz+1][iy][ix] );
        STKUNIT_ASSERT( elem_nodes[4].entity() == fixture.m_nodes[iz][iy+1][ix] );
        STKUNIT_ASSERT( elem_nodes[5].entity() == fixture.m_nodes[iz][iy+1][ix+1] );
        STKUNIT_ASSERT( elem_nodes[6].entity() == fixture.m_nodes[iz+1][iy+1][ix+1] );
        STKUNIT_ASSERT( elem_nodes[7].entity() == fixture.m_nodes[iz+1][iy+1][ix] );
      }
    }
  }
  }
  }
}

