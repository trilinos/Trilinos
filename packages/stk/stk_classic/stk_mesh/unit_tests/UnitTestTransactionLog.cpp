/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/** \TODO The design of the transaction logging functionality must
 *        be revisited / revised to clean up entity creation / deletion
 *        where two or more entities could exist with the same key.
 *        Such a condition is severely erroneous.  Until the design
 *        revisit / revision is carried out the transaction capability
 *        is disabled.
 */
#if 0

#include <sstream>
#include <stdexcept>

#include <unit_tests/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetBuckets.hpp>


#include <stk_mesh/fixtures/BoxFixture.hpp>

/** \brief When a mesh bulk data is created, the transaction is placed
 * in the bulk transaction state by default.
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyBulkOnCreate)
{
  stk_classic::mesh::MetaData meta ( stk_classic::mesh::fem_entity_rank_names() );
  stk_classic::mesh::Part  &new_part = meta.declare_part ( "another part" );
  meta.commit ();

  stk_classic::ParallelMachine comm(MPI_COMM_WORLD);
  stk_classic::mesh::BulkData bulk ( meta , comm , 100 );
  std::vector<stk_classic::mesh::Part *>  add_part;
  add_part.push_back ( &new_part );

  int  size , rank;
  rank = stk_classic::parallel_machine_rank( comm );
  size = stk_classic::parallel_machine_size( comm );

  for ( int i = 0 ; i != 100 ; i++ )
  {
    int new_id = size*i+rank;
    bulk.declare_entity ( 0 , new_id+1 , add_part );
  }
  bulk.modification_end();

  // If something shows up in the insert incremental log, then
  // not in bulk state
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_modified_buckets(0).size() == 0 );
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_deleted_buckets(0).size() == 0 );

  // Verify that things are inserted in the bulk transaction
  stk_classic::mesh::PartVector  inserted_parts;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  STKUNIT_ASSERT ( inserted_parts.size() > 0 );
}

/** \brief When the transaction log is in bulk transaction mode,
 * inserting entities in the mesh will add the parts of that entity to
 * the inserted parts set.
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyBulkInsert)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part ();
  stk_classic::mesh::PartVector                  add_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  // Verify the entity did not go into the log explicitly
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );

  // Check to see if the new_part is in the inserted parts;
  stk_classic::mesh::PartVector                  inserted_parts;
  bool  found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != inserted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( found );

  // Verify there is nothing in the modified_parts set
  stk_classic::mesh::PartVector  modified_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != modified_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( !found );

  // Verify there is nothing in the deleted_parts set
  stk_classic::mesh::PartVector  deleted_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( !found );
}


/** \brief When the transaction log is in bulk transaction mode,
 * modifying the state of entities in the mesh will add the parts
 * of that entity to the inserted parts set.
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyBulkModify)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData    &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part        &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin ();
  stk_classic::mesh::Entity &new_entity = bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.change_entity_parts ( new_entity , blank_part , add_part );
  bulk.modification_end ();

  STKUNIT_ASSERT ( bulk.get_transaction_log().get_modified_buckets(0).size() == 0 );

  stk_classic::mesh::PartVector  inserted_parts;
  bool  found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != inserted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( !found );

  stk_classic::mesh::PartVector  modified_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != modified_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( found );

  stk_classic::mesh::PartVector  deleted_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      found = true;
  }
  STKUNIT_ASSERT ( !found );
}

/** \brief When the transaction log is in bulk transaction mode,
 * modifying the relations of entities in the mesh will add the parts
 * of that entity to the inserted parts set.
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyBulkAddRelation)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part, buffer_vec;
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin();
  stk_classic::mesh::Entity  &new_node = bulk.declare_entity ( 0 , 123456789 , blank_part );
  stk_classic::mesh::Entity  &existing_cell = *bulk.buckets(3)[0]->begin();
  bulk.declare_relation ( existing_cell, new_node , 10 );
  bulk.modification_end();

  // Verify that nodes were inserted.
  log.get_parts_with_inserted_entities ( buffer_vec );
  STKUNIT_ASSERT ( buffer_vec.size() > 0u );

  // Verify that the element is modified
  buffer_vec.clear();
  log.get_parts_with_modified_entities ( buffer_vec );
  STKUNIT_ASSERT ( buffer_vec.size() > 0u );

}


/** \brief When the transaction log is in bulk transaction mode,
 * deleting entities in the mesh will add the parts of that entity to the
 * deleted parts set.
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyBulkDelete)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );


  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin ();
  stk_classic::mesh::Entity *new_entity = &bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  bulk.reset_transaction ( stk_classic::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.destroy_entity ( new_entity );
  bulk.modification_end ();

  STKUNIT_ASSERT ( bulk.get_transaction_log().get_deleted_buckets(0).size() == 0 );

  stk_classic::mesh::PartVector                  inserted_parts;
  stk_classic::mesh::PartVector                  modified_parts;
  stk_classic::mesh::PartVector                  deleted_parts;
  bool  inserted_found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      inserted_found = true;
  }
  STKUNIT_ASSERT ( !inserted_found );

  bool modified_found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      modified_found = true;
  }
  STKUNIT_ASSERT ( !modified_found );

  bool deleted_found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      deleted_found = true;
  }
  STKUNIT_ASSERT ( deleted_found );
}

/** \brief It is possible span transactions across modifications.
 * reset_transaction must be called to remove entities from a
 * transcation,
 */

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyTransactionSpanningModifications)
{
  //  HCE 3/4/10:
  //  For transactions to span multiple modifications destroyed
  //  mesh entities would have to be retained across multiple
  //  transactions.  This creates a problem where the transaction
  //  has to take ownership of the mesh entities away from the
  //  creating bulk data.
  //  This capability needs to be re-thought.

  return ;


  stk_classic::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;


  // Here are two modifications.  The first adds an edge to the mesh,
  // the second changes the state of a node
  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  bulk.declare_entity ( 1 , 10001 , blank_part );
  bulk.modification_end();

  bulk.modification_begin();
  stk_classic::mesh::Entity &n = *(*bulk.buckets(0).begin())->begin();
  bulk.change_entity_parts ( n , add_part );
  bulk.modification_end();


  // Verify both changes are logged
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_modified_buckets(0).size() == 1 );
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_inserted_buckets(1).size() == 1 );

  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  // Verify the log is cleared
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );
  STKUNIT_ASSERT ( bulk.get_transaction_log().get_inserted_buckets(1).size() == 0 );


  // Cannot end a transaction while the mesh is modifiable
  // Even though the transaction can span modifications, it cannot be
  // reset in the middle of a modification
  bulk.modification_begin();
  STKUNIT_ASSERT_THROW ( bulk.reset_transaction () , std::runtime_error );
  bulk.modification_end();

}


/** \brief During an incremental
 */


STKUNIT_UNIT_TEST(UnitTestTransaction, verifyIncrementalInsert)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  // Add 4 entities to the mesh
  stk_classic::mesh::Entity *entities[4];
  entities[0] = &bulk.declare_entity ( 0 , 123456789 , blank_part );
  entities[1] = &bulk.declare_entity ( 1 , 123456789 , blank_part );
  entities[2] = &bulk.declare_entity ( 2 , 123456789 , blank_part );
  entities[3] = &bulk.declare_entity ( 3 , 123456789 , blank_part );

  // Modify one entity to ensure modification does not appear in log
  bulk.change_entity_parts ( *entities[1] , add_part );

  // Delete one entity to ensure the entity disappears from log
  bulk.destroy_entity ( entities[3] );
  bulk.modification_end();

  // The first three entities should exist in the insert buckets in
  // the transaction log
  for ( unsigned i = 0 ; i != 3 ; i++ )
  {
    // Make sure there is only one bucket
    STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(i).size() , 1u );
    // Make sure the entity is the only thing in the bucket
    STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(i)[0]->size() , 1u );

    stk_classic::mesh::Entity &new_entity = *((*log.get_inserted_buckets(i).begin())->begin());
    // Make sure we find the right entity
    STKUNIT_ASSERT_EQUAL ( &new_entity , entities[i] );
    // Verify nothing happend to modified and deleted
    STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(i).size() , 0u );
    STKUNIT_ASSERT_EQUAL ( log.get_deleted_buckets(i).size() , 0u );
  }

  // Verify entities[3] disappeared from the log
  STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(3).size() , 0u );
  STKUNIT_ASSERT_EQUAL ( log.get_deleted_buckets(3).size() , 0u );
  STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(3).size() , 0u );
}


STKUNIT_UNIT_TEST(UnitTestTransaction, verifyIncrementalModify)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  // Modify the state of a node and entity in the mesh
  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk_classic::mesh::Entity *entities[2];
  entities[0] = &*bulk.buckets(0)[0]->begin();
  entities[1] = &*bulk.buckets(3)[0]->begin();
  bulk.change_entity_parts ( *entities[0] , add_part );
  bulk.change_entity_parts ( *entities[1] , add_part );
  bulk.modification_end();

  for ( unsigned i = 0 ; i != 2 ; i++ )
  {
    unsigned enttype = i*3;
    // Make sure there is only one bucket
    STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(enttype).size() , 1u );
    // Make sure the entity is the only thing in the bucket
    STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(enttype)[0]->size() , 1u );
    stk_classic::mesh::Entity &mod_entity = *log.get_modified_buckets(enttype)[0]->begin();
    // Make sure we find the right entity
    STKUNIT_ASSERT_EQUAL ( &mod_entity , entities[i] );
    // Verify nothing happend to modified and deleted
    STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(enttype).size() , 0u );
    STKUNIT_ASSERT_EQUAL ( log.get_deleted_buckets(enttype).size() , 0u );

    // Verify the transaction recorded the modification accurately
    //  1)  Make sure the new part is not part of the previous parts
    //  2)  Make sure the previous parts are in the new parts
    STKUNIT_ASSERT ( mod_entity.transaction_bucket() != 0 );
    STKUNIT_ASSERT ( !mod_entity.transaction_bucket()->member ( new_part ) );
    stk_classic::mesh::PartVector  modified_bucket_parts;
    mod_entity.transaction_bucket()->supersets ( modified_bucket_parts );
    STKUNIT_ASSERT ( mod_entity.bucket().member_all ( modified_bucket_parts ));
  }
}


STKUNIT_UNIT_TEST(UnitTestTransaction, verifyIncrementalAddRelation)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk_classic::mesh::Entity  &new_node = bulk.declare_entity ( 0 , 123456789 , blank_part );
  stk_classic::mesh::Entity  &existing_cell = *bulk.buckets(3)[0]->begin();
  bulk.declare_relation ( existing_cell, new_node , 10 );
  bulk.modification_end();

  // Verify that no nodes were modified, only inserted.
  STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(0).size() , 1u );
  STKUNIT_ASSERT_EQUAL ( log.get_inserted_buckets(0)[0]->size() , 1u );
  STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(0).size() , 0u );
  STKUNIT_ASSERT_EQUAL ( &*log.get_inserted_buckets(0)[0]->begin() , &new_node );

  // Verify that the element is modified
  STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(3).size() , 1u );
  STKUNIT_ASSERT_EQUAL ( log.get_modified_buckets(3)[0]->size() , 1u );
  STKUNIT_ASSERT_EQUAL ( &*log.get_modified_buckets(3)[0]->begin() , &existing_cell );

  // Make sure the parts have not changed for the existing cell
  stk_classic::mesh::PartVector old_parts , new_parts;
  STKUNIT_ASSERT ( existing_cell.transaction_bucket() != 0 );
  existing_cell.transaction_bucket()->supersets ( old_parts );
  STKUNIT_ASSERT ( existing_cell.bucket().member_all ( old_parts ) );
  existing_cell.bucket().supersets ( new_parts );
  STKUNIT_ASSERT ( existing_cell.transaction_bucket()->member_all ( new_parts ) );

}


STKUNIT_UNIT_TEST(UnitTestTransaction, verifyIncrementalDelete)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,old_parts;
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  // destroy does not delete.  element will not be deleted until next
  // transaction reset
  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk_classic::mesh::Entity     *deleted_cell = &*bulk.buckets(3)[0]->begin();

  // Record the old parts for testing later
  deleted_cell->bucket().supersets ( old_parts );
  stk_classic::mesh::EntityId  deleted_cell_id = deleted_cell->identifier();
  bulk.destroy_entity ( deleted_cell );
  bulk.modification_end();

  // Verify that the element is deleted
  STKUNIT_ASSERT_EQUAL ( log.get_deleted_buckets(3).size() , 1u );
  STKUNIT_ASSERT_EQUAL ( log.get_deleted_buckets(3)[0]->size() , 1u );
  STKUNIT_ASSERT_EQUAL ( (*log.get_deleted_buckets(3)[0]->begin()).identifier() , deleted_cell_id );

  // Check for the old parts
  deleted_cell = &*log.get_deleted_buckets(3)[0]->begin();
  STKUNIT_ASSERT ( deleted_cell->transaction_bucket() != 0 );
  STKUNIT_ASSERT ( deleted_cell->transaction_bucket()->member_all ( old_parts ) );
  stk_classic::mesh::PartVector  old_in_trans;
  deleted_cell->transaction_bucket()->supersets ( old_in_trans );
  STKUNIT_ASSERT_EQUAL ( old_in_trans.size() , old_parts.size() );
}

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyParallelChangeOwnership)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  stk_classic::mesh::PartVector   add_part,blank_part;
//  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test needs four processes to work
  if ( fixture.comm_size() < 4 ) return;

  bulk.modification_begin ();
  stk_classic::mesh::Entity  *entity = 0;
  bulk.declare_entity ( 0 , fixture.comm_rank()+1 , blank_part );
  if ( fixture.comm_rank() < 3 )
    entity = &bulk.declare_entity ( 0 , 1234 , blank_part );
  bulk.modification_end();

  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  std::vector <stk_classic::mesh::EntityProc> change_owner;
  if ( entity )
    if ( fixture.comm_rank() == entity->owner_rank() )
    {
      int other_rank = fixture.comm_rank()==0?1:0;
      change_owner.push_back ( std::make_pair ( entity , other_rank ) );
    }
  bulk.modification_begin();
  bulk.change_entity_owner ( change_owner );
  bulk.modification_end();

  /********* This needs to be fixed:  we need to know what correct
   * behavior should be
  if ( entity )
  {
    if ( fixture.comm_rank() < 3 )
    {
      STKUNIT_ASSERT ( entity->transaction_bucket()->transaction_state() == stk_classic::mesh::Transaction::MODIFIED );
    }
  }
  ******************/
}

STKUNIT_UNIT_TEST(UnitTestTransaction, verifyParallelResolutionModify)
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();

  stk_classic::mesh::Part                       &new_part = fixture.get_test_part();
  const stk_classic::mesh::MetaData             &meta = fixture.meta_data();
  const stk_classic::mesh::Transaction   &log = bulk.get_transaction_log();
  stk_classic::mesh::PartVector   add_part,old_parts;
  add_part.push_back ( &new_part );

  // This test need only run in parallel
  if ( fixture.comm_size() == 1 ) return;
  fixture.generate_boxes ();


  // Find a node to alter, preferable one that is shared
  const std::vector<stk_classic::mesh::EntityProc> &shared_entities = bulk.shared_entities();
  stk_classic::mesh::Entity  *node_to_modify = 0;
  for ( unsigned i = 0 ; i != shared_entities.size() ;i++ )
  {
    if ( shared_entities[i].first->entity_rank() == 0 )
      if ( shared_entities[i].first->bucket().member ( meta.locally_owned_part () ) )
      {
        node_to_modify = shared_entities[i].first;
        break;
      }
  }

  // Once found, tell all processes which one.  If not found, tell
  // them that as well
  int       *found_node_list     = new int [ bulk.parallel_size() ];
  stk_classic::mesh::EntityId  *found_node_id_list  = new stk_classic::mesh::EntityId [ bulk.parallel_size() ];

#ifdef STK_HAS_MPI
  stk_classic::mesh::EntityId node_id = node_to_modify ? node_to_modify->identifier() : 0;
  int found_a_node = node_to_modify ? 1 : 0;

  MPI_Allgather ( &found_a_node , 1 , MPI_INT , found_node_list , 1 , MPI_INT , bulk.parallel() );
  MPI_Allgather ( &node_id , 1 , MPI_INT , found_node_id_list , 1 , MPI_INT , bulk.parallel() );
#endif

  // Modify the node
  bulk.reset_transaction ( stk_classic::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin ();
  if ( node_to_modify )
    bulk.change_entity_parts ( *node_to_modify , add_part );
  bulk.modification_end ();

  // Verify parallel consistent modification
  // First, loop over everythin in the modified buckets
  std::vector<stk_classic::mesh::Bucket *>::const_iterator  cur_modified_node_bucket = log.get_modified_buckets(0).begin();
  while ( cur_modified_node_bucket != log.get_modified_buckets(0).end() )
  {
    stk_classic::mesh::BucketIterator  cur_modified_node = (*cur_modified_node_bucket)->begin();
    while ( cur_modified_node != (*cur_modified_node_bucket)->begin() )
    {
      // For everything located in the buckets, verify it was changed
      // by another process
      bool valid_change = false;
      for ( unsigned i = 0 ; i != bulk.parallel_size() ; i++ )
        if ( found_node_list[i] == 1 )
          if ( cur_modified_node->identifier() == found_node_id_list[i] )
            valid_change = true;
      STKUNIT_ASSERT ( valid_change );
      ++cur_modified_node;
    }
    ++cur_modified_node_bucket;
  }

  delete [] found_node_list;
  delete [] found_node_id_list;
}

#endif

