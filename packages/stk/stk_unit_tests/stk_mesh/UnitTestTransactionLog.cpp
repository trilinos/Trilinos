// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetBuckets.hpp>


#include <stk_mesh/fixtures/BoxFixture.hpp>

/** \brief When a mesh bulk data is created, the transaction is placed
 * in the bulk transaction state by default.
 */

TEST(UnitTestTransaction, verifyBulkOnCreate)
{
  stk::mesh::MetaData meta ( stk::mesh::fem_entity_rank_names() );
  stk::mesh::Part  &new_part = meta.declare_part ( "another part" );
  meta.commit ();

  stk::ParallelMachine comm(MPI_COMM_WORLD);
  stk::mesh::BulkData bulk ( meta , comm , 100 );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &new_part );

  int  size , rank;
  rank = stk::parallel_machine_rank( comm );
  size = stk::parallel_machine_size( comm );

  for ( int i = 0 ; i != 100 ; i++ )
  {
    int new_id = size*i+rank;
    bulk.declare_entity ( 0 , new_id+1 , add_part );
  }
  bulk.modification_end();

  // If something shows up in the insert incremental log, then
  // not in bulk state
  ASSERT_TRUE ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );
  ASSERT_TRUE ( bulk.get_transaction_log().get_modified_buckets(0).size() == 0 );
  ASSERT_TRUE ( bulk.get_transaction_log().get_deleted_buckets(0).size() == 0 );

  // Verify that things are inserted in the bulk transaction
  stk::mesh::PartVector  inserted_parts;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  ASSERT_TRUE ( inserted_parts.size() > 0 );
}

/** \brief When the transaction log is in bulk transaction mode,
 * inserting entities in the mesh will add the parts of that entity to
 * the inserted parts set.
 */

TEST(UnitTestTransaction, verifyBulkInsert)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part ();
  stk::mesh::PartVector                  add_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  // Verify the entity did not go into the log explicitly
  ASSERT_TRUE ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );

  // Check to see if the new_part is in the inserted parts;
  stk::mesh::PartVector                  inserted_parts;
  bool  found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != inserted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( found );

  // Verify there is nothing in the modified_parts set
  stk::mesh::PartVector  modified_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != modified_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( !found );

  // Verify there is nothing in the deleted_parts set
  stk::mesh::PartVector  deleted_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( !found );
}


/** \brief When the transaction log is in bulk transaction mode,
 * modifying the state of entities in the mesh will add the parts
 * of that entity to the inserted parts set.
 */

TEST(UnitTestTransaction, verifyBulkModify)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData    &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part        &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin ();
  stk::mesh::Entity new_entity = bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.change_entity_parts ( new_entity , blank_part , add_part );
  bulk.modification_end ();

  ASSERT_TRUE ( bulk.get_transaction_log().get_modified_buckets(0).size() == 0 );

  stk::mesh::PartVector  inserted_parts;
  bool  found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != inserted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( !found );

  stk::mesh::PartVector  modified_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != modified_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( found );

  stk::mesh::PartVector  deleted_parts;
  found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      found = true;
  }
  ASSERT_TRUE ( !found );
}

/** \brief When the transaction log is in bulk transaction mode,
 * modifying the relations of entities in the mesh will add the parts
 * of that entity to the inserted parts set.
 */

TEST(UnitTestTransaction, verifyBulkAddRelation)
{
  stk::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part, buffer_vec;
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin();
  stk::mesh::Entity new_node = bulk.declare_entity ( 0 , 123456789 , blank_part );
  stk::mesh::Entity existing_cell = *bulk.buckets(3)[0]->begin();
  bulk.declare_relation ( existing_cell, new_node , 10 );
  bulk.modification_end();

  // Verify that nodes were inserted.
  log.get_parts_with_inserted_entities ( buffer_vec );
  ASSERT_TRUE ( buffer_vec.size() > 0u );

  // Verify that the element is modified
  buffer_vec.clear();
  log.get_parts_with_modified_entities ( buffer_vec );
  ASSERT_TRUE ( buffer_vec.size() > 0u );

}


/** \brief When the transaction log is in bulk transaction mode,
 * deleting entities in the mesh will add the parts of that entity to the
 * deleted parts set.
 */

TEST(UnitTestTransaction, verifyBulkDelete)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );


  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin ();
  stk::mesh::Entity new_entity = &bulk.declare_entity ( 0 , fixture.comm_size()*1000 + fixture.comm_rank() , add_part );
  bulk.modification_end ();

  bulk.reset_transaction ( stk::mesh::Transaction::BULK );
  bulk.modification_begin ();
  bulk.destroy_entity ( *new_entity );
  bulk.modification_end ();

  ASSERT_TRUE ( bulk.get_transaction_log().get_deleted_buckets(0).size() == 0 );

  stk::mesh::PartVector                  inserted_parts;
  stk::mesh::PartVector                  modified_parts;
  stk::mesh::PartVector                  deleted_parts;
  bool  inserted_found = false;
  bulk.get_transaction_log().get_parts_with_inserted_entities ( inserted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( inserted_parts[i] == &new_part )
      inserted_found = true;
  }
  ASSERT_TRUE ( !inserted_found );

  bool modified_found = false;
  bulk.get_transaction_log().get_parts_with_modified_entities ( modified_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( modified_parts[i] == &new_part )
      modified_found = true;
  }
  ASSERT_TRUE ( !modified_found );

  bool deleted_found = false;
  bulk.get_transaction_log().get_parts_with_deleted_entities ( deleted_parts );
  for ( size_t i = 0 ; i != deleted_parts.size() ; i++ )
  {
    if ( deleted_parts[i] == &new_part )
      deleted_found = true;
  }
  ASSERT_TRUE ( deleted_found );
}

/** \brief It is possible span transactions across modifications.
 * reset_transaction must be called to remove entities from a
 * transcation,
 */

TEST(UnitTestTransaction, verifyTransactionSpanningModifications)
{
  //  HCE 3/4/10:
  //  For transactions to span multiple modifications destroyed
  //  mesh entities would have to be retained across multiple
  //  transactions.  This creates a problem where the transaction
  //  has to take ownership of the mesh entities away from the
  //  creating bulk data.
  //  This capability needs to be re-thought.

  return ;


  stk::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;


  // Here are two modifications.  The first adds an edge to the mesh,
  // the second changes the state of a node
  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  bulk.declare_entity ( 1 , 10001 , blank_part );
  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Entity n = *(*bulk.buckets(0).begin())->begin();
  bulk.change_entity_parts ( n , add_part );
  bulk.modification_end();


  // Verify both changes are logged
  ASSERT_TRUE ( bulk.get_transaction_log().get_modified_buckets(0).size() == 1 );
  ASSERT_TRUE ( bulk.get_transaction_log().get_inserted_buckets(1).size() == 1 );

  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  // Verify the log is cleared
  ASSERT_TRUE ( bulk.get_transaction_log().get_inserted_buckets(0).size() == 0 );
  ASSERT_TRUE ( bulk.get_transaction_log().get_inserted_buckets(1).size() == 0 );


  // Cannot end a transaction while the mesh is modifiable
  // Even though the transaction can span modifications, it cannot be
  // reset in the middle of a modification
  bulk.modification_begin();
  ASSERT_THROW ( bulk.reset_transaction () , std::runtime_error );
  bulk.modification_end();

}


/** \brief During an incremental
 */


TEST(UnitTestTransaction, verifyIncrementalInsert)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData             &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  // Add 4 entities to the mesh
  stk::mesh::Entity entities[4];
  entities[0] = &bulk.declare_entity ( 0 , 123456789 , blank_part );
  entities[1] = &bulk.declare_entity ( 1 , 123456789 , blank_part );
  entities[2] = &bulk.declare_entity ( 2 , 123456789 , blank_part );
  entities[3] = &bulk.declare_entity ( 3 , 123456789 , blank_part );

  // Modify one entity to ensure modification does not appear in log
  bulk.change_entity_parts ( *entities[1] , add_part );

  // Delete one entity to ensure the entity disappears from log
  bulk.destroy_entity ( *entities[3] );
  bulk.modification_end();

  // The first three entities should exist in the insert buckets in
  // the transaction log
  for ( unsigned i = 0 ; i != 3 ; i++ )
  {
    // Make sure there is only one bucket
    ASSERT_EQ ( log.get_inserted_buckets(i).size() , 1u );
    // Make sure the entity is the only thing in the bucket
    ASSERT_EQ ( log.get_inserted_buckets(i)[0]->size() , 1u );

    stk::mesh::Entity new_entity = *((*log.get_inserted_buckets(i).begin())->begin());
    // Make sure we find the right entity
    ASSERT_EQ ( &new_entity , entities[i] );
    // Verify nothing happend to modified and deleted
    ASSERT_EQ ( log.get_modified_buckets(i).size() , 0u );
    ASSERT_EQ ( log.get_deleted_buckets(i).size() , 0u );
  }

  // Verify entities[3] disappeared from the log
  ASSERT_EQ ( log.get_modified_buckets(3).size() , 0u );
  ASSERT_EQ ( log.get_deleted_buckets(3).size() , 0u );
  ASSERT_EQ ( log.get_inserted_buckets(3).size() , 0u );
}


TEST(UnitTestTransaction, verifyIncrementalModify)
{
  stk::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  // Modify the state of a node and entity in the mesh
  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk::mesh::Entity entities[2];
  entities[0] = &*bulk.buckets(0)[0]->begin();
  entities[1] = &*bulk.buckets(3)[0]->begin();
  bulk.change_entity_parts ( *entities[0] , add_part );
  bulk.change_entity_parts ( *entities[1] , add_part );
  bulk.modification_end();

  for ( unsigned i = 0 ; i != 2 ; i++ )
  {
    unsigned enttype = i*3;
    // Make sure there is only one bucket
    ASSERT_EQ ( log.get_modified_buckets(enttype).size() , 1u );
    // Make sure the entity is the only thing in the bucket
    ASSERT_EQ ( log.get_modified_buckets(enttype)[0]->size() , 1u );
    stk::mesh::Entity mod_entity = *log.get_modified_buckets(enttype)[0]->begin();
    // Make sure we find the right entity
    ASSERT_EQ ( &mod_entity , entities[i] );
    // Verify nothing happend to modified and deleted
    ASSERT_EQ ( log.get_inserted_buckets(enttype).size() , 0u );
    ASSERT_EQ ( log.get_deleted_buckets(enttype).size() , 0u );

    // Verify the transaction recorded the modification accurately
    //  1)  Make sure the new part is not part of the previous parts
    //  2)  Make sure the previous parts are in the new parts
    ASSERT_TRUE ( mod_entity.transaction_bucket() != 0 );
    ASSERT_TRUE ( !mod_entity.transaction_bucket()->member ( new_part ) );
    stk::mesh::PartVector  modified_bucket_parts;
    mod_entity.transaction_bucket()->supersets ( modified_bucket_parts );
    ASSERT_TRUE ( mod_entity.bucket().member_all ( modified_bucket_parts ));
  }
}


TEST(UnitTestTransaction, verifyIncrementalAddRelation)
{
  stk::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk::mesh::Entity new_node = bulk.declare_entity ( 0 , 123456789 , blank_part );
  stk::mesh::Entity existing_cell = *bulk.buckets(3)[0]->begin();
  bulk.declare_relation ( existing_cell, new_node , 10 );
  bulk.modification_end();

  // Verify that no nodes were modified, only inserted.
  ASSERT_EQ ( log.get_inserted_buckets(0).size() , 1u );
  ASSERT_EQ ( log.get_inserted_buckets(0)[0]->size() , 1u );
  ASSERT_EQ ( log.get_modified_buckets(0).size() , 0u );
  ASSERT_EQ ( &*log.get_inserted_buckets(0)[0]->begin() , &new_node );

  // Verify that the element is modified
  ASSERT_EQ ( log.get_modified_buckets(3).size() , 1u );
  ASSERT_EQ ( log.get_modified_buckets(3)[0]->size() , 1u );
  ASSERT_EQ ( &*log.get_modified_buckets(3)[0]->begin() , &existing_cell );

  // Make sure the parts have not changed for the existing cell
  stk::mesh::PartVector old_parts , new_parts;
  ASSERT_TRUE ( existing_cell.transaction_bucket() != 0 );
  existing_cell.transaction_bucket()->supersets ( old_parts );
  ASSERT_TRUE ( existing_cell.bucket().member_all ( old_parts ) );
  existing_cell.bucket().supersets ( new_parts );
  ASSERT_TRUE ( existing_cell.transaction_bucket()->member_all ( new_parts ) );

}


TEST(UnitTestTransaction, verifyIncrementalDelete)
{
  stk::mesh::fixtures::BoxFixture fixture;
  fixture.generate_boxes();
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,old_parts;
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test need only run in serial
  if ( fixture.comm_size() > 1 ) return;

  // destroy does not delete.  element will not be deleted until next
  // transaction reset
  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin();
  stk::mesh::Entity deleted_cell = &*bulk.buckets(3)[0]->begin();

  // Record the old parts for testing later
  deleted_cell->bucket().supersets ( old_parts );
  stk::mesh::EntityId  deleted_cell_id = deleted_cell->identifier();
  bulk.destroy_entity ( *deleted_cell );
  bulk.modification_end();

  // Verify that the element is deleted
  ASSERT_EQ ( log.get_deleted_buckets(3).size() , 1u );
  ASSERT_EQ ( log.get_deleted_buckets(3)[0]->size() , 1u );
  ASSERT_EQ ( (*log.get_deleted_buckets(3)[0]->begin()).identifier() , deleted_cell_id );

  // Check for the old parts
  deleted_cell = &*log.get_deleted_buckets(3)[0]->begin();
  ASSERT_TRUE ( deleted_cell->transaction_bucket() != 0 );
  ASSERT_TRUE ( deleted_cell->transaction_bucket()->member_all ( old_parts ) );
  stk::mesh::PartVector  old_in_trans;
  deleted_cell->transaction_bucket()->supersets ( old_in_trans );
  ASSERT_EQ ( old_in_trans.size() , old_parts.size() );
}

TEST(UnitTestTransaction, verifyParallelChangeOwnership)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();  // Comes out of fixture in MODIFIABLE

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  stk::mesh::PartVector   add_part,blank_part;
//  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  add_part.push_back ( &new_part );

  // This test needs four processes to work
  if ( fixture.comm_size() < 4 ) return;

  bulk.modification_begin ();
  stk::mesh::Entity entity = 0;
  bulk.declare_entity ( 0 , fixture.comm_rank()+1 , blank_part );
  if ( fixture.comm_rank() < 3 )
    entity = &bulk.declare_entity ( 0 , 1234 , blank_part );
  bulk.modification_end();

  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  std::vector <stk::mesh::EntityProc> change_owner;
  if ( entity )
    if ( fixture.comm_rank() == bulk.parallel_owner_rank(entity) )
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
      ASSERT_TRUE ( entity->transaction_bucket()->transaction_state() == stk::mesh::Transaction::MODIFIED );
    }
  }
  ******************/
}

TEST(UnitTestTransaction, verifyParallelResolutionModify)
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData                   &bulk = fixture.bulk_data();
  bulk.modification_end();

  stk::mesh::Part                       &new_part = fixture.get_test_part();
  const stk::mesh::MetaData             &meta = fixture.meta_data();
  const stk::mesh::Transaction   &log = bulk.get_transaction_log();
  stk::mesh::PartVector   add_part,old_parts;
  add_part.push_back ( &new_part );

  // This test need only run in parallel
  if ( fixture.comm_size() == 1 ) return;
  fixture.generate_boxes ();


  // Find a node to alter, preferable one that is shared
  const std::vector<stk::mesh::EntityProc> &shared_entities = bulk.shared_entities();
  stk::mesh::Entity node_to_modify = 0;
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
  stk::mesh::EntityId  *found_node_id_list  = new stk::mesh::EntityId [ bulk.parallel_size() ];

#ifdef STK_HAS_MPI
  stk::mesh::EntityId node_id = node_to_modify ? node_to_modify->identifier() : 0;
  int found_a_node = node_to_modify ? 1 : 0;

  MPI_Allgather ( &found_a_node , 1 , MPI_INT , found_node_list , 1 , MPI_INT , bulk.parallel() );
  MPI_Allgather ( &node_id , 1 , MPI_INT , found_node_id_list , 1 , MPI_INT , bulk.parallel() );
#endif

  // Modify the node
  bulk.reset_transaction ( stk::mesh::Transaction::INCREMENTAL );
  bulk.modification_begin ();
  if ( node_to_modify )
    bulk.change_entity_parts ( *node_to_modify , add_part );
  bulk.modification_end ();

  // Verify parallel consistent modification
  // First, loop over everythin in the modified buckets
  stk::mesh::BucketVector::const_iterator  cur_modified_node_bucket = log.get_modified_buckets(0).begin();
  while ( cur_modified_node_bucket != log.get_modified_buckets(0).end() )
  {
    stk::mesh::BucketIterator  cur_modified_node = (*cur_modified_node_bucket)->begin();
    while ( cur_modified_node != (*cur_modified_node_bucket)->begin() )
    {
      // For everything located in the buckets, verify it was changed
      // by another process
      bool valid_change = false;
      for ( unsigned i = 0 ; i != bulk.parallel_size() ; i++ )
        if ( found_node_list[i] == 1 )
          if ( cur_modified_node->identifier() == found_node_id_list[i] )
            valid_change = true;
      ASSERT_TRUE ( valid_change );
      ++cur_modified_node;
    }
    ++cur_modified_node_bucket;
  }

  delete [] found_node_list;
  delete [] found_node_id_list;
}

#endif

