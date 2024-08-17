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

#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::iterator
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityCommListInfo.hpp"
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, operator<<
#include "stk_mesh/base/EntityLess.hpp"  // for EntityLess
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting, operator<<
#include "stk_mesh/base/Types.hpp"      // for EntityProc, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology::num_nodes
#include "stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp"  // for BoxFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <iosfwd>                       // for ostringstream, ostream
#include <set>                          // for set, etc
#include <stddef.h>                     // for size_t, NULL
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/EntityCommDatabase.hpp>  // for pack_entity_info, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field_on_mesh
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <string>                       // for string
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector, etc
namespace stk { namespace mesh { class Part; } }


using namespace stk::mesh;

// UnitTestBulkData_new is the beginnings of a refactoring of the bulk
// data unit test.  It relies on a customized BoxFixture to rapidly
// create a mesh for testing.

namespace {

/**
 * The customized box fixture used in this file for testing. This fixture
 * is similar to the BoxFixture it inherits from, with the only difference
 * being the extra parts that this fixture declares for testing purposes.
 */
struct TestBoxFixture : public fixtures::BoxFixture
{
  TestBoxFixture(stk::ParallelMachine pm = MPI_COMM_WORLD)
    : BoxFixture(pm, stk::mesh::BulkData::AUTO_AURA),
      m_test_part ( m_fem_meta.declare_part ( "Test Part" ) ),
      m_cell_part ( m_fem_meta.declare_part ( "Cell list" , stk::topology::ELEM_RANK ) ),
      m_part_A_0 ( m_fem_meta.declare_part ( "Part A 0", stk::topology::NODE_RANK ) ),
      m_part_A_1 ( m_fem_meta.declare_part ( "Part A 1", stk::topology::EDGE_RANK ) ),
      m_part_A_2 ( m_fem_meta.declare_part ( "Part A 2", stk::topology::FACE_RANK ) ),
      m_part_A_3 ( m_fem_meta.declare_part ( "Part A 3", stk::topology::ELEM_RANK ) ),
      m_part_A_superset ( m_fem_meta.declare_part ( "Part A superset" ) ),
      m_part_B_0 ( m_fem_meta.declare_part ( "Part B 0", stk::topology::NODE_RANK ) ),
      m_part_B_1 ( m_fem_meta.declare_part ( "Part B 1", stk::topology::EDGE_RANK ) ),
      m_part_B_2 ( m_fem_meta.declare_part ( "Part B 2", stk::topology::FACE_RANK ) ),
      m_part_B_3 ( m_fem_meta.declare_part ( "Part B 3", stk::topology::ELEM_RANK ) ),
      m_part_B_superset ( m_fem_meta.declare_part ( "Part B superset" ) )
  {
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_0 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_1 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_2 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_3 );

    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_0 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_1 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_2 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_3 );

    // None of the tests currently need to make any addtional changes
    // to MetaData; if this changes, the line below will have to be
    // removed.
    m_fem_meta.commit();
  }

  Part     & m_test_part;   // A simple part
  Part     & m_cell_part;   // A part to put cells in

  Part     & m_part_A_0;
  Part     & m_part_A_1;
  Part     & m_part_A_2;
  Part     & m_part_A_3;

  Part     & m_part_A_superset;

  Part     & m_part_B_0;
  Part     & m_part_B_1;
  Part     & m_part_B_2;
  Part     & m_part_B_3;

  Part     & m_part_B_superset;
};

}

TEST ( UnitTestBulkData_new , verifyAssertOwnerDeletedEntity )
{
  TestBoxFixture fixture;

  BulkData         &bulk = fixture.bulk_data();
  Part             &new_part = fixture.m_test_part;
  PartVector        add_part;
  add_part.push_back ( &new_part );

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  // Find a cell owned by this process
  Entity cell_to_delete = Entity();
  BucketVector::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) == fixture.comm_rank() )
      {
        cell_to_delete = *cur_entity;
        break;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  ASSERT_TRUE ( bulk.is_valid(cell_to_delete) );
  bulk.modification_begin();
  bulk.destroy_entity ( cell_to_delete );
  bulk.modification_end();
}


TEST ( UnitTestBulkData_new , verifyDetectsBadKey )
{
  TestBoxFixture fixture;

  BulkData         &bulk = fixture.bulk_data();
  Part             &new_part = fixture.m_test_part;
  PartVector        add_part, empty_vector;
  add_part.push_back ( &new_part );

  EntityKey bad_key1 ( static_cast<EntityRank>(45) , 1 );  // Bad entity rank
  EntityKey bad_key2 ( stk::topology::EDGE_RANK , 0 );   // Bad id

  ASSERT_THROW ( bulk.declare_entity(bad_key1.rank(),
                                     bad_key1.id(),
                                     empty_vector),
                 std::logic_error );
  ASSERT_THROW ( bulk.declare_edge(bad_key2.id(), empty_vector),
                 std::logic_error );
}

TEST ( UnitTestBulkData_new , verifyDetectsNonOwnerChange )
{
  // Set up a mesh where there are shared nodes. Take one of the nodes, and
  // have the non-owning processes try to make a change to that node; this
  // should cause an exception.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);
  int p_rank = stk::parallel_machine_rank(pm);

  fixtures::QuadFixture fixture(pm, 1 /*nx*/, p_size /*ny*/);
  fixture.m_meta.commit();
  fixture.generate_mesh();
  BulkData & bulk = fixture.m_bulk_data;

  PartVector empty_vector;

  Entity shared_node = fixture.node(1 /*x*/, 1 /*y*/);
  // Assert that this node is shared
  if ( p_size > 1 && bulk.is_valid(shared_node) && (p_rank == 0 || p_rank == 1) ) {
    std::vector<int> shared_procs;
    bulk.comm_shared_procs(bulk.entity_key(shared_node),shared_procs);
    ASSERT_GE(shared_procs.size(), 1u);
  }

  bulk.modification_begin();

  // Non-owners of shared_node will attempt to make a change to it; this should
  // cause an exception
  if (bulk.is_valid(shared_node) && p_rank != bulk.parallel_owner_rank(shared_node)) {
    ASSERT_THROW(bulk.change_entity_parts(shared_node,
                                          empty_vector,  //add parts
                                          empty_vector), //rem parts
                 std::logic_error);
  }

  bulk.modification_end();
}

TEST ( UnitTestBulkData_new , verifyExplicitAddInducedPart )
{
  TestBoxFixture fixture;
  BulkData     &bulk = fixture.bulk_data ();
  PartVector    empty_vector;
  PartVector    cell_part_vector;

  bulk.modification_begin();

  Entity new_cell = bulk.declare_element(fixture.comm_rank()+1 , empty_vector );
  Entity new_node = bulk.declare_node(fixture.comm_rank()+1 , empty_vector );

  bulk.declare_relation ( new_cell , new_node , 1 );

  cell_part_vector.push_back ( &fixture.m_cell_part );
  bulk.change_entity_parts ( new_cell , cell_part_vector );
  bulk.change_entity_parts ( new_node , cell_part_vector );
  EXPECT_TRUE(bulk.bucket(new_node).member(fixture.m_cell_part));
}

TEST ( UnitTestBulkData_new , verifyDefaultPartAddition )
{
  TestBoxFixture fixture;
  BulkData            &bulk = fixture.bulk_data ();

  bulk.modification_begin();
  Entity new_cell = bulk.declare_element(1, stk::mesh::ConstPartVector{&fixture.get_elem_part()});
  unsigned cell_num_nodes = fixture.get_elem_topology().num_nodes();
  for (unsigned i = 0; i < cell_num_nodes; ++i)
  {
    Entity new_node = bulk.declare_node(i+1);
    bulk.declare_relation(new_cell, new_node, i);
  }
  bulk.modification_end();

  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.fem_meta().universal_part() ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.fem_meta().locally_owned_part() ) );
}

TEST ( UnitTestBulkData_new , verifyChangePartsSerial )
{
  TestBoxFixture fixture;
  BulkData            &bulk = fixture.bulk_data ();
  PartVector           create_parts , remove_parts , add_parts, empty_parts;

  create_parts.push_back ( &fixture.m_test_part );
  create_parts.push_back ( &fixture.m_part_A_3 );
  remove_parts.push_back ( &fixture.m_part_A_3 );
  add_parts.push_back ( &fixture.m_part_B_superset );
  add_parts.push_back ( &fixture.m_cell_part );

  bulk.modification_begin();
  Entity new_cell = bulk.declare_element(1, stk::mesh::ConstPartVector{&fixture.get_elem_part()});
  unsigned cell_num_nodes = fixture.get_elem_topology().num_nodes();
  for (unsigned i = 0; i < cell_num_nodes; ++i)
  {
    Entity new_node = bulk.declare_node(i+1);
    bulk.declare_relation(new_cell, new_node, i);
  }
  bulk.change_entity_parts ( new_cell , create_parts , empty_parts );
  bulk.modification_end();
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , add_parts , remove_parts );
  bulk.modification_end();
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , empty_parts , add_parts );
  bulk.modification_end();
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  ASSERT_TRUE ( !bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  //Verify still a member of default parts
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.fem_meta().universal_part() ) );
  ASSERT_TRUE ( bulk.bucket(new_cell).member ( fixture.fem_meta().locally_owned_part() ) );
}

TEST ( UnitTestBulkData_new , verifyParallelAddParts )
{
  TestBoxFixture fixture;
  stk::unit_test_util::BulkDataTester &bulk = fixture.bulk_data ();
  ConstPartVector            add_part;

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  add_part.push_back ( &fixture.m_part_A_0 );

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  bulk.modification_begin();

  for ( EntityCommListInfoVector::const_iterator
        i =  bulk.my_internal_comm_list().begin();
        i != bulk.my_internal_comm_list().end() ; ++i ) {
    if ( i->key.rank() == 0 ) {
      if ( bulk.parallel_owner_rank(i->entity) == fixture.comm_rank() ) {
        bulk.change_entity_parts ( i->entity, add_part, ConstPartVector() );
      }
    }
  }

  bulk.modification_end();

  for ( EntityCommListInfoVector::const_iterator
        i =  bulk.my_internal_comm_list().begin();
        i != bulk.my_internal_comm_list().end() ; ++i ) {
    if ( i->key.rank() == 0 ) {
      ASSERT_TRUE ( bulk.bucket(i->entity).member ( fixture.m_part_A_0 ) );
    }
  }
}

TEST ( UnitTestBulkData_new , verifyInducedMembership )
{
  TestBoxFixture fixture;
  BulkData             &bulk = fixture.bulk_data ();
  PartVector            create_node_parts , create_cell_parts , empty_parts;

  create_node_parts.push_back ( &fixture.m_part_A_0 );
  create_cell_parts.push_back ( &fixture.m_cell_part );

  bulk.modification_begin();

  Entity node0 = bulk.declare_node(2);
  Entity node = bulk.declare_node(1);
  Entity cell = bulk.declare_element(1, stk::mesh::ConstPartVector{&fixture.get_elem_part()});
  bulk.change_entity_parts ( node , create_node_parts , PartVector () );
  bulk.change_entity_parts ( cell , create_cell_parts , PartVector () );
  // Add node to cell part
  RelationIdentifier cell_node_rel_id = 0;
  bulk.declare_relation ( cell , node0 , cell_node_rel_id );
  cell_node_rel_id = 1;
  bulk.declare_relation ( cell , node , cell_node_rel_id );

  unsigned cell_num_nodes = fixture.get_elem_topology().num_nodes();
  for (unsigned i = 2; i < cell_num_nodes; ++i)
  {
    Entity another_node = bulk.declare_node(i+1);
    bulk.declare_relation(cell, another_node, i);
  }

  bulk.modification_end();

  ASSERT_TRUE ( bulk.bucket(node).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.destroy_relation ( cell , node, cell_node_rel_id );
  Entity another_node = bulk.declare_node(cell_num_nodes);
  bulk.declare_relation(cell, another_node, cell_node_rel_id);
  bulk.modification_end();

  ASSERT_TRUE ( !bulk.bucket(node).member ( fixture.m_cell_part ) );
}

TEST ( UnitTestBulkData_new , verifyCanRemoveFromSetWithDifferentRankSubset )
{
  TestBoxFixture fixture;
  BulkData           &bulk = fixture.bulk_data ();
  PartVector          add_parts , add_elem_parts, remove_parts, empty_parts;

  add_parts.push_back ( &fixture.m_part_B_3 );
  add_parts.push_back ( &fixture.m_part_A_superset );
  add_elem_parts = add_parts;
  Part &elem_part = fixture.get_elem_part();
  add_elem_parts.push_back(&elem_part);

  remove_parts.push_back ( &fixture.m_part_A_superset );

  bulk.modification_begin();

  Entity e = bulk.declare_element(fixture.comm_rank()+1, add_elem_parts);
  Entity n = bulk.declare_node(fixture.comm_rank()+1, add_parts);
  bulk.declare_relation(e, n, 0);
  unsigned elem_num_nodes = fixture.get_elem_topology().num_nodes();
  for (unsigned i = 1; i < elem_num_nodes; ++i)
  {
    Entity another_node = bulk.declare_node(i+1);
    bulk.declare_relation(e, another_node, i);
  }
  bulk.modification_end();

  bulk.modification_begin();
  bulk.change_entity_parts ( e , empty_parts , remove_parts );
  bulk.modification_end();

  ASSERT_TRUE ( bulk.bucket(e).member ( fixture.m_part_B_3 ) );
  ASSERT_TRUE ( !bulk.bucket(e).member ( fixture.m_part_A_superset ) );
}


TEST ( UnitTestBulkData_new , verifyCommonGhostingName )
{
#ifndef NDEBUG
  //the throw for inconsistent name only happens in debug mode.
  TestBoxFixture fixture;
  BulkData          &bulk = fixture.bulk_data ();

  bulk.modification_begin();

  if ( fixture.comm_size() == 1 ) return;

  if ( fixture.comm_rank() == 0 )
  {
    ASSERT_THROW ( bulk.create_ghosting ( "Name 1" ) , std::runtime_error );
  }
  else
  {
    ASSERT_THROW ( bulk.create_ghosting ( "Name 2" ) , std::runtime_error );
  }
#endif
}


TEST ( UnitTestBulkData_new , verifyTrivialDestroyAllGhostings )
{
  TestBoxFixture fixture;

  if ( fixture.comm_size() == 1 ) return;

  BulkData  &bulk = fixture.bulk_data();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  bulk.modification_begin();

  Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );

  // Find a cell owned by this process
  BucketVector::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;

  std::vector<EntityProc>  to_send;
  std::vector<EntityKey>   empty_vector;
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) == fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        send_rank++;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }
  bulk.change_ghosting ( ghosting , to_send , empty_vector );
  bulk.modification_end();


  {
    std::vector<EntityProc> send_list ;
    std::vector<EntityKey>  recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    ASSERT_TRUE ( ! send_list.empty()  );
    ASSERT_TRUE ( ! recv_list.empty() );
  }

  // Usage of operator << in Ghosting.cpp
  std::ostringstream oss;
  oss << ghosting;

  bulk.modification_begin();
  bulk.destroy_all_ghosting ();
  bulk.modification_end();

  {
    std::vector<EntityProc> send_list ;
    std::vector<EntityKey>  recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    ASSERT_TRUE ( send_list.empty() );
    ASSERT_TRUE ( recv_list.empty() );
  }
}


TEST ( UnitTestBulkData_new , verifyChangeGhostingGuards )
{
  TestBoxFixture fixture1, fixture2;
  BulkData & bulk1 = fixture1.bulk_data ();
  BulkData & bulk2 = fixture2.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box1[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };
  int local_box2[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk1.modification_begin();
  fixture1.generate_boxes( root_box, local_box1 );
  ASSERT_TRUE(bulk1.modification_end());

  bulk2.modification_begin();
  fixture2.generate_boxes( root_box, local_box2 );
  ASSERT_TRUE(bulk2.modification_end());

  bulk1.modification_begin();
  bulk2.modification_begin();

  std::vector<EntityProc>  to_send;
  std::vector<EntityKey>   empty_vector;
  BucketVector::const_iterator cur_bucket = bulk1.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;
  while ( cur_bucket != bulk1.buckets(stk::topology::ELEM_RANK).end() )
  {
    Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk1.parallel_owner_rank(*cur_entity) == fixture1.comm_rank() )
      {
        if ( send_rank == fixture1.comm_size() ) send_rank = 0;
        if ( send_rank != fixture1.comm_rank() )
          to_send.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        ++send_rank;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  Ghosting &ghosting = bulk1.create_ghosting ( "Ghost 1" );
  ASSERT_THROW ( bulk1.change_ghosting ( bulk1.aura_ghosting() , to_send , empty_vector ) , std::runtime_error );

  ghosting.receive_list(empty_vector);
  ghosting.send_list(to_send);

  bulk1.modification_end();
  bulk2.modification_end();
}


TEST ( UnitTestBulkData_new , verifyOtherGhostingGuards )
{
  TestBoxFixture fixture;
  BulkData          &bulk = fixture.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  bulk.modification_begin();

  std::vector<EntityProc>  to_send_unowned;
  std::vector<EntityProc>  empty_send;
  std::vector<EntityKey>   to_remove_not_ghosted;
  std::vector<EntityKey>   empty_remove;
  BucketVector::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) != fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send_unowned.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        ++send_rank;
      }
      else
      {
        to_remove_not_ghosted.push_back ( bulk.entity_key(*cur_entity) );
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );
  if ( to_send_unowned.size() > 0 )
  {
    ASSERT_THROW ( bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove );
  }

  if ( to_remove_not_ghosted.size() > 0 )
  {
    ASSERT_NO_THROW( bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted ) );
  }
  else
  {
    bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted );
  }
  bulk.modification_end();
}


TEST ( UnitTestBulkData_new , verifyPartsOnCreate )
{
  TestBoxFixture fixture;
  BulkData           & bulk = fixture.bulk_data ();
  Part               & part_a = fixture.m_part_A_0;
  Part               & part_b = fixture.m_part_B_0;

  PartVector           create_vector;
  create_vector.push_back ( &part_a );

  bulk.modification_begin();

  Entity node = bulk.declare_node(fixture.comm_rank()+1 ,create_vector );
  bulk.modification_end();

  ASSERT_TRUE ( bulk.bucket(node).member ( part_a ) );

  bulk.modification_begin();
  create_vector.push_back ( &part_b );
  Entity node2 = bulk.declare_node(fixture.comm_size() + fixture.comm_rank() + 1 , create_vector );
  bulk.modification_end();

  ASSERT_TRUE ( bulk.bucket(node2).member ( part_a ) );
  ASSERT_TRUE ( bulk.bucket(node2).member ( part_b ) );
}

//----------------------------------------------------------------------

TEST ( UnitTestBulkData_new , verifyBoxGhosting )
{
  const int p_size = stk::parallel_machine_size( MPI_COMM_WORLD );
  if ( 8 < p_size ) { return ; }

  fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 2, 2 );
  fixture.m_meta.commit();
  fixture.generate_mesh();
  const BulkData& mesh = fixture.m_bulk_data;

  for ( size_t iz = 0 ; iz < 3 ; ++iz ) {
    for ( size_t iy = 0 ; iy < 3 ; ++iy ) {
      for ( size_t ix = 0 ; ix < 3 ; ++ix ) {
        Entity const node = fixture.node(ix,iy,iz);
        ASSERT_TRUE( mesh.is_valid(node) );

        ASSERT_TRUE( fixture.node_id(ix,iy,iz) == mesh.identifier(node) );
        fixtures::HexFixture::Scalar * const node_coord = stk::mesh::field_data(*fixture.m_coord_field, node);
        ASSERT_TRUE( node_coord != NULL );
      }
    }
  }

  for ( size_t iz = 0 ; iz < 2 ; ++iz ) {
    for ( size_t iy = 0 ; iy < 2 ; ++iy ) {
      for ( size_t ix = 0 ; ix < 2 ; ++ix ) {
        Entity const elem = fixture.elem(ix,iy,iz);
        ASSERT_TRUE( mesh.is_valid(elem) );
        size_t num_elem_nodes = mesh.num_nodes(elem);
        ASSERT_EQ( 8u , num_elem_nodes );
        Entity const *elem_nodes = mesh.begin_nodes(elem);
        // ConnectivityOrdinal const *elem_node_ords = mesh.begin_node_ordinals(elem);
        if ( 8u == num_elem_nodes ) {
          ASSERT_TRUE( elem_nodes[0] == fixture.node(ix,iy,iz));
          ASSERT_TRUE( elem_nodes[1] == fixture.node(ix+1,iy,iz));
          ASSERT_TRUE( elem_nodes[2] == fixture.node(ix+1,iy+1,iz));
          ASSERT_TRUE( elem_nodes[3] == fixture.node(ix,iy+1,iz));
          ASSERT_TRUE( elem_nodes[4] == fixture.node(ix,iy,iz+1));
          ASSERT_TRUE( elem_nodes[5] == fixture.node(ix+1,iy,iz+1));
          ASSERT_TRUE( elem_nodes[6] == fixture.node(ix+1,iy+1,iz+1));
          ASSERT_TRUE( elem_nodes[7] == fixture.node(ix,iy+1,iz+1));
        }
        // Now check access to field data via the fast rank functions.
        // Node const *eph_elem_nodes = mesh.begin_nodes(elem);
        // for ( size_t j = 0 ; j < num_elem_nodes ; ++j )
        // {
        //   fixtures::HexFixture::Scalar * const node_coord =
        //     stk::mesh::field_data( fixture.m_coord_field , eph_elem_nodes[j]);
        //   EXPECT_EQ( node_coord, elem_node_coord[ elem_node_ords[j] ] );
        // }

      }
    }
  }
}

TEST ( UnitTestBulkData_new , testUninitializedMetaData )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  std::shared_ptr<BulkData> bulk = stk::mesh::MeshBuilder(pm).create();
  MetaData& meta = bulk->mesh_meta_data();

  meta.initialize(2);

  meta.commit();

  bulk->modification_begin();

  ASSERT_THROW( bulk->declare_node(1, PartVector()), std::logic_error);
}

TEST ( UnitTestBulkData_new , testGhostHandleRemainsValidAfterRefresh )
{
  //
  // Testing if a handle to a ghosted entity remains valid before and after a
  // modification cycle in which the ghost is refreshed.
  //
  // To test this, we focus on a single node shared on 2 procs, ghosted on others
  //
  //
  // 1D Mesh (node,owner)--[elem,owner]---(...)
  //
  // <---(50,0)--[100,0]--(21,1)--[201,1]---(32,2)---[302,2]---(50,0)--->
  //

  // elem, node0, node1, owner
  EntityId elems_0[][4] = { {100, 21, 50, 0}, {201, 21, 32, 1}, {302, 32, 50, 2} };
  // node, owner
  EntityId nodes_0[][2] = { {21,1}, {50,0}, {32, 2} };

  const unsigned nelems = sizeof(elems_0)/4/sizeof(EntityId);
  const unsigned nnodes = sizeof(nodes_0)/2/sizeof(EntityId);

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 1;

  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("NODE_RANK");
  entity_rank_names.push_back("EDGE_RANK");
  entity_rank_names.push_back("FACE_RANK");
  entity_rank_names.push_back("ELEM_RANK");
  entity_rank_names.push_back("FAMILY_TREE");

  stk::mesh::MeshBuilder builder(pm);
  builder.set_spatial_dimension(spatial_dim);
  builder.set_entity_rank_names(entity_rank_names);
  std::shared_ptr<BulkData> bulkPtr = builder.create();
  MetaData& meta_data = bulkPtr->mesh_meta_data();
  Part & elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::LINE_2_1D);
  Part & node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);

  meta_data.commit();
  BulkData& mesh = *bulkPtr;
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if (p_size != 3) return;

  // Build map for node sharing
  stk::mesh::fixtures::NodeToProcsMMap nodes_to_procs;
  {
    for (unsigned ielem=0; ielem < nelems; ielem++) {
      int e_owner = static_cast<int>(elems_0[ielem][3]);
      stk::mesh::fixtures::AddToNodeProcsMMap(nodes_to_procs, elems_0[ielem][2], e_owner);
      stk::mesh::fixtures::AddToNodeProcsMMap(nodes_to_procs, elems_0[ielem][1], e_owner);
    }
  }

  //
  // Begin modification cycle so we can create the entities and relations
  //
  {
    // Create elements
    Entity elem = Entity();

    mesh.modification_begin();

    for (unsigned ielem=0; ielem < nelems; ielem++) {
      int owner = static_cast<int>(elems_0[ielem][3]);
      if (owner == p_rank) {
        elem = mesh.declare_element(elems_0[ielem][0], stk::mesh::ConstPartVector{&elem_part});

        EntityVector nodes;
        // Create node on all procs
        nodes.push_back( mesh.declare_node(elems_0[ielem][2], stk::mesh::ConstPartVector{&node_part}) );
        nodes.push_back( mesh.declare_node(elems_0[ielem][1], stk::mesh::ConstPartVector{&node_part}) );

        // Add relations to nodes
        mesh.declare_relation( elem, nodes[0], 0 );
        mesh.declare_relation( elem, nodes[1], 1 );

        // Node sharing
        stk::mesh::fixtures::DoAddNodeSharings(mesh, nodes_to_procs, mesh.identifier(nodes[0]), nodes[0]);
        stk::mesh::fixtures::DoAddNodeSharings(mesh, nodes_to_procs, mesh.identifier(nodes[1]), nodes[1]);
      }
    }

    mesh.modification_end();
  }

  // change node owners
  {
    std::vector<EntityProc> change;

    for (unsigned inode=0; inode < nnodes; inode++) {
      Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodes_0[inode][0]);
      if (mesh.is_valid(node) && mesh.parallel_owner_rank(node) == p_rank) {
        int dest = nodes_0[inode][1];
        EntityProc eproc(node, dest);
        change.push_back(eproc);
      }
    }

    mesh.change_entity_owner( change );
  }

  // The real test is here
  {
    Entity node_21_handle_before_elem_deletion = mesh.get_entity(stk::topology::NODE_RANK, 21);
    if (p_rank == 2) {
      ASSERT_TRUE(mesh.in_receive_ghost(mesh.entity_key(node_21_handle_before_elem_deletion)));
    }

    // attempt to delete a node and its elems but on a ghosted proc
    mesh.modification_begin();

    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 100);
    if (mesh.is_valid(elem)) mesh.destroy_entity(elem);

    mesh.modification_end();

    // Key check is here
    if (p_rank == 2) {
      // mesh still thinks handle is valid and refers to the same entity with the same EntityKey!
      ASSERT_TRUE(mesh.is_valid(node_21_handle_before_elem_deletion));
      ASSERT_EQ(mesh.entity_key(node_21_handle_before_elem_deletion), EntityKey(stk::topology::NODE_RANK, 21));
    }
  }
}

TEST ( UnitTestBulkData_new , testCustomBucketCapacity )
{
  const int spatial_dimension = 3;

  const unsigned non_standard_bucket_capacity = 42;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatial_dimension);
  builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
  builder.set_add_fmwk_data(true);
  builder.set_bucket_capacity(non_standard_bucket_capacity);
  std::shared_ptr<BulkData> bulk = builder.create();
  MetaData& meta = bulk->mesh_meta_data();


  Part & node_part = meta.declare_part_with_topology("node_part", stk::topology::NODE);

  meta.commit();

  PartVector    create_vector;
  create_vector.push_back ( &node_part );

  bulk->modification_begin();

  EntityId nodeID = bulk->parallel_rank()+1;
  Entity node = bulk->declare_node(nodeID, create_vector);
  bulk->modification_end();

  EXPECT_EQ( bulk->bucket(node).capacity(), non_standard_bucket_capacity );
}

TEST(BulkData, bucket_debug_ThrowsWithInvalidEntity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

#ifndef NDEBUG
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  stk::mesh::Entity invalidEntity;
  EXPECT_TRUE(nullptr == bulkPtr->bucket_ptr(invalidEntity));
  EXPECT_ANY_THROW(bulkPtr->bucket(invalidEntity));
#endif
}

