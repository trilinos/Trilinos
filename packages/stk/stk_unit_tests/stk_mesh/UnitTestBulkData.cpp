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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>  // for test_shift_ring
#include <unit_tests/Setup8Quad4ProcMesh.hpp>
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>

namespace stk { namespace mesh { class FieldBase; } }

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::BaseEntityRank;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace {

const EntityRank NODE_RANK = stk::topology::NODE_RANK;

void donate_one_element( BulkData & mesh , bool aura )
{
  const int p_rank = mesh.parallel_rank();

  Selector select_owned( MetaData::get(mesh).locally_owned_part() );

  std::vector<unsigned> before_count ;
  std::vector<unsigned> after_count ;

  count_entities( select_owned , mesh , before_count );

  // Change owner of an element on a process boundary
  // from P0 to P1, and then recount to confirm ownership change

  std::vector<EntityProc> change ;

  // A shared node:
  EntityKey node_key;
  Entity elem = Entity();

  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    if ( mesh.in_shared( i->key ) && i->key.rank() == BaseEntityRank ) {
      node_key = i->key;
      break;
    }
  }

  ASSERT_TRUE( node_key.is_valid() );

  Entity node = mesh.get_entity(node_key);
  ASSERT_TRUE( mesh.is_valid(node) );

  Entity const *node_elems_i = mesh.begin_elements(node);
  Entity const *node_elems_e = mesh.end_elements(node);
  for ( ; (node_elems_i != node_elems_e) && !mesh.is_valid(elem); ++node_elems_i)
  {
    elem = *node_elems_i;
    if ( mesh.parallel_owner_rank(elem) != p_rank ) { elem = Entity(); }
  }

  ASSERT_TRUE( mesh.is_valid(elem) );

  unsigned donated_nodes = 0 ;

  // Only process #0 donates an element and its owned nodes:
  if ( 0 == p_rank ) {
    EntityProc entry ;
    entry.first = elem ;
    entry.second = mesh.entity_comm_map_shared(mesh.entity_key(node))[0].proc;
    change.push_back( entry );

    Entity const *elem_nodes_i = mesh.begin_nodes(elem);
    Entity const *elem_nodes_e = mesh.end_nodes(elem);
    for ( ; elem_nodes_i != elem_nodes_e; ++elem_nodes_i)
    {
      if ( mesh.parallel_owner_rank(*elem_nodes_i) == p_rank ) {
        entry.first = *elem_nodes_i;
        change.push_back( entry );
        ++donated_nodes ;
      }
    }
  }

  mesh.change_entity_owner( change, aura, BulkData::MOD_END_COMPRESS_AND_SORT );

  count_entities( select_owned , mesh , after_count );

  if ( 0 == p_rank ) {
    ASSERT_EQ( before_count[3] - 1 , after_count[3] );
    ASSERT_EQ( before_count[0] - donated_nodes, after_count[0] );
  }
}

void donate_all_shared_nodes( BulkData & mesh , bool aura )
{
  const int p_rank = mesh.parallel_rank();

  const Selector select_used = MetaData::get(mesh).locally_owned_part() |
                               MetaData::get(mesh).globally_shared_part() ;

  std::vector<unsigned> before_count ;
  std::vector<unsigned> after_count ;

  count_entities( select_used , mesh , before_count );

  // Donate owned shared nodes to first sharing process.

  const stk::mesh::EntityCommListInfoVector & entity_comm = mesh.comm_list();

  ASSERT_TRUE( ! entity_comm.empty() );

  std::vector<EntityProc> change ;

  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i =  entity_comm.begin() ;
        i != entity_comm.end() &&
        i->key.rank() == BaseEntityRank ; ++i ) {
    Entity const node = i->entity;
    const stk::mesh::PairIterEntityComm ec = mesh.entity_comm_map_shared(i->key);

    if ( mesh.parallel_owner_rank(node) == p_rank && ! ec.empty() ) {
      change.push_back( EntityProc( node , ec->proc ) );
    }
  }

  mesh.change_entity_owner( change, aura, BulkData::MOD_END_COMPRESS_AND_SORT );

  count_entities( select_used , mesh , after_count );

  ASSERT_TRUE( 3 <= after_count.size() );
  ASSERT_EQ( before_count[0] , after_count[0] );
  ASSERT_EQ( before_count[1] , after_count[1] );
  ASSERT_EQ( before_count[2] , after_count[2] );
  ASSERT_EQ( before_count[3] , after_count[3] );
}


//----------------------------------------------------------------------
// Testing for mesh entities without relations

TEST(UnitTestingOfBulkData, testChangeOwner_nodes)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );
  BulkData bulk( meta , pm , 100 );

  const PartVector no_parts ;

  meta.commit();
  bulk.modification_begin();

  // Ids for all entities (all entities have type 0):

  std::vector<EntityId> ids( id_total );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    ids[i] = i + 1;
  }

  // Declare just those entities in my range of ids:

  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    bulk.declare_entity( stk::topology::NODE_RANK , ids[i] , no_parts );
  }

  ASSERT_TRUE( bulk.modification_end() );

  // Verify that I only have entities in my range:

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity e = bulk.get_entity( stk::topology::NODE_RANK , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      ASSERT_TRUE( bulk.is_valid(e) );
    }
    else {
      ASSERT_TRUE( !bulk.is_valid(e) );
    }
  }

  // Test change owner no-op first:

  std::vector<EntityProc> change ;

  bulk.change_entity_owner( change );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity e = bulk.get_entity( stk::topology::NODE_RANK , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      ASSERT_TRUE( bulk.is_valid(e) );
    }
    else {
      ASSERT_TRUE( !bulk.is_valid(e) );
    }
  }

  // Can only test changing owner in parallel.

  if ( 1 < p_size ) {
    // Give my last two ids to the next process
    // Get the previous process' last two ids

    const int p_give = ( p_rank + 1 ) % p_size ;
    const unsigned id_give = id_end - 2 ;
    const unsigned id_get  = ( id_begin + id_total - 2 ) % id_total ;

    ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_give] )) );
    ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_give+1] )) );
    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get] )) );
    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get+1] )) );

    change.resize(2);
    change[0].first = bulk.get_entity( stk::topology::NODE_RANK , ids[id_give] );
    change[0].second = p_give ;
    change[1].first = bulk.get_entity( stk::topology::NODE_RANK , ids[id_give+1] );
    change[1].second = p_give ;

    bulk.change_entity_owner( change );

    ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get] )) );
    ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get+1] )) );

    // Entities given away are destroyed until the next modification cycle
    {
      Entity const e0 = bulk.get_entity( stk::topology::NODE_RANK , ids[id_give] );
      Entity const e1 = bulk.get_entity( stk::topology::NODE_RANK , ids[id_give+1] );
      ASSERT_TRUE( !bulk.is_valid(e0) );
      ASSERT_TRUE( !bulk.is_valid(e1) );
    }

    ASSERT_TRUE( bulk.modification_begin() );
    ASSERT_TRUE( bulk.modification_end() );

    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_give] )) );
    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_give+1] )) );
  }
}

//----------------------------------------------------------------------
// Testing for creating existing mesh entities without relations

TEST(UnitTestingOfBulkData, testCreateMore)
{
  std::cout << "testCreateMore needs to be replaced with corresponding test that makes use of node sharing API enforcement" << std::endl;
  /*
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };

  const int p_size = stk::parallel_machine_size( pm );
  const int p_rank = stk::parallel_machine_rank( pm );

  if ( 1 < p_size ) {

    const unsigned id_total = nPerProc * p_size ;
    const unsigned id_begin = nPerProc * p_rank ;
    const unsigned id_end   = nPerProc * ( p_rank + 1 );

    const int spatial_dimension = 3;
    MetaData meta( spatial_dimension );

    const PartVector no_parts ;

    meta.commit();

    BulkData bulk( meta , pm , 100 );

    bulk.modification_begin();

    // Ids for all entities (all entities have type 0):

    std::vector<EntityId> ids( id_total );

    for ( unsigned i = 0 ; i < id_total ; ++i ) { ids[i] = i + 1; }

    // Declare just those entities in my range of ids:

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      bulk.declare_entity( stk::topology::NODE_RANK , ids[i] , no_parts );
    }

    ASSERT_TRUE( bulk.modification_end() );

    // Only one process create entities with previous process' last two ids

    const unsigned id_get  = ( id_begin + id_total - 2 ) % id_total ;

    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get] )) );
    ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get+1] )) );

    ASSERT_TRUE( bulk.modification_begin() );

    if ( 1 == p_rank ) {
      // These declarations create entities that already exist,
      // which will be an error.  Must create an owned entity
      // to use them, thus they become shared.

      Entity e0 = bulk.declare_entity( stk::topology::NODE_RANK , ids[ id_get ] , no_parts );
      Entity e1 = bulk.declare_entity( stk::topology::NODE_RANK , ids[ id_get + 1 ] , no_parts );

      Entity eU = bulk.declare_entity( stk::topology::EDGE_RANK , 1 , no_parts );

      bulk.declare_relation( eU , e0 , 0 );
      bulk.declare_relation( eU , e1 , 1 );
    }

    bulk.modification_end();

    if ( 1 == p_rank ) {
      Entity e0 = bulk.get_entity( stk::topology::NODE_RANK , ids[id_get] );
      Entity e1 = bulk.get_entity( stk::topology::NODE_RANK , ids[id_get+1] );
      ASSERT_TRUE( bulk.is_valid(e0) );
      ASSERT_TRUE( bulk.is_valid(e1) );
      ASSERT_TRUE( 0 == bulk.parallel_owner_rank(e0) );
      ASSERT_TRUE( 0 == bulk.parallel_owner_rank(e1) );
    }

    // Now test tripping the error condition

    bulk.modification_begin();

    if ( 0 == p_rank ) {
      bulk.declare_entity( stk::topology::NODE_RANK , ids[ id_get ] , no_parts );
      bulk.declare_entity( stk::topology::NODE_RANK , ids[ id_get + 1 ] , no_parts );
    }

    ASSERT_THROW( bulk.modification_end() , std::runtime_error );
  }
  */
}

//----------------------------------------------------------------------
TEST(UnitTestingOfBulkData, testBulkDataRankBeginEnd)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const int p_size = stk::parallel_machine_size( pm );
  if (p_size != 1) {
    return;
  }

  const size_t spatial_dim = 3;
  MetaData meta(spatial_dim, stk::mesh::entity_rank_names());
  BulkData bulk(meta, pm);
  bulk.modification_begin();
  BulkData::const_entity_iterator iter = bulk.begin_entities(stk::topology::NODE_RANK);
  BulkData::const_entity_iterator end = bulk.end_entities(stk::topology::NODE_RANK);

  ASSERT_TRUE(iter == end);

  EntityId node_id = 1;
  bulk.declare_entity(stk::topology::NODE_RANK, node_id);

  iter = bulk.begin_entities(stk::topology::NODE_RANK);
  end = bulk.end_entities(stk::topology::NODE_RANK);

  //insist that there is 1 node:
  ASSERT_TRUE(iter != end);
  ASSERT_TRUE(std::distance(iter,end) == 1u);

  //now declare an edge...
  EntityId edge_id = 1;
  bulk.declare_entity(stk::topology::EDGE_RANK, edge_id);

  iter = bulk.begin_entities(stk::topology::NODE_RANK);
  end = bulk.end_entities(stk::topology::NODE_RANK);

  //insist that there is still 1 node:
  ASSERT_TRUE(iter != end);
  ASSERT_TRUE(std::distance(iter,end) == 1u);

  iter = bulk.begin_entities(stk::topology::EDGE_RANK);
  end = bulk.end_entities(stk::topology::EDGE_RANK);

  //insist that there is 1 edge:
  ASSERT_TRUE(iter != end);
  ASSERT_TRUE(std::distance(iter,end) == 1u);

  node_id = 2;
  bulk.declare_entity(stk::topology::NODE_RANK, node_id);

  iter = bulk.begin_entities(stk::topology::NODE_RANK);
  end = bulk.end_entities(stk::topology::NODE_RANK);

  //insist that there are 2 nodes:
  ASSERT_TRUE(iter != end);
  ASSERT_TRUE(std::distance(iter,end) == 2u);

  iter = bulk.begin_entities(stk::topology::EDGE_RANK);
  end = bulk.end_entities(stk::topology::EDGE_RANK);

  //insist that there is still 1 edge:
  ASSERT_TRUE(iter != end);
  ASSERT_TRUE(std::distance(iter,end) == 1u);

  iter = bulk.begin_entities(stk::topology::FACE_RANK);
  end = bulk.end_entities(stk::topology::FACE_RANK);

  //insist that there are no faces:
  ASSERT_TRUE(iter == end);
}
//----------------------------------------------------------------------

TEST(UnitTestingOfBulkData, testChangeOwner_ring)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;

  std::vector<unsigned> local_count ;

  //------------------------------
  {
    bool aura = false;
    RingFixture ring_mesh( pm , nPerProc , false /* no element parts */ );
    BulkData & bulk = ring_mesh.m_bulk_data;
    ring_mesh.m_meta_data.commit();

    bulk.modification_begin();
    ring_mesh.generate_mesh( );
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, aura));

    ring_mesh.fixup_node_ownership(aura, BulkData::MOD_END_COMPRESS_AND_SORT);

    const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() |
                                 ring_mesh.m_meta_data.globally_shared_part() ;
    const Selector select_all = ring_mesh.m_meta_data.universal_part() ;

    count_entities( select_used , ring_mesh.m_bulk_data , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

    count_entities( select_all , ring_mesh.m_bulk_data , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

    if ( 1 < p_size ) {
      // Shift ring by two nodes and elements.

      stk::unit_test::test_shift_ring( ring_mesh, false /* no aura */ );

      count_entities( select_used , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement );

      count_entities( select_all , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement );
    }
  }

  //------------------------------
  // Test shift starting with ghosting but not regenerated ghosting.
  {
    RingFixture ring_mesh( pm , nPerProc , false /* no element parts */ );
    BulkData& bulk = ring_mesh.m_bulk_data;
    ring_mesh.m_meta_data.commit();

    bulk.modification_begin();
    ring_mesh.generate_mesh( );
    ASSERT_TRUE(bulk.modification_end());

    ring_mesh.fixup_node_ownership( );

    const Selector select_owned( ring_mesh.m_meta_data.locally_owned_part() );
    const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() |
                                 ring_mesh.m_meta_data.globally_shared_part() ;
    const Selector select_all(   ring_mesh.m_meta_data.universal_part() );

    count_entities( select_used , ring_mesh.m_bulk_data , local_count );
    ASSERT_EQ( local_count[stk::topology::NODE_RANK] , nLocalNode );
    ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK] , nLocalElement );

    count_entities( select_all , ring_mesh.m_bulk_data , local_count );
    const unsigned n_extra = 1 < p_size ? 2 : 0 ;
    ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode + n_extra );
    ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement + n_extra );

    if ( 1 < p_size ) {
      stk::unit_test::test_shift_ring( ring_mesh, false /* no aura */ );

      count_entities( select_owned , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nPerProc );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nPerProc );

      count_entities( select_used , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement );

      // All of my ghosts were disrupted and therefore deleted:
      count_entities( select_all , ring_mesh.m_bulk_data , local_count );
      ASSERT_EQ( nLocalElement , local_count[stk::topology::ELEMENT_RANK] );
      ASSERT_EQ( nLocalNode , local_count[stk::topology::NODE_RANK] );
    }
  }
  //------------------------------
  // Test shift starting with ghosting and regenerating ghosting.
  {
    RingFixture ring_mesh( pm , nPerProc , false /* no element parts */ );
    BulkData& bulk = ring_mesh.m_bulk_data;
    ring_mesh.m_meta_data.commit();

    bulk.modification_begin();
    ring_mesh.generate_mesh( );
    ASSERT_TRUE(bulk.modification_end());

    ring_mesh.fixup_node_ownership( );

    const Selector select_owned( ring_mesh.m_meta_data.locally_owned_part() );
    const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() |
                                 ring_mesh.m_meta_data.globally_shared_part() ;
    const Selector select_all(   ring_mesh.m_meta_data.universal_part() );

    count_entities( select_used , ring_mesh.m_bulk_data , local_count );
    ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode );
    ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement );

    count_entities( select_all , ring_mesh.m_bulk_data , local_count );
    const unsigned n_extra = 1 < p_size ? 2 : 0 ;
    ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode + n_extra );
    ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement + n_extra );

    if ( 1 < p_size ) {
      stk::unit_test::test_shift_ring( ring_mesh, true /* with aura */ );

      count_entities( select_owned , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nPerProc );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nPerProc );

      count_entities( select_used , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement );

      // All of my ghosts were regenerated:
      count_entities( select_all , ring_mesh.m_bulk_data , local_count );
      ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode + n_extra );
      ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement + n_extra );
    }
  }
  //------------------------------
  // Test bad owner change catching:
  if ( 1 < p_size ) {
    RingFixture ring_mesh( pm , nPerProc , false /* no element parts */ );
    BulkData& bulk = ring_mesh.m_bulk_data;
    ring_mesh.m_meta_data.commit();

    bulk.modification_begin();
    ring_mesh.generate_mesh( );
    ASSERT_TRUE(bulk.modification_end());

    ring_mesh.fixup_node_ownership( );

    std::vector<EntityProc> change ;

    if ( 0 == p_rank ) {
      change.resize(4);
      // Error to change to bad owner:
      change[0].first = ring_mesh.m_bulk_data.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[1] );
      change[0].second = p_size ;
      // Error to change a ghost:
      for ( stk::mesh::EntityCommListInfoVector::const_iterator
            ec =  ring_mesh.m_bulk_data.comm_list().begin() ;
            ec != ring_mesh.m_bulk_data.comm_list().end() ; ++ec ) {
        if ( bulk.in_receive_ghost( ec->key ) ) {
          change[1].first = ec->entity;
          break ;
        }
      }
      change[1].second = p_rank ;
      // Error to change to multiple owners:
      change[2].first = ring_mesh.m_bulk_data.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[1] );
      change[2].second = ( p_rank + 1 ) % p_size ;
      change[3].first = change[2].first ;
      change[3].second = ( p_rank + 2 ) % p_size ;
    }

    ASSERT_THROW( ring_mesh.m_bulk_data.change_entity_owner( change ),
                          std::runtime_error );
  }
  //------------------------------
  // Test move one element with initial ghosting but not regenerated ghosting:
  // last processor give its shared node to P0
  if ( 1 < p_size ) {
    RingFixture ring_mesh( pm , nPerProc , false /* no element parts */ );
    BulkData& bulk = ring_mesh.m_bulk_data;
    ring_mesh.m_meta_data.commit();

    bulk.modification_begin();
    ring_mesh.generate_mesh( );
    ASSERT_TRUE(bulk.modification_end());

    ring_mesh.fixup_node_ownership( );

    const Selector select_owned( ring_mesh.m_meta_data.locally_owned_part() );
    const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() |
                                 ring_mesh.m_meta_data.globally_shared_part() ;
    const Selector select_all(   ring_mesh.m_meta_data.universal_part() );

    std::vector<EntityProc> change ;

    if ( p_rank + 1 == p_size ) {
      EntityProc entry ;
      entry.first = ring_mesh.m_bulk_data.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[0] );
      entry.second = 0 ;
      ASSERT_EQ( p_rank , bulk.parallel_owner_rank(entry.first) );
      change.push_back( entry );
    }

    ring_mesh.m_bulk_data.change_entity_owner( change, false /* regenerate_aura */, BulkData::MOD_END_COMPRESS_AND_SORT );

    count_entities( select_owned , ring_mesh.m_bulk_data , local_count );
    const unsigned n_node = p_rank == 0          ? nPerProc + 1 : (
                            p_rank + 1 == p_size ? nPerProc - 1 :
                                                   nPerProc );

    ASSERT_EQ( n_node , local_count[stk::topology::NODE_RANK] );
    ASSERT_EQ( static_cast<unsigned>(nPerProc) , local_count[stk::topology::ELEMENT_RANK] );

    count_entities( select_used , ring_mesh.m_bulk_data , local_count );
    ASSERT_EQ( nLocalNode , local_count[stk::topology::NODE_RANK] );
    ASSERT_EQ( nLocalElement , local_count[stk::topology::ELEMENT_RANK] );

    // Moving the node disrupted ghosting on first and last process
    count_entities( select_all , ring_mesh.m_bulk_data , local_count );
    const unsigned n_extra = p_rank + 1 == p_size || p_rank == 0 ? 1 : 2 ;
    ASSERT_EQ( nLocalNode + n_extra , local_count[stk::topology::NODE_RANK] );
    ASSERT_EQ( nLocalElement + n_extra , local_count[stk::topology::ELEMENT_RANK] );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing for collection of boxes

TEST(UnitTestingOfBulkData, testChangeOwner_box)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };

  const int p_size = stk::parallel_machine_size( pm );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );

  meta.commit();

  //------------------------------
  {
    bool aura = false;
    BoxFixture fixture( pm, 100 );
    fixture.fem_meta().commit();
    BulkData & bulk = fixture.bulk_data();
    int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    bulk.modification_begin();
    fixture.generate_boxes( root_box, local_box );
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, aura));

    if ( 1 < p_size ) {
      donate_one_element( bulk , aura );
    }
  }

  if ( 1 < p_size ) {
    bool aura = false;
    BoxFixture fixture( pm, 100 );
    fixture.fem_meta().commit();
    BulkData & bulk = fixture.bulk_data();
    int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    bulk.modification_begin();
    fixture.generate_boxes( root_box, local_box );
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, aura));

    donate_all_shared_nodes( bulk , aura );
  }
  //------------------------------
  if ( 1 < p_size ) {
    bool aura = false;
    BoxFixture fixture( pm, 100 );
    fixture.fem_meta().commit();
    BulkData & bulk = fixture.bulk_data();
    int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    bulk.modification_begin();
    fixture.generate_boxes( root_box, local_box );
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, aura));

    donate_one_element( bulk , false /* no aura */ );
  }
  //------------------------------
  // Introduce ghosts:
  if ( 1 < p_size ) {
    BoxFixture fixture( pm, 100 );
    BulkData & bulk = fixture.bulk_data();
    MetaData & box_meta = fixture.fem_meta();
    box_meta.commit();
    int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    bulk.modification_begin();
    fixture.generate_boxes( root_box, local_box );
    ASSERT_TRUE(bulk.modification_end());

    std::vector<unsigned> used_count ;
    std::vector<unsigned> all_count ;

    const Selector select_owned( box_meta.locally_owned_part() );
    const Selector select_used = box_meta.locally_owned_part() |
                                 box_meta.globally_shared_part() ;
    const Selector select_all(   box_meta.universal_part() );

    count_entities( select_all , bulk , all_count );
    count_entities( select_used , bulk , used_count );

    ASSERT_TRUE( used_count[0] < all_count[0] );
    ASSERT_TRUE( used_count[3] < all_count[3] );

    donate_all_shared_nodes( bulk , false /* don't regenerate aura */ );

    count_entities( select_all , bulk , all_count );
    count_entities( select_used , bulk , used_count );

    ASSERT_EQ( used_count[0] , all_count[0] );
    ASSERT_EQ( used_count[3] , all_count[3] );
  }
}

TEST(UnitTestingOfBulkData, testModifyPropagation)
{
  // Our new modification model makes it so the modified status
  // of an entity is propagated up to higher-ranked entities
  // that have relations to the modified entity. We test this
  // by grabbing a node off of a ring mesh, modifying it, and
  // checking that its element also gets marked as modified.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const unsigned nPerProc = 2;
  const int p_size = stk::parallel_machine_size( pm );

  // this test only needs to be run w/ one processor
  if (p_size > 1) return;

  // Make a ring_mesh and add an extra part
  RingFixture ring_mesh( pm , nPerProc, false /* don't use element parts */);
  stk::mesh::Part& special_part =
    ring_mesh.m_meta_data.declare_part("special_node_part", stk::mesh::BaseEntityRank );
  ring_mesh.m_meta_data.commit();
  BulkData& bulk = ring_mesh.m_bulk_data;

  bulk.modification_begin();
  ring_mesh.generate_mesh( );
  ASSERT_TRUE(bulk.modification_end());

  ring_mesh.fixup_node_ownership( );

  // grab the first element
  EntityVector elements;
  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;
  stk::mesh::get_entities( ring_mesh.m_bulk_data, element_rank, elements );
  stk::mesh::Entity element = elements.front();

  // get one of the nodes related to this element
  ASSERT_TRUE(bulk.num_nodes(element) > 0);
  stk::mesh::Entity node = *bulk.begin_nodes(element);
  ASSERT_EQ( bulk.entity_rank(node), static_cast<unsigned>(stk::mesh::BaseEntityRank) );

  // make a modification to the node by changing its parts
  ring_mesh.m_bulk_data.modification_begin();
  stk::mesh::PartVector parts;
  parts.push_back( &special_part );
  bulk.change_entity_parts( node, parts );

  // check that the node AND it's element are marked as modified
  ASSERT_EQ ( bulk.state(node), stk::mesh::Modified );
  ASSERT_EQ ( bulk.state(element), stk::mesh::Modified );

  ASSERT_TRUE ( bulk.modification_end() );
}

TEST(UnitTestingOfBulkData, testChangeEntityOwnerFromSelfToSelf)
{
  // It should be legal to "change" entity ownership from yourself to yourself.
  //
  // 1---3---5
  // | 1 | 2 |
  // 2---4---6
  //
  // To test this, we use the mesh above, with elem 1 going on rank 0 and
  // elem 2 going on rank 1. Nodes 3,4 are shared. After the mesh is set up
  // we change the ownership of a few nodes to the same proc that already
  // owns them.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Bail if we only have one proc
  if (p_size == 1) {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  EntityVector nodes;
  const unsigned nodes_per_elem = 4, nodes_per_side = 2;

  if (p_rank < 2) {
    // We're just going to add everything to the universal part
    stk::mesh::PartVector node_parts, elem_parts;
    stk::mesh::Part &node_part = meta_data.get_topology_root_part(stk::topology::NODE);
    node_parts.push_back(&node_part);
    stk::mesh::Part &quad4_part = meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
    elem_parts.push_back(&quad4_part);

    // Create element
    const EntityRank elem_rank = stk::topology::ELEMENT_RANK;
    Entity elem = mesh.declare_entity(elem_rank,
                                        p_rank+1, //elem_id
                                        elem_parts);

    // Create nodes
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
      nodes.push_back(mesh.declare_entity(NODE_RANK, id, node_parts));
    }

    // Add relations to nodes
    unsigned rel_id = 0;
    for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
      mesh.declare_relation( elem, *itr, rel_id );
    }
  }

  mesh.modification_end();

  std::vector<EntityProc> change ;
  if (p_rank < 2) {
    // Change ownership of some nodes to the same proc that owns them

    // Add a non-shared node to change list
    if ( p_rank == 0 ) {
      EntityProc entry( nodes.front(), p_rank ) ;
      change.push_back( entry );
    }
    else {
      EntityProc entry( nodes.back(), p_rank ) ;
      change.push_back( entry );
    }

    // Add a shared node to change list
    Entity shared_node = nodes[p_rank == 0 ? nodes_per_side : 0];
    EntityId expected_id = 3;
    Part& shared_part = meta_data.globally_shared_part();
    ASSERT_TRUE( has_superset(mesh.bucket(shared_node), shared_part) );
    ASSERT_EQ(mesh.identifier(shared_node), expected_id);
    if (mesh.parallel_owner_rank(shared_node) == p_rank) {
      EntityProc entry( shared_node, p_rank );
      change.push_back( entry );
    }
  }

  mesh.change_entity_owner(change);
}

void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem, int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move)
{
    const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
    for(unsigned i=0; i<bulk.num_nodes(elem); ++i) {
        if (bulk.parallel_owner_rank(nodes[i]) == bulk.parallel_rank()) {
            entities_to_move.push_back(stk::mesh::EntityProc(nodes[i], dest_proc));
        }
    }
}

TEST(UnitTestingOfBulkData, change_entity_owner_8quads_4procs)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if (numProcs != 4) {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    setup8Quad4ProcMesh2D(bulk);
    // setup now, now change
    //  move 6 from p1 to p3, move 3 from p2 to p0
    //     p0   p1   p2   p3               p0   /p       /p   p3
    //  11---12---13---14---15          11---12------13-----14---15
    //   | 5  | 6  | 7  | 8  |    ->     | 5  | 6/3  | 7/2  | 8  |
    //   6----7----8----9---10           6----7------8------9---10
    //   | 1  | 2  | 3  | 4  |           | 1  | 2/1  | 3/0  | 4  |
    //   1----2----3----4----5           1----2------3------4----5
    // also moves any owned nodes of elem3 or elem6 along with the element
    std::vector<stk::mesh::EntityProc> entities_to_move;
    if (bulk.parallel_rank() == 1) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 6);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(bulk, elem, dest_proc, entities_to_move);
    }
    if (bulk.parallel_rank() == 2) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 3);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        add_nodes_to_move(bulk, elem, dest_proc, entities_to_move);
    }
    bulk.change_entity_owner(entities_to_move);

    std::vector<unsigned> counts(meta.entity_rank_count());
    stk::mesh::Selector owned_or_shared = meta.locally_owned_part() | meta.globally_shared_part();
    stk::mesh::count_entities(owned_or_shared, bulk, counts);

    unsigned numNodes = counts[stk::topology::NODE_RANK];
    unsigned numElems = counts[stk::topology::ELEM_RANK];

    unsigned expectedNumNodes = 10;
    unsigned expectedNumElems = 3;
    if(bulk.parallel_rank() == 1 || bulk.parallel_rank() == 2) {
        expectedNumElems = 1;
        expectedNumNodes = 4;
    }
    EXPECT_EQ(expectedNumNodes, numNodes);
    EXPECT_EQ(expectedNumElems, numElems);

    //center center point
    stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, 8);
    EXPECT_TRUE(bulk.is_valid(node8));
    int expectedOwner = 3;
    EXPECT_EQ(expectedOwner, bulk.parallel_owner_rank(node8));

    stk::mesh::PairIterEntityComm sharedNodes = bulk.entity_comm_map_shared(bulk.entity_key(node8));
    unsigned expectNumShared = 3;
    EXPECT_EQ(expectNumShared, sharedNodes.size());

    //top center point
    stk::mesh::Entity node13 = bulk.get_entity(stk::topology::NODE_RANK, 13);
    if (bulk.parallel_rank() == 3 || bulk.parallel_rank() == 2) {
        EXPECT_TRUE(bulk.in_shared(bulk.entity_key(node13)));
    } else {
        EXPECT_FALSE( bulk.in_shared(bulk.entity_key(node13)) );
    }
    expectedOwner = 3;
    EXPECT_EQ(expectedOwner, bulk.parallel_owner_rank(node13));

    //bottom center point
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    if (bulk.parallel_rank() == 0 || bulk.parallel_rank() == 1) {
        EXPECT_TRUE(bulk.in_shared(bulk.entity_key(node3)));
    } else {
        EXPECT_FALSE( bulk.in_shared(bulk.entity_key(node3)) );
    }
    expectedOwner = 1;
    EXPECT_EQ(expectedOwner, bulk.parallel_owner_rank(node3));
}

TEST(UnitTestingOfBulkData, change_entity_owner_internal_get_current_sharing_of_owned_nodes)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if (numProcs != 4) {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    setup8Quad4ProcMesh2D(bulk);
    // setup now, now change
    //  move 6 from p1 to p3, move 3 from p2 to p0
    //     p0   p1   p2   p3               p0    /p     p2   p3
    //  11---12---13---14---15          11---12------13---14---15
    //   | 5  | 6  | 7  | 8  |    ->     | 5  | 6/3  | 7  | 8  |
    //   6----7----8----9---10           6----7------8----9---10
    //   | 1  | 2  | 3  | 4  |           | 1  | 2/1  | 3  | 4  |
    //   1----2----3----4----5           1----2------3----4----5
    // also moves any owned nodes of elem3 or elem6 along with the element
    if (bulk.parallel_rank() == 1) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 6);
        int dest_proc = 3;
        std::map<Entity, std::set<int> > owned_node_sharing_map;
        stk::mesh::internal_get_current_sharing_of_owned_nodes(bulk,stk::mesh::EntityProc(elem, dest_proc),
                                                    owned_node_sharing_map);
        std::map<Entity, std::set<int> > expected_map;
        expected_map[bulk.get_entity(stk::topology::NODE_RANK, 8)].insert(2);
        expected_map[bulk.get_entity(stk::topology::NODE_RANK, 13)].insert(2);
//        typedef std::map<Entity, std::set<int> > map_of_sets;
//        for (map_of_sets::iterator i=owned_node_sharing_map.begin(); i != owned_node_sharing_map.end(); ++i)
//        {
//            std::cout << "node id=" << bulk.identifier(i->first) << "\n";
//            std::set<int> & this_set = i->second;
//            for (std::set<int>::iterator j=this_set.begin() ; j != this_set.end(); ++j)
//            {
//                std::cout << "sharing proc=" << *j << "\n";
//            }
//        }
        bool result = expected_map == owned_node_sharing_map;
        EXPECT_TRUE(result);

    }
}

TEST(UnitTestingOfBulkData, testChangeEntityOwnerOfShared)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //  proc 0 owns element 1
  //  ...
  //  proc n owns element n+1
  //
  //   1---3---5---7---9---11--13
  //   | 1 | 2 | 3 | 4 | 5 | 6 | ...
  //   2---4---6---8---10--12--14
  //     this edge moves --^
  //     element 5 moves to proc 0.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  Part& edge_part = meta_data.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();
  const EntityRank node_rank = stk::topology::NODE_RANK;
  const EntityRank edge_rank = stk::topology::EDGE_RANK;
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;

  // Bail if we have fewer than 3 procs
  if (p_size < 3) {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  const unsigned nodes_per_elem = 4, nodes_per_side = 2;
  EntityKey elem_key_chg_own(elem_rank, p_size - 1 /*id*/);
  EntityKey edge_key_chg_own(edge_rank, 1 /*id*/);
  EntityKey node_A_key_chg_own(node_rank, p_size*2-1 /*id*/);
  EntityKey node_B_key_chg_own(node_rank, p_size*2 /*id*/);
  EntityKey node_C_key(node_rank, p_size*2-1-2 /*id*/);
  EntityKey node_D_key(node_rank, p_size*2-2 /*id*/);

  // Create element
  Entity elem = mesh.declare_entity(elem_rank,
                                      p_rank+1, //elem_id
                                      elem_part);

  // If it is 2nd to last element, it is the one changing
  if (p_rank == p_size - 2) {
    EXPECT_TRUE(elem_key_chg_own == mesh.entity_key(elem));
  }

  // Create nodes
  EntityVector nodes;
  const unsigned starting_node_id = p_rank * nodes_per_side + 1;
  for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, id, node_part));
  }
  // proc 0:  1,2,3,4
  // proc 1:  3,4,5,6
  // proc 2:  5,6,7,8

  // Add element relations to nodes
  unsigned rel_id = 0;
  for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
    mesh.declare_relation( elem, *itr, rel_id );
  }

  // Create edge on last two procs

  if (p_rank >= p_size - 2) {
    Entity edge = mesh.declare_entity(edge_rank,
                                       1, // id
                                       edge_part);
    EXPECT_TRUE(mesh.entity_key(edge) == edge_key_chg_own);

    // Add element relation to edge
    mesh.declare_relation( elem, edge, 1 /*rel-id*/);

    // Add edge relations to nodes
    unsigned start_idx = p_rank == p_size - 1 ? 0 : nodes_per_side;
    unsigned end_idx = start_idx + nodes_per_side;
    rel_id = 0;
    for (unsigned idx = start_idx ;
         start_idx < end_idx;
         ++start_idx, ++rel_id) {
      mesh.declare_relation( edge, nodes[idx], rel_id );
    }
  }
  if (p_rank == 0) {
      mesh.add_node_sharing(nodes[2],1);
      mesh.add_node_sharing(nodes[3],1);
  }
  else if (p_rank == p_size-1) {
      mesh.add_node_sharing(nodes[0],p_rank-1);
      mesh.add_node_sharing(nodes[1],p_rank-1);
  }
  else {
      mesh.add_node_sharing(nodes[0],p_rank-1);
      mesh.add_node_sharing(nodes[1],p_rank-1);
      mesh.add_node_sharing(nodes[2],p_rank+1);
      mesh.add_node_sharing(nodes[3],p_rank+1);
  }

  mesh.modification_end();

  // Changing elem and edge should be ghosted or unknown on proc 0
  if (p_rank == 0) {
    // Get the two entities
    Entity changing_elem = mesh.get_entity(elem_key_chg_own);
    Entity changing_edge = mesh.get_entity(edge_key_chg_own);
    Entity changing_node_A = mesh.get_entity(node_A_key_chg_own);
    Entity changing_node_B = mesh.get_entity(node_B_key_chg_own);

    if (p_size == 3) {
      // Should be ghosted
      ASSERT_TRUE(mesh.is_valid(changing_elem));
      ASSERT_TRUE(mesh.is_valid(changing_edge));
      ASSERT_TRUE(mesh.is_valid(changing_node_A));
      ASSERT_TRUE(mesh.is_valid(changing_node_B));

      // Verify that the entities are ghosted
      Part& owned = meta_data.locally_owned_part();
      Part& shared = meta_data.globally_shared_part();
      EXPECT_TRUE(!(mesh.bucket(changing_elem).member(owned) ||
                       mesh.bucket(changing_elem).member(shared)));
      EXPECT_TRUE(!(mesh.bucket(changing_edge).member(owned) ||
                       mesh.bucket(changing_edge).member(shared)));
      EXPECT_TRUE(!(mesh.bucket(changing_node_A).member(owned) ||
                       mesh.bucket(changing_node_A).member(shared)));
      EXPECT_TRUE(!(mesh.bucket(changing_node_B).member(owned) ||
                       mesh.bucket(changing_node_B).member(shared)));

    }
    else {
      // Should be invalid
      EXPECT_TRUE(!mesh.is_valid(changing_elem));
      EXPECT_TRUE(!mesh.is_valid(changing_edge));
      EXPECT_TRUE(!mesh.is_valid(changing_node_A));
      EXPECT_TRUE(!mesh.is_valid(changing_node_B));
    }
  }


  std::vector<EntityProc> change ;
  if (p_rank == (p_size-2) ) {
    // Change ownership of changing elem and all entities in it's closure that
    // we own to proc 0.

    Entity changing_elem = mesh.get_entity(elem_key_chg_own);
    ASSERT_TRUE( mesh.is_valid(changing_elem) );
    EntityProc eproc(changing_elem, 0 /*new owner*/);
    change.push_back(eproc);

    const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
    {
      stk::mesh::Entity const *to_i = mesh.begin(changing_elem, irank);
      stk::mesh::Entity const *to_e = mesh.end(changing_elem, irank);
      for ( ; to_i != to_e; ++to_i)
      {
        if (mesh.parallel_owner_rank(*to_i) == p_rank)
        {
          EntityProc eproc(*to_i, 0 /*new owner*/);
          change.push_back(eproc);
        }
      }
    }
  }

  // What state are we in now?
  if (p_rank == p_size-1) {
      EXPECT_TRUE( mesh.in_shared(node_A_key_chg_own, p_rank-1) );
      EXPECT_TRUE( mesh.in_shared(node_B_key_chg_own, p_rank-1) );
      EXPECT_TRUE( mesh.in_shared(edge_key_chg_own, p_rank-1) );
  }
  else if (p_rank == p_size-2 ) {
      EXPECT_TRUE( mesh.in_shared(node_A_key_chg_own, p_rank+1) );
      EXPECT_TRUE( mesh.in_shared(node_B_key_chg_own, p_rank+1) );
      EXPECT_TRUE( mesh.in_shared(edge_key_chg_own, p_rank+1) );
      EXPECT_TRUE( mesh.in_shared(node_C_key, p_rank-1) );
      EXPECT_TRUE( mesh.in_shared(node_D_key, p_rank-1) );
  }
  else if ( (p_rank == p_size-3) && (p_size>3) ) {
      EXPECT_TRUE( mesh.in_shared(node_C_key, p_rank+1) );
      EXPECT_TRUE( mesh.in_shared(node_D_key, p_rank+1) );
  }
  else if (p_rank == 0) {
      if (p_size>3) {
          EXPECT_FALSE( mesh.in_shared(node_C_key, p_rank-3) );
          EXPECT_FALSE( mesh.in_shared(node_D_key, p_rank-3) );
      }
      EXPECT_FALSE( mesh.in_shared(node_A_key_chg_own, p_size-1) );
      EXPECT_FALSE( mesh.in_shared(node_B_key_chg_own, p_size-1) );
      EXPECT_FALSE( mesh.in_shared(edge_key_chg_own, p_size-1) );
  }

  mesh.change_entity_owner(change);

  Entity changing_elem = mesh.get_entity(elem_key_chg_own);
  Entity changing_edge = mesh.get_entity(edge_key_chg_own);
  Entity changing_node_A = mesh.get_entity(node_A_key_chg_own);
  Entity changing_node_B = mesh.get_entity(node_B_key_chg_own);

  // Changing elem and edge should now be owned by proc 0
  if (p_rank == 0) {
    // Get the two ghosted entities, check that they were found
    ASSERT_TRUE(mesh.is_valid(changing_elem));
    ASSERT_TRUE(mesh.is_valid(changing_edge));
    ASSERT_TRUE(mesh.is_valid(changing_node_A));
    ASSERT_TRUE(mesh.is_valid(changing_node_B));

    // Verify that the entities are ghosted
    Part& owned = meta_data.locally_owned_part();
    EXPECT_TRUE( mesh.bucket(changing_elem).member(owned) );
    EXPECT_TRUE( mesh.bucket(changing_edge).member(owned) );
    EXPECT_TRUE( mesh.bucket(changing_node_A).member(owned) );
    EXPECT_TRUE( mesh.bucket(changing_node_B).member(owned) );

    EXPECT_TRUE( mesh.in_shared( edge_key_chg_own, p_size-1 ) );
    EXPECT_TRUE( mesh.in_shared( node_A_key_chg_own, p_size-1 ) );
    EXPECT_TRUE( mesh.in_shared( node_B_key_chg_own, p_size-1 ) );
  }
  if ((p_rank == p_size-3) && (p_size > 3)) {
      EXPECT_TRUE( mesh.in_shared( node_C_key, 0 ) );
      EXPECT_TRUE( mesh.in_shared( node_D_key, 0 ) );
  }
  if (p_rank == p_size-2) {
      EXPECT_TRUE(!mesh.is_valid(changing_elem));
      EXPECT_TRUE(!mesh.is_valid(changing_edge));
      EXPECT_TRUE(!mesh.is_valid(changing_node_A));
      EXPECT_TRUE(!mesh.is_valid(changing_node_B));
  }
  if (p_rank == p_size-1) {
      EXPECT_TRUE( mesh.in_shared( edge_key_chg_own, 0 ) );
      EXPECT_TRUE( mesh.in_shared( node_A_key_chg_own, 0 ) );
      EXPECT_TRUE( mesh.in_shared( node_B_key_chg_own, 0 ) );
  }
}

TEST(UnitTestingOfBulkData, testFamilyTreeGhosting)
{
  // A family tree is a higher-rank entity (rank = element_rank() + 1) that
  // has down-relations to elements used, for example, to hold parent/child
  // relations in an adapted mesh.
  //
  // 1---3---5---7
  // | 1 | 2 | 3 | ...
  // 2---4---6---8
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc.
  // After the mesh is set up we add rank-3 (family tree) entities and have them point down to
  // just the single rank-2 elements.  Then we check that they are properly
  // ghosted after modification_end.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;

  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  entity_rank_names.push_back("FAMILY_TREE");

  MetaData meta_data(spatial_dim, entity_rank_names);
  const unsigned nodes_per_elem = 4, nodes_per_side = 2;
  Part &elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4 );
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  Part& owned  = meta_data.locally_owned_part();
  Part& shared = meta_data.globally_shared_part();

  //
  // Begin modification cycle so we can create the entities and relations
  //

  mesh.modification_begin();

  EntityVector nodes;
  const EntityRank family_tree_rank = static_cast<EntityRank>(stk::topology::ELEMENT_RANK + 1);
  const EntityId my_family_tree_id = p_rank+1;

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&elem_part);

  // Create element
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  Entity elem = mesh.declare_entity(elem_rank,
                                    p_rank+1, //elem_id
                                    elem_parts);

  // Create nodes
  const unsigned starting_node_id = p_rank * nodes_per_side + 1;
  for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
    nodes.push_back(mesh.declare_entity(NODE_RANK, id, empty_parts));
  }
  if (p_rank > 0) {
    mesh.add_node_sharing(nodes[0], p_rank - 1);
    mesh.add_node_sharing(nodes[1], p_rank - 1);
  }
  if (p_rank < (p_size - 1)) {
    mesh.add_node_sharing(nodes[2], p_rank + 1);
    mesh.add_node_sharing(nodes[3], p_rank + 1);
  }

  // Add relations to nodes
  unsigned rel_id = 0;
  for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
    mesh.declare_relation( elem, *itr, rel_id );
  }

  // Create family tree
  Entity family_tree = mesh.declare_entity(family_tree_rank,
                                           my_family_tree_id,
                                           empty_parts);
  // Add relation to element
  unsigned downward_ordinal = 0; // we only have 1 down relation, it has ordinal 0
  mesh.declare_relation( family_tree, elem, downward_ordinal);

  mesh.modification_end();

  //
  // Test correctness of ghosting: Check that adjacent family-trees are ghosted on this proc
  //

  // Compute and store ids of adjacent family-trees
  std::vector<EntityId> family_tree_ghost_ids;
  if (p_rank > 0) {
    family_tree_ghost_ids.push_back(my_family_tree_id - 1);
  }
  if (p_rank < p_size - 1) {
    family_tree_ghost_ids.push_back(my_family_tree_id + 1);
  }

  // Check that my_family_tree exists and I own it
  Entity my_family_tree = mesh.get_entity(family_tree_rank, my_family_tree_id);
  ASSERT_TRUE(mesh.is_valid(my_family_tree));
  ASSERT_TRUE( (p_rank) == mesh.parallel_owner_rank(my_family_tree));

  // Check that adjacent family-trees exist and are ghosted
  for (std::vector<EntityId>::const_iterator
       itr = family_tree_ghost_ids.begin(); itr != family_tree_ghost_ids.end(); ++itr) {
    int expected_ghosted_family_tree_id = *itr;

    Entity expected_ghosted_family_tree = mesh.get_entity(family_tree_rank, expected_ghosted_family_tree_id);
    ASSERT_TRUE(mesh.is_valid(expected_ghosted_family_tree));
    ASSERT_TRUE(expected_ghosted_family_tree_id - 1 == mesh.parallel_owner_rank(expected_ghosted_family_tree));

    stk::mesh::Bucket& bucket = mesh.bucket(expected_ghosted_family_tree);
    ASSERT_TRUE(!bucket.member(owned) && !bucket.member(shared));
  }
}

TEST(UnitTestingOfBulkData, testChangeEntityPartsOfShared)
{
  //
  // This unit-test is designed to test what happens when a shared entity
  // is moved on one processor during the same modification cycle in which
  // it was declared.
  //
  // 1---3---5
  // | 1 | 2 |
  // 2---4---6
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. Node 3 is the node we'll be testing.
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  const EntityRank node_rank = stk::topology::NODE_RANK;
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;

  stk::mesh::Part& extra_node_part = meta_data.declare_part("extra_node_part", node_rank);
  meta_data.commit();

  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Bail unless in parallel
  if (p_size == 1) {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  if (p_rank < 2) {
    mesh.modification_begin();

    const unsigned nodes_per_elem = 4, nodes_per_side = 2;
    EntityKey node_key_to_move(node_rank, 3 /*id*/);

    // We're just going to add everything to the universal part
    stk::mesh::PartVector empty_parts, elem_parts;
    stk::mesh::Part &quad4_part = meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
    elem_parts.push_back(&quad4_part);

    // Create element
    Entity elem = mesh.declare_entity(elem_rank,
                                        p_rank+1, //elem_id
                                        elem_parts);

    // Create nodes
    EntityVector nodes;
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
      nodes.push_back(mesh.declare_entity(NODE_RANK, id, empty_parts));
    }

    // Add relations to nodes
    unsigned rel_id = 0;
    for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
      mesh.declare_relation( elem, *itr, rel_id );
    }

    // On the processor that does *not* end up as the owner of the node, change its parts
    Entity changing_node = mesh.get_entity(node_key_to_move);
    if (p_rank == 1) {
      PartVector add_parts(1, &extra_node_part);
      mesh.change_entity_parts(changing_node, add_parts);
    }

    mesh.modification_end();

    // Expect that this is a shared node
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(changing_node)).empty());

    // Expect that part change had no impact since it was on the proc that did not end
    // up as the owner
    EXPECT_FALSE(mesh.bucket(changing_node).member(extra_node_part));

    mesh.modification_begin();

    // On the processor that owns the node, change its parts
    if (p_rank == 0) {
      PartVector add_parts(1, &extra_node_part);
      mesh.change_entity_parts(changing_node, add_parts);
    }

    mesh.modification_end();

    // Expect that the part change *did* have an impact
    EXPECT_TRUE(mesh.bucket(changing_node).member(extra_node_part));
  }
  else {
    // On extra procs, do bare minimum
    mesh.modification_begin();
    mesh.modification_end();
    mesh.modification_begin();
    mesh.modification_end();
  }
}

TEST(UnitTestingOfBulkData, testParallelSideCreation)
{
  //
  // This unit-test is designed to test what happens when a shared sides are created on
  // both processors that share the side with different global ids.  Then synchonization
  // is done as a second step.
  //
  // 1---3---5
  // | 1 | 2 |
  // 2---4---6
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. Node 3 is the node we'll be testing.
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  const EntityRank node_rank = stk::topology::NODE_RANK;
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const EntityRank side_rank = stk::topology::EDGE_RANK;

  stk::mesh::Part& side_part = meta_data.declare_part_with_topology("side_part", stk::topology::LINE_2);
  stk::mesh::Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4);

  meta_data.commit();

  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Bail unless in parallel
  if (p_size == 1) {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  if (p_rank < 2) {
    mesh.modification_begin();

    const unsigned nodes_per_elem = 4, nodes_per_side = 2;
    EntityKey node_key_to_move(node_rank, 3 /*id*/);

    // We're just going to add everything to the universal part
    stk::mesh::PartVector empty_parts;

    // Create nodes
    EntityVector nodes;
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
      nodes.push_back(mesh.declare_entity(NODE_RANK, id, empty_parts));
    }

    // Create element
    const EntityId elem_id = p_rank+1;
    //Entity elem = mesh.declare_entity(elem_rank, elem_id, empty_parts);
    Entity elem = mesh.declare_entity(elem_rank, elem_id, elem_part);

    // Create local version of side (with different, artificial id on each processor)
    const EntityId tmp_side_id = p_rank+51;
    //Entity side = mesh.declare_entity(side_rank, tmp_side_id, empty_parts);
    Entity side = mesh.declare_entity(side_rank, tmp_side_id, side_part);

    // Add element relations to nodes
    unsigned elem_rel_id = 0;
    for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++elem_rel_id) {
      mesh.declare_relation( elem, *itr, elem_rel_id );
    }

    // Add side relations to nodes and element
    EntityVector side_nodes;
    side_nodes.push_back( mesh.get_entity( NODE_RANK, 3 ) );
    side_nodes.push_back( mesh.get_entity( NODE_RANK, 4 ) );
    unsigned side_rel_id = 0;
    for (EntityVector::iterator itr = side_nodes.begin(); itr != side_nodes.end(); ++itr, ++side_rel_id) {
      mesh.declare_relation( side, *itr, side_rel_id );
    }
    mesh.declare_relation( elem, side, 0 );

    mesh.modification_end();

    // Expect that the side is not shared, but the nodes of side are shared
    EXPECT_TRUE(mesh.entity_comm_map_shared(mesh.entity_key(side)).empty());
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(side_nodes[0])).empty());
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(side_nodes[1])).empty());

    // Now "detect" that there is a duplicate aura side using the side nodes
    EntityVector sides;
    get_entities_through_relations(mesh, side_nodes, side_rank, sides);
    EXPECT_EQ(2u, sides.size());

    mesh.modification_begin();

    // Delete the local side and create new, shared side
    mesh.destroy_relation(elem, side, 0);
    bool successfully_destroyed = mesh.destroy_entity(side);
    EXPECT_TRUE(successfully_destroyed);

    const EntityId side_id = 1;
    side = mesh.declare_entity(side_rank, side_id, side_part);

    // Add side relations to nodes and element
    side_rel_id = 0;
    for (EntityVector::iterator itr = side_nodes.begin(); itr != side_nodes.end(); ++itr, ++side_rel_id) {
      mesh.declare_relation( side, *itr, side_rel_id );
    }
    mesh.declare_relation( elem, side, 0 );

    mesh.modification_end();

    // Expect that the side is shared, and nodes of side are shared
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(side)).empty());
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(side_nodes[0])).empty());
    EXPECT_FALSE(mesh.entity_comm_map_shared(mesh.entity_key(side_nodes[1])).empty());

    // Check that there is only a single side using the side nodes
    get_entities_through_relations(mesh, side_nodes, side_rank, sides);
    EXPECT_EQ(1u, sides.size());
  }
  else {
    // On extra procs, do bare minimum
    mesh.modification_begin();
    mesh.modification_end();
    mesh.modification_begin();
    mesh.modification_end();
  }
}

TEST(UnitTestingOfBulkData, test_final_modification_end)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();

  BulkData mesh(meta_data, pm);

  mesh.modification_begin();
  mesh.final_modification_end();

  ASSERT_THROW(mesh.modification_begin(), std::logic_error );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of field_data_footprint(.)
TEST( UnitTestingOfBulkData, test_total_field_data_footprint )
{
  // Test 3x1x1 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  const stk::mesh::BulkData &mesh = hf.m_bulk_data;

  // Call function we're testing
  size_t field_data_footprint = mesh.total_field_data_footprint(stk::topology::NODE_RANK);

  // Alternative computation explicitly gathers buckets.
  size_t node_fields_footprint = 0;
  const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
  for (size_t i = 0; i < node_buckets.size(); ++i)
  {
    node_fields_footprint +=
        node_buckets[i]->capacity() * field_bytes_per_entity(hf.m_coord_field, *node_buckets[i]);
  }

  EXPECT_EQ(node_fields_footprint, field_data_footprint);
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of get_buckets and get_entities functions
TEST( UnitTestingOfBulkData, test_get_entities )
{
  // Test 3x4x4 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 4;
  const unsigned NZ = 40;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

  hf.m_meta.commit();
  hf.generate_mesh();
  const stk::mesh::BulkData &mesh = hf.m_bulk_data;

  Selector select_owned( MetaData::get(mesh).locally_owned_part() );
  const stk::mesh::BucketVector &bucket_ptrs = mesh.get_buckets(stk::topology::NODE_RANK, select_owned);
  stk::mesh::EntityVector entities;
  mesh.get_entities(stk::topology::NODE_RANK, select_owned, entities);
  //
  //  Confirm that the number of entities exracted by either bucket or entity access is identical
  //
  int numBucketEntities=0;
  std::map<stk::mesh::Entity, int> entityMap;
  for(unsigned int ibucket=0; ibucket<bucket_ptrs.size(); ++ibucket) {
    numBucketEntities += bucket_ptrs[ibucket]->size();
    for(unsigned iobj=0; iobj< bucket_ptrs[ibucket]->size(); ++iobj) {
      entityMap[(*bucket_ptrs[ibucket])[iobj]] = 1;
    }
  }
  int numEntities = entities.size();
  EXPECT_EQ(numBucketEntities, numEntities);
  //
  //  Confirm that the actual contents of the entity lists are identical
  //
  for(unsigned int iobj=0; iobj<entities.size(); ++iobj) {
    ASSERT_TRUE(entityMap.find(entities[iobj]) != entityMap.end());
  }  
  //
  //  confirm the total number of entities is the expected (41*5*4) = 820, the total number of unique nodes in the mesh
  //
  int globalNumEntities = numEntities;
  stk::all_reduce(MPI_COMM_WORLD, stk::ReduceSum<1>(&globalNumEntities));

  EXPECT_EQ(globalNumEntities, 820);
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of communicate_field_data

typedef stk::mesh::Field<int>  PressureFieldType ;

static void test_sync_1(stk::mesh::BulkData& eMesh, PressureFieldType& pressure_field, bool sync_shared, bool sync_aura)
{
  unsigned p_rank = eMesh.parallel_rank();
  unsigned p_size = eMesh.parallel_size();
  static_cast<void>(p_size);

  const stk::mesh::BucketVector & buckets = eMesh.buckets( stk::topology::NODE_RANK );

  enum Type{ Owned, Shared, Ghost };

  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
  {
    {
      stk::mesh::Bucket & bucket = **k ;

      const unsigned num_elements_in_bucket = bucket.size();

      for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
      {
        stk::mesh::Entity entity = bucket[iEntity];
        int * const p = stk::mesh::field_data( pressure_field , entity );
        stk::mesh::EntityId id=eMesh.identifier(entity);

        if (bucket.owned())
        {
          p[0] = (p_rank+1)*100+id;
        }
        else if (bucket.shared())
        {
          p[0] = -((eMesh.parallel_owner_rank(entity)+1)*100+id);
        }
        else
        {
          p[0] = ((p_rank+1)*1000 + id);
        }

      }
    }
  }

  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(&pressure_field);

    // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
    if (sync_aura) stk::mesh::communicate_field_data(eMesh.aura_ghosting(), fields);

    // the shared part (just the shared boundary)
    if (sync_shared) stk::mesh::copy_owned_to_shared(eMesh, fields);
  }

  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
  {
    {
      stk::mesh::Bucket & bucket = **k ;

      const unsigned num_elements_in_bucket = bucket.size();

      for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
      {
        stk::mesh::Entity entity = bucket[iEntity];
        stk::mesh::EntityId id = eMesh.identifier(entity);
        int * const p = stk::mesh::field_data( pressure_field , entity );
        double p_e = (p_rank+1)*100+id;
        if (bucket.owned())
        {
          ASSERT_EQ(p[0], p_e);
        }
        else if (bucket.shared())
        {
          p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
          if (sync_shared)
          {
            ASSERT_EQ(p[0], p_e);
          }
        }
        else
        {
          p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
          if (sync_aura)
            ASSERT_EQ(p[0], p_e);
        }
      }
    }
  }
}

TEST(UnitTestingOfBulkData, testFieldComm)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // run this with exercise_field_sync_bug = true, and 3 <= nprocs <= 4 to show the possible bug
  bool exercise_field_sync_bug = true;

  const unsigned p_size = stk::parallel_machine_size( pm );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );

  meta.commit();

  //------------------------------
  // 3d seems to be fine...
  if (p_size <= 4)
  {
    const int root_box[3][2] = { { 0 , 2 } , { 0 , 2 } , { 0 , 1 } };  // emulate 2d box

    BoxFixture fixture( pm, 100 );
    PressureFieldType& p_field = fixture.fem_meta().declare_field<PressureFieldType>(stk::topology::NODE_RANK, "p");
    stk::mesh::put_field( p_field , fixture.fem_meta().universal_part());
    fixture.fem_meta().commit();
    BulkData & bulk = fixture.bulk_data();
    int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    bulk.modification_begin();
    fixture.generate_boxes( root_box, local_box );
    bulk.modification_end();

    {
      bool shared_aura = false;
      bool shared = false;
      test_sync_1(bulk, p_field, shared, shared_aura);
      test_sync_1(bulk, p_field, false, true);
      if (exercise_field_sync_bug || p_size <=2)
        {
          test_sync_1(bulk, p_field, true, false);
          test_sync_1(bulk, p_field, true, true);
        }
    }
  }

  //------------------------------
  // 2d, not so much
  if (p_size <= 4)
  {
    stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
    PressureFieldType& p_field = fixture.m_meta.declare_field<PressureFieldType>(stk::topology::NODE_RANK, "p");
    stk::mesh::put_field( p_field , fixture.m_meta.universal_part());
    fixture.m_meta.commit();
    fixture.generate_mesh();
    stk::mesh::BulkData & bulk = fixture.m_bulk_data;

    {
      bool shared_aura = false;
      bool shared = false;
      test_sync_1(bulk, p_field, shared, shared_aura);
      test_sync_1(bulk, p_field, false, true);
      if (exercise_field_sync_bug || p_size <=2)
        {
          test_sync_1(bulk, p_field, true, false);
          test_sync_1(bulk, p_field, true, true);
        }
    }
  }
}

// testing comm lists and custom ghosting

TEST(UnitTestingOfBulkData, testCommList)
{
  /**
   * This is a boiled-down version of a stk_adapt situation that is failing
   *   where a custom ghosted node is later shared because it becomes part
   *   of an element (this occurs during hanging-node refinement).  It is
   *   a definite edge case, but exposed an assumption in the comm_mesh_verify_parallel_consistency()
   *   function, that a comm list can't have both a shared and custom ghosted node.
   * This test currently just verifies the issue exists.  When comm_mesh_verify_parallel_consistency()
   *   is fixed, this test should be modified to enforce no failure under debug mode.
   *
   *  Mesh
   *    7---8---9  P0 owns nodes 1,2,4,5; P, elem 1
   *    | 3 | 4 |  P1 : 3,6, elem 2
   *    4---5---6  P2 : 7,8, elem 3
   *    | 1 | 2 |  P3 : 9,   elem 4
   *    1---2---3
   *
   *  node 5 ghosted to proc 3, node 9 ghosted to proc 0
   * 
   *  Comm list looks like this after the ghosting operations (obtained from print_comm_list()):
   *    legend: (ghost_id, proc)
   *    P3: NODE[9] owner(3)     (1,0) (1,1) (1,2) (2,0)
   *    P0: NODE[9] owner(3)     (1,3) (2,3)
   *
   *  elem 1 modified to replace relation elem[1] -> node[1] with node[9] (previously ghosted)
   *
   *  This last step induces a comm list like this 
   *
   *    P3: NODE[9] owner(3) mod (0,0) (1,1) (1,2) (2,0)
   *    P0: NODE[9] owner(3) mod (0,3) (2,3)
   *
   *  Note the repetition of proc 3 in both the shared (ghost_id=0) and ghosted (id=2) list
   *  for P0 (and repetition of proc 0 in P3's list).  This causes the
   *  comm_mesh_verify_parallel_consistency() test to fail, although
   *  this should be a valid mesh, e.g., for periodic b.c.'s we might want this situation.
   *
   */

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const unsigned p_size = stk::parallel_machine_size( pm );
  const unsigned p_rank = stk::parallel_machine_rank( pm );

  //------------------------------
  if (p_size != 4)
    return;

  //------------------------------
  // test begin/end pair
  {
    stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
    fixture.m_meta.commit();
    fixture.generate_mesh();
    stk::mesh::BulkData & bulk = fixture.m_bulk_data;
    bulk.modification_begin();
    bulk.modification_end();
  }

  //------------------------------
  // test begin/end pair with mesh mods
  {
    stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
    fixture.m_meta.commit();
    fixture.generate_mesh();
    stk::mesh::BulkData & bulk = fixture.m_bulk_data;
    bulk.modification_begin();
    bulk.modification_end();

    bulk.modification_begin();

    // add some custom ghosting
    stk::mesh::Ghosting & ghosting = bulk.create_ghosting( std::string("new_nodes") );
    std::vector<stk::mesh::EntityKey> receive;
    ghosting.receive_list( receive );

    std::vector<stk::mesh::EntityProc> nodes_to_ghost;

    if (p_rank == 0)
      {
        stk::mesh::Entity node_5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
        EXPECT_TRUE(bulk.is_valid(node_5));
        nodes_to_ghost.push_back( stk::mesh::EntityProc(node_5, 3) );
      }

    if (p_rank == 3)
      {
        stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 9);
        EXPECT_TRUE(bulk.is_valid(node));
        nodes_to_ghost.push_back( stk::mesh::EntityProc(node, 0) );
      }

    bulk.change_ghosting( ghosting, nodes_to_ghost, receive);

    bulk.modification_end();

    bulk.modification_begin();

    if (p_rank == 0)
      {
        stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 9);
        EXPECT_TRUE(bulk.is_valid(node));
        stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
        EXPECT_TRUE(bulk.is_valid(node1));
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
        EXPECT_TRUE(bulk.is_valid(elem));
        bulk.destroy_relation( elem, node1, 0);
        bulk.declare_relation( elem , node , 0 );
      }

    bool failed=false;
    try {
      bulk.modification_end();
    }
    catch ( const std::exception & X ) {
      failed = true;
    }
    EXPECT_FALSE(failed);

  }

}

std::string printGhostData(stk::mesh::BulkData & bulkData, stk::mesh::Entity entity)
{
    std::ostringstream oss;
    std::vector<stk::mesh::EntityGhostData> egd;
    stk::mesh::get_ghost_data(bulkData, entity, egd);
    for (size_t z=0 ; z<egd.size() ; ++z) {
        oss << "P" << bulkData.parallel_rank() << ":  " << egd[z] << std::endl;
    }
    return oss.str();
}

std::string printGhostDataByRank(stk::mesh::BulkData & bulkData, stk::topology::rank_t rank)
{
    const stk::mesh::BucketVector & buckets = bulkData.buckets(rank);
    std::ostringstream oss;
    oss << "P" << bulkData.parallel_rank() << ":  rank=" << rank << std::endl;
    for (size_t k=0 ; k<buckets.size() ; ++k) {
        const stk::mesh::Bucket::iterator begin = buckets[k]->begin();
        const stk::mesh::Bucket::iterator end = buckets[k]->end();
        for (stk::mesh::Bucket::iterator it=begin ; it != end ; ++it) {
            oss << printGhostData(bulkData,*it);
        }
    }
    return oss.str();
}


TEST(UnitTestOfEntityGhostData, EntityGhostData)
{
    std::string gold_result = "(Entity_lid=0, direction=SEND, processor=128, ghosting level=LOCALLY_OWNED)";
    stk::mesh::EntityGhostData data;
    data.direction = stk::mesh::EntityGhostData::SEND;
    data.ghostingLevel = stk::mesh::EntityGhostData::LOCALLY_OWNED;
    data.processor = 128;
    std::ostringstream oss;
    oss << data;
    EXPECT_EQ( gold_result, oss.str() );
}


TEST(UnitTestOfEntityGhostData, get_ghost_data)
{
    using std::string;
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);

    if (psize == 3) { // Skip unless we're on 3 processors
        stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
        const string generatedMeshSpecification = "generated:1x1x3";
        stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

        if (stkMeshBulkData.parallel_rank() == 0) {
            std::ostringstream oss;
            for (stk::topology::rank_t rank=stk::topology::NODE_RANK ; rank <= stk::topology::ELEMENT_RANK ; ++rank) {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            string gold_result =
                    string("P0:  rank=NODE_RANK\n") +
                    string("P0:  (Entity_gid=1, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=1, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=2, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=2, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=3, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=3, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=4, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=4, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=9, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=10, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=11, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=12, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=5, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=5, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n") +
                    string("P0:  (Entity_gid=5, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=6, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=6, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n") +
                    string("P0:  (Entity_gid=6, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=7, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=7, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n") +
                    string("P0:  (Entity_gid=7, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=8, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=8, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n") +
                    string("P0:  (Entity_gid=8, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n") +
                    string("P0:  rank=EDGE_RANK\n") +
                    string("P0:  rank=FACE_RANK\n") +
                    string("P0:  rank=ELEMENT_RANK\n") +
                    string("P0:  (Entity_gid=1, rank=3, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n") +
                    string("P0:  (Entity_gid=1, rank=3, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P0:  (Entity_gid=2, rank=3, direction=RECEIVE, processor=1, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str() );
        }
        else if (stkMeshBulkData.parallel_rank() == 1) {
            std::ostringstream oss;
            for (stk::topology::rank_t rank=stk::topology::NODE_RANK ; rank <= stk::topology::ELEMENT_RANK ; ++rank) {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            std::string gold_result =
                    string("P1:  rank=NODE_RANK\n") +
                    string("P1:  (Entity_gid=5, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=6, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=7, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=8, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=1, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=2, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=3, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=4, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=13, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=14, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=15, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=16, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=9, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n") +
                    string("P1:  (Entity_gid=9, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=9, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=10, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n") +
                    string("P1:  (Entity_gid=10, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=10, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=11, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n") +
                    string("P1:  (Entity_gid=11, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=11, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=12, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n") +
                    string("P1:  (Entity_gid=12, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n") +
                    string("P1:  (Entity_gid=12, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n") +
                    string("P1:  rank=EDGE_RANK\n") +
                    string("P1:  rank=FACE_RANK\n") +
                    string("P1:  rank=ELEMENT_RANK\n") +
                    string("P1:  (Entity_gid=2, rank=3, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n") +
                    string("P1:  (Entity_gid=2, rank=3, direction=SEND, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=2, rank=3, direction=SEND, processor=2, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=1, rank=3, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P1:  (Entity_gid=3, rank=3, direction=RECEIVE, processor=2, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str() );
        }
        else { // if (stkMeshBulkData.parallel_rank() == 2)
            std::ostringstream oss;
            for (stk::topology::rank_t rank=stk::topology::NODE_RANK ; rank <= stk::topology::ELEMENT_RANK ; ++rank) {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            std::string gold_result =
                    string("P2:  rank=NODE_RANK\n") +
                    string("P2:  (Entity_gid=13, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n") +
                    string("P2:  (Entity_gid=13, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=14, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n") +
                    string("P2:  (Entity_gid=14, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=15, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n") +
                    string("P2:  (Entity_gid=15, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=16, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n") +
                    string("P2:  (Entity_gid=16, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=9, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n") +
                    string("P2:  (Entity_gid=10, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n") +
                    string("P2:  (Entity_gid=11, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n") +
                    string("P2:  (Entity_gid=12, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n") +
                    string("P2:  (Entity_gid=5, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=6, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=7, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=8, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n") +
                    string("P2:  rank=EDGE_RANK\n") +
                    string("P2:  rank=FACE_RANK\n") +
                    string("P2:  rank=ELEMENT_RANK\n") +
                    string("P2:  (Entity_gid=3, rank=3, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n") +
                    string("P2:  (Entity_gid=3, rank=3, direction=SEND, processor=1, ghosting level=AURA)\n") +
                    string("P2:  (Entity_gid=2, rank=3, direction=RECEIVE, processor=1, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str() );
        }
    }
}

TEST(DocTestBulkData, ChangeSharedOwner)
{
  // This test verifies that you can change the owner of shared nodes to
  // explicitly make the higher parallel-rank processor the owner.
  // This is also tested in UnitTestingOfBulkData.testChangeOwner_box, but it is a bit too complex.
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 3) {
    return;
  }

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x3";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  int target_proc = stkMeshBulkData.parallel_rank()==2 ? 2 : stkMeshBulkData.parallel_rank()+1;
  stk::mesh::EntityProcVec vec;
  stk::mesh::Selector locally_owned_and_globally_shared_selector = stkMeshMetaData.locally_owned_part() & stkMeshMetaData.globally_shared_part();
  {
      const stk::mesh::BucketVector& locally_owned_and_globally_shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, locally_owned_and_globally_shared_selector);
      int num_locally_owned_and_shared_nodes = 0;
      for (size_t k=0 ; k<locally_owned_and_globally_shared_buckets.size() ; ++k) {
          stk::mesh::Bucket & b = *locally_owned_and_globally_shared_buckets[k];
          for (size_t i=0 ; i<b.size() ; ++i) {
              stk::mesh::EntityProc tmp( b[i], target_proc );
              vec.push_back(tmp);
              ++num_locally_owned_and_shared_nodes;
          }
      }
      if (stkMeshBulkData.parallel_rank() == 0) {
          EXPECT_EQ( 4, num_locally_owned_and_shared_nodes );
      }
      else if (stkMeshBulkData.parallel_rank() == 1) {
          EXPECT_EQ( 4, num_locally_owned_and_shared_nodes );
      }
      else {
          EXPECT_EQ( 0, num_locally_owned_and_shared_nodes );
      }
  }
  stkMeshBulkData.change_entity_owner(vec);

  {
      const stk::mesh::BucketVector& locally_owned_and_globally_shared_buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, locally_owned_and_globally_shared_selector);
      int num_locally_owned_and_shared_nodes = 0;
      for (size_t k=0 ; k<locally_owned_and_globally_shared_buckets.size() ; ++k) {
          stk::mesh::Bucket & b = *locally_owned_and_globally_shared_buckets[k];
          for (size_t i=0 ; i<b.size() ; ++i) {
              stk::mesh::EntityProc tmp( b[i], target_proc );
              ++num_locally_owned_and_shared_nodes;
          }
      }
      if (stkMeshBulkData.parallel_rank() == 0) {
          EXPECT_EQ( 0, num_locally_owned_and_shared_nodes );
      }
      else if (stkMeshBulkData.parallel_rank() == 1) {
          EXPECT_EQ( 4, num_locally_owned_and_shared_nodes );
      }
      else {
          EXPECT_EQ( 4, num_locally_owned_and_shared_nodes );
      }
      // We also expect the previously locally owned & shared nodes to not be locally owned now.
      for (size_t i=0 ; i<vec.size() ; ++i) {
          EXPECT_TRUE( stkMeshBulkData.parallel_rank() != stkMeshBulkData.entity_comm_map_owner(stkMeshBulkData.entity_key(vec[i].first)) );
          EXPECT_EQ( vec[i].second, stkMeshBulkData.entity_comm_map_owner(stkMeshBulkData.entity_key(vec[i].first)) );
      }
  }

}

// Add to documentation tests for modification_end FAQ:
// 1.  How does change owner work?  Can I change the owner of a shared entity when I'm not the owner?
//     A:  Only the owner can give ownership to another processor
// 2.  Can I delete a shared and not-locally owned entity just on one processor?
// 3.  Can I delete a shared and locally owned entity just on one processor?
// 4.  Can I change parts on a shared entity differently on different sharing processors and will modification_end figure out the right parts?
//     A:  Only the owner can change parts on the entity.
// 5.  Can I add relations to shared entities differently on different sharing processors?
// 6.  Is there a difference in any of these between "creation" and "modification" cycles of the modification?
//     For example, can I create shared entities with different parts and it works?  vs.  Modifying shared entities to produce different parts?
TEST(DocTestBulkData, onlyTheOwnerCanChangeEntityParts)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 2) {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Part & myPart = stkMeshMetaData.declare_part("new_part");
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    stk::mesh::Selector allSharedSelector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> allSharedNodes;
    stk::mesh::get_selected_entities(allSharedSelector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), allSharedNodes);
    stkMeshBulkData.modification_begin();
    stk::mesh::PartVector addParts;
    addParts.push_back(&myPart);
    const int myRank = stk::parallel_machine_rank(communicator);
    for (size_t i=0 ; i<allSharedNodes.size() ; ++i) {
        if (stkMeshBulkData.parallel_owner_rank(allSharedNodes[i]) == myRank)
        {
            stkMeshBulkData.change_entity_parts(allSharedNodes[i],addParts);
        }
        else
        {
            EXPECT_THROW(stkMeshBulkData.change_entity_parts(allSharedNodes[i],addParts),std::logic_error);
        }
    }
    EXPECT_NO_THROW( stkMeshBulkData.modification_end() );

    // Verify parts are correct on all processors.
    stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> shared_nodes;
    stk::mesh::get_selected_entities(shared_selector,stkMeshBulkData.buckets(stk::topology::NODE_RANK),shared_nodes);
    for (size_t i=0 ; i<shared_nodes.size() ; ++i) {
        EXPECT_TRUE(stkMeshBulkData.bucket(shared_nodes[i]).member(myPart));
    }
}

TEST(DocTestBulkData, onlyKeepTheOwnersParts)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 2) {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA",stk::topology::NODE_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB",stk::topology::NODE_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::EntityKey node0Key;
    const stk::mesh::EntityId node_id = 2000;
    stkMeshBulkData.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;
    if (myRank == 0)
    {
       stk::mesh::Entity element0 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id_0, elem_parts);
       stk::mesh::Entity node0 = stkMeshBulkData.declare_entity(stk::topology::NODE_RANK, node_id, partA);
       node0Key = stkMeshBulkData.entity_key(node0);
       for (stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
       {
         stkMeshBulkData.declare_relation(element0,node0,node_rel_id);
       }
    }
    else  // myRank == 1
    {
       stk::mesh::Entity element1 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id_1, elem_parts);
       stk::mesh::Entity node0 = stkMeshBulkData.declare_entity(stk::topology::NODE_RANK, node_id, partB);
       node0Key = stkMeshBulkData.entity_key(node0);
       for (stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
       {
         stkMeshBulkData.declare_relation(element1,node0,node_rel_id);
       }
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Entity node0 = stkMeshBulkData.get_entity(node0Key);
    stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(node0);

    const int nodeOwner = 0;
    EXPECT_EQ( nodeOwner, stkMeshBulkData.parallel_owner_rank(node0) );
    EXPECT_TRUE( nodeBucket.member(stkMeshMetaData.globally_shared_part()) );

    EXPECT_TRUE( nodeBucket.member(partA) );
    EXPECT_FALSE( nodeBucket.member(partB) );
    EXPECT_TRUE( stk::mesh::has_superset(nodeBucket,partA) );
    EXPECT_FALSE( stk::mesh::has_superset(nodeBucket,partB) );

    stk::mesh::Entity element0 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_0);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_1);
    const int element0Owner = 0;
    const int element1Owner = 1;
    EXPECT_EQ( element0Owner, stkMeshBulkData.parallel_owner_rank(element0) );
    EXPECT_EQ( element1Owner, stkMeshBulkData.parallel_owner_rank(element1) );
}

TEST(DocTestBulkData, newSharedNodeGetMergedPartsFromElements)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 2) {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA",stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB",stk::topology::ELEMENT_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::EntityKey node0Key;
    const stk::mesh::EntityId node_id = 2000;
    stkMeshBulkData.modification_begin();
    if (myRank == 0)
    {
       const stk::mesh::EntityId element_id = 1000;
       stk::mesh::Entity element0 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id, partA);
       stk::mesh::Entity node0 = stkMeshBulkData.declare_entity(stk::topology::NODE_RANK, node_id);
       node0Key = stkMeshBulkData.entity_key(node0);
       stkMeshBulkData.change_entity_parts(element0, elem_parts, empty_parts);
       for (stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
       {
         stkMeshBulkData.declare_relation(element0,node0,node_rel_id);
       }
    }
    else  // myRank == 1
    {
       const stk::mesh::EntityId element_id = 1001;
       stk::mesh::Entity element1 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id, partB);
       stk::mesh::Entity node0 = stkMeshBulkData.declare_entity(stk::topology::NODE_RANK, node_id);
       node0Key = stkMeshBulkData.entity_key(node0);
       stkMeshBulkData.change_entity_parts(element1, elem_parts, empty_parts);
       for (stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
       {
         stkMeshBulkData.declare_relation(element1,node0,node_rel_id);
       }
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Entity node0 = stkMeshBulkData.get_entity(node0Key);
    stk::mesh::Bucket & node0Bucket = stkMeshBulkData.bucket(node0);

    const int node0Owner = 0;
    EXPECT_EQ( node0Owner, stkMeshBulkData.parallel_owner_rank(node0) );
    EXPECT_TRUE( node0Bucket.member(stkMeshMetaData.globally_shared_part()) );

    EXPECT_TRUE( node0Bucket.member(partA) );
    EXPECT_TRUE( node0Bucket.member(partB) );
    EXPECT_TRUE( stk::mesh::has_superset(node0Bucket,partA) );
    EXPECT_TRUE( stk::mesh::has_superset(node0Bucket,partB) );
}

TEST(DocTestBulkData, mayCreateRelationsToNodesDifferently)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 2) {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA",stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB",stk::topology::ELEMENT_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::Selector allSharedSelector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> allSharedNodes;
    stk::mesh::get_selected_entities(allSharedSelector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), allSharedNodes);

    std::sort(allSharedNodes.begin(),allSharedNodes.end(),stk::mesh::EntityLess(stkMeshBulkData));
    stk::mesh::Entity sharedNode0 = allSharedNodes[0];
    stk::mesh::Entity sharedNode1 = allSharedNodes[1];

    stkMeshBulkData.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;
    stk::mesh::Entity filler_node = stkMeshBulkData.declare_entity(stk::topology::NODE_RANK, 123456);
    if (myRank == 0)
    {
       stk::mesh::Entity element0 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id_0, partA);
       stkMeshBulkData.change_entity_parts(element0, elem_parts, empty_parts);
       stk::mesh::RelationIdentifier node_rel_id = 0;
       stkMeshBulkData.declare_relation(element0,sharedNode0,node_rel_id);
       for (stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
       {
         stkMeshBulkData.declare_relation(element0,filler_node,filler_rel_id);
       }
    }
    else  // myRank == 1
    {
       stk::mesh::Entity element1 = stkMeshBulkData.declare_entity(stk::topology::ELEMENT_RANK, element_id_1, partB);
       stkMeshBulkData.change_entity_parts(element1, elem_parts, empty_parts);
       const stk::mesh::RelationIdentifier node_rel_id = 0;
       stkMeshBulkData.declare_relation(element1,sharedNode0,node_rel_id);
       for (stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
       {
         stkMeshBulkData.declare_relation(element1,filler_node,filler_rel_id);
       }
    }
    EXPECT_NO_THROW( stkMeshBulkData.modification_end() );
    {
       stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(sharedNode0);
       EXPECT_TRUE( nodeBucket.member(partA) );
       EXPECT_TRUE( nodeBucket.member(partB) );
       EXPECT_EQ( 4u, stkMeshBulkData.num_elements(sharedNode0));
    }
    stk::mesh::Entity element0 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_0);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_1);
    stkMeshBulkData.modification_begin();
    if (myRank == 0)
    {
       stk::mesh::RelationIdentifier node_rel_id = 1;
       stkMeshBulkData.declare_relation(element0,sharedNode1,node_rel_id);
    }
    else // myRank == 1
    {
       stk::mesh::RelationIdentifier node_rel_id = 1;
       stkMeshBulkData.declare_relation(element1,sharedNode1,node_rel_id);

    }
    EXPECT_NO_THROW( stkMeshBulkData.modification_end() );

    {
       stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(sharedNode1);
       EXPECT_TRUE( nodeBucket.member(partA) );
       EXPECT_TRUE( nodeBucket.member(partB) );
       EXPECT_EQ( 4u, stkMeshBulkData.num_elements(sharedNode1));
    }
}


stk::mesh::PartVector setupFixture(stk::io::StkMeshIoBroker & io)
{
    const std::string generatedMeshSpecification = "generated:1x1x2";
    io.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    io.create_input_mesh();
    stk::mesh::Part & partA = io.meta_data().declare_part("PartA",stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = io.meta_data().declare_part("PartB",stk::topology::ELEMENT_RANK);
    io.populate_bulk_data();
    stk::mesh::PartVector pv;
    pv.push_back(&partA);
    pv.push_back(&partB);
    return pv;
}

stk::mesh::EntityVector getSortedNodes(stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), nodes);
    std::sort(nodes.begin(),nodes.end(),stk::mesh::EntityLess(bulk));
    return nodes;
}

bool no_aura_declare_relation( stk::mesh::BulkData &bulk,
                               stk::mesh::Entity e_from ,
                               stk::mesh::Entity e_to ,
                               const stk::mesh::RelationIdentifier local_id)
{
  stk::mesh::Bucket const *b_from       = bulk.bucket_ptr(e_from);
  stk::mesh::Bucket const *b_to         = bulk.bucket_ptr(e_to);
  const bool ownership_ok    = (b_from && b_to && !b_from->in_aura() && b_to->in_aura());

  if (ownership_ok)
  {
    bulk.declare_relation(e_from, e_to, local_id);
    return true;
  }
  else
  {
    return false;
  }
}


TEST(DocTestBulkData, inducedPartMembershipIgnoredForNonOwnedHigherRankedEntities)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 2) {
        return;
    }
    const int myRank = stk::parallel_machine_rank(communicator);

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    stk::mesh::PartVector pv = setupFixture(stkMeshIoBroker);
    stk::mesh::Part & partA = *pv[0];
    stk::mesh::Part & partB = *pv[1];

    stk::mesh::MetaData &meta = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &bulk = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (meta.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = meta.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::EntityVector ev = getSortedNodes(bulk, stk::mesh::Selector(meta.globally_shared_part()));
    ASSERT_TRUE( ev.size() >= 2 );
    stk::mesh::Entity sharedNodeA( ev[0] );
    stk::mesh::Entity sharedNodeB( ev[1] );

    stk::mesh::Bucket & nodeABucket = bulk.bucket(sharedNodeA);
    EXPECT_FALSE( nodeABucket.member(partA) );
    EXPECT_FALSE( nodeABucket.member(partB) );
    EXPECT_EQ( 2u, bulk.num_elements(sharedNodeA));
    EXPECT_EQ( 2u, bulk.num_elements(sharedNodeB));

    // First we create off processor elements
    bulk.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;

    const stk::mesh::EntityId filler_node_id = 123456;
    stk::mesh::Entity filler_node = bulk.declare_entity(stk::topology::NODE_RANK, filler_node_id);

    bulk.add_node_sharing(filler_node, (myRank == 1 ? 0 : 1));

    if (myRank == 0)
    {
       const stk::mesh::RelationIdentifier node_rel_id = 0;
       stk::mesh::Entity element = bulk.declare_entity(stk::topology::ELEMENT_RANK, element_id_0, partA);
       bulk.change_entity_parts(element, elem_parts, empty_parts);
       bulk.declare_relation(element,sharedNodeA,node_rel_id);
       for (stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
       {
         bulk.declare_relation(element,filler_node,filler_rel_id);
       }
    }
    else  // myRank == 1
    {
       const stk::mesh::RelationIdentifier node_rel_id = 0;
       stk::mesh::Entity element = bulk.declare_entity(stk::topology::ELEMENT_RANK, element_id_1, partB);
       bulk.change_entity_parts(element, elem_parts, empty_parts);
       bulk.declare_relation(element,sharedNodeA,node_rel_id);
       for (stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
       {
         bulk.declare_relation(element,filler_node,filler_rel_id);
       }
    }

    EXPECT_NO_THROW( bulk.modification_end() );
    {
       stk::mesh::Bucket & nodeABkt = bulk.bucket(sharedNodeA);
       EXPECT_TRUE( nodeABkt.member(partA) );
       EXPECT_TRUE( nodeABkt.member(partB) );
       EXPECT_EQ( 4u, bulk.num_elements(sharedNodeA));
    }

    stk::mesh::Entity element0 = bulk.get_entity(stk::topology::ELEMENT_RANK, element_id_0 );
    stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, element_id_1 );

    stk::mesh::Part &locally_owned_part = meta.locally_owned_part();
    stk::mesh::Bucket &element0_bucket = bulk.bucket(element0);
    stk::mesh::Bucket &element1_bucket = bulk.bucket(element1);

    if (myRank == 0)
    {
        EXPECT_TRUE(element0_bucket.member(locally_owned_part));
        EXPECT_FALSE(element1_bucket.member(locally_owned_part));
    }
    else // myRank == 1
    {
        EXPECT_FALSE(element0_bucket.member(locally_owned_part));
        EXPECT_TRUE(element1_bucket.member(locally_owned_part));
    }

    // Disallowing declaring relations to/from an entity in the aura can make it
    // easier to keep the mesh consistent.
    bulk.modification_begin();
    if (myRank == 0)
    {
      const stk::mesh::RelationIdentifier node_rel_id = 1;
      EXPECT_FALSE(no_aura_declare_relation(bulk, element1, sharedNodeB, node_rel_id) );
    }
    else // myRank == 1
    {
      const stk::mesh::RelationIdentifier node_rel_id = 1;
      EXPECT_FALSE(no_aura_declare_relation(bulk, element0, sharedNodeB, node_rel_id) );
    }
    EXPECT_NO_THROW( bulk.modification_end() );

    {
      stk::mesh::Bucket & nodeBBucket = bulk.bucket(sharedNodeB);
      EXPECT_FALSE( nodeBBucket.member(partA) );
      EXPECT_FALSE( nodeBBucket.member(partB) );
      EXPECT_EQ( 2u, bulk.num_elements(sharedNodeB));

      const int element0Owner = 0;
      EXPECT_EQ( element0Owner, bulk.parallel_owner_rank(element0) );
      const int element1Owner = 1;
      EXPECT_EQ( element1Owner, bulk.parallel_owner_rank(element1) );

      if (myRank == 0)
      {
        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);

        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
        EXPECT_TRUE( filler_node == bulk.begin_nodes(element1)[1]);
      }
      else // myRank == 1
      {
        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);
        EXPECT_TRUE( filler_node == bulk.begin_nodes(element0)[1]);

        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
      }
    }

    // Again we try to create relations from ghost element to shared node B.
    // BulkData::declare_relation(..) will do this, but the owner relations
    // will reset the ghost's relations in modification_end().
    //
    // Thus, more work needs to be done to enforce the relationship reciprocity
    // aspect of mesh consistency.
    bulk.modification_begin();
    if (myRank == 0)
    {
      const stk::mesh::RelationIdentifier node_rel_id = 1;
      bulk.declare_relation(element1, sharedNodeB, node_rel_id);
      EXPECT_TRUE( sharedNodeB == bulk.begin_nodes(element1)[1]);
    }
    else // myRank == 1
    {
      const stk::mesh::RelationIdentifier node_rel_id = 1;
      bulk.declare_relation(element0, sharedNodeB, node_rel_id);
      EXPECT_TRUE( sharedNodeB == bulk.begin_nodes(element0)[1]);
    }
    EXPECT_NO_THROW( bulk.modification_end() );

    {
      stk::mesh::Bucket & nodeBBucket = bulk.bucket(sharedNodeB);
      EXPECT_FALSE( nodeBBucket.member(partA) );
      EXPECT_FALSE( nodeBBucket.member(partB) );
      EXPECT_EQ( 3u, bulk.num_elements(sharedNodeB));

      if (myRank == 0)
      {
        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);

        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
        EXPECT_TRUE( filler_node == bulk.begin_nodes(element1)[1]);
      }
      else // myRank == 1
      {
        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);
        EXPECT_TRUE( filler_node == bulk.begin_nodes(element0)[1]);

        EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
      }
    }
}


} // empty namespace

//====================
extern int gl_argc;
extern char** gl_argv;

inline std::string getOption(const std::string& option, const std::string defaultString="no")
{
    std::string returnValue = defaultString;
    if ( gl_argv != 0 )
    {
        for (int i=0;i<gl_argc;i++)
        {
            std::string input_argv(gl_argv[i]);
            if ( option == input_argv )
            {
                if ( (i+1) < gl_argc )
                {
                    returnValue = std::string(gl_argv[i+1]);
                }
                break;
            }
        }
    }
    return returnValue;
}

namespace
{

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
        stk::mesh::BulkData(mesh_meta_data, comm){}
    virtual ~BulkDataTester() {}

    void my_internal_resolve_shared_modify_delete()
    {
        this->internal_resolve_shared_modify_delete();
    }

    void reset_closure_count(Entity entity)
    {
        m_closure_count[entity.local_offset()] = 0;
    }

    uint16_t closure_count(Entity entity)
    {
        return m_closure_count[entity.local_offset()];
    }
};

TEST(BulkData, ModificationEnd)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);

    if(numProcs == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester *stkMeshBulkData = new BulkDataTester(stkMeshMetaData, communicator);

        std::string exodusFileName = getOption("-i", "generated:1x1x4");

        // STK IO module will be described in separate chapter.
        // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
        // The order of the following lines in {} are important
        {
          stk::io::StkMeshIoBroker exodusFileReader(communicator);

          // Inform STK IO which STK Mesh objects to populate later
          exodusFileReader.set_bulk_data(*stkMeshBulkData);

          exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

          // Populate the MetaData which has the descriptions of the Parts and Fields.
          exodusFileReader.create_input_mesh();

          // Populate entities in STK Mesh from Exodus file
          exodusFileReader.populate_bulk_data();
        }

        int elementToMove = 3;
        int nodeToCheck = 9;

        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,nodeToCheck);
        stk::mesh::EntityKey entityToMoveKey(stk::topology::ELEMENT_RANK, elementToMove);

        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData->comm_list().begin(), stkMeshBulkData->comm_list().end(), nodeEntityKey);

        ASSERT_TRUE(iter != stkMeshBulkData->comm_list().end());
        EXPECT_EQ(nodeEntityKey,iter->key);
        EXPECT_TRUE(stkMeshBulkData->is_valid(iter->entity));

        stkMeshBulkData->modification_begin();

        ASSERT_TRUE ( stkMeshBulkData->is_valid(stkMeshBulkData->get_entity(entityToMoveKey)) );

        if ( stkMeshBulkData->parallel_rank() == 1 )
        {
            stkMeshBulkData->destroy_entity( stkMeshBulkData->get_entity(entityToMoveKey) );
        }

        // Really testing destroy_entity
        stk::mesh::delete_shared_entities_which_are_no_longer_in_owned_closure( *stkMeshBulkData );

        iter = std::lower_bound(stkMeshBulkData->comm_list().begin(), stkMeshBulkData->comm_list().end(), nodeEntityKey);

        ASSERT_TRUE(iter != stkMeshBulkData->comm_list().end());
        EXPECT_EQ(nodeEntityKey,iter->key);

        if ( stkMeshBulkData->parallel_rank() == 0 )
        {
            EXPECT_TRUE(stkMeshBulkData->is_valid(iter->entity));
        }
        else
        {
            EXPECT_FALSE(stkMeshBulkData->is_valid(iter->entity));
        }

    //    stkMeshBulkData->my_internal_resolve_shared_modify_delete();

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(*stkMeshBulkData, globalCounts);

        delete stkMeshBulkData;
    }
}

TEST(BulkData, verify_closure_count_is_correct)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    const int myRank   = stk::parallel_machine_rank(communicator);

    if(numProcs == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester *stkMeshBulkData = new BulkDataTester(stkMeshMetaData, communicator);

        std::string exodusFileName = getOption("-i", "generated:1x1x2");

        // STK IO module will be described in separate chapter.
        // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
        // The order of the following lines in {} are important
        {
          stk::io::StkMeshIoBroker exodusFileReader(communicator);

          // Inform STK IO which STK Mesh objects to populate later
          exodusFileReader.set_bulk_data(*stkMeshBulkData);

          exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

          // Populate the MetaData which has the descriptions of the Parts and Fields.
          exodusFileReader.create_input_mesh();

          // Populate entities in STK Mesh from Exodus file
          exodusFileReader.populate_bulk_data();
        }

        stk::mesh::EntityKey element_1_key(stk::topology::ELEMENT_RANK, 1);
        stk::mesh::EntityKey element_2_key(stk::topology::ELEMENT_RANK, 2);

        stk::mesh::EntityKey node_1_key(stk::topology::NODE_RANK, 1);
        stk::mesh::EntityKey node_2_key(stk::topology::NODE_RANK, 2);
        stk::mesh::EntityKey node_3_key(stk::topology::NODE_RANK, 3);
        stk::mesh::EntityKey node_4_key(stk::topology::NODE_RANK, 4);
        stk::mesh::EntityKey node_5_key(stk::topology::NODE_RANK, 5);
        stk::mesh::EntityKey node_6_key(stk::topology::NODE_RANK, 6);
        stk::mesh::EntityKey node_7_key(stk::topology::NODE_RANK, 7);
        stk::mesh::EntityKey node_8_key(stk::topology::NODE_RANK, 8);
        stk::mesh::EntityKey node_9_key(stk::topology::NODE_RANK, 9);
        stk::mesh::EntityKey node_10_key(stk::topology::NODE_RANK, 10);
        stk::mesh::EntityKey node_11_key(stk::topology::NODE_RANK, 11);
        stk::mesh::EntityKey node_12_key(stk::topology::NODE_RANK, 12);

        if (myRank == 0)
        {
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)) );

            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)) );
        }
        else // myRank == 1
        {
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)) );

            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)) );
        }

        std::vector<EntityProc> change_owner_vector;
        if (myRank == 1)
        {
            const int target_processor = 0;
            stk::mesh::Selector my_entities_selector = stkMeshBulkData->mesh_meta_data().locally_owned_part();
            std::vector<Entity> nodes;
            stk::mesh::get_selected_entities(my_entities_selector, stkMeshBulkData->buckets(stk::topology::NODE_RANK), nodes);
            for (size_t i=0 ; i<nodes.size() ; ++i) {
                change_owner_vector.push_back(EntityProc(nodes[i],target_processor));
            }
        }
        stkMeshBulkData->change_entity_owner(change_owner_vector);

        if (myRank == 0)
        {
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)) );

            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)) );
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)) );
        }
        else // myRank == 1
        {
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)) );

            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)) );
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)) );
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)) );
        }

        delete stkMeshBulkData;
    }

}

}
