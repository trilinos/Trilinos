// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc





using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::EntityId;
using stk::mesh::fixtures::RingFixture;

//----------------------------------------------------------------------
// Testing for mesh entities without relations

TEST(UnitTestingOfBulkData, testDestroy_nodes)
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

  const PartVector no_parts ;

  meta.commit();

  BulkData bulk( meta , pm );

  // Ids for all entities (all entities have type 0):

  std::vector<EntityId> ids( id_total );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    ids[i] = i + 1;
  }

  // Declare just those entities in my range of ids:

  ASSERT_TRUE( bulk.modification_begin() );
  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    bulk.declare_entity( stk::topology::NODE_RANK , ids[i] , no_parts );
  }
  ASSERT_TRUE( bulk.modification_end() );

  // Verify that I only have entities in my range:

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity e = bulk.get_entity( stk::topology::NODE_RANK , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      ASSERT_TRUE( bulk.is_valid(e) );
      ASSERT_TRUE( p_rank == bulk.parallel_owner_rank(e) );
    }
    else {
      ASSERT_TRUE( !bulk.is_valid(e) );
    }
  }

  // Delete one entity at a time.

  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    Entity e = bulk.get_entity( stk::topology::NODE_RANK , ids[ i ] );

    ASSERT_TRUE( bulk.is_valid(e) );

    bulk.modification_begin();
    ASSERT_TRUE( bulk.destroy_entity( e ) );
    bulk.modification_end();

    if ( id_begin < i ) {
      ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[ i - 1 ] )) );
    }

    e = bulk.get_entity( stk::topology::NODE_RANK , ids[ i ] );
    ASSERT_TRUE( !bulk.is_valid(e) );
  }
}

//----------------------------------------------------------------------

void assert_is_destroyed(const BulkData& mesh, const Entity entity )
{
  ASSERT_TRUE( !mesh.is_valid(entity) || mesh.bucket(entity).capacity() == 0 );
}

TEST(UnitTestingOfBulkData, testDestroy_ring)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );
  // const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used = meta.locally_owned_part() | meta.globally_shared_part();
  Selector select_all(  meta.universal_part() );

  PartVector no_parts ;

  std::vector<unsigned> local_count ;

  //------------------------------
  { // No ghosting
    const bool aura_flag = false ;

    RingFixture mesh( pm , nPerProc , false /* No element parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk,
                                                           false /*no aura*/));

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk,
                                                           false /*no aura*/));

    // This process' first element in the loop
    // if a parallel mesh has a shared node
    Entity element = bulk.get_entity( stk::topology::ELEMENT_RANK , mesh.m_element_ids[ nLocalElement * p_rank ] );
    Entity const *element_nodes = bulk.begin_nodes(element);
    Entity node0 = element_nodes[0];
    Entity node1 = element_nodes[1];

    const size_t node0_elements = bulk.count_relations(node0);
    const size_t node1_elements = bulk.count_relations(node1);

    ASSERT_TRUE( 1 <= node0_elements && node0_elements <= 2 );
    ASSERT_TRUE( 1 <= node1_elements && node1_elements <= 2 );

    ASSERT_TRUE( bulk.begin_elements(node0)[0] == element ||
                    bulk.begin_elements(node0)[1] == element );

    ASSERT_TRUE( bulk.begin_elements(node1)[0] == element ||
                    bulk.begin_elements(node1)[1] == element );

    bulk.modification_begin();

    // Destroy the element:
    bool result = bulk.destroy_entity( element );
    element = Entity();
    ASSERT_TRUE( true == result );

    // Destroy orphanned node:
    if ( bulk.has_no_relations(node0) ) {
      ASSERT_TRUE( bulk.destroy_entity( node0 ) );
      node0 = Entity();
    }
    if ( bulk.has_no_relations(node1) ) {
      ASSERT_TRUE( bulk.destroy_entity( node1 ) );
      node1 = Entity();
    }
    ASSERT_TRUE( stk::unit_test::modification_end_wrapper(bulk, aura_flag) );

    if ( bulk.is_valid(node0) ) {
      ASSERT_EQ( node0_elements - 1 , bulk.count_relations(node0) );
    }
    if ( bulk.is_valid(node1) ) {
      ASSERT_EQ( node1_elements - 1 , bulk.count_relations(node1) );
    }
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No element parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    ASSERT_TRUE( bulk.modification_end() );

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    ASSERT_TRUE( bulk.modification_end() );

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity node = bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nNotOwned ] );

    ASSERT_TRUE( bulk.is_valid(node) );
    ASSERT_NE( p_rank , bulk.parallel_owner_rank(node) );

    ASSERT_EQ( size_t(1) , bulk.entity_comm_map_shared(bulk.entity_key(node)).size() );
    ASSERT_EQ( size_t(2) , bulk.count_relations(node) );

    EntityId node_element_ids[2] ;
    Entity const *node_elems = bulk.begin_elements(node);
    node_element_ids[0] = bulk.identifier(node_elems[0]);
    node_element_ids[1] = bulk.identifier(node_elems[1]);

    bulk.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached elements.
    // One will be owned and the other ghosted

    const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; )
    {
      --irank;
      if (bulk.connectivity_map().valid(stk::topology::NODE_RANK, irank)) {
        stk::mesh::Entity const *to_b = bulk.begin(node, irank);
        stk::mesh::Entity const *to_e = bulk.end(node, irank);
        for (; to_b != to_e;
             to_b = bulk.begin(node, irank), to_e = bulk.end(node, irank))
        {
          ASSERT_TRUE( bulk.destroy_entity(*(to_e -1)) );
        }
      }
    }
    ASSERT_TRUE( bulk.destroy_entity( node ) );

    ASSERT_TRUE( bulk.modification_end() );

    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::NODE_RANK, mesh.m_node_ids[nNotOwned] ) );
    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::ELEMENT_RANK, node_element_ids[0] ) );
    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::ELEMENT_RANK, node_element_ids[1] ) );

    // assert that no entities are shared or ghosted
    ASSERT_TRUE( bulk.comm_list().empty() );
  }
  //------------------------------
  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No element parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh( );
    ASSERT_TRUE( bulk.modification_end() );

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    ASSERT_TRUE( bulk.modification_end() );

    // The owned shared entity:
    const unsigned nOwned = ( nPerProc * ( p_rank + 1 ) ) % mesh.m_node_ids.size();
    const unsigned nNotOwned = nPerProc * p_rank ;

    Entity node_owned = bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nOwned ] );
    Entity node_not_owned = bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nNotOwned ] );

    ASSERT_TRUE( bulk.is_valid(node_owned) );
    ASSERT_TRUE( bulk.is_valid(node_not_owned));
    ASSERT_NE( p_rank , bulk.parallel_owner_rank(node_not_owned) );
    ASSERT_EQ( p_rank , bulk.parallel_owner_rank(node_owned) );

    ASSERT_EQ( 1u , bulk.entity_comm_map_shared(bulk.entity_key(node_owned)).size() );
    ASSERT_EQ( 1u , bulk.entity_comm_map_shared(bulk.entity_key(node_not_owned)).size() );
    ASSERT_EQ( 2u , bulk.count_relations(node_owned) );

    EntityId node_element_ids[2] ;
    stk::mesh::Entity const *node_elems = bulk.begin_elements(node_owned);
    node_element_ids[0] = bulk.identifier(node_elems[0]);
    node_element_ids[1] = bulk.identifier(node_elems[1]);

    bulk.modification_begin();

    // This process' first node in the loop is shared, destroy it
    // First have to destroy attached elements.
    // One will be owned and the other ghosted
    const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; )
    {
      --irank;
      if (bulk.connectivity_map().valid(stk::topology::NODE_RANK, irank)) {
        stk::mesh::Entity const *to_b = bulk.begin(node_owned, irank);
        stk::mesh::Entity const *to_e = bulk.end(node_owned, irank);
        for ( ; to_b != to_e;
              bulk.begin(node_owned, irank), to_e = bulk.end(node_owned, irank))
        {
          ASSERT_TRUE( bulk.destroy_entity(*(to_e - 1)) );
        }
      }
    }
    ASSERT_TRUE( bulk.destroy_entity( node_owned ) );

    ASSERT_TRUE( bulk.modification_end() );

    // Ownership of the other process' owned, shared, and destroyed node
    // has been transferred to this process.

    ASSERT_EQ( p_rank , bulk.parallel_owner_rank(node_not_owned) );
    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::NODE_RANK, mesh.m_node_ids[ nOwned ] ) );
    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::ELEMENT_RANK, node_element_ids[0] ) );
    assert_is_destroyed(bulk, bulk.get_entity(stk::topology::ELEMENT_RANK, node_element_ids[1] ) );

    // assert that no entities are shared or ghosted
    ASSERT_TRUE( bulk.comm_list().empty() );
  }
}

