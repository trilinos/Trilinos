/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>
#include <iostream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/fixtures/BoxFixture.hpp>
#include <stk_mesh/fixtures/RingFixture.hpp>

#include <unit_tests/UnitTestModificationEndWrapper.hpp>

#include <Shards_BasicTopologies.hpp>

using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityVector;
using stk::mesh::Part;
using stk::mesh::Relation;
using stk::mesh::Selector;
using stk::mesh::EntityId;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Ghosting;
using stk::mesh::fixtures::BoxFixture;
using stk::mesh::fixtures::HexFixture;
using stk::mesh::fixtures::RingFixture;

namespace {

const EntityRank NODE_RANK = MetaData::NODE_RANK;

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testRelationCoverage)
{
  // Test some relation error cases

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier ( MPI_COMM_WORLD );

  // Just use any random fixture for convenience
  HexFixture fixture(pm, 3 /*x*/, 3 /*y*/, 3 /*z*/);
  MetaData& meta  = fixture.m_fem_meta;
  BulkData& bulk  = fixture.m_bulk_data;

  meta.commit();

  fixture.generate_mesh();

  Entity node0 = (*bulk.buckets(MetaData::NODE_RANK)[0])[0];
  Entity node1 = (*bulk.buckets(MetaData::NODE_RANK)[0])[1];

  bulk.modification_begin();
  std::vector<Part*> empty_parts;
  stk::mesh::EntityId  new_id = bulk.parallel_rank() + 1;
  Entity edge = bulk.declare_entity( MetaData::EDGE_RANK , new_id , empty_parts );

  // Cannot declare a back relation
  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( node0 , edge , 0 /*rel ord*/) , std::runtime_error );

  // Test that you cannot declare relations between entities that are from
  // different meshes.
  {
    HexFixture fixture2(pm, 3 /*x*/, 3 /*y*/, 3 /*z*/);
    MetaData& meta2  = fixture2.m_fem_meta;
    BulkData& bulk2  = fixture2.m_bulk_data;

    meta2.commit();

    fixture2.generate_mesh();

    bulk2.modification_begin();

    Entity edge_from_other_mesh = bulk2.declare_entity( MetaData::EDGE_RANK , new_id , empty_parts );

    // Actual test is here
    STKUNIT_ASSERT_THROW ( bulk.declare_relation ( edge_from_other_mesh , node0 , 0 /*rel ord*/) , std::runtime_error );

    bulk2.modification_end();
  }

  // Test that redeclaration of relation to different node but same ordinal throws
  bulk.declare_relation ( edge , node0 , 1 );
  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( edge , node1 , 1 /*rel ord*/) , std::runtime_error );

  bulk.modification_end();
}

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testRelationNoGhosting)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  const unsigned nPerProc = 10;

  const bool aura_flag = false;
  RingFixture mesh( pm , nPerProc , false /* No element parts */ );
  mesh.m_meta_data.commit();

  BulkData& ring_bulk = mesh.m_bulk_data;

  ring_bulk.modification_begin();
  mesh.generate_mesh( );
  STKUNIT_ASSERT(stk::unit_test::modification_end_wrapper(ring_bulk,
                                                          aura_flag));
  ring_bulk.modification_begin();
  mesh.fixup_node_ownership( );
  STKUNIT_ASSERT(stk::unit_test::modification_end_wrapper(ring_bulk,
                                                          aura_flag));

  // This process' first element in the loop
  // if a parallel mesh has a shared node

  Entity elementnew = ring_bulk.get_entity( MetaData::ELEMENT_RANK , mesh.m_element_ids[ nPerProc * p_rank ] );

  ring_bulk.modification_begin();
  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    if ( p != p_rank ) {
      STKUNIT_ASSERT_EQUAL( ring_bulk.in_shared( elementnew.key() , p ), false );
      STKUNIT_ASSERT_EQUAL( ring_bulk.in_send_ghost( elementnew.key() , p ), false );
    }
  }

  elementnew = ring_bulk.get_entity( MetaData::ELEMENT_RANK , mesh.m_element_ids[ nPerProc * p_rank ] );
  STKUNIT_ASSERT_EQUAL( ring_bulk.in_send_ghost( elementnew.key() , p_rank+100 ), false );

  Entity node = ring_bulk.get_entity( MetaData::NODE_RANK , mesh.m_node_ids[ nPerProc * p_rank ] );
  STKUNIT_ASSERT_EQUAL( ring_bulk.in_shared( node.key() , p_rank+100 ), false );
}

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testRelationWithGhosting)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  const unsigned nPerProc = 10;

  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No element parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh();
    STKUNIT_ASSERT(bulk.modification_end());

    bulk.modification_begin();
    mesh.fixup_node_ownership();
    STKUNIT_ASSERT(bulk.modification_end());

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity node = bulk.get_entity( MetaData::NODE_RANK , mesh.m_node_ids[ nNotOwned ] );

    bulk.modification_begin();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      if ( p != p_rank ) {
        //FIXME for Carol the check below did not pass for -np 3 or 4
        //STKUNIT_ASSERT_EQUAL( in_shared( *node3 , p ), true );
        STKUNIT_ASSERT_EQUAL( bulk.in_send_ghost( node.key() , p ), false );
      }
    }

    //not owned and not shared
    Entity node2 = bulk.get_entity( MetaData::NODE_RANK , mesh.m_node_ids[ nPerProc * p_rank ] );

    STKUNIT_ASSERT_EQUAL( bulk.in_shared( node2.key() , p_rank+100 ), false );
    STKUNIT_ASSERT_EQUAL( bulk.in_send_ghost( node.key() , p_rank+100 ), false );
  }
}

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testDegenerateRelation)
{
  // Test that, if you set up degenerate relations, only of the relations
  // is deleted when you destroy one of the degenerate relations.
  // BulkData::destroy_relation has been changed to take a relation-id so
  // that it can work this way.
  //
  // To test this, we set up an element that has several relations
  // to the same node and then delete them one by one.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create element
  const EntityRank entity_rank = MetaData::ELEMENT_RANK;
  Entity elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);

  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, p_rank+1 /*node_id*/, empty_parts);

  // Add degenerate relations
  const unsigned nodes_per_elem = 4;
  for (unsigned i = 0; i < nodes_per_elem; ++i) {
    mesh.declare_relation( elem, node, i );
  }

  // Elem should have nodes-per-elem relations
  STKUNIT_ASSERT_EQUAL( nodes_per_elem, elem.relations().size() );

  // Destroy relation one-by-one, always checking that appropriate number
  // of relations remain.
  for (unsigned i = 0; i < nodes_per_elem; ++i) {
    mesh.destroy_relation( elem, node, i );
    STKUNIT_ASSERT_EQUAL( nodes_per_elem - (i+1), elem.relations().size() );
  }

  mesh.modification_end();
}

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testRelationAttribute)
{
  // Test relation attribute

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create element
  const EntityRank entity_rank = MetaData::ELEMENT_RANK;
  Entity elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);

  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, p_rank+1 /*node_id*/, empty_parts);

  mesh.declare_relation( elem, node, 0 );

  const Relation & my_relation = *(elem.relations(NODE_RANK).begin());
  my_relation.set_attribute(6u);

  STKUNIT_ASSERT_EQUAL( my_relation.attribute(), 6u);

  mesh.modification_end();
}

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testDoubleDeclareOfRelation)
{
  // It should be legal to declare the same relation between shared
  // entities on two procs.
  //
  // 1---3---5
  // | 1 | 2 |
  // 2---4---6
  //
  // To test this, we use the mesh above, with elem 1 going on rank 0 and
  // elem 2 going on rank 1. Nodes 3,4 are shared along with the edge between
  // nodes 3 and 4. On both procs we declare relations from the shared edge
  // to the shared nodes on both procs.
  //
  // TODO: If we change how declare_relation works, not requiring all
  // sharers to declare the same relations, but instead allowing just
  // the owner to declare relations, that should be tested here.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();
  unsigned p_size = mesh.parallel_size();

  // Bail if we only have one proc
  if (p_size == 1) {
    return;
  }

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  Entity edge = Entity();
  EntityVector nodes;
  const unsigned nodes_per_elem = 4, nodes_per_side = 2;

  if (p_rank < 2) {
    // We're just going to add everything to the universal part
    stk::mesh::PartVector empty_parts;

    // Create element
    const EntityRank entity_rank = MetaData::ELEMENT_RANK;
    Entity elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);

    // Create nodes
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
      nodes.push_back(mesh.declare_entity(NODE_RANK, id, empty_parts));
    }

    // Add relations to nodes
    unsigned rel_id = 0;
    for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
      mesh.declare_relation( elem, *itr, rel_id );
    }

    // Create edge
    const EntityRank edge_rank = meta_data.side_rank();
    edge = mesh.declare_entity(edge_rank, 1 /*id*/, empty_parts);

    // Set up relation from elem to edge
    mesh.declare_relation( elem, edge, 0 /*rel-id*/ );
  }

  mesh.modification_end();

  mesh.modification_begin();

  if (p_rank < 2) {
    // Set up relations from edge to nodes
    unsigned rel_id = 0;
    const unsigned starting_node_idx = (1 - p_rank) * nodes_per_side;
    for (unsigned node_idx = starting_node_idx;
         node_idx < starting_node_idx + nodes_per_side;
         ++node_idx, ++rel_id) {
      mesh.declare_relation( edge, nodes[node_idx], rel_id );
    }
  }

  mesh.modification_end();
}

}
