/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

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
using stk::mesh::DefaultFEM;
using stk::mesh::Ghosting;
using stk::mesh::fixtures::BoxFixture;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fem::NODE_RANK;

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testRelation)
{
  // Unit test the Part functionality in isolation:

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier ( MPI_COMM_WORLD );

  typedef stk::mesh::Field<double>  ScalarFieldType;
 // static const char method[] = "stk::mesh::UnitTestRelation" ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  unsigned max_bucket_size = 4;

  BoxFixture fixture1(pm , max_bucket_size, entity_names),
             fixture2(pm , max_bucket_size, entity_names);

  MetaData& meta  = fixture1.meta_data();
  MetaData& meta2 = fixture2.meta_data();
  const int spatial_dimension = 3;
  const EntityRank element_rank = stk::mesh::fem::element_rank(fixture1.fem());

  BulkData& bulk  = fixture1.bulk_data();
  BulkData& bulk2 = fixture2.bulk_data();

  ScalarFieldType & temperature =
    meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume =
    meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  ScalarFieldType & temperature2 =
    meta2.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume2 =
    meta2.declare_field < ScalarFieldType > ( "volume" , 4 );

  Part & universal  = meta.universal_part ();
  Part & universal2 = meta2.universal_part ();
  Part & owned      = meta.locally_owned_part ();

  stk::mesh::put_field ( temperature , NODE_RANK , universal );
  stk::mesh::put_field ( volume , element_rank  , universal );
  meta.commit();
  stk::mesh::put_field ( temperature2 , NODE_RANK , universal2 );
  stk::mesh::put_field ( volume2 , element_rank  , universal2 );

  meta2.commit();

  bulk.modification_begin();
  bulk2.modification_begin();

  const int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box1[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  int local_box2[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };

  {
    bulk.modification_begin();
    fixture1.generate_boxes(root_box, local_box1);

    const Ghosting & gg = bulk.create_ghosting( std::string("shared") );

    // Test for coverage of comm_procs in EntityComm.cpp
    EntityVector nodes;
    stk::mesh::get_entities(bulk, NODE_RANK, nodes);
    std::vector<unsigned> procs ;
    STKUNIT_ASSERT(!nodes.empty());
    stk::mesh::comm_procs( gg, *nodes.front() , procs );

    STKUNIT_ASSERT(bulk.modification_end());

    bulk.modification_begin();
    bulk.destroy_all_ghosting();
    STKUNIT_ASSERT(bulk.modification_end());
  }

  {
    bulk2.modification_begin();
    fixture2.generate_boxes(root_box, local_box2);

    bulk2.create_ghosting( std::string("shared") );

    STKUNIT_ASSERT(bulk2.modification_end());

    bulk2.modification_begin();
    bulk2.destroy_all_ghosting();
    STKUNIT_ASSERT(bulk2.modification_end());
  }

  Entity &cell = *(bulk.buckets (3)[0]->begin());
  Entity &node = bulk.buckets (0)[0]-> operator [] ( 0 );
  Entity &nodeb = bulk.buckets (0)[0]-> operator [] ( 2 );

  std::vector<Part *> parts;
  parts.push_back ( &universal );
  parts.push_back ( &owned );
  bulk.modification_begin();
  stk::mesh::EntityId  new_id = bulk.parallel_rank() + 1;
  Entity &edge = bulk.declare_entity ( 1 , new_id , parts );

  Entity &cell2 = *(bulk2.buckets (3)[0]->begin());
  Entity &node2 = *(bulk2.buckets (0)[0]->begin());

  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( node , cell , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( cell , node2 , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( cell2 , node , 0 ) , std::runtime_error );

  bulk.declare_relation ( edge , node , 1 );
  STKUNIT_ASSERT_THROW ( bulk.declare_relation ( edge , nodeb , 1 ) , std::runtime_error );
  bulk.declare_relation ( edge , nodeb , 2 );

  std::stringstream s;
  s << *edge.relations().first ;

  bulk.modification_end();

  //Testing on in_send_ghost and in_shared in EntityComm.cpp
  enum { nPerProc = 10 };
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  const unsigned nLocalEdge = nPerProc ;
  MetaData meta3( stk::mesh::fem::entity_rank_names(spatial_dimension) );

  meta3.commit();

  Selector select_owned( meta3.locally_owned_part() );
  Selector select_used = meta3.locally_owned_part() ;
  Selector select_all(  meta3.universal_part() );

  stk::mesh::PartVector no_parts ;

  std::vector<unsigned> local_count ;

  //------------------------------
  { // No ghosting
    bool aura_flag = false;
    RingFixture mesh2( pm , nPerProc , false /* No edge parts */ );
    mesh2.m_meta_data.commit();

    mesh2.m_bulk_data.modification_begin();
    mesh2.generate_mesh( );
    STKUNIT_ASSERT(stk::unit_test::modification_end_wrapper(mesh2.m_bulk_data,
                                                            aura_flag));
    mesh2.m_bulk_data.modification_begin();
    mesh2.fixup_node_ownership( );
    STKUNIT_ASSERT(stk::unit_test::modification_end_wrapper(mesh2.m_bulk_data,
                                                            aura_flag));

    // This process' first element in the loop
    // if a parallel mesh has a shared node

    Entity * edgenew = mesh2.m_bulk_data.get_entity( 1 , mesh2.m_edge_ids[ nLocalEdge * p_rank ] );

    mesh2.m_bulk_data.modification_begin();
    for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
      STKUNIT_ASSERT_EQUAL( in_shared( *edgenew , p ), false );
      STKUNIT_ASSERT_EQUAL( in_send_ghost( *edgenew , p ), false );
    }

    Entity * edgenew2 = mesh2.m_bulk_data.get_entity( 1 , mesh2.m_edge_ids[ nLocalEdge * p_rank ] );
    STKUNIT_ASSERT_EQUAL( in_send_ghost( *edgenew2 , p_rank+100 ), false );

    Entity * node3 = mesh2.m_bulk_data.get_entity( 0 , mesh2.m_node_ids[ nLocalEdge * p_rank ] );
    STKUNIT_ASSERT_EQUAL( in_shared( *node3 , p_rank+100 ), false );
  }

  { //ghosting

  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh3( pm , nPerProc , false /* No edge parts */ );
    mesh3.m_meta_data.commit();

    mesh3.m_bulk_data.modification_begin();
    mesh3.generate_mesh();
    STKUNIT_ASSERT(mesh3.m_bulk_data.modification_end());

    mesh3.m_bulk_data.modification_begin();
    mesh3.fixup_node_ownership();
    STKUNIT_ASSERT(mesh3.m_bulk_data.modification_end());

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity * node3 = mesh3.m_bulk_data.get_entity( 0 , mesh3.m_node_ids[ nNotOwned ] );
    Entity * node4 = mesh3.m_bulk_data.get_entity( 0 , mesh3.m_node_ids[ nNotOwned ] );


    EntityId node_edge_ids[2] ;
    node_edge_ids[0] = node3->relations()[0].entity()->identifier();
    node_edge_ids[1] = node3->relations()[1].entity()->identifier();

    mesh3.m_bulk_data.modification_begin();

    for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
      //FIXME for Carol the check below did not pass for -np 3 or 4
      //STKUNIT_ASSERT_EQUAL( in_shared( *node3 , p ), true );
      STKUNIT_ASSERT_EQUAL( in_send_ghost( *node3 , p ), false );
    }

    //not owned and not shared
    Entity * node5 = mesh3.m_bulk_data.get_entity( 0 , mesh3.m_node_ids[ nLocalEdge * p_rank ] );

    node_edge_ids[0] = node5->relations()[0].entity()->identifier();
    node_edge_ids[1] = node5->relations()[1].entity()->identifier();

    STKUNIT_ASSERT_EQUAL( in_shared( *node5 , p_rank+100 ), false );
    STKUNIT_ASSERT_EQUAL( in_send_ghost( *node4 , p_rank+100 ), false );
  }

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
  MetaData meta_data(stk::mesh::fem::entity_rank_names(spatial_dim));
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create element
  const EntityRank entity_rank = stk::mesh::fem::element_rank(spatial_dim);
  Entity & elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);

  // Create node
  Entity & node = mesh.declare_entity(NODE_RANK, p_rank+1 /*node_id*/, empty_parts);

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
  MetaData meta_data(stk::mesh::fem::entity_rank_names(spatial_dim));
  meta_data.commit();
  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create element
  const EntityRank entity_rank = stk::mesh::fem::element_rank(spatial_dim);
  Entity & elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);

  // Create node
  Entity & node = mesh.declare_entity(NODE_RANK, p_rank+1 /*node_id*/, empty_parts);

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
  MetaData meta_data(stk::mesh::fem::entity_rank_names(spatial_dim));
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

  Entity* elem_ptr = NULL;
  Entity* edge_ptr = NULL;
  EntityVector nodes;
  const unsigned nodes_per_elem = 4, nodes_per_side = 2;

  if (p_rank < 2) {
    // We're just going to add everything to the universal part
    stk::mesh::PartVector empty_parts;

    // Create element
    const EntityRank entity_rank = stk::mesh::fem::element_rank(spatial_dim);
    Entity & elem = mesh.declare_entity(entity_rank, p_rank+1 /*elem_id*/, empty_parts);
    elem_ptr = &elem;

    // Create nodes
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
      nodes.push_back(&mesh.declare_entity(NODE_RANK, id, empty_parts));
    }

    // Add relations to nodes
    unsigned rel_id = 0;
    for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
      mesh.declare_relation( elem, **itr, rel_id );
    }

    // Create edge
    const EntityRank edge_rank = stk::mesh::fem::side_rank(spatial_dim);
    Entity & edge = mesh.declare_entity(edge_rank, 1 /*id*/, empty_parts);
    edge_ptr = &edge;

    // Set up relation from elem to edge
    mesh.declare_relation( *elem_ptr, *edge_ptr, 0 /*rel-id*/ );
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
      mesh.declare_relation( *edge_ptr, *nodes[node_idx], rel_id );
    }
  }

  mesh.modification_end();
}
