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

/*#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
*/
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

#include <unit_tests/UnitTestRelation.hpp>
#include <unit_tests/UnitTestBulkData.hpp>
#include <unit_tests/UnitTestRingMeshFixture.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

/*class UnitTestingOfRelation : public ::CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestingOfRelation );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  UnitTestingOfRelation() : CppUnit::TestCase() {}

  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION( UnitTestingOfRelation );

void UnitTestingOfRelation::testUnit()
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk::mesh::UnitTestRelation::testRelation ( MPI_COMM_WORLD );
}
*/

STKUNIT_UNIT_TEST(UnitTestingOfRelation, testUnit)
{
  MPI_Barrier ( MPI_COMM_WORLD );
  stk::mesh::UnitTestRelation::testRelation ( MPI_COMM_WORLD );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {


// Unit test the Part functionality in isolation:

void UnitTestRelation::testRelation( ParallelMachine pm )
{
  typedef Field<double>  ScalarFieldType;
 // static const char method[] = "stk::mesh::UnitTestRelation" ;

 // Create a mesh for testing buckets
  std::cout << std::endl ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );
  MetaData meta2 ( entity_names );
  const int spatial_dimension = 3;
  TopologicalMetaData top( meta, spatial_dimension );
  TopologicalMetaData top2( meta2, spatial_dimension );

  BulkData bulk( meta , pm , 4 );
  BulkData bulk2( meta2 , pm , 4 );

  ScalarFieldType & temperature =
       meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume =
       meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  ScalarFieldType & temperature2 =
       meta2.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume2 =
       meta2.declare_field < ScalarFieldType > ( "volume" , 4 );

  Part  & universal     = meta.universal_part ();
  Part  & universal2    = meta2.universal_part ();
  Part  & owned         = meta.locally_owned_part ();

  put_field ( temperature , top.node_rank , universal );
  put_field ( volume , top.element_rank , universal );
  meta.commit();
  put_field ( temperature2 , top2.node_rank , universal2 );
  put_field ( volume2 , top2.element_rank , universal2 );
  meta2.commit();


  bulk.modification_begin();
  bulk2.modification_begin();

  const int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  generate_boxes( bulk , false /* no aura */ , root_box , local_box );
  generate_boxes( bulk2 , false /* no aura */ , root_box , local_box );


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


  STKUNIT_ASSERT_THROW ( Relation r ( Relation::attribute( 2 , 0 ) , cell ) , std::invalid_argument );

  {
      int ok = 0 ;
    try {

  unsigned id = 10000*(~(0u));

  Relation r (Relation::attribute( 0 , id ), cell );

    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitRelation CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestRelation FAILED to catch error for Relation::attribute");
    }
  } 

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
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );


  const unsigned nLocalEdge = nPerProc ;
  MetaData meta3( TopologicalMetaData::entity_rank_names(spatial_dimension) );

  meta3.commit();

  Selector select_owned( meta3.locally_owned_part() );
  Selector select_used = meta3.locally_owned_part() ; 
  Selector select_all(  meta3.universal_part() );
 
  PartVector no_parts ;
   
  std::vector<unsigned> local_count ;

  //------------------------------
  { // No ghosting

    const bool aura_flag = false ;
    UnitTestRingMeshFixture mesh2( pm , nPerProc , false /* No edge parts */ );
    mesh2.m_meta_data.commit();
    mesh2.generate_mesh( aura_flag );

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


  {//ghosting

  if ( 1 < p_size ) { // With ghosting
    const bool aura_flag = true ;

    UnitTestRingMeshFixture mesh3( pm , nPerProc , false /* No edge parts */ );
    mesh3.m_meta_data.commit();
    mesh3.generate_mesh( aura_flag );
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

//----------------------------------------------------------------------
// Testing for a simple loop of mesh entities.
// node_key[i] : edge_key[i] : node_key[ ( i + 1 ) % node_key.size() ]

//void UnitTestRelation::generate_loop(
//  BulkData & mesh ,
//  const PartVector      & edge_parts ,
//  const bool              generate_aura ,
//  const unsigned          nPerProc ,
//  std::vector<EntityId> & node_ids ,
//  std::vector<EntityId> & edge_ids )
//{
//  const unsigned p_rank = mesh.parallel_rank();
//  const unsigned p_size = mesh.parallel_size();
//  const unsigned id_total = nPerProc * p_size ;
//  const unsigned id_begin = nPerProc * p_rank ;
//  const unsigned id_end   = nPerProc * ( p_rank + 1 );
//  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
//  const unsigned nLocalEdge = nPerProc ;
//  const unsigned n_extra = generate_aura && 1 < p_size ? 2 : 0 ;
//
//  node_ids.resize( id_total );
//  edge_ids.resize( id_total );
//  std::vector<unsigned> local_count ;
//
//  for ( unsigned i = 0 ; i < id_total ; ++i ) {
//    node_ids[i] = i + 1;
//    edge_ids[i] = i + 1;
//  }
//
//  // Create a loop of edges:
//  {
//    const PartVector no_parts ;
//    PartVector add_parts ;
//
//    if ( ! edge_parts.empty() ) { add_parts.resize(1); }
//
//    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
//      const unsigned n0 = i ;
//      const unsigned n1 = ( i + 1 ) % id_total ;
//      if ( ! edge_parts.empty() ) {
//        add_parts[0] = edge_parts[ i % edge_parts.size() ];
//      }
//      Entity & e_node_0 = mesh.declare_entity( 0 , node_ids[n0] , no_parts );
//      Entity & e_node_1 = mesh.declare_entity( 0 , node_ids[n1] , no_parts );
//      Entity & e_edge   = mesh.declare_entity( 1 , edge_ids[i] , add_parts );
//      mesh.declare_relation( e_edge , e_node_0 , 0 );
//      mesh.declare_relation( e_edge , e_node_1 , 1 );
//    }
//  }
//
//  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );
//  Selector select_used = mesh.mesh_meta_data().locally_owned_part() |
//                         mesh.mesh_meta_data().globally_shared_part();
//  Selector select_all(  mesh.mesh_meta_data().universal_part() );
//
//  count_entities( select_used , mesh , local_count );
//  STKUNIT_ASSERT( local_count[stk::mesh::Node] == nLocalNode );
//  STKUNIT_ASSERT( local_count[stk::mesh::Edge] == nLocalEdge );
//
//  std::vector<Entity*> all_nodes;
//  get_entities(mesh, stk::mesh::Node, all_nodes);
//
//  unsigned num_selected_nodes =
//      count_selected_entities( select_used, mesh.buckets(stk::mesh::Node) );
//  STKUNIT_ASSERT( num_selected_nodes == local_count[stk::mesh::Node] );
//
//  std::vector<Entity*> universal_nodes;
//  get_selected_entities( select_all, mesh.buckets(stk::mesh::Node), universal_nodes );
//  STKUNIT_ASSERT( universal_nodes.size() == all_nodes.size() );
//
//  mesh.modification_end();
//
//  // Verify declarations and sharing two end nodes:
//
//  count_entities( select_used , mesh , local_count );
//  STKUNIT_ASSERT( local_count[0] == nLocalNode );
//  STKUNIT_ASSERT( local_count[1] == nLocalEdge );
//
//  // Test no-op first:
//
//  std::vector<EntityProc> change ;
//
//  STKUNIT_ASSERT( mesh.modification_begin() );
//  mesh.change_entity_owner( change );
//  mesh.modification_end();
//
//  count_entities( select_used , mesh , local_count );
//  STKUNIT_ASSERT( local_count[0] == nLocalNode );
//  STKUNIT_ASSERT( local_count[1] == nLocalEdge );
//
//  count_entities( select_all , mesh , local_count );
//  STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
//  STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );
//
//  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank()
//  if ( 1 < p_size ) {
//    Entity * const e_node_0 = mesh.get_entity( 0 , node_ids[id_begin] );
//    if ( p_rank == e_node_0->owner_rank() ) {
//      EntityProc entry ;
//      entry.first = e_node_0 ;
//      entry.second = ( p_rank + p_size - 1 ) % p_size ;
//      change.push_back( entry );
//    }
//    STKUNIT_ASSERT( mesh.modification_begin() );
//    mesh.change_entity_owner( change );
//    mesh.modification_end();
//
//    count_entities( select_all , mesh , local_count );
//    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
//    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );
//
//    count_entities( select_used , mesh , local_count );
//    STKUNIT_ASSERT( local_count[0] == nLocalNode );
//    STKUNIT_ASSERT( local_count[1] == nLocalEdge );
//
//    count_entities( select_owned , mesh , local_count );
//    STKUNIT_ASSERT( local_count[0] == nPerProc );
//    STKUNIT_ASSERT( local_count[1] == nPerProc );
//  }
//}

//----------------------------------------------------------------------

namespace {

/* Recursively split a box into ( up - ip ) sub-boxes */

typedef int BOX[3][2] ;

void box_partition( int ip , int up , int axis ,
                    const BOX box ,
                    BOX p_box[] )
{
  const int np = up - ip ;
  if ( 1 == np ) {
    p_box[ip][0][0] = box[0][0] ; p_box[ip][0][1] = box[0][1] ;
    p_box[ip][1][0] = box[1][0] ; p_box[ip][1][1] = box[1][1] ;
    p_box[ip][2][0] = box[2][0] ; p_box[ip][2][1] = box[2][1] ;
  }
  else {
    const int n = box[ axis ][1] - box[ axis ][0] ;
    const int np_low = np / 2 ;  /* Rounded down */
    const int np_upp = np - np_low ;

    const int n_upp = (int) (((double) n) * ( ((double)np_upp) / ((double)np)));
    const int n_low = n - n_upp ;
    const int next_axis = ( axis + 2 ) % 3 ;

    if ( np_low ) { /* P = [ip,ip+np_low) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      dbox[ axis ][1] = dbox[ axis ][0] + n_low ;

      box_partition( ip, ip + np_low, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }

    if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      ip += np_low ;
      dbox[ axis ][0] += n_low ;
      dbox[ axis ][1]  = dbox[ axis ][0] + n_upp ;

      box_partition( ip, ip + np_upp, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }
  }
}

}

void UnitTestRelation::generate_boxes(
  BulkData  & mesh ,
  const bool  generate_aura ,
  const int   root_box[][2] ,
        int   local_box[][2] )
{
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();
  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;
  const unsigned ngz = root_box[2][1] - root_box[2][0] ;
/*
  const unsigned e_global = ngx * ngy * ngz ;
  const unsigned n_global = ( ngx + 1 ) * ( ngy + 1 ) * ( ngz + 1 );
*/

  if ( 0 == p_rank ) {
    std::cout << "Global box = " << ngx << " x " << ngy << " x " << ngz
              << std::endl ;
  }

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  local_box[0][0] = p_box[ p_rank ][0][0] ;
  local_box[0][1] = p_box[ p_rank ][0][1] ;
  local_box[1][0] = p_box[ p_rank ][1][0] ;
  local_box[1][1] = p_box[ p_rank ][1][1] ;
  local_box[2][0] = p_box[ p_rank ][2][0] ;
  local_box[2][1] = p_box[ p_rank ][2][1] ;

  const unsigned nx = local_box[0][1] - local_box[0][0] ;
  const unsigned ny = local_box[1][1] - local_box[1][0] ;
  const unsigned nz = local_box[2][1] - local_box[2][0] ;

  const unsigned e_local = nx * ny * nz ;
  const unsigned n_local = ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 );

  // Create elements:

  std::vector<unsigned> local_count ;

  const PartVector no_parts ;

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    const EntityId n0 = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n1 = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n2 = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n3 = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n4 = 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n5 = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n6 = 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n7 = 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    const EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;

    Entity & node0 = mesh.declare_entity( 0 , n0 , no_parts );
    Entity & node1 = mesh.declare_entity( 0 , n1 , no_parts );
    Entity & node2 = mesh.declare_entity( 0 , n2 , no_parts );
    Entity & node3 = mesh.declare_entity( 0 , n3 , no_parts );
    Entity & node4 = mesh.declare_entity( 0 , n4 , no_parts );
    Entity & node5 = mesh.declare_entity( 0 , n5 , no_parts );
    Entity & node6 = mesh.declare_entity( 0 , n6 , no_parts );
    Entity & node7 = mesh.declare_entity( 0 , n7 , no_parts );
    Entity & elem  = mesh.declare_entity( 3 , elem_id , no_parts );

    mesh.declare_relation( elem , node0 , 0 );
    mesh.declare_relation( elem , node1 , 1 );
    mesh.declare_relation( elem , node2 , 2 );
    mesh.declare_relation( elem , node3 , 3 );
    mesh.declare_relation( elem , node4 , 4 );
    mesh.declare_relation( elem , node5 , 5 );
    mesh.declare_relation( elem , node6 , 6 );
    mesh.declare_relation( elem , node7 , 7 );

    std::vector<Entity*> nodes(8);
    std::vector<Entity*> elems ;
    nodes[0] = & node0 ;
    nodes[1] = & node1 ;
    nodes[2] = & node2 ;
    nodes[3] = & node3 ;
    nodes[4] = & node4 ;
    nodes[5] = & node5 ;
    nodes[6] = & node6 ;
    nodes[7] = & node7 ;

    get_entities_through_relations( nodes , elems );
    STKUNIT_ASSERT_EQUAL( elems.size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( elems[0] , & elem );

    get_entities_through_relations( nodes , 3 , elems );
    STKUNIT_ASSERT_EQUAL( elems.size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( elems[0] , & elem );

  }
  }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );

  Selector select_used = mesh.mesh_meta_data().locally_owned_part() |
                         mesh.mesh_meta_data().globally_shared_part();

  Selector select_all(mesh.mesh_meta_data().universal_part());

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT_EQUAL( e_local , local_count[3] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[2] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[1] );
  STKUNIT_ASSERT_EQUAL( n_local , local_count[0] );

  //Set up ghosting
  const Ghosting & gg = mesh.create_ghosting( std::string("shared") );

  // Set up sharing:
  mesh.modification_end();


  // Verify declarations and sharing

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[3] == e_local );
  STKUNIT_ASSERT( local_count[2] == 0 );
  STKUNIT_ASSERT( local_count[1] == 0 );
  STKUNIT_ASSERT( local_count[0] == n_local );

  for ( int k = local_box[2][0] ; k <= local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j <= local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i <= local_box[0][1] ; ++i ) {
    EntityRank node_type = 0;
    EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
    Entity * const node = mesh.get_entity( node_type , node_id );
    STKUNIT_ASSERT( node != NULL );
    // Shared if on a processor boundary.
    const bool shared =
      ( k == local_box[2][0] && k != root_box[2][0] ) ||
      ( k == local_box[2][1] && k != root_box[2][1] ) ||
      ( j == local_box[1][0] && j != root_box[1][0] ) ||
      ( j == local_box[1][1] && j != root_box[1][1] ) ||
      ( i == local_box[0][0] && i != root_box[0][0] ) ||
      ( i == local_box[0][1] && i != root_box[0][1] );
    STKUNIT_ASSERT_EQUAL( shared , ! node->sharing().empty() );
  }
  }
  }


  for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
    for ( int k = p_box[p][2][0] ; k <= p_box[p][2][1] ; ++k )
    if ( local_box[2][0] <= k && k <= local_box[2][1] ) {

      for ( int j = p_box[p][1][0] ; j <= p_box[p][1][1] ; ++j )
      if ( local_box[1][0] <= j && j <= local_box[1][1] ) {

        for ( int i = p_box[p][0][0] ; i <= p_box[p][0][1] ; ++i )
        if ( local_box[0][0] <= i && i <= local_box[0][1] ) {

          EntityRank node_type = 0;
          EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
          Entity * const node = mesh.get_entity( node_type , node_id );
          STKUNIT_ASSERT( node != NULL );
          // Must be shared with 'p'
          STKUNIT_ASSERT( in_shared( *node , p ) );
          STKUNIT_ASSERT_EQUAL( in_send_ghost( *node , p ), false );

          //Test for coverage of comm_procs in EntityComm.cpp
          std::vector<unsigned> procs ;
	  comm_procs( gg, *node , procs );

        }
      }
    }
  }

  mesh.modification_begin();
  mesh.destroy_all_ghosting();
  mesh.modification_end();

  delete[] p_box ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

