/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <unit_tests/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <unit_tests/UnitTestBucket.hpp>
#include <unit_tests/UnitTestMesh.hpp>
#include <unit_tests/UnitTestBoxMeshFixture.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Transaction.hpp>
#include <stk_mesh/baseImpl/BucketImpl.hpp>


STKUNIT_UNIT_TEST(UnitTestingOfBucket, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk::mesh::UnitTestBucket::testBucket ( MPI_COMM_WORLD );
  stk::mesh::UnitTestBucket::testTopologyHelpers( MPI_COMM_WORLD );
  stk::mesh::UnitTestBucket::test_get_involved_parts( MPI_COMM_WORLD );
  stk::mesh::UnitTestBucket::testBucket2( MPI_COMM_WORLD );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {


// Unit test the Part functionality in isolation:

void UnitTestBucket::testBucket( ParallelMachine pm )
{
  typedef Field<double>  ScalarFieldType;
 // static const char method[] = "stk::mesh::UnitTestBucket" ;

 // Create a mesh for testing buckets
  std::cout << std::endl ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );
  BulkData bulk( meta , pm , 4 );

  ScalarFieldType & temperature =
       meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume =
       meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  Part  & universal     = meta.universal_part ();
  put_field ( temperature , Node , universal );
  put_field ( volume , Element , universal );
  meta.commit();


  const int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };

  bulk.modification_begin();
  generate_boxes( bulk , false /* no aura */ , root_box , local_box );
  bulk.modification_end();


 //  First, test for streaming IO;
  {
    std::string gold1;
 // Parallel and Serial runs have different part intersections for the first
 // bucket
    if ( bulk.parallel_size() == 1 )
       gold1 = "Bucket( EntityRank0 : {UNIVERSAL} {OWNS} )";
    else
       gold1 = "Bucket( EntityRank0 : {UNIVERSAL} )";
    Bucket *b1 = bulk.buckets(0)[0];
    std::stringstream  out1_str;
    out1_str << (*b1);
    bool result = out1_str.str() == gold1;
    STKUNIT_ASSERT_EQUAL ( result , true );

 // Need to validate print against parallel gold string
    std::stringstream  gold2;
    gold2 << gold1 << "\n";
  }

 // Second, update state of bucket until circular cue is filled
  {
   /* Need to set some data in state, rotate look for it, rotate 3 more times
      and look for it again */
    for ( size_t i = 0 ; i != 10 ; ++i )
      bulk.update_field_data_states ();
  }

 // Third, checking field_data_valid (...)
  {

    const std::vector< FieldBase * > &field_bases = meta.get_fields();
    STKUNIT_ASSERT_THROW(field_data_valid ( *field_bases[0] , *bulk.buckets(3)[0] , 1 , "error" ) , std::runtime_error);
    STKUNIT_ASSERT_EQUAL(field_data_valid ( *field_bases[0] , *bulk.buckets(0)[0] , 1 , "no_error" ) , true);
    STKUNIT_ASSERT_THROW(field_data_valid ( *field_bases[0] , *bulk.buckets(3)[0] , 99 , "error" ) , std::runtime_error);


    stk::mesh::MetaData meta2 ( entity_names );
    BulkData bulk2( meta2 , pm , 4 );

    ScalarFieldType & temperature2 =
       meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
    ScalarFieldType & volume2 =
       meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
    Part  & universal     = meta2.universal_part ();
    put_field ( temperature2 , Node , universal );
    put_field ( volume2 , Element , universal );
    meta2.commit();

    //Cover line containing messsage for wrong MetaData used
    const std::vector< FieldBase * > &field_bases2 = meta2.get_fields();
    STKUNIT_ASSERT_THROW(field_data_valid ( *field_bases2[0] , *bulk.buckets(0)[0] , 1 , "error" ) , std::runtime_error);

  }

 // Fourth, check has_superset (...) and membership functions
  {
    PartVector tmp(2) ;
    tmp[0] = & meta.universal_part();
    tmp[1] = & meta.locally_owned_part();
    STKUNIT_ASSERT_EQUAL ( has_superset ( *bulk.buckets(0)[0] , tmp ) , bulk.parallel_size() == 1 );
    STKUNIT_ASSERT ( bulk.buckets(0)[0]->member_any ( tmp ) );
    STKUNIT_ASSERT_EQUAL ( bulk.buckets(0)[0]->member_all ( tmp ) , bulk.parallel_size() == 1 );
    STKUNIT_ASSERT ( bulk.buckets(0)[0]->member ( **meta.get_parts().begin() ) );
  }

 // Fifth, check throw_field_data_array (...)
  {
    STKUNIT_ASSERT_THROW ( throw_field_data_array ( *meta.get_fields()[0] , 10 ) , std::runtime_error );
  }
}

void UnitTestBucket::testTopologyHelpers( ParallelMachine pm )
{
  // Tests to complete coverage of TopologyHelpers.cpp - C.Brickley - 12 May 2010

  stk::mesh::MetaData meta ( stk::unit_test::get_entity_rank_names ( 3 ) );

  stk::mesh::Part &new_part4 = meta.declare_part ( "another part");

  stk::mesh::Part & partLeft_1 = meta.declare_part( "block_left_1", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_1 );

  stk::mesh::Part & partLeft_2 = meta.declare_part( "block_left_2", Element );
  stk::mesh::set_cell_topology< shards::Wedge<15>  >( partLeft_2 );

  stk::mesh::Part & partLeft_3 = meta.declare_part( "block_left_3", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_3 );

  stk::mesh::Part & part_A_3 = meta.declare_part( "A_3", Element);
  stk::mesh::set_cell_topology< shards::Quadrilateral<4>  >( part_A_3 );

  stk::mesh::Part & part_B_3 = meta.declare_part( "B_3", Element);

  meta.commit();

  stk::mesh::BulkData bulk ( meta , pm , 100 );
  std::vector<stk::mesh::Part *>  add_part4;
  add_part4.push_back ( &new_part4 );

  int  size , rank;
  rank = stk::parallel_machine_rank( pm );
  size = stk::parallel_machine_size( pm );
  PartVector tmp(1);

  bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 93 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( 0 , new_id+1 , add_part4 );
  }

  int new_id = size * (++id_base) + rank;
  Entity & elem  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & elem1  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & elem2  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & elem3  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & elem4  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & element2  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  new_id = size * (++id_base) + rank;
  Entity & element3  = bulk.declare_entity( 3 , new_id+1 , add_part4 );

  // More coverage obtained for assert_entity_owner in BulkData.cpp:
  /*
  stk::mesh::Entity &n = *bulk.buckets(1)[0]->begin();
  {
    int ok = 0 ;
    try {
      bulk.change_entity_parts ( n , add_part4 );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for change_entity_parts");
    }
  }
  */

  // Coverage for get_cell_topology in TopologyHelpers.cpp

  tmp[0] = & part_A_3;
  bulk.change_entity_parts ( elem1 , tmp );
  bulk.change_entity_parts ( elem1 , tmp );
  tmp[0] = & part_B_3;
  bulk.change_entity_parts ( elem1 , tmp );

  const CellTopologyData * elem1_topology = stk::mesh::get_cell_topology( elem1 );

  {
    bool result = true;
    if (elem1_topology != shards::getCellTopologyData< shards::Wedge<15> >()) {
      result = false;
    }

    STKUNIT_ASSERT_EQUAL(result, false);
  }

  // Coverage for get_cell_topology in TopologyHelpers.cpp; (FAILED WITH MULTIPLE LOCAL TOPOLOGIES)

  {
    int ok = 0 ;
    try {
      //assign 3 different parts with different topologies:

      tmp[0] = & part_A_3;
      bulk.change_entity_parts ( elem , tmp );

      tmp[0] = & partLeft_1;
      bulk.change_entity_parts ( elem , tmp );

      tmp[0] = & partLeft_2;
      bulk.change_entity_parts ( elem , tmp );

      const CellTopologyData * elem_topology = stk::mesh::get_cell_topology( elem );

      {
        bool result = true;
        if (elem_topology != shards::getCellTopologyData< shards::Triangle<3> >()) {
          result = false;
        }

        STKUNIT_ASSERT_EQUAL(result, false);
      }
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for Element, elem");
    }
  }

  // Coverage for get_adjacent_entities in TopologyHelpers.cpp (get_adjacent_entities(..))
  // Element, elem2, has NULL topology

  std::vector<EntitySideComponent> adjacent_entities;
  get_adjacent_entities( elem2 , 3, 1, adjacent_entities);

  // More Coverage for get_adjacent_entities in TopologyHelpers.cpp
  // Element has Quadrilateral topology

  {
    int ok = 0 ;
    try {
      //assign a Quadrilateral part
      tmp[0] = & part_A_3;
      bulk.change_entity_parts ( elem3 , tmp );

      std::vector<EntitySideComponent> adjacent_entities2;

      //use 4 as subcell_rank - an incorrect value
      get_adjacent_entities( elem3 , 4, 0, adjacent_entities2);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for declare_element_side, elem2");
    }
  }

  // More Coverage for get_adjacent_entities in TopologyHelpers.cpp
  // Element has Quadrilateral topology

  {
    int ok = 0 ;
    try {
      std::vector<EntitySideComponent> adjacent_entities3;

      //use 8 as subcell_identifier - an incorrect value
      get_adjacent_entities( elem3 , 1, 8, adjacent_entities3);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for declare_element_side, elem2");
    }
  }

  // Coverage for declare_element_side - line 129 in TopologyHelpers.cpp - "Cannot discern element topology"

  {
    int ok = 0 ;
    try {
      new_id = size * (++id_base) + rank;
      stk::mesh::Entity &face = stk::mesh::declare_element_side( bulk, 3, elem4, new_id+1, &partLeft_2);
      stk::mesh::PairIterRelation rel = face.relations(stk::mesh::Node);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for declare_element_side, elem4");
    }
  }

  // Go all way the through declare_element_side - use new element

  stk::mesh::EntityId elem_node[4];

  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3 ;
  elem_node[3] = 4 ;

  stk::mesh::EntityId elem_id(size * (++id_base) + rank + 1);
  stk::mesh::Entity& element = stk::mesh::declare_element(bulk, part_A_3, elem_id, elem_node );

  new_id = size * (++id_base) + rank;
  stk::mesh::Entity &face2 = stk::mesh::declare_element_side( bulk, new_id+1, element, 0);

  stk::mesh::PairIterRelation rel2 = face2.relations(stk::mesh::Node);

  // Coverage of element_side_polarity in TopologyHelpers.cpp 168-181 and 200-215

  bool state = false;

  {
    bool result = true;
    if (state != element_side_polarity( element, face2, 0) ) {
      result = false;
    }

    STKUNIT_ASSERT_EQUAL(result, false);
  }

  // Coverage of element_side_polarity in TopologyHelpers.cpp

  {
    int ok = 0 ;
    try {
      tmp[0] = & part_A_3;
      bulk.change_entity_parts ( element2 , tmp );

      element_side_polarity( element2, face2, -1);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for element_side_polarity");
    }
  }

  // Coverage of element_side_polarity in TopologyHelpers.cpp - NULL = elem_top

  {
    int ok = 0 ;
    try {
      tmp[0] = & part_B_3;
      bulk.change_entity_parts ( element3 , tmp );

      element_side_polarity( element3, face2, 0);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for get_cell_topology for element_side_polarity");
    }
  }

  //Coverage of TopologyHelpers element_local_side_id(const Entity & elem, const Entity & side )

  tmp[0] = & part_A_3;
  bulk.change_entity_parts ( element3 , tmp );

  int localSide;

  localSide = element_local_side_id ( element3, face2);

  STKUNIT_ASSERT_EQUAL(localSide, -1);

  //Checks line 310-311 in TopologyHelpers element_local_side_id(const Entity & elem, const std::vector<Entity*>& entity_nodes)

  const CellTopologyData* celltopology = get_cell_topology(element);

  unsigned subcell_rank = 1;
  unsigned subcell_identifier = celltopology->subcell_count[subcell_rank] - 1;

  unsigned num_nodes =  celltopology->subcell[subcell_rank][subcell_identifier].topology->node_count;

  PairIterRelation relations = face2.relations(Node);

  std::vector<Entity*> node_entities;

  for (unsigned itr = 0; itr < num_nodes; ++itr) {
    node_entities.push_back(relations[itr].entity());
  }

  /*
  localSide = element_local_side_id ( face2, node_entities); //uninitialized memory access within node_entities!!

  STKUNIT_ASSERT_EQUAL(localSide, -1);
  */

  //Checks line 337 in TopologyHelpers for element_local_side_id(const Entity & elem, const std::vector<Entity*>& entity_nodes)

  std::vector<Entity*> node_entities2;

  localSide = element_local_side_id ( element, node_entities2);

  STKUNIT_ASSERT_EQUAL(localSide, -1);

  bulk.modification_end();
}

//----------------------------------------------------------------------
// Testing for a simple loop of mesh entities.
// node_key[i] : edge_key[i] : node_key[ ( i + 1 ) % node_key.size() ]

void UnitTestBucket::generate_loop(
  BulkData & mesh ,
  const PartVector      & edge_parts ,
  const bool              generate_aura ,
  const unsigned          nPerProc ,
  std::vector<EntityId> & node_ids ,
  std::vector<EntityId> & edge_ids )
{
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;
  const unsigned n_extra = generate_aura && 1 < p_size ? 2 : 0 ;

  node_ids.resize( id_total );
  edge_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    node_ids[i] = i + 1;
    edge_ids[i] = i + 1;
  }

  // Create a loop of edges:
  {
    const PartVector no_parts ;
    PartVector add_parts ;

    if ( ! edge_parts.empty() ) { add_parts.resize(1); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! edge_parts.empty() ) {
        add_parts[0] = edge_parts[ i % edge_parts.size() ];
      }
      Entity & e_node_0 = mesh.declare_entity( 0 , node_ids[n0] , no_parts );
      Entity & e_node_1 = mesh.declare_entity( 0 , node_ids[n1] , no_parts );
      Entity & e_edge   = mesh.declare_entity( 1 , edge_ids[i] , add_parts );
      mesh.declare_relation( e_edge , e_node_0 , 0 );
      mesh.declare_relation( e_edge , e_node_1 , 1 );
    }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );
  Selector select_used = select_owned | mesh.mesh_meta_data().globally_shared_part();
  Selector select_all(  mesh.mesh_meta_data().universal_part() );

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[stk::mesh::Node] == nLocalNode );
  STKUNIT_ASSERT( local_count[stk::mesh::Edge] == nLocalEdge );

  std::vector<Entity*> all_nodes;
  get_entities(mesh, stk::mesh::Node, all_nodes);

  unsigned num_selected_nodes =
      count_selected_entities( select_used, mesh.buckets(stk::mesh::Node) );
  STKUNIT_ASSERT( num_selected_nodes == local_count[stk::mesh::Node] );

  std::vector<Entity*> universal_nodes;
  get_selected_entities( select_all, mesh.buckets(stk::mesh::Node), universal_nodes );
  STKUNIT_ASSERT( universal_nodes.size() == all_nodes.size() );

  mesh.modification_end();

  // Verify declarations and sharing two end nodes:

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  if ( 1 < p_size ) {
    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    Entity * const node0 = mesh.get_entity( Node , node_ids[n0] );
    Entity * const node1 = mesh.get_entity( Node , node_ids[n1] );

    STKUNIT_ASSERT( node0 != NULL );
    STKUNIT_ASSERT( node1 != NULL );

    STKUNIT_ASSERT_EQUAL( node0->sharing().size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( node1->sharing().size() , size_t(1) );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( mesh.modification_begin() );
  mesh.change_entity_owner( change );
  mesh.modification_end();

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  count_entities( select_all , mesh , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank()
  if ( 1 < p_size ) {
    Entity * const e_node_0 = mesh.get_entity( 0 , node_ids[id_begin] );
    if ( p_rank == e_node_0->owner_rank() ) {
      EntityProc entry ;
      entry.first = e_node_0 ;
      entry.second = ( p_rank + p_size - 1 ) % p_size ;
      change.push_back( entry );
    }
    STKUNIT_ASSERT( mesh.modification_begin() );
    mesh.change_entity_owner( change );
    mesh.modification_end();

    count_entities( select_all , mesh , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    count_entities( select_used , mesh , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_owned , mesh , local_count );
    STKUNIT_ASSERT( local_count[0] == nPerProc );
    STKUNIT_ASSERT( local_count[1] == nPerProc );
  }
}

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

void UnitTestBucket::generate_boxes(
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
  }
  }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );
  Selector select_used = select_owned | mesh.mesh_meta_data().globally_shared_part();
  Selector select_all(  mesh.mesh_meta_data().universal_part() );

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT_EQUAL( e_local , local_count[3] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[2] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[1] );
  STKUNIT_ASSERT_EQUAL( n_local , local_count[0] );

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
    if (mesh.parallel_size() > 1) {
      STKUNIT_ASSERT_EQUAL( shared , ! node->sharing().empty() );
    }
  }
  }
  }

  size_t count_shared_node_pairs = 0 ;
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
          PairIterEntityComm iter = node->sharing();
          for ( ; ! iter.empty() && iter->proc != p ; ++iter );
          STKUNIT_ASSERT( ! iter.empty() );

          ++count_shared_node_pairs ;
        }
      }
    }
  }

  size_t count_shared_entities = 0 ;
  for ( std::vector<Entity*>::const_iterator
        i = mesh.entity_comm().begin() ; i != mesh.entity_comm().end() ; ++i ) {
    const PairIterEntityComm ec = (**i).sharing();
    count_shared_entities += ec.size();
  }
  STKUNIT_ASSERT_EQUAL( count_shared_entities , count_shared_node_pairs );

  delete[] p_box ;
}


void UnitTestBucket::test_get_involved_parts(ParallelMachine pm)
{

  // Tests to cover get_involved_parts for GetBuckets.cpp - C.Brickley - 12 May 2010

  stk::mesh::MetaData meta ( stk::unit_test::get_entity_rank_names ( 3 ) );

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  stk::mesh::Part & partLeft_1 = meta.declare_part( "block_left_1", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_1 );

  stk::mesh::Part & partLeft_2 = meta.declare_part( "block_left_2", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_2 );

  stk::mesh::Part & partLeft_3 = meta.declare_part( "block_left_3", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_3 );

  meta.commit();

  PartVector union_parts;
  union_parts.push_back(&partLeft_1);
  union_parts.push_back(&partLeft_2);

  BulkData bulk( meta , pm , 100 );
  PartVector add_part4, no_part;
  add_part4.push_back ( &partLeft_1 );

  bulk.modification_begin();
  int  size , rank;
  rank = stk::parallel_machine_rank( pm );
  size = stk::parallel_machine_size( pm );

  for ( int id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank + 1;
    bulk.declare_entity( 3 , new_id , add_part4 );
    bulk.declare_entity( Node , new_id , no_part );
  }

  bulk.modification_end();

  const std::vector<Bucket*> & buckets = bulk.buckets( Element );

  std::vector<Bucket*>::const_iterator k;

  k = buckets.begin();

  //test 1 covers aecond section of "if" statement in while loop
  stk::mesh::get_involved_parts( union_parts, **k, involved_parts);

  //test 2 covers union_parts.size() = 0
  PartVector union_parts2(0) ;
  stk::mesh::get_involved_parts( union_parts2, **k, involved_parts);

  //test 3 covers first section of "if" statement in while loop
  const std::vector<Bucket*> & buckets2 = bulk.buckets( Node );
  std::vector<Bucket*>::const_iterator k2;

  k2 = buckets2.begin();
  stk::mesh::get_involved_parts( union_parts, **k2, involved_parts);

  // tests on throw_error and BucketIterator in bucket.cpp/hpp

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }
  typedef Field<double>  ScalarFieldType;

  stk::mesh::MetaData meta2 ( entity_names );
  BulkData bulk2( meta2 , pm , 4 );

  ScalarFieldType & temperature2 =
     meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
  ScalarFieldType & volume2 =
     meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
  Part  & universal     = meta2.universal_part ();
  put_field ( temperature2 , Node , universal );
  put_field ( volume2 , Element , universal );
  meta2.commit();

  bulk2.modification_begin();
  bulk2.declare_entity( Edge , rank+1 , no_part );
  bulk2.modification_end();

  const std::vector<Bucket*> & buckets3 = bulk2.buckets( Edge );

  std::vector<Bucket*>::const_iterator k3;

  k3 = buckets3.begin();

  stk::mesh::Bucket& b3 = **k3;
  stk::mesh::BucketIterator bitr3 = b3.begin();

  stk::mesh::Bucket& b2 = **k2;
  stk::mesh::BucketIterator bitr2 = b2.begin();

  //tests operator != given iterator from different bucket - bucket.hpp

  {
    int ok = 0 ;
    try {

      if ( bitr2  !=  bitr3 ){
        // bitr3.throw_error("is NULL") ;
      }

    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for throw_error");
    }
  }

  //tests operator - given iterator from different bucket - bucket.hpp
  {
    int ok = 0 ;
    try {

    const ptrdiff_t n = bitr2 - bitr3 ;

      if ( n  !=  0 ){
        // bitr3.throw_error("is NULL") ;
      }
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for iterator from a different bucket");
    }
  }

}


void UnitTestBucket::testBucket2(ParallelMachine pm)
{

  // Tests to cover print, has_superset and BucketLess::operator() for Buckets.cpp - C.Brickley - 2nd June 2010

  stk::mesh::MetaData meta ( stk::unit_test::get_entity_rank_names ( 3 ) );

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  stk::mesh::Part & partLeft_1 = meta.declare_part( "block_left_1", Element );

  stk::mesh::Part & partLeft_2 = meta.declare_part( "block_left_2", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_2 );

  stk::mesh::Part & partLeft_3 = meta.declare_part( "block_left_3", Element );
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( partLeft_3 );

  meta.commit();

  BulkData bulk( meta , pm , 100 );
  std::vector<stk::mesh::Part *>  add_part4;
  add_part4.push_back ( &partLeft_1 );

  bulk.modification_begin();
  int  size , rank;
  rank = stk::parallel_machine_rank( pm );
  size = stk::parallel_machine_size( pm );

  for ( int id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( 3 , new_id+1 , add_part4 );
  }

  bulk.modification_end();

  const std::vector<Bucket*> & buckets2 = bulk.buckets( Element );

  std::vector<Bucket*>::const_iterator k2;

  k2 = buckets2.begin();

  stk::mesh::Bucket& b2 = **k2;
  stk::mesh::BucketIterator bitr2 = b2.begin();

  //define a new meta and bulkdata
  std::vector<std::string> entity_names(10);

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  typedef Field<double>  ScalarFieldType;

  stk::mesh::MetaData meta2 ( entity_names );
  BulkData bulk2( meta2 , pm , 4 );

  ScalarFieldType & temperature2 =
       meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
  ScalarFieldType & volume2 =
       meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
  Part  & universal     = meta2.universal_part ();
  put_field ( temperature2 , Node , universal );
  put_field ( volume2 , Element , universal );

  typedef Field<double>  VectorFieldType;
  typedef Field<double>  ElementNodePointerFieldType;

  meta2.commit();

  //Test to cover print function in Bucket.cpp
  std::cout << std::endl << "Bucket test" << std::endl ;
  stk::mesh::print(std::cout, "  ", b2);

  //Test to cover has_superset function in Bucket.cpp
  STKUNIT_ASSERT_EQUAL ( has_superset ( b2 , partLeft_3 ) , false );

  //Test on BucketLess::operator() in Bucket.cpp/hpp

  enum { KEY_TMP_BUFFER_SIZE = 64 };

  const unsigned max = ~(0u);

  unsigned key_tmp_buffer[ KEY_TMP_BUFFER_SIZE ];

  std::vector<unsigned> key_tmp_vector ;

  const unsigned key_size = 2 + 3 ;

  unsigned * const key =
    ( key_size <= KEY_TMP_BUFFER_SIZE )
    ? key_tmp_buffer
    : ( key_tmp_vector.resize( key_size ) , & key_tmp_vector[0] );


  key[ key[0] = 3 + 1 ] = max;

  {
    unsigned * const k = key + 1 ;
    for ( unsigned i = 0 ; i < 3 ; ++i ) { k[i] = 1 ; }
  }

  /*
  impl::BucketImpl::last_bucket_in_family( *k2 );

  const unsigned * t = key;
  const Bucket * u = last_bucket;

  BucketLess Buck;

  bool res = Buck(  &t[0], &u[0] );

  STKUNIT_EXPECT_EQUAL( res, false );
  */
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk
