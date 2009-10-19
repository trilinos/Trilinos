
#include <sstream>

#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

#include <unit_tests/UnitTestBucket.hpp>

STKUNIT_UNIT_TEST(UnitTestingOfBucket, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk::mesh::UnitTestBucket::testBucket ( MPI_COMM_WORLD );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {


// Unit test the Part functionality in isolation:

void UnitTestBucket::testBucket( ParallelMachine pm )
{
 // static const char method[] = "stk::mesh::UnitTestBucket" ;

 // Create a mesh for testing buckets
  std::cout << std::endl ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityType" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );
  typedef Field<double>  ScalarFieldType;

  ScalarFieldType & temperature =
       meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume =
       meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  Part  & universal     = meta.universal_part ();
  put_field ( temperature , Node , universal );
  put_field ( volume , Element , universal );
  meta.commit();

  BulkData bulk( meta , pm , 4 );

  const int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  generate_boxes( bulk , false /* no aura */ , root_box , local_box );

  
 //  First, test for streaming IO;
  {
    std::string gold1;
 // Parallel and Serial runs have different part intersections for the first
 // bucket
    if ( bulk.parallel_size() == 1 )
       gold1 = "Bucket( EntityType0 : {UNIVERSAL} {USES} {OWNS} )";
    else
       gold1 = "Bucket( EntityType0 : {UNIVERSAL} )";
    Bucket *b1 = bulk.buckets(0)[0];
    std::stringstream  out1_str;
    out1_str << (*b1);
    STKUNIT_ASSERT ( out1_str.str() == gold1 );

 // Need to validate print against parallel gold string
    std::stringstream  gold2;
    gold2 << gold1 << "\n";
    std::cout << "----->\n";
    print ( std::cout , "- " , *b1 );
    std::cout << "<-----\n";
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
  }

 // Fourth, check has_superset (...)
  {
    STKUNIT_ASSERT_EQUAL ( has_superset ( *bulk.buckets(0)[0] , meta.get_parts() ) , bulk.parallel_size() == 1 );
  }

 // Fifth, check throw_field_data_array (...)
  {
    STKUNIT_ASSERT_THROW ( throw_field_data_array ( *meta.get_fields()[0] , 10 ) , std::runtime_error );
  }

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
      mesh.declare_relation( e_edge , e_node_0 , 0 , 0 );
      mesh.declare_relation( e_edge , e_node_1 , 1 , 0 );
    }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );
  Selector select_used( mesh.mesh_meta_data().locally_used_part() );
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
    const std::vector<EntityProc> & shared = mesh.shared_entities();

    STKUNIT_ASSERT_EQUAL( shared.size() , size_t(2) );

    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    STKUNIT_ASSERT( shared[0].first->identifier() == node_ids[n0] );
    STKUNIT_ASSERT( shared[1].first->identifier() == node_ids[n1] );
    STKUNIT_ASSERT_EQUAL( shared[0].first->sharing().size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( shared[1].first->sharing().size() , size_t(1) );
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

    mesh.declare_relation( elem , node0 , 0 , 0 );
    mesh.declare_relation( elem , node1 , 1 , 0 );
    mesh.declare_relation( elem , node2 , 2 , 0 );
    mesh.declare_relation( elem , node3 , 3 , 0 );
    mesh.declare_relation( elem , node4 , 4 , 0 );
    mesh.declare_relation( elem , node5 , 5 , 0 );
    mesh.declare_relation( elem , node6 , 6 , 0 );
    mesh.declare_relation( elem , node7 , 7 , 0 );
  }
  }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );

  stk::mesh::PartVector part_intersection;
  part_intersection.push_back( &(mesh.mesh_meta_data().locally_used_part()) );
  stk::mesh::PartVector part_union;
  part_union.push_back( &(mesh.mesh_meta_data().locally_used_part()) );

  Selector select_used( part_intersection, part_union );

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
    EntityType node_type = 0;
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

  size_t count_shared_node_pairs = 0 ;
  for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
    for ( int k = p_box[p][2][0] ; k <= p_box[p][2][1] ; ++k )
    if ( local_box[2][0] <= k && k <= local_box[2][1] ) {

      for ( int j = p_box[p][1][0] ; j <= p_box[p][1][1] ; ++j )
      if ( local_box[1][0] <= j && j <= local_box[1][1] ) {

        for ( int i = p_box[p][0][0] ; i <= p_box[p][0][1] ; ++i )
        if ( local_box[0][0] <= i && i <= local_box[0][1] ) {

          EntityType node_type = 0;
          EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
          Entity * const node = mesh.get_entity( node_type , node_id );
          STKUNIT_ASSERT( node != NULL );
          // Must be shared with 'p'
          PairIterEntityProc iter = node->sharing();
          for ( ; ! iter.empty() && iter->second != p ; ++iter );
          STKUNIT_ASSERT( ! iter.empty() );

          ++count_shared_node_pairs ;
        }
      }
    }
  }
  STKUNIT_ASSERT_EQUAL( mesh.shared_entities().size() , count_shared_node_pairs );

  delete[] p_box ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

