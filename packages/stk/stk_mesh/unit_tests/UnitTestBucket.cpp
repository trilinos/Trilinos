/*------------------------------------------------------------------------*/
/*         _        Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
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

#include <unit_tests/UnitTestBucket.hpp>
#include <unit_tests/UnitTestMesh.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Transaction.hpp>
#include <stk_mesh/baseImpl/BucketImpl.hpp>
#include <stk_mesh/base/Ghosting.hpp>

using stk::unit_test::UnitTestBucket;

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;
using stk::mesh::PairIterEntityComm;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::BucketIterator;
using stk::mesh::Selector;
using stk::mesh::Field;
using stk::mesh::FieldBase;
using stk::mesh::put_field;
using stk::mesh::TopologicalMetaData;

using stk::ParallelMachine;
using std::cout;
using std::endl;


STKUNIT_UNIT_TEST(UnitTestingOfBucket, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  UnitTestBucket::testBucket ( MPI_COMM_WORLD );
  UnitTestBucket::test_get_involved_parts( MPI_COMM_WORLD );
  UnitTestBucket::testBucket2( MPI_COMM_WORLD );
  UnitTestBucket::test_EntityComm( MPI_COMM_WORLD );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------



// Unit test the Part functionality in isolation:

void UnitTestBucket::testBucket( ParallelMachine pm )
{
  typedef Field<double>  ScalarFieldType;
 // static const char method[] = "stk::mesh::UnitTestBucket" ;

 // Create a mesh for testing buckets
  cout << endl ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );
  BulkData bulk( meta , pm , 4 );
  const int spatial_dimension = 3;
  TopologicalMetaData top(meta,spatial_dimension);

  ScalarFieldType & temperature =
       meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  ScalarFieldType & volume =
       meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  Part  & universal     = meta.universal_part ();
  put_field ( temperature , top.node_rank , universal );
  put_field ( volume , top.element_rank , universal );
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


    MetaData meta2 ( entity_names );
    BulkData bulk2( meta2 , pm , 4 );

    ScalarFieldType & temperature2 =
       meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
    ScalarFieldType & volume2 =
       meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
    Part  & universal2     = meta2.universal_part ();
    put_field ( temperature2 , top.node_rank , universal2 );
    put_field ( volume2 , top.element_rank, universal2 );
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
              << endl ;
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

  const int spatial_dimension = 3;
  MetaData meta ( TopologicalMetaData::entity_rank_names ( spatial_dimension ) );
  TopologicalMetaData top( meta, spatial_dimension );

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  Part & partLeft_1 = top.declare_part<shards::Tetrahedron<4> >( "block_left_1" );

  Part & partLeft_2 = top.declare_part<shards::Tetrahedron<4> >( "block_left_2" );

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
    bulk.declare_entity( top.node_rank , new_id , no_part );
  }

  bulk.modification_end();

  const std::vector<Bucket*> & buckets = bulk.buckets( top.element_rank );

  std::vector<Bucket*>::const_iterator k;

  k = buckets.begin();

  //test 1 covers aecond section of "if" statement in while loop
  get_involved_parts( union_parts, **k, involved_parts);

  //test 2 covers union_parts.size() = 0
  PartVector union_parts2(0) ;
  get_involved_parts( union_parts2, **k, involved_parts);

  //test 3 covers first section of "if" statement in while loop
  const std::vector<Bucket*> & buckets2 = bulk.buckets( top.node_rank );
  std::vector<Bucket*>::const_iterator k2;

  k2 = buckets2.begin();
  get_involved_parts( union_parts, **k2, involved_parts);

  // tests on throw_error and BucketIterator in bucket.cpp/hpp

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }
  typedef Field<double>  ScalarFieldType;

  MetaData meta2 ( entity_names );
  BulkData bulk2( meta2 , pm , 4 );

  ScalarFieldType & temperature2 =
     meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
  ScalarFieldType & volume2 =
     meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
  Part  & universal     = meta2.universal_part ();
  put_field ( temperature2 , top.node_rank , universal );
  put_field ( volume2 , top.element_rank , universal );
  meta2.commit();

  bulk2.modification_begin();
  bulk2.declare_entity( top.edge_rank, rank+1 , no_part );
  bulk2.modification_end();

  const std::vector<Bucket*> & buckets3 = bulk2.buckets( top.edge_rank );

  std::vector<Bucket*>::const_iterator k3;

  k3 = buckets3.begin();

  Bucket& b3 = **k3;
  BucketIterator bitr3 = b3.begin();

  Bucket& b2 = **k2;
  BucketIterator bitr2 = b2.begin();

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
      cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << endl ;
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
      cout << "UnitTestBucket CORRECTLY caught error for : "
                << x.what()
                << endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestBucket FAILED to catch error for iterator from a different bucket");
    }
  }

}


void UnitTestBucket::testBucket2(ParallelMachine pm)
{

  // Tests to cover print, has_superset and BucketLess::operator() for Buckets.cpp - C.Brickley - 2nd June 2010

  const int spatial_dimension = 3;
  MetaData meta ( TopologicalMetaData::entity_rank_names ( spatial_dimension ) );
  TopologicalMetaData top( meta, spatial_dimension );

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  Part & partLeft_1 = meta.declare_part( "block_left_1", top.element_rank );

  Part & partLeft_3 = top.declare_part<shards::Tetrahedron<4> >( "block_left_3" );

  meta.commit();

  BulkData bulk( meta , pm , 100 );
  std::vector<Part *>  add_part4;
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

  const std::vector<Bucket*> & buckets2 = bulk.buckets( top.element_rank );

  std::vector<Bucket*>::const_iterator k2;

  k2 = buckets2.begin();

  Bucket& b2 = **k2;
  BucketIterator bitr2 = b2.begin();

  //define a new meta and bulkdata
  std::vector<std::string> entity_names(10);

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  typedef Field<double>  ScalarFieldType;

  MetaData meta2 ( entity_names );
  BulkData bulk2( meta2 , pm , 4 );

  ScalarFieldType & temperature2 =
       meta2.declare_field < ScalarFieldType > ( "temperature2" , 4 );
  ScalarFieldType & volume2 =
       meta2.declare_field < ScalarFieldType > ( "volume2" , 4 );
  Part  & universal     = meta2.universal_part ();
  put_field ( temperature2 , top.node_rank , universal );
  put_field ( volume2 , top.element_rank , universal );

  typedef Field<double>  VectorFieldType;
  typedef Field<double>  ElementNodePointerFieldType;

  meta2.commit();

  //Test to cover print function in Bucket.cpp
  cout << endl << "Bucket test" << endl ;
  print(cout, "  ", b2);

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

void UnitTestBucket::test_EntityComm( ParallelMachine pm )
{
  const int spatial_dimension = 3;
  MetaData meta ( TopologicalMetaData::entity_rank_names ( spatial_dimension ) );
  TopologicalMetaData top( meta, spatial_dimension );

  BulkData bulk ( meta , pm , 100 );
  std::vector<Part *>  add_part4;

  cout << endl << "Bucket test line 0.1" << endl ;
  Part & partLeft_1 = top.declare_part<shards::Tetrahedron<4> >( "block_left_1" );
  cout << endl << "Bucket test line 0.2" << endl;
  meta.commit();

  add_part4.push_back ( &partLeft_1 );

  int  size , rank;
  rank = stk::parallel_machine_rank( pm );
  size = stk::parallel_machine_size( pm );
  PartVector tmp(1);

  bulk.modification_begin();
  cout << endl << "Bucket test line 1" << endl ;
  //int id_base = 0;
  //int new_id = size * id_base + rank;
  //  for ( id_base = 0 ; id_base < 93 ; ++id_base )
  //  {
  //   int new_id = size * id_base + rank;
  //   bulk.declare_entity( 0 , new_id+1 , add_part4 );
  //  }

  bulk.modification_end();

  /*  cout << endl << "Bucket test line 3" << endl ;
  bool result = in_shared(elem);
  if( result) {
     STKUNIT_ASSERT_EQUAL( result , true );
  }
  cout << endl << "Bucket test line 4" << endl ;

  result = in_receive_ghost(elem);
  if( result) {
     STKUNIT_ASSERT_EQUAL( result , true );
  }

    for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
      cout << endl << "in relation h and p =" << p << endl ;

      STKUNIT_ASSERT_EQUAL( in_send_ghost( *elem , p ), false );
      cout << endl << "in relation ii =" << endl
   }

  cout << endl << "Bucket test line 5" << endl ;
  result = in_send_ghost(elem);
  if( result) {
     STKUNIT_ASSERT_EQUAL( result , true );
     }

  cout << endl << "Bucket test line 6" << endl ;

  unsigned proc = rank;
  unsigned procnew = rank+10;

  result = in_shared(elem, proc);
  if( result) {
     STKUNIT_ASSERT_EQUAL( result , true );
  }
  cout << endl << "Bucket test line 7" << endl ;  */
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------


