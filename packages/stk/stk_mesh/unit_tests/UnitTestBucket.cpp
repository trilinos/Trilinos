/*------------------------------------------------------------------------*/
/*         _        Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t
#include <sstream>                      // for ostringstream, etc
#include <stk_mesh/base/Bucket.hpp>     // for has_superset, Bucket, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetBuckets.hpp>  // for get_involved_parts
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc
#include "mpi.h"                        // for MPI_Barrier, MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class Part; } }

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
using stk::mesh::BucketVector;

typedef Field<double> ScalarFieldType;

namespace {

const EntityRank NODE_RANK = stk::topology::NODE_RANK;

TEST(UnitTestingOfBucket, testBucket)
{
  // Unit test the Part functionality in isolation:

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create a mesh for testing buckets...

  // Create dummy names for entity ranks to be given to MetaData
  std::vector<std::string> entity_names(5);
  for ( size_t i = 0 ; i < 5 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  // Create MetaData, BulkData
  unsigned max_bucket_size = 4;
  stk::mesh::fixtures::BoxFixture fixture(pm, max_bucket_size, entity_names);
  MetaData& meta = fixture.fem_meta();
  BulkData& bulk = fixture.bulk_data();
  // Create two scalar fields, temperature and volume. Put temperature
  // on all the nodes and put volume on all the elements.
  unsigned number_of_states = 4;

  ScalarFieldType & temperature =
    meta.declare_field < ScalarFieldType > (stk::topology::NODE_RANK, "temperature" , number_of_states );
  ScalarFieldType & volume =

    meta.declare_field < ScalarFieldType > ( stk::topology::ELEMENT_RANK, "volume" , number_of_states );
  Part & universal     = meta.universal_part ();
  put_field ( temperature , universal );
  put_field ( volume , universal );
  meta.commit();

  // Generate the mesh
  int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  //  First, test for streaming IO;
  {
    std::string gold1;
    gold1 = "Bucket( EntityRank0 : {UNIVERSAL} {OWNS} )";
    Bucket *b1 = bulk.buckets(stk::topology::NODE_RANK)[0];
    std::stringstream  out1_str;
    out1_str << (*b1);
    bool equal = (gold1 == out1_str.str());
    ASSERT_EQ ( equal, true );
  }

  // Second, update state of bucket until circular cue is filled
  {
    /* Need to set some data in state, rotate look for it, rotate 3 more times
       and look for it again */
    for ( size_t i = 0 ; i != 10 ; ++i )
      bulk.update_field_data_states ();
  }

  // next, check has_superset (...) and membership functions
  {
    PartVector tmp(2) ;
    tmp[0] = & meta.universal_part();
    tmp[1] = & meta.locally_owned_part();
    ASSERT_EQ ( has_superset ( *bulk.buckets(stk::topology::NODE_RANK)[0] , tmp ) , bulk.parallel_size() == 1 );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_any ( tmp ) );
    ASSERT_EQ ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_all ( tmp ) , bulk.parallel_size() == 1 );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member ( **meta.get_parts().begin() ) );
  }
}

TEST(UnitTestingOfBucket, testGetInvolvedParts)
{
  // Tests to cover get_involved_parts for GetBuckets.cpp - C.Brickley - 12 May 2010

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const int spatial_dimension = 3;

  MetaData meta( spatial_dimension );
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  const EntityRank edge_rank    = stk::topology::EDGE_RANK;

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  Part & partLeft_1 = meta.declare_part_with_topology("block_left_1", stk::topology::TET_4 );

  Part & partLeft_2 = meta.declare_part_with_topology("block_left_2", stk::topology::TET_4 );

  meta.commit();

  PartVector union_parts;
  union_parts.push_back(&partLeft_1);
  union_parts.push_back(&partLeft_2);

  BulkData bulk( meta , pm , 100 );
  PartVector add_part4, no_part;
  add_part4.push_back ( &partLeft_1 );

  bulk.modification_begin();
  int rank = stk::parallel_machine_rank( pm );
  int size = stk::parallel_machine_size( pm );

  for ( int id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank + 1;
    bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id , add_part4 );
    bulk.declare_entity( NODE_RANK , new_id , no_part );
  }

  bulk.modification_end();

  const BucketVector & buckets = bulk.buckets( element_rank );

  BucketVector::const_iterator k;

  k = buckets.begin();

  //test 1 covers aecond section of "if" statement in while loop
  get_involved_parts( union_parts, **k, involved_parts);

  //test 2 covers union_parts.size() = 0
  PartVector union_parts2(0) ;
  get_involved_parts( union_parts2, **k, involved_parts);

  //test 3 covers first section of "if" statement in while loop
  const BucketVector & buckets2 = bulk.buckets( NODE_RANK );
  BucketVector::const_iterator k2;

  k2 = buckets2.begin();
  get_involved_parts( union_parts, **k2, involved_parts);

  // tests on throw_error and BucketIterator in bucket.cpp/hpp

  MetaData meta2 (spatial_dimension);
  BulkData bulk2( meta2 , pm , 4 );

  unsigned number_of_states = 4;

  ScalarFieldType & temperature2 =
    meta2.declare_field < ScalarFieldType >(stk::topology::NODE_RANK, "temperature2" , number_of_states);
  ScalarFieldType & volume2 =

    meta2.declare_field < ScalarFieldType >(stk::topology::ELEMENT_RANK, "volume2", number_of_states);
  Part  & universal = meta2.universal_part ();
  put_field ( temperature2 , universal );
  put_field ( volume2 , universal );

  meta2.commit();

  bulk2.modification_begin();
  bulk2.declare_entity( edge_rank, rank+1 , no_part );
  bulk2.modification_end();
}

TEST(UnitTestingOfBucket, testBucket2)
{
  // Tests to cover print, has_superset and BucketLess::operator() for Buckets.cpp - C.Brickley - 2nd June 2010

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;

  PartVector involved_parts(2) ;
  involved_parts[0] = & meta.universal_part();
  involved_parts[1] = & meta.locally_owned_part();

  Part & partLeft_1 = meta.declare_part("block_left_1", element_rank);

  Part & partLeft_3 = meta.declare_part_with_topology("block_left_3", stk::topology::TET_4 );

  meta.commit();

  BulkData bulk( meta , pm , 100 );
  std::vector<Part *>  add_part4;
  add_part4.push_back ( &partLeft_1 );

  bulk.modification_begin();
  int rank = stk::parallel_machine_rank( pm );
  int size = stk::parallel_machine_size( pm );

  for ( int id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id+1 , add_part4 );
  }

  bulk.modification_end();

  const BucketVector & buckets2 = bulk.buckets( element_rank );

  BucketVector::const_iterator k2;

  k2 = buckets2.begin();

  Bucket& b2 = **k2;

  //define a new meta and bulkdata
  std::vector<std::string> entity_names(10);

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta2 ( spatial_dimension, entity_names );
  BulkData bulk2( meta2 , pm , 4 );

  unsigned number_of_states = 4;

  ScalarFieldType & temperature2 =
    meta2.declare_field < ScalarFieldType >(stk::topology::NODE_RANK, "temperature2" , number_of_states);
  ScalarFieldType & volume2 =
    meta2.declare_field < ScalarFieldType >(stk::topology::ELEMENT_RANK, "volume2", number_of_states);
  Part  & universal     = meta2.universal_part ();
  put_field ( temperature2 , universal );
  put_field ( volume2 , universal );

  typedef Field<double>  VectorFieldType;

  meta2.commit();

  //Test to cover print function in Bucket.cpp
  std::ostringstream oss;
  print(oss, "  ", b2);

  //Test to cover has_superset function in Bucket.cpp
  ASSERT_EQ ( has_superset ( b2 , partLeft_3 ) , false );

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

  // FIXME: The code below needs to be fixed or removed
  /*
  Bucket::last_bucket_in_family( *k2 );

  const unsigned * t = key;
  const Bucket * u = last_bucket;

  BucketLess Buck;

  bool res = Buck(  &t[0], &u[0] );

  EXPECT_EQ( res, false );
  */
}

TEST(UnitTestingOfBucket, testEntityComm)
{
  // FIXME: With so much code commented out, this unit-test does
  // not appear to be testing anything.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );

  BulkData bulk ( meta , pm , 100 );
  std::vector<Part *>  add_part4;

  Part & partLeft_1 = meta.declare_part_with_topology( "block_left_1", stk::topology::TET_4 );
  meta.commit();

  add_part4.push_back ( &partLeft_1 );

  //int rank = stk::parallel_machine_rank( pm );
  // int size = stk::parallel_machine_size( pm );
  PartVector tmp(1);

  bulk.modification_begin();

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
     ASSERT_EQ( result , true );
  }
  cout << endl << "Bucket test line 4" << endl ;

  result = in_receive_ghost(elem);
  if( result) {
     ASSERT_EQ( result , true );
  }

    for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
      cout << endl << "in relation h and p =" << p << endl ;

      ASSERT_EQ( in_send_ghost( *elem , p ), false );
      cout << endl << "in relation ii =" << endl
   }

  cout << endl << "Bucket test line 5" << endl ;
  result = in_send_ghost(elem);
  if( result) {
     ASSERT_EQ( result , true );
     }

  cout << endl << "Bucket test line 6" << endl ;

  unsigned proc = rank;
  unsigned procnew = rank+10;

  result = in_shared(elem, proc);
  if( result) {
     ASSERT_EQ( result , true );
  }
  cout << endl << "Bucket test line 7" << endl ;  */
}

}
