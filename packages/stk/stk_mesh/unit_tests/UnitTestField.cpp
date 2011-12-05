/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <boost/range.hpp>
#include <boost/foreach.hpp>

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

typedef shards::ArrayDimTag::size_type size_type;

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( CTAG )

template< class FieldType >
void print_bucket_array( const FieldType & f , const stk::mesh::Bucket & k )
{
  typedef stk::mesh::BucketArray< FieldType > ArrayType ;
  std::ostringstream oss;

  ArrayType a( f , k.begin(), k.end() );
  ArrayType b( f , k );

  oss << "  BucketArray[" << f.name() << "](" ;

  if ( a.size() != b.size() ) {
    throw std::runtime_error("UnitTestField FAILED BucketArray dimensions not consistant with Bucket::iterator");
  }

  if ( a.size() ) {
    for ( unsigned i = 0 ; i < ArrayType::Rank ; ++i ) {
      if ( i ) { oss << "," ; }
      oss << a.dimension(i);
      if (a.dimension(i) != b.dimension(i)) {
        throw std::runtime_error("UnitTestField FAILED BucketArray dimensions not consistant with Bucket::iterator");
      }
    }
  }
  oss << ")" << std::endl ;
}

STKUNIT_UNIT_TEST(UnitTestField, testCartesian)
{
  // Test the Cartesian array dimension tag

  const stk::mesh::Cartesian& cartesian_tag = stk::mesh::Cartesian::tag();

  std::string to_str = cartesian_tag.to_string(3 /*size*/, 1 /*idx*/);
  std::string expected_str("y");
  STKUNIT_ASSERT_EQUAL( (to_str == expected_str), true);

  //should throw if we supply a size < 3:
  STKUNIT_ASSERT_THROW( cartesian_tag.to_string(2 /*size*/, 1 /*idx*/),
                        std::runtime_error );

  size_type expected_idx = 1;
  size_type idx = cartesian_tag.to_index(3 /*size*/, "y" /*dim*/);
  STKUNIT_ASSERT_EQUAL( idx, expected_idx );

  //should throw if we supply a "z" along with size==2:
  STKUNIT_ASSERT_THROW( cartesian_tag.to_index(2 /*size*/, "z" /*dim*/),
                        std::runtime_error );
}

STKUNIT_UNIT_TEST(UnitTestField, testCylindrical)
{
  // Test the Cylindrical array dimension tag

  const stk::mesh::Cylindrical& cylindrical_tag =stk::mesh::Cylindrical::tag();

  std::string to_str = cylindrical_tag.to_string(3 /*size*/, 1 /*idx*/);
  std::string expected_str("a");
  STKUNIT_ASSERT_EQUAL( (to_str == expected_str), true );

  //should throw if we supply a size < 3:
  STKUNIT_ASSERT_THROW( cylindrical_tag.to_string(2 /*size*/, 1 /*idx*/),
                        std::runtime_error );

  size_type expected_idx = 1;
  size_type idx = cylindrical_tag.to_index(3 /*size*/, "a" /*dim*/);
  STKUNIT_ASSERT_EQUAL( idx, expected_idx );

  //should throw if we supply a "z" along with size==2:
  STKUNIT_ASSERT_THROW( cylindrical_tag.to_index(2 /*size*/, "z" /*dim*/),
                        std::runtime_error );
}

STKUNIT_UNIT_TEST(UnitTestField, testFullTensor)
{
  // Test the FullTensor array dimension tag

  const stk::mesh::FullTensor&  fulltensor_tag = stk::mesh::FullTensor::tag();

  std::string to_str = fulltensor_tag.to_string(9 /*size*/, 1 /*idx*/);
  std::string expected_str("yy");
  STKUNIT_ASSERT_EQUAL( (to_str == expected_str), true );

  //should throw if we supply a size < 9:
  STKUNIT_ASSERT_THROW( fulltensor_tag.to_string(2 /*size*/, 1 /*idx*/),
                        std::runtime_error );

  size_type expected_idx = 6;
  size_type idx = fulltensor_tag.to_index(9 /*size*/, "yx" /*dim*/);
  STKUNIT_ASSERT_EQUAL( idx, expected_idx );

  //should throw if we supply a "zz" along with size==2:
  STKUNIT_ASSERT_THROW( fulltensor_tag.to_index(2 /*size*/, "zz" /*dim*/),
                        std::runtime_error );
}

STKUNIT_UNIT_TEST(UnitTestField, testSymmetricTensor)
{
  // Test the SymmetricTensor array dimension tag

  const stk::mesh::SymmetricTensor& symmetrictensor_tag =
    stk::mesh::SymmetricTensor::tag();

  std::string to_str = symmetrictensor_tag.to_string(9 /*size*/, 1 /*idx*/);
  std::string expected_str("yy");
  STKUNIT_ASSERT_EQUAL( (to_str == expected_str), true);

  //should throw if we supply a size < 9:
  STKUNIT_ASSERT_THROW( symmetrictensor_tag.to_string(2 /*size*/, 1 /*idx*/),
                        std::runtime_error );

  size_type expected_idx = 1;
  size_type idx = symmetrictensor_tag.to_index(6 /*size*/, "yy" /*dim*/);
  STKUNIT_ASSERT_EQUAL( idx, expected_idx );

  //should throw if we supply a "xz" along with size==5:
  STKUNIT_ASSERT_THROW( symmetrictensor_tag.to_index(5 /*size*/, "xz" /*dim*/),
                        std::runtime_error );
}

STKUNIT_UNIT_TEST(UnitTestField, testFieldDataArray)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for some test fields
  typedef stk::mesh::Field<double>                rank_zero_field ;
  typedef stk::mesh::Field<double,ATAG>           rank_one_field ;
  typedef stk::mesh::Field<double,ATAG,BTAG>      rank_two_field ;
  typedef stk::mesh::Field<double,ATAG,BTAG,CTAG> rank_three_field ;

  const std::string name0("test_field_0");
  const std::string name1("test_field_1");
  const std::string name2("test_field_2");
  const std::string name3("test_field_3");

  const int spatial_dimension = 3;
  stk::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( stk::mesh::fem::FEMMetaData::get_meta_data(meta_data) , pm );

  rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field >( name0 );
  rank_one_field   & f1 = meta_data.declare_field< rank_one_field >(  name1 );
  rank_three_field & f3 = meta_data.declare_field< rank_three_field >( name3 );
  rank_two_field   & f2 = meta_data.declare_field< rank_two_field >(  name2 );

  // confirm that declaring field with erroneous type throws exception
  typedef stk::mesh::Field<double,CTAG> error_type ;
  STKUNIT_ASSERT_THROW(meta_data.declare_field< error_type >( name1 ),
                       std::runtime_error);

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK );
  stk::mesh::Part & p2 = meta_data.declare_part("P2", NODE_RANK );
  stk::mesh::Part & p3 = meta_data.declare_part("P3", NODE_RANK );

  stk::mesh::put_field( f0 , NODE_RANK , p0 );
  stk::mesh::put_field( f1 , NODE_RANK , p1 , 10 );
  stk::mesh::put_field( f2 , NODE_RANK , p2 , 10 , 20 );
  stk::mesh::put_field( f3 , NODE_RANK , p3 , 10 , 20 , 30 );

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare a 10 nodes on each part

  for ( unsigned i = 1 ; i < 11 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p0 ) );
  }

  for ( unsigned i = 11 ; i < 21 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p1 ) );
  }

  for ( unsigned i = 21 ; i < 31 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p2 ) );
  }

  for ( unsigned i = 31 ; i < 41 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p3 ) );
  }

  // Go through node_buckets and print the all the fields on each bucket
  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( NODE_RANK );

  for ( std::vector< stk::mesh::Bucket *>::const_iterator
        ik = node_buckets.begin() ; ik != node_buckets.end() ; ++ik ) {
    stk::mesh::Bucket & k = **ik ;

    std::vector< stk::mesh::Part * > parts ;
    k.supersets( parts );

    for ( std::vector< stk::mesh::Part * >::iterator
          ip = parts.begin() ; ip != parts.end() ; ++ip ) {
      oss << " " << (*ip)->name();
    }

    print_bucket_array( f0 , k );
    print_bucket_array( f1 , k );
    print_bucket_array( f2 , k );
    print_bucket_array( f3 , k );
  }
}

STKUNIT_UNIT_TEST(UnitTestField, testFieldWithSelector)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for test field
  typedef stk::mesh::Field<double>    rank_zero_field ;

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( stk::mesh::fem::FEMMetaData::get_meta_data(meta_data) , pm );

  rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field >( name0 );

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK );

  stk::mesh::Selector select_p0 = p0;
  std::cout <<"select_p0: "<< select_p0 << std::endl;

  stk::mesh::put_field( f0 , NODE_RANK , select_p0 );

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare 10 nodes on each part

  for ( unsigned i = 1 ; i < 11 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p0 ) );
  }

  for ( unsigned i = 11 ; i < 21 ; ++i ) {
    bulk_data.declare_entity( NODE_RANK , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p1 ) );
  }

  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( NODE_RANK );

  unsigned num = stk::mesh::count_selected_entities(select_p0, node_buckets);
  
  STKUNIT_ASSERT_EQUAL( 10u, num );

  stk::mesh::Selector select_f0 = stk::mesh::selectField(f0);

  std::cout <<"select_f0: "<< select_f0 << std::endl;

  unsigned num_f0 = stk::mesh::count_selected_entities(select_f0, node_buckets);
  STKUNIT_ASSERT_EQUAL(10u, num_f0);

  std::vector<stk::mesh::Bucket*> f0_buckets;
  stk::mesh::get_buckets(select_p0, bulk_data.buckets(NODE_RANK), f0_buckets);
  unsigned num_buckets = f0_buckets.size();
  STKUNIT_ASSERT_EQUAL(1u, num_buckets);

  BOOST_FOREACH(stk::mesh::Bucket* b, f0_buckets) {
    unsigned f0_size = b->field_data_size(f0);
    STKUNIT_ASSERT_EQUAL(8u, f0_size);
  }
}

STKUNIT_UNIT_TEST(UnitTestField, testFieldWithSelectorAnd)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  typedef stk::mesh::Field<double,shards::ArrayDimension>           rank_one_field ;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( stk::mesh::fem::FEMMetaData::get_meta_data(meta_data) , pm );

  rank_one_field  & f0 = meta_data.declare_field< rank_one_field >( name0 );

  stk::mesh::EntityRank elem_rank = meta_data.element_rank();
  stk::mesh::Part & elements = meta_data.declare_part("Elements", elem_rank);
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );
  stk::mesh::Part & tet4s = meta_data.declare_part("Tet4", elem_rank );

  stk::mesh::Selector elem_hex_selector = elements & hex8s;
  stk::mesh::Selector elem_tet_selector = elements & tet4s;
  std::cout <<"elem_hex_selector: "<< elem_hex_selector << std::endl;
  std::cout <<"elem_tet_selector: "<< elem_tet_selector << std::endl;

  stk::mesh::put_field( f0 , elem_rank , elem_hex_selector, 8u );
  stk::mesh::put_field( f0 , elem_rank , elem_tet_selector, 4u );

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare 10 elements on each part

  stk::mesh::PartVector parts;
  parts.push_back(&elements);
  parts.push_back(&hex8s);

  for ( unsigned i = 1 ; i < 11 ; ++i ) {
    bulk_data.declare_entity( elem_rank , i , parts );
  }

  parts.clear();
  parts.push_back(&elements);
  parts.push_back(&tet4s);

  for ( unsigned i = 11 ; i < 21 ; ++i ) {
    bulk_data.declare_entity( elem_rank , i , parts );
  }

  stk::mesh::BucketVector f0_buckets;
  stk::mesh::get_buckets(elem_hex_selector, bulk_data.buckets(elem_rank), f0_buckets);

  BOOST_FOREACH(stk::mesh::Bucket* b, f0_buckets) {
    unsigned f0_size = b->field_data_size(f0);
    STKUNIT_ASSERT_EQUAL(64u, f0_size);
  }

  f0_buckets.clear();

  stk::mesh::get_buckets(elem_tet_selector, bulk_data.buckets(elem_rank), f0_buckets);

  BOOST_FOREACH(stk::mesh::Bucket* b, f0_buckets) {
    unsigned f0_size = b->field_data_size(f0);
    STKUNIT_ASSERT_EQUAL(32u, f0_size);
  }
}


STKUNIT_UNIT_TEST(UnitTestField, testFieldWithSelectorInvalid)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  typedef stk::mesh::Field<double,shards::ArrayDimension>           rank_one_field ;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( stk::mesh::fem::FEMMetaData::get_meta_data(meta_data) , pm );

  rank_one_field  & f0 = meta_data.declare_field< rank_one_field >( name0 );

  stk::mesh::EntityRank elem_rank = meta_data.element_rank();
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );

  stk::mesh::Part & universal_part = meta_data.universal_part();
  stk::mesh::Selector elem_hexA_selector = hex8s;
  stk::mesh::Selector elem_hexB_selector = universal_part & hex8s;

  std::cout <<"elem_hexA_selector: "<< elem_hexA_selector << std::endl;
  std::cout <<"elem_hexB_selector: "<< elem_hexB_selector << std::endl;

  stk::mesh::put_field( f0 , elem_rank , elem_hexA_selector, 8u );
  STKUNIT_ASSERT_THROW( 
    stk::mesh::put_field( f0 , elem_rank , elem_hexA_selector, 4u ),
    std::runtime_error
  );
  stk::mesh::put_field( f0 , elem_rank , elem_hexB_selector, 4u );

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  stk::mesh::PartVector parts;
  parts.push_back(&hex8s);
  STKUNIT_ASSERT_THROW( 
    bulk_data.declare_entity( elem_rank , 1 , parts ),
    std::runtime_error
  );

}

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( CTAG )

} //namespace <anonymous>

