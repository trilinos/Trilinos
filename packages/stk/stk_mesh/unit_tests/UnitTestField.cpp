
#include <unit_tests/stk_utest_macros.hpp>

#include <stdexcept>
#include <iostream>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>

namespace {
void testCartesian();
void testFieldDataArray( stk::ParallelMachine pm );

STKUNIT_UNIT_TEST(UnitTestField, testUnit)
{
#if defined( STK_HAS_MPI )
  stk::ParallelMachine pworld = MPI_COMM_WORLD ;
  stk::ParallelMachine pself  = MPI_COMM_SELF ;
#else
  stk::ParallelMachine pworld = parallel_machine_null();
  stk::ParallelMachine pself  = parallel_machine_null();
#endif
  if ( 0 == stk::parallel_machine_rank( pworld ) ) {
    // Nothing parallel being tested, only run on one process
    testCartesian();
    testFieldDataArray( pself );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( CTAG )

template< class FieldType >
void print_bucket_array( const FieldType & f , const stk::mesh::Bucket & k )
{
  typedef stk::mesh::BucketArray< FieldType > ArrayType ;

  try {
    ArrayType a( f , k.begin(), k.end() );
    ArrayType b( f , k );

    std::cout << "  BucketArray[" << f.name() << "](" ;

    if ( a.size() != b.size() ) {
      throw std::runtime_error("UnitTestField FAILED BucketArray dimensions not consistant with Bucket::iterator");
    }

    if ( a.size() ) {
      for ( unsigned i = 0 ; i < ArrayType::Rank ; ++i ) {
        if ( i ) { std::cout << "," ; }
        std::cout << a.dimension(i);
        if (a.dimension(i) != b.dimension(i)) {
          throw std::runtime_error("UnitTestField FAILED BucketArray dimensions not consistant with Bucket::iterator");
        }
      }
    }
    std::cout << ")" << std::endl ;
  }
  catch( const std::exception & ) {
  }
}

void testCartesian()
{
  const stk::mesh::Cartesian&  cartesian_tag = stk::mesh::Cartesian::tag();

  std::string to_str = cartesian_tag.to_string(3, 1);
  std::string expected_str("y");
  STKUNIT_ASSERT_EQUAL( to_str, expected_str );

  //should throw if we supply a size < 3:
  STKUNIT_ASSERT_THROW( cartesian_tag.to_string(2, 1), std::runtime_error );

  shards::ArrayDimTag::size_type expected_idx = 1;
  shards::ArrayDimTag::size_type idx = cartesian_tag.to_index(3, "y");

  STKUNIT_ASSERT_EQUAL( idx, expected_idx );

  //should throw if we supply a "z" along with size==2:
  STKUNIT_ASSERT_THROW( cartesian_tag.to_index(2, "z"), std::runtime_error );
}

void testFieldDataArray( stk::ParallelMachine pm )
{
  typedef stk::mesh::Field<double>                rank_zero_field ;
  typedef stk::mesh::Field<double,ATAG>           rank_one_field ;
  typedef stk::mesh::Field<double,ATAG,BTAG>      rank_two_field ;
  typedef stk::mesh::Field<double,ATAG,BTAG,CTAG> rank_three_field ;

  std::cout << "UnitTestField BEGIN:" << std::endl ;

  const std::string name0("test_field_0");
  const std::string name1("test_field_1");
  const std::string name2("test_field_2");
  const std::string name3("test_field_3");

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );

  rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field >( name0 );
  rank_one_field   & f1 = meta_data.declare_field< rank_one_field >(  name1 );
  rank_three_field & f3 = meta_data.declare_field< rank_three_field >( name3 );
  rank_two_field   & f2 = meta_data.declare_field< rank_two_field >(  name2 );

  {
    int ok = 0 ;
    try {
      typedef stk::mesh::Field<double,CTAG> error_type ;
      meta_data.declare_field< error_type >( name1 );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestField CORRECTLY caught error: "
                << x.what()
                << std::endl ;
    }
    if ( ! ok ) {
      throw std::runtime_error("UnitTestField FAILED to catch error");
    }
  }

  stk::mesh::Part & p0 = meta_data.declare_part("P0", stk::mesh::Node );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", stk::mesh::Node );
  stk::mesh::Part & p2 = meta_data.declare_part("P2", stk::mesh::Node );
  stk::mesh::Part & p3 = meta_data.declare_part("P3", stk::mesh::Node );

  stk::mesh::put_field( f0 , stk::mesh::Node , p0 );
  stk::mesh::put_field( f1 , stk::mesh::Node , p1 , 10 );
  stk::mesh::put_field( f2 , stk::mesh::Node , p2 , 10 , 20 );
  stk::mesh::put_field( f3 , stk::mesh::Node , p3 , 10 , 20 , 30 );

  stk::mesh::print( std::cout , "  " , f0 ); std::cout << std::endl ;

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data , pm );

  for ( unsigned i = 1 ; i < 11 ; ++i ) {
    bulk_data.declare_entity( stk::mesh::Node , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p0 ) );
  }

  for ( unsigned i = 11 ; i < 21 ; ++i ) {
    bulk_data.declare_entity( stk::mesh::Node , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p1 ) );
  }

  for ( unsigned i = 21 ; i < 31 ; ++i ) {
    bulk_data.declare_entity( stk::mesh::Node , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p2 ) );
  }

  for ( unsigned i = 31 ; i < 41 ; ++i ) {
    bulk_data.declare_entity( stk::mesh::Node , i ,
                              std::vector< stk::mesh::Part * >( 1 , & p3 ) );
  }

  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( stk::mesh::Node );

  for ( std::vector< stk::mesh::Bucket *>::const_iterator
        ik = node_buckets.begin() ; ik != node_buckets.end() ; ++ik ) {
    stk::mesh::Bucket & k = **ik ;

    std::vector< stk::mesh::Part * > parts ;
    k.supersets( parts );

    std::cout << "Bucket:" ;
    for ( std::vector< stk::mesh::Part * >::iterator
          ip = parts.begin() ; ip != parts.end() ; ++ip ) {
      std::cout << " " << (*ip)->name();
    }
    std::cout << std::endl ;

    print_bucket_array( f0 , k );
    print_bucket_array( f1 , k );
    print_bucket_array( f2 , k );
    print_bucket_array( f3 , k );
  }

  std::cout << "UnitTestField END" << std::endl ;
}

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( CTAG )

}//namespace <anonymous>

