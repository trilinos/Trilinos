/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

namespace stk {
namespace mesh {

class UnitTestMetaData {
public:
  UnitTestMetaData() {}

  void testEntityKey();
  void testPart();
  void testPartVector();
  void testField();
  void testFieldRestriction();
  void testMetaData();
  void testProperty();
  void testEntityRepository();
};

}//namespace mesh
}//namespace stk

namespace {

STKUNIT_UNIT_TEST(UnitTestEntityKey, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testEntityKey();
}

STKUNIT_UNIT_TEST(UnitTestPart, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testPart();
}

STKUNIT_UNIT_TEST(UnitTestPartVector, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testPartVector();
}

STKUNIT_UNIT_TEST(UnitTestField, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testField();
}

STKUNIT_UNIT_TEST(UnitTestFieldRestriction, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testFieldRestriction();
}

STKUNIT_UNIT_TEST(UnitTestMetaData, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testMetaData();
}

STKUNIT_UNIT_TEST(UnitTestProperty, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testProperty();
}

STKUNIT_UNIT_TEST(UnitTestEntityRepository, testUnit)
{
  stk::mesh::UnitTestMetaData umeta;
  umeta.testEntityRepository();
}

}//namespace <anonymous>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
void UnitTestMetaData::testEntityKey()
{
  EntityKey
    key_bad_zero  = EntityKey(),
    key_good_0_1  = EntityKey( 0 , 1 ),
    key_good_1_1  = EntityKey( 1 , 1 ),
    key_good_2_10 = EntityKey( 2 , 10 );

//   EntityId val = entity_id(key_good_2_10) ;

//   EntityKey good_key_gone_bad;
//   entity_good_key_gone_bad.key = val ;

  STKUNIT_ASSERT( ! entity_key_valid( key_bad_zero ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_0_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_1_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_2_10 ) );

  STKUNIT_ASSERT( 0  == entity_rank(key_good_0_1));
  STKUNIT_ASSERT( 1  == entity_rank(key_good_1_1) );
  STKUNIT_ASSERT( 2  == entity_rank(key_good_2_10) );
  STKUNIT_ASSERT( 1  == entity_id(key_good_0_1) );
  STKUNIT_ASSERT( 1  == entity_id(key_good_1_1) );
  STKUNIT_ASSERT( 10 == entity_id(key_good_2_10) );

  EntityKey
    key_order_1_12 = EntityKey( 1 , 12 ),
    key_order_2_10 = EntityKey( 2 , 10 );

  STKUNIT_ASSERT( key_order_1_12 < key_order_2_10);
  STKUNIT_ASSERT( !(key_order_1_12 > key_order_2_10));

//  STKUNIT_ASSERT( ! entity_key_valid( good_key_gone_bad ) );

//   std::cout.unsetf( std::ios::dec);
//   std::cout.setf( std::ios::hex);
//   std::cout << "TEST entity_key_type "
//             << ", key_good_2_10 = " << key_good_2_10.key
//             << ", good_key_gone_bad = " << good_key_gone_bad.key
//             << std::endl ;
//   std::cout.unsetf( std::ios::hex);
//   std::cout.setf( std::ios::dec);

  EntityKey key01(0,1), key_default;
  key01 = key_default;

  STKUNIT_ASSERT_THROW( EntityKey( ~0u , 1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( EntityKey( 0 , ~stk::mesh::EntityKey::raw_key_type(0) ) , std::runtime_error );
}

//----------------------------------------------------------------------
// Unit test the Part functionality in isolation:

void UnitTestMetaData::testPart()
{
  MetaData * const m = NULL ;

  impl::PartRepository partRepo(m);

  Part & universal = *partRepo.universal_part();

  STKUNIT_ASSERT( universal.supersets().empty() );
  STKUNIT_ASSERT( 1u == universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL( universal.subsets()[0] , & universal );

  //--------------------------------------------------------------------
  // Test multiple part creation

  enum { NPARTS = 100 };

  Part * parts[ NPARTS ] ;

  parts[0] = & universal ;

  for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
    std::ostringstream name ;
    name << "Part_" << i ;
    parts[i] = partRepo.declare_part( name.str() , 0 );
  }
  parts[99] = partRepo.declare_part( "Part_99" , 1 );

  STKUNIT_ASSERT( universal.supersets().empty() );
  STKUNIT_ASSERT( NPARTS == universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL( universal.subsets()[0] , & universal );

  for ( unsigned i = 1 ; i < NPARTS ; ++i ) {
    STKUNIT_ASSERT( parts[i]->subsets().empty() );
    STKUNIT_ASSERT( parts[i]->intersection_of().empty() );
    STKUNIT_ASSERT( parts[i]->mesh_meta_data_ordinal() == i );
    STKUNIT_ASSERT( 1u == parts[i]->supersets().size() );
    STKUNIT_ASSERT( & universal == parts[i]->supersets()[0] );
    STKUNIT_ASSERT_EQUAL( parts[i] , universal.subsets()[i] );
    STKUNIT_ASSERT_EQUAL( parts[i] , find( universal.subsets() , parts[i]->name() ) );
  }

  //--------------------------------------------------------------------
  // Test multiple parts and transitive subset declarations:

  partRepo.declare_subset( * parts[3], * parts[4] );
  partRepo.declare_subset( * parts[4], * parts[5] );

  partRepo.declare_subset( * parts[1], * parts[2] );
  // 1 and 2 pick up 4 and 5 via transitive relationship:
  partRepo.declare_subset( * parts[2], * parts[3] );

  STKUNIT_ASSERT( 4u == parts[1]->subsets().size() );
  STKUNIT_ASSERT( 3u == parts[2]->subsets().size() );
  STKUNIT_ASSERT( 2u == parts[3]->subsets().size() );
  STKUNIT_ASSERT( 1u == parts[4]->subsets().size() );
  STKUNIT_ASSERT( 0u == parts[5]->subsets().size() );

  STKUNIT_ASSERT( contain( parts[1]->subsets() , * parts[2] ) );
  STKUNIT_ASSERT( contain( parts[1]->subsets() , * parts[3] ) );
  STKUNIT_ASSERT( contain( parts[1]->subsets() , * parts[4] ) );
  STKUNIT_ASSERT( contain( parts[1]->subsets() , * parts[5] ) );

  STKUNIT_ASSERT( contain( parts[5]->supersets() , * parts[1] ) );
  STKUNIT_ASSERT( contain( parts[5]->supersets() , * parts[2] ) );
  STKUNIT_ASSERT( contain( parts[5]->supersets() , * parts[3] ) );
  STKUNIT_ASSERT( contain( parts[5]->supersets() , * parts[4] ) );

  //--------------------------------------------------------------------
  // Test declaration of an intersection part

  PartVector intersection ;
  intersection.push_back( parts[1] );
  intersection.push_back( parts[2] );
  intersection.push_back( parts[3] );
  intersection.push_back( parts[4] ); // Smallest subset of 1..4

  // Test filtering of trivial intersection

  STKUNIT_ASSERT( parts[4] == partRepo.declare_part( intersection ) );

  // Test non-trivial intersection:

  intersection.push_back( parts[6] );
  intersection.push_back( parts[7] );

  Part & pint_4_6_7 = * partRepo.declare_part( intersection );

  STKUNIT_ASSERT( 3u == pint_4_6_7.intersection_of().size() );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[4] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[6] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[7] ) );

  STKUNIT_ASSERT( 7u == pint_4_6_7.supersets().size() );

  // Test redeclaration of intersection, should give the same part back:

  STKUNIT_ASSERT( pint_4_6_7 == * partRepo.declare_part( intersection ) );

  partRepo.declare_subset( pint_4_6_7, * parts[8] );

  //--------------------------------------------------------------------
  // Test intersection-induced subset relationship

  partRepo.declare_subset( * parts[7], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  partRepo.declare_subset( * parts[6], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  partRepo.declare_subset( * parts[3], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  partRepo.declare_subset( * parts[4], * parts[10] );
  STKUNIT_ASSERT( contain( pint_4_6_7.subsets() , * parts[10] ) );

  // Test intersection-induced subset relationship triggered from a subset

  partRepo.declare_subset( * parts[7], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );


  partRepo.declare_subset( * parts[6], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  partRepo.declare_subset( * parts[3], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  partRepo.declare_subset( * parts[5], * parts[11] );
  STKUNIT_ASSERT( contain( pint_4_6_7.subsets() , * parts[11] ) );

  std::cout << std::endl << "Part: test intersection generated" << std::endl ;
  print( std::cout , "  " , pint_4_6_7 );
  print( std::cout , "  " , * parts[10] );
  print( std::cout , "  " , * parts[11] );

  std::cout << std::endl ;

  //Test to cover assert_contain in PartRepository - Part is not a subset
  {
    const int spatial_dimension = 1;
    std::vector< std::string > dummy_names(spatial_dimension);
    dummy_names[0].assign("dummy");

    stk::mesh::MetaData meta( dummy_names );

    stk::mesh::Part &new_part4 = meta.declare_part ( "another part");
    meta.commit();

    PartVector intersection2 ;
    intersection2.push_back( &new_part4 );

    impl::PartRepository partRepo2(m);
    STKUNIT_ASSERT_THROW(
        partRepo2.declare_part(intersection2), 
        std::runtime_error
        );
  }

  //Test to cover assert_same_universe in PartRepository - Part is not in the same universe
  {
    int ok = 0 ;
    try {
      impl::PartRepository partRepo2(m);

      PartVector intersection2 ;
      std::vector< std::string > dummy_names(1);
      dummy_names[0].assign("dummy");

      stk::mesh::MetaData meta2( dummy_names );

      stk::mesh::Part &new_part4 = meta2.declare_part ( "another part");
      meta2.commit();

      intersection2.push_back( &new_part4 );

      PartVector::const_iterator i = intersection2.begin() ;
      Part * const p = *i ;
      partRepo2.declare_subset(*parts[5], *p);
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for assert_same_universe in PartRepository");
    }
  }

  //Test to cover assert_not_superset in PartRepository - parts[11] is not a superset of parts[7]
  {
    int ok = 0 ;
    try {
      partRepo.declare_subset(*parts[11], *parts[7] );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for assert_not_superset in PartRepository");
    }
  }

  //Test to cover assert_not_superset in PartRepository - parts[99] is not the same rank of parts[11]
  {
    int ok = 0 ;
    try {
      partRepo.declare_subset(*parts[11], *parts[99] );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for assert_rank_ordering in PartRepository");
    }
  }

  //Test to cover declare_part(arg_name, arg_rank) - Part_99 of rank 1 already exists..
  {
    int ok = 0 ;
    try {
      partRepo.declare_part("Part_99", 0 );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for declare_part in PartRepository");
    }
  }

  //Test to cover declare_part(const PartVector & part_intersect) in PartRepository - failed from malicious abuse
  {
    int ok = 0 ;
    try {
      //create part with intersection

      // Test declaration of an intersection part

      impl::PartRepository partRepo2(m);

      enum { NPARTS = 100 };

      Part * parts2[ NPARTS ] ;

      parts2[0] = & universal ;

      for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
        std::ostringstream name ;
        name << "Part_" << i ;
        parts2[i] = partRepo2.declare_part( name.str() , 0 );
      }

      partRepo2.declare_subset( * parts2[3], * parts2[4] );
      partRepo2.declare_subset( * parts2[4], * parts2[5] );

      partRepo2.declare_subset( * parts2[1], * parts2[2] );
      // 1 and 2 pick up 4 and 5 via transitive relationship:
      partRepo2.declare_subset( * parts2[2], * parts2[3] );

      PartVector intersection2 ;
      intersection2.push_back( parts2[1] );
      intersection2.push_back( parts2[2] );
      intersection2.push_back( parts2[3] );
      intersection2.push_back( parts2[4] ); // Smallest subset of 1..4
      intersection2.push_back( parts2[6] );
      intersection2.push_back( parts2[7] );

      Part & partB = * partRepo2.declare_part("{Part_4^Part_6^Part_7}", 0);

      std::cout << "UnitTestMetaData name of partB is " << partB.name()  << std::endl ;

      Part & pintersect_4_6_7 = * partRepo2.declare_part( intersection2 );

      std::cout << "UnitTestMetaData name of intersection part, pintersect_4_6_7 is " << pintersect_4_6_7.name()  << std::endl ;
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      ///     throw std::runtime_error("UnitTestMetaData FAILED to catch error for declare_part in PartRepository");
    }
  }

  //--------------------------------------------------------------------
  // Test error trapping:

  //bool generated_exception = false ;

  // Can only declare parts on a universal part:

//  generated_exception = false ;
//  try { parts[1]->declare_part( "error" , 0 ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

  // Cannot be circular:

//  generated_exception = false ;
//  try { parts[1]->declare_part( intersection ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

  // Universal part cannot be a subset:

//  generated_exception = false ;
//  try { parts[1]->declare_subset( universal ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

  // Cannot have self-subset

//  generated_exception = false ;
//  try { parts[1]->declare_subset( * parts[1] ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

  // Cannot have circular-subset

//  generated_exception = false ;
//  try { parts[5]->declare_subset( * parts[1] ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

//  Part & part_rank_1 = universal.declare_part( std::string("PartRank1") , 1 );
//  try { parts[1]->declare_subset( part_rank_1 ); }
//  catch( const std::exception & x ) {
//    std::cout << "Part: correctly caught " << x.what() << std::endl ;
//    generated_exception = true ;
//  }
//  STKUNIT_ASSERT( generated_exception );

  //--------------------------------------------------------------------

  return ;
}

//----------------------------------------------------------------------

void UnitTestMetaData::testPartVector()
{
  std::vector< std::string > dummy_names(1);
  dummy_names[0].assign("dummy");

  stk::mesh::MetaData m( dummy_names );

  impl::PartRepository partRepo(&m);

  Part * const pa = partRepo.declare_part( std::string("a") , 0 );
  Part * const pb = partRepo.declare_part( std::string("b") , 0 );
  Part * const pc = partRepo.declare_part( std::string("c") , 0 );
  Part * const pd = partRepo.declare_part( std::string("d") , 0 );
  Part * const pe = partRepo.declare_part( std::string("e") , 0 );
  Part * const pf = partRepo.declare_part( std::string("f") , 0 );

  STKUNIT_ASSERT( ! intersect( *pa , *pb ) );
  STKUNIT_ASSERT( ! intersect( *pb , *pc ) );
  STKUNIT_ASSERT( ! intersect( *pc , *pd ) );
  STKUNIT_ASSERT( ! intersect( *pd , *pe ) );
  STKUNIT_ASSERT( ! intersect( *pe , *pf ) );
  STKUNIT_ASSERT( ! intersect( *pf , *pa ) );

  PartVector vabc , vbcd , vdef , vresult ;

  vabc.push_back( pa );
  vabc.push_back( pb );
  vabc.push_back( pc );

  vbcd.push_back( pb );
  vbcd.push_back( pc );
  vbcd.push_back( pd );

  vdef.push_back( pd );
  vdef.push_back( pe );
  vdef.push_back( pf );

  order( vabc );
  order( vbcd );
  order( vdef );

  vresult.clear();
  STKUNIT_ASSERT_EQUAL( size_t(2) , intersect( vabc , vbcd ) );
  size_t intersect_size = intersect( vabc , vbcd , vresult );
  STKUNIT_ASSERT_EQUAL( size_t(2) , intersect_size );
  STKUNIT_ASSERT_EQUAL( pb , vresult[0] );
  STKUNIT_ASSERT_EQUAL( pc , vresult[1] );

  vresult.clear();
  STKUNIT_ASSERT_EQUAL( size_t(1) , intersect( vdef , vbcd ) );
  intersect_size = intersect( vdef , vbcd , vresult );
  STKUNIT_ASSERT_EQUAL( size_t(1) , intersect_size );
  STKUNIT_ASSERT_EQUAL( pd , vresult[0] );

  vresult.clear();
  STKUNIT_ASSERT_EQUAL( size_t(0) , intersect( vdef , vabc ) );
  intersect_size = intersect( vdef , vabc , vresult );
  STKUNIT_ASSERT_EQUAL( size_t(0) , intersect_size );
  STKUNIT_ASSERT_EQUAL( size_t(0) , vresult.size() );

  Part * const pabc = partRepo.declare_part( std::string("abc") , 0 );
  Part * const pbcd = partRepo.declare_part( std::string("bcd") , 0 );
  Part * const pdef = partRepo.declare_part( std::string("def") , 0 );

  partRepo.declare_subset( * pabc, *pa );
  partRepo.declare_subset( * pabc, *pb );
  partRepo.declare_subset( * pabc, *pc );

  partRepo.declare_subset( * pbcd, *pb );
  partRepo.declare_subset( * pbcd, *pc );
  partRepo.declare_subset( * pbcd, *pd );

  partRepo.declare_subset( * pdef, *pd );
  partRepo.declare_subset( * pdef, *pe );
  partRepo.declare_subset( * pdef, *pf );

  STKUNIT_ASSERT( intersect( *pabc , *pa ) );
  STKUNIT_ASSERT( intersect( *pabc , *pb ) );
  STKUNIT_ASSERT( intersect( *pabc , *pc ) );
  STKUNIT_ASSERT( intersect( *pa , *pabc ) );
  STKUNIT_ASSERT( intersect( *pb , *pabc ) );
  STKUNIT_ASSERT( intersect( *pc , *pabc ) );

  STKUNIT_ASSERT( intersect( *pabc , *pbcd ) );
  STKUNIT_ASSERT( ! intersect( *pabc , *pdef ) );
}

//----------------------------------------------------------------------

namespace {

// Simple tags for testing field dimensions

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( CTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( DTAG )

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( CTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( DTAG )

}

void UnitTestMetaData::testField()
{
  MetaData * const meta_null = NULL ;

  // Declaration of a field allocates one field object
  // per state of the field.  These fields are inserted
  // into a vector of fields of the base class.

  impl::FieldRepository field_repo;
  const FieldVector  & allocated_fields = field_repo.get_fields();

  //------------------------------
  // Declare a double precision scalar field of one state.
  // Not an array; therefore, is rank zero.
  // Test the query methods for accuracy.
  FieldBase * const fA =
    field_repo.declare_field( std::string("A"),
                              data_traits<double>() ,
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              1     /* # States */ ,
                              meta_null );

  STKUNIT_ASSERT( allocated_fields.size() == 1 );
  STKUNIT_ASSERT( fA != NULL );
  STKUNIT_ASSERT( fA == allocated_fields[0] );
  STKUNIT_ASSERT( fA->name() == std::string("A") );
  STKUNIT_ASSERT( fA->type_is<double>() );
  STKUNIT_ASSERT( fA->state() == StateNone );
  STKUNIT_ASSERT( fA->rank()  == 0 );

  //------------------------------
  // Declare a field with an invalid suffix in its name.
  // Suffixes corresponding to "OLD" "N" "NM1" "NM2" "NM3" "NM4"
  // are not allowed as these are automatically appended to
  // the declared variable name for multistate fields.
  {
    FieldBase * tmp = NULL ;
    bool caught = false ;
    try {
      tmp = field_repo.declare_field( "A_OLD" ,
                                      data_traits<double>() ,
                                      0     /* # Ranks */ ,
                                      NULL  /* dimension tags */ ,
                                      1     /* # States */ ,
                                      meta_null );
    }
    catch( const std::exception & x ) {
      std::cout << "Field: Correctly caught " << x.what() << std::endl ;
      caught = true ;
    }
    STKUNIT_ASSERT( caught );
    STKUNIT_ASSERT( tmp == NULL );
    STKUNIT_ASSERT( allocated_fields.size() == 1 );
  }

  //------------------------------
  // Declare a double precision scalar field of two states.
  // Not an array; therefore, is rank zero.
  // Test that two fields: "B" and "B_OLD" were created.
  // Test the query methods for accuracy.

  FieldBase * const fB =
    field_repo.declare_field( std::string("B"),
                              data_traits<int>(),
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              2     /* # States */ ,
                              meta_null );

  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( fB != NULL );
  STKUNIT_ASSERT( fB == allocated_fields[1] );
  STKUNIT_ASSERT( fB->name() == std::string("B") );
  STKUNIT_ASSERT( fB->type_is<int>() );
  STKUNIT_ASSERT( fB->state() == StateNew );
  STKUNIT_ASSERT( fB->rank() == 0 );

  const FieldBase * const fB_old = allocated_fields[2] ;
  STKUNIT_ASSERT( fB_old->name() == std::string("B_OLD") );
  STKUNIT_ASSERT( fB_old->type_is<int>() );
  STKUNIT_ASSERT( fB_old->state() == StateOld );
  STKUNIT_ASSERT( fB_old->rank() == 0 );

  //------------------------------
  // Redeclare field must give back the previous field:

  FieldBase * const fB_redundant =
    field_repo.declare_field( std::string("B"),
                              data_traits<int>(),
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              2     /* # States */ ,
                              meta_null );

  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( fB == fB_redundant );

  //------------------------------
  // Declare a double precision array field of four states.
  // Test that four fields: were created.
  // Test the query methods for accuracy.

  const shards::ArrayDimTag * dim_tags[] =
    { & ATAG::tag() , & BTAG::tag() , & CTAG::tag() , & DTAG::tag() };

  FieldBase * const fC =
    field_repo.declare_field( std::string("C"),
                              data_traits<double>(),
                              3         /* # Ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              4         /* # States */ ,
                              meta_null );

  STKUNIT_ASSERT( allocated_fields.size() == 7 );
  STKUNIT_ASSERT( fC != NULL );
  STKUNIT_ASSERT( fC == allocated_fields[3] );
  STKUNIT_ASSERT( fC->name() == std::string("C") );
  STKUNIT_ASSERT( fC->type_is<double>() );
  STKUNIT_ASSERT( fC->state() == StateNew );
  STKUNIT_ASSERT( fC->rank() == 3 );

  const FieldBase * const fC_n = allocated_fields[4] ;
  const FieldBase * const fC_nm1 = allocated_fields[5] ;
  const FieldBase * const fC_nm2 = allocated_fields[6] ;

  STKUNIT_ASSERT( fC     == fC->field_state( StateNP1 ) );
  STKUNIT_ASSERT( fC_n   == fC->field_state( StateN ) );
  STKUNIT_ASSERT( fC_nm1 == fC->field_state( StateNM1 ) );
  STKUNIT_ASSERT( fC_nm2 == fC->field_state( StateNM2 ) );

  STKUNIT_ASSERT( fC     == fC_n->field_state( StateNP1 ) );
  STKUNIT_ASSERT( fC_n   == fC_n->field_state( StateN ) );
  STKUNIT_ASSERT( fC_nm1 == fC_n->field_state( StateNM1 ) );
  STKUNIT_ASSERT( fC_nm2 == fC_n->field_state( StateNM2 ) );

  STKUNIT_ASSERT( fC     == fC_nm1->field_state( StateNP1 ) );
  STKUNIT_ASSERT( fC_n   == fC_nm1->field_state( StateN ) );
  STKUNIT_ASSERT( fC_nm1 == fC_nm1->field_state( StateNM1 ) );
  STKUNIT_ASSERT( fC_nm2 == fC_nm1->field_state( StateNM2 ) );

  STKUNIT_ASSERT( fC     == fC_nm2->field_state( StateNP1 ) );
  STKUNIT_ASSERT( fC_n   == fC_nm2->field_state( StateN ) );
  STKUNIT_ASSERT( fC_nm1 == fC_nm2->field_state( StateNM1 ) );
  STKUNIT_ASSERT( fC_nm2 == fC_nm2->field_state( StateNM2 ) );

  STKUNIT_ASSERT( fC_n->name() == std::string("C_N") );
  STKUNIT_ASSERT( fC_n->type_is<double>() );
  STKUNIT_ASSERT( fC_n->state() == StateN );
  STKUNIT_ASSERT( fC_n->rank() == 3 );

  STKUNIT_ASSERT( fC_nm1->name() == std::string("C_NM1") );
  STKUNIT_ASSERT( fC_nm1->type_is<double>() );
  STKUNIT_ASSERT( fC_nm1->state() == StateNM1 );
  STKUNIT_ASSERT( fC_nm1->rank() == 3 );

  STKUNIT_ASSERT( fC_nm2->name() == std::string("C_NM2") );
  STKUNIT_ASSERT( fC_nm2->type_is<double>() );
  STKUNIT_ASSERT( fC_nm2->state() == StateNM2 );
  STKUNIT_ASSERT( fC_nm2->rank() == 3 );

  //------------------------------
  // Redeclare field must give back the previous field:
  //------------------------------

  for ( unsigned i = 0 ; i < allocated_fields.size() ; ++i ) {
    FieldBase * const f = allocated_fields[i] ;
    STKUNIT_ASSERT( f->mesh_meta_data_ordinal() == i );
  }

  std::cout << std::endl ;
}

//----------------------------------------------------------------------
// Test field restrictions: the mapping of ( field , part ) -> dimensions

void UnitTestMetaData::testFieldRestriction()
{
  const char method[] = "UnitTestMetaData::testFieldRestriction" ;

  // Arrays for array dimension tags and dimension values
  const shards::ArrayDimTag * dim_tags[] =
    { & ATAG::tag() , & BTAG::tag() , & CTAG::tag() , & DTAG::tag() ,
      & ATAG::tag() , & BTAG::tag() , & CTAG::tag() , & DTAG::tag() };

  unsigned stride[8] ;

  stride[0] = 10 ;
  for ( unsigned i = 1 ; i < 8 ; ++i ) {
    stride[i] = ( i + 1 ) * stride[i-1] ;
  }

  MetaData * const meta_null = NULL ;

  impl::FieldRepository field_repo;
  const FieldVector  & allocated_fields = field_repo.get_fields();

  //------------------------------
  // Declare a rank two and one state:

  FieldBase * const f2 =
    field_repo.declare_field( std::string("F2"),
                              data_traits<int>(),
                              2         /* # ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              1         /* # states */ ,
                              meta_null );

  //------------------------------
  // Declare a rank three and two states:

  FieldBase * const f3 =
    field_repo.declare_field( std::string("F3"),
                              data_traits<int>(),
                              3         /* # ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              2         /* # states */ ,
                              meta_null );

  FieldBase * const f3_old = f3->m_impl.field_state( StateOld ) ;

  //------------------------------
  // Test for correctness of vector of declared fields.
  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( f2 == allocated_fields[0] );
  STKUNIT_ASSERT( f3 == allocated_fields[1] );

  //------------------------------
  // Test for correctness of field internal state access:

  STKUNIT_ASSERT( f2     == f2->field_state( StateNone ) );
  STKUNIT_ASSERT( NULL   == f2->field_state( StateOld ) );
  STKUNIT_ASSERT( f3     == f3->field_state( StateNew ) );
  STKUNIT_ASSERT( f3_old == f3->field_state( StateOld ) );
  STKUNIT_ASSERT( NULL   == f3->field_state( StateNM1 ) );
  STKUNIT_ASSERT( f3     == f3_old->field_state( StateNew ) );
  STKUNIT_ASSERT( f3_old == f3_old->field_state( StateOld ) );
  STKUNIT_ASSERT( NULL   == f3_old->field_state( StateNM1 ) );

  STKUNIT_ASSERT( f2->rank() == 2 );
  STKUNIT_ASSERT( f3->rank() == 3 );
  STKUNIT_ASSERT( f3_old->rank() == 3 );

  //------------------------------
  // Declare some parts for restrictions:

  std::vector< std::string > dummy_names(1);
  dummy_names[0].assign("dummy");

  stk::mesh::MetaData m( dummy_names );

  impl::PartRepository partRepo( &m );
  Part & universal = * partRepo.universal_part();

  Part & pA = * partRepo.declare_part( std::string("A") , 0 );
  Part & pB = * partRepo.declare_part( std::string("B") , 0 );
  Part & pC = * partRepo.declare_part( std::string("C") , 0 );
  Part & pD = * partRepo.declare_part( std::string("D") , 0 );

  // Declare three restrictions:

  f3->m_impl.insert_restriction( method , 0 , pA , stride );
  f3->m_impl.insert_restriction( method , 1 , pB , stride + 1 );
  f3->m_impl.insert_restriction( method , 2 , pC , stride + 2 );

  // Check for correctness of restrictions:

  STKUNIT_ASSERT( f3->restrictions().size() == 3 );
  STKUNIT_ASSERT( f3->restrictions()[0].key ==
                  EntityKey( 0 , pA.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[1].key ==
                  EntityKey( 1 , pB.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[2].key ==
                  EntityKey( 2 , pC.mesh_meta_data_ordinal() ) );

  f3->m_impl.insert_restriction( method , 0 , pB , stride + 1 );

  STKUNIT_ASSERT_EQUAL( f3->max_size( 0 ) , stride[3] );

  //------------------------------
  // Check for error detection of bad stride:
  {
    bool caught = false ;
    try {
      unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
      f3->m_impl.insert_restriction( method , 0 , pA , bad_stride );
    }
    catch( const std::exception & x ) {
      caught = true ;
      std::cout << "Field: Correctly caught: " << x.what() << std::endl ;
    }
    STKUNIT_ASSERT( caught );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Check for error detection in re-declaring an incompatible
  // field restriction.
  {
    bool caught = false ;
    try {
      f3->m_impl.insert_restriction( method , 0 , pA , stride + 1 );
    }
    catch( const std::exception & x ) {
      caught = true ;
      std::cout << "Field: Correctly caught: " << x.what() << std::endl ;
    }
    STKUNIT_ASSERT( caught );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Verify and clean out any redundant restructions:

  f2->m_impl.verify_and_clean_restrictions( method , universal.subsets() );
  f3->m_impl.verify_and_clean_restrictions( method , universal.subsets() );

  STKUNIT_ASSERT( f3->restrictions().size() == 4 );

  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.

  partRepo.declare_subset( pD, pA );
  f2->m_impl.insert_restriction( method , 0 , pA , stride );
  f2->m_impl.insert_restriction( method , 0 , pD , stride );

  STKUNIT_ASSERT( f2->restrictions().size() == 2 );

  f2->m_impl.verify_and_clean_restrictions( method , universal.subsets() );

  STKUNIT_ASSERT( f2->restrictions().size() == 1 );

  {
    const FieldBase::Restriction & rA = f2->restriction( 0 , pA );
    const FieldBase::Restriction & rD = f2->restriction( 0 , pD );
    STKUNIT_ASSERT( & rA == & rD );
    STKUNIT_ASSERT( entity_id( rA.key ) == pD.mesh_meta_data_ordinal() );
  }

  //------------------------------
  // Introduce a new restriction, then introduce a
  // subset-superset relationship that renders the new restriction
  // redundant and incompatible.
  // Check that the verify_and_clean_restrictions method detects
  // this error condition.
  {
    f2->m_impl.insert_restriction( method , 0 , pB , stride + 1 );
    f2->m_impl.verify_and_clean_restrictions( method , universal.subsets() );
    partRepo.declare_subset( pD, pB );
    bool caught = false ;
    try {
      f2->m_impl.verify_and_clean_restrictions( method , universal.subsets() );
    }
    catch( const std::exception & x ) {
      caught = true ;
      std::cout << "Field: Correctly caught: " << x.what() << std::endl ;
    }
    STKUNIT_ASSERT( caught );
  }

  //Test to cover print function in FieldBaseImpl.cpp
  {
    //Create a new field with MetaData m and two restrictions

    FieldBase * const f4 =
      field_repo.declare_field( std::string("F4"),
                                data_traits<int>(),
                                3         /* # ranks */ ,
                                dim_tags  /* dimension tags */ ,
                                2         /* # states */ ,
                                &m );

    partRepo.declare_subset( pD, pA );
    partRepo.declare_subset( pC, pB );

    f4->m_impl.insert_restriction( method , 0 , pA , stride );
    f4->m_impl.insert_restriction( method , 1 , pB , stride + 1 );
    stk::mesh::impl::print(std::cout, "Field f4", *f4);
  }

  //Coverage for error from print_restriction in FieldBaseImpl.cpp when there is no stride (!stride[i])
  //Call print_restriction from insert_restriction
  {
    int ok = 0 ;
    try {
      unsigned arg_stride[2];

      arg_stride[0] = 1;
      arg_stride[1] = 0;

      f2->m_impl.insert_restriction(method, 0, pA, arg_stride);

    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for iterator from print_restriction inFieldBaseImpl.cpp");
    }
  }

  std::cout << std::endl ;
}

//----------------------------------------------------------------------

void UnitTestMetaData::testProperty()
{
  std::vector< std::string > dummy_names(1);
  dummy_names[0].assign("dummy");

  stk::mesh::MetaData meta_data( dummy_names );
  stk::mesh::MetaData meta_data2( dummy_names );

  stk::mesh::Property<int> & pi = meta_data.declare_property<int>("my_i");
  stk::mesh::Property<double> & pf = meta_data.declare_property<double>("my_d");
  stk::mesh::Property<double> & px = meta_data.declare_property<double>("my_x");
  stk::mesh::Property<double> & pa = meta_data.declare_property<double>("my_array",5);

  STKUNIT_ASSERT( pi.type_is<int>() );
  STKUNIT_ASSERT( pf.type_is<double>() );
  STKUNIT_ASSERT( px.type_is<double>() );
  STKUNIT_ASSERT( pa.type_is<double>() );

  STKUNIT_ASSERT( ! pi.type_is<double>() );
  STKUNIT_ASSERT( ! pa.type_is<int>() );

  STKUNIT_ASSERT_EQUAL( pi.size() , 1u );
  STKUNIT_ASSERT_EQUAL( pf.size() , 1u );
  STKUNIT_ASSERT_EQUAL( px.size() , 1u );
  STKUNIT_ASSERT_EQUAL( pa.size() , 5u );

  meta_data.put_property( pi , meta_data.locally_owned_part() );

  STKUNIT_ASSERT( stk::mesh::property_data( pi , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT( ! stk::mesh::property_data( px , meta_data.locally_owned_part() ) );

  stk::mesh::Property<int> * property();

  const stk::mesh::PropertyBase * p = NULL;
  property_data_throw( *p, meta_data.locally_owned_part());

  ///////stk::mesh::Part & part_1 = meta_data.declare_part( "a_part_1", 0 );

  const std::string& str = "my_i";
  const std::string& str2 = "my_d";

  stk::mesh::Property<int> * const pProp = meta_data.get_property<int>( str );
  STKUNIT_ASSERT( (*pProp).type_is<int>() );

  {
    int ok = 0 ;
    try {
      meta_data.get_property<double>( str );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for get_property");
    }
  }

  {
    int ok = 0 ;
    try {
      meta_data.get_property<int>( str2 );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for get_property");
    }
  }

  //Final coverage of MetaData.hpp - declare_property
  stk::mesh::Property<double> & pb = meta_data.declare_property<double>("my_array",0);
  STKUNIT_ASSERT( (pb).type_is<double>() );

  //More coverage of Property.hpp
  STKUNIT_ASSERT( stk::mesh::property_data( pi , meta_data2.locally_owned_part() ));
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {
int stencil_test_function( unsigned  from_type ,
                           unsigned  to_type ,
                           unsigned  identifier )
{
  return 0;
}

}

namespace stk {
namespace mesh {

void UnitTestMetaData::testMetaData()
{
  const int spatial_dimension = 3;
  const std::vector<std::string> & type_names = TopologicalMetaData::entity_rank_names(spatial_dimension);

  MetaData metadata2( type_names );
  MetaData metadata3( type_names );
  MetaData metadata4( type_names );

  STKUNIT_ASSERT_THROW(metadata3.assert_committed("test throw"),std::logic_error);
  metadata3.commit();
  STKUNIT_ASSERT_THROW(metadata3.assert_not_committed("test throw"),std::logic_error);
  STKUNIT_ASSERT_THROW(metadata3.assert_same_mesh_meta_data("test throw",metadata2),std::logic_error);
  std::string test_string = "this_part_does_not_exist";
  STKUNIT_ASSERT_THROW(metadata3.get_part(test_string,"test_throw"),std::runtime_error);

  //test declare part
  Part * const pa = & metadata4.declare_part( std::string("a") , 0 );
  Part * const pb = & metadata4.declare_part( std::string("b") , 0 );
  Part * const pc = & metadata4.declare_part( std::string("c") , 0 );
  Part * const pd = & metadata4.declare_part( std::string("d") , 0 );
  Part * const pe = & metadata4.declare_part( std::string("e") , 0 );
  Part * const pf = & metadata4.declare_part( std::string("f") , 0 );
  Part * const pg = & metadata4.declare_part( std::string("g") , 0 );
  Part * const ph = & metadata4.declare_part( std::string("h") , 0 );

  metadata4.declare_part_relation(*pe,stencil_test_function,*pg);

  PartVector part_vector;

  part_vector.push_back(pa);
  part_vector.push_back(pb);
  part_vector.push_back(pc);
  part_vector.push_back(pd);

  //Part * const intersection_part = &
  metadata4.declare_part(part_vector);

  STKUNIT_ASSERT_THROW( metadata4.declare_part_subset(*pe,*pe), std::runtime_error);
  metadata4.declare_part_subset(*pd,*pf);

  STKUNIT_ASSERT_THROW( metadata4.declare_part_relation(*pg,stencil_test_function,*ph), std::runtime_error);
  STKUNIT_ASSERT_THROW( metadata4.declare_part_relation(*pe,NULL,*pe), std::runtime_error);

  {
    bool caught_throw = false;
      try {
        std::vector<std::string> empty_names;
        MetaData metadata5(empty_names);
      }
    catch(...) {
      caught_throw = true;
    }
    STKUNIT_ASSERT_EQUAL(caught_throw, true);
  }

  int i = 2;

  const std::string& i_name2 = metadata2.entity_rank_name( i );

  STKUNIT_ASSERT( i_name2 == type_names[i] );

  EntityRank one_rank_higher_than_defined = type_names.size();
  STKUNIT_ASSERT_THROW( 
    metadata2.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
    );
}

void UnitTestMetaData::testEntityRepository()
{
  //Test Entity repository - covering EntityRepository.cpp/hpp

  stk::mesh::MetaData meta ( stk::mesh::fem_entity_rank_names() );
  stk::mesh::Part & part = meta.declare_part( "another part");

  meta.commit();

  stk::mesh::BulkData bulk ( meta , MPI_COMM_WORLD , 100 );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &part );

  int  size , rank;
  rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  size = stk::parallel_machine_size( MPI_COMM_WORLD );
  PartVector tmp(1);

  bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( 0 , new_id+1 , add_part );
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity & elem  = bulk.declare_entity( 3 , new_id+1 , add_part );

  //new_id = size * (++id_base) + rank;
  // stk::mesh::Entity & elem2  = bulk.declare_entity( 3 , new_id+1 , add_part );

  stk::mesh::impl::EntityRepository e;

  e.comm_clear( elem );

  e.comm_clear_ghosting( elem );

  const stk::mesh::Ghosting & ghost = bulk.shared_aura();

  bulk.modification_end();

  STKUNIT_ASSERT_FALSE(e.erase_ghosting(elem, ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  STKUNIT_ASSERT_FALSE(e.erase_comm_info(elem, comm_info));

  STKUNIT_ASSERT(e.insert_comm_info(elem, comm_info));

  //Checking internal_create_entity

  e.internal_create_entity( stk::mesh::EntityKey( 3, 2 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 5 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 7 ));

  //Checking get_entity with invalid key - no rank or id
  {
    int ok = 0 ;
    try {

      stk::mesh::Entity * elem3 = e.get_entity(stk::mesh::EntityKey());
      if(elem3){
        // CAROL FIXME
      }

    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }
    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for get_entity - invalid key");
    }
  }
}

//----------------------------------------------------------------------

}
}

