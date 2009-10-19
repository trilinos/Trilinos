
#include <sstream>

#include <unit_tests/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

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

}//namespace <anonymous>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
void UnitTestMetaData::testEntityKey()
{
  enum { OK = StaticAssert<sizeof(EntityKey) == sizeof(EntityKeyValue)>::OK };

  EntityKey
    key_bad_zero  = EntityKey(),
    key_bad_id    = EntityKey( 0 , 0 ),
    key_good_0_1  = EntityKey( 0 , 1 ),
    key_good_1_1  = EntityKey( 1 , 1 ),
    key_good_2_10 = EntityKey( 2 , 10 );

//   EntityId val = entity_id(key_good_2_10) ;

//   EntityKey good_key_gone_bad;
//   entity_good_key_gone_bad.key = val ;

  STKUNIT_ASSERT( ! entity_key_valid( key_bad_zero ) );
  STKUNIT_ASSERT( ! entity_key_valid( key_bad_id ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_0_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_1_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( key_good_2_10 ) );

  STKUNIT_ASSERT( 0  == entity_type(key_good_0_1));
  STKUNIT_ASSERT( 1  == entity_type(key_good_1_1) );
  STKUNIT_ASSERT( 2  == entity_type(key_good_2_10) );
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
//             << ": key_bad_id = "    << key_bad_id
//             << ", key_good_2_10 = " << key_good_2_10.key
//             << ", good_key_gone_bad = " << good_key_gone_bad.key
//             << std::endl ;
//   std::cout.unsetf( std::ios::hex);
//   std::cout.setf( std::ios::dec);

  EntityKey key01(0,1), key_default;
  key01 = key_default;
}

//----------------------------------------------------------------------
// Unit test the Part functionality in isolation:

void UnitTestMetaData::testPart()
{
  MetaData * const m = NULL ;

  Part universal( m );

  STKUNIT_ASSERT( universal.supersets().empty() );
  STKUNIT_ASSERT( 1u == universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL( universal.subsets()[0] , & universal );

  //--------------------------------------------------------------------
  // Test multiple part creation

  enum { NPARTS = 100 };

  Part * parts[ NPARTS ] ;

  parts[0] = & universal ;

  for ( int i = 1 ; i < NPARTS ; ++i ) {
    std::ostringstream name ;
    name << "Part_" << i ;
    parts[i] = & universal.declare_part( name.str() , 0 );
  }

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

  parts[3]->declare_subset( * parts[4] );
  parts[4]->declare_subset( * parts[5] );

  parts[1]->declare_subset( * parts[2] );
  // 1 and 2 pick up 4 and 5 via transitive relationship:
  parts[2]->declare_subset( * parts[3] );

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

  STKUNIT_ASSERT( * parts[4] == universal.declare_part( intersection ) );

  // Test non-trivial intersection:

  intersection.push_back( parts[6] );
  intersection.push_back( parts[7] );

  Part & pint_4_6_7 = universal.declare_part( intersection );

  STKUNIT_ASSERT( 3u == pint_4_6_7.intersection_of().size() );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[4] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[6] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[7] ) );

  STKUNIT_ASSERT( 7u == pint_4_6_7.supersets().size() );

  // Test redeclaration of intersection, should give the same part back:

  STKUNIT_ASSERT( pint_4_6_7 == universal.declare_part( intersection ) );

  pint_4_6_7.declare_subset( * parts[8] );

  //--------------------------------------------------------------------
  // Test intersection-induced subset relationship

  parts[7]->declare_subset( * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  parts[6]->declare_subset( * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  parts[3]->declare_subset( * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  parts[4]->declare_subset( * parts[10] );
  STKUNIT_ASSERT( contain( pint_4_6_7.subsets() , * parts[10] ) );

  // Test intersection-induced subset relationship triggered from a subset

  parts[7]->declare_subset( * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  parts[6]->declare_subset( * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  parts[3]->declare_subset( * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  parts[5]->declare_subset( * parts[11] );
  STKUNIT_ASSERT( contain( pint_4_6_7.subsets() , * parts[11] ) );

  std::cout << std::endl << "Part: test intersection generated" << std::endl ;
  print( std::cout , "  " , pint_4_6_7 );
  print( std::cout , "  " , * parts[10] );
  print( std::cout , "  " , * parts[11] );

  //--------------------------------------------------------------------
  // Test error trapping:

  bool generated_exception = false ;

  // Can only declare parts on a universal part:

  generated_exception = false ;
  try { parts[1]->declare_part( "error" , 0 ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  // Cannot be circular:

  generated_exception = false ;
  try { parts[1]->declare_part( intersection ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  // Universal part cannot be a subset:

  generated_exception = false ;
  try { parts[1]->declare_subset( universal ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  // Cannot have self-subset

  generated_exception = false ;
  try { parts[1]->declare_subset( * parts[1] ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  // Cannot have circular-subset

  generated_exception = false ;
  try { parts[5]->declare_subset( * parts[1] ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  Part & part_rank_1 = universal.declare_part( std::string("PartRank1") , 1 );
  try { parts[1]->declare_subset( part_rank_1 ); }
  catch( const std::exception & x ) {
    std::cout << "Part: correctly caught " << x.what() << std::endl ;
    generated_exception = true ;
  }
  STKUNIT_ASSERT( generated_exception );

  //--------------------------------------------------------------------

  {
    Part bad_universal( m );
    Part & badA = bad_universal.declare_part( std::string("badA") , 0 );

    generated_exception = false ;
    const size_t part_0_nsub = parts[1]->subsets().size();
    try {
      parts[1]->declare_subset( badA );
    }
    catch( const std::exception & x ) {
      std::cout << "Part: correctly caught " << x.what() << std::endl ;
      generated_exception = true ;
    }
    STKUNIT_ASSERT( generated_exception );
    STKUNIT_ASSERT_EQUAL( part_0_nsub , parts[1]->subsets().size() );

    generated_exception = false ;
    try {
      bad_universal.declare_part( intersection );
    }
    catch( const std::exception & x ) {
      std::cout << "Part: correctly caught " << x.what() << std::endl ;
      generated_exception = true ;
    }
    STKUNIT_ASSERT( generated_exception );
  }

  std::cout << std::endl ;

  //--------------------------------------------------------------------

  return ;
}

//----------------------------------------------------------------------

void UnitTestMetaData::testPartVector()
{
  MetaData * const m = NULL ;

  Part universal( m );

  Part * const pa = & universal.declare_part( std::string("a") , 0 );
  Part * const pb = & universal.declare_part( std::string("b") , 0 );
  Part * const pc = & universal.declare_part( std::string("c") , 0 );
  Part * const pd = & universal.declare_part( std::string("d") , 0 );
  Part * const pe = & universal.declare_part( std::string("e") , 0 );
  Part * const pf = & universal.declare_part( std::string("f") , 0 );

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
  STKUNIT_ASSERT_EQUAL( size_t(2) , intersect( vabc , vbcd , vresult ) );
  STKUNIT_ASSERT_EQUAL( pb , vresult[0] );
  STKUNIT_ASSERT_EQUAL( pc , vresult[1] );

  vresult.clear();
  STKUNIT_ASSERT_EQUAL( size_t(1) , intersect( vdef , vbcd ) );
  STKUNIT_ASSERT_EQUAL( size_t(1) , intersect( vdef , vbcd , vresult ) );
  STKUNIT_ASSERT_EQUAL( pd , vresult[0] );

  vresult.clear();
  STKUNIT_ASSERT_EQUAL( size_t(0) , intersect( vdef , vabc ) );
  STKUNIT_ASSERT_EQUAL( size_t(0) , intersect( vdef , vabc , vresult ) );
  STKUNIT_ASSERT_EQUAL( size_t(0) , vresult.size() );

  Part * const pabc = & universal.declare_part( std::string("abc") , 0 );
  Part * const pbcd = & universal.declare_part( std::string("bcd") , 0 );
  Part * const pdef = & universal.declare_part( std::string("def") , 0 );

  pabc->declare_subset( *pa );
  pabc->declare_subset( *pb );
  pabc->declare_subset( *pc );

  pbcd->declare_subset( *pb );
  pbcd->declare_subset( *pc );
  pbcd->declare_subset( *pd );

  pdef->declare_subset( *pd );
  pdef->declare_subset( *pe );
  pdef->declare_subset( *pf );

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

  std::vector<FieldBase*> allocated_fields ;

  //------------------------------
  // Declare a double precision scalar field of one state.
  // Not an array; therefore, is rank zero.
  // Test the query methods for accuracy.
  FieldBase * const fA =
    FieldBase::declare_field( std::string("A"),
                              data_traits<double>() ,
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              1     /* # States */ ,
                              meta_null , allocated_fields );

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
      tmp = FieldBase::declare_field( "A_OLD" ,
                                      data_traits<double>() ,
                                      0     /* # Ranks */ ,
                                      NULL  /* dimension tags */ ,
                                      1     /* # States */ ,
                                      meta_null , allocated_fields );
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
    FieldBase::declare_field( std::string("B"),
                              data_traits<int>(),
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              2     /* # States */ ,
                              meta_null , allocated_fields );

  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( fB != NULL );
  STKUNIT_ASSERT( fB == allocated_fields[1] );
  STKUNIT_ASSERT( fB->name() == std::string("B") );
  STKUNIT_ASSERT( fB->type_is<int>() );
  STKUNIT_ASSERT( fB->state() == StateNew );
  STKUNIT_ASSERT( fB->rank() == 0 );

  FieldBase * const fB_old = allocated_fields[2] ;
  STKUNIT_ASSERT( fB_old->name() == std::string("B_OLD") );
  STKUNIT_ASSERT( fB_old->type_is<int>() );
  STKUNIT_ASSERT( fB_old->state() == StateOld );
  STKUNIT_ASSERT( fB_old->rank() == 0 );

  //------------------------------
  // Redeclare field must give back the previous field:

  FieldBase * const fB_redundant =
    FieldBase::declare_field( std::string("B"),
                              data_traits<int>(),
                              0     /* # Ranks */ ,
                              NULL  /* dimension tags */ ,
                              2     /* # States */ ,
                              meta_null , allocated_fields );

  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( fB == fB_redundant );

  //------------------------------
  // Declare a double precision array field of four states.
  // Test that four fields: were created.
  // Test the query methods for accuracy.

  const shards::ArrayDimTag * dim_tags[] =
    { & ATAG::tag() , & BTAG::tag() , & CTAG::tag() , & DTAG::tag() };

  FieldBase * const fC =
    FieldBase::declare_field( std::string("C"),
                              data_traits<double>(),
                              3         /* # Ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              4         /* # States */ ,
                              meta_null , allocated_fields );

  STKUNIT_ASSERT( allocated_fields.size() == 7 );
  STKUNIT_ASSERT( fC != NULL );
  STKUNIT_ASSERT( fC == allocated_fields[3] );
  STKUNIT_ASSERT( fC->name() == std::string("C") );
  STKUNIT_ASSERT( fC->type_is<double>() );
  STKUNIT_ASSERT( fC->state() == StateNew );
  STKUNIT_ASSERT( fC->rank() == 3 );

  FieldBase * const fC_n = allocated_fields[4] ;
  FieldBase * const fC_nm1 = allocated_fields[5] ;
  FieldBase * const fC_nm2 = allocated_fields[6] ;

  STKUNIT_ASSERT( fC     == fC->m_field_states[ StateNP1 ] );
  STKUNIT_ASSERT( fC_n   == fC->m_field_states[ StateN ] );
  STKUNIT_ASSERT( fC_nm1 == fC->m_field_states[ StateNM1 ] );
  STKUNIT_ASSERT( fC_nm2 == fC->m_field_states[ StateNM2 ] );

  STKUNIT_ASSERT( fC     == fC_n->m_field_states[ StateNP1 ] );
  STKUNIT_ASSERT( fC_n   == fC_n->m_field_states[ StateN ] );
  STKUNIT_ASSERT( fC_nm1 == fC_n->m_field_states[ StateNM1 ] );
  STKUNIT_ASSERT( fC_nm2 == fC_n->m_field_states[ StateNM2 ] );

  STKUNIT_ASSERT( fC     == fC_nm1->m_field_states[ StateNP1 ] );
  STKUNIT_ASSERT( fC_n   == fC_nm1->m_field_states[ StateN ] );
  STKUNIT_ASSERT( fC_nm1 == fC_nm1->m_field_states[ StateNM1 ] );
  STKUNIT_ASSERT( fC_nm2 == fC_nm1->m_field_states[ StateNM2 ] );

  STKUNIT_ASSERT( fC     == fC_nm2->m_field_states[ StateNP1 ] );
  STKUNIT_ASSERT( fC_n   == fC_nm2->m_field_states[ StateN ] );
  STKUNIT_ASSERT( fC_nm1 == fC_nm2->m_field_states[ StateNM1 ] );
  STKUNIT_ASSERT( fC_nm2 == fC_nm2->m_field_states[ StateNM2 ] );

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
  // Delete the allocated fields.

  for ( unsigned i = 0 ; i < allocated_fields.size() ; ++i ) {
    FieldBase * const f = allocated_fields[i] ;
    STKUNIT_ASSERT( f->mesh_meta_data_ordinal() == i );
    delete f ;
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

  std::vector<FieldBase*> allocated_fields ;

  //------------------------------
  // Declare a rank two and one state:

  FieldBase * const f2 =
    FieldBase::declare_field( std::string("F2"),
                              data_traits<int>(),
                              2         /* # ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              1         /* # states */ ,
                              meta_null , allocated_fields );

  //------------------------------
  // Declare a rank three and two states:

  FieldBase * const f3 =
    FieldBase::declare_field( std::string("F3"),
                              data_traits<int>(),
                              3         /* # ranks */ ,
                              dim_tags  /* dimension tags */ ,
                              2         /* # states */ ,
                              meta_null , allocated_fields );

  FieldBase * const f3_old = f3->m_field_states[ StateOld ] ;

  //------------------------------
  // Test for correctness of vector of declared fields.
  STKUNIT_ASSERT( allocated_fields.size() == 3 );
  STKUNIT_ASSERT( f2 == allocated_fields[0] );
  STKUNIT_ASSERT( f3 == allocated_fields[1] );

  //------------------------------
  // Test for correctness of field internal state access:

  STKUNIT_ASSERT( f2     == f2->m_field_states[ StateNone ] );
  STKUNIT_ASSERT( NULL   == f2->m_field_states[ StateOld ] );
  STKUNIT_ASSERT( f3     == f3->m_field_states[ StateNew ] );
  STKUNIT_ASSERT( f3_old == f3->m_field_states[ StateOld ] );
  STKUNIT_ASSERT( NULL   == f3->m_field_states[ StateNM1 ] );
  STKUNIT_ASSERT( f3     == f3_old->m_field_states[ StateNew ] );
  STKUNIT_ASSERT( f3_old == f3_old->m_field_states[ StateOld ] );
  STKUNIT_ASSERT( NULL   == f3_old->m_field_states[ StateNM1 ] );

  STKUNIT_ASSERT( f2->rank() == 2 );
  STKUNIT_ASSERT( f3->rank() == 3 );
  STKUNIT_ASSERT( f3_old->rank() == 3 );

  //------------------------------
  // Declare some parts for restrictions:

  Part universal( meta_null );

  Part & pA = universal.declare_part( std::string("A") , 0 );
  Part & pB = universal.declare_part( std::string("B") , 0 );
  Part & pC = universal.declare_part( std::string("C") , 0 );
  Part & pD = universal.declare_part( std::string("D") , 0 );

  // Declare three restrictions:

  f3->insert_restriction( method , 0 , pA , stride );
  f3->insert_restriction( method , 1 , pB , stride + 1 );
  f3->insert_restriction( method , 2 , pC , stride + 2 );

  // Check for correctness of restrictions:

  STKUNIT_ASSERT( f3->restrictions().size() == 3 );
  STKUNIT_ASSERT( f3->restrictions()[0].key ==
                  EntityKey( 0 , pA.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[1].key ==
                  EntityKey( 1 , pB.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[2].key ==
                  EntityKey( 2 , pC.mesh_meta_data_ordinal() ) );

  f3->insert_restriction( method , 0 , pB , stride + 1 );

  STKUNIT_ASSERT_EQUAL( f3->max_size( 0 ) , stride[3] );

  //------------------------------
  // Check for error detection of bad stride:
  {
    bool caught = false ;
    try {
      unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
      f3->insert_restriction( method , 0 , pA , bad_stride );
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
      f3->insert_restriction( method , 0 , pA , stride + 1 );
    }
    catch( const std::exception & x ) {
      caught = true ;
      std::cout << "Field: Correctly caught: " << x.what() << std::endl ;
    }
    STKUNIT_ASSERT( caught );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Verify and clean out any redundant restructions:

  f2->verify_and_clean_restrictions( method , universal.subsets() );
  f3->verify_and_clean_restrictions( method , universal.subsets() );

  STKUNIT_ASSERT( f3->restrictions().size() == 4 );

  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.

  pD.declare_subset( pA );
  f2->insert_restriction( method , 0 , pA , stride );
  f2->insert_restriction( method , 0 , pD , stride );

  STKUNIT_ASSERT( f2->restrictions().size() == 2 );

  f2->verify_and_clean_restrictions( method , universal.subsets() );

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
    f2->insert_restriction( method , 0 , pB , stride + 1 );
    f2->verify_and_clean_restrictions( method , universal.subsets() );
    pD.declare_subset( pB );
    bool caught = false ;
    try {
      f2->verify_and_clean_restrictions( method , universal.subsets() );
    }
    catch( const std::exception & x ) {
      caught = true ;
      std::cout << "Field: Correctly caught: " << x.what() << std::endl ;
    }
    STKUNIT_ASSERT( caught );
  }

  //------------------------------
  // Delete allocated fields:

  for ( unsigned i = 0 ; i < allocated_fields.size() ; ++i ) {
    delete allocated_fields[i] ;
  }

  std::cout << std::endl ;
}

//----------------------------------------------------------------------

void UnitTestMetaData::testProperty()
{
  std::vector< std::string > dummy_names(1);
  dummy_names[0].assign("dummy");

  stk::mesh::MetaData meta_data( dummy_names );

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
  meta_data.put_property( px , meta_data.locally_used_part() );

  STKUNIT_ASSERT( stk::mesh::property_data( pi , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT( stk::mesh::property_data( px , meta_data.locally_used_part() ) != NULL);
  STKUNIT_ASSERT( ! stk::mesh::property_data( pi , meta_data.locally_used_part() ) );
  STKUNIT_ASSERT( ! stk::mesh::property_data( px , meta_data.locally_owned_part() ) );
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace{
int stencil_test_function( unsigned  from_type ,
                                        unsigned  to_type ,
                                        unsigned  identifier ,
                                        unsigned  kind )
{
  return 0;
}

}


namespace stk {
namespace mesh {

void UnitTestMetaData::testMetaData()
{
  const std::vector<std::string> & type_names = fem_entity_type_names();

  STKUNIT_ASSERT_EQUAL( type_names.size() , (size_t) EntityTypeEnd );

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

  const std::string& i_name2 = metadata2.entity_type_name( i );

  STKUNIT_ASSERT( i_name2 == type_names[i] );

  i = EntityTypeEnd;
  bool caught_throw = false;
  try {
    metadata2.entity_type_name( i );
  }
  catch(...) {
    caught_throw = true;
  }

  STKUNIT_ASSERT_EQUAL(caught_throw, true);
}

//----------------------------------------------------------------------

}
}

