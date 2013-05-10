/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.               */
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
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {

class UnitTestFieldImpl {
public:
  UnitTestFieldImpl() {}

  void testFieldRestriction();

};

}//namespace mesh
}//namespace stk

namespace {

STKUNIT_UNIT_TEST(UnitTestFieldRestriction, testUnit)
{
  stk::mesh::UnitTestFieldImpl ufield;
  ufield.testFieldRestriction();
}

}//namespace <anonymous>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

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

//----------------------------------------------------------------------
// Test field restrictions: the mapping of ( field , part ) -> dimensions

void UnitTestFieldImpl::testFieldRestriction()
{
  unsigned stride[8] ;

  stride[0] = 10 ;
  for ( unsigned i = 1 ; i < 8 ; ++i ) {
    stride[i] = ( i + 1 ) * stride[i-1] ;
  }

  std::vector< std::string > dummy_names(4, "dummy");

  MetaData meta_data(0 /*dim*/,dummy_names);

  const FieldVector  & allocated_fields = meta_data.get_fields();

  //------------------------------

  typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorField;

  FieldBase * const f2 =
    &meta_data.declare_field<VectorField>( std::string("F2"), 1/* # states */ );

  //------------------------------

  FieldBase * const f3 =
    &meta_data.declare_field<VectorField>( std::string("F3"), 2/* #states*/);

  FieldBase * const f3_old = f3->field_state( StateOld ) ;

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

  //------------------------------
  // Declare some parts for restrictions:

  Part & pA = meta_data.declare_part( std::string("A") , 0 );
  Part & pB = meta_data.declare_part( std::string("B") , 0 );
  Part & pC = meta_data.declare_part( std::string("C") , 0 );
  Part & pD = meta_data.declare_part( std::string("D") , 0 );

  // Declare three restrictions:

  meta_data.declare_field_restriction(*f3, 0 , pA , stride );
  meta_data.declare_field_restriction(*f3, 1 , pB , stride + 1 );
  meta_data.declare_field_restriction(*f3, 2 , pC , stride + 2 );

  // Check for correctness of restrictions:

  STKUNIT_ASSERT( f3->restrictions().size() == 3 );
  STKUNIT_ASSERT( f3->restrictions()[0] ==
                  FieldRestriction( 0 , pA.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[1] ==
                  FieldRestriction( 1 , pB.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[2] ==
                  FieldRestriction( 2 , pC.mesh_meta_data_ordinal() ) );

  meta_data.declare_field_restriction(*f3, 0 , pB , stride + 1 );

  STKUNIT_ASSERT_EQUAL( f3->max_size( 0 ) , 20u );

  //------------------------------
  // Check for error detection of bad stride:
  {
    unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
    STKUNIT_ASSERT_THROW(
      meta_data.declare_field_restriction(*f3, 0 , pA , bad_stride ),
      std::runtime_error
    );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Check for error detection in re-declaring an incompatible
  // field restriction.
  {
    STKUNIT_ASSERT_THROW(
      meta_data.declare_field_restriction(*f3, 0 , pA , stride + 1 ),
      std::runtime_error
    );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Verify and clean out any redundant restructions:

  STKUNIT_ASSERT( f3->restrictions().size() == 4 );

  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.

  std::cout<<"pA ord: "<<pA.mesh_meta_data_ordinal()<<", pD ord: "<<pD.mesh_meta_data_ordinal()<<std::endl;
  meta_data.declare_part_subset( pD, pA );
  meta_data.declare_field_restriction(*f2, 0 , pA , stride );
  meta_data.declare_field_restriction(*f2, 0 , pD , stride );

  STKUNIT_ASSERT( f2->restrictions().size() == 1 );

  {
    const FieldBase::Restriction & rA = f2->restriction( 0 , pA );
    const FieldBase::Restriction & rD = f2->restriction( 0 , pD );
    STKUNIT_ASSERT( & rA == & rD );
    STKUNIT_ASSERT( rA.part_ordinal() == pD.mesh_meta_data_ordinal() );
  }

  //------------------------------
  // Introduce a new restriction, then introduce a
  // subset-superset relationship that renders the new restriction
  // redundant and incompatible.
  // Check that the verify_and_clean_restrictions method detects
  // this error condition.
  {
    meta_data.declare_field_restriction(*f2, 0 , pB , stride + 1 );
    STKUNIT_ASSERT_THROW(
      meta_data.declare_part_subset( pD, pB ),
      std::runtime_error
    );
  }

  //Test to cover print function in FieldBaseImpl.cpp and FieldBase.cpp
  {
    //Create a new field with MetaData m and two restrictions

    FieldBase * const f4 =
      &meta_data.declare_field<VectorField>( std::string("F4"),
                                2         /* # states */ );

    meta_data.declare_part_subset( pD, pA );
    meta_data.declare_part_subset( pC, pB );

    meta_data.declare_field_restriction(*f4, 0 , pA , stride );
    meta_data.declare_field_restriction(*f4, 1 , pB , stride + 1 );
    stk::mesh::impl::print(std::cout, "Field f4", *f4);

    //test stride[i] / stride[i-1] section of else-if
    stk::mesh::print(std::cout, "Field f4", *f4);
  }

  //Further tests to cover print function in FieldBase.cpp
  {
    //test stride[i] % stride[i-1] section of else-if

    //Create a new field with MetaData m and two restrictions

    FieldBase * const f5 =
      &meta_data.declare_field<VectorField>( std::string("F5"),
                                2         /* # states */ );

    unsigned stride2[8] ;
    stride2[0] = 10 ;
    for ( unsigned i = 1 ; i < 3 ; ++i ) {
      stride2[i] = stride[i-1];
    }
    for ( unsigned i = 3 ; i < 8 ; ++i ) {
      stride2[i] = 0;
    }
    meta_data.declare_field_restriction(*f5, 0 , pA, stride2 );

    stk::mesh::print(std::cout, "Field f5", *f5);

  }

  //Coverage for error from print_restriction in FieldBaseImpl.cpp when there is no stride (!stride[i])
  //Call print_restriction from insert_restriction
  {
    unsigned arg_no_stride[2];

    arg_no_stride[0] = 1;
    arg_no_stride[1] = 0;

    STKUNIT_ASSERT_THROW(
      meta_data.declare_field_restriction(*f2, 0, pA, arg_no_stride),
      std::runtime_error
    );
  }

  //Coverage of ordinal in FieldRestriction.hpp:
  {
    const FieldRestrictionVector & rMap = f3->restrictions();
    const FieldRestrictionVector::const_iterator ie = rMap.end() ;
          FieldRestrictionVector::const_iterator i = rMap.begin();

    EntityId entity_id = 0;
    unsigned max = 0 ;

    for ( ; i != ie ; ++i ) {
      if ( i->part_ordinal() == entity_id ) {
	const unsigned len = pA.mesh_meta_data_ordinal() ? i->stride( pA.mesh_meta_data_ordinal() - 1 ) : 1 ;
        if ( max < len ) { max = len ; }
      }
    }
  }
}


}
}

