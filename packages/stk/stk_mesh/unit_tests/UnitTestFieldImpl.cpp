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
#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/PartRelation.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>

#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {

class UnitTestFieldImpl {
public:
  UnitTestFieldImpl() {}

  void testField();
  void testFieldRestriction();
  void testFieldRelation();

};

}//namespace mesh
}//namespace stk

namespace {

STKUNIT_UNIT_TEST(UnitTestField, testUnit)
{
  stk::mesh::UnitTestFieldImpl ufield;
  ufield.testField();
}

STKUNIT_UNIT_TEST(UnitTestFieldRestriction, testUnit)
{
  stk::mesh::UnitTestFieldImpl ufield;
  ufield.testFieldRestriction();
}

STKUNIT_UNIT_TEST(UnitTestFieldRelation, testUnit)
{
  stk::mesh::UnitTestFieldImpl ufield;
  ufield.testFieldRelation();
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

void UnitTestFieldImpl::testField()
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
    STKUNIT_ASSERT_THROW(
      field_repo.declare_field( "A_STKFS_OLD" ,
                                data_traits<double>() ,
                                0     /* # Ranks */ ,
                                NULL  /* dimension tags */ ,
                                1     /* # States */ ,
                                meta_null ),
      std::runtime_error);
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
  STKUNIT_ASSERT( fB_old->name() == std::string("B_STKFS_OLD") );
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

  STKUNIT_ASSERT( fC_n->name() == std::string("C_STKFS_N") );
  STKUNIT_ASSERT( fC_n->type_is<double>() );
  STKUNIT_ASSERT( fC_n->state() == StateN );
  STKUNIT_ASSERT( fC_n->rank() == 3 );

  STKUNIT_ASSERT( fC_nm1->name() == std::string("C_STKFS_NM1") );
  STKUNIT_ASSERT( fC_nm1->type_is<double>() );
  STKUNIT_ASSERT( fC_nm1->state() == StateNM1 );
  STKUNIT_ASSERT( fC_nm1->rank() == 3 );

  STKUNIT_ASSERT( fC_nm2->name() == std::string("C_STKFS_NM2") );
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

  //Coverage of EntityDimension::name in FieldData.cpp
  {
    const stk::mesh::EntityDimension&  entity_dimension_tag = stk::mesh::EntityDimension::tag();
    // static const char * name();

    entity_dimension_tag.name();
  }
}


//----------------------------------------------------------------------
// Test field restrictions: the mapping of ( field , part ) -> dimensions

void UnitTestFieldImpl::testFieldRestriction()
{
  const char method[] = "UnitTestFieldImpl::testFieldRestriction" ;

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
  STKUNIT_ASSERT( f3->restrictions()[0] ==
                  FieldRestriction( 0 , pA.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[1] ==
                  FieldRestriction( 1 , pB.mesh_meta_data_ordinal() ) );
  STKUNIT_ASSERT( f3->restrictions()[2] ==
                  FieldRestriction( 2 , pC.mesh_meta_data_ordinal() ) );

  f3->m_impl.insert_restriction( method , 0 , pB , stride + 1 );

  STKUNIT_ASSERT_EQUAL( f3->max_size( 0 ) , stride[3] );

  //------------------------------
  // Check for error detection of bad stride:
  {
    unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
    STKUNIT_ASSERT_THROW(
      f3->m_impl.insert_restriction( method , 0 , pA , bad_stride ),
      std::runtime_error
    );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Check for error detection in re-declaring an incompatible
  // field restriction.
  {
    STKUNIT_ASSERT_THROW(
      f3->m_impl.insert_restriction( method , 0 , pA , stride + 1 ),
      std::runtime_error
    );
    STKUNIT_ASSERT( f3->restrictions().size() == 4 );
  }

  // Verify and clean out any redundant restructions:

  f2->m_impl.verify_and_clean_restrictions( method , partRepo.get_all_parts() );
  f3->m_impl.verify_and_clean_restrictions( method , partRepo.get_all_parts() );

  STKUNIT_ASSERT( f3->restrictions().size() == 4 );

  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.

  partRepo.declare_subset( pD, pA );
  f2->m_impl.insert_restriction( method , 0 , pA , stride );
  f2->m_impl.insert_restriction( method , 0 , pD , stride );

  STKUNIT_ASSERT( f2->restrictions().size() == 2 );

  f2->m_impl.verify_and_clean_restrictions( method , partRepo.get_all_parts() );

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
    f2->m_impl.insert_restriction( method , 0 , pB , stride + 1 );
    f2->m_impl.verify_and_clean_restrictions( method , partRepo.get_all_parts() );
    partRepo.declare_subset( pD, pB );
    STKUNIT_ASSERT_THROW(
      f2->m_impl.verify_and_clean_restrictions( method , partRepo.get_all_parts() ),
      std::runtime_error
    );
  }

  //Test to cover print function in FieldBaseImpl.cpp and FieldBase.cpp
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

    //test stride[i] / stride[i-1] section of else-if
    stk::mesh::print(std::cout, "Field f4", *f4);
  }

  //Further tests to cover print function in FieldBase.cpp
  {
    //test stride[i] % stride[i-1] section of else-if

    //Create a new field with MetaData m and two restrictions

    FieldBase * const f5 =
      field_repo.declare_field( std::string("F5"),
                                data_traits<int>(),
                                3         /* # ranks */ ,
                                dim_tags  /* dimension tags */ ,
                                2         /* # states */ ,
                                &m );

    unsigned stride2[8] ;
    stride2[0] = 10 ;
    for ( unsigned i = 1 ; i < 3 ; ++i ) {
      stride2[i] = stride[i-1];
    }
    for ( unsigned i = 3 ; i < 8 ; ++i ) {
      stride2[i] = 0;
    }
    //  STKUNIT_ASSERT_THROW(
    f5->m_impl.insert_restriction( method , 0 , pA, stride2 );
    //  std::runtime_error
    //  );

    // f5->m_impl.insert_restriction( method , 1 , pB ,stride2 );

    stk::mesh::print(std::cout, "Field f5", *f5);

  }

  //Coverage for error from print_restriction in FieldBaseImpl.cpp when there is no stride (!stride[i])
  //Call print_restriction from insert_restriction
  {
    unsigned arg_no_stride[2];

    arg_no_stride[0] = 1;
    arg_no_stride[1] = 0;

    STKUNIT_ASSERT_THROW(
      f2->m_impl.insert_restriction(method, 0, pA, arg_no_stride),
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

// Unit test the FieldRelation copy constructor:
void UnitTestFieldImpl::testFieldRelation()
{

  FieldRelation rA;
  FieldRelation rB(rA);

  rA = rB;

}

}
}

