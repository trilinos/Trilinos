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

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/PartRelation.hpp>

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::impl::PartRepository;

namespace {

STKUNIT_UNIT_TEST(UnitTestPart, testUnit)
{
  const int spatial_dimension = 3;
  MetaData m(stk::mesh::fem::entity_rank_names(spatial_dimension));
  PartRepository partRepo(&m);
  PartRepository partRepo2(&m);
  PartRepository partRepo3(&m);
  PartRepository partRepo4(&m);
  Part & universal = *partRepo.universal_part();
  PartVector intersection;
  PartVector intersection2;
  PartVector intersection4;
  PartVector intersection5;
  m.commit();

  STKUNIT_ASSERT(  universal.supersets().empty() );
  STKUNIT_ASSERT( 1u ==  universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL(  universal.subsets()[0] , &  universal );

  //--------------------------------------------------------------------
  // Test multiple part creation

  enum { NPARTS = 100 };

  Part * parts[NPARTS];

  parts[0] = &  universal ;

  for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
    std::ostringstream name ;
    name << "Part_" << i ;
    parts[i] =  partRepo.declare_part( name.str() , 0 );
  }
  parts[99] =  partRepo.declare_part( "Part_99" , 1 );

  STKUNIT_ASSERT(  universal.supersets().empty() );
  STKUNIT_ASSERT( NPARTS ==  universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL(  universal.subsets()[0] , &  universal );

  for ( unsigned i = 1 ; i < NPARTS ; ++i ) {
    STKUNIT_ASSERT( parts[i]->subsets().empty() );
    STKUNIT_ASSERT( parts[i]->intersection_of().empty() );
    STKUNIT_ASSERT( parts[i]->mesh_meta_data_ordinal() == i );
    STKUNIT_ASSERT( 1u == parts[i]->supersets().size() );
    STKUNIT_ASSERT( &  universal == parts[i]->supersets()[0] );
    STKUNIT_ASSERT_EQUAL( parts[i] ,  universal.subsets()[i] );
    STKUNIT_ASSERT_EQUAL( parts[i] , find(  universal.subsets() , parts[i]->name() ) );
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

  intersection.push_back( parts[1] );
  intersection.push_back( parts[2] );
  intersection.push_back( parts[3] );
  intersection.push_back( parts[4] ); // Smallest subset of 1..4

  // Test filtering of trivial intersection

  STKUNIT_ASSERT( parts[4] ==  partRepo.declare_part(  intersection ) );

  // Test non-trivial intersection:

  intersection.push_back( parts[6] );
  intersection.push_back( parts[7] );

  Part & pint_4_6_7 = *  partRepo.declare_part(  intersection );

  STKUNIT_ASSERT( 3u == pint_4_6_7.intersection_of().size() );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[4] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[6] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[7] ) );

  STKUNIT_ASSERT( 7u == pint_4_6_7.supersets().size() );

  // Test redeclaration of intersection, should give the same part back:

  STKUNIT_ASSERT( pint_4_6_7 == *  partRepo.declare_part(  intersection ) );

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

  //Test to cover assert_same_universe in PartRepository - Part is not in the same universe
  {
      std::vector< std::string > dummy_names(1);
      dummy_names[0].assign("dummy");

      stk::mesh::MetaData meta2( dummy_names );

      stk::mesh::Part &part_not_in_same_universe = meta2.declare_part ( "part_not_in_same_universe");
      meta2.commit();

       intersection4.push_back( &part_not_in_same_universe );

      PartVector::const_iterator i =  intersection4.begin() ;
      Part * const p = *i ;
      STKUNIT_ASSERT_THROW(
           partRepo3.declare_subset(*parts[5], *p),
           std::runtime_error
           );
  }

  //Test to cover assert_not_superset in PartRepository - parts[11] is not a superset of parts[7]
  {
      STKUNIT_ASSERT_THROW(
	   partRepo.declare_subset(*parts[11], *parts[7] ),
          std::runtime_error
	  );
  }

  //Test to cover assert_not_superset in PartRepository - parts[99] is not the same rank of parts[11]
  {
      STKUNIT_ASSERT_THROW(
	   partRepo.declare_subset(*parts[11], *parts[99] ),
          std::runtime_error
	  );
  }

  //Test to cover declare_part(arg_name, arg_rank) - Part_99 of rank 1 already exists..
  {
      STKUNIT_ASSERT_THROW(
	   partRepo.declare_part("Part_99", 0 ),
          std::runtime_error
	  );
  }

  //Test to cover declare_part(const PartVector & part_intersect) in PartRepository - failed from malicious abuse
  {
    int ok = 0 ;
    try {
      //create part with intersection

      // Test declaration of an intersection part

      enum { NPARTS = 100 };

      Part * parts2[ NPARTS ] ;

      parts2[0] = &  universal ;

      for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
        std::ostringstream name ;
        name << "Part_" << i ;
        parts2[i] =  partRepo4.declare_part( name.str() , 0 );
      }

       partRepo4.declare_subset( * parts2[3], * parts2[4] );
       partRepo4.declare_subset( * parts2[4], * parts2[5] );

       partRepo4.declare_subset( * parts2[1], * parts2[2] );
      // 1 and 2 pick up 4 and 5 via transitive relationship:
       partRepo4.declare_subset( * parts2[2], * parts2[3] );

       intersection5.push_back( parts2[1] );
       intersection5.push_back( parts2[2] );
       intersection5.push_back( parts2[3] );
       intersection5.push_back( parts2[4] ); // Smallest subset of 1..4
       intersection5.push_back( parts2[6] );
       intersection5.push_back( parts2[7] );

      Part & partB = *  partRepo2.declare_part("{Part_4^Part_6^Part_7}", 0);

      std::cout << "UnitTestPart name of partB is " << partB.name()  << std::endl ;

      Part & pintersect_4_6_7 = *  partRepo2.declare_part(   intersection5 );

      std::cout << "UnitTestPart name of intersection part, pintersect_4_6_7 is " << pintersect_4_6_7.name()  << std::endl ;
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestPart CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }

    if ( ! ok ) {
          throw std::runtime_error("UnitTestPart FAILED to catch error for declare_part in PartRepository");
    }
  }
}

STKUNIT_UNIT_TEST(UnitTestPart, testPartNotaSubset)
{
  //Test to cover assert_contain in PartRepository - Part is not a subset
  const int spatial_dimension = 3;
  MetaData m(stk::mesh::fem::entity_rank_names(spatial_dimension));
  PartVector intersection;
  PartRepository partRepo(&m);

  stk::mesh::Part &part_not_a_subset =  m.declare_part ( "part_not_a_subset");
  m.commit();

  intersection.push_back( &part_not_a_subset );

  STKUNIT_ASSERT_THROW(
      partRepo.declare_part(intersection),
      std::runtime_error
      );
}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestPart, testPartVector)
{
  const int spatial_dimension = 3;
  MetaData m(stk::mesh::fem::entity_rank_names(spatial_dimension));
  PartRepository partRepo(&m);

  Part * const pa =  partRepo.declare_part( std::string("a") , 0 );
  Part * const pb =  partRepo.declare_part( std::string("b") , 0 );
  Part * const pc =  partRepo.declare_part( std::string("c") , 0 );
  Part * const pd =  partRepo.declare_part( std::string("d") , 0 );
  Part * const pe =  partRepo.declare_part( std::string("e") , 0 );
  Part * const pf =  partRepo.declare_part( std::string("f") , 0 );

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

  Part * const pabc =  partRepo.declare_part( std::string("abc") , 0 );
  Part * const pbcd =  partRepo.declare_part( std::string("bcd") , 0 );
  Part * const pdef =  partRepo.declare_part( std::string("def") , 0 );

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

// Unit test the PartRelation copy constructor:
STKUNIT_UNIT_TEST(UnitTestPart, testPartRelation)
{
  PartRelation rA;
  PartRelation rB( rA);

  rA = rB;
}

} // empty namespace
