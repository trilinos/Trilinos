
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
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/PartRelation.hpp>

using stk::mesh::MetaData;
using stk::mesh::TopologicalMetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::impl::PartRepository;

class UnitTestPart {
public:
  UnitTestPart();
  ~UnitTestPart() {}

  const int spatial_dimension;
  MetaData m;
  MetaData m2;
  PartRepository partRepo;
  PartRepository partRepo2;
  PartRepository partRepo3;
  PartRepository partRepo4;
  PartRepository partRepo5;
  PartRepository partRepo6;
  Part & universal;
  PartVector intersection;
  PartVector intersection2;
  PartVector intersection3;
  PartVector intersection4;
  PartVector intersection5;
  PartRelation rA;

};

UnitTestPart::UnitTestPart()
  :spatial_dimension(3)
  ,m(TopologicalMetaData::entity_rank_names(spatial_dimension))
  ,m2(TopologicalMetaData::entity_rank_names(spatial_dimension))
  ,partRepo(&m)
  ,partRepo2(&m)
  ,partRepo3(&m)
  ,partRepo4(&m)
  ,partRepo5(&m)
  ,partRepo6(&m2)
  ,universal(*partRepo.universal_part())
  ,intersection()
  ,intersection2()
  ,intersection3()
  ,intersection4()
  ,intersection5()
  ,rA()
{
  m.commit();
}

namespace {

STKUNIT_UNIT_TEST(UnitTestPart, testUnit)
{
  UnitTestPart upart;


  STKUNIT_ASSERT( upart.universal.supersets().empty() );
  STKUNIT_ASSERT( 1u == upart.universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL( upart.universal.subsets()[0] , & upart.universal );

  //--------------------------------------------------------------------
  // Test multiple part creation

  enum { NPARTS = 100 };

  Part * parts[NPARTS];

  parts[0] = & upart.universal ;

  for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
    std::ostringstream name ;
    name << "Part_" << i ;
    parts[i] = upart.partRepo.declare_part( name.str() , 0 );
  }
  parts[99] = upart.partRepo.declare_part( "Part_99" , 1 );


  STKUNIT_ASSERT( upart.universal.supersets().empty() );
  STKUNIT_ASSERT( NPARTS == upart.universal.subsets().size() );
  STKUNIT_ASSERT_EQUAL( upart.universal.subsets()[0] , & upart.universal );

  for ( unsigned i = 1 ; i < NPARTS ; ++i ) {
    STKUNIT_ASSERT( parts[i]->subsets().empty() );
    STKUNIT_ASSERT( parts[i]->intersection_of().empty() );
    STKUNIT_ASSERT( parts[i]->mesh_meta_data_ordinal() == i );
    STKUNIT_ASSERT( 1u == parts[i]->supersets().size() );
    STKUNIT_ASSERT( & upart.universal == parts[i]->supersets()[0] );
    STKUNIT_ASSERT_EQUAL( parts[i] , upart.universal.subsets()[i] );
    STKUNIT_ASSERT_EQUAL( parts[i] , find(upart. universal.subsets() , parts[i]->name() ) );
  }

  //--------------------------------------------------------------------
  // Test multiple parts and transitive subset declarations:

  upart.partRepo.declare_subset( * parts[3], * parts[4] );
  upart.partRepo.declare_subset( * parts[4], * parts[5] );

  upart.partRepo.declare_subset( * parts[1], * parts[2] );
  // 1 and 2 pick up 4 and 5 via transitive relationship:
  upart.partRepo.declare_subset( * parts[2], * parts[3] );

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

  
  upart.intersection.push_back( parts[1] );
  upart.intersection.push_back( parts[2] );
  upart.intersection.push_back( parts[3] );
  upart.intersection.push_back( parts[4] ); // Smallest subset of 1..4

  // Test filtering of trivial intersection

  STKUNIT_ASSERT( parts[4] == upart.partRepo.declare_part( upart.intersection ) );

  // Test non-trivial intersection:

  upart.intersection.push_back( parts[6] );
  upart.intersection.push_back( parts[7] );

  Part & pint_4_6_7 = * upart.partRepo.declare_part( upart.intersection );

  STKUNIT_ASSERT( 3u == pint_4_6_7.intersection_of().size() );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[4] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[6] ) );
  STKUNIT_ASSERT( contain( pint_4_6_7.intersection_of() , * parts[7] ) );

  STKUNIT_ASSERT( 7u == pint_4_6_7.supersets().size() );

 // Test redeclaration of intersection, should give the same part back:

  STKUNIT_ASSERT( pint_4_6_7 == * upart.partRepo.declare_part( upart.intersection ) );

  upart.partRepo.declare_subset( pint_4_6_7, * parts[8] );

  //--------------------------------------------------------------------
  // Test intersection-induced subset relationship

  upart.partRepo.declare_subset( * parts[7], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  upart.partRepo.declare_subset( * parts[6], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  upart.partRepo.declare_subset( * parts[3], * parts[10] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[10] ) );

  upart.partRepo.declare_subset( * parts[4], * parts[10] );
  STKUNIT_ASSERT( contain( pint_4_6_7.subsets() , * parts[10] ) );

  // Test intersection-induced subset relationship triggered from a subset

  upart.partRepo.declare_subset( * parts[7], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );


  upart.partRepo.declare_subset( * parts[6], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  upart.partRepo.declare_subset( * parts[3], * parts[11] );
  STKUNIT_ASSERT( ! contain( pint_4_6_7.subsets() , * parts[11] ) );

  upart.partRepo.declare_subset( * parts[5], * parts[11] );
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

      upart.intersection4.push_back( &part_not_in_same_universe );

      PartVector::const_iterator i = upart.intersection4.begin() ;
      Part * const p = *i ;
      STKUNIT_ASSERT_THROW(
          upart.partRepo3.declare_subset(*parts[5], *p),
	  std::runtime_error
          );
    }

  //Test to cover assert_not_superset in PartRepository - parts[11] is not a superset of parts[7]
  {
      STKUNIT_ASSERT_THROW(
	  upart.partRepo.declare_subset(*parts[11], *parts[7] ),
          std::runtime_error
	  );
    }

  //Test to cover assert_not_superset in PartRepository - parts[99] is not the same rank of parts[11]
  {
      STKUNIT_ASSERT_THROW(
	  upart.partRepo.declare_subset(*parts[11], *parts[99] ),
          std::runtime_error
	  );
    }

  //Test to cover declare_part(arg_name, arg_rank) - Part_99 of rank 1 already exists..
  {
      STKUNIT_ASSERT_THROW(
	  upart.partRepo.declare_part("Part_99", 0 ),
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

      parts2[0] = & upart.universal ;

      for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
        std::ostringstream name ;
        name << "Part_" << i ;
        parts2[i] = upart.partRepo4.declare_part( name.str() , 0 );
      }

      upart.partRepo4.declare_subset( * parts2[3], * parts2[4] );
      upart.partRepo4.declare_subset( * parts2[4], * parts2[5] );

      upart.partRepo4.declare_subset( * parts2[1], * parts2[2] );
      // 1 and 2 pick up 4 and 5 via transitive relationship:
      upart.partRepo4.declare_subset( * parts2[2], * parts2[3] );


      upart.intersection5.push_back( parts2[1] );
      upart.intersection5.push_back( parts2[2] );
      upart.intersection5.push_back( parts2[3] );
      upart.intersection5.push_back( parts2[4] ); // Smallest subset of 1..4
      upart.intersection5.push_back( parts2[6] );
      upart.intersection5.push_back( parts2[7] );

      Part & partB = * upart.partRepo2.declare_part("{Part_4^Part_6^Part_7}", 0);

      std::cout << "UnitTestPart name of partB is " << partB.name()  << std::endl ;

      Part & pintersect_4_6_7 = * upart.partRepo2.declare_part(  upart.intersection5 );

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
  UnitTestPart upart;

  stk::mesh::Part &part_not_a_subset = upart.m2.declare_part ( "part_not_a_subset");
  upart.m2.commit();

  upart.intersection3.push_back( &part_not_a_subset );

  STKUNIT_ASSERT_THROW(
      upart.partRepo6.declare_part(upart.intersection3), 
      std::runtime_error
      );

}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestPart, testPartVector)
{
  UnitTestPart upart;

  Part * const pa = upart.partRepo5.declare_part( std::string("a") , 0 );
  Part * const pb = upart.partRepo5.declare_part( std::string("b") , 0 );
  Part * const pc = upart.partRepo5.declare_part( std::string("c") , 0 );
  Part * const pd = upart.partRepo5.declare_part( std::string("d") , 0 );
  Part * const pe = upart.partRepo5.declare_part( std::string("e") , 0 );
  Part * const pf = upart.partRepo5.declare_part( std::string("f") , 0 );

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

  Part * const pabc = upart.partRepo5.declare_part( std::string("abc") , 0 );
  Part * const pbcd = upart.partRepo5.declare_part( std::string("bcd") , 0 );
  Part * const pdef = upart.partRepo5.declare_part( std::string("def") , 0 );

  upart.partRepo5.declare_subset( * pabc, *pa );
  upart.partRepo5.declare_subset( * pabc, *pb );
  upart.partRepo5.declare_subset( * pabc, *pc );

  upart.partRepo5.declare_subset( * pbcd, *pb );
  upart.partRepo5.declare_subset( * pbcd, *pc );
  upart.partRepo5.declare_subset( * pbcd, *pd );

  upart.partRepo5.declare_subset( * pdef, *pd );
  upart.partRepo5.declare_subset( * pdef, *pe );
  upart.partRepo5.declare_subset( * pdef, *pf );

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
  UnitTestPart upart;
  PartRelation rB(upart.rA);

  upart.rA = rB;

}


}
