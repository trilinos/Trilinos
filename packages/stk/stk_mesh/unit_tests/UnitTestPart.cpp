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
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/base/FieldRelation.hpp>

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::impl::PartRepository;

namespace {

STKUNIT_UNIT_TEST(UnitTestPart, testUnit)
{
  const int spatial_dimension = 3;
  MetaData m(spatial_dimension);
  PartRepository partRepo(&m);
  PartRepository partRepo2(&m);
  PartRepository partRepo3(&m);
  PartRepository partRepo4(&m);
  Part & universal = *partRepo.universal_part();
  m.commit();

  STKUNIT_ASSERT(  universal.supersets().empty() );
  STKUNIT_ASSERT( 1u ==  partRepo.get_all_parts().size() );

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
  STKUNIT_ASSERT( NPARTS ==  partRepo.get_all_parts().size() );
  STKUNIT_ASSERT_EQUAL(  partRepo.get_all_parts()[0] , &  universal );

  for ( unsigned i = 1 ; i < NPARTS ; ++i ) {
    STKUNIT_ASSERT( parts[i]->subsets().empty() );
    STKUNIT_ASSERT( parts[i]->mesh_meta_data_ordinal() == i );
    STKUNIT_ASSERT( 1u == parts[i]->supersets().size() );
    STKUNIT_ASSERT( &  universal == parts[i]->supersets()[0] );
    STKUNIT_ASSERT_EQUAL( parts[i] ,  universal.subsets()[i-1] );
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

}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestPart, testPartVector)
{
  const int spatial_dimension = 3;
  MetaData m(spatial_dimension);
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

} // empty namespace
