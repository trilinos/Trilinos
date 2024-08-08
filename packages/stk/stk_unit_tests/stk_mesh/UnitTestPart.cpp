// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stddef.h>                     // for size_t
#include <sstream>                      // for ostringstream, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Part.hpp>       // for PartLess, intersect, Part, etc
#include <stk_mesh/baseImpl/PartRepository.hpp>  // for PartRepository
#include <string>                       // for string, char_traits
#include <vector>                       // for vector
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/SortAndUnique.hpp"  // for sort_and_unique

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::impl::PartRepository;

namespace {

TEST(UnitTestPart, testUnit)
{
  const int spatial_dimension = 3;
  MetaData m(spatial_dimension);
  PartRepository partRepo(&m);
  PartRepository partRepo2(&m);
  PartRepository partRepo3(&m);
  PartRepository partRepo4(&m);
  Part & universal = *partRepo.universal_part();
  m.commit();

  ASSERT_TRUE(  universal.supersets().empty() );
  ASSERT_TRUE( 1u ==  partRepo.get_all_parts().size() );

  //--------------------------------------------------------------------
  // Test multiple part creation

  enum { NPARTS = 100 };

  Part * parts[NPARTS];

  parts[0] = &  universal ;

  for ( int i = 1 ; i < NPARTS-1 ; ++i ) {
    std::ostringstream name ;
    name << "Part_" << i ;
    parts[i] =  partRepo.declare_part( name.str() , stk::topology::NODE_RANK );
  }
  parts[99] =  partRepo.declare_part( "Part_99" , stk::topology::EDGE_RANK );

  ASSERT_TRUE(  universal.supersets().empty() );
  ASSERT_TRUE( NPARTS ==  partRepo.get_all_parts().size() );
  ASSERT_EQ(  partRepo.get_all_parts()[0] , &  universal );

  for ( unsigned i = 1 ; i < NPARTS ; ++i ) {
    ASSERT_TRUE( parts[i]->subsets().empty() );
    ASSERT_TRUE( parts[i]->mesh_meta_data_ordinal() == i );
    ASSERT_TRUE( 1u == parts[i]->supersets().size() );
    ASSERT_TRUE( &  universal == parts[i]->supersets()[0] );
    ASSERT_EQ( parts[i] ,  universal.subsets()[i-1] );
    ASSERT_EQ( parts[i] , find(  universal.subsets() , parts[i]->name() ) );
  }

  //--------------------------------------------------------------------
  // Test multiple parts and transitive subset declarations:

  partRepo.declare_subset( * parts[3], * parts[4] );
  partRepo.declare_subset( * parts[4], * parts[5] );

  partRepo.declare_subset( * parts[1], * parts[2] );
  // 1 and 2 pick up 4 and 5 via transitive relationship:
  partRepo.declare_subset( * parts[2], * parts[3] );

  ASSERT_TRUE( 4u == parts[1]->subsets().size() );
  ASSERT_TRUE( 3u == parts[2]->subsets().size() );
  ASSERT_TRUE( 2u == parts[3]->subsets().size() );
  ASSERT_TRUE( 1u == parts[4]->subsets().size() );
  ASSERT_TRUE( 0u == parts[5]->subsets().size() );

  ASSERT_TRUE( contain( parts[1]->subsets() , * parts[2] ) );
  ASSERT_TRUE( contain( parts[1]->subsets() , * parts[3] ) );
  ASSERT_TRUE( contain( parts[1]->subsets() , * parts[4] ) );
  ASSERT_TRUE( contain( parts[1]->subsets() , * parts[5] ) );

  ASSERT_TRUE( contain( parts[5]->supersets() , * parts[1] ) );
  ASSERT_TRUE( contain( parts[5]->supersets() , * parts[2] ) );
  ASSERT_TRUE( contain( parts[5]->supersets() , * parts[3] ) );
  ASSERT_TRUE( contain( parts[5]->supersets() , * parts[4] ) );

}

//----------------------------------------------------------------------

TEST(UnitTestPart, testPartVector)
{
  const int spatial_dimension = 3;
  MetaData m(spatial_dimension);
  PartRepository partRepo(&m);

  Part * const pa =  partRepo.declare_part( std::string("a") , stk::topology::NODE_RANK );
  Part * const pb =  partRepo.declare_part( std::string("b") , stk::topology::NODE_RANK );
  Part * const pc =  partRepo.declare_part( std::string("c") , stk::topology::NODE_RANK );
  Part * const pd =  partRepo.declare_part( std::string("d") , stk::topology::NODE_RANK );
  Part * const pe =  partRepo.declare_part( std::string("e") , stk::topology::NODE_RANK );
  Part * const pf =  partRepo.declare_part( std::string("f") , stk::topology::NODE_RANK );

  ASSERT_TRUE( ! intersect( *pa , *pb ) );
  ASSERT_TRUE( ! intersect( *pb , *pc ) );
  ASSERT_TRUE( ! intersect( *pc , *pd ) );
  ASSERT_TRUE( ! intersect( *pd , *pe ) );
  ASSERT_TRUE( ! intersect( *pe , *pf ) );
  ASSERT_TRUE( ! intersect( *pf , *pa ) );

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

  stk::util::sort_and_unique( vabc, stk::mesh::PartLess() );
  stk::util::sort_and_unique( vbcd, stk::mesh::PartLess() );
  stk::util::sort_and_unique( vdef, stk::mesh::PartLess() );

  vresult.clear();
  ASSERT_EQ( size_t(2) , intersect( vabc , vbcd ) );
  size_t intersect_size = intersect( vabc , vbcd , vresult );
  ASSERT_EQ( size_t(2) , intersect_size );
  ASSERT_EQ( pb , vresult[0] );
  ASSERT_EQ( pc , vresult[1] );

  vresult.clear();
  ASSERT_EQ( size_t(1) , intersect( vdef , vbcd ) );
  intersect_size = intersect( vdef , vbcd , vresult );
  ASSERT_EQ( size_t(1) , intersect_size );
  ASSERT_EQ( pd , vresult[0] );

  vresult.clear();
  ASSERT_EQ( size_t(0) , intersect( vdef , vabc ) );
  intersect_size = intersect( vdef , vabc , vresult );
  ASSERT_EQ( size_t(0) , intersect_size );
  ASSERT_EQ( size_t(0) , vresult.size() );

  Part * const pabc =  partRepo.declare_part( std::string("abc") , stk::topology::NODE_RANK );
  Part * const pbcd =  partRepo.declare_part( std::string("bcd") , stk::topology::NODE_RANK );
  Part * const pdef =  partRepo.declare_part( std::string("def") , stk::topology::NODE_RANK );

  partRepo.declare_subset( * pabc, *pa );
  partRepo.declare_subset( * pabc, *pb );
  partRepo.declare_subset( * pabc, *pc );

  partRepo.declare_subset( * pbcd, *pb );
  partRepo.declare_subset( * pbcd, *pc );
  partRepo.declare_subset( * pbcd, *pd );

  partRepo.declare_subset( * pdef, *pd );
  partRepo.declare_subset( * pdef, *pe );
  partRepo.declare_subset( * pdef, *pf );

  ASSERT_TRUE( intersect( *pabc , *pa ) );
  ASSERT_TRUE( intersect( *pabc , *pb ) );
  ASSERT_TRUE( intersect( *pabc , *pc ) );
  ASSERT_TRUE( intersect( *pa , *pabc ) );
  ASSERT_TRUE( intersect( *pb , *pabc ) );
  ASSERT_TRUE( intersect( *pc , *pabc ) );

  ASSERT_TRUE( intersect( *pabc , *pbcd ) );
  ASSERT_TRUE( ! intersect( *pabc , *pdef ) );
}

} // empty namespace
