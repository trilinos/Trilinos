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

#include <stddef.h>                     // for NULL
#include <iostream>                     // for ostream, operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateOld, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"  // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {

using namespace stk::mesh;

class UnitTestFieldImpl  : public ::testing::Test
{
public:
  UnitTestFieldImpl()
  : dummy_names(4, "dummy"),
    meta_data(3 /*dim*/, dummy_names),
    pA(meta_data.declare_part( std::string("A") , stk::topology::NODE_RANK)),
    pB(meta_data.declare_part( std::string("B") , stk::topology::NODE_RANK)),
    pC(meta_data.declare_part( std::string("C") , stk::topology::NODE_RANK)),
    pD(meta_data.declare_part( std::string("D") , stk::topology::NODE_RANK))
  {
  }

  void testFieldRestriction();
  void testRestrictions3parts();
  void testRestrictions3parts_union();
  void testRestrictions3parts_differentNumScalars();
  void testRestrictions3parts_supersetSubset();
  void testRestrictions4parts_1union_1intersection();
  void testRestrictions4parts_2unions();

  void print_restrictions(const FieldRestrictionVector& restrs)
  {
    for(const FieldRestriction& restr : restrs) {
      std::cout << restr.selector() << " : " << restr.num_scalars_per_entity() << std::endl;
    }
  }

private:
  std::vector<std::string> dummy_names;
  MetaData meta_data;

  Part& pA;
  Part& pB;
  Part& pC;
  Part& pD;
};

void UnitTestFieldImpl::testFieldRestriction()
{
  unsigned stride[8] ;

  stride[0] = 10 ;
  for ( unsigned i = 1 ; i < 8 ; ++i ) {
    stride[i] = ( i + 1 ) * stride[i-1] ;
  }

  const FieldVector  & allocated_fields = meta_data.get_fields();

  //------------------------------

  FieldBase * const f2 = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("F2"), 1/* # states */ );

  //------------------------------

  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);
  FieldBase * const edgeField = &meta_data.declare_field<double>( stk::topology::EDGE_RANK, std::string("edgeField"), 2/* #states*/);
  FieldBase * const faceField = &meta_data.declare_field<double>( stk::topology::FACE_RANK, std::string("faceField"), 2/* #states*/);

  FieldBase * const f3_old = nodeField->field_state( StateOld ) ;

  //------------------------------
  // Test for correctness of vector of declared fields.
  ASSERT_EQ(7u,  allocated_fields.size());
  ASSERT_TRUE( f2 == allocated_fields[0] );
  ASSERT_TRUE( nodeField == allocated_fields[1] );

  //------------------------------
  // Test for correctness of field internal state access:

  ASSERT_TRUE( f2     == f2->field_state( StateNone ) );
  ASSERT_TRUE( NULL   == f2->field_state( StateOld ) );
  ASSERT_TRUE( nodeField     == nodeField->field_state( StateNew ) );
  ASSERT_TRUE( f3_old == nodeField->field_state( StateOld ) );
  ASSERT_TRUE( NULL   == nodeField->field_state( StateNM1 ) );
  ASSERT_TRUE( nodeField     == f3_old->field_state( StateNew ) );
  ASSERT_TRUE( f3_old == f3_old->field_state( StateOld ) );
  ASSERT_TRUE( NULL   == f3_old->field_state( StateNM1 ) );

  // Declare three restrictions:

  meta_data.declare_field_restriction(*nodeField, pA, stride[0], stride[0]);
  meta_data.declare_field_restriction(*edgeField, pB, stride[1], stride[1]);
  meta_data.declare_field_restriction(*faceField, pC, stride[2], stride[2]);

  // Check for correctness of restrictions:

  ASSERT_TRUE( nodeField->restrictions().size() == 1 );
  ASSERT_TRUE( nodeField->restrictions()[0] == FieldRestriction( pA ) );
  ASSERT_TRUE( edgeField->restrictions()[0] == FieldRestriction( pB ) );
  ASSERT_TRUE( faceField->restrictions()[0] == FieldRestriction( pC ) );

  meta_data.declare_field_restriction(*nodeField, pB, stride[1], stride[1]);

  ASSERT_EQ( nodeField->max_size() , 20u );

  //------------------------------
  // Check for error detection of bad stride:
  {
    unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
    Selector selA(pA);
    ASSERT_THROW(meta_data.declare_field_restriction(*nodeField , selA, 5*4*6 , bad_stride[0] ), std::runtime_error);
    ASSERT_THROW(meta_data.declare_field_restriction(*nodeField , pA, 5*4*6 , bad_stride[0] ), std::runtime_error);
    ASSERT_EQ(2u, nodeField->restrictions().size());
  }

  // Check for error detection in re-declaring an incompatible
  // field restriction.
  {
    Selector selA(pA);
    ASSERT_THROW(meta_data.declare_field_restriction(*nodeField, selA, stride[1], stride[1]), std::runtime_error);
    ASSERT_THROW(meta_data.declare_field_restriction(*nodeField, pA, stride[1], stride[1]), std::runtime_error);
    ASSERT_EQ(2u, nodeField->restrictions().size());
  }

  // Verify and clean out any redundant restructions:

  ASSERT_TRUE( nodeField->restrictions().size() == 2 );

#ifndef NDEBUG
  //The following checking/cleaning of restrictions can be expensive, so it
  //is only done in debug mode.
  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.
  Part & pD = meta_data.declare_part( std::string("D") , stk::topology::NODE_RANK );

  meta_data.declare_part_subset( pD, pA );
  meta_data.declare_field_restriction(*f2, pA, stride[0], stride[0]);
  meta_data.declare_field_restriction(*f2, pD, stride[0], stride[0]);

  unsigned expected = 1;
  ASSERT_TRUE( f2->restrictions().size() == expected );

  {
    const FieldBase::Restriction & rA = stk::mesh::find_restriction(*f2, stk::topology::NODE_RANK, pA );
    const FieldBase::Restriction & rD = stk::mesh::find_restriction(*f2, stk::topology::NODE_RANK, pD );
    ASSERT_TRUE( & rA == & rD );
    ASSERT_TRUE( rA.selector() == pD );
  }

  //------------------------------
  // Introduce a new restriction, then introduce a
  // subset-superset relationship that renders the new restriction
  // redundant and incompatible.
  // Check that the verify_and_clean_restrictions method detects
  // this error condition.
  {
    meta_data.declare_field_restriction(*f2, pB, stride[1], stride[1]);
    ASSERT_THROW(meta_data.declare_part_subset(pD, pB), std::runtime_error);
  }

  //Coverage for error from print_restriction in FieldBaseImpl.cpp when there is no stride (!stride[i])
  //Call print_restriction from insert_restriction
  {
    unsigned arg_no_stride[2];

    arg_no_stride[0] = 1;
    arg_no_stride[1] = 0;

    ASSERT_ANY_THROW(meta_data.declare_field_restriction(*f2, pA, 0, arg_no_stride[0]));
  }
#endif
}

void UnitTestFieldImpl::testRestrictions3parts()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  const unsigned numScalarsPerEntity = 3;

  meta_data.declare_field_restriction(*nodeField, pA, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pB, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pC, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 1;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

void UnitTestFieldImpl::testRestrictions3parts_union()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  const unsigned numScalarsPerEntity = 3;

  Selector abcUnion = pA | pB | pC;
  meta_data.declare_field_restriction(*nodeField, abcUnion, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pC, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 1;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

void UnitTestFieldImpl::testRestrictions3parts_differentNumScalars()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  unsigned numScalarsPerEntity = 3;

  meta_data.declare_field_restriction(*nodeField, pA, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pB, numScalarsPerEntity, 1);

  numScalarsPerEntity = 1;
  meta_data.declare_field_restriction(*nodeField, pC, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 2;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

void UnitTestFieldImpl::testRestrictions3parts_supersetSubset()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  meta_data.declare_part_subset(pB, pC);

  const unsigned numScalarsPerEntity = 3;

  meta_data.declare_field_restriction(*nodeField, pA, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pB, numScalarsPerEntity, 1);
  meta_data.declare_field_restriction(*nodeField, pC, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 1;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

void UnitTestFieldImpl::testRestrictions4parts_1union_1intersection()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  const unsigned numScalarsPerEntity = 3;

  Selector unionAB = pA | pB;
  meta_data.declare_field_restriction(*nodeField, unionAB, numScalarsPerEntity, 1);

  Selector intersectionCD = pC & pD;
  meta_data.declare_field_restriction(*nodeField, intersectionCD, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 1;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

void UnitTestFieldImpl::testRestrictions4parts_2unions()
{
  FieldBase * const nodeField = &meta_data.declare_field<double>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);

  const unsigned numScalarsPerEntity = 3;

  Selector unionAB = pA | pB;
  meta_data.declare_field_restriction(*nodeField, unionAB, numScalarsPerEntity, 1);

  Selector unionCD = pC | pD;
  meta_data.declare_field_restriction(*nodeField, unionCD, numScalarsPerEntity, 1);

  const unsigned expectedNumRestrictions = 1;
  EXPECT_EQ(expectedNumRestrictions, nodeField->restrictions().size());
  print_restrictions(nodeField->restrictions());
}

TEST_F(UnitTestFieldImpl, testRestriction)
{
  testFieldRestriction();
}

TEST_F(UnitTestFieldImpl, testRestrictions3parts)
{
  testRestrictions3parts();
}

TEST_F(UnitTestFieldImpl, testRestrictions3parts_union)
{
  testRestrictions3parts_union();
}

TEST_F(UnitTestFieldImpl, testRestrictions3parts_differentNumScalars)
{
  testRestrictions3parts_differentNumScalars();
}

TEST_F(UnitTestFieldImpl, testRestrictions3parts_supersetSubset)
{
  testRestrictions3parts_supersetSubset();
}

TEST_F(UnitTestFieldImpl, testRestrictions4parts_1union_1intersection)
{
  testRestrictions4parts_1union_1intersection();
}

TEST_F(UnitTestFieldImpl, testRestrictions4parts_2unions)
{
  testRestrictions4parts_2unions();
}

}//namespace <anonymous>

