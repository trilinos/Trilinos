// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector
#include "Shards_Array.hpp"
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateOld, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"  // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc






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

TEST(UnitTestFieldRestriction, testUnit)
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
    &meta_data.declare_field<VectorField>( stk::topology::NODE_RANK, std::string("F2"), 1/* # states */ );

  //------------------------------

//  FieldBase * const f3 = &meta_data.declare_field<VectorField>( std::string("F3"), 2/* #states*/);
  FieldBase * const nodeField = &meta_data.declare_field<VectorField>( stk::topology::NODE_RANK, std::string("nodeField"), 2/* #states*/);
  FieldBase * const edgeField = &meta_data.declare_field<VectorField>( stk::topology::EDGE_RANK, std::string("edgeField"), 2/* #states*/);
  FieldBase * const faceField = &meta_data.declare_field<VectorField>( stk::topology::FACE_RANK, std::string("faceField"), 2/* #states*/);

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

  //------------------------------
  // Declare some parts for restrictions:

  Part & pA = meta_data.declare_part( std::string("A") , stk::topology::NODE_RANK );
  Part & pB = meta_data.declare_part( std::string("B") , stk::topology::NODE_RANK );
  Part & pC = meta_data.declare_part( std::string("C") , stk::topology::NODE_RANK );
  Part & pD = meta_data.declare_part( std::string("D") , stk::topology::NODE_RANK );

  // Declare three restrictions:

  meta_data.declare_field_restriction(*nodeField, pA , stride[nodeField->field_array_rank()-1], stride[0] );
  meta_data.declare_field_restriction(*edgeField, pB , stride[edgeField->field_array_rank()], stride[1] );
  meta_data.declare_field_restriction(*faceField, pC , stride[faceField->field_array_rank()+1], stride[2] );

  // Check for correctness of restrictions:

  ASSERT_TRUE( nodeField->restrictions().size() == 1 );
  ASSERT_TRUE( nodeField->restrictions()[0] ==
                  FieldRestriction( pA ) );
  ASSERT_TRUE( edgeField->restrictions()[0] ==
                  FieldRestriction( pB ) );
  ASSERT_TRUE( faceField->restrictions()[0] ==
                  FieldRestriction( pC ) );

  meta_data.declare_field_restriction(*nodeField , pB , stride[nodeField->field_array_rank()], stride[1] );

  ASSERT_EQ( nodeField->max_size( stk::topology::NODE_RANK ) , 20u );

  //------------------------------
  // Check for error detection of bad stride:
  {
    unsigned bad_stride[4] = { 5 , 4 , 6 , 3 };
    ASSERT_THROW(
      meta_data.declare_field_restriction(*nodeField , pA, 5*4*6 , bad_stride[0] ),
      std::runtime_error
    );
    ASSERT_EQ(2u, nodeField->restrictions().size());
  }

  // Check for error detection in re-declaring an incompatible
  // field restriction.
  {
    ASSERT_THROW(
      meta_data.declare_field_restriction(*nodeField , pA , stride[nodeField->field_array_rank()], stride[1] ),
      std::runtime_error
    );
    ASSERT_EQ(2u, nodeField->restrictions().size());
  }

  // Verify and clean out any redundant restructions:

  ASSERT_TRUE( nodeField->restrictions().size() == 2 );

  //------------------------------
  // Introduce a redundant restriction, clean it, and
  // check that it was cleaned.

  std::cout<<"pA ord: "<<pA.mesh_meta_data_ordinal()<<", pD ord: "<<pD.mesh_meta_data_ordinal()<<std::endl;
  meta_data.declare_part_subset( pD, pA );
  meta_data.declare_field_restriction(*f2 , pA , stride[f2->field_array_rank()-1], stride[0] );
  meta_data.declare_field_restriction(*f2 , pD , stride[f2->field_array_rank()-1], stride[0] );

  ASSERT_TRUE( f2->restrictions().size() == 1 );

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
    meta_data.declare_field_restriction(*f2 , pB , stride[f2->field_array_rank()], stride[1] );
    ASSERT_THROW(
      meta_data.declare_part_subset( pD, pB ),
      std::runtime_error
    );
  }

  //Coverage for error from print_restriction in FieldBaseImpl.cpp when there is no stride (!stride[i])
  //Call print_restriction from insert_restriction
  {
    unsigned arg_no_stride[2];

    arg_no_stride[0] = 1;
    arg_no_stride[1] = 0;

    ASSERT_THROW(
      meta_data.declare_field_restriction(*f2, pA, 0, arg_no_stride[0]),
      std::runtime_error
    );
  }
}


}
}

