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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/FieldRestriction.hpp>  // for FieldRestriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator<<
namespace stk { namespace mesh { class Part; } }

namespace {

TEST( UnitTestFieldRestriction, defaultConstruct )
{
  stk::mesh::FieldRestriction fr;
  EXPECT_EQ( fr.selector(), stk::mesh::Selector() );
  EXPECT_EQ( 0, fr.num_scalars_per_entity() );
  EXPECT_EQ( fr.dimension(), 0 );
}


TEST( UnitTestFieldRestriction, construct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);
  EXPECT_EQ( numScalarsPerEntity, fr.num_scalars_per_entity() );

  const int dim = 1;
  fr.set_dimension(dim);
  EXPECT_EQ( fr.selector(), stk::mesh::Selector(part_a) );
  EXPECT_EQ(dim, fr.dimension());
}


TEST( UnitTestFieldRestriction, copyConstruct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);
  const int dim = 1;
  fr.set_dimension(dim);

  stk::mesh::FieldRestriction tmpfr(fr);
  EXPECT_EQ( numScalarsPerEntity, tmpfr.num_scalars_per_entity() );
  EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
  EXPECT_EQ( dim, tmpfr.dimension() );
}

TEST( UnitTestFieldRestriction, selects_part)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  stk::mesh::FieldRestriction fr(part_a);
  EXPECT_TRUE(fr.selects(part_a));
  EXPECT_FALSE(fr.selects(part_b));
}

TEST( UnitTestFieldRestriction, union_selects_part)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  stk::mesh::FieldRestriction fr(part_a);
  fr.add_union(part_b);
  EXPECT_TRUE(fr.selects(part_a));
  EXPECT_TRUE(fr.selects(part_b));
}

TEST( UnitTestFieldRestriction, operatorEqual )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  stk::mesh::FieldRestriction fr(part_a);
  stk::mesh::FieldRestriction tmpfr(part_b);

  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);
  tmpfr.set_num_scalars_per_entity(numScalarsPerEntity+10);
  const int dim = 1;
  fr.set_dimension(dim);

  tmpfr = fr;
  EXPECT_EQ( numScalarsPerEntity, tmpfr.num_scalars_per_entity() );
  EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
  EXPECT_EQ( dim, tmpfr.dimension() );
}


TEST( UnitTestFieldRestriction, operatorLess )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_b);
    EXPECT_EQ( frA < frB, true );
    EXPECT_EQ( frB < frA, false );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_a);
    EXPECT_EQ( frA == frB, true );
  }
}


TEST( UnitTestFieldRestriction, operatorLessInvalid )
{
  stk::mesh::Selector emptySelector;
  stk::mesh::FieldRestriction frA(emptySelector);
  stk::mesh::FieldRestriction frB;
  EXPECT_EQ( frA < frB, false );
  EXPECT_EQ( frB < frA, false );
}


TEST( UnitTestFieldRestriction, operatorEqualEqual_and_NotEqual )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.set_num_scalars_per_entity(25);
    stk::mesh::FieldRestriction frB(part_a);
    frB.set_num_scalars_per_entity(10);
    EXPECT_EQ( frA == frB, true );
    EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.set_num_scalars_per_entity(3);
    stk::mesh::FieldRestriction frB(part_b);
    frB.set_num_scalars_per_entity(3);
    EXPECT_EQ( frA == frB, false );
    EXPECT_EQ( frA != frB, true );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_a);
    EXPECT_EQ( frA == frB, true );
    EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA;
    stk::mesh::FieldRestriction frB;
    EXPECT_EQ( frA == frB, true );
    EXPECT_EQ( frA != frB, false );
  }
}

} //namespace <anonymous>

