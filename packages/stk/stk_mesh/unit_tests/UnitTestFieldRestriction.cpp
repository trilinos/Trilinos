/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Types.hpp>

namespace {

STKUNIT_UNIT_TEST( UnitTestFieldRestriction, defaultConstruct )
{
  stk::mesh::FieldRestriction fr;
  STKUNIT_EXPECT_EQ( fr.selector(), stk::mesh::Selector() );
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), 0 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 0 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, construct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  STKUNIT_EXPECT_EQ( fr.selector(), stk::mesh::Selector(part_a) );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, copyConstruct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  stk::mesh::FieldRestriction tmpfr(fr);
  STKUNIT_EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqual )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  stk::mesh::FieldRestriction fr(part_a);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }
  stk::mesh::FieldRestriction tmpfr(part_b);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    tmpfr.stride(i) = i+10;
  }

  tmpfr = fr;
  STKUNIT_EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLess )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_b);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_a);
    STKUNIT_EXPECT_EQ( frA == frB, true );
  }
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLessInvalid )
{
  stk::mesh::Selector emptySelector;
  stk::mesh::FieldRestriction frA(emptySelector);
  stk::mesh::FieldRestriction frB;
  STKUNIT_EXPECT_EQ( frA < frB, false );
  STKUNIT_EXPECT_EQ( frB < frA, false );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqualEqual_and_NotEqual )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.stride(0) = 25;
    stk::mesh::FieldRestriction frB(part_a);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.stride(0) = 3;
    stk::mesh::FieldRestriction frB(part_b);
    frB.stride(0) = 3;
    STKUNIT_EXPECT_EQ( frA == frB, false );
    STKUNIT_EXPECT_EQ( frA != frB, true );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    stk::mesh::FieldRestriction frB(part_a);
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA;
    stk::mesh::FieldRestriction frB;
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, not_equal_stride )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.stride(0) = 25;
    stk::mesh::FieldRestriction frB(part_a);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    for (stk::mesh::Ordinal i=0 ; i < stk::mesh::MaximumFieldDimension ; ++i ) {
      frA.stride(i) = i+1;
    }
    stk::mesh::FieldRestriction frB(part_a);
    for (stk::mesh::Ordinal i=0 ; i < stk::mesh::MaximumFieldDimension ; ++i ) {
      frB.stride(i) = i+1;
    }
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), false );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), false );
    frB.stride(stk::mesh::MaximumFieldDimension-1) = 1;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
}


} //namespace <anonymous>

