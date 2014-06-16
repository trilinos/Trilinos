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

#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Types.hpp>

namespace {


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, defaultConstruct )
{
  stk_classic::mesh::FieldRestriction fr;
  STKUNIT_EXPECT_EQ( fr.part_ordinal(), stk_classic::mesh::InvalidPartOrdinal );
  for (stk_classic::mesh::Ordinal i = 0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), 0 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 0 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, construct )
{
  stk_classic::mesh::FieldRestriction fr(1,2);
  for (stk_classic::mesh::Ordinal i = 0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  STKUNIT_EXPECT_EQ( fr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( fr.part_ordinal(), 2u );
  const int max_field_dimension = stk_classic::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, copyConstruct )
{
  stk_classic::mesh::FieldRestriction fr(1,2);
  for (stk_classic::mesh::Ordinal i = 0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  stk_classic::mesh::FieldRestriction tmpfr(fr);
  STKUNIT_EXPECT_EQ( tmpfr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( tmpfr.part_ordinal(), 2u );
  const int max_field_dimension = stk_classic::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqual )
{
  stk_classic::mesh::FieldRestriction fr(1,2);
  for (stk_classic::mesh::Ordinal i = 0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }
  stk_classic::mesh::FieldRestriction tmpfr(3,4);
  for (stk_classic::mesh::Ordinal i = 0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i) {
    tmpfr.stride(i) = i+10;
  }

  tmpfr = fr;
  STKUNIT_EXPECT_EQ( tmpfr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( tmpfr.part_ordinal(), 2u );
  const int max_field_dimension = stk_classic::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLess )
{
  {
    stk_classic::mesh::FieldRestriction frA(1,1);
    stk_classic::mesh::FieldRestriction frB(1,2);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
  {
    stk_classic::mesh::FieldRestriction frA(1,1);
    stk_classic::mesh::FieldRestriction frB(2,1);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
  {
    stk_classic::mesh::FieldRestriction frA(0,1);
    stk_classic::mesh::FieldRestriction frB(1,2);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLessInvalid )
{
  stk_classic::mesh::FieldRestriction frA(1,2);
  stk_classic::mesh::FieldRestriction frB;
  STKUNIT_EXPECT_EQ( frA < frB, true );
  STKUNIT_EXPECT_EQ( frB < frA, false );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqualEqual_and_NotEqual )
{
  {
    stk_classic::mesh::FieldRestriction frA(1,2);
    frA.stride(0) = 25;
    stk_classic::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
  {
    stk_classic::mesh::FieldRestriction frA(1,1);
    frA.stride(0) = 3;
    stk_classic::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 3;
    STKUNIT_EXPECT_EQ( frA == frB, false );
    STKUNIT_EXPECT_EQ( frA != frB, true );
  }
  {
    stk_classic::mesh::FieldRestriction frA(1,2);
    stk_classic::mesh::FieldRestriction frB(2,2);
    STKUNIT_EXPECT_EQ( frA == frB, false );
    STKUNIT_EXPECT_EQ( frA != frB, true );
  }
  {
    stk_classic::mesh::FieldRestriction frA;
    stk_classic::mesh::FieldRestriction frB;
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, not_equal_stride )
{
  {
    stk_classic::mesh::FieldRestriction frA(1,2);
    frA.stride(0) = 25;
    stk_classic::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
  {
    stk_classic::mesh::FieldRestriction frA(1,2);
    for (stk_classic::mesh::Ordinal i=0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i ) {
      frA.stride(0) = i+1;
    }
    stk_classic::mesh::FieldRestriction frB(1,2);
    for (stk_classic::mesh::Ordinal i=0 ; i < stk_classic::mesh::MaximumFieldDimension ; ++i ) {
      frB.stride(0) = i+1;
    }
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), false );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), false );
    frB.stride(stk_classic::mesh::MaximumFieldDimension-1) = 1;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
}


} //namespace <anonymous>

