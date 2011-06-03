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
  stk::mesh::FieldRestriction fr;
  //STKUNIT_EXPECT_EQ( fr.entity_rank(), stk::mesh::InvalidEntityRank );
  STKUNIT_EXPECT_EQ( fr.part_ordinal(), stk::mesh::InvalidPartOrdinal );
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), 0 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 0 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, construct )
{
  stk::mesh::FieldRestriction fr(1,2);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  STKUNIT_EXPECT_EQ( fr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( fr.part_ordinal(), 2u );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( fr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( fr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, copyConstruct )
{
  stk::mesh::FieldRestriction fr(1,2);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }

  stk::mesh::FieldRestriction tmpfr(fr);
  STKUNIT_EXPECT_EQ( tmpfr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( tmpfr.part_ordinal(), 2u );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqual )
{
  stk::mesh::FieldRestriction fr(1,2);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    fr.stride(i) = i+1;
  }
  stk::mesh::FieldRestriction tmpfr(3,4);
  for (stk::mesh::Ordinal i = 0 ; i < stk::mesh::MaximumFieldDimension ; ++i) {
    tmpfr.stride(i) = i+10;
  }

  tmpfr = fr;
  STKUNIT_EXPECT_EQ( tmpfr.entity_rank(), 1u );
  STKUNIT_EXPECT_EQ( tmpfr.part_ordinal(), 2u );
  const int max_field_dimension = stk::mesh::MaximumFieldDimension;
  for (int i = 0 ; i < max_field_dimension ; ++i) {
    STKUNIT_EXPECT_EQ( tmpfr.stride(i), i+1 );
  }
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLess )
{
  {
    stk::mesh::FieldRestriction frA(1,1);
    stk::mesh::FieldRestriction frB(1,2);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
  {
    stk::mesh::FieldRestriction frA(1,1);
    stk::mesh::FieldRestriction frB(2,1);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
  {
    stk::mesh::FieldRestriction frA(0,1);
    stk::mesh::FieldRestriction frB(1,2);
    STKUNIT_EXPECT_EQ( frA < frB, true );
    STKUNIT_EXPECT_EQ( frB < frA, false );
  }
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorLessInvalid )
{
  stk::mesh::FieldRestriction frA(1,2);
  stk::mesh::FieldRestriction frB;
  STKUNIT_EXPECT_EQ( frA < frB, true );
  STKUNIT_EXPECT_EQ( frB < frA, false );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqualEqual_and_NotEqual )
{
  {
    stk::mesh::FieldRestriction frA(1,2);
    frA.stride(0) = 25;
    stk::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA(1,1);
    frA.stride(0) = 3;
    stk::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 3;
    STKUNIT_EXPECT_EQ( frA == frB, false );
    STKUNIT_EXPECT_EQ( frA != frB, true );
  }
  {
    stk::mesh::FieldRestriction frA(1,2);
    stk::mesh::FieldRestriction frB(2,2);
    STKUNIT_EXPECT_EQ( frA == frB, false );
    STKUNIT_EXPECT_EQ( frA != frB, true );
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
  {
    stk::mesh::FieldRestriction frA(1,2);
    frA.stride(0) = 25;
    stk::mesh::FieldRestriction frB(1,2);
    frB.stride(0) = 10;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
  {
    stk::mesh::FieldRestriction frA(1,2);
    for (stk::mesh::Ordinal i=0 ; i < stk::mesh::MaximumFieldDimension ; ++i ) {
      frA.stride(0) = i+1;
    }
    stk::mesh::FieldRestriction frB(1,2);
    for (stk::mesh::Ordinal i=0 ; i < stk::mesh::MaximumFieldDimension ; ++i ) {
      frB.stride(0) = i+1;
    }
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), false );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), false );
    frB.stride(stk::mesh::MaximumFieldDimension-1) = 1;
    STKUNIT_EXPECT_EQ( frA.not_equal_stride(frB), true );
    STKUNIT_EXPECT_EQ( frB.not_equal_stride(frA), true );
  }
}


} //namespace <anonymous>

