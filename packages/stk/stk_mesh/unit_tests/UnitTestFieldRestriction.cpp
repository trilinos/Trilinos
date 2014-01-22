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
  STKUNIT_EXPECT_EQ( 0, fr.num_scalars_per_entity() );
  STKUNIT_EXPECT_EQ( fr.dimension(), 0 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, construct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);
  STKUNIT_EXPECT_EQ( numScalarsPerEntity, fr.num_scalars_per_entity() );

  STKUNIT_EXPECT_EQ( fr.selector(), stk::mesh::Selector(part_a) );
  STKUNIT_EXPECT_EQ( fr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, copyConstruct )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");

  stk::mesh::FieldRestriction fr(part_a);
  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);

  stk::mesh::FieldRestriction tmpfr(fr);
  STKUNIT_EXPECT_EQ( numScalarsPerEntity, tmpfr.num_scalars_per_entity() );
  STKUNIT_EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
  STKUNIT_EXPECT_EQ( tmpfr.dimension(), 1 );
}


STKUNIT_UNIT_TEST( UnitTestFieldRestriction, operatorEqual )
{
  stk::mesh::MetaData meta(3);
  stk::mesh::Part& part_a = meta.declare_part("a");
  stk::mesh::Part& part_b = meta.declare_part("b");

  stk::mesh::FieldRestriction fr(part_a);
  stk::mesh::FieldRestriction tmpfr(part_b);

  const int numScalarsPerEntity = 2;
  fr.set_num_scalars_per_entity(numScalarsPerEntity);
  tmpfr.set_num_scalars_per_entity(numScalarsPerEntity+10);

  tmpfr = fr;
  STKUNIT_EXPECT_EQ( numScalarsPerEntity, tmpfr.num_scalars_per_entity() );
  STKUNIT_EXPECT_EQ( tmpfr.selector(), stk::mesh::Selector(part_a) );
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
    frA.set_num_scalars_per_entity(25);
    stk::mesh::FieldRestriction frB(part_a);
    frB.set_num_scalars_per_entity(10);
    STKUNIT_EXPECT_EQ( frA == frB, true );
    STKUNIT_EXPECT_EQ( frA != frB, false );
  }
  {
    stk::mesh::FieldRestriction frA(part_a);
    frA.set_num_scalars_per_entity(3);
    stk::mesh::FieldRestriction frB(part_b);
    frB.set_num_scalars_per_entity(3);
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

} //namespace <anonymous>

