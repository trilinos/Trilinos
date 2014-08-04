/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/FieldRestriction.hpp>  // for FieldRestriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
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

