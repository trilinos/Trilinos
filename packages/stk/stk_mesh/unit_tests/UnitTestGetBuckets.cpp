/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>


#include <stk_mesh/fixtures/SelectorFixture.hpp>

namespace {

  using stk::mesh::fixtures::SelectorFixture ;

STKUNIT_UNIT_TEST( UnitTestGetBuckets, ExampleFixture )
{
  SelectorFixture fix ;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    STKUNIT_ASSERT_TRUE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

}

} // namespace

