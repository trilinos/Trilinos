/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>                    // for sort
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/fixtures/SelectorFixture.hpp>  // for SelectorFixture
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }






namespace {

using stk::mesh::Bucket;
using stk::mesh::Part;
using stk::mesh::Selector;
using stk::mesh::BulkData;
using stk::mesh::Entity;

// Just make copies so we don't have to worry about constness
template <class T>
void sort_and_compare_eq(std::vector<T*> results,
                         std::vector<T*> expected_results)
{
  std::sort(expected_results.begin(), expected_results.end());
  std::sort(results.begin(), results.end());
  STKUNIT_ASSERT_EQ(results.size(), expected_results.size());
  for (unsigned i = 0; i < results.size(); ++i) {
    STKUNIT_ASSERT(results[i] == expected_results[i]);
  }
}

STKUNIT_UNIT_TEST( UnitTestGetBuckets, ExampleFixture )
{
  // Using the SelectorFixture, test for correct bucket membership and
  // correct results from get_buckets.

  // Generate mesh

  stk::mesh::fixtures::SelectorFixture fix ;
  fix.m_meta_data.commit();

  fix.m_bulk_data.modification_begin();
  fix.generate_mesh();
  STKUNIT_ASSERT(fix.m_bulk_data.modification_end());
  const stk::mesh::BulkData& mesh = fix.m_bulk_data;

  // Check bucket membership correctness

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity1);
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity2);
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity3);
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity4);
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity5);
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

}

} // empty namespace
