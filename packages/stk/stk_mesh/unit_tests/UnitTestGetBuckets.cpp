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

#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

#include <algorithm>

namespace {

using stk::mesh::Bucket;
using stk::mesh::Part;
using stk::mesh::Selector;

const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

void sort_and_compare_eq(std::vector<Bucket*>& results,
                         std::vector<Bucket*>& expected_results)
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

  // Check bucket membership correctness

  {
    const Bucket & bucket = fix.m_entity1->bucket();
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = fix.m_entity2->bucket();
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = fix.m_entity3->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = fix.m_entity4->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = fix.m_entity5->bucket();
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partA ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partB ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partC ) );
    STKUNIT_ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  // Check get_buckets correctness

  const std::vector<Bucket*> & node_buckets = fix.m_bulk_data.buckets(NODE_RANK);

  {
    std::vector<Bucket*> get_buckets, expected_results;

    // get buckets selected by partA
    Selector selector(fix.m_partA);
    stk::mesh::get_buckets(selector, node_buckets, get_buckets);

    expected_results.push_back(&(fix.m_entity1->bucket()));
    expected_results.push_back(&(fix.m_entity2->bucket()));

    sort_and_compare_eq(get_buckets, expected_results);
  }

  {
    std::vector<Bucket*> get_buckets, expected_results;

    // get buckets selected by partB
    Selector selector(fix.m_partB);
    stk::mesh::get_buckets(selector, node_buckets, get_buckets);

    expected_results.push_back(&(fix.m_entity2->bucket()));
    expected_results.push_back(&(fix.m_entity3->bucket()));

    sort_and_compare_eq(get_buckets, expected_results);
  }

  {
    std::vector<Bucket*> get_buckets, expected_results;

    // get buckets selected by partC
    Selector selector(fix.m_partC);
    stk::mesh::get_buckets(selector, node_buckets, get_buckets);

    expected_results.push_back(&(fix.m_entity3->bucket()));
    expected_results.push_back(&(fix.m_entity4->bucket()));

    sort_and_compare_eq(get_buckets, expected_results);
  }
}

} // namespace

