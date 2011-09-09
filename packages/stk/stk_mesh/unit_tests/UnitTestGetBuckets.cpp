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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

#include <algorithm>

#include <boost/foreach.hpp>

namespace {

using stk::mesh::Bucket;
using stk::mesh::Part;
using stk::mesh::Selector;
using stk::mesh::BulkData;
using stk::mesh::Entity;

const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

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

void check_selected_buckets(const Selector& selector,
                            const BulkData& mesh,
                            const std::vector<Bucket*>& node_buckets,
                            const std::vector<Entity*>& expected_selected_entities)
{
  //
  // Check buckets
  //

  std::vector<Bucket*> get_buckets, get_buckets_range, get_buckets_alt_range, expected_buckets;
  std::set<Bucket*> selected_bucket_set;

  // Compute vector of buckets that should be selected
  BOOST_FOREACH( Entity* entity, expected_selected_entities ) {
    selected_bucket_set.insert(&(entity->bucket()));
  }
  std::copy(selected_bucket_set.begin(), selected_bucket_set.end(), std::back_inserter(expected_buckets));

  // get buckets selected by selector
  stk::mesh::get_buckets(selector, node_buckets, get_buckets);

  // get buckets selected by selector using range API
  // all buckets are node buckets, so this should work
  stk::mesh::AllSelectedBucketsRange buckets_range = stk::mesh::get_buckets(selector, mesh);
  for (stk::mesh::AllSelectedBucketsIterator buckets_itr = boost::begin(buckets_range), buckets_end = boost::end(buckets_range);
       buckets_itr != buckets_end;
       ++buckets_itr) {
    get_buckets_range.push_back(*buckets_itr);
  }
  //BOOST_FOREACH( Bucket* bucket, stk::mesh::get_buckets(selector, mesh) ) {
  //  get_buckets_range.push_back(bucket);
  //}

  // get buckets selected by selector using alteranate range API
  // all buckets are node buckets, so this should work
  stk::mesh::AllBucketsRange all_buckets = stk::mesh::get_buckets(mesh);
  buckets_range = stk::mesh::get_buckets(selector, all_buckets);
  for (stk::mesh::AllSelectedBucketsIterator buckets_itr = boost::begin(buckets_range), buckets_end = boost::end(buckets_range);
       buckets_itr != buckets_end;
       ++buckets_itr) {
    get_buckets_alt_range.push_back(*buckets_itr);
  }
  //BOOST_FOREACH( Bucket* bucket, stk::mesh::get_buckets(selector, all_buckets) ) {
  //  get_buckets_alt_range.push_back(bucket);
  //}

  sort_and_compare_eq(get_buckets,           expected_buckets);
  sort_and_compare_eq(get_buckets_range,     expected_buckets);
  sort_and_compare_eq(get_buckets_alt_range, expected_buckets);

  //
  // Check entities
  //

  std::vector<Entity*> get_entities_range;

  stk::mesh::SelectedBucketRangeEntityIteratorRange selected_entity_range = stk::mesh::get_selected_entities( selector, all_buckets );
  for (stk::mesh::SelectedBucketRangeEntityIterator selected_entity_itr = boost::begin(selected_entity_range),
                                                    selected_entity_end = boost::end(selected_entity_range);
       selected_entity_itr != selected_entity_end;
       ++selected_entity_itr) {
    get_entities_range.push_back(*selected_entity_itr);
  }
  // TODO: Figure out why BOOST_FOREACH does not work well with selected entity iterators
  // BOOST_FOREACH(Entity* entity, stk::mesh::get_selected_entities(selector, all_buckets) ) {
  //   get_entities_range.push_back(entity);
  // }

  sort_and_compare_eq(get_entities_range, expected_selected_entities);
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
    std::vector<Bucket*> get_buckets_range;

    // Get all node buckets using range API

    stk::mesh::AllBucketsRange all_node_buckets = stk::mesh::get_buckets(NODE_RANK, fix.m_bulk_data);
    for (stk::mesh::AllBucketsIterator buckets_itr = boost::begin(all_node_buckets), buckets_end = boost::end(all_node_buckets);
         buckets_itr != buckets_end;
         ++buckets_itr) {
      get_buckets_range.push_back(*buckets_itr);
    }
    // BOOST_FOREACH( Bucket* bucket, stk::mesh::get_buckets(NODE_RANK, fix.m_bulk_data) ) {
    //   get_buckets_range.push_back(bucket);
    // }

    sort_and_compare_eq(get_buckets_range, node_buckets);
  }

  {
    std::vector<Entity*> expected_selected_entities;
    expected_selected_entities.push_back(fix.m_entity1);
    expected_selected_entities.push_back(fix.m_entity2);

    check_selected_buckets(Selector(fix.m_partA), fix.m_bulk_data, node_buckets, expected_selected_entities);
  }

  {
    std::vector<Entity*> expected_selected_entities;
    expected_selected_entities.push_back(fix.m_entity2);
    expected_selected_entities.push_back(fix.m_entity3);

    check_selected_buckets(Selector(fix.m_partB), fix.m_bulk_data, node_buckets, expected_selected_entities);
  }

  {
    std::vector<Entity*> expected_selected_entities;
    expected_selected_entities.push_back(fix.m_entity3);
    expected_selected_entities.push_back(fix.m_entity4);

    check_selected_buckets(Selector(fix.m_partC), fix.m_bulk_data, node_buckets, expected_selected_entities);
  }

  // Check get_entities correctness

  {
    std::vector<Entity*> all_nodes_expected, all_nodes, all_nodes_range;

    all_nodes_expected.push_back(fix.m_entity1);
    all_nodes_expected.push_back(fix.m_entity2);
    all_nodes_expected.push_back(fix.m_entity3);
    all_nodes_expected.push_back(fix.m_entity4);
    all_nodes_expected.push_back(fix.m_entity5);

    stk::mesh::get_entities( fix.m_bulk_data, NODE_RANK, all_nodes );
    sort_and_compare_eq(all_nodes, all_nodes_expected);

    stk::mesh::BucketVectorEntityIteratorRange entity_range = stk::mesh::get_entities(NODE_RANK, fix.m_bulk_data);
    for (stk::mesh::BucketVectorEntityIterator entity_itr = boost::begin(entity_range), entity_end = boost::end(entity_range);
         entity_itr != entity_end;
         ++entity_itr) {
      all_nodes_range.push_back(*entity_itr);
    }
    // BOOST_FOREACH( Entity* const entity, stk::mesh::get_entities( NODE_RANK, fix.m_bulk_data ) ) {
    //   all_nodes_range.push_back(entity);
    // }
    sort_and_compare_eq(all_nodes_range, all_nodes_expected);
  }
}

} // empty namespace
