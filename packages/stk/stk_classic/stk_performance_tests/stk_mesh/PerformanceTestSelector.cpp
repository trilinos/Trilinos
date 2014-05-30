#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

namespace {

using stk_classic::mesh::fixtures::VariableSelectorFixture ;

}

STKUNIT_UNIT_TEST( PerformanceTestSelector, start)
{
  size_t N = 5;
  VariableSelectorFixture fix(N);
  stk_classic::mesh::Selector selectUnion;
  for (size_t part_i = 0 ; part_i<N ; ++part_i) {
    selectUnion |= *fix.m_declared_part_vector[part_i];
  }
  std::vector<stk_classic::mesh::Bucket*> buckets_out;
  unsigned entity_rank = 0;
  get_buckets(selectUnion, fix.m_BulkData.buckets(entity_rank), buckets_out);
  STKUNIT_ASSERT_EQUAL( buckets_out.size(), N );
  // Construct once for large N
  // Graph time for get_buckets against 1..N
}

STKUNIT_UNIT_TEST( PerformanceTestSelector, timings)
{
  // Construction

  // If we are running with STL in debug mode we shrink the problem
  // down in order to keep things running in a reasonable amount of
  // time.
#ifdef _GLIBCXX_DEBUG
  size_t N = 1000;
#else
  size_t N = 10000;
#endif

  VariableSelectorFixture fix(N);

  std::vector<double> selector_creation(N/2);
  std::vector<double> get_buckets_usage(N/2);
  double start_time = stk_classic::wall_time();
  size_t timing_index = 0;
  for (size_t n = 1 ; n<N; n*=2) {
    // Selector creation
    stk_classic::mesh::Selector selectUnion;
    for (size_t part_i = 0 ; part_i<n ; ++part_i) {
      selectUnion |= *fix.m_declared_part_vector[part_i];
    }
    selector_creation[timing_index] = stk_classic::wall_dtime(start_time);

    // Selector usage:
    std::vector<stk_classic::mesh::Bucket*> buckets_out;
    unsigned entity_rank = 0;
    get_buckets(selectUnion, fix.m_BulkData.buckets(entity_rank), buckets_out);
    get_buckets_usage[timing_index] = stk_classic::wall_dtime(start_time);
    ++timing_index;
  }

  // Print out table
  std::cout << "\"N\" \"selector_creation_time\" \"get_buckets_time\" " << std::endl;
  timing_index = 0;
  for (size_t n = 1 ; n<N; n*=2) {
    std::cout << n << " " << selector_creation[timing_index] << " " << get_buckets_usage[timing_index] << std::endl;
    ++timing_index;
  }
  STKUNIT_EXPECT_TRUE(true);
}
