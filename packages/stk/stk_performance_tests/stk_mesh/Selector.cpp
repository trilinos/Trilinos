#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fixtures/SelectorFixture.hpp>

#include <stdexcept>

namespace {

using stk::mesh::fixtures::VariableSelectorFixture;

}

STKUNIT_UNIT_TEST(selector_timings, selector_timings)
{
  // A scaling test for selector and bucket operations

  // If we are running with STL in debug mode we shrink the problem
  // down in order to keep things running in a reasonable amount of
  // time.
#ifdef _GLIBCXX_DEBUG
  size_t N = 1000;
#else
  size_t N = 100000;
#endif

  VariableSelectorFixture fix(N);

  std::vector<double> selector_creation(N/2);
  std::vector<double> get_buckets_usage(N/2);
  size_t total_buckets_grabbed = 0;

  double total_selector_time = 0.0, total_bucket_time = 0.0;
  size_t timing_index = 0;
  for (size_t n = 1 ; n<N; n*=2, ++timing_index) {
    // Selector creation
    double start_time = stk::wall_time();
    stk::mesh::Selector selectUnion;
    for (size_t part_i = 0 ; part_i<n ; ++part_i) {
      selectUnion |= *fix.m_declared_part_vector[part_i];
    }
    selector_creation[timing_index] = stk::wall_dtime(start_time);
    total_selector_time += selector_creation[timing_index];

    // Selector usage:
    start_time = stk::wall_time();
    unsigned entity_rank = 0;
    stk::mesh::BucketVector const& buckets_out =  fix.m_BulkData.get_buckets(entity_rank, selectUnion);
    total_buckets_grabbed += buckets_out.size();
    get_buckets_usage[timing_index] = stk::wall_dtime(start_time);
    total_bucket_time += get_buckets_usage[timing_index];
  }

  // Print out table
  std::cout << "total_buckets_grabbed: " << total_buckets_grabbed << std::endl;
  std::cout << "\"N\" \"selector_creation_time\" \"get_buckets_time\" " << std::endl;
  timing_index = 0;
  for (size_t n = 1 ; n<N; n*=2) {
    std::cout << n << " " << selector_creation[timing_index] << " " << get_buckets_usage[timing_index] << std::endl;
    ++timing_index;
  }

  const double total_time = total_selector_time + total_bucket_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {total_selector_time, total_bucket_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Selector unions", "Get buckets", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}
