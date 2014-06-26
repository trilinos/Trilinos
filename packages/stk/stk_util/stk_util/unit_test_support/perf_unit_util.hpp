#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <gtest/gtest.h>

#include <valgrind/callgrind.h>

/**
 * Useful utilities for performance unit testing.
 */

namespace stk {

#define PERFORMANCE_TEST_PREAMBLE(expected_np)          \
  CALLGRIND_START_INSTRUMENTATION;                      \
                                                        \
  stk::check_valgrind_version();                        \
                                                        \
  stk::ParallelMachine pm = MPI_COMM_WORLD;             \
                                                        \
  const int p_size = stk::parallel_machine_size(pm);    \
                                                        \
  ASSERT_EQ(p_size, expected_np)


#define PERFORMANCE_TEST_POSTAMBLE()              \
  stk::print_memory_sum_all_procs(pm);            \
                                                  \
  stk::print_debug_skip(pm)


inline void print_memory_sum_all_procs(stk::ParallelMachine pm)
{
  const size_t p_rank = stk::parallel_machine_rank(pm);
  size_t my_peak = stk::allocator_memory_usage<void>::peak_memory();
  size_t peak_sum = 0;
  int err = MPI_Reduce(&my_peak, &peak_sum, 1 /*size*/, MPI_LONG_LONG, MPI_SUM, 0 /*root*/, pm);
  ASSERT_EQ(err, MPI_SUCCESS);

  if (p_rank == 0) {
    std::cout << "\nSTKPERF peak memory sum: " << peak_sum << std::endl;
  }
}

inline void check_valgrind_version()
{
  ASSERT_EQ(__VALGRIND_MAJOR__, 3);
  ASSERT_EQ(__VALGRIND_MINOR__, 8);
}

inline void print_debug_skip(stk::ParallelMachine pm)
{
#ifndef NDEBUG
  // We're in debug; need to tell test script not to validate cycle count
  const size_t p_rank = stk::parallel_machine_rank(pm);
  if (p_rank == 0) {
    std::cout << "\nSTKPERF SKIP VALIDATION" << std::endl;
  }
#endif
}

} // namespace stk
