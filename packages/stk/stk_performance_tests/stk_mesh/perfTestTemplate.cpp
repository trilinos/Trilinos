//BEGIN
#include <stk_util/unit_test_support/perf_unit_util.hpp>

#include <gtest/gtest.h>

namespace {

TEST(stkMeshHowTo, makePerformanceUnitTest)
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  // Set up problem, this is not timed

  CALLGRIND_TOGGLE_COLLECT;

  // Do thing you want to test

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  // Clean up if needed, this is not timed

  PERFORMANCE_TEST_POSTAMBLE();
}

}
//END
