//BEGIN
#include <gtest/gtest.h>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/util/perf_util.hpp>

namespace stk {

TEST(stkMeshHowTo, makeNightlyTrackedTest)
{
  double start_time = stk::cpu_time();

  // <do-something>

  double first_timer = stk::cpu_time() - start_time;

  start_time = stk::cpu_time();

  // <do-something-else>

  double second_timer = stk::cpu_time() - start_time;

  double total_time = first_timer + second_timer;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] =
    {first_timer, second_timer, total_time};
  const char* timer_names[NUM_TIMERS] =
    {"First Timer Name", "Second Timer Name", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);
}

}
//END
