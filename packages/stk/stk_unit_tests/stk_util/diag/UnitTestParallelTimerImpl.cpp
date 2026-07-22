#include "gtest/gtest.h"
#include "stk_util/diag/ParallelTimerImpl.hpp"
#include "stk_util/diag/Timer.hpp"
#include "stk_util/diag/TimerImpl.hpp"
#include "stk_util/diag/TimerMetricTraits.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace {
stk::diag::impl::ParallelTimer create_timer(const std::string& name, double val)
{
  stk::diag::impl::ParallelTimer timer;
  timer.m_name            = name;
  timer.m_cpuTime.m_value = val;
  timer.m_cpuTime.m_sum   = val;
  timer.m_cpuTime.m_min   = val;
  timer.m_cpuTime.m_max   = val;

  return timer;
}
}

namespace stk::diag {

class TimerTester
{
  public:
    TimerTester(Timer& timer) :
      m_timer(timer)
    {}

    double getCPUTime() const
    {
      return m_timer.getMetric<CPUTime>().m_accumulatedLap;
    }

    void setCPUTime(double val)
    {
      m_timer.m_timerImpl->m_cpuTime.m_accumulatedLap = val;
    }

  private:
    Timer& m_timer;
};
}

TEST(ParallelTimer, MergeSingleLevelTimers)
{
  double val1 = 1.0, val2 = 2.0;
  stk::diag::impl::ParallelTimer t1 = create_timer("timer1", val1);
  stk::diag::impl::ParallelTimer t2 = create_timer("timer2", val2);

  stk::diag::impl::merge_parallel_timer(t1, t2, false);

  EXPECT_EQ(t1.m_cpuTime.m_value, val1);
  EXPECT_EQ(t1.m_cpuTime.m_sum,   val1 + val2);
  EXPECT_EQ(t1.m_cpuTime.m_min,   val1);
  EXPECT_EQ(t1.m_cpuTime.m_max,   val2);
}

TEST(ParallelTimer, MergeTwoLevelTimers)
{
  double val1 = 1.0, val2 = 2.0, val3 = 3.0, val4 = 4.0;
  stk::diag::impl::ParallelTimer t1 = create_timer("timer1", val1);
  stk::diag::impl::ParallelTimer t2 = create_timer("timer2", val2);
  stk::diag::impl::ParallelTimer t3 = create_timer("timer3", val3);
  stk::diag::impl::ParallelTimer t4 = create_timer("timer2", val4);

  t1.m_subtimerList.push_back(t2);
  t3.m_subtimerList.push_back(t4);

  stk::diag::impl::merge_parallel_timer(t1, t3, false);

  EXPECT_EQ(t1.m_cpuTime.m_value, val1);
  EXPECT_EQ(t1.m_cpuTime.m_sum,   val1 + val3);
  EXPECT_EQ(t1.m_cpuTime.m_min,   val1);
  EXPECT_EQ(t1.m_cpuTime.m_max,   val3);
  EXPECT_EQ(t1.m_subtimerList.size(), 1U);

  stk::diag::impl::ParallelTimer t2Merged = t1.m_subtimerList.front();
  EXPECT_EQ(t2Merged.m_cpuTime.m_value, val2);
  EXPECT_EQ(t2Merged.m_cpuTime.m_sum,   val2 + val4);
  EXPECT_EQ(t2Merged.m_cpuTime.m_min,   val2);
  EXPECT_EQ(t2Merged.m_cpuTime.m_max,   val4);
}

TEST(ParallelTimer, MergeTwoLevelTimersDifferentNames)
{
  double val1 = 1.0, val2 = 2.0, val3 = 3.0, val4 = 4.0;
  stk::diag::impl::ParallelTimer t1 = create_timer("timer1", val1);
  stk::diag::impl::ParallelTimer t2 = create_timer("timer2", val2);
  stk::diag::impl::ParallelTimer t3 = create_timer("timer3", val3);
  stk::diag::impl::ParallelTimer t4 = create_timer("timer4", val4);

  t1.m_subtimerList.push_back(t2);
  t3.m_subtimerList.push_back(t4);

  stk::diag::impl::merge_parallel_timer(t1, t3, false);

  EXPECT_EQ(t1.m_cpuTime.m_value, val1);
  EXPECT_EQ(t1.m_cpuTime.m_sum,   val1 + val3);
  EXPECT_EQ(t1.m_cpuTime.m_min,   val1);
  EXPECT_EQ(t1.m_cpuTime.m_max,   val3);
  EXPECT_EQ(t1.m_subtimerList.size(), 2U);

  stk::diag::impl::ParallelTimer t2Copy = t1.m_subtimerList.front();
  EXPECT_EQ(t2Copy.m_cpuTime.m_value, val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_sum,   val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_min,   val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_max,   val2);

  stk::diag::impl::ParallelTimer t4Copy = *(++t1.m_subtimerList.begin());
  EXPECT_EQ(t4Copy.m_cpuTime.m_value, val4);
  EXPECT_EQ(t4Copy.m_cpuTime.m_sum,   val4);
  EXPECT_EQ(t4Copy.m_cpuTime.m_min,   val4);
  EXPECT_EQ(t4Copy.m_cpuTime.m_max,   val4);
}

TEST(ParallelTimer, MergeThreeLevelTimers)
{
  double val1 = 1.0, val2 = 2.0, val3 = 3.0, val4 = 4.0;
  double val5 = 5.0, val6 = 6.0;
  stk::diag::impl::ParallelTimer t1 = create_timer("timer1", val1);
  stk::diag::impl::ParallelTimer t2 = create_timer("timer2", val2);
  stk::diag::impl::ParallelTimer t3 = create_timer("timer3", val3);
  stk::diag::impl::ParallelTimer t4 = create_timer("timer1", val4);
  stk::diag::impl::ParallelTimer t5 = create_timer("timer2", val5);
  stk::diag::impl::ParallelTimer t6 = create_timer("timer3", val6);

  t2.m_subtimerList.push_back(t3);
  t1.m_subtimerList.push_back(t2);

  t5.m_subtimerList.push_back(t6);
  t4.m_subtimerList.push_back(t5);

  stk::diag::impl::merge_parallel_timer(t1, t4, false);

  EXPECT_EQ(t1.m_cpuTime.m_value, val1);
  EXPECT_EQ(t1.m_cpuTime.m_sum,   val1 + val4);
  EXPECT_EQ(t1.m_cpuTime.m_min,   val1);
  EXPECT_EQ(t1.m_cpuTime.m_max,   val4);
  EXPECT_EQ(t1.m_subtimerList.size(), 1U);

  stk::diag::impl::ParallelTimer t2Merged = t1.m_subtimerList.front();
  EXPECT_EQ(t2Merged.m_cpuTime.m_value, val2);
  EXPECT_EQ(t2Merged.m_cpuTime.m_sum,   val2 + val5);
  EXPECT_EQ(t2Merged.m_cpuTime.m_min,   val2);
  EXPECT_EQ(t2Merged.m_cpuTime.m_max,   val5);
  EXPECT_EQ(t2Merged.m_subtimerList.size(), 1U);

  stk::diag::impl::ParallelTimer t3Merged = t2Merged.m_subtimerList.front();
  EXPECT_EQ(t3Merged.m_cpuTime.m_value, val3);
  EXPECT_EQ(t3Merged.m_cpuTime.m_sum,   val3 + val6);
  EXPECT_EQ(t3Merged.m_cpuTime.m_min,   val3);
  EXPECT_EQ(t3Merged.m_cpuTime.m_max,   val6);
  EXPECT_EQ(t3Merged.m_subtimerList.size(), 0U);
}

TEST(ParallelTimer, MergeThreeLevelTimersDifferentNames)
{
  double val1 = 1.0, val2 = 2.0, val3 = 3.0, val4 = 4.0;
  double val5 = 5.0, val6 = 6.0;
  stk::diag::impl::ParallelTimer t1 = create_timer("timer1", val1);
  stk::diag::impl::ParallelTimer t2 = create_timer("timer2", val2);
  stk::diag::impl::ParallelTimer t3 = create_timer("timer3", val3);
  stk::diag::impl::ParallelTimer t4 = create_timer("timer4", val4);
  stk::diag::impl::ParallelTimer t5 = create_timer("timer5", val5);
  stk::diag::impl::ParallelTimer t6 = create_timer("timer6", val6);

  t2.m_subtimerList.push_back(t3);
  t1.m_subtimerList.push_back(t2);

  t5.m_subtimerList.push_back(t6);
  t4.m_subtimerList.push_back(t5);

  stk::diag::impl::merge_parallel_timer(t1, t4, false);

  EXPECT_EQ(t1.m_cpuTime.m_value, val1);
  EXPECT_EQ(t1.m_cpuTime.m_sum,   val1 + val4);
  EXPECT_EQ(t1.m_cpuTime.m_min,   val1);
  EXPECT_EQ(t1.m_cpuTime.m_max,   val4);
  EXPECT_EQ(t1.m_subtimerList.size(), 2U);

  stk::diag::impl::ParallelTimer t2Copy = t1.m_subtimerList.front();
  EXPECT_EQ(t2Copy.m_cpuTime.m_value, val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_sum,   val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_min,   val2);
  EXPECT_EQ(t2Copy.m_cpuTime.m_max,   val2);
  EXPECT_EQ(t2Copy.m_subtimerList.size(), 1U);

  stk::diag::impl::ParallelTimer t5Copy = t1.m_subtimerList.back();
  EXPECT_EQ(t5Copy.m_cpuTime.m_value, val5);
  EXPECT_EQ(t5Copy.m_cpuTime.m_sum,   val5);
  EXPECT_EQ(t5Copy.m_cpuTime.m_min,   val5);
  EXPECT_EQ(t5Copy.m_cpuTime.m_max,   val5);
  EXPECT_EQ(t5Copy.m_subtimerList.size(), 1U);

  stk::diag::impl::ParallelTimer t3Copy = t2Copy.m_subtimerList.front();
  EXPECT_EQ(t3Copy.m_cpuTime.m_value, val3);
  EXPECT_EQ(t3Copy.m_cpuTime.m_sum,   val3);
  EXPECT_EQ(t3Copy.m_cpuTime.m_min,   val3);
  EXPECT_EQ(t3Copy.m_cpuTime.m_max,   val3);
  EXPECT_EQ(t3Copy.m_subtimerList.size(), 0U);

  stk::diag::impl::ParallelTimer t6Copy = t5Copy.m_subtimerList.front();
  EXPECT_EQ(t6Copy.m_cpuTime.m_value, val6);
  EXPECT_EQ(t6Copy.m_cpuTime.m_sum,   val6);
  EXPECT_EQ(t6Copy.m_cpuTime.m_min,   val6);
  EXPECT_EQ(t6Copy.m_cpuTime.m_max,   val6);
  EXPECT_EQ(t6Copy.m_subtimerList.size(), 0U);
}


TEST(ParallelTimer, CollectTimersChunkSize1)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  int commSize = stk::parallel_machine_size(comm);
  int commRank = stk::parallel_machine_rank(comm);
  double cpuTimeVal = commRank;
  stk::diag::Timer rootTimer = stk::diag::createRootTimer("rootTimer", stk::diag::TimerSet(stk::diag::getEnabledTimerMetricsMask()));
  stk::diag::TimerTester(rootTimer).setCPUTime(cpuTimeVal);

  const int maxProcsPerGather = 1;
  stk::diag::impl::ParallelTimer parallelTimer = stk::diag::impl::collect_timers(rootTimer, false, comm, maxProcsPerGather);

  if (commRank == 0)
  {
    EXPECT_EQ(parallelTimer.m_cpuTime.m_min, 0.0);
    EXPECT_EQ(parallelTimer.m_cpuTime.m_max, commSize - 1);
    EXPECT_EQ(parallelTimer.m_cpuTime.m_sum, commSize * (commSize - 1) / 2.0);
  }

  stk::diag::deleteRootTimer(rootTimer);
}

TEST(ParallelTimer, CollectTimersChunkSize2)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  int commSize = stk::parallel_machine_size(comm);
  int commRank = stk::parallel_machine_rank(comm);
  double cpuTimeVal = commRank + 1;
  stk::diag::Timer rootTimer = stk::diag::createRootTimer("rootTimer", stk::diag::TimerSet(stk::diag::getEnabledTimerMetricsMask()));
  stk::diag::TimerTester(rootTimer).setCPUTime(cpuTimeVal);

  const int maxProcsPerGather = 2;
  stk::diag::impl::ParallelTimer parallelTimer = stk::diag::impl::collect_timers(rootTimer, false, comm, maxProcsPerGather);

  if (commRank == 0)
  {
    EXPECT_EQ(parallelTimer.m_cpuTime.m_min, 1.0);
    EXPECT_EQ(parallelTimer.m_cpuTime.m_max, commSize);
    EXPECT_EQ(parallelTimer.m_cpuTime.m_sum, commSize * (1 + commSize) / 2.0);
  }

  stk::diag::deleteRootTimer(rootTimer);
}