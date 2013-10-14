#include <gtest/gtest.h>
#include <stk_util/util/Scheduler.hpp>

namespace
{
TEST(StkUtilTestForDocumentation, SchedulerWithTimeInterval)
{
    stk::util::Scheduler scheduler;

    const stk::util::Time startTime = 0.0;
    const stk::util::Time timeInterval = 1.0;
    scheduler.add_interval(startTime, timeInterval);

    stk::util::Step timeStep = 0;
    EXPECT_TRUE(scheduler.is_it_time(0.0, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(0.5, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(1.0, timeStep++));
}

TEST(StkUtilTestForDocumentation, SchedulerWithTerminationTime)
{
    stk::util::Scheduler scheduler;

    const stk::util::Time startTime = 2.0;
    const stk::util::Time timeInterval = 10.0;
    scheduler.add_interval(startTime, timeInterval);

    const stk::util::Time terminationTime = 8.2;
    scheduler.set_termination_time(terminationTime);

    stk::util::Step timeStep = 0;
    EXPECT_FALSE(scheduler.is_it_time(startTime - 1.0, timeStep++));
    const stk::util::Time firstTimeAfterStartTime = terminationTime-0.1;
    EXPECT_TRUE(scheduler.is_it_time(firstTimeAfterStartTime, timeStep++));
    const stk::util::Time firstTimeAfterTerminationTime = terminationTime+0.1;
    EXPECT_TRUE(scheduler.is_it_time(firstTimeAfterTerminationTime, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(terminationTime+0.2, timeStep++));
}

TEST(StkUtilTestForDocumentation, SchedulerWithStepInterval)
{
    stk::util::Scheduler scheduler;

    const stk::util::Step startStep = 0;
    const stk::util::Step stepInterval = 2;
    scheduler.add_interval(startStep, stepInterval);

    stk::util::Time time = 0.0;
    const stk::util::Time dt = 0.1;
    EXPECT_TRUE(scheduler.is_it_time(time, 0));
    time+=dt;
    EXPECT_FALSE(scheduler.is_it_time(time, 1));
    time+=dt;
    EXPECT_TRUE(scheduler.is_it_time(time, 2));
    time+=dt;
    EXPECT_FALSE(scheduler.is_it_time(time, 3));
    time+=dt;
    EXPECT_TRUE(scheduler.is_it_time(time, 4));
}
}
