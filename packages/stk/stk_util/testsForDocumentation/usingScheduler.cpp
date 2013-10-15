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
    const stk::util::Time firstAfterTermination = terminationTime+0.1;
    EXPECT_TRUE(scheduler.is_it_time(firstAfterTermination, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(terminationTime+0.2, timeStep++));
}

TEST(StkUtilTestForDocumentation, SchedulerIntervalCheck)
{
    stk::util::Scheduler scheduler;

    const stk::util::Step startStep = 0;
    const stk::util::Step stepInterval = 4;
    scheduler.add_interval(startStep, stepInterval);

    const stk::util::Time dt = 0.1;
    for (stk::util::Step timeStep=0;timeStep<100;timeStep+=3)
    {
        stk::util::Time time = timeStep*dt;
        bool check = scheduler.is_it_time(time, timeStep);
        if ( timeStep % stepInterval == 0 )
        {
            EXPECT_TRUE(check);
        }
        else
        {
            EXPECT_FALSE(check);
        }
    }
}

TEST(StkUtilTestForDocumentation, SchedulerWithTwoTimeIntervals)
{
    stk::util::Scheduler scheduler;
    const stk::util::Time startTime1 = 0.0;
    const stk::util::Time delta1 = 0.1;
    scheduler.add_interval(startTime1, delta1);
    const stk::util::Time startTime2 = 0.9;
    const stk::util::Time delta2 = 0.3;
    scheduler.add_interval(startTime2, delta2);

    stk::util::Step timeStep = 0;
    EXPECT_TRUE(scheduler.is_it_time(0.0, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(0.07, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(0.14, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(0.62, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(0.6999999, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(0.77, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(0.9, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(0.97, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(1.04, timeStep++));
    EXPECT_FALSE(scheduler.is_it_time(1.11, timeStep++));
    EXPECT_TRUE(scheduler.is_it_time(1.27, timeStep++));
}
}
