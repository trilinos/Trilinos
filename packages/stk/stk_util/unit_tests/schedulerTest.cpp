#include <gtest/gtest.h>
#include <stk_util/util/Scheduler.hpp>

namespace
{

TEST(SchedulerTest, timeInterval)
{
    stk::util::Scheduler scheduler;

    const stk::util::Time startTime = 0.0;
    const stk::util::Time dt = 1.0;
    scheduler.add_interval(startTime, dt);

    const stk::util::Time terminationTime = 5.2;
    scheduler.set_termination_time(terminationTime);

    const stk::util::Step unusedStep = 0;
    EXPECT_TRUE(scheduler.is_it_time(0.0, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(0.5, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(1.0, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(1.1, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(2.5, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(2.6, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(2.999999999, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(3.0, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(4.0, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(4.0, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(5.0, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(terminationTime-0.1, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(terminationTime+0.1, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(terminationTime+0.2, unusedStep));
}

TEST(SchedulerTest, stepInterval)
{
    stk::util::Scheduler scheduler;

    const stk::util::Step startStep = 0;
    const stk::util::Step stepInterval = 2;
    scheduler.add_interval(startStep, stepInterval);

    stk::util::Time unusedTime = 0.0;
    const stk::util::Time dt = 0.1;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 0));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.is_it_time(unusedTime, 1));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 2));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.is_it_time(unusedTime, 3));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 4));
}

TEST(SchedulerTest, explicitStep)
{
    stk::util::Scheduler scheduler;

    const stk::util::Step startStep = 0;
    const stk::util::Step stepInterval = 2;
    scheduler.add_interval(startStep, stepInterval);
    const stk::util::Step additionalStepToOutput = 3;
    scheduler.add_explicit(additionalStepToOutput);

    stk::util::Time unusedTime = 0.0;
    const stk::util::Time dt = 0.1;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 0));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.is_it_time(unusedTime, 1));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 2));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 3));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.is_it_time(unusedTime, 4));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.is_it_time(unusedTime, 5));
}

TEST(SchedulerTest, forceWrite)
{
    stk::util::Scheduler scheduler;

    const stk::util::Time startTime = 0.0;
    const stk::util::Time dt = 1.0;
    scheduler.add_interval(startTime, dt);

    const stk::util::Step unusedStep = 0;
    EXPECT_TRUE(scheduler.is_it_time(0.0, unusedStep));
    scheduler.set_force_schedule();
    EXPECT_TRUE(scheduler.is_it_time(0.5, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(1.0, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(1.1, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(2.5, unusedStep));
    EXPECT_FALSE(scheduler.is_it_time(2.6, unusedStep));
    EXPECT_TRUE(scheduler.is_it_time(2.999999999, unusedStep));
}

TEST(SchedulerTest, emptyScheduler)
{
    stk::util::Scheduler scheduler;
    const stk::util::Time terminationTime = 5.0;
    scheduler.set_termination_time(terminationTime);
    EXPECT_FALSE(scheduler.is_it_time(0.0, 0));
    EXPECT_FALSE(scheduler.is_it_time(0.5, 1));
    EXPECT_FALSE(scheduler.is_it_time(3.5, 2));
    EXPECT_TRUE(scheduler.is_it_time(terminationTime, 2));
    EXPECT_FALSE(scheduler.is_it_time(terminationTime+0.5, 2));
}

}
