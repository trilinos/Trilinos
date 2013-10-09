#include <gtest/gtest.h>

#include <stk_io/IOScheduler.hpp>

namespace {

TEST(SchedulerTest, timeInterval)
{
    stk::io::IOScheduler scheduler;

    const stk::io::Time startTime = 0.0;
    const stk::io::Time dt = 1.0;
    scheduler.add_interval(startTime, dt);

    const stk::io::Step unusedStep = 0;
    const stk::io::Time terminationTime = 5.2;
    EXPECT_TRUE(scheduler.write_now(0.0, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(0.5, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(1.0, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(1.1, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(2.5, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(2.6, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(2.999999999, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(5.0, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(terminationTime-0.1, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(terminationTime+0.1, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(terminationTime+0.2, unusedStep, terminationTime));
}

TEST(SchedulerTest, stepInterval)
{
    stk::io::IOScheduler scheduler;

    const stk::io::Step startStep = 0;
    const stk::io::Step stepInterval = 2;
    scheduler.add_interval(startStep, stepInterval);

    stk::io::Time unusedTime = 0.0;
    const stk::io::Time terminationTime = 5.0;
    const stk::io::Time dt = 0.1;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 0, terminationTime));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.write_now(unusedTime, 1, terminationTime));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 2, terminationTime));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.write_now(unusedTime, 3, terminationTime));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 4, terminationTime));
}

TEST(SchedulerTest, explicitStep)
{
    stk::io::IOScheduler scheduler;

    const stk::io::Step startStep = 0;
    const stk::io::Step stepInterval = 2;
    scheduler.add_interval(startStep, stepInterval);
    const stk::io::Step additionalStepToOutput = 3;
    scheduler.add_explicit(additionalStepToOutput);

    stk::io::Time unusedTime = 0.0;
    const stk::io::Time terminationTime = 5.0;
    const stk::io::Time dt = 0.1;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 0, terminationTime));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.write_now(unusedTime, 1, terminationTime));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 2, terminationTime));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 3, terminationTime));
    unusedTime+=dt;
    EXPECT_TRUE(scheduler.write_now(unusedTime, 4, terminationTime));
    unusedTime+=dt;
    EXPECT_FALSE(scheduler.write_now(unusedTime, 5, terminationTime));
}

TEST(SchedulerTest, forceWrite)
{
    stk::io::IOScheduler scheduler;

    const stk::io::Time startTime = 0.0;
    const stk::io::Time dt = 1.0;
    scheduler.add_interval(startTime, dt);

    const stk::io::Step unusedStep = 0;
    const stk::io::Time terminationTime = 5.0;
    EXPECT_TRUE(scheduler.write_now(0.0, unusedStep, terminationTime));
    scheduler.set_force_write();
    EXPECT_TRUE(scheduler.write_now(0.5, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(1.0, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(1.1, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(2.5, unusedStep, terminationTime));
    EXPECT_FALSE(scheduler.write_now(2.6, unusedStep, terminationTime));
    EXPECT_TRUE(scheduler.write_now(2.999999999, unusedStep, terminationTime));
}

TEST(SchedulerTest, emptyScheduler)
{
    stk::io::IOScheduler scheduler;
    const stk::io::Time terminationTime = 5.0;
    EXPECT_FALSE(scheduler.write_now(0.0, 0, terminationTime));
    EXPECT_FALSE(scheduler.write_now(0.5, 1, terminationTime));
    EXPECT_FALSE(scheduler.write_now(3.5, 2, terminationTime));
    EXPECT_TRUE(scheduler.write_now(terminationTime, 2, terminationTime));
    EXPECT_FALSE(scheduler.write_now(terminationTime+0.5, 2, terminationTime));
}

}
