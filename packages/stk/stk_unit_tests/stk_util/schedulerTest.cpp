// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>
#include <stk_util/environment/Scheduler.hpp>

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
