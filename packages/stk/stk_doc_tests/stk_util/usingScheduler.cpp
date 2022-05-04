
// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

//BEGINSchedulerExample
#include "gtest/gtest.h"
#include "stk_util/environment/Scheduler.hpp"  // for Scheduler, Time, Step

namespace
{
TEST(StkUtilTestForDocumentation, TimeBasedScheduling)
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

TEST(StkUtilTestForDocumentation, TimeBasedSchedulingWithTerminationTime)
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

TEST(StkUtilTestForDocumentation, StepBasedScheduler)
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

TEST(StkUtilTestForDocumentation, TimeBasedSchedulerWithTwoTimeIntervals)
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
//ENDSchedulerExample
