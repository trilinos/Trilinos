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

#include "gtest/gtest.h"
#include "stk_unit_test_utils/getOption.h"      // for get_command_line_option
#include "stk_util/diag/PrintTimer.hpp"         // for printTimersTable
#include "stk_util/diag/Timer.hpp"              // for Timer, TimeBlock, Timer::Metric, createRo...
#include "stk_util/diag/TimerMetricTraits.hpp"  // for LapCount (ptr only), MetricTraits, Metric...
#include "stk_util/parallel/Parallel.hpp"       // for MPI_COMM_WORLD
#include <cstddef>                              // for size_t
#include <chrono>                               // for milliseconds
#include <cmath>                                // for sin
#include <iostream>                             // for ostringstream, operator<<, basic_ostream
#include <memory>                               // for allocator_traits<>::value_type
#include <string>                               // for string, char_traits
#include <thread>                               // for sleep_for
#include <vector>                               // for vector

using stk::unit_test_util::get_command_line_option;

enum {
  TIMER_DOMAIN		= 0x00001000,		///< Enable domain timers
  TIMER_REGION		= 0x00002000,		///< Enable region timers
  TIMER_PROCEDURE	= 0x00004000,		///< Enable procedure timers
  TIMER_MECHANICS	= 0x00008000,		///< Enable mechanics timers
  TIMER_ALGORITHM	= 0x00010000,		///< Enable algorithm timers
  TIMER_SOLVER		= 0x00020000,		///< Enable solver timers
  TIMER_CONTACT		= 0x00040000,		///< Enable contact timers
  TIMER_MATERIAL	= 0x00080000,		///< Enable material timers
  TIMER_SEARCH		= 0x00100000,		///< Enable search timers
  TIMER_TRANSFER	= 0x00200000,		///< Enable transfer timers
  TIMER_ADAPTIVITY	= 0x00400000 		///< Enable adaptivity
};

enum {
  TIMER_UNUSED_1	= 0x00001000,		///< Enable unused 1
  TIMER_PROFILE_1	= 0x00002000,		///< Enable profile 1 timers
  TIMER_PROFILE_2	= 0x00004000,		///< Enable profile 2 timers
  TIMER_PROFILE_3	= 0x00008000,		///< Enable profile 3 timers
  TIMER_PROFILE_4	= 0x00010000,		///< Enable profile 4 timers
  TIMER_APP_1		= 0x00020000,		///< Enable application defined 1
  TIMER_APP_2		= 0x00040000,		///< Enable application defined 2
  TIMER_APP_3		= 0x00080000		///< Enable application defined 3
};

namespace {

double
work(int repetitions = 1000)
{
  double x = 1.0;

  for (int i = 0; i < repetitions; ++i) {
    double di = i;
    x += std::sin(di);
  }
  return x;
}


stk::diag::TimerSet &
unitTestTimerSet()
{
  static stk::diag::TimerSet s_unitTestTimerSet(TIMER_REGION);

  return s_unitTestTimerSet;
}


stk::diag::TimerSet &
unitTestSecondTimerSet()
{
  static stk::diag::TimerSet s_unitTestSecondTimerSet(TIMER_APP_3);

  return s_unitTestSecondTimerSet;
}


stk::diag::Timer &unitTestTimer() {
  const std::string name("Unit test timer");
  static stk::diag::Timer s_unitTestTimer (stk::diag::createRootTimer(name, unitTestTimerSet()));

  return s_unitTestTimer;
}


struct RootObject
{
  RootObject()
    : m_timer("Root object", TIMER_REGION, unitTestTimer())
  {}

  stk::diag::Timer      m_timer;
};


struct Object 
{
  Object(const std::string &name, RootObject &root_object) 
    : m_id(0),
      m_name(name),
      m_timer(name, root_object.m_timer)
  {}
  
  Object(int id, const Object &parent) 
    : m_id(id),
      m_name(id_name(id)),
      m_timer(m_name, parent.m_timer)
  {}
  
  static std::string id_name(int id) {
    std::ostringstream s;
    s << "Object id " << id << " run";
    return s.str();
  }

  void run() {
    stk::diag::TimeBlock _time(m_timer);
    m_x += work();
  }
  
  int                   m_id;
  std::string           m_name;
  stk::diag::Timer      m_timer;
  double                m_x;
};

} // namespace <empty>

constexpr double millisecTolerance = 1.e-2; //10 millis for robustness

TEST(UnitTestTimer, UnitTest)
{
  stk::diag::TimeBlock root_time_block(unitTestTimer());

  std::ostringstream strout;
  
  // Create subtimer and test lap time
  {
    static stk::diag::Timer lap_timer("One second Wall time twice", unitTestTimer());
    
    stk::diag::TimeBlock _time(lap_timer);    
    double x = work(100);
    static_cast<void>(x);
    std::ostringstream oss;
    oss << x << std::endl;

    std::this_thread::sleep_for(std::chrono::milliseconds(20));

    lap_timer.lap();
    
    stk::diag::MetricTraits<stk::diag::WallTime>::Type lap_time = lap_timer.getMetric<stk::diag::WallTime>().getLap();
  
    EXPECT_GE(lap_time, (0.02 - millisecTolerance));

    std::this_thread::sleep_for(std::chrono::milliseconds(20));

    lap_timer.stop();
    
    lap_time = lap_timer.getMetric<stk::diag::WallTime>().getLap();
  
    EXPECT_GE(lap_time, (0.04 - millisecTolerance));
  }

  {
    static stk::diag::Timer run_timer("Run 100 times twice", unitTestTimer());
    
    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    ASSERT_TRUE(lap_count == 100);
  }

  // Create second timer set
  {
    static stk::diag::Timer second_timer("Second timer set", unitTestTimer(), unitTestSecondTimerSet());
    static stk::diag::Timer second_timer_on_default("On default", second_timer);
    static stk::diag::Timer second_timer_on("On", TIMER_APP_3, second_timer);
    static stk::diag::Timer second_timer_off("Off", TIMER_APP_1, second_timer);
    
    stk::diag::TimeBlock _time(second_timer);
    stk::diag::TimeBlock _time1(second_timer_on_default);
    stk::diag::TimeBlock _time2(second_timer_on);
    stk::diag::TimeBlock _time3(second_timer_off);

    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  // Grab previous subtimer and run 100 laps
  {
    static stk::diag::Timer run_timer("Run 100 times twice", unitTestTimer());
    
    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    ASSERT_TRUE(lap_count == 200);
  }

  // Create root object
  RootObject root_object;
    
  {
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = root_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    ASSERT_TRUE(lap_count == 0);
  }

  // Create object
  {
    Object time_object("One object", root_object);
    
    for (int i = 0; i < 100; ++i) {
      time_object.run();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = time_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    ASSERT_TRUE(lap_count == 100);
  }

  // Create object tree
  {
    std::vector<Object> object_vector;
    object_vector.push_back(Object("Object Tree", root_object));

    int id = 0;
    for (size_t i = 0; i < 2; ++i) {
      size_t ix = object_vector.size();
      object_vector.push_back(Object(id++, object_vector[0]));
      for (size_t j = 0; j < 2; ++j) {
        size_t jx = object_vector.size();
        object_vector.push_back(Object(id++, object_vector[ix]));
        for (int k = 0; k < 2; ++k) {    
          object_vector.push_back(Object(id++, object_vector[jx]));
        }
      }
    }
    
    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);
    
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    EXPECT_EQ(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(0));

    for (size_t j = 0; j < object_vector.size(); ++j) 
      object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);    

    //root node has 15, children node 7, grandchildren 3, greatgrandchildren 1
    // 15 + 2*7 + 4*3 + 8*1 = 49 for laptimes that are cumulative
    EXPECT_EQ(15u, object_vector[0].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap());
    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    EXPECT_EQ(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(49));

    for (size_t i = 1; i < 100; ++i) 
      for (size_t j = 0; j < object_vector.size(); ++j) 
        object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);    

    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    EXPECT_EQ(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(4900));

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, true);

    for (size_t i = 1; i < 100; ++i) 
      for (size_t j = 0; j < object_vector.size(); ++j) 
        object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, true);
  }
}

TEST(UnitTestTimer, YuugeNumberOfTimers)
{
    std::ostringstream strout;

    stk::diag::TimerSet unitTestSecondTimerSet(TIMER_APP_3);
    const std::string name("Unit Test Timer");
    stk::diag::Timer rootTimer (stk::diag::createRootTimer(name, unitTestSecondTimerSet));
    {
      stk::diag::TimeBlock root_time_block(rootTimer);

      unsigned numTimers = 100;
      numTimers = stk::unit_test_util::get_command_line_option("-numTimers", numTimers);
      {
          static std::vector<stk::diag::Timer> lap_timers;
          for (unsigned i = 0; i < numTimers; ++i)
          {
              std::ostringstream os;
              os << "Really Great Timer For Index with lots of space in between the name and index" << i;
              lap_timers.emplace_back(os.str(), rootTimer);
              stk::diag::TimeBlock _time(lap_timers[i]);
              // work(1);
          }
          for (unsigned i = 0; i < numTimers; ++i)
          {
              stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = lap_timers[i].getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
              EXPECT_EQ(1u, lap_count);
          }
      }

      stk::diag::printTimersTable(strout, rootTimer, stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);
    }
    stk::diag::deleteRootTimer(rootTimer);
}

TEST(UnitTestTimer, lapCount_2start_1stop)
{
  stk::diag::TimerSet timerSet(TIMER_APP_3);
  stk::diag::Timer rootTimer(stk::diag::createRootTimer("Root Timer", timerSet));
  rootTimer.start();
  rootTimer.start();
  rootTimer.stop();
  unsigned expectedLapCount = 1;
  EXPECT_EQ(expectedLapCount, rootTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap());
  stk::diag::deleteRootTimer(rootTimer);
}

TEST(UnitTestTimer, MultipleStarts)
{
  stk::diag::TimerSet timerSet(TIMER_APP_3);
  stk::diag::Timer rootTimer(stk::diag::createRootTimer("Root Timer", timerSet));
  {
    stk::diag::TimeBlock root_time_block(rootTimer);
    stk::diag::Timer childTimer("Child Timer", rootTimer);
  
    std::ostringstream strout;
    
    {
      stk::diag::TimeBlock childBlock(childTimer);
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    {
      childTimer.start();
      childTimer.start();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      childTimer.stop();
      childTimer.stop();
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    double expectedTime = 0.2;
    double measuredTime = childTimer.getMetric<stk::diag::WallTime>().getAccumulatedLap();
    double relativeError = (std::abs(expectedTime - measuredTime)/expectedTime);
    double relativeTolerance = 0.1;
    EXPECT_TRUE(relativeError < relativeTolerance);
    EXPECT_EQ(childTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap(), 2u);
    stk::diag::printTimersTable(strout, rootTimer, stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);
  }
  stk::diag::deleteRootTimer(rootTimer);
}

TEST(UnitTestTimer, NestedTimers)
{
  stk::diag::TimerSet timerSet(TIMER_APP_3);
  stk::diag::Timer rootTimer(stk::diag::createRootTimer("Root Timer", timerSet));
  {
    stk::diag::TimeBlock root_time_block(rootTimer);
    stk::diag::Timer childTimer("Child Timer", rootTimer);
    stk::diag::Timer grandChildTimer(" Grandchild Timer", childTimer);

    std::ostringstream strout;
    
    {
      stk::diag::TimeBlock childBlock(childTimer);
      stk::diag::TimeBlock grandChildBlock(grandChildTimer);
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    {
      stk::diag::TimeBlock lap_timer(grandChildTimer);
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    EXPECT_GE(childTimer.getMetric<stk::diag::WallTime>().getAccumulatedLap(), (0.1 - millisecTolerance));
    EXPECT_EQ(childTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap(), 2u);
    stk::diag::printTimersTable(strout, rootTimer, stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);
  }
  stk::diag::deleteRootTimer(rootTimer);
}

TEST(UnitTestTimer, ReferenceCounting)
{
  stk::diag::TimerSet timerSet(TIMER_APP_3);
  stk::diag::Timer rootTimer(stk::diag::createRootTimer("Root Timer", timerSet));
  {
    stk::diag::TimeBlock root_time_block(rootTimer);
    stk::diag::Timer childTimer("Child Timer", rootTimer);
    stk::diag::Timer grandChildTimer("Grandchild Timer", childTimer);
    stk::diag::Timer grandChildTimer2("Grandchild Timer 2", childTimer);

    std::ostringstream strout;
    
    {
      stk::diag::TimeBlock grandChildBlock(grandChildTimer);
      stk::diag::TimeBlock grandChildBlock2(grandChildTimer);
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    EXPECT_GE(childTimer.getMetric<stk::diag::WallTime>().getAccumulatedLap(), 0.01);
    EXPECT_EQ(childTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap(), 1u);

    {
      stk::diag::TimeBlock grandChildBlock(grandChildTimer);
      {
        stk::diag::TimeBlock grandChildBlock2(grandChildTimer2);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    EXPECT_GE(childTimer.getMetric<stk::diag::WallTime>().getAccumulatedLap(), 0.03);
    EXPECT_EQ(childTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap(), 2u);

    {
      stk::diag::TimeBlock childBlock(childTimer);
      {
        stk::diag::TimeBlock grandChildBlock(grandChildTimer);
        {
          stk::diag::TimeBlock grandChildBlock2(grandChildTimer2);
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    EXPECT_GE(childTimer.getMetric<stk::diag::WallTime>().getAccumulatedLap(), 0.06);
    EXPECT_EQ(childTimer.getMetric<stk::diag::LapCount>().getAccumulatedLap(), 3u);

    stk::diag::printTimersTable(strout, rootTimer, stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);
  }
  stk::diag::deleteRootTimer(rootTimer);
}
