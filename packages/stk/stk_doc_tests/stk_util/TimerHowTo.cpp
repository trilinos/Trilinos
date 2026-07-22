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

//BEGINTimerTests
#include "gtest/gtest.h"
#include "stk_unit_test_utils/stringAndNumberComparisons.hpp"  // for areStringsEqualWithToleran...
#include "stk_util/diag/PrintTimer.hpp"                        // for printTimersTable
#include "stk_util/diag/Timer.hpp"                             // for Timer, createRootTimer
#include "stk_util/diag/TimerMetricTraits.hpp"                 // for METRICS_ALL
#include <unistd.h>                                            // for usleep
#include <iosfwd>                                              // for ostringstream
#include <string>                                              // for string

namespace
{

#if defined(NDEBUG)
  const double tolerance = 0.10;
#else
  const double tolerance = 0.25;
#endif

void doWork()
{
  ::usleep(1e5);
}

TEST(StkDiagTimerHowTo, useTheRootTimer)
{
  stk::diag::TimerSet enabledTimerSet(0);
  stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);

  {
    stk::diag::TimeBlock totalTestRuntime(rootTimer);
    doWork();

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                     \
                 Timer                   Count       CPU Time              Wall Time       \
---------------------------------------- ----- --------------------- --------------------- \
totalTestRuntime                           1        SKIP  SKIP             00:00:00.100 SKIP      \
                                                                                           \
Took 0.0001 seconds to generate the table above.                                           \
    ";
    using stk::unit_test_util::areStringsEqualWithToleranceForNumbers;
    EXPECT_TRUE(areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));
  }

  stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, useChildTimers)
{
  enum {CHILDMASK1 = 1, CHILDMASK2 = 2};
  stk::diag::TimerSet enabledTimerSet(CHILDMASK1 | CHILDMASK2);
  stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
  rootTimer.start();

  stk::diag::Timer childTimer1("childTimer1", CHILDMASK1, rootTimer);
  stk::diag::Timer childTimer2("childTimer2", CHILDMASK2, rootTimer);

  {
    stk::diag::TimeBlock timeStuffInThisScope(childTimer1);
    stk::diag::TimeBlock timeStuffInThisScopeAgain(childTimer2);
    doWork();
  }

  std::ostringstream outputStream;
  bool printTimingsOnlySinceLastPrint = false;
  stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

  {
    stk::diag::TimeBlock timeStuffInThisScope(childTimer1);
    doWork();
  }

  stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

  std::string expectedOutput = "                                                         \
                 Timer                   Count       CPU Time              Wall Time       \
---------------------------------------- ----- --------------------- --------------------- \
totalTestRuntime                             1        SKIP   SKIP        00:00:00.100 SKIP        \
  childTimer1                                1        SKIP   SKIP        00:00:00.100 SKIP        \
  childTimer2                                1        SKIP   SKIP        00:00:00.100 SKIP        \
                                                                                           \
Took 0.0001 seconds to generate the table above.                                           \
                 Timer                   Count       CPU Time              Wall Time       \
---------------------------------------- ----- --------------------- --------------------- \
totalTestRuntime                             1        SKIP   SKIP        00:00:00.200 SKIP        \
  childTimer1                                2        SKIP   SKIP        00:00:00.200 SKIP        \
  childTimer2                                1        SKIP   SKIP        00:00:00.100 SKIP        \
                                                                                           \
Took 0.0001 seconds to generate the table above.                                           \
            ";
  using stk::unit_test_util::areStringsEqualWithToleranceForNumbers;
  EXPECT_TRUE(areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));

  stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, disableChildTimers)
{
  enum {CHILDMASK1 = 1, CHILDMASK2 = 2};
  stk::diag::TimerSet enabledTimerSet(CHILDMASK2);
  stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
  rootTimer.start();

  stk::diag::Timer disabledTimer("disabledTimer", CHILDMASK1, rootTimer);
  stk::diag::Timer enabledTimer("enabledTimer", CHILDMASK2, rootTimer);

  {
    stk::diag::TimeBlock timeStuffInThisScope(disabledTimer);
    stk::diag::TimeBlock timeStuffInThisScopeAgain(enabledTimer);
    doWork();
  }

  std::ostringstream outputStream;
  bool printTimingsOnlySinceLastPrint = false;
  stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

  {
    stk::diag::TimeBlock timeStuffInThisScope(disabledTimer);
    doWork();
  }

  stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

  std::string expectedOutput = "                                                         \
                 Timer                   Count       CPU Time              Wall Time       \
---------------------------------------- ----- --------------------- --------------------- \
totalTestRuntime                             1        SKIP   SKIP        00:00:00.100 SKIP        \
  enabledTimer                               1        SKIP   SKIP        00:00:00.100 SKIP        \
                                                                                           \
Took 0.0001 seconds to generate the table above.                                           \
                 Timer                   Count       CPU Time              Wall Time       \
---------------------------------------- ----- --------------------- --------------------- \
totalTestRuntime                             1        SKIP   SKIP        00:00:00.200 SKIP        \
  enabledTimer                               1        SKIP   SKIP        00:00:00.100 SKIP        \
                                                                                           \
Took 0.0001 seconds to generate the table above.                                           \
            ";
  using stk::unit_test_util::areStringsEqualWithToleranceForNumbers;
  EXPECT_TRUE(areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));

  stk::diag::deleteRootTimer(rootTimer);
}

}
//ENDTimerTests
