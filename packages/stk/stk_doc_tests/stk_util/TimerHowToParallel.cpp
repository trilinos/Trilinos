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

//BEGINTimerParallelTests
#include "gtest/gtest.h"
#include "stk_unit_test_utils/stringAndNumberComparisons.hpp"  // for areStringsEqualWithToleran...
#include "stk_util/diag/PrintTimer.hpp"                        // for printTimersTable
#include "stk_util/diag/Timer.hpp"                             // for Timer, createRootTimer
#include "stk_util/diag/TimerMetricTraits.hpp"                 // for METRICS_ALL
#include "stk_util/parallel/Parallel.hpp"                      // for MPI_Comm_rank, MPI_Comm_size
#include <unistd.h>                                            // for usleep
#include <iosfwd>                                              // for ostringstream
#include <string>                                              // for string

namespace
{

const double tolerance = 0.10;

void doWork()
{
  ::usleep(1e5);
}

TEST(StkDiagTimerHowTo, useTimersInParallel)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs == 2)
  {
    enum {CHILDMASK1 = 1};
    stk::diag::TimerSet enabledTimerSet(CHILDMASK1);
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
    rootTimer.start();

    stk::diag::Timer childTimer1("childTimer1", CHILDMASK1, rootTimer);

    {
      stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(childTimer1, communicator);
      doWork();
    }

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint, communicator);

    int procId = stk::parallel_machine_rank(communicator);
    if(procId == 0)
    {
      std::string expectedOutput = "                                                  \
                                CPU Time            CPU Time            CPU Time            \
                               Wall Time           Wall Time           Wall Time            \
     Timer      Count      Sum (% of System)   Min (% of System)   Max (% of System)        \
                           Sum (% of System)   Min (% of System)   Max (% of System)        \
SKIP ----- --------------------- --------------------- SKIP SKIP SKIP --------------------- \
totalTestRuntime  2            SKIP  SKIP          SKIP  SKIP          SKIP  SKIP           \
                               00:00:00.200 SKIP          00:00:00.100 SKIP          00:00:00.100 SKIP           \
  childTimer1     2            SKIP  SKIP          SKIP  SKIP          SKIP  SKIP           \
                               00:00:00.200 SKIP          00:00:00.100 SKIP          00:00:00.100 SKIP           \
                                                                                            \
Took SKIP seconds to generate the table above.                                            \
                    ";
std::cerr<<expectedOutput<<" : "<<outputStream.str()<<std::endl;
      using stk::unit_test_util::areStringsEqualWithToleranceForNumbers;
      EXPECT_TRUE(areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));
    }

    stk::diag::deleteRootTimer(rootTimer);
  }
}

}
//ENDTimerParallelTests
