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

#include <string>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include <stk_unit_test_utils/timer.hpp>
#include <vector>
#include <iostream>
#include <chrono>
#include <thread>

TEST(PrintTimersTable, performance)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  const unsigned NUM_BATCHES = 5;
  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for(unsigned b=0; b<NUM_BATCHES; ++b) {
    batchTimer.start_batch_timer();

    for(unsigned r=0; r<NUM_RUNS; ++r) {
      constexpr int CHILDMASK1 = 1;
      using namespace std::chrono_literals;

      stk::diag::TimerSet enabledTimerSet(CHILDMASK1);
      stk::diag::Timer rootTimer(stk::diag::createRootTimer("root-timer", enabledTimerSet));
      stk::diag::Timer child1("child1-timer", CHILDMASK1, rootTimer);
      stk::diag::Timer child2("child2-timer", CHILDMASK1, rootTimer);
      stk::diag::Timer child3("child3-timer", CHILDMASK1, rootTimer);
      stk::diag::Timer child4("child4-timer", CHILDMASK1, rootTimer);
      stk::diag::Timer child5("child5-timer", CHILDMASK1, rootTimer);

      rootTimer.start();
      std::this_thread::sleep_for(20ms);

      {
        stk::diag::TimeBlock timeBlock1(child1, comm);
        std::this_thread::sleep_for(100ms);

        {
          stk::diag::TimeBlock timeBlock2(child2, comm);
          std::this_thread::sleep_for(100ms);

          {
            stk::diag::TimeBlock timeBlock3(child3, comm);
            std::this_thread::sleep_for(100ms);
    
            {
              stk::diag::TimeBlock timeBlock4(child4, comm);
              std::this_thread::sleep_for(100ms);
    
              {
                stk::diag::TimeBlock timeBlock5(child5, comm);
                std::this_thread::sleep_for(100ms);
              }
            }
          }
        }
      }

      std::ostringstream oss;
      stk::diag::printTimersTable(oss, rootTimer, stk::diag::METRICS_ALL, false, comm);

      stk::diag::deleteRootTimer(rootTimer);
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_RUNS);
}

