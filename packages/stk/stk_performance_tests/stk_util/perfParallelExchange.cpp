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
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <vector>
#include <iostream>
#include <chrono>
#include <thread>

void fill_comm_data(int numProcs, int myProc,
                    unsigned BANDWIDTH, unsigned NUM_VALUES,
                    std::vector<int>& sendOffsets,
                    std::vector<double>& sendData,
                    std::vector<int>& recvOffsets,
                    std::vector<double>& recvData)
{
  sendOffsets.assign(numProcs+1, 0);

  const int halfBandwidth = BANDWIDTH/2;
  const int rawFirstDestProc = myProc - halfBandwidth;
  const int rawLastDestProc = myProc + halfBandwidth;
  int firstDestProc = std::max(0, rawFirstDestProc);
  if (firstDestProc == myProc) {
    firstDestProc += 1;
  }
 
  int lastDestProc = std::min(numProcs-1, rawLastDestProc);
  if (lastDestProc == myProc) {
    lastDestProc -= 1;
  }

  //first, put send-sizes in sendOffsets:
  for(int p=firstDestProc; p<=lastDestProc; ++p) {
    if (p != myProc) {
      sendOffsets[p] = NUM_VALUES;
    }
  }

  recvOffsets = sendOffsets;

  //the send/recv procs are 'symmetric', but data-sizes will be non-symmetric:
  //add some variation to the send-sizes
  for(int p=firstDestProc; p<=lastDestProc; ++p) {
    if (p != myProc) {
      sendOffsets[p] += p;
      recvOffsets[p] += myProc;
    }
  }

  //now, convert sizes to offsets:
  unsigned sendOffset = 0;
  unsigned recvOffset = 0;
  for(unsigned p=0; p<sendOffsets.size(); ++p) {
    unsigned size = sendOffsets[p];
    sendOffsets[p] = sendOffset;
    sendOffset += size;

    size = recvOffsets[p];
    recvOffsets[p] = recvOffset;
    recvOffset += size;
  }

  sendData.assign(sendOffset, 42.0);
  recvData.assign(recvOffset, 0);
}

TEST(PerfPllDataExchange, nonsym_knownsizes)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if (numProcs < 3) { GTEST_SKIP(); }

  const unsigned NUM_BATCHES = 5;
  const unsigned NUM_RUNS = stk::unit_test_util::get_command_line_option("-r", 5);
  const unsigned defaultBandwidth = std::max(2, numProcs/5);
  const unsigned BANDWIDTH = stk::unit_test_util::get_command_line_option("-b", defaultBandwidth);
  const unsigned NUM_VALUES = stk::unit_test_util::get_command_line_option("-v", 20);

  const int myProc = stk::parallel_machine_rank(comm);

  if (myProc == 0) {
    std::cout << "Using NUM_RUNS=" << NUM_RUNS
              << ", BANDWIDTH=" << BANDWIDTH
              << ", NUM_VALUES=" << NUM_VALUES
              << std::endl;
  }
  
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  std::vector<int> sendOffsets, recvOffsets;
  std::vector<double> sendData, recvData;

  fill_comm_data(numProcs, myProc, BANDWIDTH, NUM_VALUES, sendOffsets, sendData, recvOffsets, recvData);

  for(unsigned b=0; b<NUM_BATCHES; ++b) {
    batchTimer.start_batch_timer();

    for(unsigned r=0; r<NUM_RUNS; ++r) {

      constexpr bool checkInput = false;
      stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendData.data(),
                                                     recvOffsets.data(), recvData.data(),
                                                     comm, checkInput);
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_RUNS);
}

