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
#include "stk_util/parallel/Parallel.hpp"        // for parallel_machine_rank, parallel_machine_...
#include "stk_util/parallel/ParallelReduce.hpp"  // for all_reduce_max, all_reduce_maxloc, all_r...
#include "stk_util/stk_config.h"                 // for STK_HAS_MPI
#include <cstdint>                               // for uint64_t, uint32_t
#include <limits>                                // for numeric_limits
#include <vector>                                // for vector


#if defined ( STK_HAS_MPI )
TEST(ParallelComm, AllReduceLoc)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    int myProcId = stk::parallel_machine_rank(comm);
    int numProcs = stk::parallel_machine_size(comm);

    int nvalues=5;
    uint64_t limit_32bit_integer = std::numeric_limits<uint32_t>::max();

    std::vector<double> values(nvalues);
    std::vector<uint64_t> locations(nvalues);
    std::vector<double> globalValues(nvalues);
    std::vector<uint64_t> globalLocations(nvalues);

    for(int n = 0; n < nvalues; n++) {
      values[n] = myProcId*numProcs+n;
      locations[n] = limit_32bit_integer + myProcId*numProcs*numProcs+n; //Want to test when outside 32-bit range for index.
    }

    stk::all_reduce_maxloc(comm, values.data(), locations.data(), globalValues.data(), globalLocations.data(), nvalues);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ((numProcs-1)*numProcs+n, globalValues[n]);
      EXPECT_EQ(limit_32bit_integer + (numProcs-1)*numProcs*numProcs+n, globalLocations[n]);
    }

    stk::all_reduce_minloc(comm, values.data(), locations.data(), globalValues.data(), globalLocations.data(), nvalues);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n, globalValues[n]);
      EXPECT_EQ(limit_32bit_integer+n, globalLocations[n]);
    }
}

TEST(ParallelComm, AllReduce)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    int myProcId = stk::parallel_machine_rank(comm);
    int numProcs = stk::parallel_machine_size(comm);

    int nvalues=5;

    std::vector<double> values(nvalues);
    std::vector<double> globalValues(nvalues);

    for(int n = 0; n < nvalues; n++) {
      values[n] = myProcId*numProcs+n;
    }

    stk::all_reduce_max(comm, values.data(), globalValues.data(), nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ((numProcs-1)*numProcs+n, globalValues[n]);
    }

    stk::all_reduce_min(comm, values.data(), globalValues.data(), nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n, globalValues[n]);
    }

    int alpha = 0;
    for(int n = 0; n < numProcs; ++n) {
      alpha += n*numProcs;
    }

    stk::all_reduce_sum(comm, values.data(), globalValues.data(), nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n*numProcs + alpha, globalValues[n]);
    }
}

TEST(ParallelComm, GetGlobal)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    int myProcId = stk::parallel_machine_rank(comm);
    int numProcs = stk::parallel_machine_size(comm);

    double localValue = myProcId * numProcs;

    double globalMax = stk::get_global_max(comm, localValue);
    EXPECT_EQ(globalMax, (numProcs-1)*numProcs);

    double globalMin = stk::get_global_min(comm, localValue);
    EXPECT_EQ(globalMin, 0);

    int expectedSum = 0;
    for(int n = 0; n < numProcs; ++n) {
      expectedSum += n*numProcs;
    }

    double globalSum = stk::get_global_sum(comm, localValue);
    EXPECT_EQ(globalSum, expectedSum);
}

TEST(ParallelComm, uint64_reduce_min)
{
  // Test coverage of openmpi-4.1.4 bug with unsigned integer min reduction
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  int myProcId = stk::parallel_machine_rank(comm);

  const uint64_t goldValue = 10;

  uint64_t id = (myProcId == 0) ? goldValue : std::numeric_limits<uint64_t>::max();
  uint64_t gid;

  stk::all_reduce_min(comm, &id, &gid, 1u);

  EXPECT_EQ(goldValue, gid);
}

#endif

