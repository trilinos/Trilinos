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
#include "stk_util/environment/memory_util.hpp"

TEST(MemoryUtil, maxMinAvg)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs > 2) { return; }

  double localValue = 1.0*(stk::parallel_machine_rank(comm)+1);
  double max = 0.0;
  double min = 0.0;
  double avg = 0.0;
  stk::get_max_min_avg(comm, localValue, max, min, avg);

  double expectedMax = numProcs==1 ? 1.0 : 2.0;
  double expectedMin = 1.0;
  double expectedAvg = numProcs==1 ? 1.0 : 1.5;

  EXPECT_NEAR(expectedMax, max, 1.e-12);
  EXPECT_NEAR(expectedMin, min, 1.e-12);
  EXPECT_NEAR(expectedAvg, avg, 1.e-12);
}

TEST(MemoryUtil, maxMinAvg_null_comm)
{
  MPI_Comm comm = stk::parallel_machine_null();

  double localValue = 1.0;
  double max = 0.0;
  double min = 0.0;
  double avg = 0.0;
  stk::get_max_min_avg(comm, localValue, max, min, avg);

  EXPECT_NEAR(1.0, max, 1.e-12);
  EXPECT_NEAR(1.0, min, 1.e-12);
  EXPECT_NEAR(1.0, avg, 1.e-12);
}

