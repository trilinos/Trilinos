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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <vector>

//BEGINAllReduce
TEST(AllReduce, combinedOps)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  int myRank = stk::parallel_machine_rank(comm);

  constexpr int NumMin = 2;
  constexpr int NumMax = 3;
  constexpr int NumSum = 4;

  std::vector<int> ints(NumMin, myRank);
  std::vector<float> floats(NumMax, 1.0*myRank);
  std::vector<float> doubles(NumSum, 1.0*(myRank+1));

  stk::all_reduce(comm, stk::ReduceMin<NumMin>(ints.data())
                  & stk::ReduceMax<NumMax>(floats.data())
                  & stk::ReduceSum<NumSum>(doubles.data()));

  const int expectedMin = 0;
  const float expectedMax = 1.0;
  const double expectedSum = 3.0;

  for(int thisInt : ints) {
    EXPECT_EQ(expectedMin, thisInt);
  }

  for(float thisFloat : floats) {
    EXPECT_EQ(expectedMax, thisFloat);
  }

  for(double thisDouble : doubles) {
    EXPECT_EQ(expectedSum, thisDouble);
  }
}
//ENDAllReduce

