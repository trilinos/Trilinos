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

#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/parallel/GenerateParallelConsistentIDs.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string



TEST(UnitTestParallel, testGenerateParallelConsistentIDs) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  //
  //  Test 1, nominal element death case
  //    Each processor starts with 1000 unique densly packed ids ordered 1->N.  Add 10 new ids per processor
  //    produced by kill every 100th element creating 10 new elements per processor.
  //
  unsigned maxAllowableId = ~0U;
  std::vector<unsigned> existingIds(1000);
  unsigned firstId = mpi_rank*1000;

  std::vector<unsigned> localOrderArray;
  std::vector<unsigned> newIds;

  for(unsigned i=0; i<1000; ++i) {
    unsigned currentElemId = firstId+i;
    existingIds[i] = currentElemId;
    if(i%100 == 0) {
      //  local order array is the 'source element' which provides a unique repeatable ordering
      localOrderArray.push_back(currentElemId);
      newIds.push_back(0);
    }
  }
  int returnCount = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray, newIds, MPI_COMM_WORLD);
  EXPECT_EQ(returnCount, mpi_size*10);
}
