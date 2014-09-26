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
#include <stk_util/parallel/ParallelIndexGapFinder.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string
#include <stk_util/environment/ReportHandler.hpp>

//------------------------------------------------------------------------------
//
//  Test the generate_parallel_unique_ids function in stk_util/parallel
//

TEST(UnitTestParallel, testParallelIndexGapFinder) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  //
  //  Test 1, dense filling of initially half full space
  //    First 500 of 1000 ids in use(0->500), return the remaining 500 (500->999)
  //
  {
    std::vector<unsigned> localIdsInUse;

    int targetRank = 0;

    for(unsigned i=0; i < 500; ++i) {
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
      }
      targetRank++;
      if(targetRank >= mpi_size) {
        targetRank = 0;
      }
    }
    std::vector<unsigned> returnIds;
    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 500, returnIds);

    std::sort(returnIds.begin(), returnIds.end());


    EXPECT_EQ(returnCode, 0);
    EXPECT_EQ(returnIds.size(), (unsigned)500);

    for(int i=0; i<500; ++i) {
      EXPECT_EQ(returnIds[i], (unsigned)(500+i));
    }
  }
}
