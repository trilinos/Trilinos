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
#include <stk_util/parallel/ParallelVectorConcat.hpp>


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
  
  //
  //  Test 2, dense filling of initial half full space.
  //     Initially all even ids in use.  Note, this is likely the worst case for performance of this routine, but it should
  //     still return correct results.
  //
  
  {
    std::vector<unsigned> localIdsInUse;
    int targetRank = 0;
    for(unsigned i=0; i < 1000; ++i) {
      if(i%2 != 0) {
        targetRank += 5;
        continue;
      }
      if(targetRank >= mpi_size) {
        targetRank = 0;
      }
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
      }
      targetRank++;
    }

    std::vector<unsigned> returnIds;
    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 500, returnIds);
    
    EXPECT_EQ(returnCode, 0);
    EXPECT_EQ(returnIds.size(), (unsigned)500);

    std::sort(returnIds.begin(), returnIds.end());
    for(unsigned i=0; i<500; i++) {
      EXPECT_EQ(returnIds[i], i*2+1);
    }
    
  }  
  //
  //  Test 3,
  //  Boundary case, same as 3 but only 1 processor has any initial ids
  //
  {
    std::vector<unsigned> localIdsInUse;
    int targetRank = mpi_size/2;
    for(unsigned i=0; i < 1000; ++i) {
      if(i%2 != 0) {
        continue;
      }
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
      }
    }
    std::vector<unsigned> returnIds;
    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 500, returnIds);

    EXPECT_EQ(returnCode, 0);
    EXPECT_EQ(returnIds.size(), (unsigned)500);

    std::sort(returnIds.begin(), returnIds.end());
    for(unsigned i=0; i<500; i++) {
      EXPECT_EQ(returnIds[i], i*2+1);
    }
  }
  

  //
  //  Test 4, 
  //    Sparse filling spread id in a logarithm fashion clumped at zero
  //
  
  {
    int maxId = 1000000;
    int numId = 1024;

    int curStart = 0;

    int targetRank = mpi_size/2;

    std::vector<unsigned> localIdsInUse;


    while(numId > 0) {
      numId = numId/2;
      for(int i=0; i<numId; ++i) {
        if(mpi_rank == targetRank) {
          localIdsInUse.push_back(curStart+i);
        }
        targetRank++;
        if(targetRank >= mpi_size){
          targetRank = 0;
        }        
      }
      curStart = curStart + (maxId-curStart)/2;
    }

    //
    //  Successively fill to 10% capacity
    //
    for(int iter=0; iter<10; ++iter) {
      std::vector<unsigned> returnIds;
      int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, maxId, localIdsInUse, 10000, returnIds);

      EXPECT_EQ(returnCode, 0);
      EXPECT_EQ(returnIds.size(), (unsigned)10000);
       

      std::vector<unsigned> globalInUse;

      stk::parallel_vector_concat(MPI_COMM_WORLD, localIdsInUse, globalInUse);

      globalInUse.insert(globalInUse.end(), returnIds.begin(), returnIds.end());
       
      std::sort(globalInUse.begin(), globalInUse.end());

      for(unsigned i=0; i<globalInUse.size() -1; ++i) {
        EXPECT_NE(globalInUse[i], globalInUse[i+1]);
      }
      for(unsigned iret=0; iret<returnIds.size(); ++iret) {
        if(mpi_rank == targetRank) {
          localIdsInUse.push_back(returnIds[iret]);
        }
        targetRank++;
      }



    }    
  }  


}
