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
#include "stk_util/parallel/Parallel.hpp"                // for MPI_COMM_WORLD, MPI_Scan, MPI_SUM
#include "stk_util/parallel/ParallelIndexGapFinder.hpp"  // for parallel_index_gap_finder_global
#include "stk_util/parallel/ParallelVectorConcat.hpp"    // for parallel_vector_concat
#include <algorithm>                                     // for sort
#include <cstdint>                                       // for uint64_t
#include <vector>                                        // for vector<>::iterator, vector


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
    std::vector<uint64_t> localIdsInUse;
    uint64_t numLocalIdsRequested = 0;
    int targetRank = 0;
    for(uint64_t i=0; i < 500; ++i) {
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
        numLocalIdsRequested++;
      }
      targetRank++;
      if(targetRank >= mpi_size) {
        targetRank = 0;
      }
    }
    std::vector<uint64_t> returnIds;
    uint64_t myFirstNewId;

    MPI_Scan(&numLocalIdsRequested, &myFirstNewId, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    myFirstNewId-=numLocalIdsRequested;


    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 
                                                           numLocalIdsRequested, 500, myFirstNewId, returnIds);


    std::vector<uint64_t> globalReturnIds;
    stk::parallel_vector_concat(MPI_COMM_WORLD, returnIds, globalReturnIds);


    std::sort(globalReturnIds.begin(), globalReturnIds.end());
     
    EXPECT_EQ(returnCode, 0);

    EXPECT_EQ(globalReturnIds.size(), (uint64_t)500);

    
    for(int i=0; i<500; ++i) {
      EXPECT_EQ(globalReturnIds[i], (uint64_t)(500+i));
    }
  }
  
  
  //
  //  Test 2, dense filling of initial half full space.
  //     Initially all even ids in use.  Note, this is likely the worst case for performance of this routine, but it should
  //     still return correct results.
  //
  {
    std::vector<uint64_t> localIdsInUse;
    uint64_t numLocalIdsRequested = 0;
    int targetRank = 0;
    for(uint64_t i=0; i < 1000; ++i) {
      if(i%2 != 0) {
        targetRank += 5;
        continue;
      }
      if(targetRank >= mpi_size) {
        targetRank = 0;
      }
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
        numLocalIdsRequested++;
      }
      targetRank++;
    }

    

    std::vector<uint64_t> returnIds;




    uint64_t myFirstNewId;
    MPI_Scan(&numLocalIdsRequested, &myFirstNewId, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    myFirstNewId-=numLocalIdsRequested;
    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 
                                                           numLocalIdsRequested, 500, myFirstNewId, returnIds);
    std::vector<uint64_t> globalReturnIds;
    stk::parallel_vector_concat(MPI_COMM_WORLD, returnIds, globalReturnIds);

    EXPECT_EQ(returnCode, 0);
    EXPECT_EQ(globalReturnIds.size(), (uint64_t)500);
    std::sort(globalReturnIds.begin(), globalReturnIds.end());
    for(uint64_t i=0; i<500; i++) {
      EXPECT_EQ(globalReturnIds[i], i*2+1);
    }
  }  
  
  
  //
  //  Test 3,
  //  Boundary case, same as 3 but only 1 processor has any initial ids, only 1 processor will recieve any ids
  //
  {
    std::vector<uint64_t> localIdsInUse;
    int targetRank = mpi_size/2;
    for(uint64_t i=0; i < 1000; ++i) {
      if(i%2 != 0) {
        continue;
      }
      if(mpi_rank == targetRank) {
        localIdsInUse.push_back(i);
      }
    }
    std::vector<uint64_t> returnIds;
 
    uint64_t numLocalIdsRequested;
    if(mpi_rank == mpi_size/3) {
      numLocalIdsRequested = 500;
    } else {
      numLocalIdsRequested = 0;
    }



    uint64_t myFirstNewId;
    MPI_Scan(&numLocalIdsRequested, &myFirstNewId, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    myFirstNewId-=numLocalIdsRequested;
    int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, 1000, localIdsInUse, 
                                                           numLocalIdsRequested, 500, myFirstNewId, returnIds);
    std::vector<uint64_t> globalReturnIds;
    stk::parallel_vector_concat(MPI_COMM_WORLD, returnIds, globalReturnIds);

    
    if(mpi_rank == mpi_size/3) {
      EXPECT_EQ(returnIds.size(), (uint64_t)500);
    } else {
      EXPECT_EQ(returnIds.size(), (uint64_t)0);
    }

    EXPECT_EQ(returnCode, 0);
    EXPECT_EQ(globalReturnIds.size(), (uint64_t)500);
    std::sort(globalReturnIds.begin(), globalReturnIds.end());
    for(uint64_t i=0; i<500; i++) {
      EXPECT_EQ(globalReturnIds[i], i*2+1);
    }
  }
  
  
  //
  //  Test 4, 
  //    Sparse filling spread id in a logarithm fashion clumped at zero
  //
  {
    uint64_t maxId = 1000000;
    uint64_t numId = 1024;
    uint64_t curStart = 0;
    int targetRank = mpi_size/2;
    std::vector<uint64_t> localIdsInUse;
    while(numId > 0) {
      numId = numId/2;
      for(uint64_t i=0; i<numId; ++i) {
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


    uint64_t numLocalIdsRequested=0;
    for(uint64_t i=0; i<10000; ++i) {
      if(i%mpi_size == (uint64_t)mpi_rank) {
        numLocalIdsRequested++;
      }
    }


    //
    //  Successively fill to 10% capacity
    //
    for(int iter=0; iter<10; ++iter) {
      std::vector<uint64_t> returnIds;




      uint64_t myFirstNewId;
      MPI_Scan(&numLocalIdsRequested, &myFirstNewId, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      myFirstNewId-=numLocalIdsRequested;
      int returnCode = stk::parallel_index_gap_finder_global(MPI_COMM_WORLD, 0, maxId, localIdsInUse,
                                                      numLocalIdsRequested, 10000, myFirstNewId, returnIds);
      std::vector<uint64_t> globalReturnIds;
      stk::parallel_vector_concat(MPI_COMM_WORLD, returnIds, globalReturnIds);





      EXPECT_EQ(returnCode, 0);
      EXPECT_EQ(globalReturnIds.size(), (uint64_t)10000);
      std::vector<uint64_t> globalInUse;
      stk::parallel_vector_concat(MPI_COMM_WORLD, localIdsInUse, globalInUse);
      globalInUse.insert(globalInUse.end(), globalReturnIds.begin(), globalReturnIds.end());
      std::sort(globalInUse.begin(), globalInUse.end());
      for(uint64_t i=0; i<globalInUse.size() -1; ++i) {
        EXPECT_NE(globalInUse[i], globalInUse[i+1]);
      }
      for(uint64_t iret=0; iret<globalReturnIds.size(); ++iret) {
        if(mpi_rank == targetRank) {
          localIdsInUse.push_back(globalReturnIds[iret]);
        }
        targetRank++;
      }
    }    
  }  
  

}
