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
#include "stk_util/parallel/GenerateParallelUniqueIDs.hpp"  // for generate_parallel_unique_ids
#include "stk_util/parallel/Parallel.hpp"                   // for MPI_COMM_WORLD, MPI_Allreduce
#include "stk_util/parallel/ParallelVectorConcat.hpp"       // for parallel_vector_concat
#include <algorithm>                                        // for binary_search, sort
#include <cstdint>                                          // for uint64_t
#include <vector>                                           // for vector

//-------------------------------------------------------------------------------------------------------------------
// Verify that a result generated from the unique ids function meets all properties nessecary of unique ids
void  VerifyUniqueIds(uint64_t maxAllowableId, 
                          std::vector<uint64_t>& existingIdsLocal,
                          std::vector<uint64_t>&  newIdsLocal,
                          stk::ParallelMachine comm) {
  //
  //  Generate a global list of all ids created matched and all ids that existed perviously
  //
  std::vector<uint64_t> existingIdsGlobal;
  stk::parallel_vector_concat(comm, existingIdsLocal, existingIdsGlobal);
  std::sort(existingIdsGlobal.begin(), existingIdsGlobal.end());

  std::vector<uint64_t> newIdsGlobal;
  stk::parallel_vector_concat(comm, newIdsLocal, newIdsGlobal);

     
  //
  //  Verify properties of the returned ids
  //
  for(uint64_t i=0; i<newIdsGlobal.size(); ++i) {
    uint64_t curNewID = newIdsGlobal[i];
    //    1)  Verify all ids must be less than or equal to maxAllowableId
    EXPECT_LE(curNewID, maxAllowableId); 
    //    2)  Verify no overlap between existing and new id set
    EXPECT_FALSE(std::binary_search(existingIdsGlobal.begin(), existingIdsGlobal.end(), curNewID));
  }

}

//------------------------------------------------------------------------------
//
//  Test the generate_parallel_unique_ids function in stk_util/parallel
//
TEST(UnitTestParallel, testGenerateParallelUniqueIDs) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  //
  //  Test 1, nominal element death case
  //    Each processor starts with 1000 unique densly packed ids ordered 1->N.  Add 10 new ids per processor
  //    produced by kill every 100th element creating 10 new elements per processor.
  //
  uint64_t maxAllowableId = ~0U;

  {
    std::vector<uint64_t> newIds;
    std::vector<uint64_t> existingIds(1000);
    uint64_t firstId = mpi_rank*1000;

    for(uint64_t i=0; i<1000; ++i) {
      uint64_t currentElemId = firstId+i;
      existingIds[i] = currentElemId;
    }
    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, 10, MPI_COMM_WORLD);

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);
 
    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }
  //
  //  Test 2, same total source data as case 1, but source data scattered among processors
  //
  {
    int numRequest=0;
    std::vector<uint64_t> newIds;
    std::vector<uint64_t> existingIds;
    uint64_t destProc = mpi_size/2;    
    for(int i=0; i<1000*mpi_size; ++i) {
      destProc += i;
      if(destProc >= (uint64_t)mpi_size) {
        destProc = destProc%mpi_size;
      }
      if(destProc == (uint64_t)mpi_rank) {
        existingIds.push_back(i);
        if(i%100 == 0) {
          numRequest++;
        }
      }
    }
    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, numRequest, MPI_COMM_WORLD);

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);

    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }
  //
  //  Test 2.5, same as case one, except all new ids to be created on one processor
  //
    {
    int numRequest=0;
    std::vector<uint64_t> newIds;
    std::vector<uint64_t> existingIds;
    uint64_t destProc = mpi_size/2;    
    for(int i=0; i<1000*mpi_size; ++i) {
      destProc += i;
      if(destProc >= (uint64_t)mpi_size) {
        destProc = destProc%mpi_size;
      }
      if(destProc == (uint64_t)mpi_rank) {
        existingIds.push_back(i);
      }
    }

    if(mpi_rank == mpi_size/2) {
      numRequest = mpi_size*10;
    }


    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, numRequest, MPI_COMM_WORLD);

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);

    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }

  //
  //  Test 3, sparse fill
  //

  {
    std::vector<uint64_t> newIds;
    std::vector<uint64_t> existingIds(1000);
    uint64_t firstId = mpi_rank*1000;

    for(uint64_t i=0; i<1000; ++i) {
      uint64_t currentElemId = firstId+i;
      existingIds[i] = currentElemId;
    }

    if(mpi_rank == 0) {
      existingIds.push_back(~0U);
    }

    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, 10, MPI_COMM_WORLD);

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);
 
    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }


}
