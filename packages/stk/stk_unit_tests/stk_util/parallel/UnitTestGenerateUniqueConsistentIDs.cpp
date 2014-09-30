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
#include <stk_util/parallel/GenerateParallelUniqueIDs.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string
#include <stk_util/environment/ReportHandler.hpp>


//-------------------------------------------------------------------------------------------------------------------
// Verify that a result generated from the unique ids function meets all properties nessecary of unique ids
void  VerifyUniqueIds(unsigned maxAllowableId, 
                          std::vector<unsigned>& existingIdsLocal,
                          std::vector<unsigned>&  newIdsLocal,
                          stk::ParallelMachine comm) {
  //
  //  Generate a global list of all ids created matched and all ids that existed perviously
  //
  std::vector<unsigned> existingIdsGlobal;
  stk::parallel_vector_concat(comm, existingIdsLocal, existingIdsGlobal);
  std::sort(existingIdsGlobal.begin(), existingIdsGlobal.end());

  std::vector<unsigned> newIdsGlobal;
  stk::parallel_vector_concat(comm, newIdsLocal, newIdsGlobal);

     
  //
  //  Verify properties of the returned ids
  //
  for(unsigned i=0; i<newIdsGlobal.size(); ++i) {
    unsigned curNewID = newIdsGlobal[i];
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
  unsigned maxAllowableId = ~0U;

  {
    std::vector<unsigned> newIds;
    std::vector<unsigned> existingIds(1000);
    unsigned firstId = mpi_rank*1000;

    for(unsigned i=0; i<1000; ++i) {
      unsigned currentElemId = firstId+i;
      existingIds[i] = currentElemId;
    }
    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, 10, MPI_COMM_WORLD);

    unsigned localSize = newIds.size();
    unsigned globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (unsigned)mpi_size*10);
 
    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }
  //
  //  Test 2, same total source data as case 1, but source data scattered among processors
  //
  {
    int numRequest=0;
    std::vector<unsigned> newIds;
    std::vector<unsigned> existingIds;
    unsigned destProc = mpi_size/2;    
    for(int i=0; i<1000*mpi_size; ++i) {
      destProc += i;
      if(destProc >= (unsigned)mpi_size) {
        destProc = destProc%mpi_size;
      }
      if(destProc == (unsigned)mpi_rank) {
        existingIds.push_back(i);
        if(i%100 == 0) {
          numRequest++;
        }
      }
    }
    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, numRequest, MPI_COMM_WORLD);

    unsigned localSize = newIds.size();
    unsigned globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (unsigned)mpi_size*10);

    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }
  //
  //  Test 3, sparse fill
  //

  {
    std::vector<unsigned> newIds;
    std::vector<unsigned> existingIds(1000);
    unsigned firstId = mpi_rank*1000;

    for(unsigned i=0; i<1000; ++i) {
      unsigned currentElemId = firstId+i;
      existingIds[i] = currentElemId;
    }

    if(mpi_rank == 0) {
      existingIds.push_back(~0U);
    }

    newIds = stk::generate_parallel_unique_ids(maxAllowableId, existingIds, 10, MPI_COMM_WORLD);

    unsigned localSize = newIds.size();
    unsigned globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (unsigned)mpi_size*10);
 
    VerifyUniqueIds(maxAllowableId, existingIds, newIds, MPI_COMM_WORLD);
  }


}
