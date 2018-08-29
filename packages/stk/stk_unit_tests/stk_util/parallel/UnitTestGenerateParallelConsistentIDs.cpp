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
#include <stk_util/util/ReportHandler.hpp>
#include <stdint.h>


//-------------------------------------------------------------------------------------------------------------------
// Helper functions for test error checking, consolidate pairs of local order and resultant ids into a single
// global sorted list.

std::vector<std::pair<uint64_t, uint64_t> > GenerateGlobalPairedList(std::vector<uint64_t>& localOrderArray, std::vector<uint64_t>&  newIdsLocal, stk::ParallelMachine comm) {

  ThrowRequireMsg(localOrderArray.size() == newIdsLocal.size(),"Inconsistent Sizes in GenerateGlobalPairedList");

  std::vector<std::pair<uint64_t, uint64_t> > idOrderPairsLocal;
  std::vector<std::pair<uint64_t, uint64_t> > idOrderPairsGlobal;
  for(uint64_t i = 0; i < localOrderArray.size(); ++i) {
    idOrderPairsLocal.emplace_back(localOrderArray[i], newIdsLocal[i]);
  }

  stk::parallel_vector_concat(comm, idOrderPairsLocal, idOrderPairsGlobal);

  std::sort(idOrderPairsGlobal.begin(), idOrderPairsGlobal.end());

  return idOrderPairsGlobal;
}

//-------------------------------------------------------------------------------------------------------------------
// Verify that a result generated from the consistent ids function meets all properties nessecary of consistent ids
void  VerifyConsistentIds(uint64_t maxAllowableId, 
                          std::vector<uint64_t>& existingIdsLocal,
                          std::vector<uint64_t>& localOrderArray,
                          std::vector<uint64_t>&  newIdsLocal,
                          stk::ParallelMachine comm) {

  EXPECT_EQ(localOrderArray.size(), newIdsLocal.size());
  //
  //  Generate a global list of all ids created matched with and sorted by their order array entry
  //
  std::vector<std::pair<uint64_t, uint64_t> > idOrderPairsGlobal = GenerateGlobalPairedList(localOrderArray, newIdsLocal, comm);

  std::vector<uint64_t> existingIdsGlobal;
  stk::parallel_vector_concat(comm, existingIdsLocal, existingIdsGlobal);
  std::sort(existingIdsGlobal.begin(), existingIdsGlobal.end());
  //
  //  Verify properties of the returned ids
  //
  for(uint64_t i=0; i<idOrderPairsGlobal.size(); ++i) {
    uint64_t curNewID = idOrderPairsGlobal[i].second;
    //    1)  Verify all ids must be less than or equal to maxAllowableId
    EXPECT_LE(curNewID, maxAllowableId); 
    //    2)  Verify no overlap between existing and new id set
    EXPECT_FALSE(std::binary_search(existingIdsGlobal.begin(), existingIdsGlobal.end(), curNewID));
    //    3)  Verify ids come back in order dictated by the ordering array
    if(i > 0) {
      uint64_t oldNewID = idOrderPairsGlobal[i-1].second;
      EXPECT_LT(oldNewID, curNewID);
      EXPECT_LT(idOrderPairsGlobal[i-1].first, idOrderPairsGlobal[i].first);
    }
  }

}

//------------------------------------------------------------------------------
//
//  Test the generate_parallel_consistent_ids function in stk_util/parallel
//
TEST(UnitTestParallel, testGenerateParallelConsistentIDs) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  //
  //  Test 1, nominal element death case
  //    Each processor starts with 1000 unique densly packed ids ordered 1->N.  Add 10 new ids per processor
  //    produced by kill every 100th element creating 10 new elements per processor.
  //
  uint64_t maxAllowableId = ~0U;

  std::vector<uint64_t> newIds1;
  std::vector<uint64_t> localOrderArray1;

  {
    std::vector<uint64_t> existingIds(1000);
    uint64_t firstId = mpi_rank*1000;


    for(uint64_t i=0; i<1000; ++i) {
      uint64_t currentElemId = firstId+i;
      existingIds[i] = currentElemId;
      if(i%100 == 0) {
        //  local order array is the 'source element' which provides a unique repeatable ordering
        localOrderArray1.push_back(currentElemId);
      }
    }
    newIds1 = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray1, MPI_COMM_WORLD);

    uint64_t localSize = newIds1.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);
    VerifyConsistentIds(maxAllowableId, existingIds, localOrderArray1, newIds1, MPI_COMM_WORLD);
  }
  //
  //  Test 2, same total source data as case 1, but source data scattered among processors
  //
  std::vector<uint64_t> newIds2;
  std::vector<uint64_t> localOrderArray2;
  {
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
          localOrderArray2.push_back(i);
        }
      }
    }
    newIds2 = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray2, MPI_COMM_WORLD);

    uint64_t localSize = newIds2.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);
    VerifyConsistentIds(maxAllowableId, existingIds, localOrderArray2, newIds2, MPI_COMM_WORLD);
  }
  //
  //  Test 3, verify test 2 and test 3 generate the same ids in identical order, confirming the id generation scheme is decomposition independent
  //
  std::vector<std::pair<uint64_t, uint64_t> > idOrderPairsGlobal1 = GenerateGlobalPairedList(localOrderArray1, newIds1, MPI_COMM_WORLD);
  std::vector<std::pair<uint64_t, uint64_t> > idOrderPairsGlobal2 = GenerateGlobalPairedList(localOrderArray2, newIds2, MPI_COMM_WORLD);
  EXPECT_EQ(idOrderPairsGlobal1.size(), idOrderPairsGlobal2.size());
  for(uint64_t i=0; i < idOrderPairsGlobal1.size(); ++i) {
    EXPECT_EQ(idOrderPairsGlobal1[i].first,  idOrderPairsGlobal2[i].first);
    EXPECT_EQ(idOrderPairsGlobal1[i].second, idOrderPairsGlobal2[i].second);
  }
  //
  //  Test 4, sparse fill
  //
  {
    std::vector<uint64_t> newIds;
    std::vector<uint64_t> existingIds(1000);
    std::vector<uint64_t> localOrderArray;

    uint64_t firstId = mpi_rank*1000;

    for(uint64_t i=0; i<1000; ++i) {
      uint64_t currentElemId = firstId+i;
      existingIds[i] = currentElemId;
    }

    if(mpi_rank == 0) {
      existingIds.push_back(~0U);
    }


    for(int i=0; i<10; ++i) {
      localOrderArray.push_back(mpi_rank*10 + i);
    }


    newIds = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray, MPI_COMM_WORLD);
  

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)mpi_size*10);
    VerifyConsistentIds(maxAllowableId, existingIds, localOrderArray, newIds, MPI_COMM_WORLD);
  } 

}
