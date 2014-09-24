// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string


class TestStruct {
public:
  unsigned int data1;
  double data2;
  double data3;
};


TEST(UnitTestParallel, testParallelVectorConcat) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  //
  //  Test 1, concat integer lists, each process P creates n-1 copies of P thus on 6 processors the 
  //  expected concat list is:
  //
  //  2 3 3 4 4 4 5 5 5 5
  //
  std::vector<int> localIntList;
  for(int i=0; i<mpi_rank-1; ++i) {
    localIntList.push_back(mpi_rank);
  }
  std::vector<int> globalIntList;
  int status = stk::parallel_vector_concat(MPI_COMM_WORLD, localIntList, globalIntList);
  EXPECT_EQ(status, MPI_SUCCESS);
  std::vector<int> expectedIntList;
  for(int iproc=0; iproc<mpi_size; ++iproc) {
    for(int i=0; i<iproc-1; ++i) {
      expectedIntList.push_back(iproc);
    }
  }
  for(unsigned int i=0; i<expectedIntList.size(); ++i) {
    EXPECT_EQ(globalIntList[i], expectedIntList[i]);
  }
  //
  //  Test 2, concat a bit more complex data struct with two reals and an int, same general format
  //  as before
  //
  std::vector<TestStruct> localStructList;
  std::vector<TestStruct> expectedStructList;
  for(int iproc=0; iproc<mpi_size; ++iproc) {
    for(int i=0; i<mpi_size; ++i) {
      TestStruct ts;
      ts.data1 = iproc;
      ts.data2 = 1.23;
      ts.data3 = 1.0/i;
      if(iproc == mpi_rank) {
        localStructList.push_back(ts);
      }
      expectedStructList.push_back(ts);
    }
  }
  std::vector<TestStruct> globalStructList;
  status = stk::parallel_vector_concat(MPI_COMM_WORLD, localStructList, globalStructList);
  EXPECT_EQ(status, MPI_SUCCESS);
  for(unsigned int i=0; i<expectedStructList.size(); ++i) {
    EXPECT_EQ(globalStructList[i].data1, expectedStructList[i].data1);
    EXPECT_EQ(globalStructList[i].data2, expectedStructList[i].data2);
    EXPECT_EQ(globalStructList[i].data3, expectedStructList[i].data3);
  }
  //
  //  Test 3, boundary case, check for handling of empty list.
  //
  std::vector<char> localCharList;
  std::vector<char> globalCharList(10, 'a');
  status = stk::parallel_vector_concat(MPI_COMM_WORLD, localCharList, globalCharList);
  EXPECT_EQ(status, MPI_SUCCESS);
  EXPECT_EQ(globalCharList.size(), (unsigned)0);

}
