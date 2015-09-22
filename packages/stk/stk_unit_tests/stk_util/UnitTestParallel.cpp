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
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll, CommBuffer
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_write_string
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string
 
TEST(UnitTestParallel, testUnit)
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  std::string s;
  std::ostringstream strout;

  std::cout << "all_write_string " << std::flush;
  
//  for (size_t i = 0; i < 250000; ++i) {
  for (size_t i = 0; i < 100; ++i) {
    if (mpi_rank == 0 && i%1000 == 0)
      std::cout << "." << std::flush;
    
    stk::all_write_string(MPI_COMM_WORLD, strout, s);
  }
  
  ASSERT_LT(mpi_rank, mpi_size);
}

TEST(UnitTestParallel, testCommAll)
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  ASSERT_LT(mpi_rank, mpi_size);

  bool propagate_local_error_flags = true;
  stk::CommAll comm_all(MPI_COMM_WORLD, propagate_local_error_flags);

  ASSERT_EQ(comm_all.parallel_size(), mpi_size);
  ASSERT_EQ(comm_all.parallel_rank(), mpi_rank);

  for (int p = 0; p<mpi_size; ++p) {
    if (p - mpi_rank <= 10 || mpi_rank - p <= 10) {
      stk::CommBuffer& buf = comm_all.send_buffer(p);
      buf.skip<unsigned long>(1);
    }
  }

  bool symmetric = false;
  bool local_error = false;
  bool global_bad_input = comm_all.allocate_buffers(mpi_size / 4, symmetric, local_error);
  ASSERT_EQ(global_bad_input, false);
}

