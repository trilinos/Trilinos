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

#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll
#include <stk_util/parallel/MPI.hpp>
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include <stk_util/stk_config.h>
#include <limits>

#if defined ( STK_HAS_MPI )

TEST(ParallelComm, CommAllDestructor)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD ;
    int mpi_size = stk::parallel_machine_size(comm);
    const unsigned zero = 0 ;
    std::vector<unsigned> msg_size( mpi_size , zero );
    const unsigned * const s_size = & msg_size[0] ;
    {
        stk::CommAll sparse;
        // This will allocate both m_send and m_recv regardless of mpi_size
        sparse.allocate_symmetric_buffers( comm, s_size );
        stk::CommBuffer * recv_buffer = &sparse.recv_buffer(0);
        stk::CommBuffer * send_buffer = &sparse.send_buffer(0);
        ASSERT_TRUE( recv_buffer != send_buffer );
    }
    // This should not produce a memory leak.
}

TEST(ParallelComm, CommunicateMPILocInt)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);
    
    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    int nvalues=5;
    int64_t limit_32bit_integer = std::numeric_limits<int32_t>::max();

    sierra::MPI::Loc<int> * const vin  = new sierra::MPI::Loc<int>[nvalues] ;
    sierra::MPI::Loc<int> * const vout = new sierra::MPI::Loc<int>[nvalues] ;

    for(int n = 0; n < nvalues; n++) {
      vin[n].m_value = myProcId*numProcs+n;
      vin[n].m_loc = limit_32bit_integer + myProcId*numProcs*numProcs+n; //Want to test when outside 32-bit range for index.
    }

    MPI_Allreduce( vin, vout, nvalues,
               sierra::MPI::Datatype<sierra::MPI::Loc<int> >::type(),
               sierra::MPI::get_mpi_loc_op<int, std::greater<int>, int64_t >(),
               comm);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ((numProcs-1)*numProcs+n, vout[n].m_value);
      EXPECT_EQ(limit_32bit_integer + (numProcs-1)*numProcs*numProcs+n, vout[n].m_loc);
    }
    
    MPI_Allreduce( vin, vout, nvalues,
               sierra::MPI::Datatype<sierra::MPI::Loc<int> >::type(),
               sierra::MPI::get_mpi_loc_op<int, std::less<int>, int64_t>(),
               comm);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n, vout[n].m_value);
      EXPECT_EQ(limit_32bit_integer+n, vout[n].m_loc);
    }

    delete[] vin;
    delete[] vout;
}

TEST(ParallelComm, CommunicateMPILocDouble)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);
    
    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    int nvalues=1;
    double value_offset = 55.5;
    int64_t limit_32bit_integer = std::numeric_limits<int32_t>::max();

    sierra::MPI::Loc<double> * const vin  = new sierra::MPI::Loc<double>[nvalues] ;
    sierra::MPI::Loc<double> * const vout = new sierra::MPI::Loc<double>[nvalues] ;

    for(int n = 0; n < nvalues; n++) {
      vin[n].m_value = value_offset+myProcId*numProcs+n;
      vin[n].m_loc = limit_32bit_integer + myProcId*numProcs*numProcs+n; //Want to test when outside 32-bit range for index.
    }

    MPI_Allreduce( vin, vout, nvalues,
               sierra::MPI::Datatype<sierra::MPI::Loc<double> >::type(),
               sierra::MPI::get_mpi_loc_op<double, std::greater<double>, int64_t >(),
               comm);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(value_offset+double((numProcs-1)*numProcs+n), vout[n].m_value);
      EXPECT_EQ(limit_32bit_integer + (numProcs-1)*numProcs*numProcs+n, vout[n].m_loc);
    }
    
    MPI_Allreduce( vin, vout, nvalues,
               sierra::MPI::Datatype<sierra::MPI::Loc<double> >::type(),
               sierra::MPI::get_mpi_loc_op<double, std::less<double>, int64_t >(),
               comm);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(value_offset+n, vout[n].m_value);
      EXPECT_EQ(limit_32bit_integer+n, vout[n].m_loc);
    }

    delete[] vin;
    delete[] vout;
}

TEST(ParallelComm, comm_sizes)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
  int myProc = stk::parallel_machine_rank(comm);

  unsigned num_msg_bound = numProcs, num_msg_max = 0;
  std::vector<unsigned> send_sizes(numProcs, 0), recv_sizes(numProcs, 0);
  
  for(int p=myProc-2; p<=myProc+2; ++p) {
    if (p >= 0 && p < numProcs) {
      send_sizes[p] = myProc+1;
    }
  }

  stk::comm_sizes(comm, num_msg_bound, num_msg_max, &send_sizes[0], &recv_sizes[0]);

  for(int p=0; p<numProcs; ++p) {
    if (p < (myProc-2) || p > (myProc+2)) {
      EXPECT_EQ(0u, recv_sizes[p]);
    }
    else {
      EXPECT_EQ((unsigned)p+1, recv_sizes[p]);
    }
  }
}

#endif
