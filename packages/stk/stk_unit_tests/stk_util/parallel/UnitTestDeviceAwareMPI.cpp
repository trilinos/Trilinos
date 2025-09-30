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
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/DeviceAwareMPI.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

#ifdef STK_HAS_MPI

TEST(DeviceAwareMPI, trueIfOpenMPIAndCuda)
{
#if defined(OMPI_MAJOR_VERSION) && defined(KOKKOS_ENABLE_CUDA)
  EXPECT_TRUE(stk::have_device_aware_mpi());
#else
  GTEST_SKIP()<<"trueIfOpenMPIAndCuda";
#endif
}

TEST(DeviceAwareMPI, trueIfATS2MPIAndCuda)
{
#if defined(IBM_SPECTRUM_MPI) && defined(KOKKOS_ENABLE_CUDA)
  EXPECT_TRUE(stk::have_device_aware_mpi());
#else
  GTEST_SKIP()<<"trueIfATS2MPIAndCuda";
#endif
}

TEST(DeviceAwareMPI, falseIfOpenMpiButNoCuda)
{
#if defined(OMPI_MAJOR_VERSION) && !defined(KOKKOS_ENABLE_CUDA)
  EXPECT_FALSE(stk::have_device_aware_mpi());
#else
  GTEST_SKIP()<<"falseIfOpenMpiButNoCuda";
#endif
}

TEST(DeviceAwareMPI, falseIfIntelMpi)
{
#if defined(I_MPI_VERSION)
  EXPECT_FALSE(stk::have_device_aware_mpi());
#else
  GTEST_SKIP()<<"falseIfIntel, I_MPI_VERSION not defined";
#endif
}

void check_device_aware_mpi_send_recv(MPI_Comm comm)
{
  const int numProcs = stk::parallel_machine_size(comm);
  ASSERT_EQ(2, numProcs);
  const int myProc = stk::parallel_machine_rank(comm);
  const int otherProc = 1 - myProc;
  const int msgTag = 10101;

  using BufferViewType = Kokkos::View<double*,stk::ngp::MemSpace>;
  constexpr size_t N = 8;
  constexpr double tol = 1.e-7;
  constexpr double goldValue = 3.14159;

  if (myProc == 0) {
    BufferViewType sendBuf("sendBuf",N);
    Kokkos::deep_copy(sendBuf, goldValue);
    EXPECT_EQ(MPI_SUCCESS, MPI_Send(sendBuf.data(), N, MPI_DOUBLE, otherProc, msgTag, comm));
  }
  else {
    BufferViewType recvBuf("recvBuf",N);
    Kokkos::deep_copy(recvBuf, 0.0);
    MPI_Status status;
    EXPECT_EQ( MPI_SUCCESS, MPI_Recv(recvBuf.data(), N, MPI_DOUBLE, otherProc, msgTag, comm, &status));

    BufferViewType::host_mirror_type hostRecvBuf = Kokkos::create_mirror_view(recvBuf);
    Kokkos::deep_copy(hostRecvBuf, recvBuf);
    for(size_t i=0; i<N; ++i) {
      EXPECT_NEAR(goldValue, hostRecvBuf(i), tol);
    }
  }
}

TEST(DeviceAwareMPI, basicSendRecv_np2)
{
  if (stk::have_device_aware_mpi() && stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    check_device_aware_mpi_send_recv(MPI_COMM_WORLD);
  }
  else {
    GTEST_SKIP();
  }
}

void check_device_aware_mpi_isend_irecv_waitany(MPI_Comm comm)
{
  const int numProcs = stk::parallel_machine_size(comm);
  ASSERT_EQ(2, numProcs);
  const int myProc = stk::parallel_machine_rank(comm);
  const int otherProc = 1 - myProc;
  const int msgTag = 10101;

  using BufferViewType = Kokkos::View<double*,stk::ngp::MemSpace>;
  constexpr size_t N = 4200;
  constexpr double tol = 1.e-7;
  constexpr double goldValue = 3.14159;

  const int numCommProcs = 1;
  std::vector<MPI_Request> sendRequests(numCommProcs);
  std::vector<MPI_Request> recvRequests(numCommProcs);
  std::vector<MPI_Status> statuses(numCommProcs);

  BufferViewType sendBuf("sendBuf",N);
  BufferViewType recvBuf("recvBuf",N);

  stk::parallel_machine_barrier(comm);

  if (myProc == 1) {
    Kokkos::deep_copy(recvBuf, 0.0); //theoretically unnecessary since Views initialize by default.
    EXPECT_EQ(MPI_SUCCESS, MPI_Irecv(recvBuf.data(), N, MPI_DOUBLE, otherProc, msgTag, comm, &recvRequests[0]));
  }
  if (myProc == 0) {
    Kokkos::deep_copy(sendBuf, goldValue);
    EXPECT_EQ(MPI_SUCCESS, MPI_Isend(sendBuf.data(), N, MPI_DOUBLE, otherProc, msgTag, comm, &sendRequests[0]));
  }

  Kokkos::fence();

  if (myProc == 1) {
    int idx = 99;
    MPI_Waitany(numCommProcs, recvRequests.data(), &idx, MPI_STATUS_IGNORE);
    EXPECT_EQ(0, idx);

    BufferViewType::host_mirror_type hostRecvBuf = Kokkos::create_mirror_view(recvBuf);
    Kokkos::deep_copy(hostRecvBuf, recvBuf);
    for(size_t i=0; i<N; ++i) {
      EXPECT_NEAR(goldValue, hostRecvBuf(i), tol);
    }
  }

  if (myProc == 0) {
    MPI_Waitall(numCommProcs, sendRequests.data(), statuses.data());
  }
}

TEST(DeviceAwareMPI, basicIsendIrecv_np2)
{
  if (stk::have_device_aware_mpi() && stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    check_device_aware_mpi_isend_irecv_waitany(MPI_COMM_WORLD);
  }
  else {
    GTEST_SKIP();
  }
}

#endif

