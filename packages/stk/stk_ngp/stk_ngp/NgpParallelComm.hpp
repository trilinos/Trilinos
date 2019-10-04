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

#ifndef NGPPARALLELCOMM_HPP
#define NGPPARALLELCOMM_HPP

#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_ngp/NgpSpaces.hpp"
#include "Kokkos_Core.hpp"

namespace ngp {

using CommProcsViewType = Kokkos::View<int*, MemSpace>;

using OffsetViewType = Kokkos::View<unsigned*, MemSpace>;

template <typename T>
using BufferViewType = Kokkos::View<T*, MemSpace>;

template<typename T, typename ExchangeHandler>
void parallel_data_exchange_sym_pack_unpack(MPI_Comm mpi_communicator,
                                            const std::vector<int> & comm_procs,
                                            ExchangeHandler & exchangeHandler,
                                            bool deterministic)
{
#if defined( STK_HAS_MPI)
  const int pRank = stk::parallel_machine_rank(mpi_communicator);
  const int msgTag = 10242;
  size_t num_comm_procs = comm_procs.size();
  int dataTypeSize = sizeof(T);

  CommProcsViewType deviceCommProcs("DeviceCommProcs", num_comm_procs);
  CommProcsViewType::HostMirror hostCommProcs = Kokkos::create_mirror_view(deviceCommProcs);
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    hostCommProcs(proc) = comm_procs[proc];
  }
  Kokkos::deep_copy(deviceCommProcs, hostCommProcs);

  size_t totalSizeForAllProcs = 0;
  std::vector<size_t> messageSizes(num_comm_procs, 0);
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int iproc = comm_procs[proc];
    exchangeHandler.hostSizeMessages(iproc, messageSizes[proc]);
    totalSizeForAllProcs += messageSizes[proc];
  }

  OffsetViewType bufferOffsets = OffsetViewType(Kokkos::ViewAllocateWithoutInitializing("BufferOffsets"), messageSizes.size()+1);
  OffsetViewType::HostMirror hostBufferOffsets = Kokkos::create_mirror_view(bufferOffsets);

  BufferViewType<T> deviceSendData = BufferViewType<T>(Kokkos::ViewAllocateWithoutInitializing("BufferSendData"), totalSizeForAllProcs);
  BufferViewType<T> deviceRecvData = BufferViewType<T>(Kokkos::ViewAllocateWithoutInitializing("BufferRecvData"), totalSizeForAllProcs);

  hostBufferOffsets[0] = 0;
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    hostBufferOffsets[proc+1] = hostBufferOffsets[proc] + messageSizes[proc];
  }
  Kokkos::deep_copy(bufferOffsets, hostBufferOffsets);

  std::vector<MPI_Request> sendRequests(num_comm_procs);
  std::vector<MPI_Request> recvRequests(num_comm_procs);
  std::vector<MPI_Status> statuses(num_comm_procs);

  Kokkos::parallel_for(num_comm_procs, KOKKOS_LAMBDA(size_t iproc)
  {
    const size_t dataBegin = bufferOffsets[iproc];
    const size_t dataEnd   = bufferOffsets[iproc+1];
    BufferViewType<T> buffer =  Kokkos::subview( deviceSendData, Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
    exchangeHandler.devicePackMessage(pRank, deviceCommProcs(iproc), buffer);
  });

  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int iproc = comm_procs[proc];
    const size_t dataBegin = hostBufferOffsets[proc];
    const size_t dataEnd   = hostBufferOffsets[proc+1];
    int bufSize = (dataEnd-dataBegin) * dataTypeSize;
    MPI_Irecv((deviceRecvData.data()+dataBegin), bufSize, MPI_CHAR, iproc, msgTag, mpi_communicator, &recvRequests[proc]);
    MPI_Isend((deviceSendData.data()+dataBegin), bufSize, MPI_CHAR, iproc, msgTag, mpi_communicator, &sendRequests[proc]);
  }

  MPI_Status status;
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int idx = static_cast<int>(proc);
    if (deterministic) {
        MPI_Wait(&recvRequests[proc], &status);
    }
    else {
        MPI_Waitany(static_cast<int>(num_comm_procs), recvRequests.data(), &idx, &status);
    }

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(size_t)
    {
      const size_t dataBegin = bufferOffsets[idx];
      const size_t dataEnd   = bufferOffsets[idx+1];
      BufferViewType<T> buffer =  Kokkos::subview( deviceRecvData, Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
      exchangeHandler.deviceUnpackMessage(pRank, deviceCommProcs[idx], buffer);
    });
  }

  MPI_Waitall(static_cast<int>(num_comm_procs), sendRequests.data(), statuses.data());
#endif
}

}

#endif
