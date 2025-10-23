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

#ifndef STK_MESH_NGPPARALLELCOMM_HPP
#define STK_MESH_NGPPARALLELCOMM_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/MPI.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/NgpParallelDataExchange.hpp>
#include <stk_mesh/baseImpl/DoOp.hpp>
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {

using CommProcsViewType = Kokkos::View<int*, stk::ngp::MemSpace>;

using OffsetViewType = Kokkos::View<unsigned*, stk::ngp::MemSpace>;

template<typename T, Operation OP, typename ExchangeHandler>
void ngp_parallel_data_exchange_sym_pack_unpack(MPI_Comm mpi_communicator,
                                                ExchangeHandler & exchangeHandler,
                                                bool includeGhosts,
                                                bool deterministic)
{
#if defined( STK_HAS_MPI)
  const int msgTag = 10242;

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup");
  auto& ngpMesh = exchangeHandler.get_ngp_mesh();
  auto ngpFields = exchangeHandler.get_ngp_fields();
  auto& bulkData = ngpMesh.get_bulk_on_host();
  STK_ThrowRequireMsg((std::is_same_v<T,typename ExchangeHandler::T>), "NGP parallel_data_exchange requires that exchangeHandler is templated on the same data-type.");

  const std::vector<stk::mesh::EntityRank>& fieldRanks = exchangeHandler.get_field_ranks();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - max map extent");

  std::vector<int> comm_procs;
  size_t totalMeshIndicesOffsets = 0;

  for (int proc = 0; proc < bulkData.parallel_size(); ++proc) {
    for(stk::mesh::EntityRank fieldRank : fieldRanks) {
      auto sharedCommMap = bulkData.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(fieldRank, proc, includeGhosts);
      if (sharedCommMap.extent(0) > 0) {
        comm_procs.push_back(proc);
        totalMeshIndicesOffsets += sharedCommMap.extent(0);
      }
    }
  }
  stk::util::sort_and_unique(comm_procs);

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - all message sizing");
  size_t totalMsgSizeForAllProcs = 0;
  const size_t num_comm_procs = comm_procs.size();
  std::vector<size_t> messageSizes(num_comm_procs, 0);
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int iproc = comm_procs[proc];
    exchangeHandler.hostSizeMessages(iproc, messageSizes[proc], includeGhosts);
    totalMsgSizeForAllProcs += messageSizes[proc];
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - allocation");
  auto hostBufferOffsets = exchangeHandler.get_host_buffer_offsets();
  Kokkos::resize(Kokkos::WithoutInitializing, hostBufferOffsets, messageSizes.size()+1);

  exchangeHandler.resize_device_mpi_buffers(totalMsgSizeForAllProcs);
  auto& deviceSendData = exchangeHandler.get_device_send_data();
  auto& deviceRecvData = exchangeHandler.get_device_recv_data();

  auto& deviceMeshIndicesOffsets = exchangeHandler.get_device_mesh_indices_offsets();
  auto& hostMeshIndicesOffsets = exchangeHandler.get_host_mesh_indices_offsets();
  //we don't want to just resize these every time, because it turns out that
  //Kokkos::resize is expensive even if the size is staying the same or decreasing.
  if (deviceMeshIndicesOffsets.extent(0) < (totalMeshIndicesOffsets+num_comm_procs) ||
      deviceMeshIndicesOffsets.extent(0) != hostMeshIndicesOffsets.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshIndicesOffsets, totalMeshIndicesOffsets+num_comm_procs);
  }

  if (hostMeshIndicesOffsets.extent(0) < (totalMeshIndicesOffsets+num_comm_procs) ||
      deviceMeshIndicesOffsets.extent(0) != hostMeshIndicesOffsets.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, hostMeshIndicesOffsets, totalMeshIndicesOffsets+num_comm_procs);
  }

  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - offset inits");
  hostBufferOffsets[0] = 0;
  size_t hostMeshIndicesIdx = num_comm_procs;
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    hostBufferOffsets[proc+1] = hostBufferOffsets[proc] + messageSizes[proc];

    hostMeshIndicesOffsets(proc) = hostMeshIndicesIdx;
    unsigned baseProcOffset = hostMeshIndicesIdx;

    unsigned hostMeshIndicesCounter = 0;
    unsigned hostMeshIndicesOffsetsCounter = 0;

    for(stk::mesh::EntityRank fieldRank : fieldRanks) {
      auto sharedCommMap = bulkData.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(fieldRank, comm_procs[proc], includeGhosts);
      hostMeshIndicesIdx += sharedCommMap.extent(0);

      for (size_t i = 0; i < sharedCommMap.extent(0); ++i) {
        hostMeshIndicesOffsets(baseProcOffset+hostMeshIndicesCounter+i) = hostMeshIndicesOffsetsCounter;
        auto meshIndex = sharedCommMap(i);

        for (auto ngpField : ngpFields) {
          if (ngpField->get_rank() == fieldRank) {
            size_t numComponents = stk::mesh::field_scalars_per_entity(*(ngpField->get_field_base()), meshIndex.bucket_id);

            if (numComponents != 0) {
              hostMeshIndicesOffsetsCounter += numComponents;
            }
          }
        }
      }

      hostMeshIndicesCounter += sharedCommMap.extent(0);
    }
  }

  Kokkos::deep_copy(deviceMeshIndicesOffsets, hostMeshIndicesOffsets);
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  std::vector<MPI_Request> sendRequests(num_comm_procs);
  std::vector<MPI_Request> recvRequests(num_comm_procs);
  std::vector<MPI_Status> statuses(num_comm_procs);

  using NgpFieldType = typename ExchangeHandler::field_type;
  using ThisExecSpace = typename ExchangeHandler::mesh_type::MeshExecSpace;
  using ThisRangePolicy = stk::ngp::RangePolicy<ThisExecSpace>;

  Kokkos::Profiling::pushRegion("NGP MPI - message pack");
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    auto iProc = comm_procs[proc];
    auto dataBegin = hostBufferOffsets[proc];
    auto baseProcOffset = hostMeshIndicesOffsets(proc);
    size_t meshIndicesCounter = 0;

    for(stk::mesh::EntityRank fieldRank : fieldRanks) {
      auto hostSharedCommMap = bulkData.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(fieldRank, iProc, includeGhosts);

      Kokkos::parallel_for(ThisRangePolicy(0, hostSharedCommMap.extent(0)),
        KOKKOS_LAMBDA(size_t idx) {
          auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iProc, includeGhosts);
          auto fastMeshIndex = deviceSharedCommMap(idx);

          int sendBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset+meshIndicesCounter+idx);
          const auto& ngpFieldsOnDevice = exchangeHandler.get_ngp_fields_on_device();

          for (size_t fieldIdx = 0; fieldIdx < ngpFieldsOnDevice.extent(0); ++fieldIdx) {
            NgpFieldType const& field = ngpFieldsOnDevice(fieldIdx);
            if (field.get_rank() == fieldRank) {
              stk::mesh::EntityFieldData<T> fieldData = field(fastMeshIndex);
              unsigned numComponents = fieldData.size();
              for (unsigned comp = 0; comp < numComponents; ++comp) {
                deviceSendData(dataBegin + sendBufferStartIdx++) = fieldData[comp];
              }
            }
          }
        }
      );
      Kokkos::fence();
      meshIndicesCounter += hostSharedCommMap.extent(0);
    }
  }
  Kokkos::fence();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI - message send/recv (non-blocking)");
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int iproc = comm_procs[proc];
    const size_t dataBegin = hostBufferOffsets[proc];
    const size_t dataEnd   = hostBufferOffsets[proc+1];
    int bufSize = (dataEnd-dataBegin);
    MPI_Irecv((deviceRecvData.data()+dataBegin), bufSize, sierra::MPI::Datatype<T>::type(), iproc, msgTag, mpi_communicator, &recvRequests[proc]);
    MPI_Isend((deviceSendData.data()+dataBegin), bufSize, sierra::MPI::Datatype<T>::type(), iproc, msgTag, mpi_communicator, &sendRequests[proc]);
  }
  Kokkos::Profiling::popRegion();

  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    Kokkos::Profiling::pushRegion("NGP MPI - message waits");
    int idx = static_cast<int>(proc);
    if (deterministic) {
      MPI_Wait(&recvRequests[proc], MPI_STATUS_IGNORE);
    }
    else {
      MPI_Waitany(static_cast<int>(num_comm_procs), recvRequests.data(), &idx, MPI_STATUS_IGNORE);
    }
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("NGP MPI - message unpacking");
    impl::DoOp<T,OP> doOperation;
    auto iProc = comm_procs[idx];
    auto dataBegin = hostBufferOffsets[idx];
    auto baseProcOffset = hostMeshIndicesOffsets(idx);
    size_t meshIndicesCounter = 0;

    for(stk::mesh::EntityRank fieldRank : fieldRanks) {
      auto hostSharedCommMap = bulkData.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(fieldRank, iProc, includeGhosts);

      Kokkos::parallel_for(ThisRangePolicy(0, hostSharedCommMap.extent(0)),
        KOKKOS_LAMBDA(size_t index) {
          auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iProc, includeGhosts);
          auto fastMeshIndex = deviceSharedCommMap(index);
          int recvBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset+meshIndicesCounter+index);

          const auto& ngpFieldsOnDevice = exchangeHandler.get_ngp_fields_on_device();
          for (size_t fieldIdx = 0; fieldIdx < ngpFieldsOnDevice.extent(0); ++fieldIdx) {
            NgpFieldType const& field = ngpFieldsOnDevice(fieldIdx);
            if (field.get_rank() == fieldRank) {
              stk::mesh::EntityFieldData<T> fieldData = field(fastMeshIndex);
              unsigned numComponents = fieldData.size();

              for (unsigned comp = 0; comp < numComponents; ++comp) {
                auto rcv = deviceRecvData(dataBegin + recvBufferStartIdx++);
                fieldData[comp] = doOperation(fieldData[comp], rcv);
              }
            }
          }
        }
      );
      Kokkos::fence();
      meshIndicesCounter += hostSharedCommMap.extent(0);
    }
    Kokkos::Profiling::popRegion();
  }

  MPI_Waitall(static_cast<int>(num_comm_procs), sendRequests.data(), statuses.data());
  Kokkos::Profiling::popRegion();
#endif
}

}
}

#endif
