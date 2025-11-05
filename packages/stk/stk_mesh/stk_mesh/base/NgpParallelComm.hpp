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
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpParallelDataExchange.hpp>
#include <stk_mesh/baseImpl/DoOp.hpp>
#include <Kokkos_Core.hpp>
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/MPITag.hpp"

namespace stk::mesh {

using CommProcsViewType = Kokkos::View<int*, stk::ngp::MemSpace>;

using OffsetViewType = Kokkos::View<unsigned*, stk::ngp::MemSpace>;

template<typename Scalar, Operation OP, typename NgpSpace>
void ngp_parallel_data_excahnge_sym_pack_unpack(MPI_Comm mpi_communicator,
                                                const BulkData & mesh,
                                                const std::vector<const FieldBase*>& fields,
                                                bool includeGhosts,
                                                bool deterministic)
{
#if defined( STK_HAS_MPI)
  auto msgTag = 0;
  auto mpiTag = stk::MPITag();
  if (stk::util::get_common_coupling_version() < 19) {
    msgTag = 10242;
  } else {
    mpiTag = get_mpi_tag_manager().get_tag(mpi_communicator, 10242);
    msgTag = mpiTag;
  }

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup");

  auto& ngpMesh = get_updated_ngp_mesh(mesh);

  const auto fieldRanks = impl::assemble_rank_list(fields);
  const auto rankPerField = impl::assemble_rank_per_field<NgpSpace>(fields);
  const auto fieldDataOnDevice = impl::assemble_field_data_on_device<Scalar, NgpSpace>(fields);

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - max map extent");

  const auto comm_procs = impl::assemble_comm_procs_list<NgpSpace>(mesh, fieldRanks, includeGhosts);
  const auto num_comm_procs = comm_procs.size();
  const auto totalMeshIndicesOffsets = impl::compute_total_mesh_indices_offsets<NgpSpace>(mesh, fieldRanks, includeGhosts);

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - all message sizing");

  const auto messageSizes = impl::compute_message_sizes<NgpSpace>(mesh, fields, comm_procs, fieldRanks, includeGhosts);
  const auto totalMsgSizeForAllProcs = std::accumulate(std::cbegin(messageSizes), std::cend(messageSizes), 0UL);

  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - allocation");

  auto& hostBufferOffsets = ngpMesh.get_ngp_parallel_sum_host_buffer_offsets();
  Kokkos::resize(Kokkos::WithoutInitializing, hostBufferOffsets, messageSizes.size() + 1);

  using BufferView = Kokkos::View<Scalar*, typename NgpSpace::mem_space>;
  auto deviceSendData = BufferView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceSendData"), totalMsgSizeForAllProcs);
  auto deviceRecvData = BufferView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceRecvData"), totalMsgSizeForAllProcs);

  auto& deviceMeshIndicesOffsets = ngpMesh.get_ngp_parallel_sum_device_mesh_indices_offsets();
  auto& hostMeshIndicesOffsets = ngpMesh.get_ngp_parallel_sum_host_mesh_indices_offsets();

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

  impl::fill_host_buffer_offsets<NgpSpace>(hostBufferOffsets,
                                           hostMeshIndicesOffsets,
                                           mesh,
                                           fields,
                                           comm_procs,
                                           fieldRanks,
                                           messageSizes,
                                           includeGhosts);
  Kokkos::deep_copy(deviceMeshIndicesOffsets, hostMeshIndicesOffsets);

  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  std::vector<MPI_Request> sendRequests(num_comm_procs);
  std::vector<MPI_Request> recvRequests(num_comm_procs);
  std::vector<MPI_Status> statuses(num_comm_procs);

  Kokkos::Profiling::pushRegion("NGP MPI - message pack");

  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    auto iproc = comm_procs[proc];
    auto dataBegin = hostBufferOffsets[proc];
    auto baseProcOffset = hostMeshIndicesOffsets(proc);

    impl::fill_device_send_data<NgpSpace>(deviceSendData,
                                          fieldDataOnDevice,
                                          rankPerField,
                                          deviceMeshIndicesOffsets,
                                          ngpMesh,
                                          mesh,
                                          fieldRanks,
                                          iproc,
                                          dataBegin,
                                          baseProcOffset,
                                          includeGhosts);
  }
  Kokkos::fence();

  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI - message send/recv (non-blocking)");

  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    int iproc = comm_procs[proc];
    const size_t dataBegin = hostBufferOffsets[proc];
    const size_t dataEnd = hostBufferOffsets[proc + 1];
    int bufSize = (dataEnd-dataBegin);
    MPI_Irecv((deviceRecvData.data() + dataBegin), bufSize, sierra::MPI::Datatype<Scalar>::type(), iproc, msgTag, mpi_communicator, &recvRequests[proc]);
    MPI_Isend((deviceSendData.data() + dataBegin), bufSize, sierra::MPI::Datatype<Scalar>::type(), iproc, msgTag, mpi_communicator, &sendRequests[proc]);
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

    impl::DoOp<Scalar, OP> doOperation;
    auto iproc = comm_procs[idx];
    auto dataBegin = hostBufferOffsets[idx];
    auto baseProcOffset = hostMeshIndicesOffsets(idx);

    impl::unpack_device_recv_data<NgpSpace>(fieldDataOnDevice,
                                            deviceRecvData,
                                            rankPerField,
                                            deviceMeshIndicesOffsets,
                                            ngpMesh,
                                            mesh,
                                            fieldRanks,
                                            iproc,
                                            dataBegin,
                                            baseProcOffset,
                                            doOperation,
                                            includeGhosts);

    Kokkos::Profiling::popRegion();
  }

  MPI_Waitall(static_cast<int>(num_comm_procs), sendRequests.data(), statuses.data());

  Kokkos::Profiling::popRegion();
#endif
}
} // namespace stk::mesh

#endif
