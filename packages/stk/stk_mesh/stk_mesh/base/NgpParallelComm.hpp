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
#include "stk_mesh/base/Types.hpp"
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpParallelDataExchange.hpp>
#include <stk_mesh/baseImpl/DoOp.hpp>
#include <Kokkos_Core.hpp>
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/DeviceAwareMPI.hpp"
#include "stk_util/parallel/MPITag.hpp"

namespace stk::mesh {

template<typename Scalar, Operation OP, typename NgpSpace, typename NgpMeshType, Layout LayoutValue>
void ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh_layout(const BulkData & mesh,
                                                               NgpMeshType& ngpMesh,
                                                               const std::vector<const FieldBase*>& fields,
                                                               bool includeGhosts,
                                                               bool deterministic)
{
#if defined( STK_HAS_MPI)
  Kokkos::Profiling::pushRegion("NGP parallel_data_exchange");
  MPI_Comm mpi_communicator = mesh.parallel();
  auto msgTag = 0;
  auto mpiTag = stk::MPITag();
  if (stk::util::get_common_coupling_version() < 19) {
    msgTag = 10242;
  } else {
    mpiTag = get_mpi_tag_manager().get_tag(mpi_communicator, 10242);
    msgTag = mpiTag;
  }

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup");

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - assemble field ranks, offsets, and data");

  const std::vector<EntityRank> fieldRanks = impl::assemble_rank_list(fields);
  const std::vector<unsigned> rankFieldOffsets = impl::group_fields_by_rank(fields, fieldRanks);
  const auto fieldDataOnDevice = impl::assemble_field_data_on_device<Scalar, NgpSpace, LayoutValue>(fields, fieldRanks);

  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - procs and mesh indices offsets");

  const std::vector<int> comm_procs = impl::assemble_comm_procs_list(mesh, fieldRanks, includeGhosts);
  const size_t num_comm_procs = comm_procs.size();
  const size_t totalMeshIndicesOffsets = impl::compute_total_mesh_indices_offsets(mesh, fieldRanks, comm_procs, includeGhosts);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - offsets allocation");

  auto& hostBufferOffsets = ngpMesh.get_ngp_parallel_sum_host_buffer_offsets();
  if (hostBufferOffsets.extent(0) < (comm_procs.size()+1)) {
    Kokkos::resize(Kokkos::WithoutInitializing, hostBufferOffsets, comm_procs.size() + 1);
  }

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

  impl::fill_host_buffer_offsets(hostBufferOffsets,
                                 hostMeshIndicesOffsets,
                                 mesh,
                                 fields,
                                 comm_procs,
                                 fieldRanks,
                                 includeGhosts);
  Kokkos::deep_copy(deviceMeshIndicesOffsets, hostMeshIndicesOffsets);

  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI bookkeeping setup - buffer allocation");

  std::vector<MPI_Request> sendRequests(num_comm_procs);
  std::vector<MPI_Request> recvRequests(num_comm_procs);
  std::vector<MPI_Status> statuses(num_comm_procs);

  using BufferView = Kokkos::View<Scalar*, typename NgpSpace::mem_space>;
#if KOKKOS_VERSION >= 50000
  using BufferHostView = typename BufferView::host_mirror_type;
#else
  using BufferHostView = typename BufferView::HostMirror;
#endif
  const auto totalMsgSizeForAllProcs = hostBufferOffsets(num_comm_procs);
  auto& sendAndRecvData = ngpMesh.get_ngp_parallel_sum_device_byte_buffer();
  size_t numBytes = (deterministic) ? 3*totalMsgSizeForAllProcs*sizeof(Scalar) : 2*totalMsgSizeForAllProcs*sizeof(Scalar);
  if (sendAndRecvData.extent(0) < numBytes) {
    Kokkos::resize(Kokkos::WithoutInitializing, sendAndRecvData, numBytes);
  }
  Scalar* sendAndRecvDataPtr = reinterpret_cast<Scalar*>(sendAndRecvData.data());
  auto deviceSendData = BufferView(sendAndRecvDataPtr, totalMsgSizeForAllProcs);
  auto deviceRecvData = BufferView(sendAndRecvDataPtr+totalMsgSizeForAllProcs, totalMsgSizeForAllProcs);

  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("NGP MPI - message send/recv (non-blocking)");

  auto launchReceives = [&](auto& recvView)
  {
    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      int iproc = comm_procs[proc];
      const size_t dataBegin = hostBufferOffsets[proc];
      const size_t dataEnd = hostBufferOffsets[proc + 1];
      int bufSize = (dataEnd-dataBegin);
      MPI_Irecv((recvView.data() + dataBegin), bufSize, sierra::MPI::Datatype<Scalar>::type(), iproc, msgTag, mpi_communicator, &recvRequests[proc]);
    }
  };

  auto& hostSendAndRecvData = ngpMesh.get_ngp_parallel_sum_host_byte_buffer();

  const bool useDeviceAwareMpi = use_device_aware_mpi();
  if (!useDeviceAwareMpi) {
    if (hostSendAndRecvData.extent(0) < sendAndRecvData.extent(0)) {
      Kokkos::resize(Kokkos::WithoutInitializing, hostSendAndRecvData, sendAndRecvData.extent(0));
    }
  }
  Scalar* hostSendAndRecvDataPtr = reinterpret_cast<Scalar*>(hostSendAndRecvData.data());
  auto hostSendData = useDeviceAwareMpi ? BufferHostView("deviceSendDataHost", 0) : BufferHostView(hostSendAndRecvDataPtr, totalMsgSizeForAllProcs);
  auto hostRecvData = useDeviceAwareMpi ? BufferHostView("deviceRecvDataHost", 0) : BufferHostView(hostSendAndRecvDataPtr+totalMsgSizeForAllProcs, totalMsgSizeForAllProcs);
  if (useDeviceAwareMpi) {
    Kokkos::Profiling::pushRegion("NGP MPI - message send (device-aware-mpi)");
    launchReceives(deviceRecvData);
    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      int iproc = comm_procs[proc];
      const size_t dataBegin = hostBufferOffsets[proc];
      const size_t dataEnd = hostBufferOffsets[proc + 1];
      const size_t baseProcOffset = hostMeshIndicesOffsets(proc);

      impl::fill_device_send_data<NgpSpace>(deviceSendData,
                                            fieldDataOnDevice,
                                            deviceMeshIndicesOffsets,
                                            ngpMesh,
                                            mesh,
                                            fieldRanks,
                                            rankFieldOffsets,
                                            iproc,
                                            dataBegin,
                                            baseProcOffset,
                                            includeGhosts);
      const int bufSize = (dataEnd-dataBegin);
      Kokkos::fence();
      MPI_Isend((deviceSendData.data() + dataBegin), bufSize, sierra::MPI::Datatype<Scalar>::type(), iproc, msgTag, mpi_communicator, &sendRequests[proc]);
    }
    Kokkos::Profiling::popRegion();
  }
  else {
    Kokkos::Profiling::pushRegion("NGP MPI - message send/recv (non-device-aware-mpi)");
    launchReceives(hostRecvData);
    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      int iproc = comm_procs[proc];
      const size_t dataBegin = hostBufferOffsets[proc];
      const size_t dataEnd = hostBufferOffsets[proc + 1];
      const size_t baseProcOffset = hostMeshIndicesOffsets(proc);

      impl::fill_device_send_data<NgpSpace>(deviceSendData,
                                            fieldDataOnDevice,
                                            deviceMeshIndicesOffsets,
                                            ngpMesh,
                                            mesh,
                                            fieldRanks,
                                            rankFieldOffsets,
                                            iproc,
                                            dataBegin,
                                            baseProcOffset,
                                            includeGhosts);
      const int bufSize = (dataEnd-dataBegin);
      BufferView deviceSendDataRange(deviceSendData.data()+dataBegin,bufSize);
      BufferView hostSendDataRange(hostSendData.data()+dataBegin,bufSize);
      Kokkos::deep_copy(hostSendDataRange, deviceSendDataRange);
      MPI_Isend((hostSendData.data() + dataBegin), bufSize, sierra::MPI::Datatype<Scalar>::type(), iproc, msgTag, mpi_communicator, &sendRequests[proc]);
    }
    Kokkos::Profiling::popRegion();
  }

  Kokkos::Profiling::popRegion();

  if (deterministic) {
    auto selfSendData = BufferView(sendAndRecvDataPtr+2*totalMsgSizeForAllProcs, totalMsgSizeForAllProcs);
    impl::DoOp<Scalar, OP> doOperation;
    Kokkos::deep_copy(selfSendData, doOperation.initial_value());

    size_t my_proc = stk::parallel_machine_rank(mpi_communicator);
    auto my_proc_location = std::lower_bound(std::begin(comm_procs), std::end(comm_procs), my_proc);

    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      int idx = static_cast<int>(proc);
      auto iproc = comm_procs[idx];
      auto dataBegin = hostBufferOffsets[idx];
      auto baseProcOffset = hostMeshIndicesOffsets(idx);
      Kokkos::Profiling::pushRegion("NGP MPI - clear data");

      impl::zero_field_data<NgpSpace>(fieldDataOnDevice,
                                              selfSendData,
                                              deviceMeshIndicesOffsets,
                                              ngpMesh,
                                              mesh,
                                              fieldRanks,
                                              rankFieldOffsets,
                                              iproc,
                                              dataBegin,
                                              baseProcOffset,
                                              doOperation,
                                              includeGhosts);

      Kokkos::Profiling::popRegion();
    }

    for (auto proc_location = std::begin(comm_procs); proc_location < my_proc_location; ++proc_location) {
      Kokkos::Profiling::pushRegion("NGP MPI - message waits");

      size_t proc = std::distance(std::begin(comm_procs), proc_location);
      int idx = static_cast<int>(proc);
      auto iproc = comm_procs[idx];

      MPI_Wait(&recvRequests[proc], MPI_STATUS_IGNORE);

      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("NGP MPI - message unpacking");

      if (!useDeviceAwareMpi) {
        const size_t dataBegin = hostBufferOffsets[idx];
        const size_t dataEnd = hostBufferOffsets[idx + 1];
        BufferView deviceRecvDataRange(deviceRecvData.data()+dataBegin,dataEnd-dataBegin);
        BufferView hostRecvDataRange(hostRecvData.data()+dataBegin,dataEnd-dataBegin);
        Kokkos::deep_copy(deviceRecvDataRange, hostRecvDataRange);
      }

      auto dataBegin = hostBufferOffsets[idx];
      auto baseProcOffset = hostMeshIndicesOffsets(idx);

      impl::unpack_device_recv_data<NgpSpace>(fieldDataOnDevice,
                                              deviceRecvData,
                                              deviceMeshIndicesOffsets,
                                              ngpMesh,
                                              mesh,
                                              fieldRanks,
                                              rankFieldOffsets,
                                              iproc,
                                              dataBegin,
                                              baseProcOffset,
                                              doOperation,
                                              includeGhosts);

      Kokkos::Profiling::popRegion();
    }

    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      int idx = static_cast<int>(proc);
      auto iproc = comm_procs[idx];

      Kokkos::Profiling::pushRegion("NGP MPI - cache unpacking");
      auto dataBegin = hostBufferOffsets[idx];
      auto baseProcOffset = hostMeshIndicesOffsets(idx);

      impl::unpack_device_recv_data<NgpSpace>(fieldDataOnDevice,
                                              selfSendData,
                                              deviceMeshIndicesOffsets,
                                              ngpMesh,
                                              mesh,
                                              fieldRanks,
                                              rankFieldOffsets,
                                              iproc,
                                              dataBegin,
                                              baseProcOffset,
                                              doOperation,
                                              includeGhosts);
      Kokkos::Profiling::popRegion();
    }

    for (auto proc_location = my_proc_location; proc_location < std::end(comm_procs); ++proc_location) {
      Kokkos::Profiling::pushRegion("NGP MPI - message waits");

      size_t proc = std::distance(std::begin(comm_procs), proc_location);
      int idx = static_cast<int>(proc);
      auto iproc = comm_procs[idx];

      MPI_Wait(&recvRequests[proc], MPI_STATUS_IGNORE);

      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("NGP MPI - message unpacking");

      if (!useDeviceAwareMpi) {
        const size_t dataBegin = hostBufferOffsets[idx];
        const size_t dataEnd = hostBufferOffsets[idx + 1];
        BufferView deviceRecvDataRange(deviceRecvData.data()+dataBegin,dataEnd-dataBegin);
        BufferView hostRecvDataRange(hostRecvData.data()+dataBegin,dataEnd-dataBegin);
        Kokkos::deep_copy(deviceRecvDataRange, hostRecvDataRange);
      }

      auto dataBegin = hostBufferOffsets[idx];
      auto baseProcOffset = hostMeshIndicesOffsets(idx);

      impl::unpack_device_recv_data<NgpSpace>(fieldDataOnDevice,
                                              deviceRecvData,
                                              deviceMeshIndicesOffsets,
                                              ngpMesh,
                                              mesh,
                                              fieldRanks,
                                              rankFieldOffsets,
                                              iproc,
                                              dataBegin,
                                              baseProcOffset,
                                              doOperation,
                                              includeGhosts);

      Kokkos::Profiling::popRegion();
    }
  }
  else {
    for (size_t proc = 0; proc < num_comm_procs; ++proc) {
      Kokkos::Profiling::pushRegion("NGP MPI - message waits");

      int idx = static_cast<int>(proc);
      MPI_Waitany(static_cast<int>(num_comm_procs), recvRequests.data(), &idx, MPI_STATUS_IGNORE);

      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("NGP MPI - message unpacking");

      if (!useDeviceAwareMpi) {
        const size_t dataBegin = hostBufferOffsets[idx];
        const size_t dataEnd = hostBufferOffsets[idx + 1];
        BufferView deviceRecvDataRange(deviceRecvData.data()+dataBegin,dataEnd-dataBegin);
        BufferView hostRecvDataRange(hostRecvData.data()+dataBegin,dataEnd-dataBegin);
        Kokkos::deep_copy(deviceRecvDataRange, hostRecvDataRange);
      }

      impl::DoOp<Scalar, OP> doOperation;
      auto iproc = comm_procs[idx];
      auto dataBegin = hostBufferOffsets[idx];
      auto baseProcOffset = hostMeshIndicesOffsets(idx);

      impl::unpack_device_recv_data<NgpSpace>(fieldDataOnDevice,
                                              deviceRecvData,
                                              deviceMeshIndicesOffsets,
                                              ngpMesh,
                                              mesh,
                                              fieldRanks,
                                              rankFieldOffsets,
                                              iproc,
                                              dataBegin,
                                              baseProcOffset,
                                              doOperation,
                                              includeGhosts);

      Kokkos::Profiling::popRegion();
    }
  }

  MPI_Waitall(static_cast<int>(num_comm_procs), sendRequests.data(), statuses.data());

  Kokkos::Profiling::popRegion();
#endif
}

inline
void verify_all_same_host_data_layout(const std::vector<const FieldBase*>& fields)
{
  if (fields.size() > 1) {
    const auto hostLayoutValue = fields[0]->host_data_layout();
    for(size_t i=1; i<fields.size(); ++i) {
      STK_ThrowRequireMsg(fields[i]->host_data_layout() == hostLayoutValue,"Cannot mix fields on host with different host_layout_values. Must be either all left or all right.");
    }
  }
}

template<typename Scalar, Operation OP, typename NgpSpace, typename NgpMeshType>
void ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh(const BulkData & mesh,
                                                        NgpMeshType& ngpMesh,
                                                        const std::vector<const FieldBase*>& fields,
                                                        bool includeGhosts,
                                                        bool deterministic)
{
  STK_ThrowRequireMsg(!fields.empty(), "Shouldn't get to here with empty fields vector.");
  constexpr bool isHostSpace = std::is_same_v<NgpSpace, ngp::HostSpace>;
  if constexpr (isHostSpace) {
    verify_all_same_host_data_layout(fields);
    if (fields[0]->host_data_layout() == Layout::Right) {
      ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh_layout<Scalar, OP, NgpSpace, NgpMeshType, Layout::Right>(mesh, ngpMesh, fields, includeGhosts, deterministic);
    }
    else {
      ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh_layout<Scalar, OP, NgpSpace, NgpMeshType, Layout::Left>(mesh, ngpMesh, fields, includeGhosts, deterministic);
    }
  }
  else {
    ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh_layout<Scalar, OP, NgpSpace, NgpMeshType, Layout::Left>(mesh, ngpMesh, fields, includeGhosts, deterministic);
  }
}

template<typename Scalar, Operation OP, typename NgpSpace>
void ngp_parallel_data_exchange_sym_pack_unpack(const BulkData & mesh,
                                                const std::vector<const FieldBase*>& fields,
                                                bool includeGhosts,
                                                bool deterministic)
{
  constexpr bool operatingOnHost = std::is_same_v<NgpSpace, ngp::HostSpace>;

  if constexpr (operatingOnHost) {
    //This handles the case of a non-GPU build and STK_USE_DEVICE_MESH is not defined,
    //(i.e., DeviceSpace == HostSpace),
    //or STK_USE_DEVICE_MESH is defined but the user explicitly requests HostSpace.
    HostMesh hostNgpMesh(mesh);
    ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh<Scalar, OP, NgpSpace>(mesh, hostNgpMesh, fields, includeGhosts, deterministic);
    return;
  }

  //If operatingOnHost is false, i.e., NgpSpace == DeviceSpace, then it is correct
  //to use get_updated_ngp_mesh(mesh).
  auto& ngpMesh = get_updated_ngp_mesh(mesh);
  ngp_parallel_data_exchange_sym_pack_unpack_ngpmesh<Scalar, OP, NgpSpace>(mesh, ngpMesh, fields, includeGhosts, deterministic);
}

} // namespace stk::mesh

#endif
