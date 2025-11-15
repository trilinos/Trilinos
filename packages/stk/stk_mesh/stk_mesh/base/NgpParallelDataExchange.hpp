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

#ifndef STK_MESH_NGP_PARALLEL_DATA_EXCHANGE_HPP
#define STK_MESH_NGP_PARALLEL_DATA_EXCHANGE_HPP

#include <stk_util/util/SortAndUnique.hpp>
#include <Kokkos_Core.hpp>
#include "stk_mesh/base/BulkData.hpp"

namespace stk::mesh::impl {

inline std::vector<EntityRank> assemble_rank_list(const std::vector<const FieldBase*>& fields)
{
  auto fieldRanks = std::vector<EntityRank>{};
  fieldRanks.reserve(fields.size());
  for(const FieldBase* f : fields) {
    fieldRanks.push_back(f->entity_rank());
  }
  stk::util::sort_and_unique(fieldRanks);
  return fieldRanks;
}

template <typename NgpSpace>
EntityRankViewType<typename NgpSpace::mem_space> assemble_rank_per_field(const std::vector<const FieldBase*>& fields)
{
  auto rankPerField = EntityRankViewType<typename NgpSpace::mem_space>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "fieldRanksOnDevice"), fields.size());
  auto rankPerFieldHost = Kokkos::create_mirror_view(rankPerField);
  for(size_t fieldIdx = 0; fieldIdx < fields.size(); ++fieldIdx) {
    rankPerFieldHost(fieldIdx) = fields[fieldIdx]->entity_rank();
  }
  Kokkos::deep_copy(rankPerField, rankPerFieldHost);
  return rankPerField;
}

template <typename Scalar, typename NgpSpace>
auto assemble_field_data_on_device(const std::vector<const FieldBase*>& fields)
{
  using FieldDataType = decltype(fields.front()->template data<Scalar, ReadWrite, NgpSpace>());
  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  auto fieldDataOnDevice = FieldDataView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "fieldDataOnDevice"), fields.size());
  auto fieldDataOnDeviceHost = Kokkos::create_mirror_view(fieldDataOnDevice);
  for (size_t fieldIdx = 0; fieldIdx < fields.size(); ++fieldIdx) {
    fieldDataOnDeviceHost(fieldIdx) = fields[fieldIdx]->template data<Scalar, ReadWrite, NgpSpace>();
  }
  Kokkos::deep_copy(fieldDataOnDevice, fieldDataOnDeviceHost);
  return fieldDataOnDevice;
}

template <typename NgpSpace>
std::vector<int> assemble_comm_procs_list(const BulkData& mesh, const std::vector<EntityRank>& fieldRanks, bool includeGhosts)
{
  std::vector<int> comm_procs;
  for (int proc = 0; proc < mesh.parallel_size(); ++proc) {
    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMapSize = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, proc, includeGhosts).extent(0);
      if (sharedCommMapSize > 0) {
        comm_procs.push_back(proc);
      }
    }
  }
  stk::util::sort_and_unique(comm_procs);
  return comm_procs;
}

template <typename NgpSpace>
size_t compute_total_mesh_indices_offsets(const BulkData& mesh, const std::vector<EntityRank>& fieldRanks, bool includeGhosts)
{
  size_t totalMeshIndicesOffsets = 0;
  for (int proc = 0; proc < mesh.parallel_size(); ++proc) {
    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMapSize = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, proc, includeGhosts).extent(0);
      if (sharedCommMapSize > 0) {
        totalMeshIndicesOffsets += sharedCommMapSize;
      }
    }
  }
  return totalMeshIndicesOffsets;
}

template <typename NgpSpace>
std::vector<size_t> compute_message_sizes(const BulkData& mesh,
                                          const std::vector<const FieldBase*>& fields,
                                          const std::vector<int>& comm_procs,
                                          const std::vector<EntityRank>& fieldRanks,
                                          bool includeGhosts)
{
  std::vector<size_t> messageSizes;
  messageSizes.reserve(comm_procs.size());
  for (auto iproc : comm_procs) {
    size_t numValues = 0;
    for (EntityRank fieldRank : fieldRanks) {
      auto commMapIndices = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, iproc, includeGhosts);
      for (size_t i = 0; i < commMapIndices.extent(0); ++i) {
        const unsigned bucketId = commMapIndices(i).bucket_id;
        for (const FieldBase* field : fields) {
          if (field->entity_rank() == fieldRank) {
            numValues += field_scalars_per_entity(*field, bucketId);
          }
        }
      }
    }
    messageSizes.push_back(numValues);
  }
  return messageSizes;
}

template <typename NgpSpace, typename BufferType, typename OffsetsType>
void fill_host_buffer_offsets(const BufferType& hostBufferOffsets,
                              const OffsetsType& hostMeshIndicesOffsets,
                              const BulkData& mesh,
                              const std::vector<const FieldBase*>& fields,
                              const std::vector<int>& comm_procs,
                              const std::vector<EntityRank>& fieldRanks,
                              const std::vector<size_t>& messageSizes,
                              bool includeGhosts)
{
  hostBufferOffsets(0) = 0;
  auto num_comm_procs = comm_procs.size();
  size_t hostMeshIndicesIdx = num_comm_procs;
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    hostBufferOffsets[proc+1] = hostBufferOffsets[proc] + messageSizes[proc];

    hostMeshIndicesOffsets(proc) = hostMeshIndicesIdx;
    unsigned baseProcOffset = hostMeshIndicesIdx;

    unsigned hostMeshIndicesCounter = 0;
    unsigned hostMeshIndicesOffsetsCounter = 0;

    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMap = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, comm_procs[proc], includeGhosts);
      hostMeshIndicesIdx += sharedCommMap.extent(0);

      for (size_t i = 0; i < sharedCommMap.extent(0); ++i) {
        hostMeshIndicesOffsets(baseProcOffset + hostMeshIndicesCounter + i) = hostMeshIndicesOffsetsCounter;
        auto meshIndex = sharedCommMap(i);

        for (const FieldBase* field : fields) {
          if (field->entity_rank() == fieldRank) {
            hostMeshIndicesOffsetsCounter += field_scalars_per_entity(*field, meshIndex.bucket_id);
          }
        }
      }
      hostMeshIndicesCounter += sharedCommMap.extent(0);
    }
  }
}

template <typename NgpSpace, typename SendDataType, typename FieldDataType, typename RankPerFieldType, typename OffsetsType, typename NgpMeshType>
void fill_device_send_data(const SendDataType& deviceSendData,
                           const FieldDataType& fieldDataOnDevice,
                           const RankPerFieldType& rankPerField,
                           const OffsetsType& deviceMeshIndicesOffsets,
                           const NgpMeshType& ngpMesh,
                           const BulkData& mesh,
                           const std::vector<EntityRank>& fieldRanks,
                           int iproc,
                           int dataBegin,
                           int baseProcOffset,
                           bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (EntityRank fieldRank : fieldRanks) {
    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, iproc, includeGhosts);

    Kokkos::parallel_for(stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMap.extent(0)),
      KOKKOS_LAMBDA(size_t idx) {
        auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iproc, includeGhosts);
        if (idx >= deviceSharedCommMap.extent(0)) {
          return;
        }
        auto fastMeshIdx = deviceSharedCommMap(idx);

        int sendBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + idx);
        for (size_t fieldIdx = 0; fieldIdx < fieldDataOnDevice.extent(0); ++fieldIdx) {
          if (fieldRank == rankPerField(fieldIdx)) {
            auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(fastMeshIdx);
            for (ScalarIdx scalar : entityValues.scalars()) {
              deviceSendData(dataBegin + sendBufferStartIdx++) = entityValues(scalar);
            }
          }
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMap.extent(0);
  }
}

template <typename NgpSpace, typename FieldDataType, typename RecvDataType, typename RankPerFieldType, typename OffsetsType, typename NgpMeshType, typename OP>
void unpack_device_recv_data(const FieldDataType& fieldDataOnDevice,
                             const RecvDataType& deviceRecvData,
                             const RankPerFieldType& rankPerField,
                             const OffsetsType& deviceMeshIndicesOffsets,
                             const NgpMeshType& ngpMesh,
                             const BulkData& mesh,
                             const std::vector<EntityRank>& fieldRanks,
                             int iproc,
                             int dataBegin,
                             int baseProcOffset,
                             const OP& doOperation,
                             bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (EntityRank fieldRank : fieldRanks) {
    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename NgpSpace::mem_space>(fieldRank, iproc, includeGhosts);

    Kokkos::parallel_for(stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMap.extent(0)),
      KOKKOS_LAMBDA(size_t index) {
        auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iproc, includeGhosts);
        if (index >= deviceSharedCommMap.extent(0)) {
          return;
        }
        auto fastMeshIndex = deviceSharedCommMap(index);

        int recvBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + index);
        for (size_t fieldIdx = 0; fieldIdx < fieldDataOnDevice.extent(0); ++fieldIdx) {
          if (fieldRank == rankPerField(fieldIdx)) {
            auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(fastMeshIndex);
            for (ScalarIdx scalar : entityValues.scalars()) {
              auto rcv = deviceRecvData(dataBegin + recvBufferStartIdx++);
              entityValues(scalar) = doOperation(entityValues(scalar), rcv);
            }
          }
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMap.extent(0);
  }
}

} // namespace stk::mesh::impl

#endif
