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

inline
std::vector<EntityRank> assemble_rank_list(const std::vector<const FieldBase*>& fields)
{
  auto fieldRanks = std::vector<EntityRank>{};
  fieldRanks.reserve(fields.size());
  for(const FieldBase* f : fields) {
    fieldRanks.push_back(f->entity_rank());
  }
  stk::util::sort_and_unique(fieldRanks);
  return fieldRanks;
}

inline
std::vector<unsigned> group_fields_by_rank(const std::vector<const FieldBase*>& fields,
                                           const std::vector<EntityRank>& fieldRanks)
{
  std::vector<unsigned> rankFieldOffsets;
  rankFieldOffsets.reserve(fieldRanks.size()+1);

  unsigned fieldRankOffset = 0;

  for(EntityRank rank : fieldRanks) {
    rankFieldOffsets.push_back(fieldRankOffset);
    for(const FieldBase* field : fields) {
      if (field->entity_rank() == rank) {
        ++fieldRankOffset;
      }
    }
  }

  rankFieldOffsets.push_back(fieldRankOffset);

  return rankFieldOffsets;
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

template <typename Scalar, typename NgpSpace, Layout LayoutValue>
auto assemble_field_data_on_device(const std::vector<const FieldBase*>& fields,
                                   const std::vector<EntityRank>& fieldRanks)
{
  using FieldDataType = decltype(fields.front()->template data<Scalar, ReadWrite, NgpSpace, LayoutValue>());
  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  auto fieldDataOnDevice = FieldDataView("fieldDataOnDevice", fields.size());
  auto fieldDataOnDeviceHost = Kokkos::create_mirror_view(stk::ngp::HostPinnedSpace{}, fieldDataOnDevice);
  unsigned offset = 0;
  for(EntityRank rank : fieldRanks) {
    for (size_t fieldIdx = 0; fieldIdx < fields.size(); ++fieldIdx) {
      const auto& field = *fields[fieldIdx];
      if (field.entity_rank() == rank) {
        STK_ThrowRequireMsg(field.type_is<Scalar>(),
                            "Cannot mix Fields with different datatypes.  Field '" <<
                            field.name() <<
                            "' is of type " <<
                            field.data_traits().type_info.name() <<
                            " when the entire set of Fields must be of type " <<
                            fields[0]->data_traits().type_info.name());
        fieldDataOnDeviceHost(offset++) = field.template data<Scalar,ReadWrite,NgpSpace,LayoutValue>();
      }
    }
  }

  auto execSpace = typename NgpSpace::exec_space{};
  Kokkos::Experimental::copy(execSpace, fieldDataOnDeviceHost, fieldDataOnDevice);
  return fieldDataOnDevice;
}

inline
std::vector<int> assemble_comm_procs_list(const BulkData& mesh, const std::vector<EntityRank>& fieldRanks, bool includeGhosts)
{
  std::vector<int> comm_procs;
  for (int proc = 0; proc < mesh.parallel_size(); ++proc) {
    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMapSize = mesh.template volatile_fast_shared_comm_map_size<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, proc, includeGhosts);
      if (sharedCommMapSize > 0) {
        comm_procs.push_back(proc);
        break;
      }
    }
  }
  stk::util::sort_and_unique(comm_procs);
  return comm_procs;
}

inline
size_t compute_total_mesh_indices_offsets(const BulkData& mesh, const std::vector<EntityRank>& fieldRanks, const std::vector<int>& commProcs, bool includeGhosts)
{
  size_t totalMeshIndicesOffsets = 0;
  for (int proc : commProcs) {
    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMapSize = mesh.template volatile_fast_shared_comm_map_size<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, proc, includeGhosts);
      if (sharedCommMapSize > 0) {
        totalMeshIndicesOffsets += sharedCommMapSize;
      }
    }
  }
  return totalMeshIndicesOffsets;
}

template <typename BufferType, typename OffsetsType>
void fill_host_buffer_offsets(const BufferType& hostBufferOffsets,
                              const OffsetsType& hostMeshIndicesOffsets,
                              const BulkData& mesh,
                              const std::vector<const FieldBase*>& fields,
                              const std::vector<int>& comm_procs,
                              const std::vector<EntityRank>& fieldRanks,
                              bool includeGhosts)
{
  hostBufferOffsets(0) = 0;
  auto num_comm_procs = comm_procs.size();
  size_t hostMeshIndicesIdx = num_comm_procs;
  for (size_t proc = 0; proc < num_comm_procs; ++proc) {
    hostMeshIndicesOffsets(proc) = hostMeshIndicesIdx;
    unsigned baseProcOffset = hostMeshIndicesIdx;

    unsigned hostMeshIndicesCounter = 0;
    unsigned hostMeshIndicesOffsetsCounter = 0;

    for (EntityRank fieldRank : fieldRanks) {
      auto sharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, comm_procs[proc], includeGhosts);
      hostMeshIndicesIdx += sharedCommMap.extent(0);

      unsigned prevBucketId = InvalidOrdinal;
      unsigned prevNumBytes = 0;
      for (size_t i = 0; i < sharedCommMap.extent(0); ++i) {
        hostMeshIndicesOffsets(baseProcOffset + hostMeshIndicesCounter + i) = hostMeshIndicesOffsetsCounter;
        unsigned bucketId = sharedCommMap(i).bucket_id;

        if (bucketId != prevBucketId) {
          prevBucketId = bucketId;
          prevNumBytes = 0;
          for (const FieldBase* field : fields) {
            if (field->entity_rank() == fieldRank) {
              prevNumBytes += field_scalars_per_entity(*field, bucketId);
            }
          }
        }

        hostMeshIndicesOffsetsCounter += prevNumBytes;
      }
      hostMeshIndicesCounter += sharedCommMap.extent(0);
    }

    hostBufferOffsets(proc+1) = hostBufferOffsets(proc) + hostMeshIndicesOffsetsCounter;
  }
}

template <typename NgpSpace, typename SendDataType, typename FieldDataType, typename OffsetsType, typename NgpMeshType>
requires ngp::is_host_space<NgpSpace>
void fill_device_send_data(const SendDataType& deviceSendData,
                           const FieldDataType& fieldDataOnDevice,
                           const OffsetsType& deviceMeshIndicesOffsets,
                           const NgpMeshType&,
                           const BulkData& mesh,
                           const std::vector<EntityRank>& fieldRanks,
                           const std::vector<unsigned>& rankFieldOffsets,
                           int iproc,
                           int dataBegin,
                           int baseProcOffset,
                           bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    const int sendBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter);
    typename SendDataType::value_type* deviceSendDataPtr = deviceSendData.data()+dataBegin+sendBufferStartIdx;
    for(size_t idx=0; idx<hostSharedCommMapSize; ++idx) {
      for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
        auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(hostSharedCommMap(idx));
        const int numScalars = entityValues.num_scalars();
        for (ScalarIdx scalar(0); scalar<numScalars; ++scalar) {
          deviceSendDataPtr[scalar] = entityValues(scalar);
        }
        deviceSendDataPtr += numScalars;
      }
    }
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

template <typename NgpSpace, typename SendDataType, typename FieldDataType, typename OffsetsType, typename NgpMeshType>
void fill_device_send_data(const SendDataType& deviceSendData,
                           const FieldDataType& fieldDataOnDevice,
                           const OffsetsType& deviceMeshIndicesOffsets,
                           const NgpMeshType& ngpMesh,
                           const BulkData& mesh,
                           const std::vector<EntityRank>& fieldRanks,
                           const std::vector<unsigned>& rankFieldOffsets,
                           int iproc,
                           int dataBegin,
                           int baseProcOffset,
                           bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    Kokkos::parallel_for("fill_device_send_data", stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMapSize),
      KOKKOS_LAMBDA(size_t idx) {
          auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iproc, includeGhosts);
          if (idx >= deviceSharedCommMap.extent(0)) {
            return;
          }
          const int sendBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + idx);
          typename SendDataType::value_type* deviceSendDataPtr = deviceSendData.data()+dataBegin+sendBufferStartIdx;
          for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
            auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(deviceSharedCommMap(idx));
            const int numScalars = entityValues.num_scalars();
            for (ScalarIdx scalar(0); scalar<numScalars; ++scalar) {
              deviceSendDataPtr[scalar] = entityValues(scalar);
            }
            deviceSendDataPtr += numScalars;
          }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

template <typename NgpSpace, typename FieldDataType, typename RecvDataType, typename OffsetsType, typename NgpMeshType, typename OP>
requires ngp::is_host_space<NgpSpace>
void zero_field_data(const FieldDataType& fieldDataOnDevice,
                             const RecvDataType& selfSendData,
                             const OffsetsType& deviceMeshIndicesOffsets,
                             const NgpMeshType&,
                             const BulkData& mesh,
                             const std::vector<EntityRank>& fieldRanks,
                             const std::vector<unsigned>& rankFieldOffsets,
                             int iproc,
                             int dataBegin,
                             int baseProcOffset,
                             const OP& doOperation,
                             bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    Kokkos::parallel_for("zero_field_data", stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMapSize),
      [&](size_t index) {
        const int selfSendStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + index);
        typename RecvDataType::value_type* selfSendDataPtr = selfSendData.data()+dataBegin+selfSendStartIdx;
        for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
          auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(hostSharedCommMap(index));
          for (ScalarIdx scalar : entityValues.scalars()) {
            selfSendDataPtr[scalar] = entityValues(scalar);
            entityValues(scalar) = doOperation.initial_value();
          }
          selfSendDataPtr += entityValues.num_scalars();
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

template <typename NgpSpace, typename FieldDataType, typename RecvDataType, typename OffsetsType, typename NgpMeshType, typename OP>
void zero_field_data(const FieldDataType& fieldDataOnDevice,
                             const RecvDataType& selfSendData,
                             const OffsetsType& deviceMeshIndicesOffsets,
                             const NgpMeshType& ngpMesh,
                             const BulkData& mesh,
                             const std::vector<EntityRank>& fieldRanks,
                             const std::vector<unsigned>& rankFieldOffsets,
                             int iproc,
                             int dataBegin,
                             int baseProcOffset,
                             const OP& doOperation,
                             bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    Kokkos::parallel_for("zero_field_data", stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMapSize),
      KOKKOS_LAMBDA(size_t index) {
        auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iproc, includeGhosts);
        if (index >= deviceSharedCommMap.extent(0)) {
          return;
        }
        const int selfSendStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + index);
        typename RecvDataType::value_type* selfSendDataPtr = selfSendData.data()+dataBegin+selfSendStartIdx;
        for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
          auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(deviceSharedCommMap(index));
          for (ScalarIdx scalar : entityValues.scalars()) {
            selfSendDataPtr[scalar] = entityValues(scalar);
            entityValues(scalar) = doOperation.initial_value();
          }
          selfSendDataPtr += entityValues.num_scalars();
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

template <typename NgpSpace, typename FieldDataType, typename RecvDataType, typename OffsetsType, typename NgpMeshType, typename OP>
requires ngp::is_host_space<NgpSpace>
void unpack_device_recv_data(const FieldDataType& fieldDataOnDevice,
                             const RecvDataType& deviceRecvData,
                             const OffsetsType& deviceMeshIndicesOffsets,
                             const NgpMeshType&,
                             const BulkData& mesh,
                             const std::vector<EntityRank>& fieldRanks,
                             const std::vector<unsigned>& rankFieldOffsets,
                             int iproc,
                             int dataBegin,
                             int baseProcOffset,
                             const OP& doOperation,
                             bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    Kokkos::parallel_for("unpack_device_recv_data", stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMapSize),
      [&](size_t index) {
        const int recvBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + index);
        typename RecvDataType::value_type* deviceRecvDataPtr = deviceRecvData.data()+dataBegin+recvBufferStartIdx;
        for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
          auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(hostSharedCommMap(index));
          for (ScalarIdx scalar : entityValues.scalars()) {
            entityValues(scalar) = doOperation(entityValues(scalar), deviceRecvDataPtr[scalar]);
          }
          deviceRecvDataPtr += entityValues.num_scalars();
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

template <typename NgpSpace, typename FieldDataType, typename RecvDataType, typename OffsetsType, typename NgpMeshType, typename OP>
void unpack_device_recv_data(const FieldDataType& fieldDataOnDevice,
                             const RecvDataType& deviceRecvData,
                             const OffsetsType& deviceMeshIndicesOffsets,
                             const NgpMeshType& ngpMesh,
                             const BulkData& mesh,
                             const std::vector<EntityRank>& fieldRanks,
                             const std::vector<unsigned>& rankFieldOffsets,
                             int iproc,
                             int dataBegin,
                             int baseProcOffset,
                             const OP& doOperation,
                             bool includeGhosts)
{
  size_t meshIndicesCounter = 0;
  for (unsigned i=0; i<fieldRanks.size(); ++i) {
    EntityRank fieldRank = fieldRanks[i];
    auto fieldRange = Kokkos::pair{rankFieldOffsets[i], rankFieldOffsets[i+1]};

    auto hostSharedCommMap = mesh.template volatile_fast_shared_comm_map<typename stk::ngp::DeviceSpace::mem_space>(fieldRank, iproc, includeGhosts);
    unsigned hostSharedCommMapSize = hostSharedCommMap.extent(0);

    Kokkos::parallel_for("unpack_device_recv_data", stk::ngp::RangePolicy<typename NgpSpace::exec_space>(0, hostSharedCommMapSize),
      KOKKOS_LAMBDA(size_t index) {
        auto deviceSharedCommMap = ngpMesh.volatile_fast_shared_comm_map(fieldRank, iproc, includeGhosts);
        if (index >= deviceSharedCommMap.extent(0)) {
          return;
        }
        const int recvBufferStartIdx = deviceMeshIndicesOffsets(baseProcOffset + meshIndicesCounter + index);
        typename RecvDataType::value_type* deviceRecvDataPtr = deviceRecvData.data()+dataBegin+recvBufferStartIdx;
        for (unsigned fieldIdx = fieldRange.first; fieldIdx < fieldRange.second; ++fieldIdx) {
          auto entityValues = fieldDataOnDevice(fieldIdx).entity_values(deviceSharedCommMap(index));
          for (ScalarIdx scalar : entityValues.scalars()) {
            entityValues(scalar) = doOperation(entityValues(scalar), deviceRecvDataPtr[scalar]);
          }
          deviceRecvDataPtr += entityValues.num_scalars();
        }
      }
    );
    Kokkos::fence();
    meshIndicesCounter += hostSharedCommMapSize;
  }
}

} // namespace stk::mesh::impl

#endif
