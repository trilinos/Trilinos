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

#include <stk_util/stk_config.h>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <stk_util/util/SortAndUnique.hpp>   // for PairIter

#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/HostMesh.hpp>
#include <stk_mesh/base/NgpFieldParallel.hpp>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/MeshCommImplUtils.hpp>

#include <utility>                      // for pair
#include <sstream>                      // for basic_ostream::operator<<, etc
#include <set>

namespace stk {
namespace mesh {

bool find_proc_before_index(const EntityCommInfoVector& infovec, int proc, int index)
{
    for(int i=0; i<index; ++i) {
        if (proc == infovec[i].proc) {
            return true;
        }
    }
    return false;
}

struct InfoRankLess {
  bool operator()(const EntityCommListInfo& info, EntityRank rank) const
  { return info.key.rank() < rank; }
  bool operator()(EntityRank rank, const EntityCommListInfo& info) const
  { return rank < info.key.rank(); }
};

void compute_field_entity_ranges(const EntityCommListInfoVector& commListInfo,
                                 const std::vector<const FieldBase*>& fields,
                                 std::vector<std::pair<int,int>>& fieldRange)
{
  int numFields = fields.size();
  for (int fi = 0; fi < numFields; ++fi) {
    EntityRank fieldRank = fields[fi]->entity_rank();
    EntityCommListInfoVector::const_iterator startIter = std::lower_bound(commListInfo.begin(), commListInfo.end(),
                                                                          fieldRank, InfoRankLess());
    EntityCommListInfoVector::const_iterator endIter = std::upper_bound(startIter, commListInfo.end(), fieldRank,
                                                                        InfoRankLess());
    fieldRange[fi].first = std::distance(commListInfo.begin(), startIter);
    fieldRange[fi].second = std::distance(commListInfo.begin(), endIter);
  }
}

void do_initial_sync_to_host(const std::vector<const FieldBase*>& fields)
{
  for(const FieldBase* fptr : fields) {
    fptr->synchronize<stk::mesh::ReadWrite>();  // Sync to host and mark as modified
  }
}

template <Layout DataLayout>
void copy_ghost_entity_data_into_buffer(const ConstFieldDataBytes<stk::ngp::HostMemSpace>& fieldDataBytes,
                                        const EntityCommListInfo& ecli, PairIterEntityComm info, unsigned ghostId,
                                        const MeshIndex& meshIndex, std::byte* sendBuffer,
                                        const std::vector<int>& procOffset, std::vector<int>& entityOffset)
{
  for (; !info.empty(); ++info) {
    if (info->ghost_id == ghostId) {
      const int proc = info->proc;
      std::byte* bufferPtr = &sendBuffer[procOffset[proc]+entityOffset[proc]];
      auto entityBytes = fieldDataBytes.entity_bytes<DataLayout>(meshIndex);
      for (ByteIdx byte : entityBytes.bytes()) {
        bufferPtr[byte] = entityBytes(byte);
      }
      entityOffset[proc] += entityBytes.num_bytes();
    }
  }
}

template <Layout DataLayout>
void copy_ghost_entity_data_out_of_buffer(const FieldDataBytes<stk::ngp::HostMemSpace>& fieldDataBytes,
                                          const EntityCommListInfo& ecli, PairIterEntityComm info, unsigned ghostId,
                                          int owner, const MeshIndex& meshIndex, const std::byte* recvBuffer,
                                          const std::vector<int>& procOffset, std::vector<int>& entityOffset)
{
  for (; !info.empty(); ++info) {
    if (info->ghost_id == ghostId) {
      const std::byte* bufferPtr = &recvBuffer[procOffset[owner]+entityOffset[owner]];
      auto entityBytes = fieldDataBytes.entity_bytes<DataLayout>(meshIndex);
      for (ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = bufferPtr[byte];
      }
      entityOffset[owner] += entityBytes.num_bytes();
      break;
    }
  }
}

//----------------------------------------------------------------------
void communicate_field_data(const Ghosting& ghosts, const std::vector<const FieldBase*>& fields)
{
  if ( fields.empty()) {
    return;
  }

  const BulkData & mesh = ghosts.mesh();
  const int parallelSize = mesh.parallel_size();

  mesh.confirm_host_mesh_is_synchronized_from_device();
  do_initial_sync_to_host(fields);

  if (parallelSize == 1) {
    return;
  }

  const int numFields = fields.size();

  std::vector<int> sendSizes(parallelSize);
  std::vector<int> recvSizes(parallelSize);
  std::vector<int> sendOffsets(parallelSize + 1);
  std::vector<int> recvOffsets(parallelSize + 1);

  std::vector<std::pair<int,int>> fieldEntityRange(numFields, std::make_pair(-1, -1));

  const EntityCommListInfoVector& commInfoVec = mesh.internal_comm_list();
  compute_field_entity_ranges(commInfoVec, fields, fieldEntityRange);

  const unsigned ghostId = ghosts.ordinal();
  const EntityCommDatabase& commDB = mesh.internal_comm_db();
  std::vector<int> commProcs;

  for (const EntityCommListInfo& ecli : mesh.internal_comm_list()) {
    Entity entity = ecli.entity;
    const MeshIndex meshIdx = mesh.mesh_index(entity);
    const unsigned bucketId = meshIdx.bucket->bucket_id();
    EntityRank rank = meshIdx.bucket->entity_rank();

    unsigned bytesPerEntityAllFields = 0;
    for ( const FieldBase* fptr : fields) {
      const FieldBase& field = *fptr;

      if (is_matching_rank(field, rank)) {
        bytesPerEntityAllFields += field_bytes_per_entity(field, bucketId);
      }
    }

    if (bytesPerEntityAllFields == 0) {
      continue;
    }

    const bool owned = meshIdx.bucket->owned();
    PairIterEntityComm info = commDB.comm(ecli.entity_comm);
    if (owned) {
      for (; !info.empty(); ++info) {
        if (info->ghost_id == ghostId) {
          sendSizes[info->proc] += bytesPerEntityAllFields;
        }
      }
    }
    else {
      const int owner = mesh.parallel_owner_rank(entity);
      for (; !info.empty(); ++info) {
        if (info->ghost_id == ghostId) {
          recvSizes[owner] += bytesPerEntityAllFields;
          break;  //jump out since we know we're only receiving one message for this entity from the one-and-only owner
        }
      }
    }
  }

  for (int p = 1; p < (parallelSize + 1); ++p) {
    sendOffsets[p] = sendOffsets[p-1] + sendSizes[p-1];
    recvOffsets[p] = recvOffsets[p-1] + recvSizes[p-1];
    sendSizes[p-1] = 0;  // Done with this data; Will use as a running offset to fill buffer
    recvSizes[p-1] = 0;
  }

  const int totalSendSize = sendOffsets.back();
  const int totalRecvSize = recvOffsets.back();

  std::unique_ptr<std::byte[]> sendBuffer(new std::byte[totalSendSize]);
  std::unique_ptr<std::byte[]> recvBuffer(new std::byte[totalRecvSize]);
  const int parallelRank = mesh.parallel_rank();

  for (const EntityCommListInfo& ecli : mesh.internal_comm_list()) {
    const int owner = mesh.parallel_owner_rank(ecli.entity);
    if (owner == parallelRank) {
      Entity e = ecli.entity;
      const MeshIndex meshIdx = mesh.mesh_index(e);
      EntityRank rank = meshIdx.bucket->entity_rank();
      PairIterEntityComm info = commDB.comm(ecli.entity_comm);

      for (const FieldBase* fptr : fields) {
        const FieldBase& field = *fptr;

        if (!is_matching_rank(field, rank)) continue;

        auto& fieldDataBytes = field.data_bytes<const std::byte>();
        if (field.host_data_layout() == Layout::Right) {
          copy_ghost_entity_data_into_buffer<Layout::Right>(fieldDataBytes, ecli, info, ghostId, meshIdx,
                                                            sendBuffer.get(), sendOffsets, sendSizes);
        }
        else if (field.host_data_layout() == Layout::Left) {
          copy_ghost_entity_data_into_buffer<Layout::Left>(fieldDataBytes, ecli, info, ghostId, meshIdx,
                                                           sendBuffer.get(), sendOffsets, sendSizes);
        }
        else {
          STK_ThrowErrorMsg("Unsupported host Field layout detected: " << field.host_data_layout());
        }
      }
    }
  }

  parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendBuffer.get(),
                                              recvOffsets.data(), recvBuffer.get(), mesh.parallel(),
                                              mesh.is_mesh_consistency_check_on());

  for (const EntityCommListInfo& ecli : mesh.internal_comm_list()) {
    const int owner = mesh.parallel_owner_rank(ecli.entity);
    if (owner != parallelRank) {
      Entity e = ecli.entity;
      const MeshIndex meshIdx = mesh.mesh_index(e);
      EntityRank rank = meshIdx.bucket->entity_rank();

      for (const FieldBase* fptr : fields) {
        const FieldBase& field = *fptr;

        if (!is_matching_rank(field, rank)) continue;

        PairIterEntityComm info = commDB.comm(ecli.entity_comm);

        auto& fieldDataBytes = field.data_bytes<std::byte>();
        if (field.host_data_layout() == Layout::Right) {
          copy_ghost_entity_data_out_of_buffer<Layout::Right>(fieldDataBytes, ecli, info, ghostId, owner, meshIdx,
                                                              recvBuffer.get(), recvOffsets, recvSizes);
        }
        else if (field.host_data_layout() == Layout::Left) {
          copy_ghost_entity_data_out_of_buffer<Layout::Left>(fieldDataBytes, ecli, info, ghostId, owner, meshIdx,
                                                             recvBuffer.get(), recvOffsets, recvSizes);
        }
        else {
          STK_ThrowErrorMsg("Unsupported host Field layout detected: " << field.host_data_layout());
        }
      }
    }
  }
}


template <Layout DataLayout>
void copy_entity_data_into_buffer(const BulkData& bulk,
                                  const ConstFieldDataBytes<stk::ngp::HostMemSpace>& fieldDataBytes,
                                  const std::pair<int, int>& fieldRange, const EntityCommDatabase& commDB,
                                  const EntityCommListInfoVector& commInfoVec, std::vector<int>& commProcs,
                                  std::byte* sendBuffer, const std::vector<int>& procOffset,
                                  std::vector<int>& entityOffset)
{
  for (int i = fieldRange.first; i < fieldRange.second; ++i) {
    const MeshIndex& meshIndex = bulk.mesh_index(commInfoVec[i].entity);
    const Bucket& bucket = *meshIndex.bucket;

    if (bucket.owned()) {
      auto entityBytes = fieldDataBytes.entity_bytes<DataLayout>(meshIndex);

      if (entityBytes.is_field_defined()) {
        impl::fill_sorted_procs(commDB.comm(commInfoVec[i].entity_comm), commProcs);
        for (int proc : commProcs) {
          std::byte* bufferPtr = &sendBuffer[procOffset[proc]+entityOffset[proc]];
          for (ByteIdx byte : entityBytes.bytes()) {
            bufferPtr[byte] = entityBytes(byte);
          }
          entityOffset[proc] += entityBytes.num_bytes();
        }
      }
    }
  }
}

template <Layout DataLayout>
void copy_entity_data_out_of_buffer(const BulkData& bulk, FieldDataBytes<stk::ngp::HostMemSpace>& fieldDataBytes,
                                    const std::pair<int, int>& fieldRange, const EntityCommListInfoVector& commInfoVec,
                                    const std::byte* recvBuffer, const std::vector<int>& procOffset,
                                    std::vector<int>& entityOffset)
{
  const int myRank = bulk.parallel_rank();
  for (int i = fieldRange.first; i < fieldRange.second; ++i) {
    const int owner = bulk.parallel_owner_rank(commInfoVec[i].entity);
    if (owner == myRank) {
      continue;
    }

    const int recvSize = procOffset[owner+1] - procOffset[owner];
    if (recvSize == 0) {
      continue;
    }

    const MeshIndex& meshIndex = bulk.mesh_index(commInfoVec[i].entity);

    auto entityBytes = fieldDataBytes.entity_bytes<DataLayout>(meshIndex);
    const std::byte* bufferPtr = &recvBuffer[procOffset[owner] + entityOffset[owner]];
    for (ByteIdx byte : entityBytes.bytes()) {
      entityBytes(byte) = bufferPtr[byte];
    }
    entityOffset[owner] += entityBytes.num_bytes();
  }
}

void communicate_field_data(const BulkData& mesh, const std::vector<const FieldBase*>& fields)
{
  mesh.confirm_host_mesh_is_synchronized_from_device();
  const int parallelSize = mesh.parallel_size();
  if ( fields.empty() || parallelSize == 1) {
    return;
  }

  do_initial_sync_to_host(fields);

  const int numFields = fields.size();

  std::vector<int> sendSizes(parallelSize);
  std::vector<int> recvSizes(parallelSize);
  std::vector<int> sendOffsets(parallelSize + 1);
  std::vector<int> recvOffsets(parallelSize + 1);

  std::vector<std::pair<int,int>> fieldEntityRange(fields.size(), std::make_pair(-1, -1));

  const EntityCommListInfoVector& commInfoVec = mesh.internal_comm_list();
  compute_field_entity_ranges(commInfoVec, fields, fieldEntityRange);

  const EntityCommDatabase& commDB = mesh.internal_comm_db();
  std::vector<int> commProcs;

  for (int fi = 0; fi < numFields; ++fi) {
    const FieldBase& field = *fields[fi];
    for(int i = fieldEntityRange[fi].first; i<fieldEntityRange[fi].second; ++i) {
      const MeshIndex& meshIndex = mesh.mesh_index(commInfoVec[i].entity);
      const Bucket& bucket = *meshIndex.bucket;
      const unsigned bucketId = bucket.bucket_id();
      const unsigned bytesPerEntity = field_bytes_per_entity(field, bucketId);

      if (bytesPerEntity == 0) {
        continue;
      }

      const bool owned = bucket.owned();
      if (owned) {
        impl::fill_sorted_procs(commDB.comm(commInfoVec[i].entity_comm), commProcs);
        for (int proc : commProcs) {
          sendSizes[proc] += bytesPerEntity;
        }
      }
      else {
        const int owner = mesh.parallel_owner_rank(commInfoVec[i].entity);
        recvSizes[owner] += bytesPerEntity;
      }
    }
  }

  for (int p = 1; p < (parallelSize + 1); ++p) {
    sendOffsets[p] = sendOffsets[p-1] + sendSizes[p-1];
    recvOffsets[p] = recvOffsets[p-1] + recvSizes[p-1];
    sendSizes[p-1] = 0;  // Done with this data; Will use as a running offset to fill buffer
    recvSizes[p-1] = 0;
  }

  const int totalSendSize = sendOffsets.back();
  const int totalRecvSize = recvOffsets.back();

  std::unique_ptr<std::byte[]> sendBuffer(new std::byte[totalSendSize]);
  std::unique_ptr<std::byte[]> recvBuffer(new std::byte[totalRecvSize]);

  for (int fi = 0; fi < numFields; ++fi) {
    const FieldBase& field = *fields[fi];
    auto& fieldDataBytes = field.data_bytes<const std::byte>();

    if (field.host_data_layout() == Layout::Right) {
      copy_entity_data_into_buffer<Layout::Right>(mesh, fieldDataBytes, fieldEntityRange[fi], commDB, commInfoVec,
                                                  commProcs, sendBuffer.get(), sendOffsets, sendSizes);
    }
    else if (field.host_data_layout() == Layout::Left) {
      copy_entity_data_into_buffer<Layout::Left>(mesh, fieldDataBytes, fieldEntityRange[fi], commDB, commInfoVec,
                                                 commProcs, sendBuffer.get(), sendOffsets, sendSizes);
    }
    else {
      STK_ThrowErrorMsg("Unsupported host Field layout detected: " << field.host_data_layout());
    }
  }

  parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendBuffer.get(),
                                              recvOffsets.data(), recvBuffer.get(), mesh.parallel(),
                                              mesh.is_mesh_consistency_check_on());

  //now unpack and store the recvd data
  for (int fi = 0; fi < numFields; ++fi) {
    const FieldBase& field = *fields[fi];
    auto& fieldDataBytes = field.data_bytes<std::byte>();

    if (field.host_data_layout() == Layout::Right) {
      copy_entity_data_out_of_buffer<Layout::Right>(mesh, fieldDataBytes, fieldEntityRange[fi], commInfoVec,
                                                    recvBuffer.get(), recvOffsets, recvSizes);
    }
    else if (field.host_data_layout() == Layout::Left) {
      copy_entity_data_out_of_buffer<Layout::Left>(mesh, fieldDataBytes, fieldEntityRange[fi], commInfoVec,
                                                   recvBuffer.get(), recvOffsets, recvSizes);
    }
    else {
      STK_ThrowErrorMsg("Unsupported host Field layout detected: " << field.host_data_layout());
    }
  }
}

//----------------------------------------------------------------------

/** Sum (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */

namespace {

template <typename T, Operation OP>
struct DoOp;

template <typename T>
struct DoOp<T, Operation::SUM>
{
  T operator()(T lhs, T rhs) const
  { return lhs + rhs; }
};

template <typename T>
struct DoOp<T, Operation::MIN>
{
  T operator()(T lhs, T rhs) const
  { return lhs < rhs ? lhs : rhs; }
};

template <typename T>
struct DoOp<T, Operation::MAX>
{
  T operator()(T lhs, T rhs) const
  { return lhs > rhs ? lhs : rhs; }
};

template <typename T, Operation OP>
void parallel_op_impl(const BulkData& mesh, std::vector<const FieldBase*> fields, bool deterministic = false)
{
  if (fields.empty()) {
    return;
  }
  mesh.confirm_host_mesh_is_synchronized_from_device();

  std::vector<int> commProcs = mesh.all_sharing_procs(fields[0]->entity_rank());
  stk::mesh::EntityRank firstFieldRank = fields[0]->entity_rank();
  for (size_t i = 1; i < fields.size(); ++i) {
    const FieldBase& field = *fields[i];
    if (field.entity_rank() != firstFieldRank) {
      const std::vector<int>& sharingProcs = mesh.all_sharing_procs(field.entity_rank());
      for (int p : sharingProcs) {
        stk::util::insert_keep_sorted_and_unique(p, commProcs);
      }
    }
  }

  auto msgPacker = [&fields, &mesh](int proc, std::vector<T>& sendBuffer)
  {
    sendBuffer.clear();
    size_t totalSize = 0;
    for (size_t j = 0; j < fields.size(); ++j) {
      const FieldBase& f = *fields[j];
      STK_ThrowRequireMsg(f.type_is<T>(), "Cannot mix Fields with different datatypes.  Field '" << f.name() <<
                          "' is of type " << f.data_traits().type_info.name() << " when the entire set of Fields " <<
                          "must be of type " << fields[0]->data_traits().type_info.name());

      auto hostCommMapIndices = mesh.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(f.entity_rank(), proc);
      for (size_t i = 0; i < hostCommMapIndices.extent(0); ++i) {
        const unsigned bucketId = hostCommMapIndices(i).bucket_id;
        const unsigned scalarsPerEntity = field_bytes_per_entity(f, bucketId) / sizeof(T);
        totalSize += scalarsPerEntity;
      }
    }
    sendBuffer.reserve(totalSize);

    for (size_t j = 0; j < fields.size(); ++j) {
      const FieldBase& f = *fields[j];
      auto hostCommMapIndices = mesh.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(f.entity_rank(), proc);

      if (f.host_data_layout() == Layout::Right) {
        auto fieldData = f.data<T, ReadOnly, stk::ngp::HostMemSpace, Layout::Right>();
        for (size_t i = 0; i < hostCommMapIndices.extent(0); ++i) {
          const unsigned bucketId = hostCommMapIndices(i).bucket_id;
          const unsigned bucketOrd = hostCommMapIndices(i).bucket_ord;
          auto entityValues = fieldData.entity_values(FastMeshIndex{bucketId, bucketOrd});
          for (ScalarIdx scalar : entityValues.scalars()) {
            sendBuffer.push_back(entityValues(scalar));
          }
        }
      }
      else if (f.host_data_layout() == Layout::Left) {
        auto fieldData = f.data<T, ReadOnly, stk::ngp::HostMemSpace, Layout::Left>();
        for (size_t i = 0; i < hostCommMapIndices.extent(0); ++i) {
          const unsigned bucketId = hostCommMapIndices(i).bucket_id;
          const unsigned bucketOrd = hostCommMapIndices(i).bucket_ord;
          auto entityValues = fieldData.entity_values(FastMeshIndex{bucketId, bucketOrd});
          for (ScalarIdx scalar : entityValues.scalars()) {
            sendBuffer.push_back(entityValues(scalar));
          }
        }
      }
      else {
        STK_ThrowErrorMsg("Unsupported host Field data layout: " << f.host_data_layout());
      }
    }
  };

  auto msgUnpacker = [&fields, &mesh](int proc, std::vector<T>& recvBuffer)
  {
    DoOp<T, OP> doOp;

    unsigned recvOffset = 0;
    for (size_t j = 0; j < fields.size(); ++j) {
      const FieldBase& f = *fields[j] ;
      auto hostCommMapIndices = mesh.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(f.entity_rank(), proc);

      if (f.host_data_layout() == Layout::Right) {
        auto fieldData = f.data<T, ReadWrite, stk::ngp::HostMemSpace, Layout::Right>();
        for (size_t i = 0; i < hostCommMapIndices.extent(0); ++i) {
          const unsigned bucketId = hostCommMapIndices(i).bucket_id;
          const unsigned bucketOrd = hostCommMapIndices(i).bucket_ord;
          auto entityValues = fieldData.entity_values(FastMeshIndex{bucketId, bucketOrd});
          for (ScalarIdx scalar : entityValues.scalars()) {
            entityValues(scalar) = doOp(entityValues(scalar), recvBuffer[recvOffset++]);
          }
        }
      }
      else if (f.host_data_layout() == Layout::Left) {
        auto fieldData = f.data<T, ReadWrite, stk::ngp::HostMemSpace, Layout::Left>();
        for (size_t i = 0; i < hostCommMapIndices.extent(0); ++i) {
          const unsigned bucketId = hostCommMapIndices(i).bucket_id;
          const unsigned bucketOrd = hostCommMapIndices(i).bucket_ord;
          auto entityValues = fieldData.entity_values(FastMeshIndex{bucketId, bucketOrd});
          for (ScalarIdx scalar : entityValues.scalars()) {
            entityValues(scalar) = doOp(entityValues(scalar), recvBuffer[recvOffset++]);
          }
        }
      }
      else {
        STK_ThrowErrorMsg("Unsupported host Field data layout: " << f.host_data_layout());
      }
    }
  };

  MPI_Comm comm = mesh.parallel();
  stk::parallel_data_exchange_sym_pack_unpack<T>(comm, commProcs, msgPacker, msgUnpacker, deterministic);
}

template <Operation OP>
inline
void parallel_op(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic)
{
  if (mesh.parallel_size() == 1 || fields.empty()) return;

  if (fields[0]->type_is<double>()) {
    parallel_op_impl<double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<int>()) {
    parallel_op_impl<int, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<float>()) {
    parallel_op_impl<float, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<long double>()) {
    parallel_op_impl<long double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<unsigned long>()) {
    parallel_op_impl<unsigned long, OP>(mesh, fields, deterministic);
  }
  else {
    STK_ThrowErrorMsg("Field is of unsupported type: " << fields[0]->data_traits().type_info.name());
  }
}

}

void parallel_sum(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic)
{
  parallel_op<Operation::SUM>(mesh, fields, deterministic);
}

//----------------------------------------------------------------------

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
void parallel_max(const BulkData& mesh, const std::vector<const FieldBase*>& fields)
{
  parallel_op<Operation::MAX>(mesh, fields, false);
}

/** Communicate and take the minimum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (minimum) field values on each sharing proc.
 */
void parallel_min(const BulkData& mesh, const std::vector<const FieldBase*>& fields)
{
  parallel_op<Operation::MIN>(mesh, fields, false);
}

template<typename T,Operation OP>
void parallel_op_including_ghosts_T(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic)
{
  using HostNgpField = stk::mesh::HostField<T,stk::ngp::HostMemSpace>;
  std::vector<HostNgpField*> ngpFields(fields.size());

  for(unsigned i=0; i<fields.size(); ++i) {
    STK_ThrowRequireMsg(fields[i]->type_is<T>(), "parallel_<op>_including_ghosts requires that all fields have the same scalar data type.");
    ngpFields[i] = new HostNgpField(mesh, *fields[i]);
  }

  stk::mesh::HostMesh hostNgpMesh(mesh);
  stk::mesh::parallel_op_including_ghosts_device_mpi<OP>(hostNgpMesh, ngpFields, deterministic);

  for(HostNgpField* ngpField : ngpFields) {
    delete ngpField;
  }
}

template <Operation OP>
void parallel_op_including_ghosts_impl(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic)
{
  if ( fields.empty() ) { return; }
  mesh.confirm_host_mesh_is_synchronized_from_device();

  if (fields[0]->type_is<long double>()) {
    parallel_op_including_ghosts_T<long double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<double>()) {
    parallel_op_including_ghosts_T<double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<float>()) {
    parallel_op_including_ghosts_T<float, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<int>()) {
    parallel_op_including_ghosts_T<int, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<unsigned long>()) {
    parallel_op_including_ghosts_T<unsigned long, OP>(mesh, fields, deterministic);
  }
  else {
    STK_ThrowRequireMsg(false, "Error, parallel_op only operates on fields of type long double, double, float, int, or unsigned long.");
  }
}

void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic)
{
  parallel_op_including_ghosts_impl<Operation::SUM>(mesh, fields, deterministic);
}

void parallel_max_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic)
{
  parallel_op_including_ghosts_impl<Operation::MAX>(mesh, fields, deterministic);
}

void parallel_min_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic)
{
  parallel_op_including_ghosts_impl<Operation::MIN>(mesh, fields, deterministic);
}


} // namespace mesh
} // namespace stk
