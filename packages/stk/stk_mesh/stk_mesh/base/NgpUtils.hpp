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

#ifndef STK_MESH_NGPUTILS_HPP
#define STK_MESH_NGPUTILS_HPP

#include <stk_util/stk_config.h>
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Types.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <numeric>

namespace stk {
namespace mesh {

namespace impl {

template <typename T, template <typename> class NgpDebugger = DefaultNgpFieldSyncDebugger>
inline void internal_fence_no_sync_to_host(NgpField<T, NgpDebugger>& ngpField)
{
  ngpField.asyncCopyState.execSpace.fence();
}
}

inline void ngp_field_fence(MetaData& meta)
{
  auto fields = meta.get_fields();

  for(auto field : fields) {
    if(field->has_ngp_field()) {
      field->fence();
    }
  }
}

inline void require_ngp_mesh_rank_limit(const stk::mesh::MetaData& meta)
{
  const size_t maxNumRanks = stk::topology::NUM_RANKS;
  const size_t numRanks = meta.entity_rank_count();
  ThrowRequireMsg(numRanks <= maxNumRanks,
                "stk::mesh::NgpMesh: too many entity ranks ("<<numRanks
                <<"). Required to be less-or-equal stk::topology::NUM_RANKS");
}

inline stk::NgpVector<unsigned> get_bucket_ids(const stk::mesh::BulkData &bulk,
                                               stk::mesh::EntityRank rank,
                                               const stk::mesh::Selector &selector)
{
  const stk::mesh::BucketVector &buckets = bulk.get_buckets(rank, selector);
  stk::NgpVector<unsigned> bucketIds(buckets.size());
  for(size_t i=0; i<buckets.size(); i++)
    bucketIds[i] = buckets[i]->bucket_id();
  bucketIds.copy_host_to_device();
  return bucketIds;
}

template <typename ViewType>
void transpose_contiguous_device_data_into_buffer(stk::ngp::ExecSpace & execSpace, unsigned numEntitiesInBlock, unsigned numPerEntity,
                                                  ViewType & deviceView, ViewType & bufferView)
{
  const auto& rangePolicy = Kokkos::RangePolicy<stk::ngp::ExecSpace>(execSpace, 0, numEntitiesInBlock);

  Kokkos::parallel_for("transpose_contiguous_device_data_into_buffer", rangePolicy,
    KOKKOS_LAMBDA(const int& entityIdx) {
      for (unsigned i = 0; i < numPerEntity; i++) {
        bufferView(entityIdx, i) = deviceView(ORDER_INDICES(entityIdx, i));
      }
    }
  );
}

template <typename ViewType>
void transpose_buffer_into_contiguous_device_data(stk::ngp::ExecSpace & execSpace, unsigned numEntitiesInBlock, unsigned numPerEntity,
                                                  ViewType & bufferView, ViewType & deviceView)
{
  const auto& rangePolicy = Kokkos::RangePolicy<stk::ngp::ExecSpace>(execSpace, 0, numEntitiesInBlock);

  Kokkos::parallel_for("transpose_buffer_into_contiguous_device_data", rangePolicy,
    KOKKOS_LAMBDA(const int& entityIdx) {
      for (unsigned i = 0; i < numPerEntity; i++) {
        deviceView(ORDER_INDICES(entityIdx, i)) = bufferView(entityIdx, i);
      }
    }
  );
}

template <typename DeviceViewType, typename BufferViewType, typename DeviceUnsignedViewType, typename ExecSpaceType>
void transpose_all_device_data_into_buffer(ExecSpaceType & execSpace,
                                           const stk::mesh::FieldBase & stkField,
                                           DeviceViewType & deviceView,
                                           BufferViewType & bufferView,
                                           DeviceUnsignedViewType & bucketSizes,
                                           DeviceUnsignedViewType & fieldBucketNumComponentsPerEntity)
{
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    size_t numBuckets = bucketSizes.extent(0);

    typedef typename Kokkos::TeamPolicy<ExecSpaceType, stk::ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
    Kokkos::parallel_for("transpose_all_device_data_into_buffer", teamPolicy,
                         KOKKOS_LAMBDA(const TeamHandleType & team) {
                           const unsigned bucketIndex = team.league_rank();
                           const unsigned bucketSize = bucketSizes(bucketIndex);
                           const unsigned numComponentsPerEntity = fieldBucketNumComponentsPerEntity(bucketIndex);
                           Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [&, bucketIndex, numComponentsPerEntity](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    bufferView(bucketIndex, entityIdx, i) =
                                                        deviceView(bucketIndex, ORDER_INDICES(entityIdx, i));
                                                  }
                                                });
                         });
}

template <typename DeviceViewType, typename BufferViewType, typename DeviceUnsignedViewType, typename ExecSpaceType>
void transpose_buffer_into_all_device_data(ExecSpaceType & execSpace,
                                           const stk::mesh::FieldBase & stkField,
                                           BufferViewType & bufferView,
                                           DeviceViewType & deviceView,
                                           DeviceUnsignedViewType & bucketSizes,
                                           DeviceUnsignedViewType & fieldBucketNumComponentsPerEntity)
{
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    size_t numBuckets = bucketSizes.extent(0);

    typedef typename Kokkos::TeamPolicy<ExecSpaceType, stk::ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
    Kokkos::parallel_for("transpose_buffer_into_all_device_data", teamPolicy,
                         KOKKOS_LAMBDA(const TeamHandleType & team) {
                           const unsigned bucketIndex = team.league_rank();
                           const unsigned bucketSize = bucketSizes(bucketIndex);
                           const unsigned numComponentsPerEntity = fieldBucketNumComponentsPerEntity(bucketIndex);
                           Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [&, bucketIndex, numComponentsPerEntity](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    deviceView(bucketIndex, ORDER_INDICES(entityIdx, i)) =
                                                        bufferView(bucketIndex, entityIdx, i);
                                                  }
                                                });
                         });
}


}
}

#endif

