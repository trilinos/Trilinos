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

inline stk::NgpVector<unsigned> get_bucket_sizes(const stk::mesh::BulkData &bulk,
                                                 stk::mesh::EntityRank rank,
                                                 const stk::mesh::Selector &selector)
{
  const stk::mesh::BucketVector &buckets = bulk.get_buckets(rank, selector);
  stk::NgpVector<unsigned> bucketSizes(buckets.size());
  for (size_t i = 0; i < buckets.size(); ++i) {
    bucketSizes[i] = buckets[i]->size();
  }
  bucketSizes.copy_host_to_device();
  return bucketSizes;
}

inline stk::NgpVector<unsigned> get_num_components_per_entity(const stk::mesh::BulkData &bulk,
                                                              const stk::mesh::FieldBase & field,
                                                              const stk::mesh::Selector &selector)
{
  const stk::mesh::BucketVector &buckets = bulk.get_buckets(field.entity_rank(), selector);
  stk::NgpVector<unsigned> numComponentsPerEntity(buckets.size());
  for (size_t i = 0; i < buckets.size(); ++i) {
    numComponentsPerEntity[i] = stk::mesh::field_scalars_per_entity(field, *buckets[i]);
  }
  numComponentsPerEntity.copy_host_to_device();
  return numComponentsPerEntity;
}

inline size_t get_max_ngp_field_allocation_bytes(const MetaData & meta) {
  const FieldVector & fields = meta.get_fields();
  return std::accumulate(fields.begin(), fields.end(), 0u, [](size_t maxFieldDataBytes, const FieldBase * field) {
    return std::max(maxFieldDataBytes, get_total_ngp_field_allocation_bytes(*field));
  });
}


template <typename ViewType>
void transpose_contiguous_device_data_into_buffer(unsigned numEntitiesInBlock, unsigned numPerEntity,
                                                  ViewType & deviceView, ViewType & bufferView)
{
  Kokkos::parallel_for(numEntitiesInBlock,
    KOKKOS_LAMBDA(const int& entityIdx) {
      for (unsigned i = 0; i < numPerEntity; i++) {
        bufferView(entityIdx, i) = deviceView(ORDER_INDICES(entityIdx, i));
      }
    }
  );
}

template <typename ViewType>
void transpose_buffer_into_contiguous_device_data(unsigned numEntitiesInBlock, unsigned numPerEntity,
                                                  ViewType & bufferView, ViewType & deviceView)
{
  Kokkos::parallel_for(numEntitiesInBlock,
    KOKKOS_LAMBDA(const int& entityIdx) {
      for (unsigned i = 0; i < numPerEntity; i++) {
        deviceView(ORDER_INDICES(entityIdx, i)) = bufferView(entityIdx, i);
      }
    }
  );
}

template <typename DeviceViewType, typename BufferViewType>
void transpose_all_device_data_into_buffer(const stk::mesh::FieldBase & stkField,
                                           DeviceViewType & deviceView,
                                           BufferViewType & bufferView)
{
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    stk::NgpVector<unsigned> bucketSizes = get_bucket_sizes(stkField.get_mesh(), stkField.entity_rank(), selector);
    stk::NgpVector<unsigned> bucketNumComponentsPerEntity = stk::mesh::get_num_components_per_entity(stkField.get_mesh(),
                                                                                                     stkField, selector);
    size_t numBuckets = bucketSizes.size();

    typedef Kokkos::TeamPolicy<stk::mesh::ExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::ExecSpace>(numBuckets, Kokkos::AUTO);
    Kokkos::parallel_for(teamPolicy,
                         KOKKOS_LAMBDA(const TeamHandleType & team) {
                           const unsigned bucketIndex = team.league_rank();
                           const unsigned bucketSize = bucketSizes.device_get(bucketIndex);
                           const unsigned numComponentsPerEntity = bucketNumComponentsPerEntity.device_get(bucketIndex);
                           Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [&, bucketIndex, numComponentsPerEntity](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    bufferView(bucketIndex, entityIdx, i) =
                                                        deviceView(bucketIndex, ORDER_INDICES(entityIdx, i));
                                                  }
                                                });
                         });
}

template <typename DeviceViewType, typename BufferViewType>
void transpose_buffer_into_all_device_data(const stk::mesh::FieldBase & stkField,
                                           BufferViewType & bufferView,
                                           DeviceViewType & deviceView)
{
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    stk::NgpVector<unsigned> bucketSizes = get_bucket_sizes(stkField.get_mesh(), stkField.entity_rank(), selector);
    stk::NgpVector<unsigned> bucketNumComponentsPerEntity = stk::mesh::get_num_components_per_entity(stkField.get_mesh(),
                                                                                                     stkField, selector);
    size_t numBuckets = bucketSizes.size();

    typedef Kokkos::TeamPolicy<stk::mesh::ExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::ExecSpace>(numBuckets, Kokkos::AUTO);
    Kokkos::parallel_for(teamPolicy,
                         KOKKOS_LAMBDA(const TeamHandleType & team) {
                           const unsigned bucketIndex = team.league_rank();
                           const unsigned bucketSize = bucketSizes.device_get(bucketIndex);
                           const unsigned numComponentsPerEntity = bucketNumComponentsPerEntity.device_get(bucketIndex);
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

