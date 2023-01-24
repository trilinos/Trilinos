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

#ifndef STK_MESH_NGPFIELD_AUX_HPP
#define STK_MESH_NGPFIELD_AUX_HPP

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {
namespace impl {

enum NgpFieldSyncMode {
  INVALID = 0,
  HOST_TO_DEVICE = 1,
  HOST_TO_DEVICE_ASYNC = 2,
  DEVICE_TO_HOST = 3,
  DEVICE_TO_HOST_ASYNC = 4
};

template <typename DeviceViewType, typename DeviceUnsignedViewType, typename ExecSpaceType>
void transpose_from_pinned_and_mapped_memory(ExecSpaceType & execSpace,
                                             FieldDataPointerDeviceViewType& ptrToPinnedAndMappedMemory,
                                             DeviceViewType & deviceView,
                                             DeviceUnsignedViewType & bucketSizes,
                                             DeviceUnsignedViewType & fieldBucketNumComponentsPerEntity)
{
  using ValueType = typename DeviceViewType::value_type;
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceView.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_from_pinned_and_mapped_memory", teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType & team) {
                         const unsigned bucketIndex = team.league_rank();
                         const unsigned bucketSize = bucketSizes(bucketIndex);
                         const unsigned numComponentsPerEntity = fieldBucketNumComponentsPerEntity(bucketIndex);

                         bool isScalar = (numComponentsPerEntity == 1);

                         if(isScalar) {
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [=](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * entityIdx + i;
                                                    ValueType* data = bucketPtr + offset;

                                                    deviceView(bucketIndex, ORDER_INDICES(entityIdx, i)) = *data;
                                                  }
                                                });
                         } else {
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numComponentsPerEntity),
                                                [=](const int& componentIdx) {
                                                  for (unsigned i = 0; i < bucketSize; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * i + componentIdx;
                                                    ValueType* data = bucketPtr + offset;

                                                    deviceView(bucketIndex, ORDER_INDICES(i, componentIdx)) = *data;
                                                  }
                                                });
                         }
                       });
}

template <typename DeviceViewType, typename DeviceUnsignedViewType, typename ExecSpaceType>
void transpose_new_and_modified_buckets_to_device(ExecSpaceType & execSpace,
                                                  FieldDataPointerDeviceViewType& ptrToPinnedAndMappedMemory,
                                                  DeviceViewType & deviceView,
                                                  DeviceUnsignedViewType & bucketSizes,
                                                  DeviceUnsignedViewType & fieldBucketNumComponentsPerEntity,
                                                  DeviceUnsignedViewType & bucketsMarkedModified)
{
  using ValueType = typename DeviceViewType::value_type;
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceView.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_from_pinned_and_mapped_memory", teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType & team) {
                         const unsigned bucketIndex = team.league_rank();
                         const unsigned bucketSize = bucketSizes(bucketIndex);
                         const unsigned numComponentsPerEntity = fieldBucketNumComponentsPerEntity(bucketIndex);

                         if(bucketsMarkedModified(bucketIndex) == 0) { return; }

                         bool isScalar = (numComponentsPerEntity == 1);

                         if(isScalar) {
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [=](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * entityIdx + i;
                                                    ValueType* data = bucketPtr + offset;

                                                    deviceView(bucketIndex, ORDER_INDICES(entityIdx, i)) = *data;
                                                  }
                                                });
                         } else {
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numComponentsPerEntity),
                                                [=](const int& componentIdx) {
                                                  for (unsigned i = 0; i < bucketSize; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * i + componentIdx;
                                                    ValueType* data = bucketPtr + offset;

                                                    deviceView(bucketIndex, ORDER_INDICES(i, componentIdx)) = *data;
                                                  }
                                                });
                         }
                       });
}

template <typename DeviceViewType, typename DeviceUnsignedViewType, typename ExecSpaceType>
void transpose_to_pinned_and_mapped_memory(ExecSpaceType & execSpace,
                                             FieldDataPointerDeviceViewType& ptrToPinnedAndMappedMemory,
                                             DeviceViewType & deviceView,
                                             DeviceUnsignedViewType & bucketSizes,
                                             DeviceUnsignedViewType & fieldBucketNumComponentsPerEntity)
{
  using ValueType = typename DeviceViewType::value_type;
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceView.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_to_zero_copy_pinned_memory", teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType & team) {

                         const unsigned bucketIndex = team.league_rank();
                         const unsigned bucketSize = bucketSizes(bucketIndex);
                         const unsigned numComponentsPerEntity = fieldBucketNumComponentsPerEntity(bucketIndex);

                         bool isScalar = (numComponentsPerEntity == 1);

                         if(isScalar) { 
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
                                                [=](const int& entityIdx) {
                                                  for (unsigned i = 0; i < numComponentsPerEntity; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * entityIdx + i;
                                                    ValueType* data = bucketPtr + offset;

                                                    *data = deviceView(bucketIndex, ORDER_INDICES(entityIdx, i));
                                                  }
                                                });
                         } else {
                          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numComponentsPerEntity),
                                                [=](const int& componentIdx) {
                                                  for (unsigned i = 0; i < bucketSize; ++i) {
                                                    ValueType* bucketPtr = reinterpret_cast<ValueType*>(ptrToPinnedAndMappedMemory(bucketIndex));
                                                    uint64_t offset = numComponentsPerEntity * i + componentIdx;
                                                    ValueType* data = bucketPtr + offset;

                                                    *data = deviceView(bucketIndex, ORDER_INDICES(i, componentIdx));
                                                  }
                                                });
                         }
                       });
}

}
}
}
#endif
