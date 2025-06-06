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

#include "stk_mesh/base/NgpTypes.hpp"
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {
namespace impl {

//------------------------------------------------------------------------------
template <typename ValueType, typename ExecSpaceType, typename DeviceFieldMetaDataViewType, typename FieldOrdinalType,
          typename DeviceBucketsModifiedType>
void transpose_modified_buckets_to_device(const ExecSpaceType& execSpace,
                                          const DeviceFieldMetaDataViewType& deviceFieldMetaData,
                                          const FieldOrdinalType& fieldIndex,
                                          const DeviceBucketsModifiedType& bucketsMarkedModified,
                                          Layout hostDataLayout)
{
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceFieldMetaData.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_new_and_modified_buckets_to_device", teamPolicy,
    KOKKOS_LAMBDA(const TeamHandleType & team) {
      const unsigned bucketIndex = team.league_rank();

      const DeviceFieldMetaData& fieldMetaData = deviceFieldMetaData(bucketIndex);
      ValueType* hostBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_hostData);
      ValueType* deviceBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_data);

      if (bucketsMarkedModified(fieldIndex, bucketIndex) == 0 ||
          hostBucketPtr == nullptr || deviceBucketPtr == nullptr) { return; }

      const unsigned bucketSize = fieldMetaData.m_bucketSize;
      const unsigned bucketCapacity = fieldMetaData.m_bucketCapacity;
      const unsigned numScalarsPerEntity = fieldMetaData.m_numComponentsPerEntity*fieldMetaData.m_numCopiesPerEntity;

      bool isScalarField = (numScalarsPerEntity == 1);

      if (isScalarField) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
          [=](const unsigned& entity) {
            deviceBucketPtr[entity] = hostBucketPtr[entity];
          });
      }
      else {
        if (hostDataLayout == Layout::Right) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& component) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                const ValueType* hostData = hostBucketPtr + entity*numScalarsPerEntity + component;
                ValueType* deviceData = deviceBucketPtr + component*bucketCapacity + entity;
                *deviceData = *hostData;
              }
            }
          );
        }
        else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& component) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                const ValueType* hostData = hostBucketPtr + component*bucketCapacity + entity;
                ValueType* deviceData = deviceBucketPtr + component*bucketCapacity + entity;
                *deviceData = *hostData;
              }
            }
          );
        }
      }
    }
  );
}

template <typename ValueType, typename DeviceFieldMetaDataViewType, typename ExecSpaceType>
void transpose_from_pinned_and_mapped_memory(const ExecSpaceType& execSpace,
                                             const DeviceFieldMetaDataViewType& deviceFieldMetaData,
                                             Layout hostDataLayout)
{
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;
  size_t numBuckets = deviceFieldMetaData.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_from_pinned_and_mapped_memory", teamPolicy,
    KOKKOS_LAMBDA(const TeamHandleType & team) {
      const unsigned bucketIndex = team.league_rank();

      const DeviceFieldMetaData& fieldMetaData = deviceFieldMetaData(bucketIndex);
      ValueType* hostBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_hostData);
      ValueType* deviceBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_data);

      if (hostBucketPtr == nullptr || deviceBucketPtr == nullptr) {
        return;  // Field not on this Bucket
      }

      const unsigned bucketSize = fieldMetaData.m_bucketSize;
      const unsigned bucketCapacity = fieldMetaData.m_bucketCapacity;
      const unsigned numScalarsPerEntity = fieldMetaData.m_numComponentsPerEntity*fieldMetaData.m_numCopiesPerEntity;

      bool isScalarField = (numScalarsPerEntity == 1);

      if (isScalarField) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
          [=](const unsigned& entity) {
            deviceBucketPtr[entity] = hostBucketPtr[entity];
          }
        );
      }
      else {
        if (hostDataLayout == Layout::Right) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& scalar) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                const ValueType* hostData = hostBucketPtr + entity*numScalarsPerEntity + scalar;
                ValueType* deviceData = deviceBucketPtr + scalar*bucketCapacity + entity;
                *deviceData = *hostData;
              }
            }
          );
        }
        else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& scalar) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                const ValueType* hostData = hostBucketPtr + scalar*bucketCapacity + entity;
                ValueType* deviceData = deviceBucketPtr + scalar*bucketCapacity + entity;
                *deviceData = *hostData;
              }
            }
          );
        }
      }
    }
  );
}

template <typename ValueType, typename DeviceFieldMetaDataViewType, typename ExecSpaceType>
void transpose_to_pinned_and_mapped_memory(const ExecSpaceType& execSpace,
                                           const DeviceFieldMetaDataViewType& deviceFieldMetaData,
                                           Layout hostDataLayout)
{
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceFieldMetaData.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("transpose_to_zero_copy_pinned_memory", teamPolicy,
    KOKKOS_LAMBDA(const TeamHandleType & team) {
      const unsigned bucketIndex = team.league_rank();

      const DeviceFieldMetaData& fieldMetaData = deviceFieldMetaData(bucketIndex);
      ValueType* hostBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_hostData);
      ValueType* deviceBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_data);

      if (hostBucketPtr == nullptr || deviceBucketPtr == nullptr) {
        return;  // Field not on this Bucket
      }

      const unsigned bucketSize = fieldMetaData.m_bucketSize;
      const unsigned bucketCapacity = fieldMetaData.m_bucketCapacity;
      const unsigned numScalarsPerEntity = fieldMetaData.m_numComponentsPerEntity*fieldMetaData.m_numCopiesPerEntity;

      bool isScalarField = (numScalarsPerEntity == 1);

      if (isScalarField) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
          [=](const unsigned& entity) {
            hostBucketPtr[entity] = deviceBucketPtr[entity];
          });
      }
      else {
        if (hostDataLayout == Layout::Right) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& component) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                ValueType* hostData = hostBucketPtr + entity*numScalarsPerEntity + component;
                const ValueType* deviceData = deviceBucketPtr + component*bucketCapacity + entity;
                *hostData = *deviceData;
              }
            }
          );
        }
        else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
            [=](const int& component) {
              for (unsigned entity = 0; entity < bucketSize; ++entity) {
                ValueType* hostData = hostBucketPtr + component*bucketCapacity + entity;
                const ValueType* deviceData = deviceBucketPtr + component*bucketCapacity + entity;
                *hostData = *deviceData;
              }
            }
          );
        }
      }
    });
}

template <typename ValueType, typename DeviceFieldMetaDataViewType, typename ExecSpaceType>
void fill_field_with_value(const ExecSpaceType& execSpace,
                           const DeviceFieldMetaDataViewType& deviceFieldMetaData,
                           ValueType value)
{
  using TeamHandleType = typename stk::ngp::TeamPolicy<ExecSpaceType>::member_type;

  size_t numBuckets = deviceFieldMetaData.extent(0);

  const auto& teamPolicy = stk::ngp::TeamPolicy<ExecSpaceType>(execSpace, numBuckets, Kokkos::AUTO);
  Kokkos::parallel_for("fill_field_with_value", teamPolicy,
    KOKKOS_LAMBDA(const TeamHandleType & team) {
      const unsigned bucketIndex = team.league_rank();

      const DeviceFieldMetaData& fieldMetaData = deviceFieldMetaData(bucketIndex);
      ValueType* deviceBucketPtr = reinterpret_cast<ValueType*>(fieldMetaData.m_data);

      if (deviceBucketPtr == nullptr) {
        return;  // Field not on this Bucket
      }

      const unsigned bucketSize = fieldMetaData.m_bucketSize;
      const unsigned bucketCapacity = fieldMetaData.m_bucketCapacity;
      const unsigned numScalarsPerEntity = fieldMetaData.m_numComponentsPerEntity*fieldMetaData.m_numCopiesPerEntity;

      bool isScalarField = (numScalarsPerEntity == 1);

      if (isScalarField) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, bucketSize),
          [=](const unsigned& entity) {
            deviceBucketPtr[entity] = value;
          });
      }
      else {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numScalarsPerEntity),
          [=](const int& scalar) {
            for (unsigned entity = 0; entity < bucketSize; ++entity) {
              ValueType* deviceData = deviceBucketPtr + scalar*bucketCapacity + entity;
              *deviceData = value;
            }
          });
      }
    });
}

}
}
}
#endif
