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

#include <gtest/gtest.h>
#include <stk_mesh/base/NgpTypes.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/DeviceFieldDataManager.hpp>
#include <stk_mesh/baseImpl/NgpFieldAux.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/FieldDataAllocator.hpp>
#include <stk_ngp_test/ngp_test.hpp>

namespace {

constexpr unsigned bucketCapacity = 4;
constexpr unsigned numScalarsPerEntity = 3;
constexpr unsigned bytesPerBucket = sizeof(double)*bucketCapacity*numScalarsPerEntity;

template <typename NgpMemSpace>
using DeviceBucketRawDataType = Kokkos::View<double**, Kokkos::LayoutLeft, NgpMemSpace>;
template <typename NgpMemSpace>
using DeviceBucketRawDataCollectionType = Kokkos::View<DeviceBucketRawDataType<NgpMemSpace>*, stk::ngp::HostExecSpace>;

using HostBucketRawDataType = Kokkos::View<double**, Kokkos::LayoutLeft, stk::ngp::HostPinnedSpace>;
using HostBucketRawDataCollectionType = Kokkos::View<HostBucketRawDataType*, stk::ngp::HostExecSpace>;

class TestTranspose : public ::testing::Test
{
public:
  TestTranspose() = default;
  ~TestTranspose() = default;

  double get_value(unsigned bktIndex, unsigned entityIndex, unsigned component)
  {
    return (bktIndex+1)*1000 + (entityIndex+1)*100 + (component+1);
  }

  void fill_gold_host_field_data(unsigned numBuckets)
  {
    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      for (unsigned entityIndex = 0; entityIndex < bucketCapacity; ++entityIndex) {
        for (unsigned component = 0; component < numScalarsPerEntity; ++component) {
          goldHostFieldData(bucketId)(entityIndex, component) = get_value(bucketId, entityIndex, component);
        }
      }
    }
  }

  void fill_host_field_data(unsigned numBuckets)
  {
    rawHostBucketAllocations.resize(numBuckets);
    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      rawHostBucketAllocations[bucketId] = fieldDataAllocator.host_allocate(bytesPerBucket);
      ASAN_UNPOISON_MEMORY_REGION(rawHostBucketAllocations[bucketId].data(), bytesPerBucket);
      double* hostBucketPtr = reinterpret_cast<double*>(rawHostBucketAllocations[bucketId].data());

      for (unsigned entityIndex = 0; entityIndex < bucketCapacity; ++entityIndex) {
        for (unsigned component = 0; component < numScalarsPerEntity; ++component) {
          hostBucketPtr[entityIndex*numScalarsPerEntity + component] = get_value(bucketId, entityIndex, component);
        }
      }
    }
  }

  std::byte* get_host_bucket_pointer_for_device(unsigned bucketId) {
    return fieldDataAllocator.get_host_pointer_for_device(rawHostBucketAllocations[bucketId].data());
  }

  void fill_device_field_meta_data(unsigned numBuckets)
  {
    deviceFieldMetaData = stk::mesh::DeviceFieldMetaDataArrayType<stk::mesh::NgpMeshDefaultMemSpace>("deviceFieldMetaData",
                                                                                                     numBuckets);
    hostFieldMetaData = Kokkos::create_mirror_view(deviceFieldMetaData);

    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      auto& hostMetaData = hostFieldMetaData(bucketId);
      hostMetaData.m_data = reinterpret_cast<std::byte*>(deviceFieldData(bucketId).data());
      hostMetaData.m_hostData = get_host_bucket_pointer_for_device(bucketId);
      hostMetaData.m_numComponentsPerEntity = numScalarsPerEntity;
      hostMetaData.m_numCopiesPerEntity = 1;
      hostMetaData.m_bucketSize = bucketCapacity;
      hostMetaData.m_bucketCapacity = bucketCapacity;
    }

    Kokkos::deep_copy(deviceFieldMetaData, hostFieldMetaData);
  }

  void setup_views(unsigned numBuckets)
  {
    deviceFieldData = DeviceBucketRawDataCollectionType<stk::mesh::NgpMeshDefaultMemSpace>(
          Kokkos::view_alloc("deviceFieldDataCollection", Kokkos::SequentialHostInit), numBuckets);
    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      deviceFieldData[bucketId] = DeviceBucketRawDataType<stk::mesh::NgpMeshDefaultMemSpace>(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceFieldData"), bucketCapacity, numScalarsPerEntity);
    }

    goldHostFieldData = HostBucketRawDataCollectionType(
          Kokkos::view_alloc("goldHostFieldDataCollection", Kokkos::SequentialHostInit), numBuckets);
    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      goldHostFieldData[bucketId] = HostBucketRawDataType(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceFieldData"), bucketCapacity, numScalarsPerEntity);
    }

    fill_gold_host_field_data(numBuckets);
    fill_host_field_data(numBuckets);
    fill_device_field_meta_data(numBuckets);

    deviceBucketsMarkedModified = stk::mesh::DeviceBucketsModifiedCollectionType<stk::mesh::NgpMeshDefaultMemSpace>(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketsMarkedModified"), 1, numBuckets);
    hostBucketsMarkedModified = Kokkos::create_mirror_view(deviceBucketsMarkedModified);
    Kokkos::deep_copy(hostBucketsMarkedModified, true);
    Kokkos::deep_copy(deviceBucketsMarkedModified, hostBucketsMarkedModified);
  }

  void check_device_field_data_values(unsigned numBuckets)
  {
    Kokkos::fence();
    auto checkFieldData = Kokkos::create_mirror_view(deviceFieldData(0));

    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      Kokkos::deep_copy(checkFieldData, deviceFieldData(bucketId));

      for (unsigned entityIndex = 0; entityIndex < bucketCapacity; ++entityIndex) {
        for (unsigned component = 0; component < numScalarsPerEntity; ++component) {
          EXPECT_NEAR(goldHostFieldData(bucketId)(entityIndex, component),
                      checkFieldData(entityIndex, component), 1.e-8)
                      << bucketId << ":" << entityIndex << ":" << component << ";";
        }
      }
    }
  }

  void check_host_field_data_values(unsigned numBuckets)
  {
    Kokkos::fence();

    for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
      const double* hostBucketPtr = reinterpret_cast<const double*>(rawHostBucketAllocations[bucketId].data());
      for (unsigned entityIndex = 0; entityIndex < bucketCapacity; ++entityIndex) {
        for (unsigned component = 0; component < numScalarsPerEntity; ++component) {
          EXPECT_NEAR(goldHostFieldData(bucketId)(entityIndex, component),
                      hostBucketPtr[entityIndex*numScalarsPerEntity + component], 1.e-8);
        }
      }
    }
  }

protected:
  using AllocationType = stk::FieldDataAllocator<std::byte>::HostAllocationType;

  stk::FieldDataAllocator<std::byte> fieldDataAllocator;
  std::vector<AllocationType> rawHostBucketAllocations;

  DeviceBucketRawDataCollectionType<stk::mesh::NgpMeshDefaultMemSpace> deviceFieldData;
  HostBucketRawDataCollectionType goldHostFieldData;

  stk::mesh::DeviceFieldMetaDataArrayType<stk::mesh::NgpMeshDefaultMemSpace> deviceFieldMetaData;
  stk::mesh::DeviceFieldMetaDataArrayType<stk::mesh::NgpMeshDefaultMemSpace>::host_mirror_type hostFieldMetaData;

  stk::mesh::DeviceBucketsModifiedCollectionType<stk::mesh::NgpMeshDefaultMemSpace> deviceBucketsMarkedModified;
  stk::mesh::DeviceBucketsModifiedCollectionType<stk::mesh::NgpMeshDefaultMemSpace>::host_mirror_type hostBucketsMarkedModified;
};

TEST_F(TestTranspose, transpose_from)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 3;
  setup_views(numBuckets);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory<double>(execSpace, deviceFieldMetaData,
                                                                   stk::mesh::Layout::Right);
  check_device_field_data_values(numBuckets);
}

TEST_F(TestTranspose, transpose_modified_buckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 4;
  setup_views(numBuckets);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_modified_buckets_to_device<double>(execSpace, deviceFieldMetaData, 0,
                                                                deviceBucketsMarkedModified,
                                                                stk::mesh::Layout::Right);
  check_device_field_data_values(numBuckets);
}

TEST_F(TestTranspose, transpose_to)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 5;
  setup_views(numBuckets);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory<double>(execSpace, deviceFieldMetaData,
                                                                   stk::mesh::Layout::Right);
  stk::mesh::impl::transpose_to_pinned_and_mapped_memory<double>(execSpace, deviceFieldMetaData,
                                                                 stk::mesh::Layout::Right);
  check_host_field_data_values(numBuckets);
}

}
