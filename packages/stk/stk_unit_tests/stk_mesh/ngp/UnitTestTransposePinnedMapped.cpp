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
#include <stk_mesh/baseImpl/NgpFieldAux.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_ngp_test/ngp_test.hpp>

namespace {

constexpr unsigned bucketCapacity = 4;
constexpr unsigned numPerEntity = 3;
constexpr unsigned bytesPerBucket = sizeof(double)*bucketCapacity*numPerEntity;

class TestTranspose : public ::testing::Test
{
public:
  TestTranspose() {}

  ~TestTranspose()
  {
    for(unsigned char* ptr : rawBucketAllocations) {
      fieldDataAllocator.deallocate(ptr, bytesPerBucket);
    }
    rawBucketAllocations.clear();
  }

  double get_value(unsigned bktIndex, unsigned entityIndex, unsigned component)
  {
    return (bktIndex+1)*1000 + (entityIndex+1)*100 + (component+1);
  }

  void fill_gold_host_field_data(unsigned numBuckets)
  {
    for(unsigned bkt=0; bkt<numBuckets; ++bkt) {
      for(unsigned entityIndex=0; entityIndex<bucketCapacity; ++entityIndex) {
        for(unsigned comp=0; comp<numPerEntity; ++comp) {
          goldHostFieldData(bkt,ORDER_INDICES(entityIndex,comp)) = get_value(bkt,entityIndex,comp);
        }
      }
    }
  }

  void fill_bucket_pointers(unsigned numBuckets)
  {
    rawBucketAllocations.resize(numBuckets);
    for(unsigned i=0; i<numBuckets; ++i) {
      rawBucketAllocations[i] = fieldDataAllocator.allocate(bytesPerBucket);
      double* hostBucketPtr = reinterpret_cast<double*>(rawBucketAllocations[i]);

      for(unsigned entityIndex=0; entityIndex<bucketCapacity; ++entityIndex) {
        for(unsigned comp=0; comp<numPerEntity; ++comp) {
          hostBucketPtr[entityIndex*numPerEntity + comp] = get_value(i,entityIndex,comp);
        }
      }

      double* deviceBucketPtr = hostBucketPtr;

#ifdef KOKKOS_ENABLE_CUDA
      cudaError_t status = cudaHostGetDevicePointer((void**)&deviceBucketPtr, (void*)hostBucketPtr, 0);
      ThrowRequireMsg(status == cudaSuccess, "Something went wrong during cudaHostGetDevicePointer: " + std::string(cudaGetErrorString(status)));
#elif defined(KOKKOS_ENABLE_HIP)
      hipError_t status = hipHostGetDevicePointer((void**)&deviceBucketPtr, (void*)hostBucketPtr, 0);
      ThrowRequireMsg(status == hipSuccess, "Something went wrong during hipHostGetDevicePointer: " + std::string(hipGetErrorString(status)));
#endif

      hostBucketPtrData(i) = reinterpret_cast<uintptr_t>(deviceBucketPtr);
    }

    Kokkos::deep_copy(deviceBucketPtrData, hostBucketPtrData);
  }

  void setup_views(unsigned numBuckets, double overallocationFactor)
  {
    deviceFieldData = stk::mesh::FieldDataDeviceViewType<double>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceFieldData"), numBuckets, ORDER_INDICES(bucketCapacity, numPerEntity));
    goldHostFieldData = stk::mesh::FieldDataHostViewType<double>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "goldHostFieldData"), numBuckets, ORDER_INDICES(bucketCapacity,numPerEntity));

    fill_gold_host_field_data(numBuckets);

    const unsigned numAlloc = std::lround(numBuckets*overallocationFactor);

    deviceBucketPtrData = stk::mesh::FieldDataPointerDeviceViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBktPtrData"), numAlloc);
    hostBucketPtrData = stk::mesh::FieldDataPointerHostViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "HostBktPtrData"), numAlloc);

    fill_bucket_pointers(numBuckets);

    deviceBucketSizes = stk::mesh::UnsignedViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketSizes"), numAlloc);
    hostBucketSizes = Kokkos::create_mirror_view(deviceBucketSizes);
    Kokkos::deep_copy(hostBucketSizes, bucketCapacity);
    Kokkos::deep_copy(deviceBucketSizes, hostBucketSizes);

    deviceBucketComponentsPerEntity = stk::mesh::UnsignedViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketComponentsPerEntity"), numAlloc);
    hostBucketComponentsPerEntity = Kokkos::create_mirror_view(deviceBucketComponentsPerEntity);
    Kokkos::deep_copy(hostBucketComponentsPerEntity, numPerEntity);
    Kokkos::deep_copy(deviceBucketComponentsPerEntity, hostBucketComponentsPerEntity);

    deviceBucketsMarkedModified = stk::mesh::BoolViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketsMarkedModified"), numAlloc);
    hostBucketsMarkedModified = Kokkos::create_mirror_view(deviceBucketsMarkedModified);
    Kokkos::deep_copy(hostBucketsMarkedModified, true);
    Kokkos::deep_copy(deviceBucketsMarkedModified, hostBucketsMarkedModified);
  }

  void check_field_data_values(unsigned numBuckets)
  {
    Kokkos::fence();
    auto checkFieldData = Kokkos::create_mirror_view(deviceFieldData);
    Kokkos::deep_copy(checkFieldData, deviceFieldData);

    for(unsigned bkt=0; bkt<numBuckets; ++bkt) {
      for(unsigned entityIndex=0; entityIndex<bucketCapacity; ++entityIndex) {
        for(unsigned comp=0; comp<numPerEntity; ++comp) {
          EXPECT_NEAR(goldHostFieldData(bkt,ORDER_INDICES(entityIndex,comp)), checkFieldData(bkt,ORDER_INDICES(entityIndex,comp)), 1.e-8)<<bkt<<":"<<entityIndex<<":"<<comp<<";";
        }
      }
    }
  }

  void check_host_field_data_values(unsigned numBuckets)
  {
    Kokkos::fence();
    for(unsigned bkt=0; bkt<numBuckets; ++bkt) {
      const double* hostBucketPtr = reinterpret_cast<const double*>(hostBucketPtrData(bkt));
      for(unsigned entityIndex=0; entityIndex<bucketCapacity; ++entityIndex) {
        for(unsigned comp=0; comp<numPerEntity; ++comp) {
          EXPECT_NEAR(goldHostFieldData(bkt,ORDER_INDICES(entityIndex,comp)), hostBucketPtr[entityIndex*numPerEntity+comp], 1.e-8);
        }
      }
    }
  }

protected:
  stk::mesh::AllocatorAdaptor<stk::impl::FieldDataAllocator<unsigned char>> fieldDataAllocator;
  std::vector<unsigned char*> rawBucketAllocations;

  stk::mesh::FieldDataPointerHostViewType  hostBucketPtrData;
  stk::mesh::FieldDataPointerDeviceViewType  deviceBucketPtrData;

  stk::mesh::FieldDataDeviceViewType<double>  deviceFieldData;
  stk::mesh::FieldDataHostViewType<double>  goldHostFieldData;

  stk::mesh::UnsignedViewType  deviceBucketSizes;
  stk::mesh::UnsignedViewType::HostMirror  hostBucketSizes;
  stk::mesh::UnsignedViewType  deviceBucketComponentsPerEntity;
  stk::mesh::UnsignedViewType::HostMirror  hostBucketComponentsPerEntity;
  stk::mesh::BoolViewType  deviceBucketsMarkedModified;
  stk::mesh::BoolViewType::HostMirror  hostBucketsMarkedModified;
};

TEST_F(TestTranspose, no_overallocation_transpose_from)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 2;
  const double overAllocFactor = 1.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  check_field_data_values(numBuckets);
}

TEST_F(TestTranspose, overallocate_2_transpose_from)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 3;
  const double overAllocFactor = 2.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  check_field_data_values(numBuckets);
}

TEST_F(TestTranspose, no_overallocation_transpose_new_and_modified)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 4;
  const double overAllocFactor = 1.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_new_and_modified_buckets_to_device(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity, deviceBucketsMarkedModified);
  check_field_data_values(numBuckets);
}

TEST_F(TestTranspose, overallocate_2_transpose_new_and_modified)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 200;
  const double overAllocFactor = 2.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_new_and_modified_buckets_to_device(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity, deviceBucketsMarkedModified);
  check_field_data_values(numBuckets);
}

TEST_F(TestTranspose, no_overallocation_transpose_to)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 100;
  const double overAllocFactor = 1.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  stk::mesh::impl::transpose_to_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  check_host_field_data_values(numBuckets);
}

TEST_F(TestTranspose, overallocate_2_transpose_to)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned numBuckets = 5;
  const double overAllocFactor = 2.0;
  setup_views(numBuckets, overAllocFactor);
  stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace();
  stk::mesh::impl::transpose_from_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  stk::mesh::impl::transpose_to_pinned_and_mapped_memory(execSpace, deviceBucketPtrData, deviceFieldData, deviceBucketSizes, deviceBucketComponentsPerEntity);
  check_host_field_data_values(numBuckets);
}

}
