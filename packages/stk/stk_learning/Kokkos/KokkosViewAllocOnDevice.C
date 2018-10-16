#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_MemoryPool.hpp>
#include <stk_ngp/NgpSpaces.hpp>

#ifdef KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda Device;
  typedef Kokkos::CudaSpace MemSpace;
#else
  typedef Kokkos::Serial Device;
  typedef Kokkos::HostSpace MemSpace;
#endif

typedef Kokkos::View<double*, Kokkos::LayoutRight, MemSpace>   EntitiesViewType;
typedef Kokkos::View<double*, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryUnmanaged> UnmanagedEntitiesViewType;

struct Bucket {
  UnmanagedEntitiesViewType entities;
  EntitiesViewType::HostMirror hostEntities;
};

typedef Kokkos::View<Bucket*,Kokkos::LayoutRight, MemSpace> BucketsViewType;

struct Mesh {
  BucketsViewType buckets;
  BucketsViewType::HostMirror hostBuckets;

  Mesh(size_t totalSize, size_t sizePerAlloc, size_t entitiesPerBucket)
  : pool(ngp::ExecSpace::memory_space(), totalSize, sizePerAlloc, sizePerAlloc),
    bytesPerAlloc(sizePerAlloc),
    numEntitiesPerBucket(entitiesPerBucket)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const int& bucketIndex) const {
      buckets(bucketIndex).entities = UnmanagedEntitiesViewType(static_cast<double*>(pool.allocate(bytesPerAlloc)), numEntitiesPerBucket);
      for(size_t j=0; j<numEntitiesPerBucket; ++j) {
        buckets(bucketIndex).entities(j) = j+1;
      }
  }

  void fill_buckets(size_t numBuckets)
  {
    buckets = BucketsViewType("buckets",numBuckets);
    hostBuckets = Kokkos::create_mirror_view(buckets);
  
    Kokkos::parallel_for(numBuckets, *this);
  
    //copy outer view from device to host
    Kokkos::deep_copy(hostBuckets, buckets);
  
    //copy inner views from device to host
    for(size_t bucketIndex=0; bucketIndex<numBuckets; ++bucketIndex) {
      hostBuckets(bucketIndex).hostEntities = Kokkos::create_mirror_view(hostBuckets(bucketIndex).entities);
      Kokkos::deep_copy(hostBuckets(bucketIndex).hostEntities, hostBuckets(bucketIndex).entities);
    }
  }

  Kokkos::MemoryPool<ngp::ExecSpace> pool;
  size_t bytesPerAlloc;
  size_t numEntitiesPerBucket;
};

void kokkos_view_allocate_on_device() {

  const size_t numBuckets = 5;
  const size_t numEntitiesPerBucket = 8;
  const size_t totalSize = numBuckets * numEntitiesPerBucket*sizeof(double);
  const size_t bytesPerAlloc = numEntitiesPerBucket*sizeof(double);

  Mesh mesh(totalSize, bytesPerAlloc, numEntitiesPerBucket);

  mesh.fill_buckets(numBuckets);

  //now confirm that the device data came back to the host
  for(size_t bucketIndex=0; bucketIndex<numBuckets; ++bucketIndex) {
    EXPECT_EQ(numEntitiesPerBucket, mesh.hostBuckets(bucketIndex).entities.size());
    EXPECT_EQ(numEntitiesPerBucket, mesh.hostBuckets(bucketIndex).hostEntities.size());
    for(size_t j=0; j<numEntitiesPerBucket; ++j) {
      double expected = j+1;
      EXPECT_NEAR(expected, mesh.hostBuckets(bucketIndex).hostEntities(j), 1.e-8);
    }
  }
}

TEST(KokkosView, allocateOnDevice) {
  kokkos_view_allocate_on_device();
}

