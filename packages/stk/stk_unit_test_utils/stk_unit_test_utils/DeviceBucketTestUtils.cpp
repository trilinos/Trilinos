#include "DeviceBucketTestUtils.hpp"

namespace stk {
namespace unit_test_util {

Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace>
get_entities(const DeviceBucket deviceBucket)
{
  auto policy = stk::ngp::RangePolicy<stk::ngp::ExecSpace>(0, deviceBucket.size());

  Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace> hostEntities("host_entities", deviceBucket.size());
  auto deviceEntities = Kokkos::create_mirror_view(stk::mesh::NgpMeshDefaultMemSpace{}, hostEntities);
  auto func = KOKKOS_LAMBDA (const int idx)
  {
    deviceEntities(idx) = deviceBucket[idx];
  };

  Kokkos::parallel_for("get_entities", policy, func);
  Kokkos::deep_copy(hostEntities, deviceEntities);

  return hostEntities;
}

void check_entities(const stk::mesh::Bucket& hostBucket, const DeviceBucket& deviceBucket)
{
  EXPECT_EQ(hostBucket.size(), deviceBucket.size());
  auto entitiesFromDevice = get_entities(deviceBucket);
  for (size_t i=0; i < hostBucket.size(); ++i)
  {
    EXPECT_EQ(hostBucket[i], entitiesFromDevice[i]);
  }
}

Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace>
get_connected_entities(const DeviceBucket& deviceBucket, stk::mesh::EntityRank rank,
                       unsigned entityOrdinal)
{
  constexpr size_t MAX_NUM_ENTITIES = 100;

  auto policy = stk::ngp::RangePolicy<stk::ngp::ExecSpace>(0, MAX_NUM_ENTITIES);
  Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace> hostEntities("host_entities", MAX_NUM_ENTITIES);
  auto deviceEntities = Kokkos::create_mirror_view(stk::mesh::NgpMeshDefaultMemSpace{}, hostEntities);

  auto func = KOKKOS_LAMBDA (const unsigned int idx)
  {
    stk::mesh::ConnectedEntities connectedEntities = deviceBucket.get_connected_entities(entityOrdinal, rank);
    deviceEntities(idx) = idx < connectedEntities.size() ? connectedEntities[idx] : stk::mesh::Entity();
  };

  Kokkos::parallel_for("get_connected_entities", policy, func);
  Kokkos::deep_copy(hostEntities, deviceEntities);
  size_t numEntities = 0;
  for (size_t i=0; i < MAX_NUM_ENTITIES; ++i)
  {
    if (hostEntities(i) != stk::mesh::Entity())
    {
      numEntities++;
    }
  }

  Kokkos::resize(hostEntities, numEntities);
  return hostEntities;
}

Kokkos::View<stk::mesh::Permutation*, stk::ngp::HostMemSpace>
get_connected_permutations(const DeviceBucket deviceBucket, stk::mesh::EntityRank rank,
                           unsigned entityOrdinal)
{
  constexpr size_t MAX_NUM_ENTITIES = 100;

  auto policy = stk::ngp::RangePolicy<stk::ngp::ExecSpace>(0, MAX_NUM_ENTITIES);
  Kokkos::View<stk::mesh::Permutation*, stk::ngp::HostMemSpace> hostPerms("host_perms", MAX_NUM_ENTITIES);
  auto devicePerms = Kokkos::create_mirror_view(stk::mesh::NgpMeshDefaultMemSpace{}, hostPerms);

  auto func = KOKKOS_LAMBDA (const unsigned int idx)
  {
    auto connectedPerms = deviceBucket.get_connected_permutations(entityOrdinal, rank);
    devicePerms(idx) = idx < connectedPerms.size() ? connectedPerms[idx] : stk::mesh::Permutation::INVALID_PERMUTATION;
  };

  Kokkos::parallel_for("get_connected_perms", policy, func);
  Kokkos::deep_copy(hostPerms, devicePerms);
  size_t numEntities = 0;
  for (size_t i=0; i < MAX_NUM_ENTITIES; ++i)
  {
    if (hostPerms(i) != stk::mesh::Permutation::INVALID_PERMUTATION)
    {
      numEntities++;
    }
  }

  Kokkos::resize(hostPerms, numEntities);
  return hostPerms;
}

void check_connected_entities(const stk::mesh::Bucket& hostBucket, const DeviceBucket& deviceBucket)
{
  stk::mesh::EntityRank endRank = stk::topology::ELEM_RANK;
  for (stk::mesh::EntityRank rank = stk::topology::BEGIN_RANK; rank <= endRank; ++rank)
  {
    if (rank == hostBucket.entity_rank())
    {
      continue;
    }

    for (size_t i=0; i < hostBucket.size(); ++i)
    {
      auto entitiesFromDevice = get_connected_entities(deviceBucket, rank, i);
      stk::mesh::ConnectedEntities entitiesFromHost = hostBucket.get_connected_entities(i, rank);
      EXPECT_EQ(entitiesFromDevice.size(), entitiesFromHost.size());
      for (size_t j=0; j < entitiesFromDevice.size(); ++j)
      {
        EXPECT_EQ(entitiesFromDevice(j), entitiesFromHost[j]);
      }

      if (hostBucket.has_permutation(rank))
      {
        auto permsFromDevice = get_connected_permutations(deviceBucket, rank, i);
        stk::mesh::Permutation const* hostPerm = hostBucket.begin_permutations(i, rank);
        EXPECT_EQ(permsFromDevice.size(), size_t(hostBucket.end_permutations(i, rank) - hostPerm));
        for (size_t j=0; j < permsFromDevice.size(); ++j)
        {
          EXPECT_EQ(permsFromDevice(j), hostPerm[j]);
        }
      }
    }
  }
}

void check_repo_counts(stk::mesh::impl::BucketRepository& hostBucketRepo,
                       stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace> deviceBucketRepo)
{
  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    EXPECT_EQ(hostBucketRepo.get_partitions(rank).size(), deviceBucketRepo.get_partitions(rank).size());
    EXPECT_EQ(hostBucketRepo.buckets(rank).size(), deviceBucketRepo.num_buckets(rank));
  }
}

void check_partition_attributes(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition)
{
  EXPECT_EQ(hostPartition.get_rank(), devicePartition.get_rank());
  EXPECT_EQ(hostPartition.num_buckets(), devicePartition.num_buckets());
}

void check_partition_ordinals(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition)
{
  auto deviceOrdinals = devicePartition.superset_part_ordinals();
  auto ordinalsFromDevice = Kokkos::create_mirror_view(deviceOrdinals);
  Kokkos::deep_copy(ordinalsFromDevice, deviceOrdinals);
  const std::vector<stk::mesh::PartOrdinal> hostOrdinals = hostPartition.get_legacy_partition_id();

  EXPECT_EQ(ordinalsFromDevice.extent(0), hostOrdinals.size());
  for (size_t i=0; i < hostOrdinals.size(); ++i)
  {
    EXPECT_EQ(ordinalsFromDevice(i), hostOrdinals[i]);
  }
}

void check_partition_bucket_attributes(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition)
{
  const std::vector<stk::mesh::PartOrdinal> hostOrdinals = hostPartition.get_legacy_partition_id();

  for (size_t i=0; i < hostPartition.num_buckets(); ++i)
  {
    const stk::mesh::Bucket* hostBucket = hostPartition.get_bucket(i);
    DeviceBucketWrapper deviceBucket = devicePartition.get_bucket(i);

    EXPECT_EQ(hostBucket->ngp_mesh_bucket_id(), deviceBucket->bucket_id());
    EXPECT_EQ(hostBucket->ngp_field_bucket_id(), deviceBucket->bucket_id()); //TODO: what is this?
    EXPECT_EQ(hostBucket->size(), deviceBucket->size());
    const unsigned* bucketOrdinals = hostBucket->superset_part_ordinals().first;
    size_t numOrdinals = hostBucket->superset_part_ordinals().second - bucketOrdinals;
    EXPECT_EQ(numOrdinals, hostOrdinals.size());
    for (size_t j=0; j < numOrdinals; ++j)
    {
      EXPECT_EQ(bucketOrdinals[j], hostOrdinals[j]);
    }
  }
}

}
}