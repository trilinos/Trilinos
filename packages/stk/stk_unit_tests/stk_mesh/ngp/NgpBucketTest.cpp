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

#include "stk_ngp_test/ngp_test.hpp"
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/TestHexFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/DevicePartition.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/StkNgpVector.hpp>

#ifdef STK_USE_DEVICE_MESH
class NgpBucketRepositoryTest : public stk::mesh::fixtures::TestHexFixture
{
public:
  NgpBucketRepositoryTest()
  {}

  void build_empty_mesh(const unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                        const unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity())
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
  }

  stk::mesh::Entity create_node(stk::mesh::EntityId nodeId, 
                                stk::mesh::PartVector const& initialParts = stk::mesh::PartVector{})
  {
    m_bulk->modification_begin();
    stk::mesh::Entity newNode = m_bulk->declare_node(nodeId, initialParts);
    m_bulk->modification_end();

    return newNode;
  }

  void create_parts(unsigned numParts)
  {
    for (unsigned i = 0; i < numParts; ++i) {
      m_meta->declare_part_with_topology("part" + std::to_string(i), testTopo);
    }
  }

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> get_device_part_ordinals(unsigned numPartOrdinals, unsigned partOrdinalOffset = 0)
  {
    stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> devicePartOrdinals("", numPartOrdinals);
    auto hostPartOrdinals = Kokkos::create_mirror_view(devicePartOrdinals);
    
    for (unsigned i = 0; i < numPartOrdinals; ++i) {
      auto newPart = m_meta->get_part("part" + std::to_string(i + partOrdinalOffset));
      hostPartOrdinals(i) = newPart->mesh_meta_data_ordinal();
    }

    Kokkos::deep_copy(devicePartOrdinals, hostPartOrdinals);
    return devicePartOrdinals;
  }

  void create_device_partitions_with_one_part(unsigned numPartitions, unsigned bucketCapacity = 1)
  {
    build_empty_mesh(bucketCapacity, bucketCapacity);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    create_parts(numPartitions);
    for (unsigned i = 0; i < numPartitions; ++i) {
      auto devicePartOrdinals = get_device_part_ordinals(1, i);
      deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
    }
  }

  void create_device_partitions_with_multi_part(unsigned numPartitions, unsigned numMultiPartPartitions, unsigned bucketCapacity = 1)
  {
    ASSERT_GE(numPartitions, numMultiPartPartitions);

    build_empty_mesh(bucketCapacity, bucketCapacity);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    create_parts(numPartitions + numMultiPartPartitions);

    for (unsigned i = 0; i < numPartitions - numMultiPartPartitions; ++i) {
      auto devicePartOrdinals = get_device_part_ordinals(1, i);
      deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
    }

    for (unsigned i = 0; i < numMultiPartPartitions; ++i) {
      auto devicePartOrdinals = get_device_part_ordinals(2, i);
      deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
    }
  }

  void create_one_device_partition()
  {
    build_empty_mesh(1, 1);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    create_parts(1);
    auto devicePartOrdinals = get_device_part_ordinals(1);
    deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  }

  void create_one_bucket_host_mesh()
  {
    build_empty_mesh(5, 5);

    stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
    stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);
    create_node(1, {&part1, &part2});
  }

  template <typename PartOrdinalView>
  void check_part_ordinal_equality(PartOrdinalView const& expected, PartOrdinalView const& found)
  {
    EXPECT_EQ(expected.extent(0), found.extent(0));
    if (expected.extent(0) != found.extent(0)) { return; }

    Kokkos::View<bool> result("");
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i) {
      for (unsigned i = 0; i < expected.extent(0); ++i) {
        if (expected(i) != found(i)) {
          result() = false;
          return;
        }
      }
      result() = true;
    });

    auto hostResult = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, result);
    EXPECT_TRUE(hostResult());
  }

  void create_one_partition_sparse_device_bucket_view(unsigned numBuckets)
  {
    build_empty_mesh(5, 5);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    create_parts(1);
    auto devicePartOrdinals = get_device_part_ordinals(1);
    auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);

    for (unsigned i = 0; i < numBuckets*2; ++i) {
      auto bucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
      partition->add_bucket(bucket);
    }

    for (unsigned i = 0; i < numBuckets*2; ++i) {
      if (i % 2 == 0) {
        auto bucket = deviceBucketRepo.get_bucket(testRank, i);
        partition->remove_bucket(bucket);
      }
    }
  }

  void create_two_partition_sparse_device_bucket_view(unsigned numBuckets)
  {
    build_empty_mesh(5, 5);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    unsigned numPartition = 2;
    create_parts(numPartition);
    
    for (unsigned j = 0; j < numPartition; j++) {
      auto devicePartOrdinals = get_device_part_ordinals(1, j);
      auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);

      for (unsigned i = 0; i < numBuckets*2; ++i) {
        auto bucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
        partition->add_bucket(bucket);
      }
    }

    unsigned initialNumBuckets = deviceBucketRepo.num_buckets(testRank);
    for (unsigned i = 0; i < initialNumBuckets; ++i) {
      if (i % 2 == 0) {
        auto bucket = deviceBucketRepo.get_bucket(testRank, i);
        auto partitionId = bucket->m_owningPartitionId;
        auto partition = deviceBucketRepo.get_partition(testRank, partitionId);
        partition->remove_bucket(bucket);
      }
    }
  }

  void create_three_partitions_one_empty_partition(unsigned numBuckets)
  {
    build_empty_mesh(5, 5);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    unsigned numPartition = 3;
    create_parts(numPartition);

    for (unsigned j = 0; j < numPartition; j++) {
      auto devicePartOrdinals = get_device_part_ordinals(1, j);
      auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);

      for (unsigned i = 0; i < numBuckets; ++i) {
        auto bucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
        partition->add_bucket(bucket);
      }
    }

    auto partitionToEmpty = deviceBucketRepo.get_partition(testRank, 1);
    unsigned initialNumBuckets = partitionToEmpty->num_buckets();

    for (unsigned i = 0; i < initialNumBuckets; ++i) {
      auto bucketPtr = partitionToEmpty->m_buckets(i).bucketPtr;
      partitionToEmpty->remove_bucket(bucketPtr);
    }
  }

  void check_device_bucket_repo_has_only_valid_buckets()
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    for (unsigned i = 0; i < deviceBucketRepo.num_buckets(testRank); ++i) {
      EXPECT_NE(stk::mesh::INVALID_BUCKET_ID, deviceBucketRepo.m_buckets[testRank](i).bucket_id());
    }
  }

  void check_device_bucket_repo_has_only_valid_partitions()
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    for (unsigned i = 0; i < deviceBucketRepo.num_partitions(testRank); ++i) {
      EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, deviceBucketRepo.m_partitions[testRank](i).partition_id());
    }
  }

  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData* m_meta;
  stk::mesh::EntityRank testRank = stk::topology::NODE_RANK;
  stk::topology::topology_t testTopo = stk::topology::NODE;
};

TEST_F(NgpBucketRepositoryTest, check_get_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_one_device_partition();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto& partitions = deviceBucketRepo.get_partitions(testRank);
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, partitions.extent(0));
}

TEST_F(NgpBucketRepositoryTest, create_empty_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  
  create_device_partitions_with_one_part(2);
  
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto& partitions = deviceBucketRepo.get_partitions(testRank);
  auto numPartitions = deviceBucketRepo.num_partitions(testRank);
  EXPECT_EQ(2u, numPartitions);

  auto& partition = partitions(0);
  EXPECT_EQ(0u, partition.num_buckets());
  EXPECT_EQ(0u, partition.num_entities());
}

TEST_F(NgpBucketRepositoryTest, check_invalid_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 1;
  create_device_partitions_with_one_part(numPartitions);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_ANY_THROW(deviceBucketRepo.get_partition(testRank, stk::mesh::impl::initialPartitionViewCapacity));
}

TEST_F(NgpBucketRepositoryTest, check_create_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 3;
  create_device_partitions_with_one_part(numPartitions);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto partitions = deviceBucketRepo.get_partitions(testRank);
 
  for (unsigned i = 0; i < numPartitions; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    check_part_ordinal_equality(devicePartOrdinals, partitions(i).superset_part_ordinals());
  }
}

TEST_F(NgpBucketRepositoryTest, check_search_for_partition_with_part_ordinals)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 3;
  create_device_partitions_with_multi_part(numPartitions, numPartitions);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  auto devicePartOrdinals = get_device_part_ordinals(2, 1);
  EXPECT_EQ(2u, devicePartOrdinals.extent(0));

  auto result = deviceBucketRepo.get_partition(testRank, devicePartOrdinals);
  EXPECT_NE(result, nullptr);

  auto resultPartOrdinals = result->superset_part_ordinals();
  EXPECT_EQ(2u, resultPartOrdinals.extent(0));

  check_part_ordinal_equality(devicePartOrdinals, resultPartOrdinals);
}

TEST_F(NgpBucketRepositoryTest, check_get_or_create_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh();
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  unsigned numPartitions = 2;
  create_parts(numPartitions);
  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto result = deviceBucketRepo.get_partition(testRank, devicePartOrdinals);
  EXPECT_EQ(result, nullptr);

  for (unsigned i = 0; i < numPartitions; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  }

  auto createdPartition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);
  EXPECT_NE(createdPartition, nullptr);

  result = deviceBucketRepo.get_partition(testRank, devicePartOrdinals);
  EXPECT_EQ(result, createdPartition);

  auto resultPartOrdinals = result->superset_part_ordinals();
  EXPECT_EQ(1u, resultPartOrdinals.extent(0));
  check_part_ordinal_equality(devicePartOrdinals, resultPartOrdinals);
}

TEST_F(NgpBucketRepositoryTest, device_bucket_repo_check_invalid_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(2);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  EXPECT_ANY_THROW(deviceBucketRepo.get_bucket(testRank, stk::mesh::impl::initialBucketViewCapacity));
}

TEST_F(NgpBucketRepositoryTest, device_bucket_repo_allocate_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);

  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

  EXPECT_NO_THROW(deviceBucketRepo.get_bucket(testRank, 0));
  EXPECT_NE(deviceBucketRepo.get_bucket(testRank, 0), nullptr);
}

TEST_F(NgpBucketRepositoryTest, device_bucket_repo_deallocate_bucket_no_sync_from_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto newBucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);

  EXPECT_EQ(1u, deviceBucketRepo.m_numBuckets[testRank]);
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

  EXPECT_NO_THROW(deviceBucketRepo.get_bucket(testRank, 0));
  EXPECT_NE(deviceBucketRepo.get_bucket(testRank, 0), nullptr);

  deviceBucketRepo.deallocate_bucket(newBucket);

  EXPECT_EQ(0u, deviceBucketRepo.m_numBuckets[testRank]);
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));
  EXPECT_EQ(stk::mesh::INVALID_BUCKET_ID, newBucket->bucket_id());
}

// INCOMPLETE: New entities are not yet added to device buckets
TEST_F(NgpBucketRepositoryTest, add_entity_with_parts_to_new_device_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(5, 5);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);
  create_node(1, {&part1, &part2});

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(0u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto entities = stk::mesh::get_entities(*m_bulk, testRank);
  EXPECT_EQ(1u, entities.size()); 

  stk::mesh::Part& newPart = m_meta->declare_part_with_topology("new part", testTopo);

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> devicePartOrdinals("", 1);
  auto hostPartOrdinals = Kokkos::create_mirror_view(stk::ngp::HostMemSpace{}, devicePartOrdinals);
  hostPartOrdinals(0) = newPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devicePartOrdinals, hostPartOrdinals);

  // EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

  // TODO Query buckets when entities are added to device buckets
}

// INCOMPLETE: New entities are not yet added to device buckets
TEST_F(NgpBucketRepositoryTest, add_entity_with_parts_to_existing_device_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(5, 5);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);
  create_node(1, {&part1, &part2});

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(0u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto entities = stk::mesh::get_entities(*m_bulk, testRank);
  EXPECT_EQ(1u, entities.size()); 

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> devicePartOrdinals("", 2);
  auto hostPartOrdinals = Kokkos::create_mirror_view(stk::ngp::HostMemSpace{}, devicePartOrdinals);
  hostPartOrdinals(0) = part1.mesh_meta_data_ordinal();
  hostPartOrdinals(1) = part2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devicePartOrdinals, hostPartOrdinals);

  // EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

  // TODO Query buckets when entities are added to device buckets
}

TEST_F(NgpBucketRepositoryTest, check_sync_bucket_ids_buckets_all_in_one_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 5;
  create_one_device_partition();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

  for (unsigned i = 0; i < numBuckets; ++i) {
    auto newBucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
    partition->add_bucket(newBucket);
  }

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));

  for (unsigned i = 0; i < numBuckets; ++i) {
    auto bucket = deviceBucketRepo.get_bucket(testRank, i);
    bucket->m_bucketId = 100 - i;
  }

  auto buckets = deviceBucketRepo.m_buckets[testRank];
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, buckets));

  deviceBucketRepo.sync_and_sort_bucket_ids(testRank);
  
  using UView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  UView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.m_numBuckets[testRank]);
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_buckets[testRank]));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));
}

TEST_F(NgpBucketRepositoryTest, check_sync_bucket_ids)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 5;
  create_device_partitions_with_one_part(numBuckets);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  for (unsigned i = 0; i < numBuckets; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);
    auto newBucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
    partition->add_bucket(newBucket);
  }

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));

  for (unsigned i = 0; i < numBuckets; ++i) {
    auto bucket = deviceBucketRepo.get_bucket(testRank, i);
    bucket->m_bucketId = 100 - i;
  }

  auto buckets = deviceBucketRepo.m_buckets[testRank];
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, buckets));

  deviceBucketRepo.sync_and_sort_bucket_ids(testRank);

  using UView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  UView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.m_numBuckets[testRank]);
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_buckets[testRank]));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));
}

TEST_F(NgpBucketRepositoryTest, check_sync_from_partitions_condense_bucket_view)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_one_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));

  deviceBucketRepo.sync_from_partitions();

  using UView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  UView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.m_numBuckets[testRank]);
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_buckets[testRank]));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));
}

TEST_F(NgpBucketRepositoryTest, check_sync_from_partitions_remove_empty_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 2;
  create_three_partitions_one_empty_partition(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets*2, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));

  deviceBucketRepo.sync_from_partitions();

  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_buckets[testRank]));
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_partitions[testRank]));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(testRank));

  using BucketUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  using PartitionUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DevicePartitionUView;
  BucketUView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.m_numBuckets[testRank]);
  PartitionUView compactPartitions(deviceBucketRepo.m_partitions[testRank].data(), deviceBucketRepo.m_numPartitions[testRank]);
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitions));

  check_device_bucket_repo_has_only_valid_partitions();
  check_device_bucket_repo_has_only_valid_buckets();
}

class NgpPartitionTest : public NgpBucketRepositoryTest
{
public:
  NgpPartitionTest() = default;

  void check_device_partition_has_only_valid_buckets() {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
    auto& partitions = deviceBucketRepo.m_partitions[testRank];
    auto numPartitions = deviceBucketRepo.num_partitions(testRank);

    for (unsigned i = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      auto numBuckets = partition.num_buckets();
      for (unsigned j = 0; j < numBuckets; ++j) {
        auto& bucket = partition.m_buckets(i);
        EXPECT_NE(stk::mesh::INVALID_BUCKET_ID, bucket.bucketPtr->bucket_id());
      }
    } 
  }
};

TEST_F(NgpPartitionTest, check_add_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto devicePartOrdinals = get_device_part_ordinals(1);

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
  auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  auto newBucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
  partition->add_bucket(newBucket);

  EXPECT_EQ(1u, partition->num_buckets());
  EXPECT_EQ(newBucket->bucket_id(), partition->m_buckets(0).bucketPtr->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, partition->partition_id());
}

TEST_F(NgpPartitionTest, check_remove_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto devicePartOrdinals = get_device_part_ordinals(1);

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
  auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  auto newBucket = deviceBucketRepo.allocate_bucket(testRank, devicePartOrdinals);
  partition->add_bucket(newBucket);

  EXPECT_EQ(1u, partition->num_buckets());
  EXPECT_EQ(newBucket->bucket_id(), partition->m_buckets(0).bucketPtr->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, partition->partition_id());

  partition->remove_bucket(newBucket);

  EXPECT_EQ(0u, partition->num_buckets());
  EXPECT_EQ(stk::mesh::INVALID_BUCKET_ID, newBucket->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, stk::mesh::INVALID_PARTITION_ID);
}

TEST_F(NgpPartitionTest, check_valid_bucket_view_after_sync_from_partitions_one_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_one_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

  EXPECT_EQ(numBuckets, partition->num_buckets());
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, partition->m_buckets.extent(0));

  deviceBucketRepo.sync_from_partitions();

  EXPECT_EQ(numBuckets, partition->num_buckets());
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, partition->m_buckets.extent(0));

  auto& buckets = deviceBucketRepo.m_buckets[testRank];
  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, buckets));

  check_device_bucket_repo_has_only_valid_buckets();
  check_device_partition_has_only_valid_buckets();
}

TEST_F(NgpPartitionTest, check_valid_bucket_view_after_sync_from_partitions_two_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_two_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets*2, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].extent(0));

  for (unsigned i = 0; i < 2; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

    EXPECT_EQ(numBuckets, partition->num_buckets());
    EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, partition->m_buckets.extent(0));
  }

  deviceBucketRepo.sync_from_partitions();

  for (unsigned i = 0; i < 2; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

    EXPECT_EQ(numBuckets, partition->num_buckets());
    EXPECT_EQ(stk::mesh::impl::initialBucketViewCapacity, partition->m_buckets.extent(0));

    auto& buckets = deviceBucketRepo.m_buckets[testRank];
    EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, buckets));

    check_device_bucket_repo_has_only_valid_buckets();
    check_device_partition_has_only_valid_buckets();
  }
}

class NgpBucketSearchTest : public NgpBucketRepositoryTest
{
public:
  NgpBucketSearchTest() = default;
};

TEST_F(NgpBucketSearchTest, search_partition_no_matching_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 5;
  create_device_partitions_with_one_part(numPartitions);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(5u, deviceBucketRepo.num_partitions(testRank));

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> testDevicePartOrdinals("", 2);
  auto hostPartOrdinals = Kokkos::create_mirror_view(testDevicePartOrdinals);
  hostPartOrdinals(0) = m_meta->get_part("part1")->mesh_meta_data_ordinal();
  hostPartOrdinals(1) = m_meta->get_part("part3")->mesh_meta_data_ordinal();
  Kokkos::deep_copy(testDevicePartOrdinals, hostPartOrdinals);

  auto partitions = deviceBucketRepo.get_partitions(testRank);
  auto policy = Kokkos::RangePolicy(0, deviceBucketRepo.num_partitions(testRank));

  auto searchResult = stk::mesh::impl::search_matching_device_partitions(policy, partitions, testDevicePartOrdinals);
  EXPECT_EQ(stk::mesh::INVALID_PARTITION_ID, searchResult);
}

TEST_F(NgpBucketSearchTest, search_partition_found_matching_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 5;
  create_device_partitions_with_multi_part(numPartitions, 2);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(5u, deviceBucketRepo.num_partitions(testRank));

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> testDevicePartOrdinals("", 2);
  auto hostPartOrdinals = Kokkos::create_mirror_view(testDevicePartOrdinals);
  hostPartOrdinals(0) = m_meta->get_part("part1")->mesh_meta_data_ordinal();
  hostPartOrdinals(1) = m_meta->get_part("part2")->mesh_meta_data_ordinal();
  Kokkos::deep_copy(testDevicePartOrdinals, hostPartOrdinals);

  auto partitions = deviceBucketRepo.get_partitions(testRank);
  auto policy = Kokkos::RangePolicy(0, deviceBucketRepo.num_partitions(testRank));

  auto searchResult = stk::mesh::impl::search_matching_device_partitions(policy, partitions, testDevicePartOrdinals);
  EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, searchResult);

  auto devicePartitionPartOrdinals = deviceBucketRepo.get_partition(testRank, searchResult)->superset_part_ordinals();
  check_part_ordinal_equality(testDevicePartOrdinals, devicePartitionPartOrdinals);
}
#endif
