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

#include "Kokkos_Core.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "stk_mesh/base/Part.hpp"
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
#include <string>

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

  void create_device_partitions_with_buckets(unsigned numPartitions, unsigned numBuckets)
  {
    build_empty_mesh(5, 5);
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    create_parts(numPartitions);

    for (unsigned j = 0; j < numPartitions; j++) {
      auto devicePartOrdinals = get_device_part_ordinals(1, j);
      auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);

      for (unsigned i = 0; i < numBuckets; ++i) {
        auto bucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
        partition->add_bucket(bucket);
      }
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

  void create_one_bucket_host_mesh(const unsigned initialBucketCapacity = 5,
                                   const unsigned maximumBucketCapacity = 5)
  {
    build_empty_mesh(initialBucketCapacity, maximumBucketCapacity);

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
      auto bucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
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
        auto bucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
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
        auto bucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
        partition->add_bucket(bucket);
      }
    }

    auto partitionToEmpty = deviceBucketRepo.get_partition(testRank, 1);
    unsigned initialNumBuckets = partitionToEmpty->num_buckets();

    for (unsigned i = 0; i < initialNumBuckets; ++i) {
      auto bucketPtr = partitionToEmpty->m_buckets[i].bucketPtr;
      partitionToEmpty->remove_bucket(bucketPtr);
    }
  }

  void check_device_bucket_repo_has_only_valid_buckets()
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    for (unsigned i = 0; i < deviceBucketRepo.num_buckets(testRank); ++i) {
      EXPECT_NE(stk::mesh::INVALID_BUCKET_ID, deviceBucketRepo.m_buckets[testRank][i].bucket_id());
    }
  }

  void check_device_bucket_repo_has_only_valid_partitions()
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    for (unsigned i = 0; i < deviceBucketRepo.num_partitions(testRank); ++i) {
      EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, deviceBucketRepo.m_partitions[testRank][i].partition_id());
    }
  }

  void check_device_bucket_and_partition_destruction_for_each_entity_run_rerf()
  {
    create_one_bucket_host_mesh(1, 1);

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

    stk::mesh::Selector universal = get_meta().universal_part();
    stk::mesh::for_each_entity_run(ngpMesh, testRank, universal, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const&) {});
  }

  void check_device_bucket_and_partition_destruction_for_each_entity_run_copy()
  {
    create_one_bucket_host_mesh(1, 1);

    stk::mesh::NgpMesh ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

    stk::mesh::Selector universal = get_meta().universal_part();
    stk::mesh::for_each_entity_run(ngpMesh, testRank, universal, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const&) {});
  }

  void check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_in_order_destruction()
  {
    build_empty_mesh(1, 1);
    auto numInitBucketCapacity = 32u;
    create_parts(numInitBucketCapacity+1);

    auto part = m_meta->get_part("part0");
    create_node(1, {part});

    stk::mesh::NgpMesh ngpMesh1 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo1 = ngpMesh1.get_device_bucket_repository();
    EXPECT_EQ(1u, deviceBucketRepo1.m_buckets[testRank].size());

    stk::mesh::Selector universal = get_meta().universal_part();
    stk::mesh::for_each_entity_run(ngpMesh1, testRank, universal, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const&) {});

    for (unsigned i = 1; i <= numInitBucketCapacity; ++i) {
      auto part = m_meta->get_part("part" + std::to_string(i));
      m_bulk->modification_begin();
      m_bulk->declare_node(i+1, stk::mesh::PartVector{part});
      m_bulk->modification_end();
    }

    stk::mesh::NgpMesh ngpMesh2 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo2 = ngpMesh2.get_device_bucket_repository();
    EXPECT_EQ(numInitBucketCapacity+1, deviceBucketRepo2.m_buckets[testRank].size());

    {
      stk::mesh::NgpMesh ngpMesh3 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      stk::mesh::for_each_entity_run(ngpMesh2, testRank, universal, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const&) {});
      stk::mesh::for_each_entity_run(ngpMesh3, testRank, universal, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const&) {});
    }
  }

  void check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_reverse_order_destruction()
  {
    build_empty_mesh(1, 1);
    auto numInitBucketCapacity = 32u;
    create_parts(numInitBucketCapacity*2);

    auto part = m_meta->get_part("part0");
    create_node(1, {part});
    stk::mesh::Selector universal = get_meta().universal_part();

    {
      stk::mesh::NgpMesh ngpMesh1 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo1 = ngpMesh1.get_device_bucket_repository();
      auto rank = testRank;

      EXPECT_EQ(1u, deviceBucketRepo1.m_buckets[rank].size());
      stk::mesh::EntityRank thisRank = rank;
      stk::mesh::for_each_entity_run(ngpMesh1, thisRank, universal,
        KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const& index) {
          [[maybe_unused]] const stk::mesh::Entity entity = ngpMesh1.get_entity(thisRank, index);
        }
      );
    }

    for (unsigned i = 1; i <= 11; ++i) {
      auto part = m_meta->get_part("part" + std::to_string(i));
      m_bulk->modification_begin();
      m_bulk->declare_node(i, stk::mesh::PartVector{part});
      m_bulk->modification_end();
    }

    {
      stk::mesh::NgpMesh ngpMesh2 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo2 = ngpMesh2.get_device_bucket_repository();
      EXPECT_EQ(11u, deviceBucketRepo2.m_buckets[testRank].size());
      stk::mesh::EntityRank thisRank = testRank;

      stk::mesh::for_each_entity_run(ngpMesh2, thisRank, universal,
        KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const& index) {
          [[maybe_unused]] const stk::mesh::Entity entity = ngpMesh2.get_entity(thisRank, index);
        }
      );
    }

    for (unsigned i = 12; i <= 17; ++i) {
      auto part = m_meta->get_part("part" + std::to_string(i));
      m_bulk->modification_begin();
      m_bulk->declare_node(i, stk::mesh::PartVector{part});
      m_bulk->modification_end();
    }

    {
      stk::mesh::NgpMesh ngpMesh3 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo3 = ngpMesh3.get_device_bucket_repository();
      EXPECT_EQ(17u, deviceBucketRepo3.m_buckets[testRank].size());

      for (unsigned i = 18; i <= numInitBucketCapacity+1; ++i) {
        auto part = m_meta->get_part("part" + std::to_string(i));
        m_bulk->modification_begin();
        m_bulk->declare_node(i, stk::mesh::PartVector{part});
        m_bulk->modification_end();
      }

      stk::mesh::NgpMesh ngpMesh4 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo4 = ngpMesh4.get_device_bucket_repository();
      EXPECT_EQ(numInitBucketCapacity+1, deviceBucketRepo4.m_buckets[testRank].size());
    }

    {
      stk::mesh::NgpMesh ngpMesh5 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      stk::mesh::EntityRank thisRank = testRank;
      stk::mesh::for_each_entity_run(ngpMesh5, testRank, universal,
        KOKKOS_LAMBDA(stk::mesh::FastMeshIndex const& index) {
          [[maybe_unused]] const stk::mesh::Entity entity = ngpMesh5.get_entity(thisRank, index);
        }
      );
    }
  }

  void check_device_bucket_copied_views_after_destroyed_from_copied_buckets() {
    build_empty_mesh(1, 1);
    auto numInitBucketCapacity = 32u;
    create_parts(numInitBucketCapacity*2);

    auto part = m_meta->get_part("part0");
    create_node(1, {part});

    {
      stk::mesh::NgpMesh ngpMesh1 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo1 = ngpMesh1.get_device_bucket_repository();

      EXPECT_EQ(1u, deviceBucketRepo1.num_buckets(testRank));
      EXPECT_EQ(0u, deviceBucketRepo1.m_buckets[testRank][0].bucket_id());
      EXPECT_EQ(1u, deviceBucketRepo1.m_buckets[testRank][0].m_entities.extent(0));
      EXPECT_EQ(1, deviceBucketRepo1.m_buckets[testRank][0].m_entities.use_count());

      Kokkos::parallel_for(1,
        KOKKOS_LAMBDA(const int) {
          auto& bucket = ngpMesh1.get_bucket(stk::topology::NODE_RANK, 0u);
          NGP_EXPECT_EQ(1u, bucket.m_entities(0).local_offset());
        }
      );

      for (unsigned i = 1; i <= numInitBucketCapacity; ++i) {
        auto part = m_meta->get_part("part" + std::to_string(i));
        m_bulk->modification_begin();
        m_bulk->declare_node(i+1, stk::mesh::PartVector{part});
        m_bulk->modification_end();
      }

      stk::mesh::NgpMesh ngpMesh2 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
      auto& deviceBucketRepo2 = ngpMesh2.get_device_bucket_repository();
      EXPECT_EQ(numInitBucketCapacity+1, deviceBucketRepo2.num_buckets(testRank));
      EXPECT_EQ(0u, deviceBucketRepo2.m_buckets[testRank][0].bucket_id());
      EXPECT_EQ(1u, deviceBucketRepo2.m_buckets[testRank][0].m_entities.extent(0));
      EXPECT_EQ(2, deviceBucketRepo2.m_buckets[testRank][0].m_entities.use_count());

      Kokkos::parallel_for(1,
        KOKKOS_LAMBDA(const int) {
          for (unsigned i = 0; i < ngpMesh2.num_buckets(stk::topology::NODE_RANK); ++i) {
            auto& bucket = ngpMesh2.get_bucket(stk::topology::NODE_RANK, i);
            NGP_EXPECT_EQ(i+1, bucket.m_entities(0).local_offset());
          }
        }
      );
    }

    stk::mesh::NgpMesh ngpMesh3 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo3 = ngpMesh3.get_device_bucket_repository();
    EXPECT_EQ(numInitBucketCapacity+1, deviceBucketRepo3.num_buckets(testRank));
    EXPECT_EQ(0u, deviceBucketRepo3.m_buckets[testRank][0].bucket_id());
    EXPECT_EQ(1u, deviceBucketRepo3.m_buckets[testRank][0].m_entities.extent(0));

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < ngpMesh3.num_buckets(stk::topology::NODE_RANK); ++i) {
          auto& bucket = ngpMesh3.get_bucket(stk::topology::NODE_RANK, i);
          NGP_EXPECT_EQ(i+1, bucket.m_entities(0).local_offset());
        }
      }
    );
  }

  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData* m_meta;
  stk::mesh::EntityRank testRank = stk::topology::NODE_RANK;
  stk::topology::topology_t testTopo = stk::topology::NODE;
};

TEST(NgpEntitySorting, Comparator)
{
  stk::mesh::impl::EntityCompareInvalidAtEnd compare;
  EXPECT_TRUE( compare(stk::mesh::Entity(1), stk::mesh::Entity(2)));
  EXPECT_TRUE( compare(stk::mesh::Entity(1), stk::mesh::Entity(stk::mesh::Entity::InvalidEntity)));
  EXPECT_FALSE(compare(stk::mesh::Entity(stk::mesh::Entity::InvalidEntity), stk::mesh::Entity(1)));
  EXPECT_TRUE( compare(stk::mesh::Entity(stk::mesh::Entity::MaxEntity), stk::mesh::Entity(stk::mesh::Entity::InvalidEntity)));
  EXPECT_FALSE(compare(stk::mesh::Entity(stk::mesh::Entity::InvalidEntity), stk::mesh::Entity(stk::mesh::Entity::MaxEntity)));
}

NGP_TEST_F(NgpBucketRepositoryTest, check_get_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_one_device_partition();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_NO_THROW(deviceBucketRepo.get_partitions(testRank));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(testRank));
}

NGP_TEST_F(NgpBucketRepositoryTest, create_empty_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  
  create_device_partitions_with_one_part(2);
  
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto& partitions = deviceBucketRepo.get_partitions(testRank);
  auto numPartitions = deviceBucketRepo.num_partitions(testRank);
  EXPECT_EQ(2u, numPartitions);

  auto& partition = partitions[0];
  EXPECT_EQ(0u, partition.num_buckets());
  EXPECT_EQ(0u, partition.num_entities());
}

NGP_TEST_F(NgpBucketRepositoryTest, check_invalid_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 1;
  create_device_partitions_with_one_part(numPartitions);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto invalidIdx = deviceBucketRepo.num_partitions(testRank)+1;
  EXPECT_ANY_THROW(deviceBucketRepo.get_partition(testRank, invalidIdx));
}

NGP_TEST_F(NgpBucketRepositoryTest, check_create_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numPartitions = 3;
  create_device_partitions_with_one_part(numPartitions);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto partitions = deviceBucketRepo.get_partitions(testRank);
 
  for (unsigned i = 0; i < numPartitions; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    check_part_ordinal_equality(devicePartOrdinals, partitions[i].superset_part_ordinals());
  }
}

NGP_TEST_F(NgpBucketRepositoryTest, check_search_for_partition_with_part_ordinals)
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

NGP_TEST_F(NgpBucketRepositoryTest, check_get_or_create_partition)
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

NGP_TEST_F(NgpBucketRepositoryTest, device_bucket_repo_check_invalid_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(2);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].size());
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto invalidIdx = deviceBucketRepo.num_buckets(testRank)+1;
  EXPECT_ANY_THROW(deviceBucketRepo.get_bucket(testRank, invalidIdx));
}

NGP_TEST_F(NgpBucketRepositoryTest, device_bucket_repo_construct_new_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].size());
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));
  EXPECT_NO_THROW(deviceBucketRepo.get_bucket(testRank, 0));
  EXPECT_NE(deviceBucketRepo.get_bucket(testRank, 0), nullptr);
}

NGP_TEST_F(NgpBucketRepositoryTest, device_bucket_repo_invalidate_bucket_no_sync_from_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank].size());
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto newBucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));
  EXPECT_NO_THROW(deviceBucketRepo.get_bucket(testRank, 0));
  EXPECT_NE(deviceBucketRepo.get_bucket(testRank, 0), nullptr);

  EXPECT_ANY_THROW(deviceBucketRepo.invalidate_bucket(newBucket));

  EXPECT_EQ(stk::mesh::INVALID_BUCKET_ID, newBucket->bucket_id());
}

// INCOMPLETE: New entities are not yet added to device buckets
NGP_TEST_F(NgpBucketRepositoryTest, add_entity_with_parts_to_new_device_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(5, 5);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);
  create_node(1, {&part1, &part2});

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

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
NGP_TEST_F(NgpBucketRepositoryTest, add_entity_with_parts_to_existing_device_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(5, 5);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);
  create_node(1, {&part1, &part2});

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(testRank));

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

NGP_TEST_F(NgpBucketRepositoryTest, check_sync_from_partitions_condense_bucket_view)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_one_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));

  for (unsigned i = 0; i < deviceBucketRepo.num_buckets(testRank); ++i) {
    EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank][i].size());
  }

  deviceBucketRepo.sync_from_partitions();

  using UView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  UView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.num_buckets(testRank));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
}

NGP_TEST_F(NgpBucketRepositoryTest, check_sync_from_partitions_remove_empty_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 2;
  create_three_partitions_one_empty_partition(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets*2, deviceBucketRepo.num_buckets(testRank));

  for (unsigned i = 0; i < deviceBucketRepo.num_buckets(testRank); ++i) {
    EXPECT_EQ(0u, deviceBucketRepo.m_buckets[testRank][i].size());
  }

  deviceBucketRepo.sync_from_partitions();
  EXPECT_EQ(0u, deviceBucketRepo.num_partitions(testRank));
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));

  using BucketUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  using PartitionUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DevicePartitionUView;
  BucketUView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.num_buckets(testRank));
  PartitionUView compactPartitions(deviceBucketRepo.m_partitions[testRank].data(), deviceBucketRepo.num_partitions(testRank));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitions));

  check_device_bucket_repo_has_only_valid_partitions();
  check_device_bucket_repo_has_only_valid_buckets();
}

void test_copy_to_device(stk::mesh::NgpMesh ngpMesh, stk::mesh::EntityRank testRank)
{
  EXPECT_EQ(3, ngpMesh.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {});
  EXPECT_EQ(3, ngpMesh.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
}

NGP_TEST_F(NgpBucketRepositoryTest, check_bucket_repo_ref_counts)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  stk::mesh::NgpMesh& ngpMeshInit = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
    EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
    EXPECT_EQ(0, deviceBucketRepo.m_buckets[testRank].get_view().use_count());
  }
  EXPECT_EQ(0, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());

  unsigned numBucketsToCreate = 33;
  for (unsigned i = 1; i <= numBucketsToCreate; ++i) {
    stk::mesh::Part& part = m_meta->declare_part_with_topology("testPart" + std::to_string(i), testTopo);

    m_bulk->modification_begin();
    m_bulk->declare_node(i, stk::mesh::PartVector{&part});
    m_bulk->modification_end();
  }
  EXPECT_EQ(0, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());

  // reference
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
    EXPECT_EQ(numBucketsToCreate, deviceBucketRepo.num_buckets(testRank));
    EXPECT_EQ(1, deviceBucketRepo.m_buckets[testRank].get_view().use_count());
    EXPECT_EQ(1, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
  }
  EXPECT_EQ(1, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());

  // copy
  {
    auto ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
    EXPECT_EQ(numBucketsToCreate, deviceBucketRepo.num_buckets(testRank));
    EXPECT_EQ(2, deviceBucketRepo.m_buckets[testRank].get_view().use_count());
    EXPECT_EQ(2, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
  }
  EXPECT_EQ(1, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());

  // copy into device lambda
  {
    EXPECT_EQ(1, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
    auto ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
    EXPECT_EQ(numBucketsToCreate, deviceBucketRepo.num_buckets(testRank));
    EXPECT_EQ(2, deviceBucketRepo.m_buckets[testRank].get_view().use_count());
    EXPECT_EQ(2, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());

    test_copy_to_device(ngpMesh, testRank);
    EXPECT_EQ(2, deviceBucketRepo.m_buckets[testRank].get_view().use_count());
  }
  EXPECT_EQ(1, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].get_view().use_count());
}

NGP_TEST_F(NgpBucketRepositoryTest, check_multiple_copies_of_ngp_mesh_updated_copies_requiring_bucket_resize)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);

  stk::mesh::NgpMesh ngpMeshInit = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMeshInit.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
  EXPECT_EQ(0u, ngpMeshInit.get_device_bucket_repository().num_buckets(testRank));

  unsigned numBucketsToCreate = 32;
  for (unsigned i = 1; i <= numBucketsToCreate; ++i) {
    stk::mesh::Part& part = m_meta->declare_part_with_topology("testPart" + std::to_string(i), testTopo);

    m_bulk->modification_begin();
    m_bulk->declare_node(i, stk::mesh::PartVector{&part});
    m_bulk->modification_end();
  }
  stk::mesh::NgpMesh ngpMesh2 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_EQ(0u, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].size());
  EXPECT_EQ(numBucketsToCreate, ngpMesh2.get_device_bucket_repository().m_buckets[testRank].size());

  EXPECT_EQ(0u, ngpMeshInit.get_device_bucket_repository().num_buckets(testRank));
  EXPECT_EQ(numBucketsToCreate, ngpMesh2.get_device_bucket_repository().num_buckets(testRank));

  {
    stk::mesh::Part& part = m_meta->declare_part_with_topology("testPart" + std::to_string(numBucketsToCreate+1), testTopo);
    m_bulk->modification_begin();
    m_bulk->declare_node(numBucketsToCreate+1, stk::mesh::PartVector{&part});
    m_bulk->modification_end();
  }

  stk::mesh::NgpMesh ngpMesh3 = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_EQ(0u, ngpMeshInit.get_device_bucket_repository().m_buckets[testRank].size());
  EXPECT_EQ(numBucketsToCreate, ngpMesh2.get_device_bucket_repository().m_buckets[testRank].size());
  EXPECT_EQ(numBucketsToCreate+1, ngpMesh3.get_device_bucket_repository().m_buckets[testRank].size());

  EXPECT_EQ(0u, ngpMeshInit.get_device_bucket_repository().num_buckets(testRank));
  EXPECT_EQ(numBucketsToCreate, ngpMesh2.get_device_bucket_repository().num_buckets(testRank));
  EXPECT_EQ(numBucketsToCreate+1, ngpMesh3.get_device_bucket_repository().num_buckets(testRank));
}


// NGP_TEST_F(NgpBucketRepositoryTest, check_updated_mesh_indices_after_first_bucket_removal)
// {
//   if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

//   unsigned numPartitions = 3;
//   unsigned numBuckets = 3;
//   create_device_partitions_with_buckets(numPartitions, numBuckets);

//   stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
//   auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

//   EXPECT_EQ(numBuckets*numPartitions, deviceBucketRepo.num_buckets(testRank));
//   EXPECT_EQ(numPartitions, deviceBucketRepo.num_partitions(testRank));

//   deviceBucketRepo.sync_from_partitions();

//   check_device_buckets_have_correct_fast_mesh_indices();

//   auto delBucket1 = deviceBucketRepo.get_bucket(testRank, 0);
//   auto delBucket2 = deviceBucketRepo.get_bucket(testRank, 3);

//   deviceBucketRepo.invalidate_bucket(delBucket1);
//   deviceBucketRepo.invalidate_bucket(delBucket2);

//   EXPECT_EQ(numBuckets*numPartitions-2, deviceBucketRepo.num_buckets(testRank));
//   EXPECT_EQ(numPartitions, deviceBucketRepo.num_partitions(testRank));
// }

// NGP_TEST_F(NgpBucketRepositoryTest, check_buckets_sorted_by_first_entity_identifier)
// {
//   if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

//   unsigned numBuckets = 2;
//   create_three_partitions_one_empty_partition(numBuckets);

//   stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
//   auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

//   EXPECT_EQ(numBuckets*2, deviceBucketRepo.num_buckets(testRank));
//   deviceBucketRepo.sync_from_partitions();

//   EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_buckets[testRank]));
//   EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, deviceBucketRepo.m_partitions[testRank]));
//   EXPECT_EQ(2u, deviceBucketRepo.num_partitions(testRank));

//   using BucketUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
//   using PartitionUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DevicePartitionUView;
//   BucketUView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.num_buckets(testRank));
//   PartitionUView compactPartitions(deviceBucketRepo.m_partitions[testRank].data(), deviceBucketRepo.num_partitions(testRank]);
//   EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));
//   EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitions));

//   check_device_bucket_repo_has_only_valid_partitions();
//   check_device_bucket_repo_has_only_valid_buckets();
// }

NGP_TEST_F(NgpBucketRepositoryTest, add_buckets_more_than_doubling_partition_capacity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  create_one_bucket_host_mesh(1, 1);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& hostBuckets = m_bulk->buckets(testRank);
  EXPECT_EQ(hostBuckets.size(), ngpMesh.num_buckets(testRank));

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto initDevBucketViewCapacity = 32u;

  EXPECT_EQ(1u, deviceBucketRepo.m_buckets[testRank].size());

  unsigned numNodeToAdd = initDevBucketViewCapacity-1;
  auto& parts = hostBuckets[0]->supersets();
  for (unsigned i = 1; i <= numNodeToAdd; ++i) {
    create_node(hostBuckets.size()+i, parts);
  }

  ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_EQ(initDevBucketViewCapacity, ngpMesh.num_buckets(testRank));
  EXPECT_EQ(initDevBucketViewCapacity, deviceBucketRepo.m_buckets[testRank].size());

  numNodeToAdd = initDevBucketViewCapacity+1;
  for (unsigned i = 1; i <= numNodeToAdd; ++i) {
    create_node(hostBuckets.size()+i, parts);
  }

  ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_EQ(initDevBucketViewCapacity*2+1, ngpMesh.num_buckets(testRank));
  EXPECT_EQ(initDevBucketViewCapacity*2+1, deviceBucketRepo.m_buckets[testRank].size());
}

NGP_TEST_F(NgpBucketRepositoryTest, add_host_buckets_more_than_initial_device_bucket_view_capacity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", testTopo);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", testTopo);

  unsigned numBucketsToAdd = 33u;
  for (unsigned i = 1; i <= numBucketsToAdd; ++i) {
    create_node(i, {&part1, &part2});
  }

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_EQ(numBucketsToAdd, ngpMesh.num_buckets(testRank));
}

NGP_TEST_F(NgpBucketRepositoryTest, check_device_bucket_and_partition_destruction_for_each_entity_run_ref)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  check_device_bucket_and_partition_destruction_for_each_entity_run_rerf();
}

NGP_TEST_F(NgpBucketRepositoryTest, check_device_bucket_and_partition_destruction_for_each_entity_run_copy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  check_device_bucket_and_partition_destruction_for_each_entity_run_copy();
}

NGP_TEST_F(NgpBucketRepositoryTest, check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_order1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_in_order_destruction();
}

NGP_TEST_F(NgpBucketRepositoryTest, check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_order2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  check_device_bucket_and_partition_destruction_for_each_entity_run_multiple_copies_reverse_order_destruction();
}

NGP_TEST_F(NgpBucketRepositoryTest, check_device_bucket_copied_views_after_destroyed_from_copied_buckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  check_device_bucket_copied_views_after_destroyed_from_copied_buckets();
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
      auto& partition = partitions[i];
      auto numBuckets = partition.num_buckets();
      for (unsigned j = 0; j < numBuckets; ++j) {
        auto& bucket = partition.m_buckets[i];
        EXPECT_NE(stk::mesh::INVALID_BUCKET_ID, bucket->bucket_id());
      }
    } 
  }
};

NGP_TEST_F(NgpPartitionTest, check_add_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto devicePartOrdinals = get_device_part_ordinals(1);

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
  auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  auto newBucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
  partition->add_bucket(newBucket);

  EXPECT_EQ(1u, partition->num_buckets());
  EXPECT_EQ(newBucket->bucket_id(), partition->m_buckets[0]->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, partition->partition_id());
}

NGP_TEST_F(NgpPartitionTest, check_remove_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_device_partitions_with_one_part(1);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto devicePartOrdinals = get_device_part_ordinals(1);

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(testRank));
  auto partition = deviceBucketRepo.create_partition(testRank, devicePartOrdinals);
  auto newBucket = deviceBucketRepo.construct_new_bucket(testRank, devicePartOrdinals);
  partition->add_bucket(newBucket);

  EXPECT_EQ(1u, partition->num_buckets());
  EXPECT_EQ(newBucket->bucket_id(), partition->m_buckets[0]->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, partition->partition_id());

  partition->remove_bucket(newBucket);

  EXPECT_EQ(0u, partition->num_buckets());
  EXPECT_EQ(stk::mesh::INVALID_BUCKET_ID, newBucket->bucket_id());
  EXPECT_EQ(newBucket->m_owningPartitionId, stk::mesh::INVALID_PARTITION_ID);
}

NGP_TEST_F(NgpPartitionTest, check_valid_bucket_view_after_sync_from_partitions_one_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_one_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets, deviceBucketRepo.num_buckets(testRank));

  auto devicePartOrdinals = get_device_part_ordinals(1);
  auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

  EXPECT_EQ(numBuckets, partition->num_buckets());

  deviceBucketRepo.sync_from_partitions();

  EXPECT_EQ(0u, partition->num_buckets());

  using BucketUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
  BucketUView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.num_buckets(testRank));
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));

  check_device_bucket_repo_has_only_valid_buckets();
  check_device_partition_has_only_valid_buckets();
}

NGP_TEST_F(NgpPartitionTest, check_valid_bucket_view_after_sync_from_partitions_two_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  unsigned numBuckets = 3;
  create_two_partition_sparse_device_bucket_view(numBuckets);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(numBuckets*2, deviceBucketRepo.num_buckets(testRank));

  for (unsigned i = 0; i < 2; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

    EXPECT_EQ(numBuckets, partition->num_buckets());
  }

  deviceBucketRepo.sync_from_partitions();

  for (unsigned i = 0; i < 2; ++i) {
    auto devicePartOrdinals = get_device_part_ordinals(1, i);
    auto partition = deviceBucketRepo.get_or_create_partition(testRank, devicePartOrdinals);

    EXPECT_EQ(0u, partition->num_buckets());
    EXPECT_EQ(0u, partition->m_buckets.size());

    using BucketUView = typename stk::mesh::impl::DeviceBucketRepository<stk::ngp::MemSpace>::DeviceBucketUView;
    BucketUView compactBuckets(deviceBucketRepo.m_buckets[testRank].data(), deviceBucketRepo.num_buckets(testRank));
    EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBuckets));

    check_device_bucket_repo_has_only_valid_buckets();
    check_device_partition_has_only_valid_buckets();
  }
}

class NgpBucketSearchTest : public NgpBucketRepositoryTest
{
public:
  NgpBucketSearchTest() = default;
};

NGP_TEST_F(NgpBucketSearchTest, search_partition_no_matching_partition)
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

NGP_TEST_F(NgpBucketSearchTest, search_partition_found_matching_partition)
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

class NgpBucketRepoBatchOpTest : public NgpBucketRepositoryTest
{
public:
  using PartOrdinalsProxyViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;

  NgpBucketRepoBatchOpTest() = default;

  template <typename DevicePartOrdinalView>
  PartOrdinalsProxyViewType create_part_ordinal_proxy_separate_partitions(unsigned numEntities, unsigned numPartsPerEntity, DevicePartOrdinalView const& devicePartOrdinalView)
  {
    auto rank = testRank;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), numEntities);
    
    Kokkos::parallel_for(numEntities, KOKKOS_LAMBDA(const int i) {
      partOrdinalsProxy(i).rank = rank;
      partOrdinalsProxy(i).startPtr = devicePartOrdinalView.data() + (numPartsPerEntity * i);
      partOrdinalsProxy(i).length = numPartsPerEntity;
    });

    return partOrdinalsProxy;
  }

  template <typename DevicePartOrdinalView>
  PartOrdinalsProxyViewType create_part_ordinal_proxy_same_partition(unsigned numEntities, unsigned numPartsPerEntity, DevicePartOrdinalView const& devicePartOrdinalView)
  {
    auto rank = testRank;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), numEntities);
    
    Kokkos::parallel_for(numEntities, KOKKOS_LAMBDA(const int i) {
      partOrdinalsProxy(i).rank = rank;
      partOrdinalsProxy(i).startPtr = devicePartOrdinalView.data();
      partOrdinalsProxy(i).length = numPartsPerEntity;
    });

    return partOrdinalsProxy;
  }

  template <typename DeviceEntityViewType>
  void fill_device_entities(DeviceEntityViewType& deviceEntityView)
  {
    [[maybe_unused]] stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    stk::mesh::EntityVector entities;
    m_bulk->get_entities(testRank, m_meta->universal_part(), entities);

    auto hostEntityView = Kokkos::create_mirror_view(deviceEntityView);

    for (unsigned i = 0; i < entities.size(); ++i) {
      hostEntityView(i) = entities[i];
    }
    Kokkos::deep_copy(deviceEntityView, hostEntityView);
  }

  stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> get_aggregated_part_ordinals_add_one_part(unsigned numEntities, unsigned numPartsPerEntity, stk::mesh::PartVector const& partToAdd)
  {
    stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace> devicePartOrdinals("", numEntities * (numPartsPerEntity + partToAdd.size()));
    auto hostPartOrdinals = Kokkos::create_mirror_view(devicePartOrdinals);

    int offset = 0;
    for (unsigned i = 0; i < numEntities; ++i) {
      for (unsigned j = 0; j < numPartsPerEntity; ++j) {
        auto part = m_meta->get_part("part" + std::to_string(numPartsPerEntity+j));
        hostPartOrdinals(offset++) = part->mesh_meta_data_ordinal();
      }
      for (unsigned j = 0; j < partToAdd.size(); j++) {
        hostPartOrdinals(offset++) = partToAdd[j]->mesh_meta_data_ordinal();
      }
    }

    Kokkos::deep_copy(devicePartOrdinals, hostPartOrdinals);
    return devicePartOrdinals;
  }

  void check_partition_part_ordinal_count(unsigned numPartsPerEntity)
  {
    auto rank = testRank;
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      auto& partitions = deviceBucketRepo.m_partitions[rank];

      for (unsigned i = 0; i < deviceBucketRepo.num_partitions(rank); ++i) {
        NGP_EXPECT_EQ(numPartsPerEntity, partitions[i].m_partOrdinals.extent(0));
      }
    });
  }
};

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_create_partition_separate_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh();

  auto numEntities = 2;
  auto numPartsPerEntity = 2;
  create_parts(numEntities * numPartsPerEntity);

  auto devicePartOrdinalView = get_device_part_ordinals(numEntities * numPartsPerEntity);
  auto devicePartOrdinalsProxy = create_part_ordinal_proxy_separate_partitions(numEntities, numPartsPerEntity, devicePartOrdinalView);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  deviceBucketRepo.batch_create_partitions(devicePartOrdinalsProxy);

  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(testRank));
  check_partition_part_ordinal_count(numPartsPerEntity);
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_create_partition_same_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh();

  auto numEntities = 2;
  auto numPartsPerEntity = 4;
  create_parts(numPartsPerEntity);

  auto devicePartOrdinalView = get_device_part_ordinals(numPartsPerEntity);
  auto devicePartOrdinalsProxy = create_part_ordinal_proxy_same_partition(numEntities, numPartsPerEntity, devicePartOrdinalView);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  deviceBucketRepo.batch_create_partitions(devicePartOrdinalsProxy);

  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(testRank));
  check_partition_part_ordinal_count(numPartsPerEntity);
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_get_partition_empty_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh();

  auto numEntities = 2;
  auto numPartsPerEntity = 1;
  create_parts(numEntities + numPartsPerEntity);

  stk::mesh::PartVector partToAdd1{m_meta->get_part("part0")};
  stk::mesh::PartVector partToAdd2{m_meta->get_part("part1")};
  stk::mesh::PartVector partToAdd3{m_meta->get_part("part2")};
  create_node(1, partToAdd1);
  create_node(2, partToAdd2);

  using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
  DeviceEntitiesType entities("deviceEntities", numEntities);

  fill_device_entities(entities);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  auto devicePartOrdinalView = get_aggregated_part_ordinals_add_one_part(numEntities, numPartsPerEntity, partToAdd3);
  auto devicePartOrdinalProxy = create_part_ordinal_proxy_separate_partitions(numEntities, numPartsPerEntity+1, devicePartOrdinalView);

  stk::mesh::UnsignedPairViewType<stk::ngp::MemSpace> srcDestPartitionsPerEntity("srcDestPartitionIdPerEntity", entities.size());
  Kokkos::Experimental::fill(stk::ngp::ExecSpace{}, srcDestPartitionsPerEntity, Kokkos::pair{stk::mesh::INVALID_PARTITION_ID, stk::mesh::INVALID_PARTITION_ID});
  deviceBucketRepo.batch_get_partitions(entities, devicePartOrdinalProxy, srcDestPartitionsPerEntity);

  auto hostSrcDestPartitions = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, srcDestPartitionsPerEntity);
  for (unsigned i = 0; i < hostSrcDestPartitions.extent(0); ++i) {
    EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, hostSrcDestPartitions(i).first);
    EXPECT_EQ(stk::mesh::INVALID_PARTITION_ID, hostSrcDestPartitions(i).second);
  }
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_get_partition_found_partitions)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  using DevicePartOrdinalView = stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace>;
  using HostPartOrdinalView = typename DevicePartOrdinalView::host_mirror_type;

  build_empty_mesh();

  auto numEntities = 2;
  auto numPartsPerEntity = 1;
  create_parts(numEntities + numPartsPerEntity);

  stk::mesh::PartVector partToAdd1{m_meta->get_part("part0")};
  stk::mesh::PartVector partToAdd2{m_meta->get_part("part1")};
  stk::mesh::PartVector partToAdd3{m_meta->get_part("part2")};
  create_node(1, partToAdd1);
  create_node(2, partToAdd2);

  using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
  DeviceEntitiesType entities("deviceEntities", numEntities);

  fill_device_entities(entities);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  HostPartOrdinalView hostPartOrdinal1("", 2);
  hostPartOrdinal1(0) = m_meta->get_part("part0")->mesh_meta_data_ordinal();
  hostPartOrdinal1(1) = m_meta->get_part("part2")->mesh_meta_data_ordinal(); 
  auto devicePartOrdinal1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPartOrdinal1);
  deviceBucketRepo.create_partition(testRank, devicePartOrdinal1);

  HostPartOrdinalView hostPartOrdinal2("", 2);
  hostPartOrdinal2(0) = m_meta->get_part("part1")->mesh_meta_data_ordinal();
  hostPartOrdinal2(1) = m_meta->get_part("part2")->mesh_meta_data_ordinal(); 
  auto devicePartOrdinal2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPartOrdinal2);
  deviceBucketRepo.create_partition(testRank, devicePartOrdinal2);

  auto devicePartOrdinalView = get_aggregated_part_ordinals_add_one_part(numEntities, numPartsPerEntity, partToAdd3);
  auto devicePartOrdinalProxy = create_part_ordinal_proxy_separate_partitions(numEntities, numPartsPerEntity+1, devicePartOrdinalView);

  stk::mesh::UnsignedPairViewType<stk::ngp::MemSpace> srcDestPartitionsPerEntity("srcDestPartitionIdPerEntity", entities.size());
  Kokkos::Experimental::fill(stk::ngp::ExecSpace{}, srcDestPartitionsPerEntity, Kokkos::pair{stk::mesh::INVALID_PARTITION_ID, stk::mesh::INVALID_PARTITION_ID});
  deviceBucketRepo.batch_get_partitions(entities, devicePartOrdinalProxy, srcDestPartitionsPerEntity);

  auto hostSrcDestPartitions = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, srcDestPartitionsPerEntity);
  for (unsigned i = 0; i < hostSrcDestPartitions.extent(0); ++i) {
    EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, hostSrcDestPartitions(i).first);
    EXPECT_NE(stk::mesh::INVALID_PARTITION_ID, hostSrcDestPartitions(i).second);
  }
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_move_entities_same_partition)
{
  using DevicePartOrdinalView = stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace>;
  using HostPartOrdinalView = typename DevicePartOrdinalView::host_mirror_type;
  using DevicePairViewType = stk::mesh::UnsignedPairViewType<stk::ngp::MemSpace>;
  using HostPairViewType = typename DevicePairViewType::host_mirror_type;
  using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;

  build_empty_mesh();

  create_parts(2);
  auto part0 = m_meta->get_part("part0");
  auto part1 = m_meta->get_part("part1");
  create_node(1, stk::mesh::PartVector{part0});
  create_node(2, stk::mesh::PartVector{part1});
  
  DeviceEntitiesType entities("deviceEntities", 2);
  fill_device_entities(entities);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  HostPartOrdinalView hostPartOrdinal("", 2);
  hostPartOrdinal(0) = part0->mesh_meta_data_ordinal();
  hostPartOrdinal(1) = part1->mesh_meta_data_ordinal();

  auto devicePartOrdinal = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPartOrdinal);
  deviceBucketRepo.create_partition(testRank, devicePartOrdinal);
  EXPECT_EQ(3u, deviceBucketRepo.num_partitions(testRank));

  HostPairViewType hostPairs("", 2);
  hostPairs(0) = {0, 2};
  hostPairs(1) = {1, 2};

  auto srcDestPartitions = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPairs);

  deviceBucketRepo.batch_move_entities(entities, srcDestPartitions);

  // empty buckets (with stale data) have not been cleared out from partition 0 and 1 at this point
  auto partition0 = deviceBucketRepo.get_partition(testRank, 0);
  EXPECT_EQ(0u, partition0->m_numEntities);
  EXPECT_EQ(1u, partition0->num_buckets());
  auto partition1 = deviceBucketRepo.get_partition(testRank, 1);
  EXPECT_EQ(0u, partition1->m_numEntities);
  EXPECT_EQ(1u, partition1->num_buckets());
  auto partition2 = deviceBucketRepo.get_partition(testRank, 2);
  EXPECT_EQ(2u, partition2->m_numEntities);
  EXPECT_EQ(1u, partition2->num_buckets());
  EXPECT_EQ(2u, partition2->m_buckets[0]->m_bucketSize);
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_move_entities_separate_partitions)
{
  using DevicePartOrdinalView = stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace>;
  using HostPartOrdinalView = typename DevicePartOrdinalView::host_mirror_type;
  using DevicePairViewType = stk::mesh::UnsignedPairViewType<stk::ngp::MemSpace>;
  using HostPairViewType = typename DevicePairViewType::host_mirror_type;
  using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;

  build_empty_mesh();

  create_parts(3);
  auto part0 = m_meta->get_part("part0");
  auto part1 = m_meta->get_part("part1");
  auto part2 = m_meta->get_part("part2");
  create_node(1, stk::mesh::PartVector{part0});
  create_node(2, stk::mesh::PartVector{part1});
  
  DeviceEntitiesType entities("deviceEntities", 2);
  fill_device_entities(entities);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  HostPartOrdinalView hostPartOrdinal1("", 2);
  hostPartOrdinal1(0) = part0->mesh_meta_data_ordinal();
  hostPartOrdinal1(1) = part2->mesh_meta_data_ordinal();

  auto devicePartOrdinal1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPartOrdinal1);
  deviceBucketRepo.create_partition(testRank, devicePartOrdinal1);

  HostPartOrdinalView hostPartOrdinal2("", 2);
  hostPartOrdinal2(0) = part1->mesh_meta_data_ordinal();
  hostPartOrdinal2(1) = part2->mesh_meta_data_ordinal();

  auto devicePartOrdinal2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPartOrdinal2);
  deviceBucketRepo.create_partition(testRank, devicePartOrdinal2);

  EXPECT_EQ(4u, deviceBucketRepo.num_partitions(testRank));

  HostPairViewType hostPairs("", 2);
  hostPairs(0) = {0, 2};
  hostPairs(1) = {1, 3};

  auto srcDestPartitions = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPairs);

  deviceBucketRepo.batch_move_entities(entities, srcDestPartitions);

  // empty buckets (with stale data) have not been cleared out from partition 0 and 1 at this point
  auto partition0 = deviceBucketRepo.get_partition(testRank, 0);
  EXPECT_EQ(0u, partition0->m_numEntities);
  EXPECT_EQ(1u, partition0->num_buckets());
  auto partition1 = deviceBucketRepo.get_partition(testRank, 1);
  EXPECT_EQ(0u, partition1->m_numEntities);
  EXPECT_EQ(1u, partition1->num_buckets());
  auto partition2 = deviceBucketRepo.get_partition(testRank, 2);
  EXPECT_EQ(1u, partition2->m_numEntities);
  EXPECT_EQ(1u, partition2->num_buckets());
  EXPECT_EQ(1u, partition2->m_buckets[0]->m_bucketSize);
  auto partition3 = deviceBucketRepo.get_partition(testRank, 2);
  EXPECT_EQ(1u, partition3->m_numEntities);
  EXPECT_EQ(1u, partition3->num_buckets());
  EXPECT_EQ(1u, partition3->m_buckets[0]->m_bucketSize);
}

NGP_TEST_F(NgpBucketRepoBatchOpTest, batch_move_entities_exchange_partitions)
{
  using DevicePartOrdinalView = stk::mesh::PartOrdinalViewType<stk::ngp::MemSpace>;
  using HostPartOrdinalView = typename DevicePartOrdinalView::host_mirror_type;
  using DevicePairViewType = stk::mesh::UnsignedPairViewType<stk::ngp::MemSpace>;
  using HostPairViewType = typename DevicePairViewType::host_mirror_type;
  using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;

  build_empty_mesh();

  create_parts(2);
  auto part0 = m_meta->get_part("part0");
  auto part1 = m_meta->get_part("part1");
  create_node(1, stk::mesh::PartVector{part0});
  create_node(2, stk::mesh::PartVector{part1});
  
  DeviceEntitiesType entities("deviceEntities", 2);
  fill_device_entities(entities);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(testRank));

  HostPairViewType hostPairs("", 2);
  hostPairs(0) = {0, 1};
  hostPairs(1) = {1, 0};

  auto srcDestPartitions = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostPairs);

  deviceBucketRepo.batch_move_entities(entities, srcDestPartitions);

  // empty buckets (with stale data) have not been cleared out from partition 0 and 1 at this point
  auto partition0 = deviceBucketRepo.get_partition(testRank, 0);
  EXPECT_EQ(1u, partition0->m_numEntities);
  EXPECT_EQ(1u, partition0->num_buckets());
  auto partition1 = deviceBucketRepo.get_partition(testRank, 1);
  EXPECT_EQ(1u, partition1->m_numEntities);
  EXPECT_EQ(1u, partition1->num_buckets());
}

#endif
