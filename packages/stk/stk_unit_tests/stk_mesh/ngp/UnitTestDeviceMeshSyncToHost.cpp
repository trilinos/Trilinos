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

#include "stk_mesh/base/Ngp.hpp"

#ifdef STK_USE_DEVICE_MESH

#include <gtest/gtest.h>
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_topology/topology.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/NgpFieldBLAS.hpp"
#include "stk_unit_test_utils/DeviceBucketTestUtils.hpp"
#include <stk_ngp_test/ngp_test.hpp>

namespace {

using Memspace = stk::mesh::NgpMeshDefaultMemSpace;
using DevicePartition = stk::mesh::impl::DevicePartition<Memspace>;
using DeviceBucket = stk::mesh::DeviceBucketT<Memspace>;
using DeviceBucketWrapper = stk::mesh::impl::DeviceBucketPtrWrapper<Memspace>;

class DeviceMeshSyncTester : public ::ngp_testing::Test
{
  public:
    DeviceMeshSyncTester() :
      meta(3),
      bulk(meta,
           stk::parallel_machine_world(),
           stk::mesh::BulkData::AUTO_AURA,
           false,
           std::unique_ptr<stk::mesh::FieldDataManager>(),
           maximumBucketCapacity,
           maximumBucketCapacity),
      block2(&meta.declare_part("block_2", stk::topology::ELEM_RANK)),
      field(meta.declare_field<double>(stk::topology::ELEM_RANK, "my_field"))
    {
      stk::mesh::put_field_on_mesh(field, meta.universal_part(), &field_init_value);
    }

    void setup_mesh(const std::string& meshStr)
    {
      stk::io::fill_mesh(meshStr, bulk);
      block1 = bulk.mesh_meta_data().get_part("block_1");
    }


    void move_entity_on_host()
    {
      stk::mesh::Bucket* bucket = bulk.buckets(stk::topology::ELEM_RANK)[0];
      stk::mesh::Entity elem2 = (*bucket)[1];
      std::vector<stk::mesh::Entity> entities    = {elem2};
      std::vector<stk::mesh::Part*> add_parts    = {block2};
      std::vector<stk::mesh::Part*> remove_parts = {block1};
      bulk.batch_change_entity_parts(entities, add_parts, remove_parts);
    }

    void move_entity_on_device(const std::vector<unsigned>& elemOrdinals)
    {
      move_entity_on_device(0, elemOrdinals);
    }

    void move_entity_on_device(unsigned bucketId, const std::vector<unsigned>& elemOrdinals)
    {
      stk::mesh::Bucket* bucket = bulk.buckets(stk::topology::ELEM_RANK)[0];
      Kokkos::View<stk::mesh::Entity*, Memspace> entities("entities", elemOrdinals.size());
      Kokkos::View<stk::mesh::PartOrdinal*, Memspace> add_part_ordinals("add_part_ords", 1);
      Kokkos::View<stk::mesh::PartOrdinal*, Memspace> remove_part_ordinals("remove_part_rods", 1);

      auto entities_host             = Kokkos::create_mirror_view(entities);
      auto add_part_ordinals_host    = Kokkos::create_mirror_view(add_part_ordinals);
      auto remove_part_ordinals_host = Kokkos::create_mirror_view(remove_part_ordinals);

      for (unsigned i=0; i < elemOrdinals.size(); ++i)
      {
        entities_host(i) = (*bucket)[elemOrdinals[i]];
      }
      add_part_ordinals_host(0)    = block2->mesh_meta_data_ordinal();
      remove_part_ordinals_host(0) = block1->mesh_meta_data_ordinal();

      Kokkos::deep_copy(entities, entities_host);
      Kokkos::deep_copy(add_part_ordinals,    add_part_ordinals_host);
      Kokkos::deep_copy(remove_part_ordinals, remove_part_ordinals_host);

      stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
      deviceMesh.impl_batch_change_entity_parts(entities, add_part_ordinals, remove_part_ordinals);
    }

    void set_field_on_device()
    {
      double offset = field_new_value;
      stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_ngp_mesh(bulk);
      stk::ngp::RangePolicy<stk::ngp::DeviceSpace::exec_space> policy(0, deviceMesh.num_buckets(stk::topology::ELEM_RANK));
      auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
      auto func = KOKKOS_LAMBDA(unsigned bucket_id)
      {
        const DeviceBucket& bucket = deviceMesh.get_bucket(stk::topology::ELEM_RANK, bucket_id);
        for (unsigned i=0; i < bucket.size(); ++i)
        {
          stk::mesh::Entity entity = bucket[i];
          auto entityValues = fieldData.entity_values(entity);
          entityValues(0_comp) = offset + entity.local_offset();
        }
      };

      Kokkos::parallel_for("set_field_on_device", policy, func);
    }

    void set_field_on_host()
    {
      auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
      for (stk::mesh::Bucket* bucket : bulk.buckets(stk::topology::ELEM_RANK))
      {
        for (stk::mesh::Entity entity : *bucket)
        {
          auto entityValues = fieldData.entity_values(entity);
          entityValues(0_comp) = field_new_value + entity.local_offset();
        }
      }
    }

    void test_field()
    {
      auto field_data = field.data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
      for (stk::mesh::Bucket* bucket : bulk.buckets(stk::topology::ELEM_RANK))
      {
        for (stk::mesh::Entity entity : *bucket)
        {
          auto entity_values = field_data.entity_values(entity);
          EXPECT_EQ(entity_values(0_comp), field_new_value + entity.local_offset());
        }
      }
    }

    void test_field_on_device()
    {
      double offset = field_new_value;
      stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_ngp_mesh(bulk);
      stk::ngp::RangePolicy<stk::ngp::DeviceSpace::exec_space> policy(0, deviceMesh.num_buckets(stk::topology::ELEM_RANK));
      auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
      auto func = KOKKOS_LAMBDA(unsigned bucket_id)
      {
        const DeviceBucket& bucket = deviceMesh.get_bucket(stk::topology::ELEM_RANK, bucket_id);
        for (unsigned i=0; i < bucket.size(); ++i)
        {
          stk::mesh::Entity entity = bucket[i];
          auto entityValues = fieldData.entity_values(entity);
          NGP_EXPECT_EQ(entityValues(0_comp), offset + entity.local_offset());
        }
      };

      Kokkos::parallel_for("set_field_on_device", policy, func);
    }

    stk::mesh::MetaData meta;
    stk::unit_test_util::BulkDataTester bulk;
    stk::mesh::Part* block1;
    stk::mesh::Part* block2;
    stk::mesh::Field<double>& field;

    static constexpr unsigned maximumBucketCapacity = 2;
    static constexpr double field_init_value = 1;
    static constexpr double field_new_value = 2;
};


void test_mesh_indices(const stk::mesh::BulkData& bulk)
{
  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank != stk::topology::END_RANK; ++rank)
  {
    for (stk::mesh::Bucket* bucket : bulk.buckets(rank))
    {
      for (size_t i=0; i < bucket->size(); ++i)
      {
        stk::mesh::Entity entity = (*bucket)[i];
        stk::mesh::MeshIndex meshIndex = bulk.mesh_index(entity);
        EXPECT_EQ(meshIndex.bucket, bucket);
        EXPECT_EQ(meshIndex.bucket_ordinal, i);
        EXPECT_EQ(bulk.bucket_ptr(entity), bucket);

        stk::mesh::EntityKey key = bulk.entity_key(entity);
        EXPECT_TRUE(stk::mesh::EntityKey::is_valid_id(key.id()));
        EXPECT_EQ(key.rank(), bucket->entity_rank());
        EXPECT_EQ(bulk.get_entity(key), entity);
      }
    }
  }
}

void test_selector(const stk::mesh::BulkData& bulk)
{
  stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part("block_1");
  stk::mesh::Part& block2 = *bulk.mesh_meta_data().get_part("block_2");
  stk::mesh::Selector common_selector = block1 & block2;
  std::vector<stk::mesh::Entity> common_nodes;
  stk::mesh::get_selected_entities(common_selector, bulk.buckets(stk::topology::NODE_RANK), common_nodes);

  std::set<stk::mesh::Entity> nodes_unique(common_nodes.begin(), common_nodes.end());
  EXPECT_EQ(common_nodes.size(), nodes_unique.size());

  const stk::mesh::FieldBase* coord_field = bulk.mesh_meta_data().coordinate_field();
  auto coord_field_data = coord_field->data<double, stk::mesh::ReadOnly>();

  for (stk::mesh::Entity node : common_nodes)
  {
    auto entity_values = coord_field_data.entity_values(node);
    EXPECT_EQ(entity_values(2_comp), 1.0);
  }
}

void test_all_buckets(const stk::mesh::BulkData& bulk, stk::mesh::NgpMesh deviceMesh)
{
  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    for (size_t i=0; i < bulk.buckets(rank).size(); ++i)
    {
      stk::mesh::Bucket& hostBucket = *bulk.buckets(rank)[i];
      const DeviceBucket& deviceBucket = deviceMesh.get_bucket(rank, i);

      stk::unit_test_util::check_entities(hostBucket, deviceBucket);
      stk::unit_test_util::check_connected_entities(hostBucket, deviceBucket);
    }
  }
}

void test_partition_buckets(stk::unit_test_util::BulkDataTester& bulk, stk::mesh::NgpMesh& deviceMesh)
{
  stk::mesh::impl::BucketRepository& hostBucketRepo = bulk.my_get_bucket_repository();
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  for (stk::mesh::EntityRank rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    for (size_t i=0; i < deviceBucketRepo.num_partitions(rank); ++i)
    {
      stk::mesh::impl::Partition& hostPartition = *(hostBucketRepo.get_partition(rank, i));
      stk::mesh::impl::DevicePartition<Memspace>& devicePartition = *(deviceBucketRepo.get_partition(rank, i));

      stk::unit_test_util::check_partition_attributes(hostPartition, devicePartition);
      stk::unit_test_util::check_partition_ordinals(hostPartition, devicePartition);
      stk::unit_test_util::check_partition_bucket_attributes(hostPartition, devicePartition);

      for (size_t j=0; j < hostPartition.num_buckets(); ++j)
      {
        stk::mesh::Bucket* hostBucket = *(hostPartition.begin() + j);
        DeviceBucketWrapper deviceBucket = devicePartition.get_bucket(j);

        stk::unit_test_util::check_entities(*hostBucket, *deviceBucket);
        stk::unit_test_util::check_connected_entities(*hostBucket, *deviceBucket);
      }
    }
  }
}


}

TEST_F(DeviceMeshSyncTester, CopyToHostMoveFirstElement)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_device({0});
  deviceMesh.update_bulk_data();

  test_mesh_indices(bulk);
  test_selector(bulk);
  test_all_buckets(bulk, deviceMesh);
  test_partition_buckets(bulk, deviceMesh);
}

TEST_F(DeviceMeshSyncTester, CopyToHostMoveSecondElement)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_device({1});
  deviceMesh.update_bulk_data();

  test_mesh_indices(bulk);
  test_selector(bulk);
  test_all_buckets(bulk, deviceMesh);
  test_partition_buckets(bulk, deviceMesh);
}

TEST_F(DeviceMeshSyncTester, CopyToHostMoveFirstElementWithSidesets)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2|sideset:xyzXYZ");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_device({0});
  deviceMesh.update_bulk_data();

  test_mesh_indices(bulk);
  test_selector(bulk);
  test_all_buckets(bulk, deviceMesh);
  test_partition_buckets(bulk, deviceMesh);
}

TEST_F(DeviceMeshSyncTester, CopyToHostLeavePartitionEmpty)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_device({0, 1});
  deviceMesh.update_bulk_data();

  test_mesh_indices(bulk);
  test_selector(bulk);
  test_all_buckets(bulk, deviceMesh);
  test_partition_buckets(bulk, deviceMesh);
}

TEST_F(DeviceMeshSyncTester, CopyToHostLeaveBucketEmpty)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x4");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_device({0, 1});
  deviceMesh.update_bulk_data();

  test_mesh_indices(bulk);
  test_selector(bulk);
  test_all_buckets(bulk, deviceMesh);
  test_partition_buckets(bulk, deviceMesh);
}

TEST_F(DeviceMeshSyncTester, CopyToHostSelectorMemoized)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  test_selector(bulk);
  move_entity_on_device({0});
  deviceMesh.update_bulk_data();
  test_selector(bulk);
}


TEST_F(DeviceMeshSyncTester, CopyToHostField)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  set_field_on_device();
  move_entity_on_device({0});
  std::cout << "\nfirst test on device" << std::endl;
  test_field_on_device();

  deviceMesh.update_bulk_data();
  //std::cout << "second test on device" << std::endl;
  //test_field_on_device();

  std::cout << "testing on host" << std::endl;
  test_field();
}

TEST_F(DeviceMeshSyncTester, CopyToHostFieldModifiedOnHost)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  set_field_on_host();
  move_entity_on_device({0});

  EXPECT_NO_THROW(deviceMesh.update_bulk_data());
  EXPECT_TRUE(field.need_sync_to_host());
  test_field();
}

TEST_F(DeviceMeshSyncTester, CopyToHostFieldErrorModifyOnHostAfterMeshmod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  set_field_on_host();
  move_entity_on_device({0});

  EXPECT_ANY_THROW(set_field_on_host());
}

TEST_F(DeviceMeshSyncTester, CopyToHostFieldModifyOnDeviceAfterMeshmod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  set_field_on_host();
  move_entity_on_device({0});

  set_field_on_device();
  EXPECT_NO_THROW(deviceMesh.update_bulk_data());
  test_field();
}

TEST_F(DeviceMeshSyncTester, CopyToHostFieldModifiedOnHostMoveIntoExistingBucket)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  set_field_on_host();
  EXPECT_TRUE(field.need_sync_to_device());
  EXPECT_FALSE(field.need_sync_to_host());
  move_entity_on_device({0});
  EXPECT_NO_THROW(deviceMesh.update_bulk_data());
  EXPECT_TRUE(field.need_sync_to_host());

  // move the remaining element in block_1 to block_2 (which already has a half full bucket)
  set_field_on_host();
  move_entity_on_device(1, {0});
  EXPECT_NO_THROW(deviceMesh.update_bulk_data());
  EXPECT_TRUE(field.need_sync_to_host());
  EXPECT_FALSE(field.need_sync_to_device());
  test_field();
}

TEST_F(DeviceMeshSyncTester, CopyToHostModifiedOnHost)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_host();
  EXPECT_ANY_THROW(deviceMesh.update_bulk_data());
}

TEST_F(DeviceMeshSyncTester, CopyToHostDuringModCycle)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  bulk.modification_begin();
  EXPECT_ANY_THROW(deviceMesh.update_bulk_data());
}

TEST_F(DeviceMeshSyncTester, CopyToHostModifiedOnHostSyncToDevice)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");
  stk::mesh::NgpMesh* deviceMesh = &stk::mesh::get_updated_ngp_mesh(bulk);

  move_entity_on_host();
  deviceMesh = &stk::mesh::get_updated_ngp_mesh(bulk);
  EXPECT_NO_THROW(deviceMesh->update_bulk_data());
  EXPECT_EQ(deviceMesh->synchronized_count(), bulk.synchronized_count());
}

#endif