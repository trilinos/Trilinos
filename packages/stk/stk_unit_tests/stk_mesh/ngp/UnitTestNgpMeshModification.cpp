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
#include "Kokkos_Core.hpp"
#include "ngp/NgpUnitTestUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SkinMesh.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"

#ifdef STK_USE_DEVICE_MESH
namespace
{
using ngp_unit_test_utils::check_bucket_layout;
using ngp_unit_test_utils::check_entity_parts_on_device;

class NgpBatchChangeEntityParts : public ::ngp_testing::Test
{
public:
  NgpBatchChangeEntityParts()
  {
  }

  void build_empty_mesh(unsigned initialBucketCapacity, unsigned maximumBucketCapacity)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
  }

protected:
  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData * m_meta;
};

using Memspace = stk::mesh::NgpMeshDefaultMemSpace;
using DeviceBucket = stk::mesh::DeviceBucketT<Memspace>;

class NgpBatchChangeEntityPartsTwoBlocks : public ::ngp_testing::Test
{
  public:
    NgpBatchChangeEntityPartsTwoBlocks() :
      meta(3),
      bulk(meta,
           stk::parallel_machine_world(),
           stk::mesh::BulkData::AUTO_AURA,
           false,
           std::unique_ptr<stk::mesh::FieldDataManager>(),
           maximumBucketCapacity,
           maximumBucketCapacity),
      block2(&meta.declare_part("block_2", stk::topology::ELEM_RANK))
    {
      if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
      {
        return;
      }

      stk::io::fill_mesh("generated:1x1x2", bulk);
      block1 = bulk.mesh_meta_data().get_part("block_1");
    }

    void move_entity_on_device(const std::vector<unsigned>& elemOrdinals)
    {
      stk::mesh::Bucket* bucket = bulk.buckets(stk::topology::ELEM_RANK)[0];
      Kokkos::View<stk::mesh::Entity*, Memspace> entities("entities", elemOrdinals.size());
      Kokkos::View<stk::mesh::PartOrdinal*, Memspace> add_part_ordinals("add_part_ords", 1);
      Kokkos::View<stk::mesh::PartOrdinal*, Memspace> remove_part_ordinals("remove_part_ords", 1);

      auto entities_host             = Kokkos::create_mirror_view(entities);
      auto add_part_ordinals_host    = Kokkos::create_mirror_view(add_part_ordinals);
      auto remove_part_ordinals_host = Kokkos::create_mirror_view(remove_part_ordinals);

      for (unsigned i=0; i < elemOrdinals.size(); ++i)
      {
        entities_host(i)             = (*bucket)[elemOrdinals[i]];
      }
      add_part_ordinals_host(0)    = block2->mesh_meta_data_ordinal();
      remove_part_ordinals_host(0) = block1->mesh_meta_data_ordinal();

      Kokkos::deep_copy(entities, entities_host);
      Kokkos::deep_copy(add_part_ordinals,    add_part_ordinals_host);
      Kokkos::deep_copy(remove_part_ordinals, remove_part_ordinals_host);

      stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
      deviceMesh.impl_batch_change_entity_parts(entities, add_part_ordinals, remove_part_ordinals);
    }

    stk::mesh::MetaData meta;
    stk::unit_test_util::BulkDataTester bulk;
    stk::mesh::Part* block1;
    stk::mesh::Part* block2;

    static constexpr unsigned maximumBucketCapacity = 2;
};



stk::mesh::Entity create_node(stk::mesh::BulkData& bulk, stk::mesh::EntityId nodeId,
                              const stk::mesh::PartVector& initialParts = stk::mesh::PartVector())
{
  bulk.modification_begin();
  stk::mesh::Entity newNode = bulk.declare_node(nodeId, initialParts);
  bulk.modification_end();

  return newNode;
}

template <typename MeshType>
void confirm_host_mesh_is_not_synchronized_from_device(const MeshType& ngpMesh)
{
  if constexpr (std::is_same_v<MeshType, stk::mesh::DeviceMesh>) {
    EXPECT_TRUE(ngpMesh.need_update_bulk_data());
  }
  else {
    EXPECT_FALSE(ngpMesh.need_update_bulk_data());  // If host build, HostMesh can't ever be stale
  }
}

template <typename MeshType>
void confirm_host_mesh_is_synchronized_from_device(const MeshType& ngpMesh)
{
  EXPECT_FALSE(ngpMesh.need_update_bulk_data());
}

using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
using DevicePartOrdinalsType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;

using HostEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostExecSpace>;
using HostPartOrdinalsType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::HostExecSpace>;

void fill_device_views_add_remove_part_from_node(DeviceEntitiesType& entities, DevicePartOrdinalsType& addPartOrdinals,
                                                 DevicePartOrdinalsType& removePartOrdinals, stk::mesh::NgpMesh& ngpMesh,
                                                 stk::mesh::Entity node, stk::mesh::Part* addPart,
                                                 stk::mesh::Part* removePart)
{
  const stk::mesh::BulkData& bulk = ngpMesh.get_bulk_on_host();
  stk::mesh::EntityId nodeId = bulk.identifier(node);
  const stk::mesh::PartOrdinal addPartOrdinal = (addPart) ? addPart->mesh_meta_data_ordinal()
                                                          : stk::mesh::InvalidPartOrdinal;
  const stk::mesh::PartOrdinal removePartOrdinal = (removePart) ? removePart->mesh_meta_data_ordinal()
                                                                : stk::mesh::InvalidPartOrdinal;

  Kokkos::parallel_for("Fill Device Views for Part Addition", stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(size_t /*index*/) {
      STK_NGP_ThrowRequireMsg(ngpMesh.identifier(node) == nodeId, "Unexpected node found on device");
      entities(0) = node;

      if (addPartOrdinal != stk::mesh::InvalidPartOrdinal) {
        addPartOrdinals(0) = addPartOrdinal;
      }

      if (removePartOrdinal != stk::mesh::InvalidPartOrdinal) {
        removePartOrdinals(0) = removePartOrdinal;
      }
    });
  Kokkos::fence();
}

DevicePartOrdinalsType create_device_part_ordinal(stk::mesh::PartVector const& vector)
{
  HostPartOrdinalsType hostPartOrdinals("", vector.size());
  for (unsigned i = 0; i < vector.size(); ++i) {
    hostPartOrdinals(i) = vector[i]->mesh_meta_data_ordinal();
  }
  DevicePartOrdinalsType devicePartOrdinals = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace{}, hostPartOrdinals);
  Kokkos::fence();
  return devicePartOrdinals;
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_part_ordinal_match(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank, EntitiesViewType const& entities, DevicePartOrdinalsType const& expectedPartOrdinals)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);

      auto bucket = deviceBucketRepo.get_bucket(rank, fastMeshIndex.bucket_id);
      auto& bucketPartOrdinals = bucket->get_part_ordinals();

      for (int i = expectedPartOrdinals.extent(0)-1, j = bucketPartOrdinals.extent(0)-1; i >= 0; --i, --j) {
        NGP_EXPECT_EQ(bucketPartOrdinals(j), expectedPartOrdinals(i));
      }

      auto partition = deviceBucketRepo.get_partition(rank, bucket->partition_id());
      auto partitionPartOrdinals = partition->superset_part_ordinals();

      for (int i = expectedPartOrdinals.extent(0)-1, j = partitionPartOrdinals.extent(0)-1; i >= 0; --i, --j) {
        NGP_EXPECT_EQ(partitionPartOrdinals(j), expectedPartOrdinals(i));
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_has_parts(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank, EntitiesViewType const& entities, DevicePartOrdinalsType const& expectedPartOrdinals)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);

      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIndex.bucket_id);

      for (unsigned i = 0; i<expectedPartOrdinals.extent(0); ++i) {
        NGP_EXPECT_TRUE(bucket.member(expectedPartOrdinals(i)));
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType>
void check_device_mesh_indices(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(deviceBucketRepo.num_buckets(rank),
    KOKKOS_LAMBDA(const int idx) {
      auto& buckets = deviceBucketRepo.m_buckets[rank];
      auto& bucket = buckets[idx];

      if (!bucket.is_active()) { return; }

      for (unsigned i = 0; i < bucket.size(); ++i) {
        auto entity = bucket[i];

        if (!entity.is_local_offset_valid()) { continue; }

        auto fastMeshIndex = ngpMesh.fast_mesh_index(entity);
        NGP_EXPECT_EQ(fastMeshIndex.bucket_id, bucket.bucket_id());
        NGP_EXPECT_EQ(fastMeshIndex.bucket_ord, i);
      }
    }
  );
  Kokkos::fence();
}

template<typename EntitiesHostViewType>
void init_host_field_data(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& field, EntitiesHostViewType entities)
{
  auto fieldData = field.data<stk::mesh::ReadWrite,stk::ngp::HostSpace>();
  for(unsigned idx=0; idx<entities.extent(0); ++idx) {
    auto entity = entities(idx);
    stk::mesh::EntityId id = mesh.identifier(entity);
    auto fieldEntityValues = fieldData.entity_values(entity);
    for(stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
      fieldEntityValues(i) = id;
    }
  }
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_field_data_is_id(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank,
                                    stk::mesh::Field<double>& field, EntitiesViewType const& entities)
{
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      stk::mesh::EntityId id = ngpMesh.identifier(entity);
      auto fieldEntityValues = fieldData.entity_values(fastMeshIndex);

      for (stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
        stk::mesh::EntityId fieldValue = static_cast<stk::mesh::EntityId>(fieldEntityValues(i));
        NGP_EXPECT_EQ(id, fieldValue);
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_field_data(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank,
                                    stk::mesh::Field<double>& field, EntitiesViewType const& entities,
                                    double expectedValue)
{
  constexpr double tol = 1.e-12;
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      auto fieldEntityValues = fieldData.entity_values(fastMeshIndex);

      for (stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
        double fieldValue = fieldEntityValues(i);
        NGP_EXPECT_NEAR(expectedValue, fieldValue, tol);
      }
    }
  );
  Kokkos::fence();
}

NGP_TEST_F(NgpBatchChangeEntityParts, copy_entity_field_bytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Field<double>& field1 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field1");
  stk::mesh::Field<double>& field2 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field2");
  stk::mesh::Field<double>& field3 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field3");

  const unsigned numFieldComponents = 1;
  const double initVal_field2 = 3.14;
  stk::mesh::put_field_on_mesh(field1, part1, numFieldComponents, nullptr);
  stk::mesh::put_field_on_mesh(field2, part2, numFieldComponents, &initVal_field2);
  stk::mesh::put_field_on_mesh(field3, m_meta->universal_part(), numFieldComponents, nullptr);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = node1;
  hostNodes(1) = node2;
  HostEntitiesType justNode1("justNode1", 1);
  justNode1(0) = node1;
  HostEntitiesType justNode2("justNode2", 1);
  justNode2(0) = node2;

  init_host_field_data(*m_bulk, field3, hostNodes);
  init_host_field_data(*m_bulk, field1, justNode1);
  init_host_field_data(*m_bulk, field2, justNode2);

  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);
  auto devNode1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, justNode1);
  auto devNode2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, justNode2);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field3, devNodes);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field1, devNode1);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field2, devNode2);

//  std::vector<stk::mesh::FieldBase*> field1vec = {&field1};
//  std::vector<stk::mesh::FieldBase*> field2vec = {&field2};
  std::vector<stk::mesh::FieldBase*> fieldsvec = {&field1, &field2, &field3};

  const stk::mesh::MeshIndex& node1HostIndex = m_bulk->mesh_index(node1);
  const stk::mesh::MeshIndex& node2HostIndex = m_bulk->mesh_index(node2);
  stk::mesh::FastMeshIndex node1Index{node1HostIndex.bucket->bucket_id(), node1HostIndex.bucket_ordinal};
  stk::mesh::FastMeshIndex node2Index{node2HostIndex.bucket->bucket_id(), node2HostIndex.bucket_ordinal};

  //The above calls to init_host_field_data marked field-data modified-on-host, and
  //the above calls to check_device_entity_field_data_is_id synchronized field-data to device.
  //Now we're about to modify field2 and field3 field-data on host again, but copy_entity_field_bytes
  //doesn't mark modified. so we need to do that manually for the fields that we will modify.

  field2.synchronize<stk::mesh::ReadWrite>();
  field3.synchronize<stk::mesh::ReadWrite>();

  stk::mesh::copy_entity_field_bytes<stk::ngp::HostSpace>(fieldsvec, node1Index, node2Index);

  double expected = 1;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field3, devNode2, expected);
  expected = initVal_field2;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field2, devNode2, expected);

  init_host_field_data(*m_bulk, field3, hostNodes);
  init_host_field_data(*m_bulk, field1, justNode1);
  init_host_field_data(*m_bulk, field2, justNode2);

  //The above calls to init_host_field_data did synchronize<ReadWrite> on host again, so now we
  //need to synchronize to device in preparation for the copy_entity_field_bytes function to
  //operate on up-to-date data on device.

  field1.synchronize<stk::mesh::ReadOnly,stk::ngp::DeviceSpace>();
  field2.synchronize<stk::mesh::ReadWrite,stk::ngp::DeviceSpace>();
  field3.synchronize<stk::mesh::ReadWrite,stk::ngp::DeviceSpace>();

  stk::mesh::copy_entity_field_bytes<stk::ngp::DeviceSpace>(fieldsvec, node1Index, node2Index);

  expected = 1;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field3, devNode2, expected);
  expected = initVal_field2;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field2, devNode2, expected);
}

void call_copy_entity_field_bytes_on_device(const std::vector<stk::mesh::FieldBase*>& fieldsVec,
                                            stk::mesh::FastMeshIndex node1Index,
                                            stk::mesh::FastMeshIndex node2Index)
{
  for(stk::mesh::FieldBase* field : fieldsVec) {
    field->synchronize<stk::mesh::ReadWrite,stk::ngp::DeviceSpace>();
  }

  using FieldDataBytesType = stk::mesh::FieldDataBytes<stk::ngp::DeviceSpace>;
  using FieldDataBytesViewType = Kokkos::View<FieldDataBytesType*, stk::ngp::DeviceSpace::mem_space>;

  FieldDataBytesViewType deviceFieldDataBytes[stk::topology::NUM_RANKS] = {
    FieldDataBytesViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,"FieldDataBytes"),0),
    FieldDataBytesViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,"FieldDataBytes"),0),
    FieldDataBytesViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,"FieldDataBytes"),0),
    FieldDataBytesViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,"FieldDataBytes"),0),
    FieldDataBytesViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,"FieldDataBytes"),0)
  };
  stk::mesh::assemble_field_data_bytes_on_device<stk::ngp::DeviceSpace>(fieldsVec, deviceFieldDataBytes[stk::topology::NODE_RANK]);

  using InitValsPtrViewType = Kokkos::View<stk::mesh::BytePtr*, stk::ngp::HostPinnedSpace>;
  InitValsPtrViewType initValsPtrsView("init-vals",fieldsVec.size());
  stk::mesh::assemble_field_init_vals_on_device(fieldsVec, initValsPtrsView);

  auto notParallel = Kokkos::RangePolicy<stk::ngp::DeviceSpace::exec_space>(0, 1);
  Kokkos::parallel_for("call_copy_entity_bytes_on_device", notParallel,
    KOKKOS_LAMBDA(const int) {
      stk::mesh::copy_entity_bytes_kernel<stk::mesh::Layout::Left>(deviceFieldDataBytes[stk::topology::NODE_RANK], node1Index, node2Index, initValsPtrsView);
    }
  );
  Kokkos::fence();
}

NGP_TEST_F(NgpBatchChangeEntityParts, copy_entity_field_bytes_call_on_device)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Field<double>& field1 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field1");
  stk::mesh::Field<double>& field2 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field2");
  stk::mesh::Field<double>& field3 = m_meta->declare_field<double>(stk::topology::NODE_RANK, "field3");

  const unsigned numFieldComponents = 1;
  const double initVal_field2 = 3.14;
  stk::mesh::put_field_on_mesh(field1, part1, numFieldComponents, nullptr);
  stk::mesh::put_field_on_mesh(field2, part2, numFieldComponents, &initVal_field2);
  stk::mesh::put_field_on_mesh(field3, m_meta->universal_part(), numFieldComponents, nullptr);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = node1;
  hostNodes(1) = node2;
  HostEntitiesType justNode1("justNode1", 1);
  justNode1(0) = node1;
  HostEntitiesType justNode2("justNode2", 1);
  justNode2(0) = node2;

  init_host_field_data(*m_bulk, field3, hostNodes);
  init_host_field_data(*m_bulk, field1, justNode1);
  init_host_field_data(*m_bulk, field2, justNode2);

  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);
  auto devNode1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, justNode1);
  auto devNode2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, justNode2);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field3, devNodes);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field1, devNode1);
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, field2, devNode2);

  std::vector<stk::mesh::FieldBase*> fieldsvec = {&field1, &field2, &field3};

  const stk::mesh::MeshIndex& node1HostIndex = m_bulk->mesh_index(node1);
  const stk::mesh::MeshIndex& node2HostIndex = m_bulk->mesh_index(node2);
  stk::mesh::FastMeshIndex node1Index{node1HostIndex.bucket->bucket_id(), node1HostIndex.bucket_ordinal};
  stk::mesh::FastMeshIndex node2Index{node2HostIndex.bucket->bucket_id(), node2HostIndex.bucket_ordinal};

  //The above calls to init_host_field_data marked field-data modified-on-host, and
  //the above calls to check_device_entity_field_data_is_id synchronized field-data to device.
  //Now we're about to modify field2 and field3 field-data on host again, but copy_entity_field_bytes
  //doesn't mark modified. so we need to do that manually for the fields that we will modify.

  field2.synchronize<stk::mesh::ReadWrite>();
  field3.synchronize<stk::mesh::ReadWrite>();

  stk::mesh::copy_entity_field_bytes<stk::ngp::HostSpace>(fieldsvec, node1Index, node2Index);

  double expected = 1;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field3, devNode2, expected);
  expected = initVal_field2;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field2, devNode2, expected);

  init_host_field_data(*m_bulk, field3, hostNodes);
  init_host_field_data(*m_bulk, field1, justNode1);
  init_host_field_data(*m_bulk, field2, justNode2);

  //The above calls to init_host_field_data did synchronize<ReadWrite> on host again, so now we
  //need to synchronize to device in preparation for the copy_entity_field_bytes function to
  //operate on up-to-date data on device.

  call_copy_entity_field_bytes_on_device(fieldsvec, node1Index, node2Index);

  expected = 1;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field3, devNode2, expected);
  expected = initVal_field2;
  check_device_entity_field_data(ngpMesh, stk::topology::NODE_RANK, field2, devNode2, expected);
}

NGP_TEST_F(NgpBatchChangeEntityParts, addPartToNode_host)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  std::vector<stk::mesh::Entity> entities {node1};
  std::vector<stk::mesh::Part*> addParts {&part2};
  std::vector<stk::mesh::Part*> removeParts;

  m_bulk->batch_change_entity_parts(entities, addParts, removeParts);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, addPartToNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  HostEntitiesType entities("hostEntities", 1);
  HostPartOrdinalsType addPartOrdinals("hostAddParts", 1);
  HostPartOrdinalsType removePartOrdinals("hostRemoveParts", 0);

  entities(0) = node1;
  addPartOrdinals(0) = part2.mesh_meta_data_ordinal();

  stk::mesh::HostMesh hostMesh(*m_bulk);
  hostMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  hostMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, addPartToNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, nullptr);

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  check_entity_parts_on_device(ngpMesh, entities, addPartOrdinals, removePartOrdinals, stk::topology::NODE_RANK);

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}


NGP_TEST_F(NgpBatchChangeEntityParts, removePartFromNode_host)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1, &part2});
  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);

  std::vector<stk::mesh::Entity> entities {node1};
  std::vector<stk::mesh::Part*> addParts;
  std::vector<stk::mesh::Part*> removeParts {&part1};

  m_bulk->batch_change_entity_parts(entities, addParts, removeParts);
  stk::mesh::get_updated_ngp_mesh(*m_bulk);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, removePartFromNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1, &part2});
  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);

  HostEntitiesType entities("hostEntities", 1);
  HostPartOrdinalsType addPartOrdinals("hostAddParts", 0);
  HostPartOrdinalsType removePartOrdinals("hostRemoveParts", 1);

  entities(0) = node1;
  removePartOrdinals(0) = part1.mesh_meta_data_ordinal();

  stk::mesh::HostMesh hostMesh(*m_bulk);
  hostMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  hostMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, removePartFromNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1, &part2});
  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, nullptr, &part1);

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}


NGP_TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_host)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  std::vector<stk::mesh::Entity> entities {node1};
  std::vector<stk::mesh::Part*> addParts {&part2};
  std::vector<stk::mesh::Part*> removeParts {&part1};

  m_bulk->batch_change_entity_parts(entities, addParts, removeParts);
  stk::mesh::get_updated_ngp_mesh(*m_bulk);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  HostEntitiesType entities("hostEntities", 1);
  HostPartOrdinalsType addPartOrdinals("hostAddParts", 1);
  HostPartOrdinalsType removePartOrdinals("hostRemoveParts", 1);

  entities(0) = node1;
  addPartOrdinals(0) = part2.mesh_meta_data_ordinal();
  removePartOrdinals(0) = part1.mesh_meta_data_ordinal();

  stk::mesh::HostMesh hostMesh(*m_bulk);
  hostMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  hostMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, &part1);

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

#ifdef USE_IMPL_DEVICE_MESH_MOD
  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part2});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, entities, expectedDevicePartOrdinal);
#endif

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, multipleDeviceMeshMods)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, nullptr);

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part3, nullptr);

  Kokkos::fence();

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

#ifdef USE_IMPL_DEVICE_MESH_MOD
  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1, &part2, &part3});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, entities, expectedDevicePartOrdinal);
#endif

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, failedHostAccessAfterDeviceMeshMod)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, nullptr);

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  if constexpr (std::is_same_v<stk::mesh::NgpMesh, stk::mesh::DeviceMesh>) {
    EXPECT_ANY_THROW(m_bulk->buckets(stk::topology::NODE_RANK));
    EXPECT_ANY_THROW(m_bulk->get_buckets(stk::topology::NODE_RANK, m_meta->universal_part()));
    EXPECT_ANY_THROW(m_bulk->modification_begin());
    EXPECT_ANY_THROW(m_bulk->batch_change_entity_parts(stk::mesh::EntityVector{node1}, stk::mesh::PartVector{},
                                                       stk::mesh::PartVector{}));
    EXPECT_ANY_THROW(stk::mesh::skin_mesh(*m_bulk, stk::mesh::PartVector{&part1}));
  }
}

NGP_TEST_F(NgpBatchChangeEntityParts, addPartToNode_impl_set_new_part_lists)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, nullptr);

  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh,entities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();
  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*maxNewNumPartsPerEntity);

  using PartOrdinalsIndicesViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsIndicesViewType partOrdinalsIndices("partOrdinalsBeginEndPerEntity", entities.size());

  stk::mesh::impl::set_new_part_list_per_entity(ngpMesh, entities, addPartOrdinals, removePartOrdinals,
                                                maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsIndices);

  const stk::mesh::PartVector& newParts = m_bulk->bucket(node1).supersets();

  auto newPartOrdinalsHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, newPartOrdinalsPerEntity);

  unsigned expectedSize = 1+newParts.size();
  EXPECT_EQ(expectedSize, newPartOrdinalsHost.size());
  for(unsigned i=0; i<newParts.size(); ++i) {
    EXPECT_EQ(newParts[i]->mesh_meta_data_ordinal(), newPartOrdinalsHost(i));
  }

  auto hostPartOrdinalProxy = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, partOrdinalsIndices);
  for(unsigned i = 0; i < hostPartOrdinalProxy.extent(0); ++i) {
    EXPECT_EQ(stk::topology::NODE_RANK, hostPartOrdinalProxy(i).rank);
    EXPECT_EQ(newPartOrdinalsPerEntity.data(), hostPartOrdinalProxy(i).startPtr);
    EXPECT_EQ(expectedSize, hostPartOrdinalProxy(i).length);
  }
}

NGP_TEST_F(NgpBatchChangeEntityParts, removePartToNode_impl_set_new_part_lists)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1, &part2, &part3});
  check_bucket_layout(*m_bulk, {{{"part1", "part2", "part3"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, nullptr, &part3);

  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh,entities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();
  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*maxNewNumPartsPerEntity);

  using PartOrdinalsIndicesViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsIndicesViewType partOrdinalsIndices("partOrdinalsBeginEndPerEntity", entities.size());

  stk::mesh::impl::set_new_part_list_per_entity(ngpMesh, entities, addPartOrdinals, removePartOrdinals,
                                                maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsIndices);

  const stk::mesh::PartVector& newParts = m_bulk->bucket(node1).supersets();

  auto newPartOrdinalsHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, newPartOrdinalsPerEntity);

  unsigned expectedSize = newParts.size()-1;
  EXPECT_EQ(expectedSize, newPartOrdinalsHost.size()-1);
  EXPECT_EQ(maxNewNumPartsPerEntity, newPartOrdinalsHost.size());
  for(unsigned i=0; i < expectedSize; ++i) {
    EXPECT_EQ(newParts[i]->mesh_meta_data_ordinal(), newPartOrdinalsHost(i));
  }

  auto hostPartOrdinalProxy = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, partOrdinalsIndices);
  for(unsigned i = 0; i < hostPartOrdinalProxy.extent(0); ++i) {
    EXPECT_EQ(stk::topology::NODE_RANK, hostPartOrdinalProxy(i).rank);
    EXPECT_EQ(newPartOrdinalsPerEntity.data(), hostPartOrdinalProxy(i).startPtr);
    EXPECT_EQ(expectedSize, hostPartOrdinalProxy(i).length);
  }
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_addPartToNode)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1});
  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, &part2, nullptr);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals));

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1, &part2});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, entities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_removePartFromNode)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  const unsigned nodeId = 1;
  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId, {&part1, &part2});
  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);

  DeviceEntitiesType entities("deviceEntities", 1);
  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  fill_device_views_add_remove_part_from_node(entities, addPartOrdinals, removePartOrdinals, ngpMesh,
                                              node1, nullptr, &part1);
  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals));

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part2});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, entities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_remove_parts_create_multiple_buckets_in_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1, &part2});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1, &part3});

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DeviceEntitiesType devEntities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  hostEntities(1) = node2;
  Kokkos::deep_copy(devEntities, hostEntities);
  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 2);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = part2.mesh_meta_data_ordinal();
  hostRemoveParts(1) = part3.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts));

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_move_nodes_to_another_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1, &part2});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1, &part2});

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DeviceEntitiesType devEntities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  hostEntities(1) = node2;
  Kokkos::deep_copy(devEntities, hostEntities);
  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = part2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts));

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_add_parts_create_multiple_buckets_in_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1, &part2});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1, &part3});

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DeviceEntitiesType devEntities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  hostEntities(1) = node2;
  Kokkos::deep_copy(devEntities, hostEntities);
  DevicePartOrdinalsType devAddParts("", 2);
  DevicePartOrdinalsType devRemoveParts("", 0);
  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = part2.mesh_meta_data_ordinal();
  hostAddParts(1) = part3.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts));

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1, &part2, &part3});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_remove_parts_create_single_bucket_in_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(24, 24);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1, &part2});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1, &part3});

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DeviceEntitiesType devEntities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  hostEntities(1) = node2;
  Kokkos::deep_copy(devEntities, hostEntities);
  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 2);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = part2.mesh_meta_data_ordinal();
  hostRemoveParts(1) = part3.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts));

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchChangeEntityParts, impl_batch_change_entity_parts_add_parts_create_single_buckets_in_partition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(24, 24);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  stk::mesh::Part & part3 = m_meta->declare_part_with_topology("part3", stk::topology::NODE);

  const unsigned nodeId1 = 1;
  const unsigned nodeId2 = 2;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1, &part2});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1, &part3});

  DeviceEntitiesType devEntities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  hostEntities(1) = node2;
  Kokkos::deep_copy(devEntities, hostEntities);
  DevicePartOrdinalsType devAddParts("", 2);
  DevicePartOrdinalsType devRemoveParts("", 0);
  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = part2.mesh_meta_data_ordinal();
  hostAddParts(1) = part3.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1});

  check_device_entity_has_parts(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(2u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts));

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&part1, &part2, &part3});
  check_device_entity_has_parts(ngpMesh, stk::topology::NODE_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_mesh_indices(ngpMesh, stk::topology::NODE_RANK);
}

using NgpBatchChangeEntityPartsDeathTest = NgpBatchChangeEntityParts;

TEST_F(NgpBatchChangeEntityPartsDeathTest, DISABLED_check_impl_batch_change_entity_parts_detect_internal_part)
{
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  auto& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  const unsigned nodeId1 = 1;

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  [[maybe_unused]] auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = node1;
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);
  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = m_meta->get_part(1).mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

#ifndef NDEBUG
  ASSERT_DEATH(ngpMesh.impl_batch_change_entity_parts(devEntities, devAddParts, devRemoveParts), "");
  ASSERT_DEATH(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts), "");
#endif
}

template<typename ViewType>
void init_sorted(ViewType& view)
{
  auto devicePolicy = stk::ngp::DeviceRangePolicy(0,1);
  Kokkos::parallel_for("init", devicePolicy, KOKKOS_LAMBDA(const int& /*idx*/) {
    view(0) = 0;
    view(1) = 2;
    view(2) = 2;
    view(3) = 4;
    view(4) = 9;
  });
}

template<typename ViewType>
void init_unsorted(ViewType& view)
{
  auto devicePolicy = stk::ngp::DeviceRangePolicy(0,1);
  Kokkos::parallel_for("init", devicePolicy, KOKKOS_LAMBDA(const int& /*idx*/) {
    view(0) = 0;
    view(1) = 2;
    view(2) = 1;
    view(3) = 4;
    view(4) = 9;
  });
}

TEST(NgpMeshImpl, is_sorted)
{
  Kokkos::View<int*,stk::ngp::ExecSpace> emptyView("emptyView", 0);
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},emptyView));

  Kokkos::View<int*,stk::ngp::ExecSpace> sortedView("sortedView", 5);
  init_sorted(sortedView);

  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},sortedView));

  Kokkos::View<int*,stk::ngp::ExecSpace> unsortedView("unsortedView", 5);
  init_unsorted(unsortedView);

  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},unsortedView));
}

TEST(NgpMeshImpl, get_sorted_view)
{
  Kokkos::View<int*,stk::ngp::ExecSpace> unsortedView("unsortedView", 5);
  init_unsorted(unsortedView);

  Kokkos::View<int*,stk::ngp::ExecSpace> sortedView = stk::mesh::impl::get_sorted_view(unsortedView);
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},sortedView));
}

TEST(NgpMeshConstructionTest, prevent_update_during_mesh_mod)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(3);
  auto bulk = builder.create();

  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));

  bulk->modification_begin();
  EXPECT_ANY_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));
  bulk->modification_end();

  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));
}

class NgpBatchChangeEntityPartsInducedPartMembership : public NgpBatchChangeEntityParts
{
public:
  NgpBatchChangeEntityPartsInducedPartMembership()
  {}

  template <typename PartOrdinalsProxyViewType>
  void check_part_ordinals_induced_from_one_elem(stk::mesh::PartVector const& partsForElem, stk::mesh::PartVector const& partsForNodes,
                                                 PartOrdinalsProxyViewType const& partOrdinalsProxy)
  {
    HostPartOrdinalsType hostElemPartOrdinals("hostElemPartOrdinals", partsForElem.size());
    HostPartOrdinalsType hostNodesPartOrdinals("hostNodesPartOrdinals", partsForNodes.size());

    auto copyFromPartVector = [&](stk::mesh::PartVector const& src, HostPartOrdinalsType& dest) {
                                for (unsigned i = 0; i < src.size(); ++i) {
                                  dest(i) = src[i]->mesh_meta_data_ordinal();
                                }
                              };

    copyFromPartVector(partsForElem, hostElemPartOrdinals);
    copyFromPartVector(partsForNodes, hostNodesPartOrdinals);

    auto deviceElemPartOrdinals = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostElemPartOrdinals);
    auto deviceNodesPartOrdinals = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodesPartOrdinals);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      NGP_EXPECT_EQ(stk::topology::ELEM_RANK, partOrdinalsProxy(0).rank);

      for (int i = deviceElemPartOrdinals.extent(0)-1, j = partOrdinalsProxy(0).length-1; i >= 0; --i, --j) {
        auto expected = deviceElemPartOrdinals(i);
        auto partOrdinalInProxy = partOrdinalsProxy(0).startPtr + j;
        NGP_EXPECT_EQ(expected, *partOrdinalInProxy);
      }

      for (unsigned i = 1; i < partOrdinalsProxy.extent(0); ++i) {
        for (int j = deviceNodesPartOrdinals.extent(0)-1, k = partOrdinalsProxy(i).length-1; j >= 0; --j, --k) {
          auto expected = deviceNodesPartOrdinals(j);
          auto partOrdinalInProxy = partOrdinalsProxy(i).startPtr + k;
          NGP_EXPECT_EQ(expected, *partOrdinalInProxy);
        }
      }
    });
  }

  template <typename DeviceMeshType, typename EntityViewType>
  void check_device_buckets_part_ordinal_does_not_contain_part(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank, EntityViewType const& entities, DevicePartOrdinalsType const& invalidPartOrdinals)
  {
    auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

    Kokkos::parallel_for(entities.extent(0),
      KOKKOS_LAMBDA(const int idx) {
        auto entity = entities(idx);
        auto fastMeshIndex = ngpMesh.device_mesh_index(entity);

        auto bucket = deviceBucketRepo.get_bucket(rank, fastMeshIndex.bucket_id);

        for (unsigned i = 0; i < invalidPartOrdinals.extent(0); ++i) {
          auto partOrdinal = invalidPartOrdinals(i);
          NGP_EXPECT_FALSE(bucket->member(partOrdinal));
        }
      }
    );
    Kokkos::fence();
  }
};

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_populate_all_downward_connected_entities_one_element)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(entities, hostEntities);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  EXPECT_EQ(8, maxNumDownwardConnectedEntities);

  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);
  EXPECT_EQ(9u, maxNumEntitiesForInducingParts);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);

  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  auto hostWrappedEntities = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, wrappedEntities);
  EXPECT_EQ(9u, hostWrappedEntities.extent(0));

  for (unsigned i = 0; i < hostWrappedEntities.extent(0); ++i) {
    auto wrappedEntity = hostWrappedEntities(i);
    auto rank = m_bulk->entity_rank(wrappedEntity);

    if (rank == stk::topology::ELEM_RANK) {
      EXPECT_FALSE(wrappedEntity.isForPartInduction);
    } else {
      EXPECT_TRUE(wrappedEntity.isForPartInduction);
    }
    EXPECT_EQ(i+1, static_cast<stk::mesh::Entity>(wrappedEntity).m_value);
  }
}

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_populate_all_downward_connected_entities_one_element_and_one_node)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 2);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;

  auto node = m_bulk->get_entity(stk::topology::NODE_RANK, 1);
  hostEntities(1) = node;
  Kokkos::deep_copy(entities, hostEntities);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  EXPECT_EQ(8, maxNumDownwardConnectedEntities);

  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);
  EXPECT_EQ(18u, maxNumEntitiesForInducingParts);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);

  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  EXPECT_EQ(9u, wrappedEntities.size());

  auto hostWrappedEntities = Kokkos::create_mirror_view_and_copy(stk::ngp::HostExecSpace{}, wrappedEntities);
  EXPECT_EQ(9u, hostWrappedEntities.extent(0));

  for (unsigned i = 0; i < hostWrappedEntities.extent(0); ++i) {
    auto wrappedEntity = hostWrappedEntities(i);
    auto rank = m_bulk->entity_rank(wrappedEntity);

    if (rank == stk::topology::ELEM_RANK) {
      EXPECT_FALSE(wrappedEntity.isForPartInduction);
    } else {
      EXPECT_TRUE(wrappedEntity.isForPartInduction);
    }
    EXPECT_EQ(i+1, static_cast<stk::mesh::Entity>(wrappedEntity).m_value);
  }
}

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_set_new_part_lists_no_parts_to_induce_add)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& unrankedPart = m_meta->declare_part("Part1");

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(entities, hostEntities);

  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  auto hostAddPartOrdinals = Kokkos::create_mirror_view(addPartOrdinals);
  hostAddPartOrdinals(0) = unrankedPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(addPartOrdinals, hostAddPartOrdinals);

  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);
  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  // determine resulting parts per entity including inducible parts
  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh, wrappedEntities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
  DevicePartOrdinalsType sortedAddPartOrdinals = stk::mesh::impl::get_sorted_view(addPartOrdinals);

  // create part ordinals proxy, sort and unique it (and realloc)
  // A proxy indices view to part ordinals: (rank, startPtr, length)
  using PartOrdinalsProxyViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
  stk::mesh::impl::set_new_part_list_per_entity_with_induced_parts(ngpMesh, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                                   maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  EXPECT_EQ(9u, partOrdinalsProxy.extent(0));

  stk::mesh::PartVector partsForElem{&elemPart1, &unrankedPart};
  stk::mesh::PartVector partsForNodes{&elemPart1};
  check_part_ordinals_induced_from_one_elem(partsForElem, partsForNodes, partOrdinalsProxy);
}

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_set_new_part_lists_no_parts_to_induce_remove)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& unrankedPart = m_meta->declare_part("Part1");

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1, &unrankedPart};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(entities, hostEntities);

  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);
  auto hostRemovePartOrdinals = Kokkos::create_mirror_view(removePartOrdinals);
  hostRemovePartOrdinals(0) = unrankedPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(removePartOrdinals, hostRemovePartOrdinals);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);
  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  // determine resulting parts per entity including inducible parts
  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh, wrappedEntities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
  DevicePartOrdinalsType sortedAddPartOrdinals = stk::mesh::impl::get_sorted_view(addPartOrdinals);

  // create part ordinals proxy, sort and unique it (and realloc)
  // A proxy indices view to part ordinals: (rank, startPtr, length)
  using PartOrdinalsProxyViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
  stk::mesh::impl::set_new_part_list_per_entity_with_induced_parts(ngpMesh, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                                   maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  EXPECT_EQ(9u, partOrdinalsProxy.extent(0));

  stk::mesh::PartVector partsForElem{&elemPart1};
  stk::mesh::PartVector partsForNodes{&elemPart1};
  check_part_ordinals_induced_from_one_elem(partsForElem, partsForNodes, partOrdinalsProxy);
}

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_set_new_part_lists_has_parts_to_induce_add)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(entities, hostEntities);

  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 1);
  auto hostAddPartOrdinals = Kokkos::create_mirror_view(addPartOrdinals);
  hostAddPartOrdinals(0) = elemPart2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(addPartOrdinals, hostAddPartOrdinals);

  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);
  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  // determine resulting parts per entity including inducible parts
  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh, wrappedEntities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
  DevicePartOrdinalsType sortedAddPartOrdinals = stk::mesh::impl::get_sorted_view(addPartOrdinals);

  // create part ordinals proxy, sort and unique it (and realloc)
  // A proxy indices view to part ordinals: (rank, startPtr, length)
  using PartOrdinalsProxyViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
  stk::mesh::impl::set_new_part_list_per_entity_with_induced_parts(ngpMesh, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                                   maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  EXPECT_EQ(9u, partOrdinalsProxy.extent(0));

  stk::mesh::PartVector partsForElem{&elemPart1, &elemPart2};
  stk::mesh::PartVector partsForNodes{&elemPart1, &elemPart2};
  check_part_ordinals_induced_from_one_elem(partsForElem, partsForNodes, partOrdinalsProxy);
}

NGP_TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, check_impl_set_new_part_lists_has_parts_to_induce_remove)
{
  using EntityWrapperViewType = Kokkos::View<stk::mesh::impl::EntityWrapper*>;

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType entities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(entities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(entities, hostEntities);

  DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
  DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 1);
  auto hostRemovePartOrdinals = Kokkos::create_mirror_view(removePartOrdinals);
  hostRemovePartOrdinals(0) = elemPart2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(removePartOrdinals, hostRemovePartOrdinals);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  auto maxNumDownwardConnectedEntities = stk::mesh::impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  stk::mesh::impl::populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);
  stk::mesh::impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, stk::ngp::ExecSpace{});

  // determine resulting parts per entity including inducible parts
  const unsigned maxCurrentNumPartsPerEntity = stk::mesh::impl::get_max_num_parts_per_entity(ngpMesh, wrappedEntities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  DevicePartOrdinalsType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
  DevicePartOrdinalsType sortedAddPartOrdinals = stk::mesh::impl::get_sorted_view(addPartOrdinals);

  // create part ordinals proxy, sort and unique it (and realloc)
  // A proxy indices view to part ordinals: (rank, startPtr, length)
  using PartOrdinalsProxyViewType = Kokkos::View<stk::mesh::impl::PartOrdinalsProxyIndices*>;
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
  stk::mesh::impl::set_new_part_list_per_entity_with_induced_parts(ngpMesh, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                                   maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  EXPECT_EQ(9u, partOrdinalsProxy.extent(0));

  stk::mesh::PartVector partsForElem{&elemPart1};
  stk::mesh::PartVector partsForNodes{&elemPart1};
  check_part_ordinals_induced_from_one_elem(partsForElem, partsForNodes, partOrdinalsProxy);
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_connectivity(DeviceMeshType& ngpMesh, stk::mesh::EntityRank entityRank, const EntitiesViewType& entities,
                            stk::mesh::EntityRank connRank, unsigned expectedNumConnected)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      auto connectedEntities = ngpMesh.get_connected_entities(entityRank, fastMeshIndex, connRank);
      NGP_EXPECT_EQ(expectedNumConnected, connectedEntities.size());
    }
  );
  Kokkos::fence();
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, add_unranked_part_to_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& unrankedPart = m_meta->declare_part("unrankedPart");

  const stk::mesh::EntityId elemId = 1;
  const unsigned nodeId = 1;

  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);

  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = unrankedPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  auto node = m_bulk->get_entity(stk::topology::NODE_RANK, nodeId);

  DeviceEntitiesType devNodeEntities("deviceNodeEntities", 1);
  auto hostNodeEntities = Kokkos::create_mirror_view(devNodeEntities);
  hostNodeEntities(0) = node;
  Kokkos::deep_copy(devNodeEntities, hostNodeEntities);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  const unsigned expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  const unsigned expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart1, &unrankedPart});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodeEntities, expectedDevicePartOrdinal);

  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, remove_unranked_part_from_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& unrankedPart = m_meta->declare_part("unrankedPart");

  const stk::mesh::EntityId elemId = 1;
  const unsigned nodeId = 1;

  stk::mesh::PartVector parts{&elemPart1, &unrankedPart};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);

  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = unrankedPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  auto node = m_bulk->get_entity(stk::topology::NODE_RANK, nodeId);

  DeviceEntitiesType devNodeEntities("deviceNodeEntities", 1);
  auto hostNodeEntities = Kokkos::create_mirror_view(devNodeEntities);
  hostNodeEntities(0) = node;
  Kokkos::deep_copy(devNodeEntities, hostNodeEntities);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  const unsigned expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  const unsigned expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);

  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodeEntities, expectedDevicePartOrdinal);

  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, add_ranked_part_to_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  const unsigned nodeId = 1;

  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);

  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = elemPart2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  auto node = m_bulk->get_entity(stk::topology::NODE_RANK, nodeId);

  DeviceEntitiesType devNodeEntities("deviceNodeEntities", 1);
  auto hostNodeEntities = Kokkos::create_mirror_view(devNodeEntities);
  hostNodeEntities(0) = node;
  Kokkos::deep_copy(devNodeEntities, hostNodeEntities);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  const unsigned expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  const unsigned expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart1, &elemPart2});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodeEntities, expectedDevicePartOrdinal);
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodeEntities,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, remove_ranked_part_from_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;
  const unsigned nodeId = 1;

  stk::mesh::PartVector parts{&elemPart1, &elemPart2};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);

  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = elemPart1.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  auto node = m_bulk->get_entity(stk::topology::NODE_RANK, nodeId);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  DeviceEntitiesType devNodeEntities("deviceNodeEntities", 1);
  auto hostNodeEntities = Kokkos::create_mirror_view(devNodeEntities);
  hostNodeEntities(0) = node;
  Kokkos::deep_copy(devNodeEntities, hostNodeEntities);

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart2});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodeEntities, expectedDevicePartOrdinal);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, element_with_ranked_part_try_remove_ranked_part_from_lower_rank_entity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;

  stk::mesh::PartVector parts{&elemPart};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  DevicePartOrdinalsType devAddParts("", 0);
  HostPartOrdinalsType hostRemoveParts("", 1);
  hostRemoveParts(0) = elemPart.mesh_meta_data_ordinal();
  auto devRemoveParts = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostRemoveParts);

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->get_entity(stk::topology::NODE_RANK, 3);
  hostNodes(1) = m_bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devNodes, devAddParts, devRemoveParts));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, try_add_ranked_part_to_lower_rank_entity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId = 1;

  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem = stk::mesh::declare_element(*m_bulk, parts, elemId, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  m_bulk->modification_end();

  HostPartOrdinalsType hostAddParts("", 1);
  hostAddParts(0) = elemPart2.mesh_meta_data_ordinal();
  auto devAddParts = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostAddParts);
  DevicePartOrdinalsType devRemoveParts("", 0);

  DeviceEntitiesType devEntities("deviceEntities", 1);
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  hostEntities(0) = elem;
  Kokkos::deep_copy(devEntities, hostEntities);

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->get_entity(stk::topology::NODE_RANK, 3);
  hostNodes(1) = m_bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devNodes, devAddParts, devRemoveParts));

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{&elemPart1});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, two_elements_with_diff_ranked_parts_remove_ranked_part_from_one_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  m_meta->set_part_id(part1, 1);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);
  m_meta->set_part_id(part2, 2);

  stk::io::put_io_part_attribute(part1);
  stk::io::put_io_part_attribute(part2);

  const stk::mesh::EntityId elemId1 = 1;
  const stk::mesh::EntityId elemId2 = 2;

  stk::mesh::PartVector elemPart1{&part1};
  stk::mesh::PartVector elemPart2{&part2};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem1 = stk::mesh::declare_element(*m_bulk, elemPart1, elemId1, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  const stk::mesh::Entity elem2 = stk::mesh::declare_element(*m_bulk, elemPart2, elemId2, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  HostEntitiesType hostEntities1("hostEntities1", 1);
  hostEntities1(0) = elem1;
  auto devEntities1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostEntities1);

  DevicePartOrdinalsType expectedDevicePartOrdinal = create_device_part_ordinal(elemPart1);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities1, expectedDevicePartOrdinal);

  HostEntitiesType hostEntities2("hostEntities2", 1);
  hostEntities2(0) = elem2;
  auto devEntities2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostEntities2);

  expectedDevicePartOrdinal = create_device_part_ordinal(elemPart2);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities2, expectedDevicePartOrdinal);

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->begin_nodes(elem2)[0];
  hostNodes(1) = m_bulk->begin_nodes(elem2)[1];
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  expectedDevicePartOrdinal = create_device_part_ordinal(elemPart2);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);

  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = part1.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities1, devAddParts, devRemoveParts));

  expectedDevicePartOrdinal = create_device_part_ordinal(elemPart2);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities2, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, two_elements_with_same_ranked_part_remove_part_from_one_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);

  const stk::mesh::EntityId elemId1 = 1;
  const stk::mesh::EntityId elemId2 = 2;

  stk::mesh::PartVector elemPart{&part1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem1 = stk::mesh::declare_element(*m_bulk, elemPart, elemId1, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  const stk::mesh::Entity elem2 = stk::mesh::declare_element(*m_bulk, elemPart, elemId2, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  HostEntitiesType hostEntities1("hostEntities1", 1);
  hostEntities1(0) = elem1;
  auto devEntities1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostEntities1);

  HostEntitiesType hostEntities2("hostEntities2", 1);
  hostEntities2(0) = elem2;
  auto devEntities2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostEntities2);

  HostEntitiesType hostNodes1("hostNodes", 2);
  hostNodes1(0) = m_bulk->begin_nodes(elem1)[0];
  hostNodes1(1) = m_bulk->begin_nodes(elem1)[1];
  auto devNodes1 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes1);

  HostEntitiesType hostNodes2("hostNodes", 4);
  hostNodes2(0) = m_bulk->begin_nodes(elem2)[0];
  hostNodes2(1) = m_bulk->begin_nodes(elem2)[1];
  hostNodes2(2) = m_bulk->begin_nodes(elem2)[2];
  hostNodes2(3) = m_bulk->begin_nodes(elem2)[3];
  auto devNodes2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes2);

  const unsigned expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities1,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities2,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  unsigned expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes1,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
  expectedNumConnectedElems = 2;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes2,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = part1.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities1, devAddParts, devRemoveParts));

  auto unexpectedDevicePartOrdinal = create_device_part_ordinal(elemPart);
  check_device_buckets_part_ordinal_does_not_contain_part(ngpMesh, stk::topology::ELEM_RANK, devEntities1, unexpectedDevicePartOrdinal);
  check_device_buckets_part_ordinal_does_not_contain_part(ngpMesh, stk::topology::NODE_RANK, devNodes1, unexpectedDevicePartOrdinal);

  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities1,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devEntities2,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  auto expectedDevicePartOrdinal = create_device_part_ordinal(elemPart);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devEntities2, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes2, expectedDevicePartOrdinal);

  expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes1,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
  expectedNumConnectedElems = 2;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes2,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, two_elements_with_same_ranked_part_remove_ranked_part_from_one_face)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                               "|sideset:name=internal_surface; data=1,6; split=block";
  stk::io::fill_mesh(meshSpec, *m_bulk);

  auto block1Part = m_meta->get_part("block_1");
  auto internalPart = m_meta->get_part("INTERNAL_SURFACE_BLOCK_1_QUAD4");

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  stk::mesh::EntityVector entities;
  stk::mesh::get_entities(*m_bulk, stk::topology::ELEM_RANK, entities);
  auto elem1 = entities[0];
  auto elem2 = entities[1];
  auto face = m_bulk->begin_faces(elem1)[0];

  HostEntitiesType hostElems("hostElems", 2);
  hostElems(0) = elem1;
  hostElems(1) = elem2;
  auto devElems = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostElems);

  HostEntitiesType hostFaces("hostFace", 1);
  hostFaces(0) = face;
  auto devFaces = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostFaces);

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->begin_nodes(elem2)[0];
  hostNodes(1) = m_bulk->begin_nodes(elem2)[1];
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  DevicePartOrdinalsType devAddParts("", 0);
  DevicePartOrdinalsType devRemoveParts("", 1);
  auto hostRemoveParts = Kokkos::create_mirror_view(devRemoveParts);
  hostRemoveParts(0) = block1Part->mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRemoveParts, hostRemoveParts);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devFaces, devAddParts, devRemoveParts));

  auto expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{block1Part});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::ELEM_RANK, devElems, expectedDevicePartOrdinal);

  expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{internalPart});
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::FACE_RANK, devFaces, expectedDevicePartOrdinal);
  check_device_entity_part_ordinal_match(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);
}

// void print_host_entities(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank)
// {
//   const stk::mesh::BucketVector& bkts = bulk.buckets(rank);
//   for(const stk::mesh::Bucket* bptr : bkts) {
//     auto parts = bptr->supersets();
//     std::cout<<rank<<" bkt "<<bptr->bucket_id()<<" parts: ";
//     for(const auto* part : parts) {
//       std::cout<<part->name()<<" ";
//     }
//     std::cout<<std::endl<<" entities: ";
//     for(stk::mesh::Entity ent : *bptr) {
//       std::cout<<bulk.identifier(ent)<<" ";
//     }
//     std::cout<<std::endl;
//   }
// }

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, two_elements_with_same_ranked_part_add_new_ranked_part_to_one_element_and_face)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 8);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                               "|sideset:name=internal_surface; data=1,6; split=block";
  stk::io::fill_mesh(meshSpec, *m_bulk);

  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);
  stk::mesh::Part& part3 = m_meta->declare_part_with_topology("nodePart3", stk::topology::NODE);

  auto block1Part = m_meta->get_part("block_1");
  auto internalPart = m_meta->get_part("INTERNAL_SURFACE_BLOCK_1_QUAD4");

  // print_host_entities(*m_bulk, stk::topology::ELEM_RANK);
  // print_host_entities(*m_bulk, stk::topology::FACE_RANK);
  // print_host_entities(*m_bulk, stk::topology::NODE_RANK);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto elem1 = m_bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto face = m_bulk->begin_faces(elem1)[0];

  HostEntitiesType hostElem("hostElems", 1);
  hostElem(0) = elem1;
  auto devElem = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostElem);

  HostEntitiesType hostFace("hostFace", 1);
  hostFace(0) = face;
  auto devFace = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostFace);


  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->begin_nodes(elem1)[4];
  hostNodes(1) = m_bulk->begin_nodes(elem1)[5];
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  HostEntitiesType hostNodes2("hostNodes", 2);
  hostNodes2(0) = m_bulk->begin_nodes(elem1)[0];
  hostNodes2(1) = m_bulk->begin_nodes(elem1)[1];
  auto devNodes2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes2);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);
  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = part2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  unsigned expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devElem,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  expectedNumConnectedNodes = 4;
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, devFace,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  unsigned expectedNumConnectedElems = 2;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, devFace,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes2,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devElem, devAddParts, devRemoveParts));

  auto expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{block1Part, &part2, internalPart});
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, devFace, expectedDevicePartOrdinal);
  check_device_entity_has_parts(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);

  hostAddParts(0) = part3.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devNodes2, devAddParts, devRemoveParts));

  expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{block1Part, &part2});
  check_device_entity_has_parts(ngpMesh, stk::topology::ELEM_RANK, devElem, expectedDevicePartOrdinal);
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, devFace, expectedDevicePartOrdinal);
  check_device_entity_has_parts(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);

  expectedNumConnectedNodes = 8;
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, devElem,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  expectedNumConnectedNodes = 4;
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, devFace,
                            stk::topology::NODE_RANK, expectedNumConnectedNodes);

  expectedDevicePartOrdinal = create_device_part_ordinal(stk::mesh::PartVector{block1Part, &part2, internalPart});
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, devFace, expectedDevicePartOrdinal);
  check_device_entity_has_parts(ngpMesh, stk::topology::NODE_RANK, devNodes, expectedDevicePartOrdinal);

  expectedNumConnectedElems = 2;
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, devFace,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  expectedNumConnectedElems = 1;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes2,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);

  expectedNumConnectedElems = 2;
  check_device_connectivity(ngpMesh, stk::topology::NODE_RANK, devNodes,
                            stk::topology::ELEM_RANK, expectedNumConnectedElems);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, nodal_field_data_two_elements_with_same_ranked_part_add_new_ranked_part_to_one_element_and_face)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 8);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                               "|sideset:name=internal_surface; data=1,6; split=block";
  stk::io::fill_mesh(meshSpec, *m_bulk);

  m_meta->enable_late_fields();
  stk::mesh::Field<double>& nodalField = m_meta->declare_field<double>(stk::topology::NODE_RANK, "myNodalField");

  const unsigned numFieldComponents = 1;
  stk::mesh::put_field_on_mesh(nodalField, m_meta->universal_part(), numFieldComponents, nullptr);

  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  // print_host_entities(*m_bulk, stk::topology::ELEM_RANK);
  // print_host_entities(*m_bulk, stk::topology::FACE_RANK);
  // print_host_entities(*m_bulk, stk::topology::NODE_RANK);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto elem1 = m_bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto face = m_bulk->begin_faces(elem1)[0];

  HostEntitiesType hostElem("hostElems", 1);
  hostElem(0) = elem1;
  auto devElem = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostElem);

  HostEntitiesType hostFace("hostFace", 1);
  hostFace(0) = face;
  auto devFace = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostFace);

  HostEntitiesType hostNodes("hostNodes", 2);
  hostNodes(0) = m_bulk->begin_nodes(elem1)[4];
  hostNodes(1) = m_bulk->begin_nodes(elem1)[5];
  init_host_field_data(*m_bulk, nodalField, hostNodes);
  auto devNodes = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, nodalField, devNodes);

  hostNodes(0) = m_bulk->begin_nodes(elem1)[0];
  hostNodes(1) = m_bulk->begin_nodes(elem1)[1];
  init_host_field_data(*m_bulk, nodalField, hostNodes);
  auto devNodes2 = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, hostNodes);

  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, nodalField, devNodes2);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);
  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = part2.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devElem, devAddParts, devRemoveParts));

  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, nodalField, devNodes2);

  hostNodes(0) = m_bulk->begin_nodes(elem1)[4];
  hostNodes(1) = m_bulk->begin_nodes(elem1)[5];
  check_device_entity_field_data_is_id(ngpMesh, stk::topology::NODE_RANK, nodalField, devNodes);
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, test_repeated_part_addition)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  auto numElemsInZ = stk::unit_test_util::get_command_line_option("-n", 1);
  auto numIters = stk::unit_test_util::get_command_line_option("-i", 3);

  auto genMeshDesc = "generated:10x10x" + std::to_string(numElemsInZ);
  stk::io::fill_mesh(genMeshDesc, *m_bulk);

  stk::mesh::Part& partToAdd = m_meta->declare_part_with_topology("newPart", stk::topology::HEX_8);

  stk::mesh::EntityVector entityVector;
  stk::mesh::get_entities(*m_bulk, stk::topology::ELEM_RANK, stk::mesh::Selector(*m_meta->get_part("block_1")), entityVector);

  DeviceEntitiesType devEntities("deviceEntities", entityVector.size());
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  for (unsigned i = 0; i < entityVector.size(); ++i) {
    hostEntities(i) = entityVector[i];
  }
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRemoveParts("", 0);

  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = partToAdd.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  for (int i = 0; i < numIters; ++i) {
    EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, devRemoveParts));
  }
}

TEST_F(NgpBatchChangeEntityPartsInducedPartMembership, test_repeated_part_addition_and_removal)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  auto numElemsInZ = stk::unit_test_util::get_command_line_option("-n", 1);
  auto numIters = stk::unit_test_util::get_command_line_option("-i", 3);

  auto genMeshDesc = "generated:10x10x" + std::to_string(numElemsInZ);
  stk::io::fill_mesh(genMeshDesc, *m_bulk);

  stk::mesh::Part& newPart = m_meta->declare_part_with_topology("newPart", stk::topology::HEX_8);

  stk::mesh::EntityVector entityVector;
  stk::mesh::get_entities(*m_bulk, stk::topology::ELEM_RANK, stk::mesh::Selector(*m_meta->get_part("block_1")), entityVector);

  DeviceEntitiesType devEntities("deviceEntities", entityVector.size());
  auto hostEntities = Kokkos::create_mirror_view(devEntities);
  for (unsigned i = 0; i < entityVector.size(); ++i) {
    hostEntities(i) = entityVector[i];
  }
  Kokkos::deep_copy(devEntities, hostEntities);

  DevicePartOrdinalsType emptyParts("", 0);
  DevicePartOrdinalsType devAddParts("", 1);
  DevicePartOrdinalsType devRmParts("", 1);

  auto hostAddParts = Kokkos::create_mirror_view(devAddParts);
  hostAddParts(0) = newPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devAddParts, hostAddParts);

  auto hostRmParts = Kokkos::create_mirror_view(devRmParts);
  hostRmParts(0) = newPart.mesh_meta_data_ordinal();
  Kokkos::deep_copy(devRmParts, hostRmParts);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  for (int i = 0; i < numIters; ++i) {
    EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, devAddParts, emptyParts));
    EXPECT_NO_THROW(ngpMesh.impl_batch_change_entity_parts_with_inducible_parts(devEntities, emptyParts, devRmParts));
  }
}

namespace {
void check_entities_valid(stk::mesh::NgpMesh& deviceMesh)
{
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    Kokkos::parallel_for(
      "check_entities_valid" ,
      stk::ngp::RangePolicy<stk::ngp::DeviceSpace::exec_space>(0, deviceBucketRepo.num_buckets(rank)),
      KOKKOS_LAMBDA(unsigned bucket_id)
      {
        DeviceBucket& bucket = *(deviceBucketRepo.get_bucket(rank, bucket_id));
        for (unsigned i=0; i < bucket.size(); ++i)
        {
          NGP_EXPECT_NE(bucket[i], stk::mesh::Entity());
        }
      }
    );
  }
}

void check_entities_unique(stk::mesh::NgpMesh& deviceMesh)
{
  using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>, stk::ngp::DeviceSpace::exec_space>;

  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();


  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    Policy policy({0, 0}, {deviceBucketRepo.num_buckets(rank), deviceBucketRepo.num_buckets(rank)});

    Kokkos::parallel_for(
      "check_entities_unique" ,
      policy,
      KOKKOS_LAMBDA(unsigned bucketId_i, unsigned bucketId_j)
      {
        if (bucketId_i == bucketId_j)
        {
          return;
        }
        DeviceBucket& bucket_i = *(deviceBucketRepo.get_bucket(rank, bucketId_i));
        DeviceBucket& bucket_j = *(deviceBucketRepo.get_bucket(rank, bucketId_j));
        for (unsigned i=0; i < bucket_i.size(); ++i)
        {
          for (unsigned j=0; j < bucket_j.size(); ++j)
          {
            NGP_EXPECT_NE(bucket_i[i], bucket_j[j]);
          }
        }
      }
    );
  }
}

}

TEST_F(NgpBatchChangeEntityPartsTwoBlocks, MoveElementOne)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  move_entity_on_device({0});

  check_entities_valid(deviceMesh);
  check_entities_unique(deviceMesh);
  EXPECT_EQ(deviceBucketRepo.get_bucket(stk::topology::ELEM_RANK, 0)->size(), 1u);
  EXPECT_EQ(deviceBucketRepo.get_bucket(stk::topology::ELEM_RANK, 1)->size(), 1u);
}

TEST_F(NgpBatchChangeEntityPartsTwoBlocks, MoveElementTwo)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  move_entity_on_device({1});

  check_entities_valid(deviceMesh);
  check_entities_unique(deviceMesh);
  EXPECT_EQ(deviceBucketRepo.get_bucket(stk::topology::ELEM_RANK, 0)->size(), 1u);
  EXPECT_EQ(deviceBucketRepo.get_bucket(stk::topology::ELEM_RANK, 1)->size(), 1u);
}

TEST_F(NgpBatchChangeEntityPartsTwoBlocks, MoveBothElements)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  move_entity_on_device({0, 1});

  check_entities_valid(deviceMesh);
  check_entities_unique(deviceMesh);
  EXPECT_EQ(deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK), 1u);
  EXPECT_EQ(deviceBucketRepo.get_bucket(stk::topology::ELEM_RANK, 0)->size(), 2u);
}

}  // namespace
#endif
