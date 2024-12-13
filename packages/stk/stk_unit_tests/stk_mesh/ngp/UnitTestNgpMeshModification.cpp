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
#include <stk_mesh/base/NgpMesh.hpp>
#include "ngp/NgpUnitTestUtils.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SkinMesh.hpp"

namespace
{
using ngp_unit_test_utils::check_bucket_layout;

class NgpBatchChangeEntityParts : public ::testing::Test
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
    EXPECT_EQ(ngpMesh.need_sync_to_host(), true);
  }
  else {
    EXPECT_EQ(ngpMesh.need_sync_to_host(), false);  // If host build, HostMesh can't ever be stale
  }
}

template <typename MeshType>
void confirm_host_mesh_is_synchronized_from_device(const MeshType& ngpMesh)
{
  EXPECT_EQ(ngpMesh.need_sync_to_host(), false);
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
}


TEST_F(NgpBatchChangeEntityParts, addPartToNode_host)
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

TEST_F(NgpBatchChangeEntityParts, addPartToNode_ngpHost)
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

  hostMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchChangeEntityParts, addPartToNode_ngpDevice)
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

  ngpMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}


TEST_F(NgpBatchChangeEntityParts, removePartFromNode_host)
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

TEST_F(NgpBatchChangeEntityParts, removePartFromNode_ngpHost)
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

  hostMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchChangeEntityParts, removePartFromNode_ngpDevice)
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

  ngpMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}


TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_host)
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

TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_ngpHost)
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

  hostMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(hostMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchChangeEntityParts, addAndRemovePartFromNode_ngpDevice)
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

  ngpMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchChangeEntityParts, multipleDeviceMeshMods)
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

  ngpMesh.batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  ngpMesh.sync_to_host();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchChangeEntityParts, failedHostAccessAfterDeviceMeshMod)
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

TEST_F(NgpBatchChangeEntityParts, impl_addPartToNode_ngpDevice)
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

  ngpMesh.impl_batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
//  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);
//
//  ngpMesh.sync_to_host();
//  confirm_host_mesh_is_synchronized_from_device(ngpMesh);
//
//  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId}}}, stk::topology::NODE_RANK);
}

}  // namespace
