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
#include <algorithm>
#include <iterator>
#include <map>
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
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"

#include "ngp/UnitTestNgpMeshModificationUtils.hpp"

#ifdef STK_USE_DEVICE_MESH
namespace
{
class NgpBatchDestroyEntities : public NgpBatchDeclareDestroyEntities
{
public:
  NgpBatchDestroyEntities() {}

  HostEntitiesType make_host_entities(stk::topology::rank_t rank,
                                      const std::vector<unsigned>& ids)
  {
    HostEntitiesType entities("hostEntities", ids.size());
    for (size_t i = 0; i < ids.size(); ++i) {
      entities(i) = m_bulk->get_entity(rank, ids[i]);
    }
    return entities;
  }
  
  DeviceEntitiesType make_device_entities(stk::topology::rank_t rank,
                                          const std::vector<unsigned>& ids)
  {
    DeviceEntitiesType entities("deviceEntities", ids.size());
    auto hostEntities = Kokkos::create_mirror_view(entities);
    for (size_t i = 0; i < ids.size(); ++i) {
      hostEntities(i) = m_bulk->get_entity(rank, ids[i]);
    }
    Kokkos::deep_copy(entities, hostEntities);
    return entities;
  }

  void check_host_connectivity_by_id(stk::mesh::HostMesh& hostMesh,
                                     stk::topology::rank_t fromRank,
                                     const std::vector<unsigned>& fromIds,
                                     stk::topology::rank_t toRank,
                                     unsigned expectedCount)
  {
    HostEntitiesType entities = make_host_entities(fromRank, fromIds);
    check_host_connectivity(hostMesh, fromRank, entities, toRank, expectedCount);
  }

  void check_device_connectivity_by_id(stk::mesh::NgpMesh& ngpMesh,
                                       stk::topology::rank_t fromRank,
                                       const std::vector<unsigned>& fromIds,
                                       stk::topology::rank_t toRank,
                                       unsigned expectedCount)
  {
    DeviceEntitiesType entities = make_device_entities(fromRank, fromIds);
    check_device_connectivity(ngpMesh, fromRank, entities, toRank, expectedCount);
  }

  void destroy_on_host(const std::vector<stk::mesh::Entity>& entitiesToDestroy)
  {
    HostEntitiesType entities("hostEntities", entitiesToDestroy.size());
    fill_views(entities, entitiesToDestroy);

    stk::mesh::HostMesh hostMesh(*m_bulk);
    hostMesh.batch_destroy_entities(entities);
    confirm_host_mesh_is_synchronized_from_device(hostMesh);

    hostMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(hostMesh);
  }

  void destroy_on_device(const std::vector<stk::mesh::Entity>& entitiesToDestroy, [[maybe_unused]] stk::topology::rank_t rank)
  {
    DeviceEntitiesType entities("deviceEntities", entitiesToDestroy.size());
    fill_views(entities, entitiesToDestroy);

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_destroy_entities(entities);
    confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

    DevicePartOrdinalsType addPartOrdinals("deviceAddParts", 0);
    DevicePartOrdinalsType removePartOrdinals("deviceRemoveParts", 0);
    check_entity_parts_on_device(ngpMesh, entities, addPartOrdinals, removePartOrdinals, rank);

    ngpMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(ngpMesh);
  }

  std::vector<bool> destroy_on_host_with_results(const std::vector<stk::mesh::Entity>& entitiesToDestroy)
  {
    HostEntitiesType entities("hostEntities", entitiesToDestroy.size());
    fill_views(entities, entitiesToDestroy);

    Kokkos::View<bool*, stk::ngp::HostExecSpace> wasDestroyed("wasDestroyed", entitiesToDestroy.size());

    stk::mesh::HostMesh hostMesh(*m_bulk);
    hostMesh.batch_destroy_entities(entities, wasDestroyed);
    hostMesh.update_bulk_data();

    std::vector<bool> results(entitiesToDestroy.size());
    for (size_t i = 0; i < results.size(); ++i) {
      results[i] = wasDestroyed(i);
    }
    return results;
  }

  std::vector<bool> destroy_on_device_with_results(const std::vector<stk::mesh::Entity>& entitiesToDestroy)
  {
    DeviceEntitiesType entities("deviceEntities", entitiesToDestroy.size());
    fill_views(entities, entitiesToDestroy);

    Kokkos::View<bool*, stk::ngp::MemSpace> wasDestroyed("wasDestroyed", entitiesToDestroy.size());

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_destroy_entities(entities, wasDestroyed);
    ngpMesh.update_bulk_data();

    auto hostWasDestroyed = Kokkos::create_mirror_view(wasDestroyed);
    Kokkos::deep_copy(hostWasDestroyed, wasDestroyed);

    std::vector<bool> results(entitiesToDestroy.size());
    for (size_t i = 0; i < results.size(); ++i) {
      results[i] = hostWasDestroyed(i);
    }
    return results;
  }

  bool is_valid_by_id(stk::topology::rank_t rank, unsigned id)
  {
    return m_bulk->is_valid(m_bulk->get_entity(rank, id));
  }
};

class NgpBatchDestroyNodes : public NgpBatchDestroyEntities
{
public:
  NgpBatchDestroyNodes() {}

  void destroy_on_device(const std::vector<stk::mesh::Entity>& entitiesToDestroy)
  {
    NgpBatchDestroyEntities::destroy_on_device(entitiesToDestroy, stk::topology::NODE_RANK);
  }

  static constexpr unsigned nodeId1 = 1;
  static constexpr unsigned nodeId2 = 2;
  static constexpr unsigned nodeId3 = 3;
};

class NgpBatchDestroyElements : public NgpBatchDestroyEntities
{
public:
  NgpBatchDestroyElements() {}

  stk::mesh::Entity declare_hex8_element(stk::mesh::PartVector& parts,
                                         unsigned elemId,
                                         const std::vector<unsigned>& nodeIds)
  {
    return stk::mesh::declare_element(
        *m_bulk,
        parts,
        elemId,
        stk::mesh::EntityIdVector(nodeIds.begin(), nodeIds.end()));
  }

  void destroy_on_device(const std::vector<stk::mesh::Entity>& entitiesToDestroy)
  {
    NgpBatchDestroyEntities::destroy_on_device(entitiesToDestroy, stk::topology::ELEM_RANK);
  }

  stk::mesh::Entity build_single_hex()
  {
    build_empty_mesh(1, 1);

    stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
    stk::mesh::PartVector parts{&elemPart1};

    m_bulk->modification_begin();
    const stk::mesh::Entity elem = declare_hex8_element(parts, elemId1, {1,2,3,4,5,6,7,8});
    m_bulk->modification_end();

    return elem;
  }

  stk::mesh::Entity build_one_hex_with_edges_and_faces()
  {
    build_empty_mesh(1, 1);

    stk::mesh::Part& elemPart = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
    stk::mesh::Part& facePart = m_meta->declare_part_with_topology("facePart", stk::topology::QUAD_4);
    stk::mesh::Part& edgePart = m_meta->declare_part_with_topology("edgePart", stk::topology::LINE_2);
    stk::mesh::PartVector parts{&elemPart};

    m_bulk->modification_begin();
    const stk::mesh::Entity elem = declare_hex8_element(parts, elemId1, {1,2,3,4,5,6,7,8});
    m_bulk->modification_end();

    constexpr bool connectFacesToEdges = true;
    stk::mesh::create_edges(*m_bulk, m_meta->universal_part(), &edgePart);
    stk::mesh::create_all_sides(*m_bulk, m_meta->universal_part(), {&facePart}, connectFacesToEdges);

    return elem;
  }

  size_t num_entities(stk::topology::rank_t rank)
  {
    return stk::mesh::count_entities(*m_bulk, rank, m_meta->universal_part());
  }

  static constexpr unsigned elemId1 = 1;
  static constexpr unsigned elemId2 = 2;
};

template <typename DeviceMeshType>
unsigned device_fast_mesh_index_bucket_id(DeviceMeshType& ngpMesh, stk::mesh::Entity entity)
{
  Kokkos::View<unsigned*, stk::ngp::UVMMemSpace> result("fastMeshIndexBucketId", 1);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const int /*idx*/) {
      result(0) = ngpMesh.fast_mesh_index(entity).bucket_id;
    });
  Kokkos::fence();
  return result(0);
}

template <typename DeviceMeshType>
stk::mesh::EntityKey device_entity_key(DeviceMeshType& ngpMesh, stk::mesh::Entity entity)
{
  Kokkos::View<stk::mesh::EntityKey*, stk::ngp::UVMMemSpace> result("deviceEntityKey", 1);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const int /*idx*/) {
      result(0) = ngpMesh.entity_key(entity);
    });
  Kokkos::fence();
  return result(0);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyOneNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  destroy_on_host({node1});

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyOneNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  destroy_on_device({node1});

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({node1, node2});

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodes_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({node1, node2});

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyEntityInLastBucket_invalidatesFastMeshIndex_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);  // bucket capacity 1 => each node lands in its own bucket

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  DeviceEntitiesType entitiesToDestroy("entitiesToDestroy", 1);
  fill_views(entitiesToDestroy, std::vector<stk::mesh::Entity>{node2});

  ngpMesh.batch_destroy_entities(entitiesToDestroy);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

  const unsigned numBucketsAfterDestroy = ngpMesh.num_buckets(stk::topology::NODE_RANK);
  EXPECT_EQ(numBucketsAfterDestroy, 1u);  // the last bucket was invalidated/removed

  const unsigned staleBucketId = device_fast_mesh_index_bucket_id(ngpMesh, node2);

  EXPECT_EQ(staleBucketId, stk::mesh::INVALID_BUCKET_ID);

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyNode_invalidatesDeviceEntityKey_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);  // bucket capacity 1 => each node lands in its own bucket

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part2});

  const stk::mesh::EntityKey node1Key = m_bulk->entity_key(node1);
  const stk::mesh::EntityKey node2Key = m_bulk->entity_key(node2);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  EXPECT_EQ(device_entity_key(ngpMesh, node2), node2Key);

  DeviceEntitiesType entitiesToDestroy("entitiesToDestroy", 1);
  fill_views(entitiesToDestroy, std::vector<stk::mesh::Entity>{node2});

  ngpMesh.batch_destroy_entities(entitiesToDestroy);
  confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);  // device is ahead; not yet synced to host

  const stk::mesh::EntityKey destroyedKey = device_entity_key(ngpMesh, node2);
  EXPECT_FALSE(destroyedKey.is_valid());

  EXPECT_EQ(device_entity_key(ngpMesh, node1), node1Key);

  EXPECT_FALSE(ngpMesh.get_entity(node2Key).is_local_offset_valid());
  EXPECT_EQ(ngpMesh.get_entity(node1Key), node1);

  ngpMesh.update_bulk_data();
  confirm_host_mesh_is_synchronized_from_device(ngpMesh);

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesOneRemains_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 =
      m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 =
      m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part1});
  const stk::mesh::Entity node3 = create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part1"}, {nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({node1, node3});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesOneRemains_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 =
      m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 =
      m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part1});
  const stk::mesh::Entity node3 = create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part1"}, {nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({node1, node3});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyOneNodeOneRemainsBiggerBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({node1});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyOneNodeOneRemainsBiggerBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId1}},
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({node1});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part2"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesInSameBucketOneRemainsBiggerBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1});
  create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1", "part1"}, {nodeId1, nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({node1, node2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesInSameBucketOneRemainsBiggerBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  const stk::mesh::Entity node2 = create_node(*m_bulk, nodeId2, {&part1});
  create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1", "part1"}, {nodeId1, nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({node1, node2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesOneRemainsBiggerBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part1});
  const stk::mesh::Entity node3 = create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1", "part1"}, {nodeId1, nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({node1, node3});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyNodes, destroyTwoNodesOneRemainsBiggerBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(3, 3);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  const stk::mesh::Entity node1 = create_node(*m_bulk, nodeId1, {&part1});
  create_node(*m_bulk, nodeId2, {&part1});
  const stk::mesh::Entity node3 = create_node(*m_bulk, nodeId3, {&part2});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1", "part1"}, {nodeId1, nodeId2}},
                        {{"part2"}, {nodeId3}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({node1, node3});

  check_bucket_layout(*m_bulk,
                      {
                        {{"part1"}, {nodeId2}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, destroyElement_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  stk::mesh::HostMesh hostMesh(*m_bulk);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::ELEM_RANK, {elemId1},
                                stk::topology::NODE_RANK, 8);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                stk::topology::ELEM_RANK, 1);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({elem});

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {1}},
                        {{}, {2}},
                        {{}, {3}},
                        {{}, {4}},
                        {{}, {5}},
                        {{}, {6}},
                        {{}, {7}},
                        {{}, {8}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, destroyElement_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::ELEM_RANK, {elemId1},
                                  stk::topology::NODE_RANK, 8);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                  stk::topology::ELEM_RANK, 1);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({elem});

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(0u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                  stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {1}},
                        {{}, {2}},
                        {{}, {3}},
                        {{}, {4}},
                        {{}, {5}},
                        {{}, {6}},
                        {{}, {7}},
                        {{}, {8}}
                      },
                      stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, destroyOneOfTwoElements_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem1 = declare_hex8_element(parts, elemId1, {1,2,3,4,5,6,7,8});
  declare_hex8_element(parts, elemId2, {5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  stk::mesh::HostMesh hostMesh(*m_bulk);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::ELEM_RANK, {elemId1, elemId2},
                                stk::topology::NODE_RANK, 8);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4,9,10,11,12},
                                stk::topology::ELEM_RANK, 1);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {5,6,7,8},
                                stk::topology::ELEM_RANK, 2);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}},
                        {{"elemPart1"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}},
                        {{"elemPart1"}, {9}},
                        {{"elemPart1"}, {10}},
                        {{"elemPart1"}, {11}},
                        {{"elemPart1"}, {12}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({elem1});

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::ELEM_RANK, {elemId2},
                                stk::topology::NODE_RANK, 8);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {5,6,7,8,9,10,11,12},
                                stk::topology::ELEM_RANK, 1);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4},
                                stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {1}},
                        {{}, {2}},
                        {{}, {3}},
                        {{}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}},
                        {{"elemPart1"}, {9}},
                        {{"elemPart1"}, {10}},
                        {{"elemPart1"}, {11}},
                        {{"elemPart1"}, {12}}
                      },
                      stk::topology::NODE_RANK);

  EXPECT_EQ(1u, hostMesh.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(12u, hostMesh.num_buckets(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, destroyOneOfTwoElements_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 =
      m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::PartVector parts{&elemPart1};

  m_bulk->modification_begin();
  const stk::mesh::Entity elem1 = declare_hex8_element(parts, elemId1, {1,2,3,4,5,6,7,8});
  declare_hex8_element(parts, elemId2, {5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::ELEM_RANK, {elemId1, elemId2},
                                  stk::topology::NODE_RANK, 8);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4,9,10,11,12},
                                  stk::topology::ELEM_RANK, 1);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {5,6,7,8},
                                  stk::topology::ELEM_RANK, 2);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}},
                        {{"elemPart1"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}},
                        {{"elemPart1"}, {9}},
                        {{"elemPart1"}, {10}},
                        {{"elemPart1"}, {11}},
                        {{"elemPart1"}, {12}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({elem1});

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::ELEM_RANK, {elemId2},
                                  stk::topology::NODE_RANK, 8);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {5,6,7,8,9,10,11,12},
                                  stk::topology::ELEM_RANK, 1);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4},
                                  stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {1}},
                        {{}, {2}},
                        {{}, {3}},
                        {{}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}},
                        {{"elemPart1"}, {9}},
                        {{"elemPart1"}, {10}},
                        {{"elemPart1"}, {11}},
                        {{"elemPart1"}, {12}}
                      },
                      stk::topology::NODE_RANK);

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(12u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, destroyOneOfTwoElements_differentParts_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  stk::mesh::PartVector parts1{&elemPart1};
  stk::mesh::PartVector parts2{&elemPart2};

  m_bulk->modification_begin();
  declare_hex8_element(parts1, elemId1, {1,2,3,4,5,6,7,8});
  const stk::mesh::Entity elem2 = declare_hex8_element(parts2, elemId2, {5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  stk::mesh::HostMesh hostMesh(*m_bulk);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::ELEM_RANK, {elemId1, elemId2},
                                stk::topology::NODE_RANK, 8);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4,9,10,11,12},
                                stk::topology::ELEM_RANK, 1);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {5,6,7,8},
                                stk::topology::ELEM_RANK, 2);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}},
                        {{"elemPart2"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart2"}, {9}},
                        {{"elemPart2"}, {10}},
                        {{"elemPart2"}, {11}},
                        {{"elemPart2"}, {12}},
                        {{"elemPart1", "elemPart2"}, {5}},
                        {{"elemPart1", "elemPart2"}, {6}},
                        {{"elemPart1", "elemPart2"}, {7}},
                        {{"elemPart1", "elemPart2"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_host({elem2});

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::ELEM_RANK, {elemId1},
                                stk::topology::NODE_RANK, 8);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                stk::topology::ELEM_RANK, 1);

  check_host_connectivity_by_id(hostMesh,
                                stk::topology::NODE_RANK, {9,10,11,12},
                                stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {9}},
                        {{}, {10}},
                        {{}, {11}},
                        {{}, {12}},
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  EXPECT_EQ(1u, hostMesh.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(12u, hostMesh.num_buckets(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, destroyOneOfTwoElements_differentParts_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& elemPart1 = m_meta->declare_part_with_topology("elemPart1", stk::topology::HEX_8);
  stk::mesh::Part& elemPart2 = m_meta->declare_part_with_topology("elemPart2", stk::topology::HEX_8);

  stk::mesh::PartVector parts1{&elemPart1};
  stk::mesh::PartVector parts2{&elemPart2};

  m_bulk->modification_begin();
  declare_hex8_element(parts1, elemId1, {1,2,3,4,5,6,7,8});
  const stk::mesh::Entity elem2 = declare_hex8_element(parts2, elemId2, {5,6,7,8,9,10,11,12});
  m_bulk->modification_end();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::ELEM_RANK, {elemId1, elemId2},
                                  stk::topology::NODE_RANK, 8);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4,9,10,11,12},
                                  stk::topology::ELEM_RANK, 1);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {5,6,7,8},
                                  stk::topology::ELEM_RANK, 2);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}},
                        {{"elemPart2"}, {elemId2}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart2"}, {9}},
                        {{"elemPart2"}, {10}},
                        {{"elemPart2"}, {11}},
                        {{"elemPart2"}, {12}},
                        {{"elemPart1", "elemPart2"}, {5}},
                        {{"elemPart1", "elemPart2"}, {6}},
                        {{"elemPart1", "elemPart2"}, {7}},
                        {{"elemPart1", "elemPart2"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  destroy_on_device({elem2});

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::ELEM_RANK, {elemId1},
                                  stk::topology::NODE_RANK, 8);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {1,2,3,4,5,6,7,8},
                                  stk::topology::ELEM_RANK, 1);

  check_device_connectivity_by_id(ngpMesh,
                                  stk::topology::NODE_RANK, {9,10,11,12},
                                  stk::topology::ELEM_RANK, 0);

  check_bucket_layout(*m_bulk,
                      {
                        {{"elemPart1"}, {elemId1}}
                      },
                      stk::topology::ELEM_RANK);

  check_bucket_layout(*m_bulk,
                      {
                        {{}, {9}},
                        {{}, {10}},
                        {{}, {11}},
                        {{}, {12}},
                        {{"elemPart1"}, {1}},
                        {{"elemPart1"}, {2}},
                        {{"elemPart1"}, {3}},
                        {{"elemPart1"}, {4}},
                        {{"elemPart1"}, {5}},
                        {{"elemPart1"}, {6}},
                        {{"elemPart1"}, {7}},
                        {{"elemPart1"}, {8}}
                      },
                      stk::topology::NODE_RANK);

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EXPECT_EQ(1u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(1u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));

  EXPECT_EQ(12u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(2u, deviceBucketRepo.num_partitions(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, reportsDestroyResult_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const std::vector<bool> results = destroy_on_host_with_results({elem});

  EXPECT_EQ((std::vector<bool>{true}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, reportsDestroyResult_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const std::vector<bool> results = destroy_on_device_with_results({elem});

  EXPECT_EQ((std::vector<bool>{true}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, skipsInvalidEntity_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const stk::mesh::Entity invalidEntity;  // default-constructed => invalid

  const std::vector<bool> results = destroy_on_host_with_results({elem, invalidEntity});

  EXPECT_EQ((std::vector<bool>{true, false}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, skipsInvalidEntity_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const stk::mesh::Entity invalidEntity;  // default-constructed => invalid

  const std::vector<bool> results = destroy_on_device_with_results({elem, invalidEntity});

  EXPECT_EQ((std::vector<bool>{true, false}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, skipsEntityWithUpwardConnectivity_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_single_hex();

  const stk::mesh::Entity node1 = m_bulk->get_entity(stk::topology::NODE_RANK, 1);

  const std::vector<bool> results = destroy_on_host_with_results({node1});

  EXPECT_EQ((std::vector<bool>{false}), results);
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
  EXPECT_TRUE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
}

NGP_TEST_F(NgpBatchDestroyElements, skipsEntityWithUpwardConnectivity_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_single_hex();

  const stk::mesh::Entity node1 = m_bulk->get_entity(stk::topology::NODE_RANK, 1);

  const std::vector<bool> results = destroy_on_device_with_results({node1});

  EXPECT_EQ((std::vector<bool>{false}), results);
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
  EXPECT_TRUE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
}

NGP_TEST_F(NgpBatchDestroyElements, doesNotCascadeParentAndChild_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const stk::mesh::Entity node1 = m_bulk->get_entity(stk::topology::NODE_RANK, 1);

  const std::vector<bool> results = destroy_on_host_with_results({elem, node1});

  EXPECT_EQ((std::vector<bool>{true, false}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
}

NGP_TEST_F(NgpBatchDestroyElements, doesNotCascadeParentAndChild_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  const stk::mesh::Entity node1 = m_bulk->get_entity(stk::topology::NODE_RANK, 1);

  const std::vector<bool> results = destroy_on_device_with_results({elem, node1});

  EXPECT_EQ((std::vector<bool>{true, false}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
}

NGP_TEST_F(NgpBatchDestroyElements, destroysDuplicateEntityOnce_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  // The same entity listed twice must be destroyed exactly once, with every occurrence reported true.
  const std::vector<bool> results = destroy_on_host_with_results({elem, elem});

  EXPECT_EQ((std::vector<bool>{true, true}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);
}

NGP_TEST_F(NgpBatchDestroyElements, destroysDuplicateEntityOnce_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_single_hex();

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  const std::vector<bool> results = destroy_on_device_with_results({elem, elem});

  EXPECT_EQ((std::vector<bool>{true, true}), results);
  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);

  // The element must be removed exactly once; a double-remove would corrupt these counts (or crash).
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  EXPECT_EQ(0u, deviceBucketRepo.num_buckets(stk::topology::ELEM_RANK));
  EXPECT_EQ(0u, deviceBucketRepo.num_partitions(stk::topology::ELEM_RANK));
  EXPECT_EQ(8u, deviceBucketRepo.num_buckets(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, destroysElementWithMultiRankConnectivity_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_one_hex_with_edges_and_faces();

  EXPECT_TRUE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  EXPECT_EQ(1u,  num_entities(stk::topology::ELEM_RANK));
  EXPECT_EQ(6u,  num_entities(stk::topology::FACE_RANK));
  EXPECT_EQ(12u, num_entities(stk::topology::EDGE_RANK));
  EXPECT_EQ(8u,  num_entities(stk::topology::NODE_RANK));

  destroy_on_host({elem});

  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);

  EXPECT_EQ(0u,  num_entities(stk::topology::ELEM_RANK));
  EXPECT_EQ(6u,  num_entities(stk::topology::FACE_RANK));
  EXPECT_EQ(12u, num_entities(stk::topology::EDGE_RANK));
  EXPECT_EQ(8u,  num_entities(stk::topology::NODE_RANK));
}

NGP_TEST_F(NgpBatchDestroyElements, destroysElementWithMultiRankConnectivity_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const stk::mesh::Entity elem = build_one_hex_with_edges_and_faces();

  EXPECT_EQ(1u,  num_entities(stk::topology::ELEM_RANK));
  EXPECT_EQ(6u,  num_entities(stk::topology::FACE_RANK));
  EXPECT_EQ(12u, num_entities(stk::topology::EDGE_RANK));
  EXPECT_EQ(8u,  num_entities(stk::topology::NODE_RANK));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  const unsigned syncCountBefore = ngpMesh.synchronized_count();

  DeviceEntitiesType entities("deviceEntities", 1);
  fill_views(entities, std::vector<stk::mesh::Entity>{elem});
  ngpMesh.batch_destroy_entities(entities);

  EXPECT_EQ(1u, ngpMesh.synchronized_count() - syncCountBefore);

  ngpMesh.update_bulk_data();

  EXPECT_FALSE(is_valid_by_id(stk::topology::ELEM_RANK, elemId1));
  EXPECT_TRUE(is_valid_by_id(stk::topology::NODE_RANK, 1));
  check_bucket_layout(*m_bulk, {}, stk::topology::ELEM_RANK);

  EXPECT_EQ(0u,  num_entities(stk::topology::ELEM_RANK));
  EXPECT_EQ(6u,  num_entities(stk::topology::FACE_RANK));
  EXPECT_EQ(12u, num_entities(stk::topology::EDGE_RANK));
  EXPECT_EQ(8u,  num_entities(stk::topology::NODE_RANK));
}

}  // namespace
#endif
