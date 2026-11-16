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
class NgpBatchDeclareEntities : public NgpBatchDeclareDestroyEntities
{
public:
  NgpBatchDeclareEntities() {};

  std::vector<stk::mesh::PartOrdinal> extract_part_ordinals(const stk::mesh::PartVector& addParts) const {
    std::vector<stk::mesh::PartOrdinal> addPartOrdinals;
    std::ranges::transform(addParts, std::back_inserter(addPartOrdinals), [](const stk::mesh::Part* part) {
      return part->mesh_meta_data_ordinal();
    });
    return addPartOrdinals;
  }
  void declare_on_host(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    Kokkos::View<unsigned*, Kokkos::HostSpace> entityIds("entityIds", entityIdsToDeclare.size());
    fill_views(entityIds, entityIdsToDeclare);

    HostPartOrdinalsType addPartOrdinals("addPartOrdinals", addParts.size());
    fill_views(addPartOrdinals, extract_part_ordinals(addParts));

    stk::mesh::HostMesh hostMesh(*m_bulk);

    HostEntitiesType createdEntities("createdEntities", entityIdsToDeclare.size());
    hostMesh.batch_declare_entities(stk::topology::NODE_RANK, entityIds, addPartOrdinals, createdEntities);
    confirm_host_mesh_is_not_synchronized_from_device(hostMesh);

    hostMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(hostMesh);
  }

  void declare_on_device(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    Kokkos::View<unsigned*> entityIds("entityIds", entityIdsToDeclare.size());
    fill_views(entityIds, entityIdsToDeclare);

    DevicePartOrdinalsType addPartOrdinals("addPartOrdinals", addParts.size());
    fill_views(addPartOrdinals, extract_part_ordinals(addParts));

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

    DeviceEntitiesType createdEntities("createdEntities", entityIdsToDeclare.size());
    ngpMesh.batch_declare_entities(stk::topology::NODE_RANK, entityIds, addPartOrdinals, createdEntities);
    confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

    ngpMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(ngpMesh);
  }
};

class NgpBatchDeclareNodes : public NgpBatchDeclareEntities
{
public:
  NgpBatchDeclareNodes() {};

  void declare_on_host(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    NgpBatchDeclareEntities::declare_on_host(entityIdsToDeclare, addParts);
  }

  void declare_on_device(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    NgpBatchDeclareEntities::declare_on_device(entityIdsToDeclare, addParts);
  }

  static constexpr unsigned nodeId1 = 1;
  static constexpr unsigned nodeId2 = 2;
};

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_OneBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_OneBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_DifferentBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_host({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_DifferentBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_device({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_SameBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_host({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_SameBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_device({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_twoParts_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1}, {&part1, &part2});

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_twoParts_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1}, {&part1, &part2});

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

class NgpBatchDeclareEntitiesWithNodes : public NgpBatchDeclareEntities
{
public:
  using NodeIdGrid = std::vector<std::vector<stk::mesh::EntityId>>;

  void build_empty_mesh_without_aura(unsigned initialBucketCapacity, unsigned maximumBucketCapacity)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
  }

  void setup_hex_part()
  {
    m_hexPart = &m_meta->declare_part_with_topology("hexPart", stk::topology::HEX_8);
  }

  void declare_on_host(const std::vector<stk::mesh::EntityId>& elemIds, const NodeIdGrid& nodeIds,
                       const stk::mesh::PartVector& addParts)
  {
    Kokkos::View<stk::mesh::EntityId*, Kokkos::HostSpace> elemIdsView("elemIds", elemIds.size());
    fill_views(elemIdsView, elemIds);
    Kokkos::View<stk::mesh::EntityId**, Kokkos::HostSpace> nodeIdsView("nodeIds", nodeIds.size(), num_cols(nodeIds));
    for (size_t i = 0; i < nodeIds.size(); ++i) {
      for (size_t j = 0; j < nodeIds[i].size(); ++j) {
        nodeIdsView(i, j) = nodeIds[i][j];
      }
    }
    HostPartOrdinalsType addPartOrdinals("addPartOrdinals", addParts.size());
    fill_views(addPartOrdinals, extract_part_ordinals(addParts));

    stk::mesh::HostMesh hostMesh(*m_bulk);
    HostEntitiesType declaredEntities("declaredEntities", 0);
    Kokkos::View<stk::mesh::Entity**, stk::ngp::HostExecSpace> connectedNodes("connectedNodes", 0, 0);
    hostMesh.batch_declare_entities(stk::topology::ELEM_RANK, elemIdsView, nodeIdsView, addPartOrdinals,
                                    declaredEntities, connectedNodes);
    hostMesh.update_bulk_data();

    store_results(declaredEntities, connectedNodes);
  }

  void check_declared_elements(const std::vector<stk::mesh::EntityId>& elemIds, const NodeIdGrid& nodeIds,
                               unsigned expectedNumNodes)
  {
    ASSERT_EQ(elemIds.size(), m_declaredEntities.size());
    ASSERT_EQ(nodeIds.size(), m_connectedNodes.size());

    for (size_t i = 0; i < elemIds.size(); ++i) {
      stk::mesh::Entity elem = m_bulk->get_entity(stk::topology::ELEM_RANK, elemIds[i]);
      ASSERT_TRUE(m_bulk->is_valid(elem));
      EXPECT_EQ(elem, m_declaredEntities[i]);
      EXPECT_TRUE(m_bulk->bucket(elem).member(*m_hexPart));
      EXPECT_EQ(stk::topology::HEX_8, m_bulk->bucket(elem).topology());
      EXPECT_TRUE(m_bulk->bucket(elem).owned());

      ASSERT_EQ(nodeIds[i].size(), m_bulk->num_nodes(elem));
      const stk::mesh::Entity* elemNodes = m_bulk->begin_nodes(elem);
      const stk::mesh::ConnectivityOrdinal* elemNodeOrds = m_bulk->begin_node_ordinals(elem);
      for (size_t j = 0; j < nodeIds[i].size(); ++j) {
        EXPECT_EQ(nodeIds[i][j], m_bulk->identifier(elemNodes[j]));
        EXPECT_EQ(static_cast<stk::mesh::ConnectivityOrdinal>(j), elemNodeOrds[j]);
        EXPECT_EQ(elemNodes[j], m_connectedNodes[i][j]);

        // The node is induced into the element's topology part and knows about its new element.
        EXPECT_TRUE(m_bulk->bucket(elemNodes[j]).member(*m_hexPart));
        const stk::mesh::Entity* nodeElems = m_bulk->begin_elements(elemNodes[j]);
        const stk::mesh::Entity* nodeElemsEnd = m_bulk->end_elements(elemNodes[j]);
        EXPECT_NE(nodeElemsEnd, std::find(nodeElems, nodeElemsEnd, elem));
      }
    }

    EXPECT_EQ(expectedNumNodes, stk::mesh::count_entities(*m_bulk, stk::topology::NODE_RANK, m_meta->universal_part()));

    // The synchronized host mesh must survive a full modification cycle.
    EXPECT_NO_THROW(m_bulk->modification_begin());
    EXPECT_NO_THROW(m_bulk->modification_end());
  }

protected:
  static size_t num_cols(const NodeIdGrid& nodeIds) { return nodeIds.empty() ? 0 : nodeIds[0].size(); }

  template <typename EntitiesView, typename NodesView>
  void store_results(const EntitiesView& declaredEntities, const NodesView& connectedNodes)
  {
    m_declaredEntities.assign(declaredEntities.data(), declaredEntities.data() + declaredEntities.extent(0));
    m_connectedNodes.assign(connectedNodes.extent(0), std::vector<stk::mesh::Entity>(connectedNodes.extent(1)));
    for (size_t i = 0; i < connectedNodes.extent(0); ++i) {
      for (size_t j = 0; j < connectedNodes.extent(1); ++j) {
        m_connectedNodes[i][j] = connectedNodes(i, j);
      }
    }
  }

  stk::mesh::Part* m_hexPart = nullptr;
  std::vector<stk::mesh::Entity> m_declaredEntities;
  std::vector<std::vector<stk::mesh::Entity>> m_connectedNodes;
};

const std::vector<stk::mesh::EntityId> hex1Nodes = {1, 2, 3, 4, 5, 6, 7, 8};
const std::vector<stk::mesh::EntityId> hex2Nodes = {5, 6, 7, 8, 9, 10, 11, 12};  // shares a face with hex1

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, OneHexAllNewNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();

  declare_on_host({1}, {hex1Nodes}, {m_hexPart});
  check_declared_elements({1}, {hex1Nodes}, 8);
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, TwoHexesSharingNewNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(2, 2);
  setup_hex_part();

  declare_on_host({1, 2}, {hex1Nodes, hex2Nodes}, {m_hexPart});
  check_declared_elements({1, 2}, {hex1Nodes, hex2Nodes}, 12);
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, HexReusingExistingNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(2, 2);
  setup_hex_part();
  declare_on_host({1}, {hex1Nodes}, {m_hexPart});

  declare_on_host({2}, {hex2Nodes}, {m_hexPart});
  check_declared_elements({2}, {hex2Nodes}, 12);
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, HexReusingNodesFromGeneratedMesh_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(512, 512);
  stk::io::fill_mesh("generated:1x1x1", *m_bulk);
  setup_hex_part();

  // generated:1x1x1 is nodes 1..8 with the z=1 face 5..8; stack a second hex on it.
  declare_on_host({2}, {hex2Nodes}, {m_hexPart});
  check_declared_elements({2}, {hex2Nodes}, 12);
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, ExistingElementThrows_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();
  declare_on_host({1}, {hex1Nodes}, {m_hexPart});

  EXPECT_ANY_THROW(declare_on_host({1}, {hex2Nodes}, {m_hexPart}));
  EXPECT_FALSE(m_bulk->is_valid(m_bulk->get_entity(stk::topology::NODE_RANK, 9)));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, DuplicateElementIdsThrow_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();

  EXPECT_ANY_THROW(declare_on_host({1, 1}, {hex1Nodes, hex2Nodes}, {m_hexPart}));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, WrongNumberOfNodesThrows_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();

  EXPECT_ANY_THROW(declare_on_host({1}, {{1, 2, 3, 4}}, {m_hexPart}));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, NoTopologyPartThrows_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();
  stk::mesh::Part& plainPart = m_meta->declare_part("plainPart", stk::topology::ELEM_RANK);

  EXPECT_ANY_THROW(declare_on_host({1}, {hex1Nodes}, {&plainPart}));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, OverloadWithoutNodesRejectsElemRank_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_empty_mesh(1, 1);
  setup_hex_part();

  Kokkos::View<unsigned*, Kokkos::HostSpace> elemIds("elemIds", 1);
  elemIds(0) = 1;
  HostPartOrdinalsType addPartOrdinals("addPartOrdinals", 1);
  addPartOrdinals(0) = m_hexPart->mesh_meta_data_ordinal();
  HostEntitiesType declaredEntities("declaredEntities", 1);

  stk::mesh::HostMesh hostMesh(*m_bulk);
  EXPECT_ANY_THROW(hostMesh.batch_declare_entities(stk::topology::ELEM_RANK, elemIds, addPartOrdinals,
                                                   declaredEntities));
  EXPECT_FALSE(m_bulk->is_valid(m_bulk->get_entity(stk::topology::ELEM_RANK, 1)));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, ConnectingToSharedNodeThrows_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) GTEST_SKIP();
  build_empty_mesh(512, 512);
  stk::io::fill_mesh("generated:1x1x2", *m_bulk);  // nodes 5..8 are shared between the two ranks
  setup_hex_part();

  const stk::mesh::EntityId newElemId = 10 + stk::parallel_machine_rank(MPI_COMM_WORLD);
  const stk::mesh::EntityId firstNewNodeId = 100 + 10 * stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<stk::mesh::EntityId> nodes = {5, 6, 7, 8};
  for (stk::mesh::EntityId i = 0; i < 4; ++i) { nodes.push_back(firstNewNodeId + i); }

  EXPECT_ANY_THROW(declare_on_host({newElemId}, {nodes}, {m_hexPart}));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, ConnectingToAuraGhostedNodeThrows_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) GTEST_SKIP();
  build_empty_mesh(512, 512);
  stk::io::fill_mesh("generated:1x1x2", *m_bulk);
  setup_hex_part();

  // With aura, each rank's unshared face nodes are still ghosted to the other rank.
  const int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<stk::mesh::EntityId> nodes = (myRank == 0) ? std::vector<stk::mesh::EntityId>{1, 2, 3, 4}
                                                         : std::vector<stk::mesh::EntityId>{9, 10, 11, 12};
  for (stk::mesh::EntityId i = 0; i < 4; ++i) { nodes.push_back(100 + 10 * myRank + i); }

  EXPECT_ANY_THROW(declare_on_host({stk::mesh::EntityId(10 + myRank)}, {nodes}, {m_hexPart}));
}

NGP_TEST_F(NgpBatchDeclareEntitiesWithNodes, LocalHexOnEachRank_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) GTEST_SKIP();
  build_empty_mesh_without_aura(512, 512);
  stk::io::fill_mesh("generated:1x1x2", *m_bulk);
  setup_hex_part();

  // Each rank stacks a hex on the face of its own element that is not shared, using only locally
  // owned, unshared nodes and new node ids no other rank uses.
  const int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const std::vector<stk::mesh::EntityId> ownFace = (myRank == 0) ? std::vector<stk::mesh::EntityId>{1, 2, 3, 4}
                                                                 : std::vector<stk::mesh::EntityId>{9, 10, 11, 12};
  std::vector<stk::mesh::EntityId> nodes = ownFace;
  for (stk::mesh::EntityId i = 0; i < 4; ++i) { nodes.push_back(100 + 10 * myRank + i); }
  const stk::mesh::EntityId newElemId = 10 + myRank;
  const unsigned numNodesBefore = stk::mesh::count_entities(*m_bulk, stk::topology::NODE_RANK, m_meta->universal_part());

  declare_on_host({newElemId}, {nodes}, {m_hexPart});
  check_declared_elements({newElemId}, {nodes}, numNodesBefore + 4);
}

NGP_TEST_F(NgpBatchDeclareNodes, declaredNodesAreFoundByIdOnHost_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  declare_on_device({nodeId1, nodeId2}, {&part1});

  for (unsigned nodeId : {nodeId1, nodeId2}) {
    const stk::mesh::Entity node = m_bulk->get_entity(stk::topology::NODE_RANK, nodeId);
    ASSERT_TRUE(m_bulk->is_valid(node));
    EXPECT_EQ(nodeId, m_bulk->identifier(node));
  }
}

NGP_TEST_F(NgpBatchDeclareNodes, declaredNodesFitInMeshConnectivity_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(512, 512);
  stk::io::fill_mesh("generated:1x1x1", *m_bulk);

  // New entities get local offsets past the end of the device entityKeys, which is sized from the host
  // entity-key capacity, not its size.  The connectivity must cover those offsets, or nothing can be
  // connected to the new entities.
  const std::vector<stk::mesh::Entity> newNodes =
      batch_declare_entities_on_device_for_test(*m_bulk, stk::topology::NODE_RANK, {100, 101, 102}, {});

  const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_ngp_mesh(*m_bulk);
  for (const stk::mesh::Entity node : newNodes) {
    EXPECT_LT(node.local_offset(), ngpMesh.get_mesh_connectivity().get_num_entities());
  }
}

}  // namespace

std::vector<stk::mesh::Entity> batch_declare_entities_on_device_for_test(
    stk::mesh::BulkData& bulk, stk::topology::rank_t rank, const std::vector<unsigned>& entityIds,
    const std::vector<stk::mesh::PartOrdinal>& addPartOrdinals)
{
  Kokkos::View<unsigned*> entityIdsView("entityIds", entityIds.size());
  fill_views(entityIdsView, entityIds);

  DevicePartOrdinalsType addPartOrdinalsView("addPartOrdinals", addPartOrdinals.size());
  fill_views(addPartOrdinalsView, addPartOrdinals);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  DeviceEntitiesType createdEntities("createdEntities", entityIds.size());
  ngpMesh.batch_declare_entities(rank, entityIdsView, addPartOrdinalsView, createdEntities);

  auto hostCreatedEntities = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), createdEntities);
  return std::vector<stk::mesh::Entity>(hostCreatedEntities.data(), hostCreatedEntities.data() + entityIds.size());
}

#endif
