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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_tools/mesh_tools/EdgeMidNodeDetector.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "DisconnectBlocksMeshConstruction.hpp"
#include "stk_mesh/base/Types.hpp"

#include <string>
#include <algorithm>
#include <map>
#include <utility>

using stk::unit_test_util::build_mesh;

namespace {

std::string get_text_mesh_string(stk::topology stkTopology, const std::string& textMeshTopologyString)
{
  std::string meshDesc = "0,1," + textMeshTopologyString + ",";
  for(unsigned i=1; i<=stkTopology.num_nodes(); i++) {
    meshDesc += std::to_string(i);
    meshDesc += ",";
  }
  meshDesc += "block_1";
  return meshDesc;
}

stk::mesh::PartVector setup_mesh_with_one_block_and_one_element(stk::mesh::BulkData& bulk,
                                                                stk::topology::topology_t stkTopology,
                                                                const std::string& textMeshTopologyString)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stkTopology, "block_1", 1);
  std::string meshDesc = get_text_mesh_string(stk::topology(stkTopology), textMeshTopologyString);
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
  return {&block1};
}

void test_edge_mid_node_entities(const stk::mesh::BulkData& bulk,
                                 const stk::mesh::EntityVector& edgeMidNodes,
                                 const stk::mesh::EntityVector& goldNodes)
{
  EXPECT_EQ(edgeMidNodes.size(), goldNodes.size()) << "Mid nodes size: " << edgeMidNodes.size() << " != " << goldNodes.size();

  for(auto node : edgeMidNodes) {
    auto it = std::find_if(goldNodes.begin(), goldNodes.end(), [node](stk::mesh::Entity n) {
                                                                  return n == node;
                                                                });

    EXPECT_TRUE (it != goldNodes.end()) << "Node " << bulk.identifier(node) << " not found in gold list";
  }
}

void test_edge_mid_node_ids(const stk::mesh::BulkData& bulk,
                            const stk::mesh::EntityVector& edgeMidNodes,
                            const stk::mesh::EntityIdVector& goldNodeIds)
{
  stk::mesh::EntityVector goldNodes;
  goldNodes.reserve(goldNodeIds.size());

  for(stk::mesh::EntityId id : goldNodeIds) {
    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, id);

    EXPECT_TRUE(bulk.is_valid(node)) << "Invalid node with id: " << id;
    goldNodes.push_back(node);
  }

  test_edge_mid_node_entities(bulk, edgeMidNodes, goldNodes);
}

stk::mesh::EntityVector get_edge_mid_nodes(const stk::mesh::BulkData& bulk)
{
  stk::tools::EdgeMidNodeDetector midNodeDetector(bulk);
  return midNodeDetector.get_all_edge_mid_nodes();
}

}


TEST(TestEdgeMidNode, EmptyMesh)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Particle)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::PARTICLE, "PARTICLE");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Line2_1D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(1,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::LINE_2_1D, "LINE_2_1D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Line3_1D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(1,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::LINE_3_1D, "LINE_3_1D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{3u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Beam2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::BEAM_2, "BEAM_2");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Beam3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::BEAM_3, "BEAM_3");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{3u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Spring2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SPRING_2, "SPRING_2");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Spring3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SPRING_3, "SPRING_3");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{3u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellLine2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_LINE_2, "SHELL_LINE_2");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellLine3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_LINE_3, "SHELL_LINE_3");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{3u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tri3_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TRI_3_2D, "TRI_3_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tri4_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TRI_4_2D, "TRI_4_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tri6_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TRI_6_2D, "TRI_6_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{4u, 5u, 6u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellTri3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_TRI_3, "SHELL_TRI_3");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, DISABLED_ShellTri4)
{
  // shards does not define a topology for a 4-noded triangular shell
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_TRI_4, "SHELL_TRI_4");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellTri6)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_TRI_6, "SHELL_TRI_6");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{4u, 5u, 6u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Quad4_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::QUAD_4_2D, "QUAD_4_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Quad8_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::QUAD_8_2D, "QUAD_8_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Quad9_2D)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(2,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::QUAD_9_2D, "QUAD_9_2D");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellQuad4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_QUAD_4, "SHELL_QUAD_4");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellQuad8)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_QUAD_8, "SHELL_QUAD_8");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, ShellQuad9)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::SHELL_QUAD_9, "SHELL_QUAD_9");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tet4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TET_4, "TET_4");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tet8)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TET_8, "TET_8");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tet10)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TET_10, "TET_10");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u, 9u, 10u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Tet11)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::TET_11, "TET_11");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u, 9u, 10u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Pyramid5)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::PYRAMID_5, "PYRAMID_5");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Pyramid13)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::PYRAMID_13, "PYRAMID_13");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{6u, 7u, 8u, 9u, 10u, 11u, 12u, 13u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Pyramid14)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::PYRAMID_14, "PYRAMID_14");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{6u, 7u, 8u, 9u, 10u, 11u, 12u, 13u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Wedge6)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::WEDGE_6, "WEDGE_6");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Wedge12)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::WEDGE_12, "WEDGE_12");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{7u, 8u, 9u, 10u, 11u, 12u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Wedge15)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::WEDGE_15, "WEDGE_15");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Wedge18)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::WEDGE_18, "WEDGE_18");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Hex8)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::HEX_8, "HEX_8");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Hex20)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::HEX_20, "HEX_20");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{9u, 10u, 11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u, 20u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, Hex27)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);
  setup_mesh_with_one_block_and_one_element(*bulk, stk::topology::HEX_27, "HEX_27");
  stk::mesh::EntityVector edgeMidNodes = ::get_edge_mid_nodes(*bulk);
  stk::mesh::EntityIdVector goldNodeIds{9u, 10u, 11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u, 20u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);
}

TEST(TestEdgeMidNode, verifyEdgeCount)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);

  create_part(bulk->mesh_meta_data(), stk::topology::TET_10, "block_1", 1);
  create_part(bulk->mesh_meta_data(), stk::topology::TET_4,  "block_2", 2);

  std::string meshDesc = "0,1,TET_10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,block_1\n"
                         "0,2,TET_4, 11,12,13,14,block_2";

  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  stk::tools::EdgeMidNodeDetector edges(*bulk);
  EXPECT_EQ( 6u, edges.num_mid_nodes());
  EXPECT_EQ(12u, edges.num_edges());
}

TEST(TestEdgeMidNode, verifyTwoTet10MidNodeConnectivity)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3,MPI_COMM_WORLD);

  create_part(bulk->mesh_meta_data(), stk::topology::TET_10, "block_1", 1);

  std::string meshDesc = "0,1,TET_10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,block_1\n"
                         "0,2,TET_10,11,12,13,14,15,16,17,18,19, 5,block_1";

  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  stk::tools::EdgeMidNodeDetector midNodeDetector(*bulk);
  stk::mesh::EntityVector edgeMidNodes = midNodeDetector.get_all_edge_mid_nodes();
  stk::mesh::EntityIdVector goldNodeIds{5u, 6u, 7u, 8u, 9u, 10u, 15u, 16u, 17u, 18u, 19u};
  test_edge_mid_node_ids(*bulk, edgeMidNodes, goldNodeIds);

  EXPECT_FALSE(midNodeDetector.is_mid_edge_node(1u));
  EXPECT_TRUE (midNodeDetector.is_mid_edge_node(5u));

  stk::mesh::Entity node1  = bulk->get_entity(stk::topology::NODE_RANK,  1u);
  stk::mesh::Entity node2  = bulk->get_entity(stk::topology::NODE_RANK,  2u);
  stk::mesh::Entity node5  = bulk->get_entity(stk::topology::NODE_RANK,  5u);
  stk::mesh::Entity node13 = bulk->get_entity(stk::topology::NODE_RANK, 13u);
  stk::mesh::Entity node14 = bulk->get_entity(stk::topology::NODE_RANK, 14u);

  EXPECT_TRUE(bulk->is_valid(node1));
  EXPECT_TRUE(bulk->is_valid(node2));
  EXPECT_TRUE(bulk->is_valid(node5));
  EXPECT_TRUE(bulk->is_valid(node13));
  EXPECT_TRUE(bulk->is_valid(node14));

  const std::vector<stk::tools::EdgeVertices>& edges = midNodeDetector.get_edge_vertices(node5);
  EXPECT_EQ(2u, edges.size());

  std::vector<stk::tools::EdgeVertices> goldEdges{{node1,node2}, {node13,node14}};
  for(const auto& edge : edges) {
    auto it = std::find_if(goldEdges.begin(), goldEdges.end(), [edge](const stk::tools::EdgeVertices& goldEdge) {
                                                                  return goldEdge == edge;
                                                                });

    EXPECT_TRUE (it != goldEdges.end()) << "Edge {" << bulk->identifier(edge.first)
                                        << " , "    << bulk->identifier(edge.second) << "} not found in gold list";
  }
}

