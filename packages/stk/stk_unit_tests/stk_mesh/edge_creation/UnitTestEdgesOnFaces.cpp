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

#include "gtest/gtest.h"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"

namespace {

class EdgesOnFacesFixture : public ::testing::Test
{
  public:
    void setup(const std::string& meshSpec, bool enable_aura=true)
    {
      stk::mesh::MeshBuilder builder(stk::parallel_machine_world());
      builder.set_aura_option(enable_aura ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA);
      bulk = builder.create();

      stk::io::fill_mesh(meshSpec, *bulk);

      stk::mesh::get_entities(*bulk, stk::topology::FACE_RANK, bulk->mesh_meta_data().locally_owned_part(),
                              ownedFaces);

      surfaceSelector = stk::mesh::Selector(*(bulk->get_sidesets()[0]->get_part()));
    }

    std::shared_ptr<stk::mesh::BulkData> bulk;
    stk::mesh::EntityVector ownedFaces;
    stk::mesh::Selector surfaceSelector;
};

std::pair<int, int> countEdges(stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
{
  int edgesInSelector = stk::mesh::count_entities(bulk, stk::topology::EDGE_RANK,
                                                  selector & bulk.mesh_meta_data().locally_owned_part());

  int totalEdges = stk::mesh::count_entities(bulk, stk::topology::EDGE_RANK,
                                             bulk.mesh_meta_data().locally_owned_part());

  std::array<int, 2> edgeCountsLocal = {edgesInSelector, totalEdges}, edgeCountsGlobal;
  stk::all_reduce_sum(stk::parallel_machine_world(), edgeCountsLocal.data(), edgeCountsGlobal.data(), 2);

  return std::make_pair(edgeCountsGlobal[0], edgeCountsGlobal[1]);
}

void test_faces_have_edges(stk::mesh::BulkData& bulk, const std::vector<stk::mesh::Entity>& faces,
                           int numExpectedEdgesPerFace)
{
  for (stk::mesh::Entity face : faces)
  {
    stk::mesh::Entity const* beginEdges = bulk.begin_edges(face);
    stk::mesh::Entity const* endEdges = bulk.end_edges(face);
    ptrdiff_t numEdges = endEdges - beginEdges;

    EXPECT_EQ(numEdges, numExpectedEdgesPerFace);
    for (int i=0; i < numEdges; ++i)
    {
      for (int j=0; j < numEdges; ++j)
      {
        if (i != j)
        {
          EXPECT_NE(beginEdges[i], beginEdges[j]);
        }
      }
    }
  }
}

void test_third_node_is_midpoint(stk::mesh::BulkData& bulk)
{
  std::vector<stk::mesh::Entity> edgeNodes;
  std::vector<std::array<double, 3>> edgeNodeCoords;
  const stk::mesh::FieldBase* coordField = bulk.mesh_meta_data().coordinate_field();

  for (stk::mesh::Bucket* bucket : bulk.buckets(stk::topology::EDGE_RANK))
  {
    for (stk::mesh::Entity edge : *bucket)
    {
      edgeNodes.clear();
      for (stk::mesh::Entity const* node = bulk.begin_nodes(edge); node != bulk.end_nodes(edge); ++node)
      {
        edgeNodes.push_back(*node);
        double* coordData = static_cast<double*>(stk::mesh::field_data(*coordField, *node));
        edgeNodeCoords.emplace_back(std::array<double, 3>{coordData[0], coordData[1], coordData[2]});
      }

      EXPECT_EQ(edgeNodes.size(), 3U);
      for (int i=0; i < 3; ++i)
      {
        EXPECT_FLOAT_EQ(edgeNodeCoords[2][i], 0.5*(edgeNodeCoords[0][i] + edgeNodeCoords[1][i]));
      }
    }
  }

}

}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshNoAura)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", false);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 24);
  EXPECT_EQ(edgesTotal, 24);
  test_faces_have_edges(*bulk, ownedFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshWithAura)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", true);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 24);
  EXPECT_EQ(edgesTotal, 24);
  test_faces_have_edges(*bulk, ownedFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex4x4MeshWithAuraAndCustomGhosting)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  setup("generated:4x4x4|sideset:x", false);  // TODO: turn aura back on

  bulk->modification_begin();
  stk::mesh::Ghosting& ghosting = bulk->create_ghosting("my_custom_ghosting");
  bulk->modification_end();

  bulk->modification_begin();
  int otherProc = 1 - stk::parallel_machine_rank(stk::parallel_machine_world());
  std::vector<stk::mesh::EntityProc> sendGhosts;
  for (stk::mesh::Entity face : ownedFaces)
  {
    sendGhosts.emplace_back(face, otherProc);
  }
  bulk->change_ghosting(ghosting, sendGhosts);
  bulk->modification_end();

  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  // its unfortunate we have to do this, because stk does have enough information to update
  // the downward adjacencies of the ghost faces
  bulk->modification_begin();
  sendGhosts.clear();
  for (stk::mesh::Bucket* bucket : bulk->get_buckets(stk::topology::EDGE_RANK, *edgePart & bulk->mesh_meta_data().locally_owned_part()))
  {
    for (stk::mesh::Entity edge : *bucket)
    {
      sendGhosts.emplace_back(edge, otherProc);
    }
  }
  bulk->change_ghosting(ghosting, sendGhosts);
  bulk->modification_end();


  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 40);
  EXPECT_EQ(edgesTotal, 40);

  std::vector<stk::mesh::Entity> allFaces;
  stk::mesh::get_selected_entities(surfaceSelector, bulk->buckets(stk::topology::FACE_RANK), allFaces);
  EXPECT_EQ(allFaces.size(), 4*4U);
  test_faces_have_edges(*bulk, allFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshEnableAuraLater)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", false);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  bulk->set_automatic_aura_option(stk::mesh::BulkData::AUTO_AURA, true);

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 24);
  EXPECT_EQ(edgesTotal, 24);
  test_faces_have_edges(*bulk, ownedFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex6x6MeshWithAura)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 6)
  {
    GTEST_SKIP();
  }

  setup("generated:6x6x6|sideset:x", true);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 6*7*2);
  EXPECT_EQ(edgesTotal, 6*7*2);
  test_faces_have_edges(*bulk, ownedFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshProvidePart)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", true);
  stk::mesh::Part* edgePart = &(bulk->mesh_meta_data().declare_part_with_topology("face_edges_part", stk::topology::LINE_2));

  bulk->modification_begin();
  stk::mesh::Part* edgePart2 = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector, edgePart);
  bulk->modification_end();
  EXPECT_EQ(edgePart, edgePart2);

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 24);
  EXPECT_EQ(edgesTotal, 24);
  test_faces_have_edges(*bulk, ownedFaces, 4);
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshProvideIncorrectPart)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", true);
  stk::mesh::Part* edgePart = &(bulk->mesh_meta_data().declare_part_with_topology("face_edges_part", stk::topology::QUAD_4));

  EXPECT_ANY_THROW(stk::mesh::create_edges_on_faces(*bulk, surfaceSelector, edgePart));
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshNoAuraDeleteEdges)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", false);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  bulk->modification_begin();
  stk::mesh::delete_edges_on_faces(*bulk, edgePart);
  bulk->modification_end();


  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 0);
  EXPECT_EQ(edgesTotal, 0);
  test_faces_have_edges(*bulk, ownedFaces, 0);
}

TEST_F(EdgesOnFacesFixture, Hex3x3MeshWithAuraDeleteEdges)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", true);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  bulk->modification_begin();
  stk::mesh::delete_edges_on_faces(*bulk, edgePart);
  bulk->modification_end();


  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 0);
  EXPECT_EQ(edgesTotal, 0);
  test_faces_have_edges(*bulk, ownedFaces, 0);
}

TEST_F(EdgesOnFacesFixture, Hex2x1HighOrder)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 2)
  {
    GTEST_SKIP();
  }
  std::string meshDesc = "textmesh:"
                         "0,1,HEX_20,1,3,11,9,  20,22,30,28, 2,7,10,6,    14,15,18,17, 21,26,29,25\n"
                         "0,2,HEX_20,3,5,13,11, 22,24,32,30, 4,8,12,7,    15,16,19,18, 23,27,31,26\n"
                                       // bottom plane
                         "|coordinates: 0,0,0,  1,0,0,  2,0,0,  3,0,0,  4,0,0,"
                                       "0,1,0,  2,1,0,  5,1,0,"
                                       "0,2,0,  1,2,0,  2,2,0,  3,2,0,  4,2,0,"
                                       // middle plane
                                       "0,0,1,  3,0,1,  5,0,1,"
                                       "0,2,1,  3,2,1,  5,2,1,"
                                       // top plane
                                       "0,0,1,  1,0,1,  2,0,1,  3,0,1,  4,0,1,"
                                       "0,1,1,  2,1,1,  5,1,1,"
                                       "0,2,1,  1,2,1,  2,2,1,  3,2,1,  4,2,1\n"
                         "|sideset:name=surface_1; data=1,1, 2,1";
  setup(meshDesc, true);

  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, surfaceSelector);
  bulk->modification_end();

  const auto& [edgesInSelector, edgesTotal] = countEdges(*bulk, *edgePart);
  EXPECT_EQ(edgesInSelector, 7);
  EXPECT_EQ(edgesTotal, 7);
  test_third_node_is_midpoint(*bulk);
}

TEST_F(EdgesOnFacesFixture, NoFacesSelected)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  setup("generated:3x3x3|sideset:x", false);
  bulk->modification_begin();
  stk::mesh::Part* edgePart = stk::mesh::create_edges_on_faces(*bulk, stk::mesh::Selector());
  bulk->modification_end();
  EXPECT_EQ(edgePart, nullptr);
}