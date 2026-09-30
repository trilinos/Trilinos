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
TEST_F(NgpBatchDestroyConnectivities, oneElem_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), 1);
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  hostMesh.batch_destroy_relations(hexesView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

TEST_F(NgpBatchDestroyConnectivities, TwoElems_DeleteOne_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), 1);
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  for (unsigned i = 0u; i < 5u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (unsigned i = 5u; i < quads.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  hostMesh.batch_destroy_relations(hexesView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();

  for (unsigned i = 0u; i < 5u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (unsigned i = 5u; i < quads.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (unsigned i = 0u; i < 4u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 4u; i < 8u; ++i) {
    EXPECT_TRUE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 8u; i < 12u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 12u; i < edges.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

TEST_F(NgpBatchDestroyConnectivities, TwoElems_DeleteBoth_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), 2);
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  hostMesh.batch_destroy_relations(hexesView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_does_not_have_parts(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank,
                                             EntitiesViewType const& entities,
                                             DevicePartOrdinalsType const& partOrdinals)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIndex.bucket_id);

      for (unsigned i = 0; i < partOrdinals.extent(0); ++i) {
        NGP_EXPECT_FALSE(bucket.member(partOrdinals(i)));
      }
    }
  );
  Kokkos::fence();
}

class NgpBatchDestroyRelations : public NgpBatchDestroyConnectivities
{
public:
  NgpBatchDestroyRelations() {}

  stk::mesh::Entity find_face_with_num_elements(unsigned numElems)
  {
    for (stk::mesh::Entity face : quads) {
      if (m_bulk->num_connectivity(face, stk::topology::ELEM_RANK) == numElems) {
        return face;
      }
    }
    return stk::mesh::Entity();
  }

  DeviceEntitiesType make_device_entities(const std::vector<stk::mesh::Entity>& ents)
  {
    DeviceEntitiesType deviceEnts("deviceEntities", ents.size());
    auto hostEnts = Kokkos::create_mirror_view(deviceEnts);
    for (size_t i = 0; i < ents.size(); ++i) {
      hostEnts(i) = ents[i];
    }
    Kokkos::deep_copy(deviceEnts, hostEnts);
    return deviceEnts;
  }

  HostEntitiesType make_host_entities(const std::vector<stk::mesh::Entity>& ents)
  {
    HostEntitiesType hostEnts("hostEntities", ents.size());
    for (size_t i = 0; i < ents.size(); ++i) {
      hostEnts(i) = ents[i];
    }
    return hostEnts;
  }
};

NGP_TEST_F(NgpBatchDestroyRelations, basicDestroy_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face = quads[0];
  EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));

  stk::mesh::HostMesh hostMesh(*m_bulk);

  auto hexView = make_host_entities({hex});
  auto faceView = make_host_entities({face});
  check_host_connectivity(hostMesh, stk::topology::ELEM_RANK, hexView, stk::topology::FACE_RANK, 6u);
  check_host_connectivity(hostMesh, stk::topology::FACE_RANK, faceView, stk::topology::ELEM_RANK, 1u);

  hostMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);

  check_host_connectivity(hostMesh, stk::topology::ELEM_RANK, hexView, stk::topology::FACE_RANK, 5u);
  check_host_connectivity(hostMesh, stk::topology::FACE_RANK, faceView, stk::topology::ELEM_RANK, 0u);

  hostMesh.update_bulk_data();
  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
}

NGP_TEST_F(NgpBatchDestroyRelations, basicDestroy_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face = quads[0];
  EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto hexView = make_device_entities({hex});
  auto faceView = make_device_entities({face});
  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, hexView, stk::topology::FACE_RANK, 6u);
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, faceView, stk::topology::ELEM_RANK, 1u);

  ngpMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);

  check_device_connectivity(ngpMesh, stk::topology::ELEM_RANK, hexView, stk::topology::FACE_RANK, 5u);
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, faceView, stk::topology::ELEM_RANK, 0u);

  ngpMesh.update_bulk_data();
  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
}

NGP_TEST_F(NgpBatchDestroyRelations, inducedPartRetained_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex1 = m_bulk->get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity sharedFace = find_face_with_num_elements(2u);
  ASSERT_TRUE(m_bulk->is_valid(sharedFace));
  EXPECT_TRUE(m_bulk->bucket(sharedFace).member(hexPart));

  stk::mesh::HostMesh hostMesh(*m_bulk);

  auto sharedFaceView = make_host_entities({sharedFace});
  check_host_connectivity(hostMesh, stk::topology::FACE_RANK, sharedFaceView, stk::topology::ELEM_RANK, 2u);

  hostMesh.batch_destroy_relations(make_host_entities({hex1}), stk::topology::FACE_RANK);

  check_host_connectivity(hostMesh, stk::topology::FACE_RANK, sharedFaceView, stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(m_bulk->bucket(sharedFace).member(hexPart));

  hostMesh.update_bulk_data();
  EXPECT_TRUE(m_bulk->bucket(sharedFace).member(hexPart));
  EXPECT_EQ(1u, m_bulk->num_connectivity(sharedFace, stk::topology::ELEM_RANK));
}

NGP_TEST_F(NgpBatchDestroyRelations, inducedPartRetained_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex1 = m_bulk->get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity sharedFace = find_face_with_num_elements(2u);  // middle face, shared by both hexes
  ASSERT_TRUE(m_bulk->is_valid(sharedFace));
  EXPECT_TRUE(m_bulk->bucket(sharedFace).member(hexPart));

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto sharedFaceView = make_device_entities({sharedFace});
  auto hexPartView = create_device_part_ordinal({&hexPart});
  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, sharedFaceView, stk::topology::ELEM_RANK, 2u);

  ngpMesh.batch_destroy_relations(make_device_entities({hex1}), stk::topology::FACE_RANK);

  check_device_connectivity(ngpMesh, stk::topology::FACE_RANK, sharedFaceView, stk::topology::ELEM_RANK, 1u);
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, sharedFaceView, hexPartView);

  ngpMesh.update_bulk_data();
  EXPECT_TRUE(m_bulk->bucket(sharedFace).member(hexPart));
  EXPECT_EQ(1u, m_bulk->num_connectivity(sharedFace, stk::topology::ELEM_RANK));
}

NGP_TEST_F(NgpBatchDestroyRelations, batchMixedWithInvalid_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face0 = quads[0];

  stk::mesh::HostMesh hostMesh(*m_bulk);

  std::vector<stk::mesh::Entity> ents{hex, stk::mesh::Entity(), face0};
  EXPECT_ANY_THROW(hostMesh.batch_destroy_relations(make_host_entities(ents), stk::topology::FACE_RANK));
}

NGP_TEST_F(NgpBatchDestroyRelations, batchMixedWithInvalid_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face0 = quads[0];

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  std::vector<stk::mesh::Entity> ents{hex, stk::mesh::Entity(), face0};
  EXPECT_ANY_THROW(ngpMesh.batch_destroy_relations(make_device_entities(ents), stk::topology::FACE_RANK));
}

NGP_TEST_F(NgpBatchDestroyRelations, batchMixedRanks_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face0 = quads[0];

  stk::mesh::HostMesh hostMesh(*m_bulk);

  std::vector<stk::mesh::Entity> ents{hex, face0};
  hostMesh.batch_destroy_relations(make_host_entities(ents), stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDestroyRelations, batchMixedRanks_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity face0 = quads[0];

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto face0View = make_device_entities({face0});
  auto hexPartView = create_device_part_ordinal({&hexPart});

  std::vector<stk::mesh::Entity> ents{hex, face0};
  ngpMesh.batch_destroy_relations(make_device_entities(ents), stk::topology::FACE_RANK);

  check_device_entity_does_not_have_parts(ngpMesh, stk::topology::FACE_RANK, face0View, hexPartView);

  ngpMesh.update_bulk_data();
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDestroyRelations, forceNoInducePartRetained_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& noInducePart = m_meta->declare_part("noInduceElemPart", stk::topology::ELEM_RANK,
                                                       /*forceNoInduce=*/true);

  stk::mesh::Entity face = quads[0];

  m_bulk->modification_begin();
  m_bulk->change_entity_parts(stk::mesh::EntityVector{face}, stk::mesh::PartVector{&noInducePart});
  m_bulk->modification_end();

  ASSERT_TRUE(m_bulk->bucket(face).member(hexPart));
  ASSERT_TRUE(m_bulk->bucket(face).member(noInducePart));

  stk::mesh::HostMesh hostMesh(*m_bulk);
  auto faceView = make_host_entities({face});

  hostMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);
  hostMesh.update_bulk_data();

  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(face).member(noInducePart));
}

NGP_TEST_F(NgpBatchDestroyRelations, forceNoInducePartRetained_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& noInducePart = m_meta->declare_part("noInduceElemPart", stk::topology::ELEM_RANK,
                                                       /*forceNoInduce=*/true);

  stk::mesh::Entity face = quads[0];

  m_bulk->modification_begin();
  m_bulk->change_entity_parts(stk::mesh::EntityVector{face}, stk::mesh::PartVector{&noInducePart});
  m_bulk->modification_end();

  ASSERT_TRUE(m_bulk->bucket(face).member(hexPart));
  ASSERT_TRUE(m_bulk->bucket(face).member(noInducePart));

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto faceView = make_device_entities({face});
  auto hexPartView = create_device_part_ordinal({&hexPart});
  auto noInducePartView = create_device_part_ordinal({&noInducePart});

  ngpMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);

  check_device_entity_does_not_have_parts(ngpMesh, stk::topology::FACE_RANK, faceView, hexPartView);
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, faceView, noInducePartView);

  ngpMesh.update_bulk_data();
  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(face).member(noInducePart));
}

NGP_TEST_F(NgpBatchDestroyRelations, explicitInduciblePartRetained_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& otherElemPart = m_meta->declare_part("otherElemPart", stk::topology::ELEM_RANK);

  stk::mesh::Entity face = quads[0];

  m_bulk->modification_begin();
  m_bulk->change_entity_parts(stk::mesh::EntityVector{face}, stk::mesh::PartVector{&otherElemPart});
  m_bulk->modification_end();

  ASSERT_TRUE(m_bulk->bucket(face).member(hexPart));
  ASSERT_TRUE(m_bulk->bucket(face).member(otherElemPart));
  ASSERT_FALSE(m_bulk->bucket(hexes[0]).member(otherElemPart));

  stk::mesh::HostMesh hostMesh(*m_bulk);
  auto faceView = make_host_entities({face});

  hostMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);
  hostMesh.update_bulk_data();

  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(face).member(otherElemPart));
}

NGP_TEST_F(NgpBatchDestroyRelations, explicitInduciblePartRetained_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& otherElemPart = m_meta->declare_part("otherElemPart", stk::topology::ELEM_RANK);

  stk::mesh::Entity face = quads[0];

  m_bulk->modification_begin();
  m_bulk->change_entity_parts(stk::mesh::EntityVector{face}, stk::mesh::PartVector{&otherElemPart});
  m_bulk->modification_end();

  ASSERT_TRUE(m_bulk->bucket(face).member(hexPart));
  ASSERT_TRUE(m_bulk->bucket(face).member(otherElemPart));
  ASSERT_FALSE(m_bulk->bucket(hexes[0]).member(otherElemPart));

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  auto faceView = make_device_entities({face});
  auto hexPartView = create_device_part_ordinal({&hexPart});
  auto otherElemPartView = create_device_part_ordinal({&otherElemPart});

  ngpMesh.batch_destroy_relations(faceView, stk::topology::ELEM_RANK);

  check_device_entity_does_not_have_parts(ngpMesh, stk::topology::FACE_RANK, faceView, hexPartView);
  check_device_entity_has_parts(ngpMesh, stk::topology::FACE_RANK, faceView, otherElemPartView);

  ngpMesh.update_bulk_data();
  EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(face).member(otherElemPart));
}

NGP_TEST_F(NgpBatchDestroyRelations, oneElem_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0]}), stk::topology::FACE_RANK);
  ngpMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0]}), stk::topology::EDGE_RANK);
  ngpMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

NGP_TEST_F(NgpBatchDestroyRelations, TwoElems_DeleteOne_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0]}), stk::topology::FACE_RANK);
  ngpMesh.update_bulk_data();

  for (unsigned i = 0u; i < 5u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (unsigned i = 5u; i < quads.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0]}), stk::topology::EDGE_RANK);
  ngpMesh.update_bulk_data();

  for (unsigned i = 0u; i < 5u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (unsigned i = 5u; i < quads.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(quads[i]).member(hexPart));
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (unsigned i = 0u; i < 4u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 4u; i < 8u; ++i) {
    EXPECT_TRUE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 8u; i < 12u; ++i) {
    EXPECT_FALSE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (unsigned i = 12u; i < edges.size(); ++i) {
    EXPECT_TRUE(m_bulk->bucket(edges[i]).member(hexPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

NGP_TEST_F(NgpBatchDestroyRelations, TwoElems_DeleteBoth_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& nodePart = *m_meta->get_part("nodePart");
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node : nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0], hexes[1]}), stk::topology::FACE_RANK);
  ngpMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }

  ngpMesh.batch_destroy_relations(make_device_entities({hexes[0], hexes[1]}), stk::topology::EDGE_RANK);
  ngpMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
  for (auto& node: nodes) {
    EXPECT_TRUE(m_bulk->bucket(node).member(quadPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(edgePart));
    EXPECT_TRUE(m_bulk->bucket(node).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(node).member(nodePart));
  }
}

}  // namespace
#endif
