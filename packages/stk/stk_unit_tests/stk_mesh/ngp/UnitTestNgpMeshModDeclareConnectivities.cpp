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
class NgpBatchDeclareRelations : public NgpBatchDestroyConnectivities
{
public:
  using HostOffsetsType = Kokkos::View<unsigned*, stk::ngp::HostExecSpace>;
  using HostOrdinalsType = Kokkos::View<stk::mesh::RelationIdentifier*, stk::ngp::HostExecSpace>;
  using HostPermutationsType = Kokkos::View<stk::mesh::Permutation*, stk::ngp::HostExecSpace>;
  using Host2dEntitiesType = Kokkos::View<stk::mesh::Entity**, stk::ngp::HostExecSpace>;
  using Host2dPermutationsType = Kokkos::View<stk::mesh::Permutation**, stk::ngp::HostExecSpace>;

  using DeviceOffsetsType = Kokkos::View<unsigned*, stk::ngp::MemSpace>;
  using DeviceOrdinalsType = Kokkos::View<stk::mesh::RelationIdentifier*, stk::ngp::MemSpace>;
  using DevicePermutationsType = Kokkos::View<stk::mesh::Permutation*, stk::ngp::MemSpace>;
  using Device2dEntitiesType = Kokkos::View<stk::mesh::Entity**, stk::ngp::MemSpace>;
  using Device2dPermutationsType = Kokkos::View<stk::mesh::Permutation**, stk::ngp::MemSpace>;

  NgpBatchDeclareRelations() {};

  struct Connectivity {
    std::vector<stk::mesh::Entity> toEntities;
    std::vector<stk::mesh::RelationIdentifier> ordinals;
    std::vector<stk::mesh::Permutation> permutations;

    void append(const Connectivity& c) {
      toEntities.insert(toEntities.end(), c.toEntities.begin(), c.toEntities.end());
      ordinals.insert(ordinals.end(), c.ordinals.begin(), c.ordinals.end());
      permutations.insert(permutations.end(), c.permutations.begin(), c.permutations.end());
    }

    Connectivity reversed_copy() const {
      Connectivity r;
      for (size_t i = toEntities.size(); i > 0; --i) {
        r.toEntities.push_back(toEntities[i-1]);
        r.ordinals.push_back(ordinals[i-1]);
        r.permutations.push_back(permutations[i-1]);
      }
      return r;
    }
  };

  Connectivity capture_connectivity(stk::mesh::Entity from, stk::mesh::EntityRank connectedRank) {
    Connectivity conn;
    const unsigned num = m_bulk->num_connectivity(from, connectedRank);
    const stk::mesh::Entity* e = m_bulk->begin(from, connectedRank);
    const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_ordinals(from, connectedRank);
    const stk::mesh::Permutation* p = m_bulk->begin_permutations(from, connectedRank);
    for (unsigned i = 0; i < num; ++i) {
      conn.toEntities.push_back(e[i]);
      conn.ordinals.push_back(static_cast<stk::mesh::RelationIdentifier>(o[i]));
      conn.permutations.push_back(p != nullptr ? p[i] : stk::mesh::Permutation::INVALID_PERMUTATION);
    }
    return conn;
  }

  void verify_slice(const stk::mesh::Entity* e, const stk::mesh::ConnectivityOrdinal* o,
                    const stk::mesh::Permutation* p, const Connectivity& c) {
    for (size_t i = 0; i < c.toEntities.size(); ++i) {
      EXPECT_EQ(c.toEntities[i], e[i]);
      EXPECT_EQ(c.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
      EXPECT_EQ(c.permutations[i], p[i]);
    }
  }

  std::map<stk::mesh::Entity, unsigned> capture_connectivity_counts(
      const std::vector<stk::mesh::Entity>& ents, stk::mesh::EntityRank connectedRank) {
    std::map<stk::mesh::Entity, unsigned> counts;
    for (auto& ent : ents) {
      counts[ent] = m_bulk->num_connectivity(ent, connectedRank);
    }
    return counts;
  }

  void build_crs_views(const std::vector<stk::mesh::Entity>& fromEntities,
                       const std::vector<Connectivity>& conns,
                       HostEntitiesType& fromView,
                       Kokkos::View<unsigned*, stk::ngp::HostExecSpace>& offsets,
                       HostEntitiesType& toView,
                       Kokkos::View<stk::mesh::RelationIdentifier*, stk::ngp::HostExecSpace>& ordView,
                       Kokkos::View<stk::mesh::Permutation*, stk::ngp::HostExecSpace>& permView) {
    STK_ThrowRequire(fromEntities.size() == conns.size());
    unsigned totalRelations = 0;
    for (const Connectivity& c : conns) { totalRelations += c.toEntities.size(); }

    Kokkos::resize(fromView, fromEntities.size());
    Kokkos::resize(offsets, fromEntities.empty() ? 0 : fromEntities.size() + 1);
    Kokkos::resize(toView, totalRelations);
    Kokkos::resize(ordView, totalRelations);
    Kokkos::resize(permView, totalRelations);

    unsigned idx = 0;
    for (size_t f = 0; f < fromEntities.size(); ++f) {
      fromView(f) = fromEntities[f];
      offsets(f) = idx;
      for (size_t i = 0; i < conns[f].toEntities.size(); ++i, ++idx) {
        toView(idx) = conns[f].toEntities[i];
        ordView(idx) = conns[f].ordinals[i];
        permView(idx) = conns[f].permutations[i];
      }
    }
    if (!fromEntities.empty()) { offsets(fromEntities.size()) = idx; }
  }

  void build_2d_views(const std::vector<stk::mesh::Entity>& fromEntities,
                      const std::vector<Connectivity>& conns,
                      HostEntitiesType& fromView,
                      Host2dEntitiesType& toView,
                      Host2dPermutationsType& permView) {
    STK_ThrowRequire(fromEntities.size() == conns.size());
    const size_t numFrom = fromEntities.size();
    const size_t numCols = conns.empty() ? 0 : conns[0].toEntities.size();

    Kokkos::resize(fromView, numFrom);
    Kokkos::resize(toView, numFrom, numCols);
    Kokkos::resize(permView, numFrom, numCols);
    for (size_t f = 0; f < numFrom; ++f) {
      fromView(f) = fromEntities[f];
      STK_ThrowRequire(conns[f].toEntities.size() == numCols);
      for (size_t i = 0; i < numCols; ++i) {
        const stk::mesh::RelationIdentifier ord = conns[f].ordinals[i];
        STK_ThrowRequire(ord < numCols);
        toView(f, ord) = conns[f].toEntities[i];
        permView(f, ord) = conns[f].permutations[i];
      }
    }
  }

  void build_device_crs_views(const std::vector<stk::mesh::Entity>& fromEntities,
                              const std::vector<Connectivity>& conns,
                              DeviceEntitiesType& fromView,
                              DeviceOffsetsType& offsets,
                              DeviceEntitiesType& toView,
                              DeviceOrdinalsType& ordView,
                              DevicePermutationsType& permView) {
    HostEntitiesType hFrom;
    HostOffsetsType hOff;
    HostEntitiesType hTo;
    HostOrdinalsType hOrd;
    HostPermutationsType hPerm;
    build_crs_views(fromEntities, conns, hFrom, hOff, hTo, hOrd, hPerm);

    Kokkos::resize(fromView, hFrom.extent(0));
    Kokkos::resize(offsets, hOff.extent(0));
    Kokkos::resize(toView, hTo.extent(0));
    Kokkos::resize(ordView, hOrd.extent(0));
    Kokkos::resize(permView, hPerm.extent(0));
    Kokkos::deep_copy(fromView, hFrom);
    Kokkos::deep_copy(offsets, hOff);
    Kokkos::deep_copy(toView, hTo);
    Kokkos::deep_copy(ordView, hOrd);
    Kokkos::deep_copy(permView, hPerm);
  }

  void build_device_2d_views(const std::vector<stk::mesh::Entity>& fromEntities,
                             const std::vector<Connectivity>& conns,
                             DeviceEntitiesType& fromView,
                             Device2dEntitiesType& toView,
                             Device2dPermutationsType& permView) {
    HostEntitiesType hFrom;
    Host2dEntitiesType hTo;
    Host2dPermutationsType hPerm;
    build_2d_views(fromEntities, conns, hFrom, hTo, hPerm);

    Kokkos::resize(fromView, hFrom.extent(0));
    Kokkos::resize(toView, hTo.extent(0), hTo.extent(1));
    Kokkos::resize(permView, hPerm.extent(0), hPerm.extent(1));

    Kokkos::deep_copy(fromView, hFrom);
    auto mTo = Kokkos::create_mirror_view(toView);
    auto mPerm = Kokkos::create_mirror_view(permView);
    for (size_t f = 0; f < hTo.extent(0); ++f) {
      for (size_t c = 0; c < hTo.extent(1); ++c) {
        mTo(f, c) = hTo(f, c);
        mPerm(f, c) = hPerm(f, c);
      }
    }
    Kokkos::deep_copy(toView, mTo);
    Kokkos::deep_copy(permView, mPerm);
  }

  void clear_relations_via_host(const std::vector<stk::mesh::Entity>& fromEntities,
                                stk::mesh::EntityRank connectedRank) {
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
    stk::mesh::HostMesh hostMesh(*m_bulk);
    HostEntitiesType destroyView("clearView", fromEntities.size());
    for (size_t i = 0; i < fromEntities.size(); ++i) { destroyView(i) = fromEntities[i]; }
    hostMesh.batch_destroy_relations(destroyView, connectedRank);
    hostMesh.update_bulk_data();
  }
};

TEST_F(NgpBatchDeclareRelations, hexAndQuadToEdges_raggedRoundTrip_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity quad = quads[0];

  Connectivity hexConn = capture_connectivity(hex, stk::topology::EDGE_RANK);
  Connectivity quadConn = capture_connectivity(quad, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexConn.toEntities.size());
  ASSERT_EQ(4u, quadConn.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  std::vector<stk::mesh::Entity> fromEntities{hex, quad};
  HostEntitiesType destroyView(fromEntities.data(), fromEntities.size());
  hostMesh.batch_destroy_relations(destroyView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_edges(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(quad));

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex, quad}, {hexConn, quadConn}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  ASSERT_EQ(4u, m_bulk->num_edges(quad));

  auto verify = [&](stk::mesh::Entity from, const Connectivity& conn) {
    const stk::mesh::Entity* e = m_bulk->begin_edges(from);
    const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(from);
    const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(from);
    for (size_t i = 0; i < conn.toEntities.size(); ++i) {
      EXPECT_EQ(conn.toEntities[i], e[i]);
      EXPECT_EQ(conn.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
      EXPECT_EQ(conn.permutations[i], p[i]);
    }
  };
  verify(hex, hexConn);
  verify(quad, quadConn);
}

TEST_F(NgpBatchDeclareRelations, hexToEdges_noPermutationOverload_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  Connectivity hexConn = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexConn.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType destroyView(&hex, 1);
  hostMesh.batch_destroy_relations(destroyView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  const unsigned numEdges = hexConn.toEntities.size();
  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexConn}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(numEdges, m_bulk->num_edges(hex));
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < numEdges; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
}

TEST_F(NgpBatchDeclareRelations, hexToFaces_partInductionRoundTrip_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
}

TEST_F(NgpBatchDeclareRelations, twoHex_multipleFromEntities_partInduction_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  std::vector<Connectivity> hexFaceConns{capture_connectivity(hexes[0], stk::topology::FACE_RANK),
                                         capture_connectivity(hexes[1], stk::topology::FACE_RANK)};
  ASSERT_EQ(6u, hexFaceConns[0].toEntities.size());
  ASSERT_EQ(6u, hexFaceConns[1].toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), hexes.size());
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hexes[0], hexes[1]}, hexFaceConns, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, raggedWithEmptySlice_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  Connectivity emptyConn;
  ASSERT_EQ(6u, hex0Faces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), hexes.size());
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_faces(hexes[0]));
  ASSERT_EQ(0u, m_bulk->num_faces(hexes[1]));

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hexes[0], hexes[1]}, {hex0Faces, emptyConn},
                  fromView, offsets, toView, ordView, permView);
  ASSERT_EQ(3u, offsets.extent(0));
  EXPECT_EQ(offsets(1), offsets(2));

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(0u, m_bulk->num_faces(hexes[1]));
  for (const stk::mesh::Entity& face : hex0Faces.toEntities) {
    EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, emptyInput_noOp_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];
  const unsigned origNumFaces = m_bulk->num_faces(hex);
  const unsigned origNumEdges = m_bulk->num_edges(hex);

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({}, {}, fromView, offsets, toView, ordView, permView);
  ASSERT_EQ(0u, fromView.extent(0));

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(origNumFaces, m_bulk->num_faces(hex));
  EXPECT_EQ(origNumEdges, m_bulk->num_edges(hex));
}

TEST_F(NgpBatchDeclareRelations, noPermutations_inducesParts_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_edges(hex));
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexEdges}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView);
  hostMesh.update_bulk_data();

  const unsigned numEdges = m_bulk->num_edges(hex);
  ASSERT_EQ(12u, numEdges);
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < numEdges; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
}

TEST_F(NgpBatchDeclareRelations, forceNoInduceElementPart_notInduced_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Part& blockPart = m_meta->declare_part("blockPart", stk::topology::ELEM_RANK);
  stk::mesh::Part& noInducePart = m_meta->declare_part("noInducePart", stk::topology::ELEM_RANK);
  m_meta->force_no_induce(noInducePart);

  stk::mesh::Entity hex = hexes[0];
  std::vector<stk::mesh::Entity> hexVec{hex};
  m_bulk->modification_begin();
  m_bulk->change_entity_parts(hexVec, stk::mesh::PartVector{&blockPart, &noInducePart});
  m_bulk->modification_end();

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(blockPart));
    EXPECT_FALSE(m_bulk->bucket(quad).member(noInducePart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(blockPart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(blockPart));
    EXPECT_FALSE(m_bulk->bucket(quad).member(noInducePart));
  }
}

TEST_F(NgpBatchDeclareRelations, topologyRootPartInduced_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Part& hexRootPart = m_meta->get_topology_root_part(stk::topology::HEX_8);
  stk::mesh::Entity hex = hexes[0];
  ASSERT_TRUE(m_bulk->bucket(hex).member(hexRootPart));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexRootPart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexRootPart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexRootPart));
  }
}

TEST_F(NgpBatchDeclareRelations, sharedFaceMultiParentUnion_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity sharedFace = stk::mesh::Entity();
  for (auto& quad : quads) {
    if (m_bulk->num_elements(quad) == 2) { sharedFace = quad; break; }
  }
  ASSERT_TRUE(m_bulk->is_valid(sharedFace));

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hex0Faces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hex0View(&hexes[0], 1);
  hostMesh.batch_destroy_relations(hex0View, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  for (auto& face : hex0Faces.toEntities) {
    if (face == sharedFace) {
      EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
    } else {
      EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
    }
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hexes[0]}, {hex0Faces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  for (auto& face : hex0Faces.toEntities) {
    EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, redeclareExistingRelations_isNoOp_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_face_permutations(hex);
  for (size_t i = 0; i < hexFaces.toEntities.size(); ++i) {
    EXPECT_EQ(hexFaces.toEntities[i], e[i]);
    EXPECT_EQ(hexFaces.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexFaces.permutations[i], p[i]);
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, inductionDoesNotChainAcrossRanks_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
  }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexEdges}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
  }
}

TEST_F(NgpBatchDeclareRelations, hexToFaces_2dView_roundTrip_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_faces(hex));

  HostEntitiesType fromView;
  Host2dEntitiesType toView;
  Host2dPermutationsType permView;
  build_2d_views({hex}, {hexFaces}, fromView, toView, permView);

  hostMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_face_permutations(hex);
  for (size_t i = 0; i < hexFaces.toEntities.size(); ++i) {
    EXPECT_EQ(hexFaces.toEntities[i], e[i]);
    EXPECT_EQ(hexFaces.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexFaces.permutations[i], p[i]);
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, noPermutationOverload_2dView_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  HostEntitiesType fromView;
  Host2dEntitiesType toView;
  Host2dPermutationsType permView;
  build_2d_views({hex}, {hexEdges}, fromView, toView, permView);

  hostMesh.batch_declare_relations(fromView, toView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < 12u; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
}

TEST_F(NgpBatchDeclareRelations, multipleUniformFromEntities_2dView_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  std::vector<Connectivity> hexFaceConns{capture_connectivity(hexes[0], stk::topology::FACE_RANK),
                                         capture_connectivity(hexes[1], stk::topology::FACE_RANK)};
  ASSERT_EQ(6u, hexFaceConns[0].toEntities.size());
  ASSERT_EQ(6u, hexFaceConns[1].toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexesView(hexes.data(), hexes.size());
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }

  HostEntitiesType fromView;
  Host2dEntitiesType toView;
  Host2dPermutationsType permView;
  build_2d_views({hexes[0], hexes[1]}, hexFaceConns, fromView, toView, permView);

  hostMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

TEST_F(NgpBatchDeclareRelations, positionalOrdinals_2dView_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  HostEntitiesType fromView;
  Host2dEntitiesType toView;
  Host2dPermutationsType permView;
  build_2d_views({hex}, {hexFaces}, fromView, toView, permView);

  hostMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  for (unsigned i = 0; i < 6u; ++i) {
    EXPECT_EQ(i, static_cast<unsigned>(o[i]));
    EXPECT_EQ(toView(0, i), e[i]);
  }
}

#ifndef NDEBUG
TEST_F(NgpBatchDeclareRelations, wrongComplementCount_2dView_throwsInDebug_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType fromView("from", 1);
  Host2dEntitiesType toView("to", 1, 5);
  Host2dPermutationsType permView("perm", 1, 5);
  fromView(0) = hex;
  for (unsigned j = 0; j < 5u; ++j) {
    toView(0, j) = hexFaces.toEntities[j];
    permView(0, j) = hexFaces.permutations[j];
  }

  EXPECT_ANY_THROW(hostMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK));
}
#endif

TEST_F(NgpBatchDeclareRelations, mixedRanksInSingleBatch_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  Connectivity mixed;
  mixed.append(hexFaces);
  mixed.append(hexEdges);

  stk::mesh::HostMesh hostMesh(*m_bulk);

  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(hex));
  for (auto& quad : quads) { EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart)); }
  for (auto& edge : edges) { EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart)); }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {mixed}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hex));
  EXPECT_EQ(12u, m_bulk->num_edges(hex));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart)); }
}

TEST_F(NgpBatchDeclareRelations, partialRedeclare_mixNewAndExisting_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  const std::vector<unsigned> removedIdx{0, 2, 4};
  m_bulk->modification_begin();
  for (unsigned i : removedIdx) {
    m_bulk->destroy_relation(hex, hexFaces.toEntities[i], hexFaces.ordinals[i]);
  }
  m_bulk->modification_end();
  ASSERT_EQ(3u, m_bulk->num_faces(hex));
  EXPECT_FALSE(m_bulk->bucket(hexFaces.toEntities[0]).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(hexFaces.toEntities[1]).member(hexPart));

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  EXPECT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
}

TEST_F(NgpBatchDeclareRelations, bucketGrowthAndCreate_acrossTwoDeclares_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(8, 8);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  Connectivity hex1Faces = capture_connectivity(hexes[1], stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hex0Faces.toEntities.size());
  ASSERT_EQ(6u, hex1Faces.toEntities.size());

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType hexesView(hexes.data(), hexes.size());
  hostMesh.batch_destroy_relations(hexesView, stk::topology::FACE_RANK);
  hostMesh.update_bulk_data();
  for (auto& quad : quads) { EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart)); }

  {
    HostEntitiesType fromView; HostOffsetsType offsets; HostEntitiesType toView;
    HostOrdinalsType ordView; HostPermutationsType permView;
    build_crs_views({hexes[0]}, {hex0Faces}, fromView, offsets, toView, ordView, permView);
    hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    hostMesh.update_bulk_data();
  }
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));

  {
    HostEntitiesType fromView; HostOffsetsType offsets; HostEntitiesType toView;
    HostOrdinalsType ordView; HostPermutationsType permView;
    build_crs_views({hexes[1]}, {hex1Faces}, fromView, offsets, toView, ordView, permView);
    hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    hostMesh.update_bulk_data();
  }
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
}

TEST_F(NgpBatchDeclareRelations, highFanInSharedEdges_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");

  std::vector<Connectivity> faceEdgeConns;
  for (auto& quad : quads) {
    Connectivity c = capture_connectivity(quad, stk::topology::EDGE_RANK);
    ASSERT_EQ(4u, c.toEntities.size());
    faceEdgeConns.push_back(c);
  }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart)); }

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType quadsView(quads.data(), quads.size());
  hostMesh.batch_destroy_relations(quadsView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  for (auto& edge : edges) { EXPECT_FALSE(m_bulk->bucket(edge).member(quadPart)); }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views(quads, faceEdgeConns, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) { EXPECT_EQ(4u, m_bulk->num_edges(quad)); }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart)); }
}

TEST_F(NgpBatchDeclareRelations, orderPreservedNonMonotonicOrdinals_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // Feed the edges in reversed (descending-ordinal) order.  STK canonicalizes connectivity by
  // ordinal, so the stored order comes back ascending regardless of input order (asserted below).
  // This is the direct regression for the stable composite sort key: an unstable sort over the
  // non-monotonic input would corrupt the intermediate device stages before canonicalization.
  Connectivity reversed = hexEdges.reversed_copy();

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {reversed}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  // STK canonicalizes connectivity by ordinal, so the stored order is ascending-ordinal (the
  // original `hexEdges` order) regardless of the reversed input order we fed in.  The regression is
  // that the non-monotonic input is handled without corruption and the device path produces the
  // same canonical result as the (unchanged) host reference.
  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::Entity* e = m_bulk->begin_edges(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    EXPECT_EQ(hexEdges.toEntities[i], e[i]);
    EXPECT_EQ(hexEdges.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexEdges.permutations[i], p[i]);
  }
}

TEST_F(NgpBatchDeclareRelations, multiRankSingleOwner_sliceContents_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // One from-entity gaining connectivity at two connected ranks in one batch: exercises the
  // per-rank gap-shift in insert_connectivities_no_grow (the edge slice must be shifted up when the
  // face slice grows) and the exact per-owner size reduction.
  Connectivity mixed;
  mixed.append(hexFaces);
  mixed.append(hexEdges);

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType hexView(&hex, 1);
  hostMesh.batch_destroy_relations(hexView, stk::topology::FACE_RANK);
  hostMesh.batch_destroy_relations(hexView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {mixed}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  ASSERT_EQ(12u, m_bulk->num_edges(hex));

  verify_slice(m_bulk->begin_faces(hex), m_bulk->begin_face_ordinals(hex),
               m_bulk->begin_face_permutations(hex), hexFaces);
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

TEST_F(NgpBatchDeclareRelations, highFanInExactReciprocalCounts_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();

  std::vector<Connectivity> faceEdgeConns;
  for (auto& quad : quads) {
    Connectivity c = capture_connectivity(quad, stk::topology::EDGE_RANK);
    ASSERT_EQ(4u, c.toEntities.size());
    faceEdgeConns.push_back(c);
  }

  // Record the exact reciprocal (edge->face) fan-in each edge had before removal, so we can assert
  // the two-pass per-owner sizing reproduces it exactly (no over- or under-allocation).
  std::map<stk::mesh::Entity, unsigned> expectedEdgeFaceCount =
      capture_connectivity_counts(edges, stk::topology::FACE_RANK);

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType quadsView(quads.data(), quads.size());
  hostMesh.batch_destroy_relations(quadsView, stk::topology::EDGE_RANK);
  hostMesh.update_bulk_data();
  for (auto& edge : edges) { EXPECT_EQ(0u, m_bulk->num_connectivity(edge, stk::topology::FACE_RANK)); }

  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views(quads, faceEdgeConns, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  for (auto& quad : quads) { EXPECT_EQ(4u, m_bulk->num_edges(quad)); }
  for (auto& edge : edges) {
    EXPECT_EQ(expectedEdgeFaceCount[edge], m_bulk->num_connectivity(edge, stk::topology::FACE_RANK));
  }
}

TEST_F(NgpBatchDeclareRelations, partialRedeclare_sliceStaysSorted_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // Destroy every edge, then re-declare the odd ordinals first and the even ordinals second, so the
  // second half of the batch interleaves the first.
  Connectivity oddOrdinals, evenOrdinals;
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    Connectivity& dst = (hexEdges.ordinals[i] % 2 == 1) ? oddOrdinals : evenOrdinals;
    dst.toEntities.push_back(hexEdges.toEntities[i]);
    dst.ordinals.push_back(hexEdges.ordinals[i]);
    dst.permutations.push_back(hexEdges.permutations[i]);
  }

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  Connectivity oddThenEven;
  oddThenEven.append(oddOrdinals);
  oddThenEven.append(evenOrdinals);

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {oddThenEven}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(hex);
  for (unsigned i = 0; i < 12u; ++i) {
    EXPECT_EQ(i, static_cast<unsigned>(o[i]))
        << "host edge slice not in ascending-ordinal order at slot " << i;
  }
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

TEST_F(NgpBatchDeclareRelations, reversedDeclare_sliceStaysSorted_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  // Feed the relations in descending-ordinal order into an empty slice.  The stored order must
  // still come back ascending.
  Connectivity reversed = hexEdges.reversed_copy();

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {reversed}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(hex);
  for (unsigned i = 0; i < 12u; ++i) {
    EXPECT_EQ(i, static_cast<unsigned>(o[i]))
        << "host edge slice not in ascending-ordinal order at slot " << i;
  }
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

TEST_F(NgpBatchDeclareRelations, partialRedeclare_upwardSliceStaysSorted_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  Connectivity oddOrdinals, evenOrdinals;
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    Connectivity& dst = (hexEdges.ordinals[i] % 2 == 1) ? oddOrdinals : evenOrdinals;
    dst.toEntities.push_back(hexEdges.toEntities[i]);
    dst.ordinals.push_back(hexEdges.ordinals[i]);
    dst.permutations.push_back(hexEdges.permutations[i]);
  }

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);

  Connectivity oddThenEven;
  oddThenEven.append(oddOrdinals);
  oddThenEven.append(evenOrdinals);

  stk::mesh::HostMesh hostMesh(*m_bulk);
  HostEntitiesType fromView;
  HostOffsetsType offsets;
  HostEntitiesType toView;
  HostOrdinalsType ordView;
  HostPermutationsType permView;
  build_crs_views({hex}, {oddThenEven}, fromView, offsets, toView, ordView, permView);

  hostMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
  hostMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    stk::mesh::Entity edge = hexEdges.toEntities[i];
    ASSERT_EQ(1u, m_bulk->num_connectivity(edge, stk::topology::ELEM_RANK));
    const stk::mesh::ConnectivityOrdinal* upOrds =
        m_bulk->begin_ordinals(edge, stk::topology::ELEM_RANK);
    const unsigned numUp = m_bulk->num_connectivity(edge, stk::topology::ELEM_RANK);
    for (unsigned s = 1; s < numUp; ++s) {
      EXPECT_LE(static_cast<unsigned>(upOrds[s-1]), static_cast<unsigned>(upOrds[s]))
          << "host elem slice of edge " << i << " not ascending at slot " << s;
    }
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, hexAndQuadToEdges_raggedRoundTrip_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  stk::mesh::Entity quad = quads[0];

  Connectivity hexConn = capture_connectivity(hex, stk::topology::EDGE_RANK);
  Connectivity quadConn = capture_connectivity(quad, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexConn.toEntities.size());
  ASSERT_EQ(4u, quadConn.toEntities.size());

  clear_relations_via_host({hex, quad}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(quad));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex, quad}, {hexConn, quadConn}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  ASSERT_EQ(4u, m_bulk->num_edges(quad));

  auto verify = [&](stk::mesh::Entity from, const Connectivity& conn) {
    const stk::mesh::Entity* e = m_bulk->begin_edges(from);
    const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(from);
    const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(from);
    for (size_t i = 0; i < conn.toEntities.size(); ++i) {
      EXPECT_EQ(conn.toEntities[i], e[i]);
      EXPECT_EQ(conn.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
      EXPECT_EQ(conn.permutations[i], p[i]);
    }
  };
  verify(hex, hexConn);
  verify(quad, quadConn);
}

NGP_TEST_F(NgpBatchDeclareRelations, hexToEdges_noPermutationOverload_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Entity hex = hexes[0];
  Connectivity hexConn = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexConn.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  const unsigned numEdges = hexConn.toEntities.size();
  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexConn}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(numEdges, m_bulk->num_edges(hex));
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < numEdges; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, hexToFaces_partInductionRoundTrip_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(quad).member(quadPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, twoHex_multipleFromEntities_partInduction_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  std::vector<Connectivity> hexFaceConns{capture_connectivity(hexes[0], stk::topology::FACE_RANK),
                                         capture_connectivity(hexes[1], stk::topology::FACE_RANK)};
  ASSERT_EQ(6u, hexFaceConns[0].toEntities.size());
  ASSERT_EQ(6u, hexFaceConns[1].toEntities.size());

  clear_relations_via_host({hexes[0], hexes[1]}, stk::topology::FACE_RANK);
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hexes[0], hexes[1]}, hexFaceConns, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, raggedWithEmptySlice_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  Connectivity emptyConn;
  ASSERT_EQ(6u, hex0Faces.toEntities.size());

  clear_relations_via_host({hexes[0], hexes[1]}, stk::topology::FACE_RANK);
  ASSERT_EQ(0u, m_bulk->num_faces(hexes[0]));
  ASSERT_EQ(0u, m_bulk->num_faces(hexes[1]));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hexes[0], hexes[1]}, {hex0Faces, emptyConn},
                         fromView, offsets, toView, ordView, permView);
  ASSERT_EQ(3u, offsets.extent(0));

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(0u, m_bulk->num_faces(hexes[1]));
  for (const stk::mesh::Entity& face : hex0Faces.toEntities) {
    EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, emptyInput_noOp_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];
  const unsigned origNumFaces = m_bulk->num_faces(hex);
  const unsigned origNumEdges = m_bulk->num_edges(hex);

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({}, {}, fromView, offsets, toView, ordView, permView);
  ASSERT_EQ(0u, fromView.extent(0));

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(origNumFaces, m_bulk->num_faces(hex));
  EXPECT_EQ(origNumEdges, m_bulk->num_edges(hex));
}

NGP_TEST_F(NgpBatchDeclareRelations, noPermutations_inducesParts_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& edgePart = *m_meta->get_part("edgePart");
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexEdges}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView);
    ngpMesh.update_bulk_data();
  }

  const unsigned numEdges = m_bulk->num_edges(hex);
  ASSERT_EQ(12u, numEdges);
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < numEdges; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(edgePart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, forceNoInduceElementPart_notInduced_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Part& blockPart = m_meta->declare_part("blockPart", stk::topology::ELEM_RANK);
  stk::mesh::Part& noInducePart = m_meta->declare_part("noInducePart", stk::topology::ELEM_RANK);
  m_meta->force_no_induce(noInducePart);

  stk::mesh::Entity hex = hexes[0];
  std::vector<stk::mesh::Entity> hexVec{hex};
  m_bulk->modification_begin();
  m_bulk->change_entity_parts(hexVec, stk::mesh::PartVector{&blockPart, &noInducePart});
  m_bulk->modification_end();

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(blockPart));
    EXPECT_FALSE(m_bulk->bucket(quad).member(noInducePart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(blockPart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(blockPart));
    EXPECT_FALSE(m_bulk->bucket(quad).member(noInducePart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, topologyRootPartInduced_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();

  stk::mesh::Part& hexRootPart = m_meta->get_topology_root_part(stk::topology::HEX_8);
  stk::mesh::Entity hex = hexes[0];
  ASSERT_TRUE(m_bulk->bucket(hex).member(hexRootPart));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexRootPart));
  }

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexRootPart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexRootPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, sharedFaceMultiParentUnion_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");

  stk::mesh::Entity sharedFace = stk::mesh::Entity();
  for (auto& quad : quads) {
    if (m_bulk->num_elements(quad) == 2) { sharedFace = quad; break; }
  }
  ASSERT_TRUE(m_bulk->is_valid(sharedFace));

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hex0Faces.toEntities.size());

  clear_relations_via_host({hexes[0]}, stk::topology::FACE_RANK);

  for (auto& face : hex0Faces.toEntities) {
    if (face == sharedFace) {
      EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
    } else {
      EXPECT_FALSE(m_bulk->bucket(face).member(hexPart));
    }
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hexes[0]}, {hex0Faces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  for (auto& face : hex0Faces.toEntities) {
    EXPECT_TRUE(m_bulk->bucket(face).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, redeclareExistingRelations_isNoOp_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_face_permutations(hex);
  for (size_t i = 0; i < hexFaces.toEntities.size(); ++i) {
    EXPECT_EQ(hexFaces.toEntities[i], e[i]);
    EXPECT_EQ(hexFaces.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexFaces.permutations[i], p[i]);
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, inductionDoesNotChainAcrossRanks_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
  for (auto& edge : edges) {
    EXPECT_FALSE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
  }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexEdges}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  for (auto& edge : edges) {
    EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart));
    EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, hexToFaces_2dView_roundTrip_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  ASSERT_EQ(0u, m_bulk->num_faces(hex));

  DeviceEntitiesType fromView;
  Device2dEntitiesType toView;
  Device2dPermutationsType permView;
  build_device_2d_views({hex}, {hexFaces}, fromView, toView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_face_permutations(hex);
  for (size_t i = 0; i < hexFaces.toEntities.size(); ++i) {
    EXPECT_EQ(hexFaces.toEntities[i], e[i]);
    EXPECT_EQ(hexFaces.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexFaces.permutations[i], p[i]);
  }
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, noPermutationOverload_2dView_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  DeviceEntitiesType fromView;
  Device2dEntitiesType toView;
  Device2dPermutationsType permView;
  build_device_2d_views({hex}, {hexEdges}, fromView, toView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, toView, stk::topology::EDGE_RANK);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (unsigned i = 0; i < 12u; ++i) {
    EXPECT_EQ(stk::mesh::Permutation::INVALID_PERMUTATION, p[i]);
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, multipleUniformFromEntities_2dView_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  std::vector<Connectivity> hexFaceConns{capture_connectivity(hexes[0], stk::topology::FACE_RANK),
                                         capture_connectivity(hexes[1], stk::topology::FACE_RANK)};
  ASSERT_EQ(6u, hexFaceConns[0].toEntities.size());
  ASSERT_EQ(6u, hexFaceConns[1].toEntities.size());

  clear_relations_via_host({hexes[0], hexes[1]}, stk::topology::FACE_RANK);
  for (auto& quad : quads) {
    EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart));
  }

  DeviceEntitiesType fromView;
  Device2dEntitiesType toView;
  Device2dPermutationsType permView;
  build_device_2d_views({hexes[0], hexes[1]}, hexFaceConns, fromView, toView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) {
    EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, positionalOrdinals_2dView_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);

  DeviceEntitiesType fromView;
  Device2dEntitiesType toView;
  Device2dPermutationsType permView;
  build_device_2d_views({hex}, {hexFaces}, fromView, toView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  const stk::mesh::Entity* e = m_bulk->begin_faces(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_face_ordinals(hex);
  for (unsigned i = 0; i < 6u; ++i) {
    EXPECT_EQ(i, static_cast<unsigned>(o[i]));
    EXPECT_EQ(hexFaces.toEntities[i], e[i]);
  }
}

#ifndef NDEBUG
NGP_TEST_F(NgpBatchDeclareRelations, wrongComplementCount_2dView_throwsInDebug_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  DeviceEntitiesType fromView("from", 1);
  Device2dEntitiesType toView("to", 1, 5);
  Device2dPermutationsType permView("perm", 1, 5);
  auto hFrom = Kokkos::create_mirror_view(fromView);
  auto hTo = Kokkos::create_mirror_view(toView);
  auto hPerm = Kokkos::create_mirror_view(permView);
  hFrom(0) = hex;
  for (unsigned j = 0; j < 5u; ++j) {
    hTo(0, j) = hexFaces.toEntities[j];
    hPerm(0, j) = hexFaces.permutations[j];
  }
  Kokkos::deep_copy(fromView, hFrom);
  Kokkos::deep_copy(toView, hTo);
  Kokkos::deep_copy(permView, hPerm);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  EXPECT_ANY_THROW(ngpMesh.batch_declare_relations(fromView, toView, permView, stk::topology::FACE_RANK));
}
#endif

NGP_TEST_F(NgpBatchDeclareRelations, mixedRanksInSingleBatch_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // One CRS slice carrying two different to-ranks (faces and edges) for the same from-entity.
  Connectivity mixed;
  mixed.append(hexFaces);
  mixed.append(hexEdges);

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {mixed}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  EXPECT_EQ(6u, m_bulk->num_faces(hex));
  EXPECT_EQ(12u, m_bulk->num_edges(hex));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(hexPart)); }
}

NGP_TEST_F(NgpBatchDeclareRelations, partialRedeclare_mixNewAndExisting_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());

  // Build the device mesh while the full face connectivity is present so the device connectivity
  // pool is sized for all 6 faces before we remove (and later re-add) a subset.  (Same rationale
  // as clear_relations_via_host.)
  stk::mesh::get_updated_ngp_mesh(*m_bulk);

  // Remove a subset of the hex's faces on the host, leaving the rest in place.  Re-declaring all 6
  // then feeds the device a batch that mixes already-existing relations (must be deduped) with new
  // ones (must be added) -- exercising the keep-mask compaction path.
  const std::vector<unsigned> removedIdx{0, 2, 4};
  m_bulk->modification_begin();
  for (unsigned i : removedIdx) {
    m_bulk->destroy_relation(hex, hexFaces.toEntities[i], hexFaces.ordinals[i]);
  }
  m_bulk->modification_end();
  ASSERT_EQ(3u, m_bulk->num_faces(hex));
  EXPECT_FALSE(m_bulk->bucket(hexFaces.toEntities[0]).member(hexPart));
  EXPECT_TRUE(m_bulk->bucket(hexFaces.toEntities[1]).member(hexPart));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {hexFaces}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  // No duplicate relations: exactly the original 6 distinct faces, each induced with hexPart.
  EXPECT_EQ(6u, m_bulk->num_faces(hex));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
}

NGP_TEST_F(NgpBatchDeclareRelations, bucketGrowthAndCreate_acrossTwoDeclares_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  // Bucket capacity 8: the first declare fills a face bucket to 6/8; the second declare grows that
  // bucket to 8 and then creates an additional bucket for the overflow -- exercising both
  // batch_grow_buckets and batch_create_buckets through the relation-addition (induction) path.
  build_empty_mesh(8, 8);
  create_two_hex_mesh();
  stk::mesh::Part& hexPart = *m_meta->get_part("hexPart");
  ASSERT_EQ(2u, hexes.size());

  Connectivity hex0Faces = capture_connectivity(hexes[0], stk::topology::FACE_RANK);
  Connectivity hex1Faces = capture_connectivity(hexes[1], stk::topology::FACE_RANK);
  ASSERT_EQ(6u, hex0Faces.toEntities.size());
  ASSERT_EQ(6u, hex1Faces.toEntities.size());

  clear_relations_via_host({hexes[0], hexes[1]}, stk::topology::FACE_RANK);
  for (auto& quad : quads) { EXPECT_FALSE(m_bulk->bucket(quad).member(hexPart)); }

  {
    DeviceEntitiesType fromView, toView; DeviceOffsetsType offsets;
    DeviceOrdinalsType ordView; DevicePermutationsType permView;
    build_device_crs_views({hexes[0]}, {hex0Faces}, fromView, offsets, toView, ordView, permView);
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[0]));
  for (auto& face : hex0Faces.toEntities) { EXPECT_TRUE(m_bulk->bucket(face).member(hexPart)); }

  {
    DeviceEntitiesType fromView, toView; DeviceOffsetsType offsets;
    DeviceOrdinalsType ordView; DevicePermutationsType permView;
    build_device_crs_views({hexes[1]}, {hex1Faces}, fromView, offsets, toView, ordView, permView);
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }
  EXPECT_EQ(6u, m_bulk->num_faces(hexes[1]));
  for (auto& quad : quads) { EXPECT_TRUE(m_bulk->bucket(quad).member(hexPart)); }
}

NGP_TEST_F(NgpBatchDeclareRelations, highFanInSharedEdges_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  // In the two-hex mesh the internal-plane edges are each shared by 3 faces, so re-declaring all
  // face->edge relations in one batch drives the affected-entity fan-in to 3 (>2) and the
  // associated slab-stride sizing.  (face->edge, not face->node: entities above node rank are
  // required to always retain their connected nodes, so face->node relations cannot be removed.)
  build_empty_mesh(1, 1);
  create_two_hex_mesh();
  stk::mesh::Part& quadPart = *m_meta->get_part("quadPart");

  std::vector<Connectivity> faceEdgeConns;
  for (auto& quad : quads) {
    Connectivity c = capture_connectivity(quad, stk::topology::EDGE_RANK);
    ASSERT_EQ(4u, c.toEntities.size());
    faceEdgeConns.push_back(c);
  }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart)); }

  clear_relations_via_host(quads, stk::topology::EDGE_RANK);
  for (auto& edge : edges) { EXPECT_FALSE(m_bulk->bucket(edge).member(quadPart)); }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views(quads, faceEdgeConns, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  for (auto& quad : quads) { EXPECT_EQ(4u, m_bulk->num_edges(quad)); }
  for (auto& edge : edges) { EXPECT_TRUE(m_bulk->bucket(edge).member(quadPart)); }
}

NGP_TEST_F(NgpBatchDeclareRelations, orderPreservedNonMonotonicOrdinals_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // Feed the edges in reversed (descending-ordinal) order.  Connectivity is stored in insertion
  // order, so the stored order must come back reversed.  Direct regression for the stable composite
  // sort key on device: Kokkos::sort is not stable, so a missing seq tiebreaker would scramble this.
  Connectivity reversed = hexEdges.reversed_copy();

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {reversed}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  // STK canonicalizes connectivity by ordinal, so the stored order is ascending-ordinal (the
  // original `hexEdges` order) regardless of the reversed input order we fed in.  The regression is
  // that the non-monotonic input is handled without corruption and the device path produces the
  // same canonical result as the (unchanged) host reference.
  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  const stk::mesh::Entity* e = m_bulk->begin_edges(hex);
  const stk::mesh::ConnectivityOrdinal* o = m_bulk->begin_edge_ordinals(hex);
  const stk::mesh::Permutation* p = m_bulk->begin_edge_permutations(hex);
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    EXPECT_EQ(hexEdges.toEntities[i], e[i]);
    EXPECT_EQ(hexEdges.ordinals[i], static_cast<stk::mesh::RelationIdentifier>(o[i]));
    EXPECT_EQ(hexEdges.permutations[i], p[i]);
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, multiRankSingleOwner_sliceContents_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexFaces = capture_connectivity(hex, stk::topology::FACE_RANK);
  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(6u, hexFaces.toEntities.size());
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // One from-entity gaining connectivity at two connected ranks in one batch: exercises the
  // per-rank gap-shift in insert_connectivities_no_grow (the edge slice is shifted up when the face
  // slice grows) and the exact per-owner size reduction, all on device.
  Connectivity mixed;
  mixed.append(hexFaces);
  mixed.append(hexEdges);

  clear_relations_via_host({hex}, stk::topology::FACE_RANK);
  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_faces(hex));
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {mixed}, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  ASSERT_EQ(6u, m_bulk->num_faces(hex));
  ASSERT_EQ(12u, m_bulk->num_edges(hex));

  verify_slice(m_bulk->begin_faces(hex), m_bulk->begin_face_ordinals(hex),
               m_bulk->begin_face_permutations(hex), hexFaces);
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

NGP_TEST_F(NgpBatchDeclareRelations, highFanInExactReciprocalCounts_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_two_hex_mesh();

  std::vector<Connectivity> faceEdgeConns;
  for (auto& quad : quads) {
    Connectivity c = capture_connectivity(quad, stk::topology::EDGE_RANK);
    ASSERT_EQ(4u, c.toEntities.size());
    faceEdgeConns.push_back(c);
  }

  // Record the exact reciprocal (edge->face) fan-in each edge had before removal, so we can assert
  // the two-pass per-owner sizing on device reproduces it exactly (no over/under allocation): the
  // shared internal-plane edges are each connected to 3 faces.
  std::map<stk::mesh::Entity, unsigned> expectedEdgeFaceCount =
      capture_connectivity_counts(edges, stk::topology::FACE_RANK);

  clear_relations_via_host(quads, stk::topology::EDGE_RANK);
  for (auto& edge : edges) { EXPECT_EQ(0u, m_bulk->num_connectivity(edge, stk::topology::FACE_RANK)); }

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views(quads, faceEdgeConns, fromView, offsets, toView, ordView, permView);

  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);
    ngpMesh.update_bulk_data();
  }

  for (auto& quad : quads) { EXPECT_EQ(4u, m_bulk->num_edges(quad)); }
  for (auto& edge : edges) {
    EXPECT_EQ(expectedEdgeFaceCount[edge], m_bulk->num_connectivity(edge, stk::topology::FACE_RANK));
  }
}

NGP_TEST_F(NgpBatchDeclareRelations, partialRedeclare_deviceSliceStaysSorted_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  // Destroy every edge, then re-declare the odd ordinals first and the even ordinals second.  The
  // second batch interleaves the first, so a plain append would leave the device slice ordered
  // 1 3 5 7 9 11 0 2 4 6 8 10 instead of 0..11.
  Connectivity oddOrdinals, evenOrdinals;
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    Connectivity& dst = (hexEdges.ordinals[i] % 2 == 1) ? oddOrdinals : evenOrdinals;
    dst.toEntities.push_back(hexEdges.toEntities[i]);
    dst.ordinals.push_back(hexEdges.ordinals[i]);
    dst.permutations.push_back(hexEdges.permutations[i]);
  }

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  Connectivity oddThenEven;
  oddThenEven.append(oddOrdinals);
  oddThenEven.append(evenOrdinals);

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {oddThenEven}, fromView, offsets, toView, ordView, permView);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);

  {
    auto& meshConn = ngpMesh.get_mesh_connectivity();
    auto deviceOrdinals = meshConn.get_connected_ordinals(hex, stk::topology::EDGE_RANK);
    ASSERT_EQ(12u, deviceOrdinals.size());
    for (unsigned i = 0; i < deviceOrdinals.size(); ++i) {
      EXPECT_EQ(i, static_cast<unsigned>(deviceOrdinals[i]))
          << "device edge slice not in ascending-ordinal order at slot " << i;
    }
  }

  ngpMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

NGP_TEST_F(NgpBatchDeclareRelations, reversedDeclare_deviceSliceStaysSorted_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);
  ASSERT_EQ(0u, m_bulk->num_edges(hex));

  // Feed the relations in descending-ordinal order into an empty slice.  The result must still be
  // ascending on the device.
  Connectivity reversed = hexEdges.reversed_copy();

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {reversed}, fromView, offsets, toView, ordView, permView);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);

  {
    auto& meshConn = ngpMesh.get_mesh_connectivity();
    auto deviceOrdinals = meshConn.get_connected_ordinals(hex, stk::topology::EDGE_RANK);
    ASSERT_EQ(12u, deviceOrdinals.size());
    for (unsigned i = 0; i < deviceOrdinals.size(); ++i) {
      EXPECT_EQ(i, static_cast<unsigned>(deviceOrdinals[i]))
          << "device edge slice not in ascending-ordinal order at slot " << i;
    }
  }

  ngpMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  verify_slice(m_bulk->begin_edges(hex), m_bulk->begin_edge_ordinals(hex),
               m_bulk->begin_edge_permutations(hex), hexEdges);
}

NGP_TEST_F(NgpBatchDeclareRelations, partialRedeclare_upwardSliceStaysSorted_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  create_one_hex_mesh();
  stk::mesh::Entity hex = hexes[0];

  Connectivity hexEdges = capture_connectivity(hex, stk::topology::EDGE_RANK);
  ASSERT_EQ(12u, hexEdges.toEntities.size());

  Connectivity oddOrdinals, evenOrdinals;
  for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
    Connectivity& dst = (hexEdges.ordinals[i] % 2 == 1) ? oddOrdinals : evenOrdinals;
    dst.toEntities.push_back(hexEdges.toEntities[i]);
    dst.ordinals.push_back(hexEdges.ordinals[i]);
    dst.permutations.push_back(hexEdges.permutations[i]);
  }

  clear_relations_via_host({hex}, stk::topology::EDGE_RANK);

  Connectivity oddThenEven;
  oddThenEven.append(oddOrdinals);
  oddThenEven.append(evenOrdinals);

  DeviceEntitiesType fromView;
  DeviceOffsetsType offsets;
  DeviceEntitiesType toView;
  DeviceOrdinalsType ordView;
  DevicePermutationsType permView;
  build_device_crs_views({hex}, {oddThenEven}, fromView, offsets, toView, ordView, permView);

  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
  ngpMesh.batch_declare_relations(fromView, offsets, toView, ordView, permView);

  {
    auto& meshConn = ngpMesh.get_mesh_connectivity();
    for (size_t i = 0; i < hexEdges.toEntities.size(); ++i) {
      stk::mesh::Entity edge = hexEdges.toEntities[i];
      auto upOrdinals = meshConn.get_connected_ordinals(edge, stk::topology::ELEM_RANK);
      ASSERT_EQ(1u, upOrdinals.size());
      for (unsigned s = 1; s < upOrdinals.size(); ++s) {
        EXPECT_LE(static_cast<unsigned>(upOrdinals[s-1]), static_cast<unsigned>(upOrdinals[s]))
            << "device elem slice of edge " << i << " not ascending at slot " << s;
      }
    }
  }

  ngpMesh.update_bulk_data();

  ASSERT_EQ(12u, m_bulk->num_edges(hex));
  for (stk::mesh::Entity edge : hexEdges.toEntities) {
    EXPECT_EQ(1u, m_bulk->num_connectivity(edge, stk::topology::ELEM_RANK));
  }
}

}  // namespace
#endif
