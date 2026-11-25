// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

// A tet near the edge x = y = 0.5 of a cuboid: A and B are on its +x face, X and C are outside
// beyond its +y face.  Each edge from A or B to X or C goes inside and leaves through y = 0.5, so
// its ends (on the surface and outside) do not show that it crosses the surface.
struct TetWithMultipleEdgeCrossings
{
  static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
  std::vector<stk::math::Vector3d> nodeLocs{{{0.5, 0.3, 0.0}, // A
      {0.45, 0.7, 0.0},                                       // X
      {0.3, 0.6, 0.2},                                        // C
      {0.5, 0.3, 0.4}}};                                      // B
  std::vector<std::array<unsigned, 4>> allElementConn{{{0, 1, 2, 3}}};
};

double tet_volume(const std::array<stk::math::Vector3d, 4> & x)
{
  return std::abs(Dot(x[3] - x[0], Cross(x[1] - x[0], x[2] - x[0]))) / 6.;
}

// The volume of the part of the tet ABXC with y below yPlane, for A and B below and X and C above:
// the prism between the triangles A, P(A,X), P(A,C) and B, P(B,X), P(B,C).
double volume_below_plane_y(const stk::math::Vector3d & a,
    const stk::math::Vector3d & b,
    const stk::math::Vector3d & x,
    const stk::math::Vector3d & c,
    const double yPlane)
{
  auto crossing = [yPlane](const stk::math::Vector3d & p, const stk::math::Vector3d & q)
  { return p + ((yPlane - p[1]) / (q[1] - p[1])) * (q - p); };
  const std::array<stk::math::Vector3d, 3> a3{{a, crossing(a, x), crossing(a, c)}};
  const std::array<stk::math::Vector3d, 3> b3{{b, crossing(b, x), crossing(b, c)}};
  return tet_volume({{a3[0], a3[1], a3[2], b3[2]}}) + tet_volume({{a3[0], a3[1], b3[1], b3[2]}}) +
      tet_volume({{a3[0], b3[0], b3[1], b3[2]}});
}

template <typename MESHSPEC>
class MultipleEdgeCrossingsDecompositionFixture
    : public DecompositionFixture<MESHSPEC, LSPerInterfacePolicy, 1>
{
public:
  MultipleEdgeCrossingsDecompositionFixture()
  {
    if (stk::parallel_machine_size(this->mComm) > 1) return;
    this->cdfem_support().set_cdfem_edge_degeneracy_handling(
        SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
    this->cdfem_support().register_cdfem_snap_displacements_field();
    this->build_mesh_from_spec(this->meshSpec);
    this->setup_ls_fields();
  }

  double decompose_and_compute_inside_volume(const Surface & surface)
  {
    const stk::mesh::MetaData & meta = this->mMesh.mesh_meta_data();
    const Phase_Support & phaseSupport = Phase_Support::get(meta);
    const Surface_Identifier surfaceIdentifier = this->levelset_fields()[0].identifier;
    const AnalyticSurfaceInterfaceGeometry geometry({surfaceIdentifier},
        {&surface},
        this->get_aux_meta().active_part(),
        this->cdfem_support(),
        phaseSupport);
    CDMesh::decompose_mesh(this->mMesh, geometry);

    double volume = compute_volume_or_surface_area(this->mMesh,
        stk::topology::ELEMENT_RANK,
        phaseSupport.get_negative_levelset_block_selector(surfaceIdentifier) &
            meta.locally_owned_part() & this->get_aux_meta().active_part());
    all_reduce_sum(this->mMesh.parallel(), volume);
    return volume;
  }
};

typedef MultipleEdgeCrossingsDecompositionFixture<TetWithMultipleEdgeCrossings>
    MultipleEdgeCrossingsFixture;
typedef MultipleEdgeCrossingsDecompositionFixture<RightTet> RightTetMultipleEdgeCrossings;
typedef MultipleEdgeCrossingsDecompositionFixture<RightTri>
    RightTriMultipleEdgeCrossingsDecomposition;

} // namespace

// The part of the tet inside the cuboid is the part below y = 0.5 (all of the tet has x <= 0.5),
// which needs nodes where the edges from A and B to X and C leave the cuboid.
TEST_F(MultipleEdgeCrossingsFixture, edgesFromSurfaceNodesThatPassInside_insideVolumeExact)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Cuboid cuboid(stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., 1., 10.));
  const auto & v = meshSpec.nodeLocs;
  const double exactVolume = volume_below_plane_y(v[0], v[3], v[1], v[2], 0.5);
  EXPECT_NEAR(exactVolume,
      decompose_and_compute_inside_volume(cuboid),
      1.e-12 * tet_volume({{v[0], v[1], v[2], v[3]}}));
}

// The slab 0.4 <= x <= 0.6 crosses the edges of the right tet from node 1 twice each; cutting them
// makes leaf edges from the nodes on the slab to node 1, and to the nodes beyond it, which cross it
// again, and are cut in the next round.  The part of the tet in the slab is then exact:
// the integral of (1-x)^2/2 from 0.4 to 0.6.
TEST_F(RightTetMultipleEdgeCrossings, slabAcrossTet_insideVolumeExact)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Cuboid slab(stk::math::Vector3d(0.5, 0., 0.), stk::math::Vector3d(0.2, 10., 10.));
  const double exactVolume = (0.6 * 0.6 * 0.6 - 0.4 * 0.4 * 0.4) / 6.;
  EXPECT_NEAR(exactVolume, decompose_and_compute_inside_volume(slab), 1.e-12);
}

// In 2D: the band 0.4 <= x <= 0.6 crosses the edges of the right tri from node 1 twice each;
// cutting them makes leaf edges that cross it again, cut in the next round.  The part of the tri in
// the band is then exact: the integral of 1-x from 0.4 to 0.6.
TEST_F(RightTriMultipleEdgeCrossingsDecomposition, bandAcrossTri_insideAreaExact)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Cuboid band(stk::math::Vector3d(0.5, 0., 0.), stk::math::Vector3d(0.2, 10., 10.));
  const double exactArea = 0.5 * (0.6 * 0.6 - 0.4 * 0.4);
  EXPECT_NEAR(exactArea, decompose_and_compute_inside_volume(band), 1.e-12);
}

} // namespace krino
