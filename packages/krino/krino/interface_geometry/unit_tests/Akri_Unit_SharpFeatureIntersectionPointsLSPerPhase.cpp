// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <Akri_Unit_SharpFeatureIntersectionPointChecks.hpp>

namespace krino {

namespace {

constexpr double lsPerPhaseEdgeTol = 1.e-3;

// With one level set per phase, phase i where level set i is smallest.
template <typename MESHSPEC, unsigned NUM_LS>
class LSPerPhaseSharpFeaturePointsFixture
    : public SharpFeaturePointChecks<DecompositionFixture<MESHSPEC, LSPerPhasePolicy, NUM_LS>>
{
public:
  LSPerPhaseSharpFeaturePointsFixture()
  {
    this->set_valid_proc_sizes_for_test({1});
    this->build_mesh_from_spec(meshSpec);
    if (this->is_valid_proc_size_for_test())
      this->setup_ls_fields();
  }

  // Level set i is -scales[i] at the nodes in phase i, and 0 at the others.  So within a tet it is
  // -scales[i] times the sum of the barycentric weights of the nodes in phase i, and the points
  // where three or four phases meet are at weights set by the scales.
  void set_level_sets_by_node_phase(const std::vector<unsigned> & nodePhases, const std::vector<double> & scales)
  {
    for (size_t nodeIndex = 0; nodeIndex < nodePhases.size(); ++nodeIndex)
    {
      const stk::mesh::Entity node = this->get_assigned_node_for_index(nodeIndex);
      if (!this->mMesh.is_valid(node)) continue;
      for (unsigned ls = 0; ls < NUM_LS; ++ls)
        this->node_ls_field(this->get_ls_field(ls), node) = (nodePhases[nodeIndex] == ls) ? -scales[ls] : 0.;
    }
  }

  std::unique_ptr<InterfaceGeometry> build_geometry()
  {
    const stk::mesh::MetaData & meta = this->mMesh.mesh_meta_data();
    return create_levelset_geometry(meta.spatial_dimension(),
        this->get_aux_meta().active_part(),
        this->cdfem_support(),
        Phase_Support::get(meta),
        this->levelset_fields());
  }

protected:
  MESHSPEC meshSpec;
};

typedef LSPerPhaseSharpFeaturePointsFixture<TwoRightTets, 3> TwoRightTets3LSPerPhaseSharpFeaturePoints;
typedef LSPerPhaseSharpFeaturePointsFixture<RightTet, 4> RightTet4LSPerPhaseSharpFeaturePoints;

// The two tets share the face of nodes 0, 2, 4, which are in phases 0, 1, 2; their other nodes are in
// phase 0.
const std::vector<unsigned> twoRightTetsNodePhases{0, 0, 1, 0, 2};

constexpr int quadruplePointCodimension = 3;
constexpr int triplePointCodimension = 2;
}

// The three phases meet on the shared face at its centroid, and the line where they meet runs
// through each tet to the centroid of its outer face: both tets have the same point on the shared
// face.
TEST_F(TwoRightTets3LSPerPhaseSharpFeaturePoints, triplePointOnSharedFace_samePointInBothTets)
{
  if (is_valid_proc_size_for_test())
  {
    set_level_sets_by_node_phase(twoRightTetsNodePhases, {1., 1., 1.});
    const auto geometry = build_geometry();

    EXPECT_NEAR(1. / 3.,
        expect_same_point_on_in_elements(
            *geometry, {0, 1}, {0, 2, 4}, triplePointCodimension, lsPerPhaseEdgeTol),
        1.e-12);
    expect_same_point_on_in_elements(
        *geometry, {0}, {1, 2, 4}, triplePointCodimension, lsPerPhaseEdgeTol);
    expect_same_point_on_in_elements(
        *geometry, {1}, {2, 3, 4}, triplePointCodimension, lsPerPhaseEdgeTol);
    for (size_t elementIndex : {0u, 1u})
    {
      EXPECT_EQ(2u,
          element_points(*geometry, elementIndex, triplePointCodimension, lsPerPhaseEdgeTol).size())
          << "element " << elementIndex;
    }
  }
}

// Scaling phase 1 moves the points where the phases meet to within the tolerance of the edges
// opposite node 2 (weight 1/(1+2k) there): each is a point on such an edge, the one of the shared
// face on the edge of nodes 0 and 4, the same in both tets.
TEST_F(TwoRightTets3LSPerPhaseSharpFeaturePoints, triplePointsNearEdges_edgePoints)
{
  if (is_valid_proc_size_for_test())
  {
    const double k = 1666.;
    set_level_sets_by_node_phase(twoRightTetsNodePhases, {1., k, 1.});
    const auto geometry = build_geometry();

    EXPECT_NEAR(0.5,
        expect_same_point_on_in_elements(
            *geometry, {0, 1}, {0, 4}, triplePointCodimension, lsPerPhaseEdgeTol),
        1.e-12);
    expect_same_point_on_in_elements(
        *geometry, {0}, {1, 4}, triplePointCodimension, lsPerPhaseEdgeTol);
    expect_same_point_on_in_elements(
        *geometry, {1}, {3, 4}, triplePointCodimension, lsPerPhaseEdgeTol);
    for (size_t elementIndex : {0u, 1u})
    {
      EXPECT_EQ(2u,
          element_points(*geometry, elementIndex, triplePointCodimension, lsPerPhaseEdgeTol).size())
          << "element " << elementIndex;
    }
  }
}

// Each node of the tet in its own phase, phase 0 scaled so that the quadruple point has weight eta
// at node 0; the line where phases 1, 2, 3 meet runs from it to the centroid of the face opposite
// node 0, eta away.  Within the tolerance the quadruple point is on that face, at that triple
// point: that it captures it is decided when cutting, by the domains of its node.
TEST_F(
    RightTet4LSPerPhaseSharpFeaturePoints, quadruplePointWithinToleranceOfFace_quadruplePointOnFace)
{
  if (is_valid_proc_size_for_test())
  {
    for (const double eta : {1.e-4, 5.e-4})
    {
      set_level_sets_by_node_phase({0, 1, 2, 3}, {(1./eta - 1.)/3., 1., 1., 1.});
      const auto geometry = build_geometry();

      const std::vector<IntersectionPoint> quadruplePoints =
          element_points(*geometry, 0, quadruplePointCodimension, lsPerPhaseEdgeTol);
      ASSERT_EQ(1u, quadruplePoints.size()) << "eta " << eta;
      ASSERT_EQ(1u, points_on(quadruplePoints, {1, 2, 3}).size()) << "eta " << eta;
      for (const double weight : quadruplePoints[0].get_weights())
      {
        EXPECT_NEAR(1. / 3., weight, 1.e-12) << "eta " << eta;
      }

      const std::vector<IntersectionPoint> triplePoints =
          element_points(*geometry, 0, triplePointCodimension, lsPerPhaseEdgeTol);
      const std::vector<const IntersectionPoint *> triplePointsOnFace =
          points_on(triplePoints, {1, 2, 3});
      ASSERT_EQ(1u, triplePointsOnFace.size()) << "eta " << eta;
      for (const double weight : triplePointsOnFace[0]->get_weights())
      {
        EXPECT_NEAR(1./3., weight, 1.e-12) << "eta " << eta;
      }
    }
  }
}

TEST_F(RightTet4LSPerPhaseSharpFeaturePoints, quadruplePointBeyondToleranceOfFace_quadrupleAndTriplePoints)
{
  if (is_valid_proc_size_for_test())
  {
    for (const double eta : {1.5e-3, 1.e-2, 1.e-1})
    {
      set_level_sets_by_node_phase({0, 1, 2, 3}, {(1./eta - 1.)/3., 1., 1., 1.});
      const auto geometry = build_geometry();

      const std::vector<IntersectionPoint> quadruplePoints =
          element_points(*geometry, 0, quadruplePointCodimension, lsPerPhaseEdgeTol);
      ASSERT_EQ(1u, quadruplePoints.size()) << "eta " << eta;
      ASSERT_EQ(1u, points_on(quadruplePoints, {0, 1, 2, 3}).size()) << "eta " << eta;
      EXPECT_NEAR(eta, weight_of_node(quadruplePoints[0], 0), 1.e-12) << "eta " << eta;

      const std::vector<IntersectionPoint> triplePoints =
          element_points(*geometry, 0, triplePointCodimension, lsPerPhaseEdgeTol);
      const std::vector<const IntersectionPoint *> triplePointsOnFace =
          points_on(triplePoints, {1, 2, 3});
      ASSERT_EQ(1u, triplePointsOnFace.size()) << "eta " << eta;
      for (const double weight : triplePointsOnFace[0]->get_weights())
      {
        EXPECT_NEAR(1./3., weight, 1.e-12) << "eta " << eta;
      }
    }
  }
}

} // namespace krino
